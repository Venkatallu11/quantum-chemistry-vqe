#!/usr/bin/env python3
"""
task3_family_wise_zne.py — iteration 26, Task 3. THE KEY METHODOLOGICAL
CHANGE: fit the fold response PER FAMILY, not globally, and select the
model class by predicting a HELD-OUT fold (9), never by which model
yields the best-looking final energy -- the rule this project's own
history shows every prior ZNE attempt violated (iterations 2, 11, 13:
the plateau test failed precisely because the fit was implicitly chosen
to look good).
============================================================================
MODEL CLASSES compared: linear, quadratic (np.polyfit), exponential
(A*exp(-k*fold)+E_inf, iteration-13-style), rational (Pade[1/1] via
scipy.optimize.curve_fit), stretched-exponential (A*exp(-(k*fold)^beta)+
E_inf), and a genuine Gaussian Process (scikit-learn, RBF+white-noise
kernel, installed this iteration) as the GP/spline entry -- a real GP
posterior, not a spline dressed up as one.

SELECTION RULE (non-negotiable, per explicit instruction): fit each
model class on folds [1,3,5], predict fold=9, score by |predicted -
actual_measured_fold9|. The model class with the LOWEST held-out error
is selected -- NEVER the one whose fold-0 extrapolation looks best. This
is applied independently at three granularities: GLOBAL (one model class
for the whole aggregate signal), FAMILY-WISE (one model class per
family, shared across that family's (slot,label) curves), and PER-TERM
(one model class per individual (slot,label) curve -- acknowledged as
data-poor with only 4 points, disclosed not hidden).

Once a model class is selected (by holdout only), the FINAL fold->0
extrapolation refits that class on ALL 4 available points (1,3,5,9) for
a lower-variance final estimate -- selection and final-fit are
deliberately separate steps, stated explicitly so this isn't quietly
peeking at the holdout point twice.

SCOPE CAVEAT, disclosed: Task 2's dataset covers 12 of 36 K=6 slots (all
6 diagonal + 6 representative cross-term pairs). The reconstructed
"energy" here is therefore a PARTIAL comparison: the 6 measured diagonal
+ 6 measured cross-term matrix elements are extrapolated per scheme; the
remaining 24 (of 30) cross-term slots are held at their EXACT IDEAL value
in every scheme compared (not raw-noisy, not zero) -- a neutral
placeholder chosen so every scheme is compared on the SAME partial
reconstruction, isolating the effect of the MEASURED terms' extrapolation
quality. This is NOT a claim about the full real forged energy.

Run:
    python vqe/task3_family_wise_zne.py
"""
import os
import sys
import json
import warnings
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz
from qiskit.quantum_info import Statevector, Pauli

K = 6
FOLD_FIT = [1, 3, 5]
FOLD_HOLDOUT = 9
DATASET_PATH = os.path.join(os.path.dirname(__file__), "task2_fold_response_dataset_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task3_family_wise_zne_results.json")

CDR_CITED = {  # reused, not re-derived -- iterations 9 and 19's own established real-hardware finding
    "aria-1": {"note": "2.1-2.6x WORSE than raw on real hardware, confirmed independently twice", "verdict": "actively harmful"},
    "forte-1": {"note": "2.1-2.6x WORSE than raw on real hardware, confirmed independently twice", "verdict": "actively harmful"},
}


# ---------------------------------------------------------------------------
# model classes
# ---------------------------------------------------------------------------

def fit_linear(folds, vals):
    c = np.polyfit(folds, vals, 1)
    return lambda x: np.polyval(c, x), c.tolist()


def fit_quadratic(folds, vals):
    if len(folds) < 3:
        return None, None
    c = np.polyfit(folds, vals, 2)
    return lambda x: np.polyval(c, x), c.tolist()


def fit_exponential(folds, vals):
    def model(x, A, k, Einf):
        return A * np.exp(-k * np.asarray(x)) + Einf
    try:
        p0 = [vals[0] - vals[-1], 0.1, vals[-1]]
        popt, _ = curve_fit(model, folds, vals, p0=p0, maxfev=5000)
        return lambda x: model(np.asarray(x), *popt), popt.tolist()
    except Exception:
        return None, None


def fit_rational(folds, vals):
    def model(x, a, b, c):
        x = np.asarray(x)
        return (a + b * x) / (1 + c * x)
    try:
        popt, _ = curve_fit(model, folds, vals, p0=[vals[0], 0.0, 0.1], maxfev=5000)
        return lambda x: model(np.asarray(x), *popt), popt.tolist()
    except Exception:
        return None, None


def fit_stretched_exp(folds, vals):
    def model(x, A, k, beta, Einf):
        x = np.asarray(x)
        return A * np.exp(-np.power(np.clip(k * x, 0, None), beta)) + Einf
    try:
        popt, _ = curve_fit(model, folds, vals, p0=[vals[0] - vals[-1], 0.1, 1.0, vals[-1]],
                             bounds=([-np.inf, 0, 0.1, -np.inf], [np.inf, np.inf, 5.0, np.inf]), maxfev=8000)
        return lambda x: model(np.asarray(x), *popt), popt.tolist()
    except Exception:
        return None, None


def fit_gp(folds, vals):
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
        X = np.array(folds).reshape(-1, 1)
        y = np.array(vals)
        kernel = ConstantKernel(1.0) * RBF(length_scale=3.0, length_scale_bounds=(0.5, 50)) + \
            WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-6, 1.0))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=3)
            gp.fit(X, y)
        return lambda x: gp.predict(np.asarray(x).reshape(-1, 1)), "gp_fitted"
    except Exception:
        return None, None


MODEL_CLASSES = {
    "linear": fit_linear, "quadratic": fit_quadratic, "exponential": fit_exponential,
    "rational": fit_rational, "stretched_exp": fit_stretched_exp, "gp": fit_gp,
}


def select_model_class(folds_fit, vals_fit, val_holdout):
    """SELECTION RULE: fit each class on folds_fit, predict FOLD_HOLDOUT,
    score by |predicted - val_holdout|. Returns (best_name, per_class_scores)."""
    scores = {}
    for name, fitter in MODEL_CLASSES.items():
        fn, params = fitter(folds_fit, vals_fit)
        if fn is None:
            scores[name] = {"holdout_error": float("inf"), "params": None}
            continue
        try:
            pred = float(np.atleast_1d(fn([FOLD_HOLDOUT]))[0])
            if not np.isfinite(pred):
                raise ValueError("non-finite prediction")
            scores[name] = {"holdout_error": abs(pred - val_holdout), "params": params, "predicted": pred}
        except Exception:
            scores[name] = {"holdout_error": float("inf"), "params": None}
    best_name = min(scores, key=lambda n: scores[n]["holdout_error"])
    return best_name, scores


def extrapolate_to_zero(model_name, folds_all, vals_all):
    """Refit the SELECTED class on ALL available points (1,3,5,9) for the
    final fold->0 estimate -- separate step from selection, per module
    docstring."""
    fn, _ = MODEL_CLASSES[model_name](folds_all, vals_all)
    if fn is None:
        return None
    try:
        val0 = float(np.atleast_1d(fn([0]))[0])
        return float(np.clip(val0, -1.0, 1.0))
    except Exception:
        return None


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task3_family_wise_zne.py -- family-wise ZNE, validated by predicting an unseen fold")
    print("=" * 96)

    with open(DATASET_PATH) as f:
        t2 = json.load(f)
    dataset = t2["dataset"]
    print(f"  loaded {len(dataset)} datapoints from Task 2")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]

    # exact ideal values for ALL 36 slots (for the "held at exact ideal" fallback)
    from qforge import fit_all_targets
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"])
    exact_ideal_all = {}
    for name, sol in fixed_solutions.items():
        sv = Statevector.from_instruction(build_ansatz(sol["angles"]))
        exact_ideal_all[name] = {l: float(sv.expectation_value(Pauli(l)).real) for l in non_id_labels}

    # organize: curves[(model, slot, label)] = {fold: measured}
    curves = {}
    families_of = {}
    for row in dataset:
        key = (row["model"], row["slot"], row["label"])
        curves.setdefault(key, {})[row["fold"]] = row["measured_expectation"]
        families_of[key] = row["family"]

    all_model_results = {}
    for model in ["aria-1", "forte-1"]:
        print(f"\n  ===== {model} =====")
        model_keys = [k for k in curves if k[0] == model]
        families = sorted(set(families_of[k] for k in model_keys))

        # -- GLOBAL selection: pool ALL curves' holdout errors, pick the class with lowest MEAN holdout error --
        global_scores = {name: [] for name in MODEL_CLASSES}
        for key in model_keys:
            c = curves[key]
            if not all(f in c for f in FOLD_FIT + [FOLD_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_FIT]
            _, scores = select_model_class(FOLD_FIT, vals_fit, c[FOLD_HOLDOUT])
            for name in MODEL_CLASSES:
                if np.isfinite(scores[name]["holdout_error"]):
                    global_scores[name].append(scores[name]["holdout_error"])
        global_best = min(global_scores, key=lambda n: np.mean(global_scores[n]) if global_scores[n] else float("inf"))
        print(f"    GLOBAL selected model class: {global_best} "
              f"(mean holdout error={np.mean(global_scores[global_best]):.4f})")

        # -- FAMILY-WISE selection: same, but per family --
        family_best = {}
        for fam in families:
            fam_keys = [k for k in model_keys if families_of[k] == fam]
            fam_scores = {name: [] for name in MODEL_CLASSES}
            for key in fam_keys:
                c = curves[key]
                if not all(f in c for f in FOLD_FIT + [FOLD_HOLDOUT]):
                    continue
                vals_fit = [c[f] for f in FOLD_FIT]
                _, scores = select_model_class(FOLD_FIT, vals_fit, c[FOLD_HOLDOUT])
                for name in MODEL_CLASSES:
                    if np.isfinite(scores[name]["holdout_error"]):
                        fam_scores[name].append(scores[name]["holdout_error"])
            best = min(fam_scores, key=lambda n: np.mean(fam_scores[n]) if fam_scores[n] else float("inf"))
            family_best[fam] = best
            print(f"    FAMILY '{fam}' selected model class: {best} "
                  f"(mean holdout error={np.mean(fam_scores[best]):.4f}, n={len(fam_keys)})")

        # -- PER-TERM selection: one model class per individual curve --
        per_term_best = {}
        for key in model_keys:
            c = curves[key]
            if not all(f in c for f in FOLD_FIT + [FOLD_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_FIT]
            best, _ = select_model_class(FOLD_FIT, vals_fit, c[FOLD_HOLDOUT])
            per_term_best[key] = best

        # -- reconstruct energy for each scheme: raw(fold1), global, family-wise, per-term --
        def build_raw(scheme):
            raw = {name: dict(exact_ideal_all[name]) for name in fixed_solutions}
            for key in model_keys:
                _, slot, label = key
                c = curves[key]
                if scheme == "raw_fold1":
                    if 1 in c:
                        raw[slot][label] = c[1]
                    continue
                if not all(f in c for f in FOLD_FIT + [FOLD_HOLDOUT]):
                    continue
                folds_all = sorted(c.keys())
                vals_all = [c[f] for f in folds_all]
                if scheme == "global":
                    cls = global_best
                elif scheme == "family_wise":
                    cls = family_best[families_of[key]]
                elif scheme == "per_term":
                    cls = per_term_best.get(key)
                    if cls is None:
                        continue
                val0 = extrapolate_to_zero(cls, folds_all, vals_all)
                if val0 is not None:
                    raw[slot][label] = val0
            return raw

        results_model = {}
        for scheme in ["raw_fold1", "global", "family_wise", "per_term"]:
            raw = build_raw(scheme)
            _, err = energy_and_err(p, raw, K)
            results_model[scheme] = err
            print(f"    {scheme:<12}: partial-reconstruction err_vs_exact = {err:.3f} kcal/mol")

        print(f"    CDR (cited, not re-run -- {CDR_CITED[model]['verdict']}: {CDR_CITED[model]['note']})")
        print(f"    family_wise + CDR: N/A -- CDR was not independently re-run against THIS dataset "
              f"(Task 2's circuits carry no CDR calibration data); composing an untested combination "
              f"would not be honest. See ALTERNATIVES NOT TAKEN.")

        all_model_results[model] = {
            "global_best_class": global_best, "family_best_class": family_best,
            "scheme_errors_kcal": results_model,
        }

    all_results = {
        "fold_fit": FOLD_FIT, "fold_holdout": FOLD_HOLDOUT,
        "aria-1": all_model_results["aria-1"], "forte-1": all_model_results["forte-1"],
        "cdr_cited": CDR_CITED,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return all_results


if __name__ == "__main__":
    main()
