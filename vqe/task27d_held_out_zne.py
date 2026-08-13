#!/usr/bin/env python3
"""
task27d_held_out_zne.py — iteration 27, Task D. Held-out ZNE validation,
no rescue allowed. Train on folds 1,3,5, predict 7. Then train on
1,3,5,7, predict 9. Compare linear/quadratic/exponential/rational. Keep
ONLY models that predict the held-out fold. NEVER choose a model by the
final energy it produces.
============================================================================
NO CLIPPING. NO HIDDEN REGULARIZATION. Iteration 26 Task 3 caught
per-term ZNE silently clipping unphysical extrapolations (values like
-101,291 for a quantity bounded in [-1,1]) to fake a good-looking
13.81 kcal/mol result. This file does the opposite on purpose: any curve
whose CHOSEN (held-out-validated) model extrapolates to |value| > 1 at
fold=0 is EXCLUDED from the final reconstruction and COUNTED, not
clipped and hidden. Excluded terms fall back to their RAW fold=1
measured value (the actual data, not an idealized placeholder) --
disclosed explicitly as a real limitation of what "ZNE energy" means once
some terms can't be trusted.

TWO-STAGE HELD-OUT PROCEDURE, exactly as specified:
  Stage 1: fit on folds [1,3,5], predict fold=7, score by |predicted -
    actual measured fold=7|.
  Stage 2: fit on folds [1,3,5,7], predict fold=9, score the SAME way.
A model class only "passes" for a given curve if ITS holdout error is
the lowest among the four classes AND is not catastrophically large
(no arbitrary pass threshold invented here -- "lowest of the four,
reported honestly" is the selection rule exactly as specified; whether
that's "good enough" is judged at the aggregate PASS-gate level, Task
27's own criterion 4/5, not per-curve).

Run:
    python vqe/task27d_held_out_zne.py
"""
import os
import sys
import json
import warnings
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
import ef_fragment as effrag_mod
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K_TARGET = 6   # this project's current main line; K=5 handled as a secondary comparison, see main()
FOLD_STAGE1_FIT = [1, 3, 5]
FOLD_STAGE1_HOLDOUT = 7
FOLD_STAGE2_FIT = [1, 3, 5, 7]
FOLD_STAGE2_HOLDOUT = 9
SHOTS = 100_000
N_SEEDS = 8
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")


def ckpt_path_for_K(K):
    return os.path.join(CKPT_DIR, f"task27c_full_h4_folds_K{K}.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task27d_held_out_zne_results.json")


def fit_linear(folds, vals):
    c = np.polyfit(folds, vals, 1)
    return lambda x: np.polyval(c, x)


def fit_quadratic(folds, vals):
    if len(folds) < 3:
        return None
    c = np.polyfit(folds, vals, 2)
    return lambda x: np.polyval(c, x)


def fit_exponential(folds, vals):
    def model(x, A, k, Einf):
        return A * np.exp(-k * np.asarray(x)) + Einf
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _ = curve_fit(model, folds, vals, p0=[vals[0] - vals[-1], 0.1, vals[-1]], maxfev=5000)
        return lambda x: model(np.asarray(x), *popt)
    except Exception:
        return None


def fit_rational(folds, vals):
    def model(x, a, b, c):
        x = np.asarray(x)
        return (a + b * x) / (1 + c * x)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _ = curve_fit(model, folds, vals, p0=[vals[0], 0.0, 0.1], maxfev=5000)
        return lambda x: model(np.asarray(x), *popt)
    except Exception:
        return None


MODEL_CLASSES = {"linear": fit_linear, "quadratic": fit_quadratic, "exponential": fit_exponential, "rational": fit_rational}


def select_and_predict(folds_fit, vals_fit, holdout_fold, target_fold):
    """Fits every class on folds_fit, scores by |predict(holdout_fold) -
    ACTUAL measured value at holdout_fold| (passed in via vals_fit's
    caller), returns (best_class_name, best_holdout_error, prediction_at_target_fold)."""
    pass  # implemented inline in main() where actual holdout value is available


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def run_for_K(K):
    print(f"\n  ===== K={K} =====")
    with open(ckpt_path_for_K(K)) as f:
        ck = json.load(f)
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=(K == 6))
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    group_idx = {}
    for gi, g in enumerate(groups):
        for l in g:
            group_idx[l] = gi
    diag, plus, kept = kept_slots_for_K(K)

    results_by_model = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        print(f"\n    --- {model} ---")
        # exact per-(slot,label) curves at each fold (mean over 8 bootstrap seeds, no shot-level per-seed fitting --
        # ZNE is fit on the MEAN measured curve, matching standard practice and this project's own prior convention)
        curves = {}   # (slot,label) -> {fold: value}
        for fold in [1, 3, 5, 7, 9]:
            key = f"{fold}|{model}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            per_name = {}
            for (name, group), counts in zip(tags, counts_list):
                per_name.setdefault(name, {}).setdefault(tuple(group), counts)
            for name in kept:
                for group_t, counts in per_name[name].items():
                    seed_vals = []
                    for seed in range(N_SEEDS):
                        rng = np.random.default_rng(stable_seed("task27d", K, fold, model, seed, name))
                        resampled = bootstrap_counts(counts, SHOTS, rng)
                        for l in group_t:
                            curves.setdefault((name, l), {}).setdefault(fold, []).append(expectation_from_counts(resampled, l))
        # collapse per-seed lists to means
        for key2 in curves:
            for fold in curves[key2]:
                curves[key2][fold] = float(np.mean(curves[key2][fold]))

        # -- STAGE 1: fit [1,3,5], predict 7 --
        stage1_selected = {}
        stage1_errors = {}
        for key2, c in curves.items():
            if not all(f in c for f in FOLD_STAGE1_FIT + [FOLD_STAGE1_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_STAGE1_FIT]
            best_name, best_err = None, float("inf")
            for name, fitter in MODEL_CLASSES.items():
                fn = fitter(FOLD_STAGE1_FIT, vals_fit)
                if fn is None:
                    continue
                try:
                    pred = float(np.atleast_1d(fn([FOLD_STAGE1_HOLDOUT]))[0])
                    if not np.isfinite(pred):
                        continue
                    err = abs(pred - c[FOLD_STAGE1_HOLDOUT])
                    if err < best_err:
                        best_err, best_name = err, name
                except Exception:
                    continue
            if best_name is not None:
                stage1_selected[key2] = best_name
                stage1_errors[key2] = best_err
        mean_stage1_err = float(np.mean(list(stage1_errors.values()))) if stage1_errors else float("inf")
        print(f"      STAGE 1 (fit 1,3,5 -> predict 7): {len(stage1_selected)}/{len(curves)} curves fit; "
              f"mean held-out error={mean_stage1_err:.4f}")

        # -- STAGE 2: fit [1,3,5,7], predict 9 --
        stage2_selected = {}
        stage2_errors = {}
        for key2, c in curves.items():
            if not all(f in c for f in FOLD_STAGE2_FIT + [FOLD_STAGE2_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_STAGE2_FIT]
            best_name, best_err = None, float("inf")
            for name, fitter in MODEL_CLASSES.items():
                fn = fitter(FOLD_STAGE2_FIT, vals_fit)
                if fn is None:
                    continue
                try:
                    pred = float(np.atleast_1d(fn([FOLD_STAGE2_HOLDOUT]))[0])
                    if not np.isfinite(pred):
                        continue
                    err = abs(pred - c[FOLD_STAGE2_HOLDOUT])
                    if err < best_err:
                        best_err, best_name = err, name
                except Exception:
                    continue
            if best_name is not None:
                stage2_selected[key2] = best_name
                stage2_errors[key2] = best_err
        mean_stage2_err = float(np.mean(list(stage2_errors.values()))) if stage2_errors else float("inf")
        print(f"      STAGE 2 (fit 1,3,5,7 -> predict 9): {len(stage2_selected)}/{len(curves)} curves fit; "
              f"mean held-out error={mean_stage2_err:.4f}")

        # -- final extrapolation to fold=0, using stage-2 selected classes on ALL 5 points, NO CLIPPING --
        n_excluded = 0
        n_total = 0
        extrap = {}
        for key2, c in curves.items():
            if key2 not in stage2_selected:
                continue
            n_total += 1
            cls = stage2_selected[key2]
            all_folds = sorted(c.keys())
            vals_all = [c[f] for f in all_folds]
            fn = MODEL_CLASSES[cls](all_folds, vals_all)
            if fn is None:
                n_excluded += 1
                continue
            try:
                val0 = float(np.atleast_1d(fn([0]))[0])
            except Exception:
                n_excluded += 1
                continue
            if not np.isfinite(val0) or abs(val0) > 1.0:
                n_excluded += 1
                continue   # EXCLUDED, not clipped -- per explicit instruction
            extrap[key2] = val0
        print(f"      fold->0 extrapolation: {n_total - n_excluded}/{n_total} curves gave a PHYSICAL "
              f"(|value|<=1) result; {n_excluded} EXCLUDED (unphysical, not clipped)")

        # -- reconstruct energy: extrapolated value where available, else RAW fold=1 (real data, not placeholder) --
        full = {name: {} for name in kept}
        for name in kept:
            for l in non_id_labels:
                key2 = (name, l)
                if key2 in extrap:
                    full[name][l] = extrap[key2]
                elif key2 in curves and 1 in curves[key2]:
                    full[name][l] = curves[key2][1]
        # algebraic derivation of the "-" pairs, same as Task C
        full_complete = {n2: dict(full[n2]) for n2 in diag}
        for n in range(K):
            for m in range(K):
                if n >= m:
                    continue
                un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
                full_complete[pl] = dict(full[pl])
                synth_minus = {}
                for l in non_id_labels:
                    if l not in full[pl] or l not in full_complete[un] or l not in full_complete[um]:
                        continue
                    cross = full[pl][l] - (full_complete[un][l] + full_complete[um][l]) / 2
                    synth_minus[l] = full[pl][l] - 2 * cross
                full_complete[f"(u{n}-u{m})"] = synth_minus
        _, err_zne = energy_and_err(p, full_complete, K)

        # -- compare: raw fold=1 energy, for reference --
        raw1 = {name: {} for name in kept}
        for name in kept:
            for l in non_id_labels:
                key2 = (name, l)
                if key2 in curves and 1 in curves[key2]:
                    raw1[name][l] = curves[key2][1]
        raw1_complete = {n2: dict(raw1[n2]) for n2 in diag}
        for n in range(K):
            for m in range(K):
                if n >= m:
                    continue
                un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
                raw1_complete[pl] = dict(raw1[pl])
                synth_minus = {}
                for l in non_id_labels:
                    if l not in raw1[pl] or l not in raw1_complete[un] or l not in raw1_complete[um]:
                        continue
                    cross = raw1[pl][l] - (raw1_complete[un][l] + raw1_complete[um][l]) / 2
                    synth_minus[l] = raw1[pl][l] - 2 * cross
                raw1_complete[f"(u{n}-u{m})"] = synth_minus
        _, err_raw1 = energy_and_err(p, raw1_complete, K)

        improves = err_zne < err_raw1
        print(f"      raw (fold=1) energy error: {err_raw1:.2f} kcal/mol")
        print(f"      ZNE (fold->0, excluded-not-clipped) energy error: {err_zne:.2f} kcal/mol")
        print(f"      ZNE improves on raw: {improves}")

        results_by_model[model] = {
            "stage1_n_fit": len(stage1_selected), "stage1_mean_holdout_err": mean_stage1_err,
            "stage2_n_fit": len(stage2_selected), "stage2_mean_holdout_err": mean_stage2_err,
            "n_curves_total": n_total, "n_excluded_unphysical": n_excluded,
            "err_raw_fold1": err_raw1, "err_zne_fold0": err_zne, "zne_improves_on_raw": bool(improves),
            "class_selection_counts_stage2": {c: list(stage2_selected.values()).count(c) for c in MODEL_CLASSES},
        }

    return results_by_model


def main():
    print("\n" + "=" * 96)
    print("  task27d_held_out_zne.py -- held-out ZNE validation, no rescue")
    print("=" * 96)

    all_results = {}
    for K in [6, 5]:
        all_results[K] = run_for_K(K)

    print(f"\n  -- PASS-GATE CRITERIA 4 and 5 (ZNE predicts held-out folds; ZNE improves held-out error) --")
    for K, by_model in all_results.items():
        for model, r in by_model.items():
            crit4 = r["stage1_mean_holdout_err"] < 0.1 and r["stage2_mean_holdout_err"] < 0.1  # in Pauli-expectation units
            crit5 = r["zne_improves_on_raw"]
            print(f"    K={K} {model}: criterion4(predicts held-out, err<0.1)={crit4}  "
                  f"criterion5(improves raw)={crit5}  {'BOTH PASS' if crit4 and crit5 else 'FAILS at least one'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return all_results


if __name__ == "__main__":
    main()
