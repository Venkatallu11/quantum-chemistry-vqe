#!/usr/bin/env python3
"""
task29c_manifold_estimator.py -- iteration 29, Task C. THE PHYSICAL-
MANIFOLD ESTIMATOR. The forged state is not an arbitrary 4-qubit state:
for K=6 the Schmidt basis is known exactly, every slot's target lives in
the 6-dimensional weight-2 sector, and with real amplitudes and unit
normalization it has only FIVE free parameters (not 35 for a general
density matrix in the sector, not hundreds for independent Pauli
channels).
============================================================================
PER SLOT: a_hat = argmin_a sum_l w_l [m_l - a^T P_S[l] a]^2, s.t. |a|=1.
Parametrized as a = v/|v| for unconstrained v in R^K -- automatically
unit-norm for any v != 0, no explicit constraint needed. Reuses
`build_P_S` UNCHANGED from `phys_constrained_reconstruction.py` (P_S[l] =
U^dagger P U, exact classical linear algebra, K=6 real gauge). Then
E = <psi(a_hat)|H|psi(a_hat)> directly -- normalization, particle number,
Hermiticity, and |<P>|<=1 hold BY CONSTRUCTION (a^T P_S a for unit a and
bounded P_S is automatically in-range), not by clipping.

MANIFOLD EXTRAPOLATION: fit each of the K real components a_j(lambda)
across folds {1,3,5,7,9} (the SAME two-stage held-out model-selection
machinery used everywhere else in this project), evaluate a_j(0), then
RENORMALIZE the resulting K-vector to unit norm (a per-component smooth
fit does not automatically stay on the unit sphere at an extrapolated
point). SIGN GAUGE, handled explicitly: a and -a are the SAME physical
state (a^T P_S a is invariant under global sign flip), so each fold's
independently-fitted a_hat is sign-aligned to fold=1's a_hat (flip if
dot product is negative) BEFORE any component-wise curve fitting -- an
unaligned trajectory would look discontinuous for a reason that has
nothing to do with the physics.

2q-only folding ONLY (`fold_native_2q`, the one ZNE variant that has
legitimately passed its own held-out test) -- reuses iteration 27's own
checkpoint (`task27c_full_h4_folds_K6.json`), NO new circuits, NO new
real submission.

COMPARISON BASELINE, computed fresh on the SAME checkpoint data for an
apples-to-apples number: standard per-Pauli-curve two-stage held-out ZNE
(iteration 27/28's own validated method, the "14.28 kcal/mol" family).

MANDATORY IDEAL-DATA SANITY CHECK: identical pipeline run on the `ideal`
model's own fold data.

Run:
    python vqe/task29c_manifold_estimator.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K, ckpt_path_for_K
from task27d_held_out_zne import (
    MODEL_CLASSES, fit_linear, fit_quadratic,
    FOLD_STAGE1_FIT, FOLD_STAGE1_HOLDOUT, FOLD_STAGE2_FIT, FOLD_STAGE2_HOLDOUT,
)
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

LINQUAD_ONLY = {"linear": fit_linear, "quadratic": fit_quadratic}  # Task 29G found this pair is only
                                                                    # 1.4x-2.4x worse-conditioned at
                                                                    # lambda=0 than lambda=9 (vs
                                                                    # millions-to-quintillions x for
                                                                    # rational/exponential) -- tested
                                                                    # here as the natural fix for
                                                                    # component-wise manifold extrapolation

K = 6
SHOTS = 100_000
N_SEEDS = 8
FOLD_FACTORS = [1, 3, 5, 7, 9]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task29c_manifold_estimator_results.json")


def target_coeff_vector(name, K):
    """The KNOWN ideal target, in the K-dim Schmidt-coefficient basis --
    used ONLY as an optimizer starting point (a well-motivated one, since
    we are fitting small physical deviations from a known target, not a
    blind search), never as a substitute for the measured data itself."""
    a = np.zeros(K)
    if name.startswith("u_"):
        n = int(name[2:])
        a[n] = 1.0
        return a
    # "(u{n}+u{m})" or "(u{n}-u{m})"
    inner = name.strip("()")
    sign = "+" if "+" in inner else "-"
    n_str, m_str = inner.split(sign)
    n, m = int(n_str[1:]), int(m_str[1:])
    a[n] = 1.0 / np.sqrt(2)
    a[m] = (1.0 if sign == "+" else -1.0) / np.sqrt(2)
    return a


def fit_pure_state(P_S, m_dict, w_dict, K, v0, seed):
    """Nonconvex (degree-4 polynomial on a sphere) -- a single restart is
    not reproducible run-to-run (verified: an earlier version of this
    function with 2 restarts gave 0.13 vs 0.65 kcal/mol for the SAME
    ideal-model fold=1 reconstruction on two independent runs, tracked
    down to landing in different local optima). 8 restarts (the known
    target + 7 perturbations at varied scales, deterministically seeded)
    fixes this -- disclosed and verified below, not assumed fixed."""
    labels = list(m_dict.keys())
    Ps = [P_S[l] for l in labels]
    ms = np.array([m_dict[l] for l in labels])
    ws = np.array([w_dict[l] for l in labels])

    def objective(v):
        a = v / np.linalg.norm(v)
        pred = np.array([float(np.real(a @ P @ a)) for P in Ps])
        return float(np.sum(ws * (pred - ms) ** 2))

    rng = np.random.default_rng(seed)
    inits = [v0] + [v0 + rng.normal(0, scale, K) for scale in [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2]]
    best_val, best_v = float("inf"), None
    for v_init in inits:
        res = minimize(objective, v_init, method="L-BFGS-B")
        if res.fun < best_val:
            best_val, best_v = res.fun, res.x
    a_hat = best_v / np.linalg.norm(best_v)
    return a_hat, best_val


def build_full_from_a(a_by_name, P_S, diag, K, non_id_labels):
    """Reconstructed <P_l> = a^T P_S[l] a for every kept slot's fitted
    pure state, then the SAME (u_n+u_m)/(u_n-u_m) synthesis formula used
    throughout this project (28G's own precedent) to get the off-diagonal
    slots from the reconstructed diag/plus values."""
    def expect(a, l):
        return float(np.real(a @ P_S[l] @ a))

    full = {}
    for name in diag:
        a = a_by_name[name]
        full[name] = {l: expect(a, l) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            a_pl = a_by_name[pl]
            full[pl] = {l: expect(a_pl, l) for l in non_id_labels}
            synth_minus = {}
            for l in non_id_labels:
                synth_minus[l] = full[un][l] + full[um][l] - full[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def zne_scalar_curve(fold_to_val, model_classes=None):
    """Two-stage held-out selection + extrapolation on a per-fold SCALAR,
    exactly Task 28G's `zne_on_curve` -- used here for the RAW per-Pauli
    baseline comparison's final energy assembly. `model_classes` defaults
    to the full 4-class set; pass a restricted dict (e.g. linear+quadratic
    only) to test whether narrowing the class set fixes the Task A/G
    extrapolation-conditioning pathology for THIS use (component-wise
    manifold extrapolation), without touching the baseline comparison."""
    mc = model_classes if model_classes is not None else MODEL_CLASSES
    if not all(f in fold_to_val for f in FOLD_STAGE1_FIT + [FOLD_STAGE1_HOLDOUT]):
        return None
    vals1 = [fold_to_val[f] for f in FOLD_STAGE1_FIT]
    b1n, b1e = None, float("inf")
    for name, fitter in mc.items():
        fn = fitter(FOLD_STAGE1_FIT, vals1)
        if fn is None:
            continue
        try:
            pred = float(np.atleast_1d(fn([FOLD_STAGE1_HOLDOUT]))[0])
            if np.isfinite(pred):
                err = abs(pred - fold_to_val[FOLD_STAGE1_HOLDOUT])
                if err < b1e:
                    b1e, b1n = err, name
        except Exception:
            continue
    vals2 = [fold_to_val[f] for f in FOLD_STAGE2_FIT]
    b2n, b2e = None, float("inf")
    for name, fitter in mc.items():
        fn = fitter(FOLD_STAGE2_FIT, vals2)
        if fn is None:
            continue
        try:
            pred = float(np.atleast_1d(fn([FOLD_STAGE2_HOLDOUT]))[0])
            if np.isfinite(pred):
                err = abs(pred - fold_to_val[FOLD_STAGE2_HOLDOUT])
                if err < b2e:
                    b2e, b2n = err, name
        except Exception:
            continue
    if b2n is None:
        return {"stage1_err": b1e, "stage2_err": None, "class": None, "val0": None, "excluded": True}
    all_folds = sorted(fold_to_val.keys())
    fn = mc[b2n](all_folds, [fold_to_val[f] for f in all_folds])
    val0, excluded = None, True
    if fn is not None:
        try:
            val0 = float(np.atleast_1d(fn([0]))[0])
            excluded = not np.isfinite(val0)
        except Exception:
            pass
    return {"stage1_err": b1e, "stage2_err": b2e, "class": b2n, "val0": val0, "excluded": excluded}


def main():
    print("\n" + "=" * 96)
    print("  task29c_manifold_estimator.py -- 5-parameter physical-manifold estimator, 2q-only folding")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    herm_err = max(float(np.max(np.abs(P_S[l] - P_S[l].conj().T))) for l in p["alpha_labels"])
    im_err = max(float(np.max(np.abs(P_S[l].imag))) for l in p["alpha_labels"])
    print(f"  K={K}: {len(kept)} kept slots. P_S Hermiticity err={herm_err:.2e}, "
          f"max imaginary part={im_err:.2e} (should be ~0 -- real gauge)")

    with open(ckpt_path_for_K(K)) as f:
        ck = json.load(f)
    print(f"  loaded 2q-only fold checkpoint (iteration 27, no new submission): folds={FOLD_FACTORS}")

    t_start = time.time()
    results_by_model = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        t0 = time.time()
        # -- per (fold, slot): seed-average bootstrap-resampled Pauli expectations, fit ONE pure state --
        a_by_fold = {fold: {} for fold in FOLD_FACTORS}
        raw_curves = {}   # (name, label) -> {fold: seed-averaged value} -- for the RAW per-Pauli baseline
        fit_residuals = []
        for fold in FOLD_FACTORS:
            key = f"{fold}|{model}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            per_name = {}
            for (name, group), counts in zip(tags, counts_list):
                per_name.setdefault(name, {}).setdefault(tuple(group), counts)

            for name in kept:
                seed_vals = {l: [] for l in non_id_labels if any(l in g for g in per_name[name])}
                seed_vars = {l: [] for l in seed_vals}
                for seed in range(N_SEEDS):
                    rng = np.random.default_rng(stable_seed("task29c", fold, model, name, seed))
                    for group_t, counts in per_name[name].items():
                        resampled = bootstrap_counts(counts, SHOTS, rng)
                        total = sum(resampled.values())
                        for l in group_t:
                            m = expectation_from_counts(resampled, l)
                            seed_vals[l].append(m)
                            var = max(1 - m ** 2, 1e-4) / max(total, 1)
                            seed_vars[l].append(var)
                            raw_curves.setdefault((name, l), {}).setdefault(fold, []).append(m)
                m_dict = {l: float(np.mean(v)) for l, v in seed_vals.items()}
                w_dict = {l: 1.0 / max(float(np.mean(seed_vars[l])), 1e-6) for l in seed_vals}
                v0 = target_coeff_vector(name, K)
                a_hat, resid = fit_pure_state(P_S, m_dict, w_dict, K, v0,
                                               seed=stable_seed("task29c_restart", fold, model, name))
                a_by_fold[fold][name] = a_hat
                fit_residuals.append(resid)

        for key2 in raw_curves:
            for fold in raw_curves[key2]:
                raw_curves[key2][fold] = float(np.mean(raw_curves[key2][fold]))

        # -- sign-align each slot's a_hat trajectory to fold=1 before extrapolating components --
        for name in kept:
            ref = a_by_fold[1][name]
            for fold in FOLD_FACTORS:
                if np.dot(a_by_fold[fold][name], ref) < 0:
                    a_by_fold[fold][name] = -a_by_fold[fold][name]

        # -- MANIFOLD extrapolation: fit each of K components across folds, evaluate a(0), renormalize --
        # two variants: the full 4-class set (matches production exactly) and a linear+quadratic-only
        # set (Task 29G's own evidence-based candidate fix for the extrapolation-conditioning pathology)
        a0_by_name = {}
        a0_linquad_by_name = {}
        n_component_excluded = 0
        n_component_excluded_linquad = 0
        for name in kept:
            comp_curves = np.array([a_by_fold[fold][name] for fold in FOLD_FACTORS])  # (5 folds, K)
            a0 = np.zeros(K)
            a0_lq = np.zeros(K)
            for j in range(K):
                fold_to_val = {fold: comp_curves[i, j] for i, fold in enumerate(FOLD_FACTORS)}
                r = zne_scalar_curve(fold_to_val)
                if r is None or r["val0"] is None or not np.isfinite(r["val0"]):
                    a0[j] = comp_curves[0, j]  # fallback: fold=1 raw component value
                    n_component_excluded += 1
                else:
                    a0[j] = r["val0"]
                r_lq = zne_scalar_curve(fold_to_val, model_classes=LINQUAD_ONLY)
                if r_lq is None or r_lq["val0"] is None or not np.isfinite(r_lq["val0"]):
                    a0_lq[j] = comp_curves[0, j]
                    n_component_excluded_linquad += 1
                else:
                    a0_lq[j] = r_lq["val0"]
            norm = np.linalg.norm(a0)
            a0_by_name[name] = a0 / norm if norm > 1e-9 else comp_curves[0]
            norm_lq = np.linalg.norm(a0_lq)
            a0_linquad_by_name[name] = a0_lq / norm_lq if norm_lq > 1e-9 else comp_curves[0]

        full_manifold = build_full_from_a(a0_by_name, P_S, diag, K, non_id_labels)
        _, err_manifold = energy_and_err(p, full_manifold, K)

        full_manifold_linquad = build_full_from_a(a0_linquad_by_name, P_S, diag, K, non_id_labels)
        _, err_manifold_linquad = energy_and_err(p, full_manifold_linquad, K)

        full_manifold_fold1 = build_full_from_a(a_by_fold[1], P_S, diag, K, non_id_labels)
        _, err_manifold_fold1 = energy_and_err(p, full_manifold_fold1, K)

        # -- RAW baseline, computed fresh on the SAME checkpoint: standard per-Pauli-curve ZNE --
        n_excluded_raw, n_total_raw, extrap_raw = 0, 0, {}
        for key2, c in raw_curves.items():
            if not all(f in c for f in FOLD_FACTORS):
                continue
            n_total_raw += 1
            r = zne_scalar_curve(c)
            if r is None or r["val0"] is None or not np.isfinite(r["val0"]) or abs(r["val0"]) > 1.0:
                n_excluded_raw += 1
                continue
            extrap_raw[key2] = r["val0"]
        full_raw = {name: {} for name in kept}
        for name in kept:
            for l in non_id_labels:
                key2 = (name, l)
                if key2 in extrap_raw:
                    full_raw[name][l] = extrap_raw[key2]
                elif key2 in raw_curves and 1 in raw_curves[key2]:
                    full_raw[name][l] = raw_curves[key2][1]
        full_raw_complete = {n2: dict(full_raw[n2]) for n2 in diag}
        for n in range(K):
            for m in range(K):
                if n >= m:
                    continue
                un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
                full_raw_complete[pl] = dict(full_raw[pl])
                synth_minus = {}
                for l in non_id_labels:
                    if l not in full_raw[pl] or l not in full_raw_complete[un] or l not in full_raw_complete[um]:
                        continue
                    synth_minus[l] = full_raw_complete[un][l] + full_raw_complete[um][l] - full_raw[pl][l]
                full_raw_complete[f"(u{n}-u{m})"] = synth_minus
        _, err_raw_zne = energy_and_err(p, full_raw_complete, K)

        raw1 = {name: {} for name in kept}
        for name in kept:
            for l in non_id_labels:
                key2 = (name, l)
                if key2 in raw_curves and 1 in raw_curves[key2]:
                    raw1[name][l] = raw_curves[key2][1]
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
                    synth_minus[l] = raw1_complete[un][l] + raw1_complete[um][l] - raw1[pl][l]
                raw1_complete[f"(u{n}-u{m})"] = synth_minus
        _, err_raw1 = energy_and_err(p, raw1_complete, K)

        t_elapsed = time.time() - t0
        print(f"\n  {model} ({t_elapsed:.1f}s, mean pure-state fit residual={np.mean(fit_residuals):.2e}, "
              f"{n_component_excluded}/{len(kept)*K} components fell back (4-class), "
              f"{n_component_excluded_linquad}/{len(kept)*K} fell back (linquad-only)):")
        print(f"    raw(fold=1)                       = {err_raw1:.3f} kcal/mol")
        print(f"    RAW per-Pauli-curve ZNE(0)        = {err_raw_zne:.3f} kcal/mol  "
              f"(excluded {n_excluded_raw}/{n_total_raw} curves)")
        print(f"    MANIFOLD, fold=1 only              = {err_manifold_fold1:.3f} kcal/mol")
        print(f"    MANIFOLD, extrap a(0), 4-class     = {err_manifold:.3f} kcal/mol")
        print(f"    MANIFOLD, extrap a(0), linquad-only = {err_manifold_linquad:.3f} kcal/mol")

        results_by_model[model] = {
            "err_raw_fold1": err_raw1, "err_raw_pauli_curve_zne": err_raw_zne,
            "n_excluded_raw_curves": n_excluded_raw, "n_total_raw_curves": n_total_raw,
            "err_manifold_fold1_only": err_manifold_fold1,
            "err_manifold_extrapolated_4class": err_manifold,
            "err_manifold_extrapolated_linquad": err_manifold_linquad,
            "n_component_fallback_4class": n_component_excluded,
            "n_component_fallback_linquad": n_component_excluded_linquad,
            "n_component_total": len(kept) * K,
            "mean_pure_state_fit_residual": float(np.mean(fit_residuals)),
        }

    print(f"\n  -- SUMMARY (total wall clock {time.time()-t_start:.1f}s) --")
    for model in ["ideal", "aria-1", "forte-1"]:
        r = results_by_model[model]
        note = "  (SANITY CHECK: manifold must not make the ideal control worse)" if model == "ideal" else ""
        print(f"    {model}: raw(fold1)={r['err_raw_fold1']:.2f}  RAW-ZNE={r['err_raw_pauli_curve_zne']:.2f}  "
              f"MANIFOLD(fold1)={r['err_manifold_fold1_only']:.2f}  "
              f"MANIFOLD-extrap-4class(0)={r['err_manifold_extrapolated_4class']:.2f}  "
              f"MANIFOLD-extrap-linquad(0)={r['err_manifold_extrapolated_linquad']:.2f}{note}")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results_by_model, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results_by_model


if __name__ == "__main__":
    main()
