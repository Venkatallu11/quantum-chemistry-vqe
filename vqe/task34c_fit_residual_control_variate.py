#!/usr/bin/env python3
"""
task34c_fit_residual_control_variate.py -- iteration 34, Task C. Control
variate for the champion pipeline, adapted from the proposal's section 18
(use a known physical identity that should hold exactly to reduce energy
variance) -- but the LITERAL version doesn't apply here: checked directly
against task29c_manifold_estimator.py's fit_pure_state, `a_hat = best_v /
np.linalg.norm(best_v)` HARD-ENFORCES Sum(a_j^2)=1 on every replicate, by
construction, with zero exception -- so C(a)=Sum(a_j^2)-1 is identically
0 always, carrying no information. Building anything on that would be
testing a null signal and calling it a result.

THE REAL AVAILABLE SIGNAL: `fit_pure_state` ALSO returns `best_val`, the
leftover weighted-squared-residual of the best-fit pure state against the
raw (possibly non-pure, noisy) measured data BEFORE normalization is
enforced. A bootstrap replicate whose raw noisy data looks less like a
genuine pure state should show a LARGER best_val -- a real, physically
motivated, measurable-before-any-hard-constraint quantity that plausibly
correlates with how bad that replicate's noise realization was, and
therefore with the resulting energy error. This tests exactly that
correlation, honestly, before assuming it's useful.

METHOD: reuses Task 32G/33C's own bootstrap worker discipline (resample
shots AND which of the 16 real literal-twirled MC draws contribute per
replicate, refit the manifold per slot, separate OS process per
replicate) on Task 31C's real champion-pipeline data, additionally
capturing SUM(best_val across all 21 kept slots) per replicate. Then:
  1. Real Pearson correlation of total_residual vs |E-E_exact| across
     N_BOOT replicates.
  2. If correlated, the MSE-optimal control-variate coefficient
     beta = Cov(E, residual) / Var(residual), applied as
     E_corrected = E - beta*(residual - mean(residual)), and the
     resulting variance reduction, real, measured, not assumed.

Run:
    python vqe/task34c_fit_residual_control_variate.py
"""
import os
import sys
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
N_BOOT = 48
EXACT_ENERGY_TARGET = 0.5
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task34c_fit_residual_control_variate_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def _boot_worker(p, non_id_labels, diag, kept, P_S, by_name_label, n_mc, seed):
    rng = np.random.default_rng(seed)
    a_by_name = {}
    total_residual = 0.0
    for name in kept:
        blended = {}
        for l in non_id_labels:
            key = (name, l)
            if key not in by_name_label:
                continue
            entries = by_name_label[key]
            draw_idx = rng.integers(0, len(entries), size=n_mc)
            vals_signed = []
            for idx in draw_idx:
                sign, gamma, counts = entries[idx]
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m = expectation_from_counts(resampled, l)
                vals_signed.append(sign * gamma * m)
            blended[l] = max(-1.0, min(1.0, float(np.mean(vals_signed))))
        v0 = target_coeff_vector(name, K)
        a_hat, best_val = fit_pure_state(P_S, blended, {l: 1.0 for l in blended}, K, v0,
                                          seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
        total_residual += best_val
    full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err, total_residual


def main():
    print("\n" + "=" * 96)
    print("  task34c_fit_residual_control_variate.py -- fit-residual control variate, champion pipeline")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(CKPT_PATH) as f:
        ck = json.load(f)
    tags = [tuple(t) for t in ck["tags"]]
    counts_list = ck["counts"]
    gamma_per_circuit = ck["gamma_per_circuit"]
    by_name_label = {}
    for (name, group, draw, sign), counts, gamma in zip(tags, counts_list, gamma_per_circuit):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma, counts))
    n_mc = min(len(v) for v in by_name_label.values())
    print(f"  loaded Task 31C's real checkpoint: {n_mc} MC draws/label, {len(kept)} slots")

    print(f"\n  -- running {N_BOOT} bootstrap replicates, capturing per-replicate energy error AND "
          f"total fit-residual (sum of best_val across all {len(kept)} slots) --")
    with ProcessPoolExecutor(max_workers=2) as ex:
        futures = [ex.submit(_boot_worker, p, non_id_labels, diag, kept, P_S, by_name_label, n_mc, 60_000 + i)
                   for i in range(N_BOOT)]
        results = [fut.result() for fut in futures]
    errs = np.array([r[0] for r in results])
    residuals = np.array([r[1] for r in results])

    print(f"\n  raw energy error: median={np.median(errs):.4f}  std={errs.std(ddof=1):.4f}  "
          f"MSE={np.median(errs)**2 + errs.std(ddof=1)**2:.4f}")
    print(f"  fit residual (sum best_val): median={np.median(residuals):.6f}  "
          f"range=[{residuals.min():.6f}, {residuals.max():.6f}]")

    rho = float(np.corrcoef(errs, residuals)[0, 1])
    print(f"\n  -- Pearson correlation(|E-E_exact|, total_fit_residual) across {N_BOOT} replicates: rho={rho:+.4f} --")

    bias_raw = float(np.median(errs))
    std_raw = float(errs.std(ddof=1))
    mse_raw = bias_raw ** 2 + std_raw ** 2

    if abs(rho) > 0.15:
        cov = float(np.cov(errs, residuals, ddof=1)[0, 1])
        var_res = float(residuals.var(ddof=1))
        beta = cov / var_res
        err_corrected = errs - beta * (residuals - residuals.mean())
        bias_cv = float(np.median(err_corrected))
        std_cv = float(err_corrected.std(ddof=1))
        mse_cv = bias_cv ** 2 + std_cv ** 2
        print(f"\n  -- APPLYING control variate: beta={beta:.4f} --")
        print(f"     RAW:              median={bias_raw:.4f}  std={std_raw:.4f}  MSE={mse_raw:.4f}")
        print(f"     CONTROL-VARIATE:  median={bias_cv:.4f}  std={std_cv:.4f}  MSE={mse_cv:.4f}")
        pct = 100 * (mse_raw - mse_cv) / mse_raw
        print(f"     MSE change: {pct:+.1f}%")
        verdict = ("Real, meaningful reduction" if pct > 10 else
                    "Small/marginal effect" if pct > 0 else
                    "No benefit -- control variate does not help despite nonzero correlation")
        print(f"     HONEST READ: {verdict}")
    else:
        print(f"\n  -- |rho|={abs(rho):.4f} is too weak to bother applying a control-variate correction "
              f"(threshold 0.15) -- the fit residual does NOT meaningfully predict energy error on real "
              f"data. This idea does not pay off for this pipeline, reported plainly rather than forcing "
              f"a weak correction that would mostly add noise.")
        beta = None
        mse_cv = None
        std_cv = None
        bias_cv = None

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N_BOOT": N_BOOT, "errs": errs.tolist(), "residuals": residuals.tolist(),
            "rho": rho, "bias_raw": bias_raw, "std_raw": std_raw, "mse_raw": mse_raw,
            "beta": beta, "bias_cv": bias_cv, "std_cv": std_cv, "mse_cv": mse_cv,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
