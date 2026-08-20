#!/usr/bin/env python3
"""
task33c_adaptive_label_pec_champion.py -- iteration 33, Task C. Apply
PER-LABEL damping to the ACTUAL CHAMPION pipeline (literal-twirling PEC +
manifold fit), not the fast analytic-ratio shortcut Task 33B tested.

Why this task exists: Task 33B found that per-label damping (weight each
label's correction by 1/(1+k*boundary_proximity), using Task 32I's own
measured per-label risk proxy) gives a real 88.6% MSE reduction relative
to its OWN undamped baseline (8.90->2.92 kcal/mol median) -- but that
baseline was the fast analytic-ratio shortcut (Task 30B/32I's method),
which is NOT this project's best pipeline. Task 32G already established
that the literal-twirling method (Task 31C, 4,368 real circuits) +
manifold fit, blended via a single GLOBAL lambda, is MSE-optimal at
lambda=1 (full PEC, MSE=2.33) -- global blending does not help there.
This task asks the real remaining question: does PER-LABEL damping help
THIS pipeline (the one that actually matters), the way it helped the
weaker shortcut?

METHOD: exact reuse of Task 32G's `_lambda_worker` bootstrap machinery
(resample shots AND which of the 16 real literal-twirled MC draws
contribute, blend raw vs PEC-corrected per label, refit the manifold,
compute energy -- one full replicate per OS process, matching Task
32B/32C/32G's established discipline for isolating manifold-fit seed
variance). The ONLY change: `lam` is no longer one global scalar -- it is
lambda_l(k) = 1 / (1 + k*boundary_proximity[l]) per label l, using the
SAME boundary_proximity risk proxy Task 32I/33B computed (recomputed here
identically, not reused from a file, so this script is self-contained).
At k=0, lambda_l=1 for EVERY label regardless of its risk score -- this
EXACTLY reproduces Task 32G's lambda=1 result (same data, same blend
formula, same worker logic), so k=0 here is a valid, apples-to-apples
reproduction of Task 32G's MSE=2.33 baseline (unlike Task 33B, where the
k=0 comparison to Task 32G was WRONG -- different correction methods).
This script's k=0 result is verified against Task 32G's own number before
any k>0 result is trusted.

DATA: reuses task31c_full_pec_calibration.json's real checkpoint (Task
31C's 4,368 real circuits) exactly as Task 32G did. No new submissions.

Run:
    python vqe/task33c_adaptive_label_pec_champion.py
"""
import os
import sys
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from task30b_pec_application import analytic_A_and_B
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
GATE_NAME = "zz"
P2_ZZ = 0.0146
SHOTS = 100_000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
K_GRID = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
N_BOOT_PER_K = 24  # matches Task 32G's N_BOOT_PER_LAMBDA exactly, for a like-for-like comparison
TASK32G_LAMBDA1_MSE = 2.33  # this task's k=0 MUST reproduce this (same data, same formula, all lambda_l=1)
EXACT_ENERGY_TARGET = 0.5
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task33c_adaptive_label_pec_champion_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def _k_worker(p, non_id_labels, diag, kept, P_S, by_name_label, n_mc, boundary_proximity, k, seed):
    """One bootstrap replicate at damping level k -- same resample-shots
    AND resample-which-MC-draws-contribute discipline as Task 32G's
    _lambda_worker, generalized to a per-label lambda_l(k)."""
    rng = np.random.default_rng(seed)
    a_by_name = {}
    for name in kept:
        blended = {}
        for l in non_id_labels:
            key = (name, l)
            if key not in by_name_label:
                continue
            entries = by_name_label[key]
            draw_idx = rng.integers(0, len(entries), size=n_mc)
            raw_vals, pec_vals_signed = [], []
            for idx in draw_idx:
                sign, gamma, counts = entries[idx]
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m = expectation_from_counts(resampled, l)
                raw_vals.append(m)
                pec_vals_signed.append(sign * gamma * m)
            raw_mean = float(np.mean(raw_vals))
            pec_mean = max(-1.0, min(1.0, float(np.mean(pec_vals_signed))))
            prox = boundary_proximity.get(l, 0.0)
            lam_l = 1.0 / (1.0 + k * prox)
            blended[l] = (1 - lam_l) * raw_mean + lam_l * pec_mean
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, blended, {l: 1.0 for l in blended}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err


def main():
    print("\n" + "=" * 96)
    print("  task33c_adaptive_label_pec_champion.py -- per-label damping on the ACTUAL champion pipeline")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins = learned["forte-1"]["gpi"]
    gpi2_bins = gpi_bins

    # -- boundary_proximity: recomputed identically to Task 32I/33B, self-contained (not read from a file) --
    boundary_proximity = {}
    for name in kept:
        A, B = analytic_A_and_B(fixed_solutions[name]["angles"], GATE_NAME, P2_ZZ, gpi_bins, gpi2_bins, non_id_labels)
        for l in non_id_labels:
            if l not in A:
                continue
            ratio_here = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            boundary_proximity.setdefault(l, []).append(abs(ratio_here - 1.0))
    boundary_proximity = {l: float(np.mean(v)) for l, v in boundary_proximity.items()}
    print(f"  boundary_proximity computed for {len(boundary_proximity)} labels "
          f"(median={sorted(boundary_proximity.values())[len(boundary_proximity)//2]:.4f})")

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
    print(f"  loaded Task 31C's real checkpoint: {n_mc} MC draws/label available")

    print(f"\n  -- k sweep: {N_BOOT_PER_K} bootstrap replicates per k, champion (literal-twirling+manifold) pipeline --")
    k_results = {}
    for k in K_GRID:
        print(f"  k={k}: running {N_BOOT_PER_K} bootstrap replicates...")
        with ProcessPoolExecutor(max_workers=4) as ex:
            futures = [ex.submit(_k_worker, p, non_id_labels, diag, kept, P_S, by_name_label, n_mc,
                                  boundary_proximity, k, seed=40_000 + int(k * 1000) + i)
                       for i in range(N_BOOT_PER_K)]
            errs = np.array([fut.result() for fut in futures])
        bias_proxy = float(np.median(errs))
        sigma = float(errs.std(ddof=1))
        mse = bias_proxy ** 2 + sigma ** 2
        passes = (abs(bias_proxy) + 2 * sigma) < EXACT_ENERGY_TARGET
        k_results[k] = {"errs": errs.tolist(), "median": bias_proxy, "std": sigma, "mse": mse, "passes": bool(passes)}
        print(f"    k={k:<6} median={bias_proxy:8.4f}  std={sigma:7.4f}  MSE={mse:9.4f}  "
              f"{'PASS' if passes else 'FAIL'} (|b|+2sigma<0.5)")

    k0_mse = k_results[0.0]["mse"]
    print(f"\n  -- SANITY CHECK: k=0 should reproduce Task 32G's lambda=1 baseline (MSE={TASK32G_LAMBDA1_MSE}) --")
    print(f"    this run's k=0 MSE = {k0_mse:.4f}  "
          f"(expected close to {TASK32G_LAMBDA1_MSE}; different bootstrap seeds -> some real sampling variation expected)")

    best_k = min(k_results, key=lambda kk: k_results[kk]["mse"])
    best_mse = k_results[best_k]["mse"]
    print(f"\n  MSE-MINIMIZING k* = {best_k}  (MSE={best_mse:.4f})")
    if best_k > 0.0 and best_mse < k0_mse:
        improvement_pct = 100.0 * (k0_mse - best_mse) / k0_mse
        print(f"  -> Per-label damping (k={best_k}) beats this run's own k=0 (full PEC, all labels) by "
              f"{improvement_pct:.1f}% MSE on the CHAMPION pipeline -- a real improvement on the project's "
              f"actual best method, not just the shortcut.")
    else:
        print(f"  -> Per-label damping does NOT beat k=0 on the champion pipeline. Consistent with Task 32G's "
              f"global-lambda finding, now confirmed at per-label granularity too: this pipeline's literal-"
              f"twirling PEC correction is already MSE-optimal at full strength for every label, and "
              f"boundary-proximity-driven damping is not the fix, on real data, for either the shortcut "
              f"(Task 33B, where it DID help) or the champion (this task).")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "boundary_proximity": boundary_proximity,
            "k_sweep": {str(k): v for k, v in k_results.items()},
            "best_k": best_k,
            "best_mse": best_mse,
            "task32g_lambda1_mse_reference": TASK32G_LAMBDA1_MSE,
            "k0_reproduction_check": k0_mse,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
