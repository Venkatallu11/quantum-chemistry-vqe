#!/usr/bin/env python3
"""
task32g_mse_reframing.py -- iteration 32, Task G. MINIMIZE MSE, NOT BIAS.
The target is not minimum |bias| but MSE = bias^2 + variance, with the
real acceptance criterion |b| + 2*sigma < 0.5. A slightly biased, STABLE
estimator can beat an unbiased one with catastrophic variance -- exactly
what this iteration's own numbers already show (Task 31C's 0.115 kcal/mol
"best ever" point estimate turned out to be a high-variance outlier per
Task 31F/32B; the pipeline's real submission-to-submission std is 0.264).
============================================================================
PART 1 -- MSE LEAGUE TABLE: every candidate estimator this iteration
produced, ranked by MSE = bias^2 + variance (using each estimator's own
measured bias and variance, not re-derived), against the |b|+2*sigma<0.5
criterion.

PART 2 -- ROBUST QUASI-PROBABILITY DECOMPOSITION, a concrete, testable
version: a SHRINKAGE estimator m_lambda = (1-lambda)*raw + lambda*PEC,
interpolating between raw (lambda=0, gamma=1, zero sampling overhead) and
full PEC (lambda=1). Reduces the quasi-probability correction's shot-noise
AMPLIFICATION (which scales with gamma^2) at the cost of reintroducing some
of PEC's bias-correction benefit. Swept over lambda in [0,1], MSE
(estimated via the SAME per-shot bootstrap machinery as Task 32B) is
computed at each point on Task 31C's real checkpointed data -- if the
MSE-minimizing lambda* < 1, PEC's "insist on an exact inverse" premise is
suboptimal for THIS circuit/noise regime, a real, actionable finding.

Run:
    python vqe/task32g_mse_reframing.py
"""
import os
import sys
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from qforge import combine_matrices, energy_from_alpha_matrices
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
LAMBDAS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
N_BOOT_PER_LAMBDA = 24
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32g_mse_reframing_results.json")
EXACT_ENERGY_TARGET = 0.5  # chemical-accuracy-adjacent target this project's acceptance criterion uses


def mse_league_table():
    """Every candidate estimator this iteration produced, with its OWN
    measured bias (point estimate, taken as-is) and variance (measured
    where available; None where only a single realization exists, which
    is itself the point -- an unmeasured variance is not a small one)."""
    table = [
        {"name": "Task 31C literal-PEC (N_MC=16, single realization)", "bias": 0.115, "sigma": None,
         "note": "no variance ever measured for this exact number -- Task 31F/32B show real spread is large"},
        {"name": "Task 31D covariance-aware manifold (FAIL)", "bias": 16.987, "sigma": None, "note": "disqualified"},
        {"name": "Task 31D uniform-weighted manifold, 32-seed bootstrap", "bias": 0.317, "sigma": None,
         "note": "shot-noise-only bootstrap, not full pipeline variance"},
        {"name": "Task 31F mean of 8 independent real submissions", "bias": 0.438, "sigma": 0.264,
         "note": "the most honest single number this project has: real submission-to-submission spread included"},
        {"name": "Task 32C convex SDP relaxation", "bias": 1.486, "sigma": 0.0,
         "note": "PERFECTLY reproducible (sigma=0 across processes) but worse bias -- see MSE below"},
        {"name": "Task 32D multinomial MLE + PEC-ratio (median of 8 processes)", "bias": 32.953, "sigma": None,
         "note": "range 16.8-66.3 across 8 processes -- reported median, real spread enormous"},
    ]
    for row in table:
        if row["sigma"] is not None:
            row["mse"] = row["bias"] ** 2 + row["sigma"] ** 2
            row["passes_acceptance"] = (abs(row["bias"]) + 2 * row["sigma"]) < EXACT_ENERGY_TARGET
        else:
            row["mse"] = None
            row["passes_acceptance"] = None
    return table


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def _lambda_worker(p, non_id_labels, diag, kept, P_S, by_name_label, n_mc, lam, seed):
    """ONE bootstrap draw at shrinkage level lambda: resample shots AND
    which MC draws contribute (same joint bootstrap as Task 32J), blend
    raw and PEC-corrected per label by lambda, refit manifold, compute
    energy. Separate process (Task 32B/32C's established discipline)."""
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
                raw_vals.append(m)  # unsigned, ungamma'd -- the RAW (lambda=0) measurement
                pec_vals_signed.append(sign * gamma * m)
            raw_mean = float(np.mean(raw_vals))
            pec_mean = max(-1.0, min(1.0, float(np.mean(pec_vals_signed))))
            blended[l] = (1 - lam) * raw_mean + lam * pec_mean
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, blended, {l: 1.0 for l in blended}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err


def main():
    print("\n" + "=" * 96)
    print("  task32g_mse_reframing.py -- MSE league table + shrinkage-PEC lambda sweep")
    print("=" * 96)

    print("\n  -- PART 1: MSE LEAGUE TABLE (this iteration's own candidate estimators) --")
    table = mse_league_table()
    for row in sorted(table, key=lambda r: (r["mse"] is None, r["mse"] if r["mse"] is not None else 0)):
        mse_str = f"{row['mse']:.3f}" if row["mse"] is not None else "UNMEASURED (treat as disqualifying, not small)"
        pass_str = ("PASS" if row["passes_acceptance"] else "FAIL") if row["passes_acceptance"] is not None else "N/A"
        print(f"    {row['name']:<58} bias={row['bias']:<8.3f} sigma={str(row['sigma']):<6} MSE={mse_str:<10} {pass_str}")
        print(f"        note: {row['note']}")

    print("\n  -- PART 2: SHRINKAGE-PEC lambda sweep (real data, joint bootstrap MSE at each lambda) --")
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
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

    lambda_results = {}
    for lam in LAMBDAS:
        print(f"  lambda={lam}: running {N_BOOT_PER_LAMBDA} bootstrap replicates...")
        with ProcessPoolExecutor(max_workers=4) as ex:
            futures = [ex.submit(_lambda_worker, p, non_id_labels, diag, kept, P_S, by_name_label, n_mc,
                                  lam, seed=20_000 + int(lam * 1000) + i)
                       for i in range(N_BOOT_PER_LAMBDA)]
            errs = np.array([fut.result() for fut in futures])
        bias_proxy = float(np.median(errs))  # median as a robust central estimate given known heavy tails
        sigma = float(errs.std(ddof=1))
        mse = bias_proxy ** 2 + sigma ** 2
        lambda_results[lam] = {"errs": errs.tolist(), "median": bias_proxy, "std": sigma, "mse": mse}
        print(f"    lambda={lam}: median={bias_proxy:.4f}  std={sigma:.4f}  MSE={mse:.4f}")

    best_lambda = min(lambda_results, key=lambda l: lambda_results[l]["mse"])
    print(f"\n  MSE-MINIMIZING lambda* = {best_lambda} (MSE={lambda_results[best_lambda]['mse']:.4f})")
    if best_lambda < 1.0:
        print(f"  -> FULL PEC (lambda=1) is NOT MSE-optimal for this circuit/data -- a partial correction "
              f"trades some bias for materially less shot-noise amplification, a real actionable finding.")
    else:
        print(f"  -> Full PEC (lambda=1) IS MSE-optimal here -- the exact-inverse premise holds for this "
              f"circuit/noise regime, at least under this bootstrap's joint shot+MC-draw resampling.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"mse_league_table": table, "lambda_sweep": lambda_results, "best_lambda": best_lambda},
                   f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
