#!/usr/bin/env python3
"""
task32j_bootstrap_estimator.py -- iteration 32, Task J. BOOTSTRAP THE WHOLE
ESTIMATOR END TO END. Task 31D's covariance-aware manifold fit failed
because a plug-in covariance-matrix INVERSE with n_seeds=32 < n_labels=36
sits in the classic small-sample regime where the off-diagonal structure
is estimation noise, not real correlation -- more sophisticated covariance
estimation is the wrong fix for that. Instead: resample the RAW COUNTS
(shot noise) AND which of the 16 real PEC Monte Carlo draws are used
(quasi-probability sampling noise) AND re-run the manifold fit (letting its
own real cross-process degenerate-optimum noise, quantified in Task 32B, be
part of the picture) TOGETHER, per bootstrap replicate, end to end. The
OBSERVED distribution of final energies across replicates IS the
uncertainty -- no covariance matrix, no matrix inverse, anywhere.
============================================================================
Reuses Task 31C's already-collected 4,368 real circuits' raw counts
unchanged (no new submission) -- same data Task 31C/31D/32B/32C all
reprocess, now bootstrapped jointly rather than one source at a time.

Each bootstrap replicate is a SEPARATE OS process (Task 32B/32C's own
established discipline: the manifold fit's cross-process degenerate-optimum
noise cannot be seen from seed variation within one process).

Run:
    python vqe/task32j_bootstrap_estimator.py
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
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
N_BOOTSTRAP = 200
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32j_bootstrap_estimator_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def _one_bootstrap_replicate(p, non_id_labels, diag, kept, P_S, by_name_label, n_mc, seed):
    """Picklable top-level worker: ONE full bootstrap replicate --
    resample shots for every draw, resample WHICH draws (with replacement)
    contribute, refit the manifold, compute energy. Runs in a separate
    process (see module docstring for why)."""
    rng = np.random.default_rng(seed)
    a_by_name = {}
    for name in kept:
        pec_dict = {}
        for l in non_id_labels:
            key = (name, l)
            if key not in by_name_label:
                continue
            entries = by_name_label[key]
            draw_idx = rng.integers(0, len(entries), size=n_mc)  # WITH replacement -- MC resampling
            vals_signed = []
            for idx in draw_idx:
                sign, gamma, counts = entries[idx]
                resampled = bootstrap_counts(counts, SHOTS, rng)  # shot-noise resampling
                m = expectation_from_counts(resampled, l)
                vals_signed.append(sign * gamma * m)
            pec_dict[l] = max(-1.0, min(1.0, float(np.mean(vals_signed))))
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, pec_dict, {l: 1.0 for l in pec_dict}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err


def main():
    print("\n" + "=" * 96)
    print(f"  task32j_bootstrap_estimator.py -- {N_BOOTSTRAP} full end-to-end bootstrap replicates, forte-1")
    print("=" * 96)

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
    print(f"  loaded Task 31C's real checkpoint: {n_mc} MC draws/label available")

    print(f"  running {N_BOOTSTRAP} replicates across separate processes...")
    with ProcessPoolExecutor(max_workers=2) as ex:
        futures = [ex.submit(_one_bootstrap_replicate, p, non_id_labels, diag, kept, P_S, by_name_label,
                              n_mc, seed=10_000 + i)
                   for i in range(N_BOOTSTRAP)]
        errs = []
        for i, fut in enumerate(futures):
            errs.append(fut.result())
            if (i + 1) % 20 == 0:
                print(f"    {i+1}/{N_BOOTSTRAP} replicates done")
    errs = np.array(errs)

    mean_e, median_e, std_e = float(errs.mean()), float(np.median(errs)), float(errs.std(ddof=1))
    ci_lo, ci_hi = float(np.percentile(errs, 2.5)), float(np.percentile(errs, 97.5))
    q = {qq: float(np.percentile(errs, qq)) for qq in [50, 90, 95, 99]}

    print(f"\n  -- BOOTSTRAP DISTRIBUTION (N={N_BOOTSTRAP}, no covariance matrix anywhere) --")
    print(f"    mean={mean_e:.4f}  median={median_e:.4f}  std={std_e:.4f}")
    print(f"    95% CI (percentile method) = [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"    Q50={q[50]:.3f}  Q90={q[90]:.3f}  Q95={q[95]:.3f}  Q99={q[99]:.3f}")

    print(f"\n  COMPARISON TO PRIOR POINT ESTIMATES on the SAME underlying data:")
    print(f"    Task 31C (N_MC=16 fixed, uniform-weighted, single realization): 4.495 kcal/mol")
    print(f"    Task 31D (covariance-aware, FAIL): 16.987 kcal/mol")
    print(f"    Task 31D (uniform-weighted, 32-seed bootstrap of shots only): 0.317 kcal/mol")
    print(f"    THIS bootstrap (shots + MC-draw-selection + manifold-refit, jointly, {N_BOOTSTRAP} reps): "
          f"{mean_e:.3f} +/- {std_e:.3f} (median {median_e:.3f})")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "n_bootstrap": N_BOOTSTRAP, "errs": errs.tolist(), "mean": mean_e, "median": median_e,
            "std": std_e, "ci95_lo": ci_lo, "ci95_hi": ci_hi, "quantiles": q,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
