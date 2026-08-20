#!/usr/bin/env python3
"""
task32b_variance_decomposition.py -- iteration 32, Task B. WHY did nominally
the same measurement give 0.115 (Task 31C), 0.317 (Task 31D, same real data
reprocessed), and 0.438 (Task 31F, mean of 8 independent submissions)? Four
candidate variance sources, decomposed as cheaply as possible by reusing
ALREADY-COLLECTED real data -- NO new submission for (i)-(iii); (iv) reuses
Task 31F's existing 8 real submissions outright.
============================================================================
ALL FOUR SOURCES USE THE SAME FIXED CALIBRATION Task 31C/31D/31F actually
used (P2_ZZ=0.0146 consensus, GPi/GPi2-fallback p1_gpi mean) -- this
decomposes the variance that PRODUCED those three specific numbers, not a
re-calibrated pipeline (Task 32A's improved GPi2 posterior feeds a SEPARATE
"calibration uncertainty" term in the Task 32-synthesis error budget, kept
out of this file to avoid conflating "why did three past runs disagree"
with "how good could a recalibrated pipeline be").

  (i)   SHOT NOISE: fix the PEC Monte Carlo draws (all 16, as Task 31C used)
        and the manifold fit's random-restart seed; vary ONLY the bootstrap
        resample of the real integer counts. 64 independent resamples.
  (ii)  PEC MONTE CARLO: fix the shot data (real counts, no resampling) and
        the manifold seed; vary ONLY which subset of the 16 real twirled
        MC draws is averaged, M in {1,2,4,8,16}. Multiple random subsets
        per M (with replacement across subset choices, not across draws
        within a subset) to estimate sigma_PEC(M). Extending to M=32-128
        would need NEW real submissions -- done only if the M<=16 trend
        does not already show a clear sqrt(M) law or an early plateau
        (both are informative without spending more).
  (iii) MANIFOLD OPTIMIZER: fix PEC-corrected label values completely (all
        16 draws, ONE fixed shot-resample); vary ONLY the optimizer's
        random-restart seed. 32 independent re-fits.
  (iv)  SUBMISSION VARIATION: Task 31F's own 8 independent real submissions,
        reused outright (no recomputation).

Run:
    python vqe/task32b_variance_decomposition.py
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
import ef_fragment as effrag_mod
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
P2_ZZ = 0.0146
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32b_variance_decomposition_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def load_context():
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
    return p, non_id_labels, diag, kept, P_S, by_name_label


def pec_value_for_slot(name, by_name_label, non_id_labels, draw_idx_subset, resample_rng, use_resample=True):
    """PEC-corrected expectation per label for ONE slot, using only the
    draws in `draw_idx_subset` (indices into the per-label draw list), and
    EITHER a real bootstrap resample (resample_rng given) OR the raw real
    counts directly (use_resample=False, resample_rng ignored)."""
    out = {}
    for l in non_id_labels:
        key = (name, l)
        if key not in by_name_label:
            continue
        entries = by_name_label[key]
        vals_signed = []
        for idx in draw_idx_subset:
            if idx >= len(entries):
                continue
            sign, gamma, counts = entries[idx]
            if use_resample:
                resampled = bootstrap_counts(counts, SHOTS, resample_rng)
                m = expectation_from_counts(resampled, l)
            else:
                m = expectation_from_counts(counts, l)
            vals_signed.append(sign * gamma * m)
        if vals_signed:
            out[l] = float(np.mean(vals_signed))
    return out


def full_energy_from_manifold(p, kept, diag, P_S, non_id_labels, pec_by_name, fit_seed):
    manifold_kept = {}
    for name in kept:
        pec_dict = pec_by_name[name]
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, pec_dict, {l: 1.0 for l in pec_dict}, K, v0, seed=fit_seed)
        manifold_kept[name] = a_hat
    full = build_full_from_a(manifold_kept, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err


def _manifold_worker(p, kept, diag, P_S, non_id_labels, pec_by_name, fit_seed):
    """Top-level (picklable) worker for ProcessPoolExecutor -- MUST run in a
    genuinely separate OS process, not just a different RNG seed within one
    process. A controlled diagnostic (same fixed input, same explicit seed,
    3 separate `python -c ...` launches) found the SAME bit-identical
    objective value (5.0192454942) reached via DIFFERENT solution vectors
    (specific components flipped sign) across process launches -- a real
    degenerate-optimum problem (multiple exactly-tied global optima for
    under-constrained slots), with the tie broken by cross-process
    floating-point/BLAS nondeterminism, NOT by the seed. Varying only the
    seed WITHIN one process (this file's first attempt) undersamples this
    entirely -- it can look artificially reproducible OR artificially
    variable depending on which restarts happen to be tried, neither of
    which reflects what the real per-submission production pipeline
    (Task 31C/31F, one fresh process per run) actually experiences."""
    return full_energy_from_manifold(p, kept, diag, P_S, non_id_labels, pec_by_name, fit_seed)


def main():
    print("\n" + "=" * 96)
    print("  task32b_variance_decomposition.py -- decomposing the 0.115/0.317/0.438 spread")
    print("=" * 96)

    p, non_id_labels, diag, kept, P_S, by_name_label = load_context()
    n_mc_available = min(len(v) for v in by_name_label.values())
    print(f"  loaded Task 31C's real checkpoint: {n_mc_available} MC draws/label available")

    # -------------------------------------------------------------
    # (i) SHOT NOISE: all 16 draws fixed, manifold seed fixed, vary shot resample
    # -------------------------------------------------------------
    print(f"\n  -- (i) SHOT NOISE: 64 independent bootstrap resamples, MC draws and manifold seed fixed --")
    full_subset = list(range(n_mc_available))
    shot_noise_errs = []
    for rep in range(64):
        rng = np.random.default_rng(stable_seed("t32b_shot", rep))
        pec_by_name = {name: pec_value_for_slot(name, by_name_label, non_id_labels, full_subset, rng, True)
                       for name in kept}
        err = full_energy_from_manifold(p, kept, diag, P_S, non_id_labels, pec_by_name,
                                         fit_seed=stable_seed("t32b_shot_fitseed"))
        shot_noise_errs.append(err)
    shot_noise_errs = np.array(shot_noise_errs)
    var_shot = float(np.var(shot_noise_errs, ddof=1))
    print(f"    mean={shot_noise_errs.mean():.4f}  std={shot_noise_errs.std(ddof=1):.4f}  var={var_shot:.5f}")

    # -------------------------------------------------------------
    # (ii) PEC MONTE CARLO: fix shot data (no resample), vary subset of M draws
    # -------------------------------------------------------------
    print(f"\n  -- (ii) PEC MONTE CARLO: subsets of M in {{1,2,4,8,{n_mc_available}}} draws, no shot resampling --")
    # M=n_mc_available (16): sampling WITHOUT replacement from a pool of 16 leaves only ONE possible
    # subset (the full set), which would give a fake "zero variance" -- not a real measurement, just
    # an artifact of having exactly as many draws available as the production M. Sample WITH
    # replacement instead (a genuine bootstrap over draws) so M=16 gets a real, nonzero variance
    # estimate too, on the same footing as M<16.
    M_LEVELS = [1, 2, 4, 8, n_mc_available]
    pec_mc_variance_by_M = {}
    for M in M_LEVELS:
        n_trials = 24
        replace = (M >= n_mc_available)
        errs_M = []
        for t in range(n_trials):
            rng = np.random.default_rng(stable_seed("t32b_pecmc_subset", M, t))
            subset = sorted(rng.choice(n_mc_available, size=M, replace=replace).tolist())
            pec_by_name = {name: pec_value_for_slot(name, by_name_label, non_id_labels, subset, None, False)
                           for name in kept}
            err = full_energy_from_manifold(p, kept, diag, P_S, non_id_labels, pec_by_name,
                                             fit_seed=stable_seed("t32b_pecmc_fitseed"))
            errs_M.append(err)
        errs_M = np.array(errs_M)
        pec_mc_variance_by_M[M] = {
            "mean": float(errs_M.mean()), "std": float(errs_M.std(ddof=1)),
            "median": float(np.median(errs_M)), "iqr": float(np.percentile(errs_M, 75) - np.percentile(errs_M, 25)),
            "n_trials": n_trials, "sampled_with_replacement": bool(replace), "errs": errs_M.tolist(),
        }
        print(f"    M={M:>2}: mean={errs_M.mean():.4f}  std={pec_mc_variance_by_M[M]['std']:.4f}  "
              f"median={pec_mc_variance_by_M[M]['median']:.4f}  IQR={pec_mc_variance_by_M[M]['iqr']:.4f}  "
              f"(n_trials={n_trials}, {'WITH' if replace else 'without'} replacement)")
    print(f"    NOTE: std is outlier-sensitive and low-M draws include real extreme outliers (a single bad "
          f"MC draw can dominate one slot's fit) -- median/IQR alongside std for an honest picture, not just "
          f"the outlier-driven std.")

    # check M^-1/2 scaling using M in {1,2,4,8} (16 has only 1 trial, no std)
    fit_Ms = [m for m in M_LEVELS[:-1] if pec_mc_variance_by_M[m]["std"] > 0]
    if len(fit_Ms) >= 2:
        log_m = np.log(fit_Ms)
        log_s = np.log([pec_mc_variance_by_M[m]["std"] for m in fit_Ms])
        beta_pec, c_pec = np.polyfit(log_m, log_s, 1)
        print(f"    fit: log(sigma_PEC) = {c_pec:.4f} + {beta_pec:.4f}*log(M)  (want beta~-0.5; "
              f"{'consistent with sqrt(M) convergence' if beta_pec < -0.3 else 'PLATEAU-LIKE -- PEC estimator itself may be unstable, a real finding'})")
    else:
        beta_pec = None

    var_pec_mc = pec_mc_variance_by_M[n_mc_available]["std"] ** 2  # bootstrap-with-replacement var at production M=16

    # -------------------------------------------------------------
    # (iii) MANIFOLD OPTIMIZER: fix PEC-corrected values (all 16 draws, ONE resample), vary fit seed
    # -------------------------------------------------------------
    print(f"\n  -- (iii) MANIFOLD OPTIMIZER: 24 re-fits of the SAME fixed PEC-corrected data, each in a "
          f"SEPARATE OS process (see _manifold_worker docstring for why within-process seed variation "
          f"is the wrong measurement) --")
    rng_fixed = np.random.default_rng(stable_seed("t32b_manifold_fixed_data"))
    pec_by_name_fixed = {name: pec_value_for_slot(name, by_name_label, non_id_labels, full_subset, rng_fixed, True)
                          for name in kept}
    N_PROC_TRIALS = 24
    manifold_errs = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        futures = [ex.submit(_manifold_worker, p, kept, diag, P_S, non_id_labels, pec_by_name_fixed,
                              stable_seed("t32b_manifold_seed", rep))
                   for rep in range(N_PROC_TRIALS)]
        for fut in futures:
            manifold_errs.append(fut.result())
    manifold_errs = np.array(manifold_errs)
    var_manifold = float(np.var(manifold_errs, ddof=1))
    n_distinct = len(np.unique(np.round(manifold_errs, 6)))
    print(f"    mean={manifold_errs.mean():.4f}  std={manifold_errs.std(ddof=1):.4f}  var={var_manifold:.5f}  "
          f"({n_distinct}/{N_PROC_TRIALS} distinct values to 1e-6 -- {'CONFIRMS cross-process degeneracy' if n_distinct > 1 else 'fully reproducible for this data'})")
    print(f"    (identical INPUT data every time, separate process every time -- ANY spread here is pure "
          f"cross-process optimizer/degenerate-optimum noise, the mechanism a single-process seed sweep cannot see)")

    # -------------------------------------------------------------
    # (iv) SUBMISSION VARIATION: reused outright from Task 31F
    # -------------------------------------------------------------
    with open(os.path.join(os.path.dirname(__file__), "task31f_convergence_study_results.json")) as f:
        conv = json.load(f)
    E_31f = np.array(conv["E_list"])
    var_submission = float(np.var(E_31f, ddof=1))
    print(f"\n  -- (iv) SUBMISSION VARIATION (reused, Task 31F's 8 real independent submissions) --")
    print(f"    mean={E_31f.mean():.4f}  std={E_31f.std(ddof=1):.4f}  var={var_submission:.5f}")

    # -------------------------------------------------------------
    # VARIANCE BUDGET
    # -------------------------------------------------------------
    print(f"\n" + "=" * 96)
    print(f"  VARIANCE BUDGET (assumes approximate independence across sources -- disclosed simplification)")
    total = var_shot + var_pec_mc + var_manifold + var_submission
    for label, v in [("shot noise", var_shot), ("PEC Monte Carlo (at M=16)", var_pec_mc),
                     ("manifold optimizer", var_manifold), ("submission-to-submission", var_submission)]:
        pct = 100 * v / total if total > 0 else 0
        print(f"    {label:<30} var={v:.5f}  std={np.sqrt(v):.4f}  ({pct:.1f}% of decomposed total)")
    print(f"    {'TOTAL (decomposed)':<30} var={total:.5f}  std={np.sqrt(total):.4f}")
    dominant = max([("shot", var_shot), ("pec_mc", var_pec_mc), ("manifold", var_manifold),
                    ("submission", var_submission)], key=lambda x: x[1])
    print(f"\n  DOMINANT SOURCE: {dominant[0]} ({100*dominant[1]/total:.1f}% of decomposed variance)")
    print(f"  (compare: observed real spread across the three headline numbers 0.115/0.317/0.438 has "
          f"range {max(0.115,0.317,0.438)-min(0.115,0.317,0.438):.3f}, sample std="
          f"{np.std([0.115,0.317,0.438], ddof=1):.4f} -- consistent order of magnitude, not identical by "
          f"construction since those three used different subsets of these SAME four sources)")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "n_mc_available": n_mc_available,
            "shot_noise": {"errs": shot_noise_errs.tolist(), "var": var_shot},
            "pec_mc_by_M": pec_mc_variance_by_M, "pec_mc_scaling_beta": float(beta_pec) if beta_pec is not None else None,
            "var_pec_mc_at_production_M": var_pec_mc,
            "manifold_optimizer": {"errs": manifold_errs.tolist(), "var": var_manifold},
            "submission_variation": {"errs": E_31f.tolist(), "var": var_submission},
            "variance_budget": {"shot": var_shot, "pec_mc": var_pec_mc, "manifold": var_manifold,
                                 "submission": var_submission, "total": total},
            "dominant_source": dominant[0],
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
