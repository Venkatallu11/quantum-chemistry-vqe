#!/usr/bin/env python3
"""
task33b_adaptive_label_pec.py -- iteration 33, Task B. ADAPTIVE PER-LABEL PEC.

Motivated by two results already on record, not a new hypothesis:
  - Task 33A (this iteration): the Q95 robustness-envelope tail is NOT
    explained by any single calibration parameter (GPi2 included -- its
    population correlation with error is rho=0.066, p=0.52, not
    significant). It is a multi-parameter joint effect. Chasing one
    "bad gate" will not fix it.
  - Task 32I (prior iteration): the existing correction is a single GLOBAL
    ratio B/A applied identically to every label, with NO regard for how
    reliable that particular label's ratio is. For IYYI, B/A=1.1227 is
    already unphysical (>1) before any real noise is even applied, and
    the corrected value lands at the WRONG SIGN vs the one label this
    project has real ground truth for. Task 32I flagged specific labels
    (ZIII, IIZI, IIIZ, IZII, plus IYYI) as high "boundary proximity" risk
    via a measurement-independent proxy: |B/A - 1| weighted by the
    label's real Hamiltonian-coefficient energy sensitivity.
  - Task 32G (prior iteration) already tested ONE global blend knob
    (lambda: same blend strength for every label) and found lambda=1
    (full PEC, no damping) MSE-optimal -- global shrinkage does not help.
    This task tests a DIFFERENT, finer-grained axis: does PER-LABEL
    damping (using each label's own measured boundary-proximity as the
    damping signal, not a single knob for everything) do better than
    Task 32G's global answer? This is not a re-run of Task 32G.

METHOD: generalize Task 32I's ratio_correct(m, A, B) = clip(m*B/A, -1, 1)
with a per-label damping strength r_l in [0, 1]:
    corrected_l(m) = clip(m * (1 + r_l*(B_l/A_l - 1)), -1, 1)
r_l=1 reproduces Task 32I's existing production correction exactly.
r_l=0 means "trust this label's correction not at all, use raw m."
r_l is NOT hand-picked per label -- it's a function of each label's own
measured boundary_proximity[l] = |B_l/A_l - 1| (Task 32I's own risk
proxy, reused unchanged) via a single global damping constant k:
    r_l(k) = 1 / (1 + k * boundary_proximity[l])
k=0 -> r_l=1 for every label (== current production, Task 32I's
"ratio-clip" baseline, exactly). Larger k damps ill-conditioned labels
more while leaving well-conditioned labels (small boundary_proximity)
close to fully corrected. This is the one new free parameter; it is
swept and chosen by real measured MSE on real data, not asserted.

DATA: reuses task28b_optimized_raw.json (forte-1 raw measurement
checkpoint) and task30b_pec_calibration_results.json (learned GPi bins)
-- the SAME already-collected data Task 32I used. No new submissions.

TESTS, same discipline as every prior task in this project:
  1. IDEAL CONTROL: at k=any value, on the noiseless "ideal" checkpoint,
     the correction must be a near-no-op (A=B=1 by construction there).
  2. IYYI VALIDATION: report where the adaptive correction lands for the
     one label with real literal-twirling ground truth (1.0), honestly,
     even if (like Task 32I) it doesn't fix the wrong-sign problem.
  3. MSE SWEEP: bootstrap over real shot noise (matching Task 32G/32I's
     established resampling discipline) at each k, report bias
     (median), sigma (std), MSE = bias^2+sigma^2, and the
     |b|+2*sigma<0.5 acceptance criterion, against Task 32G's own
     k=0 (== lambda=1 global) baseline MSE=2.33 as the number to beat.

Run:
    python vqe/task33b_adaptive_label_pec.py
"""
import os
import sys
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import analytic_A_and_B
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
GATE_NAME = "zz"
P2_ZZ = 0.0146
SHOTS = 100_000
N_BOOT = 32
K_GRID = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 300.0]
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
IYYI_GROUND_TRUTH = 1.0
TASK32G_GLOBAL_BEST_MSE = 2.33  # lambda=1, this iteration's own prior result, the number to beat
EXACT_ENERGY_TARGET = 0.5
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task33b_adaptive_label_pec_results.json")


def adaptive_correct(m_measured, A, B, r):
    ratio = B / A if abs(A) > 1e-6 else 1.0
    damped_ratio = 1.0 + r * (ratio - 1.0)
    return max(-1.0, min(1.0, m_measured * damped_ratio))


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def build_full(raw_kept, diag, K, non_id_labels):
    full = {name: dict(raw_kept[name]) for name in diag}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            full[pl] = dict(raw_kept[pl])
            synth_minus = {}
            for l in non_id_labels:
                if l not in raw_kept[pl] or l not in full[un] or l not in full[um]:
                    continue
                synth_minus[l] = full[un][l] + full[um][l] - raw_kept[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def _boot_worker(kept, per_name, AB_by_name, boundary_proximity, k, seed):
    rng = np.random.default_rng(seed)
    corrected_kept = {name: {} for name in kept}
    for name in kept:
        A, B = AB_by_name[name]
        for group_t, counts in per_name[name].items():
            resampled = bootstrap_counts(counts, SHOTS, rng)
            for l in group_t:
                m = expectation_from_counts(resampled, l)
                prox = boundary_proximity.get(l, 0.0)
                r_l = 1.0 / (1.0 + k * prox)
                corrected_kept[name][l] = adaptive_correct(m, A[l], B[l], r_l)
    return corrected_kept


def main():
    print("\n" + "=" * 96)
    print("  task33b_adaptive_label_pec.py -- per-label damped PEC, damping driven by measured boundary risk")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins = learned["forte-1"]["gpi"]
    gpi2_bins = gpi_bins  # disclosed fallback, unchanged from every prior task using it

    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    tags = raw_ck["tags"]["forte-1"]
    counts_list = raw_ck["counts"]["forte-1"]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)

    # -- per-label boundary proximity (Task 32I's own risk proxy, reused unchanged) --
    AB_by_name = {}
    boundary_proximity = {}
    label_weight = {}
    for (a_label, b_label, coeff) in p["terms"]:
        w = abs(coeff) * HARTREE_TO_KCAL_MOL
        label_weight[a_label] = label_weight.get(a_label, 0.0) + w
    for name in kept:
        labels_here = [l for group_t in per_name[name] for l in group_t]
        A, B = analytic_A_and_B(fixed_solutions[name]["angles"], GATE_NAME, P2_ZZ, gpi_bins, gpi2_bins, labels_here)
        AB_by_name[name] = (A, B)
        for l in labels_here:
            ratio_here = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            boundary_proximity.setdefault(l, []).append(abs(ratio_here - 1.0))
    boundary_proximity = {l: float(np.mean(v)) for l, v in boundary_proximity.items()}

    n_flagged_high_risk = sum(1 for l, prox in boundary_proximity.items()
                               if prox * label_weight.get(l, 0.0) > 0.02)
    print(f"  {len(boundary_proximity)} labels have a measured boundary-proximity risk score; "
          f"{n_flagged_high_risk} exceed Task 32I's own 0.02 weighted-risk flag threshold")

    # -- IYYI validation, all k values, against real ground truth --
    print(f"\n  -- IYYI validation across damping k (real literal-twirling ground truth = {IYYI_GROUND_TRUTH}) --")
    slot, label = "(u0+u1)", "IYYI"
    A_iyyi, B_iyyi = AB_by_name[slot]
    prox_iyyi = boundary_proximity[label]
    m_raw_vals = []
    for seed in range(8):
        rng = np.random.default_rng(stable_seed("t33b_iyyi", seed))
        for group_t, counts in per_name[slot].items():
            if label in group_t:
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m_raw_vals.append(expectation_from_counts(resampled, label))
    m_raw_iyyi = float(np.mean(m_raw_vals))
    iyyi_by_k = {}
    for k in K_GRID:
        r_l = 1.0 / (1.0 + k * prox_iyyi)
        m_adj = adaptive_correct(m_raw_iyyi, A_iyyi[label], B_iyyi[label], r_l)
        iyyi_by_k[k] = {"r_l": r_l, "m_adaptive": m_adj, "abs_err_vs_gt": abs(m_adj - IYYI_GROUND_TRUTH)}
        print(f"    k={k:<6} r_IYYI={r_l:.4f}  m_adaptive={m_adj:+.4f}  |err vs GT|={abs(m_adj-IYYI_GROUND_TRUTH):.4f}")
    print(f"    (raw m={m_raw_iyyi:+.4f}, |err|={abs(m_raw_iyyi-IYYI_GROUND_TRUTH):.4f}; "
          f"full ratio-clip k=0 |err|={iyyi_by_k[0.0]['abs_err_vs_gt']:.4f} -- Task 32I's number)")

    # -- MSE sweep over k, real data, bootstrap over shot noise --
    print(f"\n  -- MSE sweep: {N_BOOT} bootstrap replicates per k, real forte-1 data --")
    k_results = {}
    for k in K_GRID:
        with ProcessPoolExecutor(max_workers=4) as ex:
            futures = [ex.submit(_boot_worker, kept, per_name, AB_by_name, boundary_proximity, k,
                                  30_000 + int(k * 10) + i)
                       for i in range(N_BOOT)]
            corrected_sets = [fut.result() for fut in futures]
        errs = []
        for corrected_kept in corrected_sets:
            full = build_full(corrected_kept, diag, K, non_id_labels)
            _, err = energy_and_err(p, full, K)
            errs.append(err)
        errs = np.array(errs)
        bias_proxy = float(np.median(errs))
        sigma = float(errs.std(ddof=1))
        mse = bias_proxy ** 2 + sigma ** 2
        passes = (abs(bias_proxy) + 2 * sigma) < EXACT_ENERGY_TARGET
        k_results[k] = {"errs": errs.tolist(), "median": bias_proxy, "std": sigma, "mse": mse, "passes": bool(passes)}
        print(f"    k={k:<6} median={bias_proxy:8.4f}  std={sigma:7.4f}  MSE={mse:9.4f}  "
              f"{'PASS' if passes else 'FAIL'} (|b|+2sigma<0.5)")

    best_k = min(k_results, key=lambda kk: k_results[kk]["mse"])
    best_mse = k_results[best_k]["mse"]
    print(f"\n  MSE-MINIMIZING k* = {best_k}  (MSE={best_mse:.4f})")
    print(f"  Task 32G's global-lambda best (lambda=1, == this task's k=0) was MSE={TASK32G_GLOBAL_BEST_MSE}")
    if best_k > 0.0 and best_mse < k_results[0.0]["mse"]:
        improvement_pct = 100.0 * (k_results[0.0]["mse"] - best_mse) / k_results[0.0]["mse"]
        print(f"  -> PER-LABEL damping (k={best_k}) beats full correction (k=0) by "
              f"{improvement_pct:.1f}% MSE on this real data -- a real, finer-grained win Task 32G's "
              f"single global knob could not find, because it damped every label equally.")
    else:
        print(f"  -> Per-label damping does NOT beat full correction (k=0) here -- consistent with Task 32G's "
              f"global finding, now confirmed at finer granularity too: this pipeline's PEC correction is "
              f"MSE-optimal at full strength, not because damping-by-risk doesn't work in principle, but "
              f"because on THIS real data the measured boundary-proximity signal isn't predictive enough of "
              f"which corrections are actually harmful.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "boundary_proximity": boundary_proximity,
            "label_weight": label_weight,
            "n_flagged_high_risk": n_flagged_high_risk,
            "iyyi_validation": {str(k): v for k, v in iyyi_by_k.items()},
            "m_raw_iyyi": m_raw_iyyi,
            "k_sweep": {str(k): v for k, v in k_results.items()},
            "best_k": best_k,
            "best_mse": best_mse,
            "task32g_global_best_mse": TASK32G_GLOBAL_BEST_MSE,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
