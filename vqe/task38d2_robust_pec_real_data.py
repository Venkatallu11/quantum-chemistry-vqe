#!/usr/bin/env python3
"""
task38d2_robust_pec_real_data.py -- iteration 38, Task D corrected retry.
`task38d_robust_pec.py`'s synthetic-simulation approach hit a real,
disclosed scale-mismatch bug (its baseline error was 5-7x larger than
Task 38C's real-data reference at the same nominal point) -- rather than
trying to fix the guessed p_gpi2_true distribution (which would require
knowing the very thing this project has repeatedly failed to calibrate
from scratch), this retry sidesteps the problem entirely: use REAL raw
data directly, on TWO genuinely different real noise realizations
(`task28b_optimized_raw.json`'s aria-1 AND forte-1 backends -- reusing
Task 38C's own real-data pipeline unchanged, just no longer hardcoded to
forte-1 only) instead of one synthetic "true theta" ensemble.

QUESTION: does a candidate GPi2 correction (a nonzero assumed p_gpi2
for the PEC-inverse step, applied via Task 37C's own verified
`analytic_A_and_B_5param`) generalize across aria-1 and forte-1's real,
independently-realized noise, or does it only help one backend while
hurting the other -- which would be a DIRECT, real demonstration of the
exact fragility Task 38's whole redirection is about, not a synthetic
guess at it. No new real submissions (avoids the quota block hit
earlier); no guessed noise-model distribution (avoids the scale-mismatch
bug); reuses Task 38C's own `energy_of_theta` function unchanged.

Run:
    PYTHONHASHSEED=0 python vqe/task38d2_robust_pec_real_data.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts, stable_seed
from task38c_sensitivity_analysis import energy_of_theta

K = 6
GPI2_GRID = [0.0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0008, 0.0015, 0.003]
ZZ_ASSUMED = 0.014593
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
SHOTS = 100_000


def load_real_blended_backend(kept, backend_name, n_seeds=8):
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    tags = raw_ck["tags"][backend_name]
    counts_list = raw_ck["counts"][backend_name]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)
    blended = {name: {} for name in kept}
    for name in kept:
        for group_t, counts in per_name[name].items():
            for seed in range(n_seeds):
                rng = np.random.default_rng(stable_seed("t38d2", backend_name, name, seed))
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    m = expectation_from_counts(resampled, l)
                    blended[name].setdefault(l, []).append(m)
    for name in kept:
        for l in blended[name]:
            blended[name][l] = float(np.mean(blended[name][l]))
    return blended


def main():
    print("\n" + "=" * 96)
    print("  task38d2_robust_pec_real_data.py -- GPi2 correction robustness on TWO real backends (aria-1, forte-1)")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    real_aria = load_real_blended_backend(kept, "aria-1")
    real_forte = load_real_blended_backend(kept, "forte-1")
    print(f"\n  loaded real raw data: aria-1 ({sum(len(real_aria[n]) for n in kept)} residuals), "
          f"forte-1 ({sum(len(real_forte[n]) for n in kept)} residuals)")

    print(f"\n  {'p_gpi2_assumed':>16}{'err_aria-1 (kcal/mol)':>24}{'err_forte-1 (kcal/mol)':>24}"
          f"{'worst-case':>14}{'|spread|':>12}")
    results = []
    for p_gpi2 in GPI2_GRID:
        theta = np.array([ZZ_ASSUMED, p_gpi2, 0.0, 0.0, 0.0])
        e_aria = energy_of_theta(theta, kept, non_id_labels, fixed_solutions, diag, real_aria, p)
        e_forte = energy_of_theta(theta, kept, non_id_labels, fixed_solutions, diag, real_forte, p)
        worst = max(abs(e_aria), abs(e_forte))
        spread = abs(e_aria - e_forte)
        results.append((p_gpi2, e_aria, e_forte, worst, spread))
        print(f"  {p_gpi2:>16.4f}{e_aria:>24.3f}{e_forte:>24.3f}{worst:>14.3f}{spread:>12.3f}")

    baseline = results[0]
    best_worst_case = min(results, key=lambda r: r[3])
    best_spread = min(results, key=lambda r: r[4])
    print(f"\n  BASELINE (p_gpi2_assumed=0): err_aria={baseline[1]:.3f}  err_forte={baseline[2]:.3f}  "
          f"worst-case={baseline[3]:.3f}  spread={baseline[4]:.3f}")
    print(f"  BEST WORST-CASE candidate:   p_gpi2_assumed={best_worst_case[0]}  worst-case={best_worst_case[3]:.3f}  "
          f"(vs baseline {baseline[3]:.3f})")
    print(f"  BEST SPREAD (most consistent across backends): p_gpi2_assumed={best_spread[0]}  "
          f"spread={best_spread[4]:.3f}  (vs baseline {baseline[4]:.3f})")

    print(f"\n  -- HONEST READ --")
    if best_worst_case[3] < baseline[3] and best_worst_case[0] != 0.0:
        pct = 100 * (baseline[3] - best_worst_case[3]) / baseline[3]
        print(f"    A nonzero GPi2 correction (p_gpi2_assumed={best_worst_case[0]}) improves the WORST-CASE error "
              f"across both real backends by {pct:.1f}% vs applying no correction at all -- real, on real data, "
              f"no synthetic noise-model guess involved. This directly supports introducing SOME GPi2 correction.")
    else:
        print(f"    No nonzero GPi2 correction improves the worst-case error across BOTH real backends vs the "
              f"zero-correction baseline -- a real, honest negative result: any single assumed p_gpi2 that helps "
              f"one backend's real data does not reliably help the other's, which is itself a direct, real "
              f"demonstration of exactly the calibration-fragility problem this whole investigation is about.")
    if best_spread[4] < baseline[4]:
        print(f"    Separately: p_gpi2_assumed={best_spread[0]} gives the most CONSISTENT (least backend-dependent) "
              f"correction across the two real backends, even if not necessarily the lowest worst-case error.")


if __name__ == "__main__":
    main()
