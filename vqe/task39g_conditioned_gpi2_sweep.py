#!/usr/bin/env python3
"""
task39g_conditioned_gpi2_sweep.py -- iteration 39, Task G. Task F fixed
the conditioning mismatch but left p_gpi2_assumed at 0 (uncorrected)
throughout. Task 38D2 found, on the OLD (no ancilla, no postselection)
real data, that ANY nonzero p_gpi2 correction made things monotonically
WORSE. Does that still hold once the data has been postselected on the
ancilla and corrected with Task 39E's conditioned formula, or does
postselection change which correction is actually best? Direct,
real-data grid sweep, mirroring Task 38D2's own methodology exactly
(same grid, same real-data-only approach, no synthetic noise guessing).

Run:
    PYTHONHASHSEED=0 python vqe/task39g_conditioned_gpi2_sweep.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import build_full
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from ionq_simulator_binding_curve import expectation_from_counts

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
GPI2_GRID = [0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001]  # bracketing aria-1's true
# minimum (still improving past 0.00065 in the previous finer sweep) while keeping forte-1's own
# minimum (~0.0004) in view, to see whether the two backends share any reasonable joint sweet spot
NEW_CKPT = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission.partial.json")

_CACHE = {}


def corrected_energy(raw_by_slot_label, p_gpi2, kept, non_id_labels, fixed_solutions, diag, p):
    raw_kept = {}
    for name in kept:
        key = (name, round(p_gpi2, 6))
        if key not in _CACHE:
            A, B, retA, retB = analytic_A_and_B_conditioned(fixed_solutions[name]["angles"], GATE_NAME,
                                                               ZZ_ASSUMED, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0,
                                                               non_id_labels)
            _CACHE[key] = (A, B)
        A, B = _CACHE[key]
        raw_kept[name] = {}
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    full = build_full(raw_kept, diag, K, non_id_labels)
    alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def load_postselected(kept, backend_name):
    with open(NEW_CKPT) as f:
        state = json.load(f)
    blended = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for group, counts in zip(entry["groups"], entry["counts"]):
            filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
            for l in group:
                blended[name][l] = expectation_from_counts(filtered, l)
    return blended


def main():
    print("\n" + "=" * 96)
    print("  task39g_conditioned_gpi2_sweep.py -- does a nonzero GPi2 correction help NOW, on postselected data?")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    print(f"\n  {'p_gpi2_assumed':>16}{'err_aria-1':>16}{'err_forte-1':>16}{'worst-case':>14}")
    aria_ps = load_postselected(kept, "aria-1")
    forte_ps = load_postselected(kept, "forte-1")
    results = []
    for p_gpi2 in GPI2_GRID:
        e_aria = corrected_energy(aria_ps, p_gpi2, kept, non_id_labels, fixed_solutions, diag, p)
        e_forte = corrected_energy(forte_ps, p_gpi2, kept, non_id_labels, fixed_solutions, diag, p)
        worst = max(abs(e_aria), abs(e_forte))
        results.append((p_gpi2, e_aria, e_forte, worst))
        print(f"  {p_gpi2:>16.4f}{e_aria:>16.3f}{e_forte:>16.3f}{worst:>14.3f}")

    baseline = results[0]
    best = min(results, key=lambda r: r[3])
    print(f"\n  BASELINE (p_gpi2_assumed=0, this task's own conditioned-postselected result): worst-case={baseline[3]:.3f}")
    print(f"  BEST candidate: p_gpi2_assumed={best[0]}  worst-case={best[3]:.3f}")
    if best[3] < baseline[3] and best[0] != 0.0:
        pct = 100 * (baseline[3] - best[3]) / baseline[3]
        print(f"\n  A nonzero GPi2 correction DOES help on the postselected data ({pct:.1f}% worst-case reduction) -- "
              f"different from Task 38D2's finding on the un-postselected data. Postselection changes which "
              f"correction is optimal, a real and interesting result if it holds.")
    else:
        print(f"\n  Consistent with Task 38D2: no nonzero GPi2 correction improves the worst-case error here either -- "
              f"the same finding holds regardless of postselection.")


if __name__ == "__main__":
    main()
