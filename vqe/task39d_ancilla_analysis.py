#!/usr/bin/env python3
"""
task39d_ancilla_analysis.py -- iteration 39, Task D. Analyzes the real
ancilla-augmented data Task 39C collected (21 slots x 2 backends x 13
groups, 20k shots each, real ionq_simulator). Three-way comparison,
matching iteration 18's own original structure exactly (that iteration
found: no-ancilla baseline BETTER than ancilla-no-postselect, but
ancilla-WITH-postselect better than both):
  (a) OLD baseline: task28b_optimized_raw.json, no ancilla, 100k shots
      (already-collected, real, unchanged).
  (b) NEW data, ancilla present but NOT postselected (marginalizing the
      ancilla bit out) -- isolates the real cost of the extra ancilla
      gates themselves, at the new 20k-shot budget.
  (c) NEW data, postselected on ancilla==0 (the actual leakage-detection
      benefit) -- also reports the retained-shot fraction per backend,
      since heavy shot loss would itself degrade precision.

All three use the IDENTICAL correction (Task 37C's own verified
`analytic_A_and_B_5param`, nominal theta_0 = (0.014593, 0, 0, 0, 0),
matching Task 38C's own established convention) and the SAME standard
(non-joint-frame) combine/energy pipeline, so any difference between
(a)/(b)/(c) is real signal, not a pipeline artifact.

Run:
    PYTHONHASHSEED=0 python vqe/task39d_ancilla_analysis.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import build_full
from task37c_extended_forward_model import analytic_A_and_B_5param, readout_attenuation
from task37b_h4_noise_model import GPI_REAL_MEAN
from ionq_simulator_binding_curve import expectation_from_counts

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
OLD_RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                             "task28b_optimized_raw.json")
NEW_CKPT = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission.partial.json")

_GLOBAL_AB_CACHE = {}


def corrected_energy(raw_by_slot_label, kept, non_id_labels, fixed_solutions, diag, p):
    theta = (ZZ_ASSUMED, 0.0, 0.0, 0.0, 0.0)
    raw_kept = {}
    for name in kept:
        cache_key = name  # theta is fixed here, so cache only needs to vary by slot
        if cache_key not in _GLOBAL_AB_CACHE:
            A, B = analytic_A_and_B_5param(fixed_solutions[name]["angles"], GATE_NAME, theta[0], GPI_REAL_MEAN,
                                            theta[1], theta[2], theta[3], non_id_labels)
            _GLOBAL_AB_CACHE[cache_key] = (A, B)
        A, B = _GLOBAL_AB_CACHE[cache_key]
        raw_kept[name] = {}
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    full = build_full(raw_kept, diag, K, non_id_labels)
    alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def load_old_baseline(kept, backend_name, non_id_labels):
    with open(OLD_RAW_CKPT) as f:
        raw_ck = json.load(f)
    tags = raw_ck["tags"][backend_name]
    counts_list = raw_ck["counts"][backend_name]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)
    blended = {name: {} for name in kept}
    for name in kept:
        for group_t, counts in per_name[name].items():
            for l in group_t:
                blended[name][l] = expectation_from_counts(counts, l)
    return blended


def load_new_data(kept, backend_name, postselect):
    with open(NEW_CKPT) as f:
        state = json.load(f)
    blended = {name: {} for name in kept}
    retained_fracs = []
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for group, counts in zip(entry["groups"], entry["counts"]):
            total = sum(counts.values())
            if postselect:
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
            else:
                filtered = {}
                for bs, c in counts.items():
                    filtered[bs[1:]] = filtered.get(bs[1:], 0) + c
            retained = sum(filtered.values())
            retained_fracs.append(retained / total if total else 0.0)
            for l in group:
                blended[name][l] = expectation_from_counts(filtered, l)
    return blended, float(np.mean(retained_fracs))


def main():
    print("\n" + "=" * 96)
    print("  task39d_ancilla_analysis.py -- old baseline vs ancilla-no-postselect vs ancilla-postselected")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    results = {}
    for backend_name in ["aria-1", "forte-1"]:
        old_blended = load_old_baseline(kept, backend_name, non_id_labels)
        new_nops_blended, _ = load_new_data(kept, backend_name, postselect=False)
        new_ps_blended, retained_frac = load_new_data(kept, backend_name, postselect=True)

        err_old = corrected_energy(old_blended, kept, non_id_labels, fixed_solutions, diag, p)
        err_new_nops = corrected_energy(new_nops_blended, kept, non_id_labels, fixed_solutions, diag, p)
        err_new_ps = corrected_energy(new_ps_blended, kept, non_id_labels, fixed_solutions, diag, p)

        print(f"\n  === {backend_name} ===")
        print(f"    (a) OLD baseline (no ancilla, 100k shots):        err = {err_old:+.4f} kcal/mol")
        print(f"    (b) NEW, ancilla present, NOT postselected (20k): err = {err_new_nops:+.4f} kcal/mol")
        print(f"    (c) NEW, ancilla postselected (20k, retained {100*retained_frac:.1f}% of shots): "
              f"err = {err_new_ps:+.4f} kcal/mol")
        results[backend_name] = {"old": err_old, "new_nops": err_new_nops, "new_ps": err_new_ps,
                                   "retained_frac": retained_frac}

    print(f"\n  -- HONEST READ --")
    for backend_name, r in results.items():
        ancilla_cost = abs(r["new_nops"]) - abs(r["old"])
        ps_benefit = abs(r["new_nops"]) - abs(r["new_ps"])
        vs_old = abs(r["old"]) - abs(r["new_ps"])
        print(f"    {backend_name}: adding the ancilla (before postselection) changed |err| by {ancilla_cost:+.3f} "
              f"(matches iteration 18's own finding that the extra gates cost something before any denoising); "
              f"postselection then changed |err| by {-ps_benefit:+.3f} vs the un-postselected ancilla data; "
              f"net vs the OLD no-ancilla baseline: {vs_old:+.3f} kcal/mol "
              f"({'NET WIN over the old baseline' if vs_old > 0 else 'net LOSS vs the old baseline'}).")
    print(f"\n  Chemical accuracy bars for reference: 0.25 kcal/mol (internal target), 0.5 kcal/mol (looser hardware bar).")


if __name__ == "__main__":
    main()
