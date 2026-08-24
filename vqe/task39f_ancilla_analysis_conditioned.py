#!/usr/bin/env python3
"""
task39f_ancilla_analysis_conditioned.py -- iteration 39, Task F. Re-runs
Task D's real 3-way comparison, but with case (c) now using Task 39E's
CONDITIONED correction (`analytic_A_and_B_conditioned`) instead of the
unconditional `analytic_A_and_B_5param` that produced Task D's negative
result. Cases (a) and (b) are unchanged (not postselected data, so the
unconditional correction is the correct one for them) -- only (c)
changes, isolating exactly the effect of fixing the conditioning
mismatch Task D's own honest read identified as the likely cause.

Run:
    PYTHONHASHSEED=0 python vqe/task39f_ancilla_analysis_conditioned.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import build_full
from task37c_extended_forward_model import analytic_A_and_B_5param
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from ionq_simulator_binding_curve import expectation_from_counts

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
OLD_RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                             "task28b_optimized_raw.json")
NEW_CKPT = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission.partial.json")

_CACHE_UNC = {}
_CACHE_COND = {}


def corrected_energy(raw_by_slot_label, kept, non_id_labels, fixed_solutions, diag, p, conditioned):
    raw_kept = {}
    for name in kept:
        cache = _CACHE_COND if conditioned else _CACHE_UNC
        if name not in cache:
            if conditioned:
                A, B, retA, retB = analytic_A_and_B_conditioned(fixed_solutions[name]["angles"], GATE_NAME,
                                                                   ZZ_ASSUMED, GPI_REAL_MEAN, 0.0, 0.0, 0.0,
                                                                   non_id_labels)
            else:
                A, B = analytic_A_and_B_5param(fixed_solutions[name]["angles"], GATE_NAME, ZZ_ASSUMED,
                                                 GPI_REAL_MEAN, 0.0, 0.0, 0.0, non_id_labels)
            cache[name] = (A, B)
        A, B = cache[name]
        raw_kept[name] = {}
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    full = build_full(raw_kept, diag, K, non_id_labels)
    alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def load_old_baseline(kept, backend_name):
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
    print("  task39f_ancilla_analysis_conditioned.py -- same 3-way comparison, CONDITIONED correction for (c)")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    for backend_name in ["aria-1", "forte-1"]:
        old_blended = load_old_baseline(kept, backend_name)
        new_nops_blended, _ = load_new_data(kept, backend_name, postselect=False)
        new_ps_blended, retained_frac = load_new_data(kept, backend_name, postselect=True)

        err_old = corrected_energy(old_blended, kept, non_id_labels, fixed_solutions, diag, p, conditioned=False)
        err_new_nops = corrected_energy(new_nops_blended, kept, non_id_labels, fixed_solutions, diag, p, conditioned=False)
        err_new_ps_unc = corrected_energy(new_ps_blended, kept, non_id_labels, fixed_solutions, diag, p, conditioned=False)
        err_new_ps_cond = corrected_energy(new_ps_blended, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)

        print(f"\n  === {backend_name} ===")
        print(f"    (a) OLD baseline (no ancilla, 100k shots):                    err = {err_old:+.4f} kcal/mol")
        print(f"    (b) NEW, ancilla, NOT postselected (20k):                     err = {err_new_nops:+.4f} kcal/mol")
        print(f"    (c-unconditional, Task D's original) postselected, WRONG correction: err = {err_new_ps_unc:+.4f} kcal/mol")
        print(f"    (c-CONDITIONED, this task's fix) postselected, correct ratio: err = {err_new_ps_cond:+.4f} kcal/mol")
        print(f"    retained fraction: {100*retained_frac:.1f}%")
        fixed_vs_wrong = err_new_ps_unc - err_new_ps_cond
        fixed_vs_old = err_old - err_new_ps_cond
        print(f"    conditioning fix changed |err| by {fixed_vs_wrong:+.3f} kcal/mol vs the unconditional (wrong) correction")
        print(f"    net vs OLD baseline: {fixed_vs_old:+.3f} kcal/mol "
              f"({'NET WIN over old baseline' if fixed_vs_old > 0 else 'still net LOSS vs old baseline'})")

    print(f"\n  Chemical accuracy bars for reference: 0.25 kcal/mol (internal target), 0.5 kcal/mol (looser hardware bar).")


if __name__ == "__main__":
    main()
