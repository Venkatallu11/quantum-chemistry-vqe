#!/usr/bin/env python3
"""
task69b_hw_plan_rehearsal_uncertainty.py -- iteration 69b. Same 56-circuit
plan rehearsal as task69, but repeated across 6 independent bootstrap
draws per backend to get a real uncertainty estimate (mean +/- std, not
a single point value) for all three locked analyses -- directly answers
IonQ's explicit request ("...with uncertainty estimates") before any
real hardware money is spent. Still zero new cost: same already-collected
real task59 free-simulator data, just resampled with different seeds.

Run:
    PYTHONHASHSEED=0 python vqe/task69b_hw_plan_rehearsal_uncertainty.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from general_commuting_measurements import build_general_commuting_measurement_plan, expectations_from_counts
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from task59_h4_gc_no_frame_fit import CKPT_PATH, GPI2_SELECTED
from task60_ionq_no_frame_h4 import fit_all_slots, build_full_from_independent_states, energy_report, variance_weights
from task69_hw_plan_rehearsal import PANEL_SLOTS, measured_groups_for_slot, corrected_observables_partial
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

K = 6
BACKENDS = ["ideal", "aria-1", "forte-1"]
PROPOSED_SHOTS = 1100
N_TRIALS = 6
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task69b_hw_plan_rehearsal_uncertainty_results.json")


def one_trial(trial, backend_name, state, kept, groups, diagonalizers, n_groups, p, non_id_labels,
              fixed_solutions, P_S, diag):
    postselected = {name: {} for name in kept}
    kept_shots = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for gi in measured_groups_for_slot(name, n_groups):
            dg = diagonalizers[gi]
            full_counts = entry["counts"][gi]
            rng = np.random.default_rng(stable_seed("task69b", backend_name, name, gi, trial))
            resampled = bootstrap_counts(full_counts, PROPOSED_SHOTS, rng)
            filtered = {bs[1:]: c for bs, c in resampled.items() if bs[0] == "0"}
            n_kept = sum(filtered.values())
            exp = expectations_from_counts(filtered, dg) if filtered else {l: 0.0 for l in groups[gi]}
            for l in groups[gi]:
                postselected[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                kept_shots[name][l] = n_kept

    weights_raw = variance_weights(postselected, kept_shots)
    fits_raw = fit_all_slots(postselected, weights_raw, P_S, kept)
    full_raw = build_full_from_independent_states(fits_raw, diag, K, P_S, non_id_labels)
    _, errs_raw = energy_report(p, full_raw, K)

    if backend_name == "ideal":
        corrected = postselected
    else:
        corrected = corrected_observables_partial(postselected, GPI2_SELECTED[backend_name], kept, fixed_solutions)
    weights_c = variance_weights(corrected, kept_shots)
    fits_nf = fit_all_slots(corrected, weights_c, P_S, kept)
    full_nf = build_full_from_independent_states(fits_nf, diag, K, P_S, non_id_labels)
    _, errs_nf = energy_report(p, full_nf, K)

    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    rng_frame = np.random.default_rng(6900 + trial)
    U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, corrected, weight_unit,
                                             rng_frame, n_restarts=4)
    full_sf = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    _, errs_sf = energy_report(p, full_sf, K)

    return errs_raw["err_vs_exact_kcal"], errs_nf["err_vs_exact_kcal"], errs_sf["err_vs_exact_kcal"], float(chi2dof)


def main():
    print("\n" + "=" * 96)
    print(f"  task69b -- {N_TRIALS}-trial uncertainty rehearsal, 56-circuit plan @ {PROPOSED_SHOTS} shots")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    n_groups = len(groups)

    with open(CKPT_PATH) as f:
        state = json.load(f)

    results = {}
    for backend_name in BACKENDS:
        raws, nfs, sfs, chis = [], [], [], []
        for trial in range(N_TRIALS):
            r, n, s, c = one_trial(trial, backend_name, state, kept, groups, diagonalizers, n_groups, p,
                                     non_id_labels, fixed_solutions, P_S, diag)
            raws.append(r); nfs.append(n); sfs.append(s); chis.append(c)
            print(f"    {backend_name} trial {trial+1}/{N_TRIALS}: raw={r:+.3f}  no_frame={n:+.3f}  "
                  f"shared_frame={s:+.3f}  chi2/dof={c:.4f}")
        def stat(arr):
            a = np.abs(arr)
            return {"mean_signed": float(np.mean(arr)), "std": float(np.std(arr, ddof=1)),
                     "mean_abs": float(np.mean(a)), "q50_abs": float(np.percentile(a, 50)),
                     "q90_abs": float(np.percentile(a, 90))}
        results[backend_name] = {"raw": stat(raws), "no_frame": stat(nfs), "shared_frame": stat(sfs),
                                   "mean_chi2_dof": float(np.mean(chis))}
        print(f"\n  == {backend_name} summary (N={N_TRIALS} trials) ==")
        for tag, arr in [("raw", raws), ("no_frame", nfs), ("shared_frame", sfs)]:
            a = np.abs(arr)
            print(f"    {tag:<14}: mean|err|={np.mean(a):.4f}  std={np.std(arr,ddof=1):.4f}  "
                  f"Q50={np.percentile(a,50):.4f}  Q90={np.percentile(a,90):.4f} kcal/mol")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
