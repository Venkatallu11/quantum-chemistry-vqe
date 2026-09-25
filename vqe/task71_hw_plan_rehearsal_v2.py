#!/usr/bin/env python3
"""
task71_hw_plan_rehearsal_v2.py -- iteration 71. Corrected rehearsal of
the real-hardware plan, fixing two real issues caught by external review
of task69/69b:

1. GC-AWARE CORRECTION (task70): task69/69b's correction step used
   analytic_A_and_B_conditioned, which treats the measurement/basis-
   change circuit as noiseless -- exact for GC0 (0 extra gates), WRONG
   for GC1/GC2/GC3 (3/3/5 real extra native 2q gates). Now every
   group's correction uses analytic_A_and_B_conditioned_gc, propagating
   the assumed noise model through that group's own native diagonalizer
   too. Regression-verified to exactly match the old function on GC0.

2. PANEL DIRECTION COVERAGE: task69's 7-slot panel (u_0,u_1,u_2,u_4,
   (u0+u2),(u2+u5),(u4+u5)) never touched Schmidt direction 3 anywhere.
   Replaced with a 6-slot panel chosen by an explicit coverage
   constraint (all 6 directions represented at least once), which turned
   out CHEAPER (6 slots instead of 7): u_0, u_1, u_2, (u2+u5), (u4+u5),
   (u1+u3). 54 total circuits (51 new + 3 reused from task49), real cost
   $2,259.00 at the precise Azure-confirmed rate ($0.0001645/1q,
   $0.001121/2q), margin $93.16 of $2,352.16.

Same method otherwise: bootstrap already-collected real task59 free-sim
data down to 1,100 shots/circuit, 6 independent trials, all three locked
analyses (raw / no-frame / shared-frame).

Run:
    PYTHONHASHSEED=0 python vqe/task71_hw_plan_rehearsal_v2.py
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
from task70_gc_aware_correction import analytic_A_and_B_conditioned_gc
from task37b_h4_noise_model import GPI_REAL_MEAN
from native_stateprep import to_native
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
BACKENDS = ["ideal", "aria-1", "forte-1"]
PROPOSED_SHOTS = 1100
N_TRIALS = 4
PANEL_SLOTS = ["u_0", "u_1", "u_2", "(u2+u5)", "(u4+u5)", "(u1+u3)"]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task71_hw_plan_rehearsal_v2_results.json")


def measured_groups_for_slot(name, n_groups):
    return list(range(2)) if name not in PANEL_SLOTS else list(range(n_groups))


_AB_CACHE = {}
def corrected_observables_gc_aware(postselected, p_gpi2, kept, fixed_solutions, diag_natives, diagonalizers,
                                     groups, measured_groups_map):
    corrected = {name: {} for name in kept}
    for name in kept:
        for gi in measured_groups_map[name]:
            key = (name, gi, round(float(p_gpi2), 7))
            if key not in _AB_CACHE:
                A, B, _, _ = analytic_A_and_B_conditioned_gc(
                    fixed_solutions[name]["angles"], GATE_NAME, ZZ_ASSUMED, GPI_REAL_MEAN, float(p_gpi2), 0.0, 0.0,
                    diag_natives[gi], diagonalizers[gi].transformed)
                _AB_CACHE[key] = (A, B)
            A, B = _AB_CACHE[key]
            for label in groups[gi]:
                raw = float(postselected[name][label])
                denom = float(A[label])
                ratio = float(B[label] / denom) if abs(denom) > 1e-6 else 1.0
                corrected[name][label] = float(np.clip(raw * ratio, -1.0, 1.0))
    return corrected


def one_trial(trial, backend_name, state, kept, groups, diagonalizers, diag_natives, n_groups, p, non_id_labels,
              fixed_solutions, P_S, diag, measured_groups_map):
    postselected = {name: {} for name in kept}
    kept_shots = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for gi in measured_groups_map[name]:
            dg = diagonalizers[gi]
            full_counts = entry["counts"][gi]
            rng = np.random.default_rng(stable_seed("task71", backend_name, name, gi, trial))
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
        corrected = corrected_observables_gc_aware(postselected, GPI2_SELECTED[backend_name], kept, fixed_solutions,
                                                      diag_natives, diagonalizers, groups, measured_groups_map)
    weights_c = variance_weights(corrected, kept_shots)
    fits_nf = fit_all_slots(corrected, weights_c, P_S, kept)
    full_nf = build_full_from_independent_states(fits_nf, diag, K, P_S, non_id_labels)
    _, errs_nf = energy_report(p, full_nf, K)

    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    rng_frame = np.random.default_rng(7100 + trial)
    U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, corrected, weight_unit,
                                             rng_frame, n_restarts=4)
    full_sf = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    _, errs_sf = energy_report(p, full_sf, K)

    return errs_raw["err_vs_exact_kcal"], errs_nf["err_vs_exact_kcal"], errs_sf["err_vs_exact_kcal"], float(chi2dof)


def load_existing():
    if os.path.exists(RESULTS_PATH):
        with open(RESULTS_PATH) as f:
            return json.load(f)
    return {}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=BACKENDS, required=True, help="run ONE backend only, low-memory mode")
    args = ap.parse_args()
    backend_name = args.backend

    print("\n" + "=" * 96)
    print(f"  task71_hw_plan_rehearsal_v2.py -- {backend_name} only, {N_TRIALS} trials @ {PROPOSED_SHOTS} shots")
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
    diag_natives = [to_native(d.to_circuit(), GATE_NAME) for d in diagonalizers]
    measured_groups_map = {name: measured_groups_for_slot(name, n_groups) for name in kept}

    n_circuits = sum(len(v) for v in measured_groups_map.values())
    print(f"  design: {n_circuits} circuits (panel: {PANEL_SLOTS})")

    with open(CKPT_PATH) as f:
        state = json.load(f)

    raws, nfs, sfs, chis = [], [], [], []
    for trial in range(N_TRIALS):
        r, n, s, c = one_trial(trial, backend_name, state, kept, groups, diagonalizers, diag_natives, n_groups,
                                 p, non_id_labels, fixed_solutions, P_S, diag, measured_groups_map)
        raws.append(r); nfs.append(n); sfs.append(s); chis.append(c)
        print(f"    {backend_name} trial {trial+1}/{N_TRIALS}: raw={r:+.3f}  no_frame={n:+.3f}  "
              f"shared_frame={s:+.3f}  chi2/dof={c:.4f}")

    def stat(arr):
        a = np.abs(arr)
        return {"mean_signed": float(np.mean(arr)), "std": float(np.std(arr, ddof=1)),
                 "mean_abs": float(np.mean(a)), "q50_abs": float(np.percentile(a, 50)),
                 "q90_abs": float(np.percentile(a, 90))}
    entry = {"raw": stat(raws), "no_frame": stat(nfs), "shared_frame": stat(sfs), "mean_chi2_dof": float(np.mean(chis))}
    print(f"\n  == {backend_name} summary (N={N_TRIALS}) ==")
    for tag, arr in [("raw", raws), ("no_frame", nfs), ("shared_frame", sfs)]:
        a = np.abs(arr)
        print(f"    {tag:<14}: mean|err|={np.mean(a):.4f}  std={np.std(arr,ddof=1):.4f}  "
              f"Q50={np.percentile(a,50):.4f}  Q90={np.percentile(a,90):.4f} kcal/mol")

    results = load_existing()
    results[backend_name] = entry
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
