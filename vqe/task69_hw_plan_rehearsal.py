#!/usr/bin/env python3
"""
task69_hw_plan_rehearsal.py -- iteration 69. Dry-run rehearsal of the
proposed 56-circuit real-hardware plan (21x GC0 + 21x GC1, all 21 slots;
GC2+GC3 completed on 7 pre-registered panel slots: u_0, u_1, u_2, u_4,
(u0+u2), (u2+u5), (u4+u5)) BEFORE spending any real money -- exactly
this project's own established discipline (validate on free simulators
first).

METHOD: task59's real, already-collected free-simulator data (all 84
circuits x 3 backends, 20,000 shots each, real submissions to
ionq_simulator) is bootstrap-resampled DOWN to 1,100 shots/circuit --
the real shot count proposed for the actual hardware run -- restricted
to exactly the 56-circuit design's (slot, group) set. This tells us,
using REAL noise-model data at the REAL proposed shot count, whether
the three locked analyses (raw / no-frame / shared-frame) behave
sensibly BEFORE any real hardware dollar is spent. Zero new cost.

THREE ANALYSES, locked in advance, matching IonQ's Sep 24 request:
  A. RAW: postselected, UNCORRECTED expectations -> independent
     per-slot 5-parameter physical-state fit (task60's
     fit_slot_data_only), no shared frame.
  B. NO-FRAME, corrected: same as A but with the conditioned PEC+GPi2
     correction applied first (task59/task60's own ratio correction).
  C. SHARED-FRAME, corrected: same corrected data, joint 15-parameter
     frame fit (task36/task66).
For the 14 slots outside the 7-slot panel, only GC0+GC1 labels are
measured; per-slot physical fits STILL produce a prediction for every
label (model-inferred from the 5 fitted state parameters), never
backfilled from simulator data -- exactly the "missing measurements"
handling promised to IonQ.

Run:
    PYTHONHASHSEED=0 python vqe/task69_hw_plan_rehearsal.py
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
from task60_ionq_no_frame_h4 import (
    fit_all_slots, build_full_from_independent_states, energy_report,
    variance_weights,
)
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

K = 6
BACKENDS = ["ideal", "aria-1", "forte-1"]
PROPOSED_SHOTS = 1100
PANEL_SLOTS = ["u_0", "u_1", "u_2", "u_4", "(u0+u2)", "(u2+u5)", "(u4+u5)"]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task69_hw_plan_rehearsal_results.json")


def measured_groups_for_slot(name, n_groups):
    """0,1 (GC0/GC1) always; 2,3 (GC2/GC3) only for the 7 panel slots."""
    return list(range(2)) if name not in PANEL_SLOTS else list(range(n_groups))


def corrected_observables_partial(postselected, p_gpi2, kept, fixed_solutions):
    """Same ratio correction as task60.corrected_observables, but only over
    the labels ACTUALLY measured per slot -- task60's own version assumes
    full 35-label coverage everywhere (true for its original QWC checkpoint),
    which KeyErrors here since 14 of our 21 slots only measured GC0+GC1's 20
    labels. No backfill for the missing ones; they stay genuinely absent,
    to be model-inferred later by the per-slot physical-state fit."""
    from task60_ionq_no_frame_h4 import _CACHE, GATE_NAME, ZZ_ASSUMED
    from task37b_h4_noise_model import GPI_REAL_MEAN
    from task39e_conditioned_correction import analytic_A_and_B_conditioned
    corrected = {name: {} for name in kept}
    for name in kept:
        measured_labels = list(postselected[name].keys())
        key = (name, round(float(p_gpi2), 7))
        if key not in _CACHE:
            A, B, _, _ = analytic_A_and_B_conditioned(
                fixed_solutions[name]["angles"], GATE_NAME, ZZ_ASSUMED, GPI_REAL_MEAN,
                float(p_gpi2), 0.0, 0.0, measured_labels)
            _CACHE[key] = (A, B)
        A, B = _CACHE[key]
        for label in measured_labels:
            raw = float(postselected[name][label])
            denom = float(A[label])
            ratio = float(B[label] / denom) if abs(denom) > 1e-6 else 1.0
            corrected[name][label] = float(np.clip(raw * ratio, -1.0, 1.0))
    return corrected


def main():
    print("\n" + "=" * 96)
    print("  task69_hw_plan_rehearsal.py -- rehearsing the 56-circuit real-hardware plan on REAL free-sim data")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    if not os.path.exists(CKPT_PATH):
        print(f"  ERROR: {CKPT_PATH} not found.")
        return

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    n_groups = len(groups)

    n_circuits = sum(len(measured_groups_for_slot(name, n_groups)) for name in kept)
    print(f"  design: {n_circuits} circuits (21x GC0 + 21x GC1 + GC2/GC3 on {len(PANEL_SLOTS)} panel slots)")
    print(f"  bootstrapping real task59 free-sim data DOWN to {PROPOSED_SHOTS} shots/circuit (the real proposed hw shot count)")

    with open(CKPT_PATH) as f:
        state = json.load(f)

    results = {}
    for backend_name in BACKENDS:
        rng_seed_base = stable_seed("task69_rehearsal", backend_name)
        postselected = {name: {} for name in kept}
        kept_shots = {name: {} for name in kept}
        for name in kept:
            entry = state["done"][f"{backend_name}|{name}"]
            for gi in measured_groups_for_slot(name, n_groups):
                dg = diagonalizers[gi]
                full_counts = entry["counts"][gi]
                rng = np.random.default_rng(stable_seed("task69", backend_name, name, gi))
                resampled = bootstrap_counts(full_counts, PROPOSED_SHOTS, rng)
                filtered = {bs[1:]: c for bs, c in resampled.items() if bs[0] == "0"}
                n_kept = sum(filtered.values())
                exp = expectations_from_counts(filtered, dg) if filtered else {l: 0.0 for l in groups[gi]}
                for l in groups[gi]:
                    postselected[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                    kept_shots[name][l] = n_kept

        mean_labels_per_slot = float(np.mean([len(v) for v in postselected.values()]))
        print(f"\n  == {backend_name} == mean labels/slot measured: {mean_labels_per_slot:.1f} (of 35 total)")

        # ---- A. RAW: no correction ----
        weights_raw = variance_weights(postselected, kept_shots)
        fits_raw = fit_all_slots(postselected, weights_raw, P_S, kept)
        full_raw = build_full_from_independent_states(fits_raw, diag, K, P_S, non_id_labels)
        E_raw, errs_raw = energy_report(p, full_raw, K)

        # ---- B. NO-FRAME, corrected ----
        if backend_name == "ideal":
            corrected = postselected
        else:
            corrected = corrected_observables_partial(postselected, GPI2_SELECTED[backend_name], kept, fixed_solutions)
        weights_c = variance_weights(corrected, kept_shots)
        fits_nf = fit_all_slots(corrected, weights_c, P_S, kept)
        full_nf = build_full_from_independent_states(fits_nf, diag, K, P_S, non_id_labels)
        E_nf, errs_nf = energy_report(p, full_nf, K)

        # ---- C. SHARED-FRAME, corrected ----
        weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
        rng_frame = np.random.default_rng(69)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, corrected, weight_unit,
                                                 rng_frame, n_restarts=4)
        full_sf = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        E_sf, errs_sf = energy_report(p, full_sf, K)

        print(f"    A. raw            : E={E_raw:.6f} Ha  err={errs_raw['err_vs_exact_kcal']:+.4f} kcal/mol")
        print(f"    B. no-frame (corr): E={E_nf:.6f} Ha  err={errs_nf['err_vs_exact_kcal']:+.4f} kcal/mol")
        print(f"    C. shared-frame   : E={E_sf:.6f} Ha  err={errs_sf['err_vs_exact_kcal']:+.4f} kcal/mol  chi2/dof={chi2dof:.3f}")

        results[backend_name] = {
            "raw": {"E": E_raw, "err_kcal": errs_raw["err_vs_exact_kcal"]},
            "no_frame": {"E": E_nf, "err_kcal": errs_nf["err_vs_exact_kcal"]},
            "shared_frame": {"E": E_sf, "err_kcal": errs_sf["err_vs_exact_kcal"], "chi2_dof": float(chi2dof)},
            "mean_labels_per_slot": mean_labels_per_slot,
        }

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
