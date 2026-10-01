#!/usr/bin/env python3
"""
task72_rehearsal_sign_and_chi2_check.py -- iteration 72. Two checks on the
task71 hardware-plan rehearsal before it goes into the plan sent for review.

1. SCHMIDT SIGN CONVENTION. Re-running task71 on numpy 2.4.6 / scipy 1.17.1
   reproduced raw and no-frame exactly but gave shared-frame errors of ~32
   kcal/mol instead of 0.23. Cause: np.linalg.svd fixes each Schmidt pair only
   up to a joint sign, and vectors 2, 3, 5 came out flipped relative to the
   convention the stored task59 data was collected under. That turns
   (u_i+u_j) slots into (u_i-u_j) in the model, and the joint frame fit (a
   rotation near identity) cannot undo it. Independent per-slot fits (raw /
   no-frame) are unaffected. Fixed in qforge.forging.setup_fragment via an
   explicit per-fragment sign fingerprint (SCHMIDT_SIGN_REFERENCE). This
   script GATES on it: the ideal-backend data must fit the model at the exact
   frame (theta=0) to shot-noise level before any analysis runs.

   This is also an exact-state dependency of the shared-frame analysis: it
   uses the FCI Schmidt vectors including their signs, not just their span.

2. REAL CHI2. task71's printed "chi2/dof" (~0.001) is the frame fit's cost
   with unit weights -- a mean squared residual, not a chi2. Here residuals
   of the same fits are scored against shot-noise sigmas of the CORRECTED
   observables (raw binomial sigma times the correction ratio B/A, i.e.
   including the noise amplification of the correction). Labels whose
   correction ratio |B/A| < FORCED_ZERO_RATIO (the correction suppresses the
   measured value >100x, forcing it to ~0 whatever was measured -- 23-24
   of ~516 labels per noisy trial, e.g. IZZI/ZIIZ on (u1+u2), ratio 6e-4)
   carry essentially no measurement and are excluded, with the count
   reported. Left in, their near-zero sigmas make chi2/dof ~1e2-1e3 from just
   6 labels. chi2/dof ~ 1 means the residual scatter is explained by shot
   noise; > 1 means systematic error the shared frame is absorbing.

   Result: chi2/dof 1.03-1.46 on ideal, aria-1 and forte-1 alike --
   residuals consistent with shot noise at 1,100 shots on this noise-model
   data. (Real hardware may not be; that is what the run tests.)

Same data, shot count, trials and seeds as task71.

Run:
    PYTHONHASHSEED=0 python vqe/task72_rehearsal_sign_and_chi2_check.py --backend ideal
    (then aria-1, forte-1; results merge into one JSON)
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task71_hw_plan_rehearsal_v2 as t71
from task36_joint_schmidt_frame import joint_residuals, slot_vector, N_THETA
from general_commuting_measurements import expectations_from_counts
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task72_rehearsal_sign_and_chi2_check_results.json")
EXACT_FRAME_GATE = 1e-3
FORCED_ZERO_RATIO = 0.01   # ideal full-shot mean squared residual at theta=0; unpinned signs give ~0.20


def exact_frame_gate(state, kept, groups, diagonalizers, P_S, K, non_id_labels):
    data = {}
    for name in kept:
        entry = state["done"][f"ideal|{name}"]
        data[name] = {}
        for gi in range(len(groups)):
            filtered = {bs[1:]: c for bs, c in entry["counts"][gi].items() if bs[0] == "0"}
            exp = expectations_from_counts(filtered, diagonalizers[gi])
            for l in groups[gi]:
                data[name][l] = exp.get(l, 0.0)
    unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    r = joint_residuals(np.zeros(N_THETA), np.eye(K), P_S, K, kept, non_id_labels, data, unit)
    return float(r @ r / len(r))


def one_trial(trial, backend_name, state, kept, groups, diagonalizers, diag_natives, p, non_id_labels,
              fixed_solutions, P_S, diag, measured_groups_map):
    K = t71.K
    postselected = {name: {} for name in kept}
    kept_shots = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for gi in measured_groups_map[name]:
            rng = np.random.default_rng(stable_seed("task71", backend_name, name, gi, trial))
            resampled = bootstrap_counts(entry["counts"][gi], t71.PROPOSED_SHOTS, rng)
            filtered = {bs[1:]: c for bs, c in resampled.items() if bs[0] == "0"}
            n_kept = sum(filtered.values())
            exp = expectations_from_counts(filtered, diagonalizers[gi]) if filtered else {l: 0.0 for l in groups[gi]}
            for l in groups[gi]:
                postselected[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                kept_shots[name][l] = n_kept

    if backend_name == "ideal":
        corrected = postselected
    else:
        corrected = t71.corrected_observables_gc_aware(postselected, t71.GPI2_SELECTED[backend_name], kept,
                                                         fixed_solutions, diag_natives, diagonalizers, groups,
                                                         measured_groups_map)

    unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    U_hat, cost, msr = t71.fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, corrected, unit,
                                           np.random.default_rng(7100 + trial), n_restarts=4)
    full_sf = t71.build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    _, errs_sf = t71.energy_report(p, full_sf, K)

    resid = []
    n_forced = 0
    for name in kept:
        v = slot_vector(U_hat, name, K)
        for l, y in corrected[name].items():
            m_raw = postselected[name][l]
            ratio = y / m_raw if abs(m_raw) > 1e-9 else 1.0
            if abs(ratio) < FORCED_ZERO_RATIO:
                # correction ratio B/A ~ 0: the corrected value is forced to ~0 regardless of the
                # measurement, so it has no meaningful shot-noise sigma; left out of chi2 (counted).
                n_forced += 1
                continue
            sigma = np.sqrt(max(1.0 - m_raw ** 2, 1e-4) / max(kept_shots[name][l], 1)) * abs(ratio)
            resid.append((float(np.real(v @ P_S[l] @ v)) - y) / max(sigma, 1e-6))
    resid = np.array(resid)
    chi2_dof = float(resid @ resid / (len(resid) - N_THETA))
    return errs_sf["err_vs_exact_kcal"], float(msr), chi2_dof, len(resid), n_forced


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=t71.BACKENDS, required=True)
    args = ap.parse_args()
    backend_name = args.backend
    K = t71.K

    p = t71.setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, _, _ = t71.fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = t71.kept_slots_for_K(K)
    P_S = t71.build_P_S(p["alpha_labels"], np.asarray(p["u_vecs"]).T)
    groups, diagonalizers = t71.build_general_commuting_measurement_plan(non_id_labels)
    diag_natives = [t71.to_native(d.to_circuit(), t71.GATE_NAME) for d in diagonalizers]
    measured_groups_map = {name: t71.measured_groups_for_slot(name, len(groups)) for name in kept}
    with open(t71.CKPT_PATH) as f:
        state = json.load(f)

    gate = exact_frame_gate(state, kept, groups, diagonalizers, P_S, K, non_id_labels)
    print(f"  exact-frame gate (ideal, full shots): mean sq residual = {gate:.2e}  (limit {EXACT_FRAME_GATE:.0e})")
    assert gate < EXACT_FRAME_GATE, "stored data does not match the model's Schmidt sign convention -- STOP"

    rows = []
    for trial in range(t71.N_TRIALS):
        err, msr, chi2, n, n_forced = one_trial(trial, backend_name, state, kept, groups, diagonalizers, diag_natives, p,
                                      non_id_labels, fixed_solutions, P_S, diag, measured_groups_map)
        rows.append({"shared_frame_err_kcal": err, "unit_weight_msr": msr, "chi2_dof": chi2, "n_resid": n,
                     "n_forced_zero_excluded": n_forced})
        print(f"    {backend_name} trial {trial+1}/{t71.N_TRIALS}: shared_frame={err:+.3f} kcal/mol  "
              f"unit-weight msr={msr:.4f}  chi2/dof={chi2:.2f}  (n={n}, {n_forced} forced-zero excluded)")

    chis = [r["chi2_dof"] for r in rows]
    errs = np.abs([r["shared_frame_err_kcal"] for r in rows])
    entry = {"exact_frame_gate_msr": gate, "trials": rows,
             "shared_frame_mean_abs_kcal": float(np.mean(errs)),
             "chi2_dof_mean": float(np.mean(chis)), "chi2_dof_min": float(np.min(chis)),
             "chi2_dof_max": float(np.max(chis))}
    print(f"  {backend_name}: shared-frame mean|err|={entry['shared_frame_mean_abs_kcal']:.4f} kcal/mol  "
          f"chi2/dof mean={entry['chi2_dof_mean']:.2f} [{entry['chi2_dof_min']:.2f}, {entry['chi2_dof_max']:.2f}]")

    results = {}
    if os.path.exists(RESULTS_PATH):
        with open(RESULTS_PATH) as f:
            results = json.load(f)
    results[backend_name] = entry
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved -> {RESULTS_PATH}")


if __name__ == "__main__":
    main()
