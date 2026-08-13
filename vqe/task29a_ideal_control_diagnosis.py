#!/usr/bin/env python3
"""
task29a_ideal_control_diagnosis.py -- iteration 29, Task A. LOCATE THE
IDEAL-CONTROL FAILURE. All-gate ZNE's ideal control went 1.83 -> 39.31
kcal/mol (Task 28D) even though the fold construction's own unitary-
equivalence check passes at ~1e-14. For an ideal (noiseless) circuit,
E_ideal(lambda) MUST be exactly constant for every valid unitary fold --
it manifestly is not, somewhere downstream of the circuit itself.

Four layers, smallest test first, EXACT statevector throughout (no shots,
no IonQ calls at all -- this isolates the bug from shot noise entirely,
since exact quantities must satisfy |E_lambda - E_1| < 1e-10 if the
pipeline is correct):

  A: state prep only, raw statevector overlap vs fold=1 base
  B: state prep + ONE Pauli observable, exact expectation via the SAME
     measurement-circuit / pauli_expectation() code path production uses,
     just fed exact probabilities instead of sampled counts
  C: full energy reconstruction -- C1 isolates ONE off-diagonal
     (u_n - u_m) synthesis identity (the "measurement/energy-
     reconstruction bookkeeping between folds" prime suspect), C2 builds
     the full K=6 multi-slot energy exactly as production's combine_matrices
     / energy_from_alpha_matrices do
  D: feeds the (near-)constant exact per-fold energies from C2 through the
     SAME two-stage held-out ZNE extrapolator production uses, to test
     whether the extrapolator itself amplifies a near-degenerate input

Reuses optimized_native_circuit / fold_all_gates UNCHANGED from
task28d_all_gate_zne.py -- no reimplementation, no risk of testing a
different circuit than the one that actually failed.

Run:
    python vqe/task29a_ideal_control_diagnosis.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
from task27d_held_out_zne import MODEL_CLASSES, FOLD_STAGE1_FIT, FOLD_STAGE1_HOLDOUT, FOLD_STAGE2_FIT, FOLD_STAGE2_HOLDOUT
import ef_fragment as effrag_mod
from ionq_run import pauli_expectation
from task28d_all_gate_zne import optimized_native_circuit, fold_all_gates
from qiskit.quantum_info import Statevector

K = 6
FOLD_FACTORS = [1, 3, 5, 7, 9]
GATE_NAME = "ms"   # GATE_BY_MODEL["ideal"] = "ms" in task28d
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task29a_ideal_control_diagnosis_results.json")
TOL = 1e-10


def exact_probs(qc):
    sv = Statevector.from_instruction(qc)
    return sv.probabilities_dict()


def main():
    print("\n" + "=" * 96)
    print("  task29a_ideal_control_diagnosis.py -- find the first layer that breaks E_ideal(lambda)=const")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)
    print(f"  K={K}: {len(kept)} kept circuits ({len(diag)} diag + {len(plus)} plus), "
          f"{len(groups)} qubit-wise-commuting measurement groups")

    results = {}
    verdict = None

    # ---------------------------------------------------------------
    # build optimizer-frozen base circuits (same base task28d folds)
    # ---------------------------------------------------------------
    base = {name: optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME) for name in kept}

    # ---------------------------------------------------------------
    # LAYER A -- state prep only, exact statevector, no basis rotation
    # ---------------------------------------------------------------
    print("\n  -- LAYER A: state prep only, raw statevector vs fold=1 --")
    sv_base = {name: np.asarray(Statevector.from_instruction(base[name])) for name in kept}
    layerA = {}
    worst_A = 0.0
    for fold in FOLD_FACTORS:
        worst_this_fold = 0.0
        for name in kept:
            folded = fold_all_gates(base[name], fold, GATE_NAME)
            sv_f = np.asarray(Statevector.from_instruction(folded))
            idx = int(np.argmax(np.abs(sv_base[name])))
            phase = sv_f[idx] / sv_base[name][idx] if abs(sv_base[name][idx]) > 1e-9 else 1.0
            err = float(np.max(np.abs(sv_f / phase - sv_base[name])))
            worst_this_fold = max(worst_this_fold, err)
        layerA[fold] = worst_this_fold
        worst_A = max(worst_A, worst_this_fold)
        print(f"    fold={fold}: worst |sv_folded/phase - sv_base| over {len(kept)} slots = {worst_this_fold:.3e}")
    results["layerA"] = layerA
    layerA_pass = worst_A < TOL
    print(f"  LAYER A: worst={worst_A:.3e}  {'PASS' if layerA_pass else 'FAIL'} (tol={TOL:.0e})")
    if not layerA_pass:
        verdict = "A"

    # ---------------------------------------------------------------
    # LAYER B -- state prep + ONE Pauli observable, exact expectation,
    # through the SAME measurement-circuit / pauli_expectation() path
    # ---------------------------------------------------------------
    print("\n  -- LAYER B: state prep + ONE Pauli observable, exact expectation --")
    probe_name = kept[0]
    probe_group = groups[0]
    probe_label = probe_group[0]
    combined = effrag_mod.combined_basis_label(probe_group)
    print(f"    probe: slot={probe_name}  group={probe_group}  label={probe_label}")
    layerB = {}
    ref_val = None
    worst_B = 0.0
    for fold in FOLD_FACTORS:
        folded = fold_all_gates(base[probe_name], fold, GATE_NAME)
        basis_qc = native_basis_change(combined, GATE_NAME)
        qc = folded.compose(basis_qc)
        probs = exact_probs(qc)
        val = pauli_expectation(probs, probe_label)
        layerB[fold] = val
        if fold == 1:
            ref_val = val
        err = abs(val - ref_val) if ref_val is not None else 0.0
        worst_B = max(worst_B, err)
        print(f"    fold={fold}: <{probe_label}> = {val:.12f}  (|diff vs fold=1| = {err:.3e})")
    results["layerB"] = layerB
    layerB_pass = worst_B < TOL
    print(f"  LAYER B: worst diff={worst_B:.3e}  {'PASS' if layerB_pass else 'FAIL'} (tol={TOL:.0e})")
    if verdict is None and not layerB_pass:
        verdict = "B"

    # ---------------------------------------------------------------
    # exact per-(name, label) expectation table, ALL folds, ALL kept
    # slots x ALL non-identity labels -- built once, reused by C and D
    # ---------------------------------------------------------------
    print("\n  -- building exact per-fold expectation table (all slots, all labels) --")
    m_by_fold = {}
    for fold in FOLD_FACTORS:
        m = {name: {} for name in kept}
        for name in kept:
            folded = fold_all_gates(base[name], fold, GATE_NAME)
            for group in groups:
                combined = effrag_mod.combined_basis_label(group)
                basis_qc = native_basis_change(combined, GATE_NAME)
                qc = folded.compose(basis_qc)
                probs = exact_probs(qc)
                for l in group:
                    m[name][l] = pauli_expectation(probs, l)
        m_by_fold[fold] = m
        print(f"    fold={fold}: built exact expectations for {len(kept)} slots x {len(non_id_labels)} labels")

    # ---------------------------------------------------------------
    # LAYER C1 -- isolate ONE (u_n - u_m) synthesis identity: does the
    # bookkeeping "(un_val + um_val) - plus_val" stay constant across
    # folds, for a single (n, m) pair and a single label?
    # ---------------------------------------------------------------
    print("\n  -- LAYER C1: isolate ONE off-diagonal synthesis identity across folds --")
    n0, m0 = 0, 1
    un_name, um_name, pl_name = f"u_{n0}", f"u_{m0}", f"(u{n0}+u{m0})"
    probe_label_c = next(l for l in non_id_labels if l in m_by_fold[1][un_name]
                          and l in m_by_fold[1][um_name] and l in m_by_fold[1][pl_name])
    print(f"    pair=(u{n0}, u{m0})  label={probe_label_c}")
    layerC1 = {}
    ref_minus = None
    worst_C1 = 0.0
    for fold in FOLD_FACTORS:
        m = m_by_fold[fold]
        un_v, um_v, pl_v = m[un_name][probe_label_c], m[um_name][probe_label_c], m[pl_name][probe_label_c]
        minus_v = (un_v + um_v) - pl_v
        layerC1[fold] = minus_v
        if fold == 1:
            ref_minus = minus_v
        err = abs(minus_v - ref_minus)
        worst_C1 = max(worst_C1, err)
        print(f"    fold={fold}: synthesized (u{n0}-u{m0})[{probe_label_c}] = {minus_v:.12f}  "
              f"(|diff vs fold=1| = {err:.3e})")
    results["layerC1"] = layerC1
    layerC1_pass = worst_C1 < TOL
    print(f"  LAYER C1: worst diff={worst_C1:.3e}  {'PASS' if layerC1_pass else 'FAIL'} (tol={TOL:.0e})")
    if verdict is None and not layerC1_pass:
        verdict = "C1 (off-diagonal synthesis bookkeeping)"

    # ---------------------------------------------------------------
    # LAYER C2 -- full K=6 multi-slot energy reconstruction, exact
    # ---------------------------------------------------------------
    print("\n  -- LAYER C2: full energy reconstruction (all 21 slots, exact) --")
    layerC2 = {}
    ref_E = None
    worst_C2 = 0.0
    for fold in FOLD_FACTORS:
        m = m_by_fold[fold]
        full = {name: dict(m[name]) for name in diag}
        for n in range(K):
            for mm in range(K):
                if n >= mm:
                    continue
                un, um, pl = f"u_{n}", f"u_{mm}", f"(u{n}+u{mm})"
                full[pl] = dict(m[pl])
                synth_minus = {}
                for l in non_id_labels:
                    if l not in m[pl] or l not in full[un] or l not in full[um]:
                        continue
                    synth_minus[l] = (full[un][l] + full[um][l]) - m[pl][l]
                full[f"(u{n}-u{mm})"] = synth_minus
        alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
        E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                              exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
        E_kcal = errs["err_vs_exact_kcal"]
        layerC2[fold] = {"E_hartree": E, "err_vs_exact_kcal": E_kcal}
        if fold == 1:
            ref_E = E
        err = abs(E - ref_E) * HARTREE_TO_KCAL_MOL
        worst_C2 = max(worst_C2, abs(E - ref_E))
        print(f"    fold={fold}: E = {E:.12f} Ha  (err_vs_exact={E_kcal:.4f} kcal/mol, "
              f"|diff vs fold=1|={err:.3e} kcal/mol)")
    results["layerC2"] = layerC2
    layerC2_pass = worst_C2 < TOL
    print(f"  LAYER C2: worst diff={worst_C2:.3e} Ha  {'PASS' if layerC2_pass else 'FAIL'} (tol={TOL:.0e} Ha)")
    if verdict is None and not layerC2_pass:
        verdict = "C2 (full multi-slot energy assembly)"

    # ---------------------------------------------------------------
    # LAYER D -- feed the (near-)constant exact energies through the
    # SAME two-stage held-out ZNE extrapolator production uses
    # ---------------------------------------------------------------
    print("\n  -- LAYER D: exact per-fold energies through the production ZNE extrapolator --")
    curve = {fold: layerC2[fold]["E_hartree"] for fold in FOLD_FACTORS}
    print(f"    input curve (exact, Hartree): {curve}")

    def best_fit(folds_fit, vals_fit, holdout_fold, holdout_val):
        best_name, best_err, best_fn = None, float("inf"), None
        for name, fitter in MODEL_CLASSES.items():
            fn = fitter(folds_fit, vals_fit)
            if fn is None:
                continue
            try:
                pred = float(np.atleast_1d(fn([holdout_fold]))[0])
                if np.isfinite(pred):
                    err = abs(pred - holdout_val)
                    if err < best_err:
                        best_err, best_name, best_fn = err, name, fn
            except Exception:
                continue
        return best_name, best_err, best_fn

    vals_s1 = [curve[f] for f in FOLD_STAGE1_FIT]
    name1, err1, fn1 = best_fit(FOLD_STAGE1_FIT, vals_s1, FOLD_STAGE1_HOLDOUT, curve[FOLD_STAGE1_HOLDOUT])
    print(f"    stage1 (fit [1,3,5], predict 7): best_class={name1}  holdout_err={err1:.3e} Ha")

    vals_s2 = [curve[f] for f in FOLD_STAGE2_FIT]
    name2, err2, fn2 = best_fit(FOLD_STAGE2_FIT, vals_s2, FOLD_STAGE2_HOLDOUT, curve[FOLD_STAGE2_HOLDOUT])
    print(f"    stage2 (fit [1,3,5,7], predict 9): best_class={name2}  holdout_err={err2:.3e} Ha")

    all_folds = sorted(curve.keys())
    vals_all = [curve[f] for f in all_folds]
    fn_final = MODEL_CLASSES[name2](all_folds, vals_all) if name2 else None
    E0 = float(np.atleast_1d(fn_final([0]))[0]) if fn_final is not None else None
    err_vs_E1_kcal = abs(E0 - curve[1]) * HARTREE_TO_KCAL_MOL if E0 is not None else None
    print(f"    ZNE(fold=0) using class={name2}: E0={E0}  |E0 - E(fold=1)| = {err_vs_E1_kcal} kcal/mol")

    results["layerD"] = {
        "stage1_class": name1, "stage1_holdout_err_ha": err1,
        "stage2_class": name2, "stage2_holdout_err_ha": err2,
        "E0_hartree": E0, "E_fold1_hartree": curve[1],
        "err_vs_fold1_kcal": err_vs_E1_kcal,
    }
    layerD_pass = (err_vs_E1_kcal is not None) and (err_vs_E1_kcal * (1.0 / HARTREE_TO_KCAL_MOL) < TOL)
    # report in kcal/mol terms too since 1e-10 Ha is an extremely tight bar for a curve_fit-based extrapolator
    print(f"  LAYER D: |E0 - E(fold=1)| = {err_vs_E1_kcal:.6f} kcal/mol  "
          f"{'PASS' if layerD_pass else 'FAIL'} (tol={TOL*HARTREE_TO_KCAL_MOL:.2e} kcal/mol)")
    if verdict is None and not layerD_pass:
        verdict = "D (extrapolator amplifying near-degenerate input)"

    print("\n" + "=" * 96)
    if verdict is None:
        print("  ALL LAYERS PASS at 1e-10 -- the bug is NOT reproduced by this exact-statevector diagnostic.")
        print("  Implication: the failure must involve shot noise, real-hardware-specific batching/retrieval,")
        print("  or the outage-recovery job reassembly -- NOT the fold/synthesis/extrapolation math itself.")
    else:
        print(f"  FIRST FAILING LAYER: {verdict}")
    print("=" * 96 + "\n")

    results["verdict_first_failing_layer"] = verdict
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
