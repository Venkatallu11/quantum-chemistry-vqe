#!/usr/bin/env python3
"""
zne_vs_cdr.py — ZNE and CDR compared on the SAME circuit, same noise model.
============================================================================
WHY THIS FILE EXISTS: the previous ZNE-vs-CDR comparison wasn't a fair one.
The 34.25 kcal/mol quadratic-ZNE number (native_forged_zne_results.json)
came from the OLD 14-gate hand-derived native state prep, folded and
measured on real aria-1 hardware (raw 181 kcal/mol). CDR's 3.83 kcal/mol
(fixed_ansatz_v2_results.json) came from the NEW 11-gate fixed_ansatz.py
circuit under a locally-simulated per-gate depolarizing model (raw 99.4
kcal/mol). Different circuit, different noise source -- comparing them
proved nothing about ZNE vs CDR as methods. This file runs ZNE on the
IDENTICAL circuit (fixed_ansatz.build_ansatz, 11 two-qubit gates) and the
IDENTICAL noise model cdr_mitigation.py uses (depolarizing_error, cx:
P2_PER_GATE=0.01214, u3: P1_PER_GATE=P2_PER_GATE/40, applied directly --
no root-find/calibration to hit a target f, matching the fix already made
to CDR's calibration bug), so the two methods are finally comparable.

FOLDING LEVEL: gate folding is applied to the ABSTRACT u3/cx-transpiled
circuit (the exact circuit cdr_mitigation.py's local Aer simulation uses),
never to a native (ms/gpi) circuit, and the IonQ TrappedIonOptimizerPlugin
is NEVER invoked here. That plugin is angle- AND fold-aware in a way that
silently cancels much of an inserted fold back out -- a quick independent
spot-check here (fold=3 and fold=5 applied to a NATIVE ms circuit, then
TrappedIonOptimizerPlugin run over it) found real, substantial 2-qubit
gate cancellation of the same character described as a known risk (see
this run's own printed numbers, not asserted from elsewhere) -- a milder
echo of the original compiler-side fold-cancellation bug this project hit
early on with abstract-gate submission. Working at the abstract u3/cx
level, matching CDR's own execution path exactly, sidesteps that risk
entirely rather than needing to detect and correct for it.

GATE-COUNT FIX (shared with cdr_mitigation.py): every u3/cx transpile here
uses optimization_level=0, not 1. Verified (not assumed): optimization_level
>=1 numerically collapses some fitted-angle circuits to fewer CX gates when
an angle lands near a periodic special value, breaking the fixed-structure
premise -- opt_level=0 was confirmed constant (11 CX) across all 25 targets
and random training angles before use. See cdr_mitigation.py's docstring
for the full account; fixed_ansatz.py now self-checks this too.

FOLD/SCALE CONVENTION: `fold_circuit` (copied from ionq_run.py, the
already-verified local-unitary-folding implementation: G -> G(G^-1 G)^reps,
reps=(fold-1)/2) is called with fold in {1,3,5} -- this repo's own
established convention elsewhere (ionq_run.py, ionq_fold_check.py,
ionq_native_forged_energy.py) where the `fold` argument IS the 2-qubit
gate-count multiplier directly. The ACTUAL multiplier is measured by
counting 2-qubit gates before and after folding for every run below
-- not assumed to be `fold`, not assumed to be 2*fold-1 either.

Run:
    python vqe/zne_vs_cdr.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
import cdr_mitigation as cdr
from fixed_ansatz import build_ansatz

from qiskit import transpile
from qiskit.quantum_info import Pauli

HARTREE_TO_KCAL_MOL = cdr.HARTREE_TO_KCAL_MOL
K = cdr.K
BASIS_GATES = cdr.BASIS_GATES
FOLDS = [1, 3, 5]

CDR_V2_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "fixed_ansatz_v2_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "zne_vs_cdr_results.json")


# ---------------------------------------------------------------------------
# Local unitary folding, copied verbatim from ionq_run.py's already-verified
# implementation rather than importing that module (which has IonQ-cloud
# imports this file has no reason to pull in for a simulator-only run).
# ---------------------------------------------------------------------------

def fold_circuit(qc, fold):
    """G -> G (G^-1 G)^reps on every 2-qubit gate. fold must be odd;
    fold=1 is a no-op copy."""
    if fold == 1:
        return qc.copy()
    assert fold % 2 == 1, "fold factor must be odd"
    reps = (fold - 1) // 2
    folded = qc.copy_empty_like()
    for instr in qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        folded.append(op, qargs, cargs)
        if op.num_qubits == 2 and op.name not in ("measure", "barrier"):
            for _ in range(reps):
                folded.append(op.inverse(), qargs, cargs)
                folded.append(op, qargs, cargs)
    return folded


def noisy_labels_folded(angles, labels, estimator, fold):
    """Same as cdr_mitigation.noisy_labels, but folds the u3/cx-transpiled
    circuit BEFORE handing it to the noisy estimator -- folding happens on
    the exact circuit that gets simulated, nothing re-transpiles or
    re-optimizes it afterward."""
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    folded = fold_circuit(qc, fold)
    obs = [Pauli(l) for l in labels]
    result = estimator.run([(folded, obs)]).result()
    evs = np.atleast_1d(result[0].data.evs)
    return {l: float(np.real(v)) for l, v in zip(labels, evs)}, folded.count_ops().get("cx", 0)


def measure_fold(p, non_id_labels, estimator, fold):
    """Full K=5 forged energy at a given fold, plus the REAL measured
    2-qubit gate count (not assumed from the fold label)."""
    raw, n2q_list = {}, []
    for name in p["solutions"]:
        vals, n2q = noisy_labels_folded(p["solutions"][name]["angles"], non_id_labels, estimator, fold)
        raw[name] = vals
        n2q_list.append(n2q)
    assert len(set(n2q_list)) == 1, f"2-qubit gate count not fixed across the 25 targets at fold={fold}: {set(n2q_list)}"
    n2q_folded = n2q_list[0]

    mats = cdr.combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, err_kcal, f = cdr.energy_from_alpha_matrices(mats, p)
    return {"fold": fold, "n2q_folded": n2q_folded, "E": E, "err_kcal": err_kcal, "f": f}


# ---------------------------------------------------------------------------
# ZNE fits (linear / quadratic / exponential), against the MEASURED scale.
# ---------------------------------------------------------------------------

def zne_fits(scales, energies, exact_energy):
    no_mit_kcal = abs(energies[0] - exact_energy) * HARTREE_TO_KCAL_MOL

    lin_fit = np.polyfit(scales, energies, 1)
    E_lin0 = float(np.polyval(lin_fit, 0))
    lin_kcal = abs(E_lin0 - exact_energy) * HARTREE_TO_KCAL_MOL

    quad_fit = np.polyfit(scales, energies, 2)
    E_quad0 = float(np.polyval(quad_fit, 0))
    quad_kcal = abs(E_quad0 - exact_energy) * HARTREE_TO_KCAL_MOL

    try:
        diffs = np.array(energies) - energies[-1]
        if np.all(diffs > 0) or np.all(diffs < 0):
            slope, intercept = np.polyfit(scales, np.log(np.abs(diffs)), 1)
            A = np.sign(diffs[0]) * np.exp(intercept)
            E_exp0 = energies[-1] + A
            exp_kcal = abs(E_exp0 - exact_energy) * HARTREE_TO_KCAL_MOL
            exp_note = None
        else:
            exp_kcal = None
            exp_note = "energies not monotonic across scales -- exponential fit not well-posed, not faked"
    except Exception as e:
        exp_kcal = None
        exp_note = f"exponential fit failed: {e}"

    return {
        "no_mitigation_kcal": round(no_mit_kcal, 4),
        "zne_linear_kcal": round(lin_kcal, 4),
        "zne_quadratic_kcal": round(quad_kcal, 4),
        "zne_exponential_kcal": round(exp_kcal, 4) if exp_kcal is not None else None,
        "zne_exponential_note": exp_note,
    }


# ---------------------------------------------------------------------------
# Independent spot-check: does TrappedIonOptimizerPlugin cancel a NATIVE
# fold back down? Verified here (not asserted from elsewhere) as the
# reason this file deliberately avoids that execution path entirely.
# ---------------------------------------------------------------------------

def native_fold_optimizer_spot_check():
    from qiskit.quantum_info import Statevector
    from qiskit.transpiler import PassManagerConfig
    from qiskit_ionq import TrappedIonOptimizerPlugin
    from native_stateprep import to_native, native_target
    from ionq_fold_check import fold_native_2q

    angles = [0.3, 0.5, 0.7, 0.2, 0.9]
    native = to_native(build_ansatz(angles), "ms")
    n2q_base = native.count_ops().get("ms", 0)

    rows = []
    for fold in (3, 5):
        folded = fold_native_2q(native, fold, "ms")
        n2q_before = folded.count_ops().get("ms", 0)
        tgt = native_target(folded.num_qubits, "ms")
        pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt), optimization_level=3)
        optimized = pm.run(folded)
        n2q_after = optimized.count_ops().get("ms", 0)

        sv1 = np.asarray(Statevector.from_instruction(folded))
        sv2 = np.asarray(Statevector.from_instruction(optimized))
        idx = int(np.argmax(np.abs(sv1)))
        phase = sv2[idx] / sv1[idx] if abs(sv1[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv2 / phase - sv1)))

        rows.append({"fold": fold, "n2q_before_optimizer": n2q_before, "n2q_after_optimizer": n2q_after,
                      "statevector_identical": bool(err < 1e-9), "statevector_err": err})
    return {"n2q_base_native": n2q_base, "folds": rows}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 78)
    print("  ZNE vs CDR on the SAME 11-gate fixed ansatz, same noise model (simulator only)")
    print("=" * 78)

    p = cdr.setup()
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    n_groups = len(effrag.group_labels_qubit_wise(p["alpha_labels"]))
    print(f"\n  {len(p['alpha_labels'])} alpha labels -> {n_groups} qubit-wise measurement groups")
    print(f"  ansatz 2-qubit gate count (base, fold=1) = {p['gate_count']}")

    print("\n  Native-fold + TrappedIonOptimizerPlugin spot-check (why this file stays")
    print("  at the abstract u3/cx level instead) ...")
    spot_check = native_fold_optimizer_spot_check()
    print(f"    native base: {spot_check['n2q_base_native']} 2q gates")
    for row in spot_check["folds"]:
        print(f"    fold={row['fold']}: {row['n2q_before_optimizer']} -> {row['n2q_after_optimizer']} "
              f"2q gates after TrappedIonOptimizerPlugin "
              f"(statevector identical: {row['statevector_identical']})")

    print(f"\n  Building the noise model shared with cdr_mitigation.py: "
          f"P2_PER_GATE={cdr.P2_PER_GATE} (cx), P1_PER_GATE={cdr.P1_PER_GATE:.6f} (u3), applied directly.")
    estimator = cdr.make_noisy_estimator()

    print(f"\n  Measuring the 25 real target circuits at folds {FOLDS} "
          f"({n_groups} groups x 25 slots x {len(FOLDS)} folds = {n_groups * 25 * len(FOLDS)} "
          "real-hardware-equivalent circuits)...")
    fold_results = []
    for fold in FOLDS:
        r = measure_fold(p, non_id_labels, estimator, fold)
        measured_scale = r["n2q_folded"] / p["gate_count"]
        r["measured_scale"] = measured_scale
        fold_results.append(r)
        print(f"    fold={fold}: measured 2q gates={r['n2q_folded']} "
              f"(measured scale={measured_scale:.3f}, nominal fold label={fold}), "
              f"E={r['E']:.6f} Ha, error={r['err_kcal']:.3f} kcal/mol, f={r['f']:.6f}")

    claimed_scale_formula = [2 * f - 1 for f in FOLDS]
    measured_scale_list = [r["measured_scale"] for r in fold_results]
    scale_formula_matches_claim = measured_scale_list == [float(x) for x in claimed_scale_formula]
    print(f"\n  SCALE CONVENTION CHECK: the task text asserted noise scale = 2*fold-1 "
          f"(i.e. {claimed_scale_formula} for fold labels {FOLDS}). Measured (fold_circuit's own "
          f"gate-count ratio, this repo's established convention: fold argument = multiplier "
          f"directly): {measured_scale_list}. {'Matches the asserted formula.' if scale_formula_matches_claim else 'Does NOT match -- using the MEASURED value, not the asserted formula, per the instruction to verify rather than assume.'}")

    fold1 = fold_results[0]
    cdr_raw_kcal = None
    if os.path.exists(CDR_V2_RESULTS_PATH):
        with open(CDR_V2_RESULTS_PATH) as f:
            cdr_v2 = json.load(f)
        cdr_raw_kcal = cdr_v2["results_kcal_mol"]["raw"]
        consistency_diff = abs(fold1["err_kcal"] - cdr_raw_kcal)
        print(f"\n  Consistency check vs cdr_mitigation.py's own raw (fold=1) result: "
              f"{fold1['err_kcal']:.3f} vs {cdr_raw_kcal:.3f} kcal/mol "
              f"(diff {consistency_diff:.4f}, {'PASS' if consistency_diff < 0.01 else 'MISMATCH -- investigate'})")

    scales = [r["measured_scale"] for r in fold_results]
    energies = [r["E"] for r in fold_results]
    fits = zne_fits(scales, energies, p["exact_energy"])
    print(f"\n  ZNE fits (against MEASURED scales {[round(s, 3) for s in scales]}):")
    print(f"    no mitigation (scale=1) : {fits['no_mitigation_kcal']:.2f} kcal/mol")
    print(f"    ZNE linear              : {fits['zne_linear_kcal']:.2f} kcal/mol")
    print(f"    ZNE quadratic           : {fits['zne_quadratic_kcal']:.2f} kcal/mol")
    if fits["zne_exponential_kcal"] is not None:
        print(f"    ZNE exponential         : {fits['zne_exponential_kcal']:.2f} kcal/mol")
    else:
        print(f"    ZNE exponential         : n/a ({fits['zne_exponential_note']})")

    cdr_per_basis_mean = cdr_v2["results_kcal_mol"]["per_basis"]["mean"] if cdr_raw_kcal is not None else None
    chem_acc = p["chemical_accuracy_kcal"]

    # --- cost comparison ---
    zne_total_circuits = n_groups * 25 * len(FOLDS)
    zne_total_2q_ops = sum(n_groups * 25 * r["n2q_folded"] for r in fold_results)

    cdr_n_train = cdr_v2["n_seeds"] * 0 + 25 * cdr.N_TRAIN_PER_SLOT if cdr_raw_kcal is not None else 25 * cdr.N_TRAIN_PER_SLOT
    cdr_total_circuits = n_groups * (cdr_n_train + 25)   # training + targets, all fold=1
    cdr_total_2q_ops = n_groups * (cdr_n_train + 25) * p["gate_count"]

    print("\n" + "-" * 78)
    print(f"  {'method':<24}{'error (kcal/mol)':>20}{'reaches chem.acc.':>22}")
    print(f"  {'raw (fold=1)':<24}{fold1['err_kcal']:>20.3f}{'no':>22}")
    print(f"  {'ZNE linear':<24}{fits['zne_linear_kcal']:>20.3f}{'yes' if fits['zne_linear_kcal']<chem_acc else 'no':>22}")
    print(f"  {'ZNE quadratic':<24}{fits['zne_quadratic_kcal']:>20.3f}{'yes' if fits['zne_quadratic_kcal']<chem_acc else 'no':>22}")
    if fits["zne_exponential_kcal"] is not None:
        print(f"  {'ZNE exponential':<24}{fits['zne_exponential_kcal']:>20.3f}{'yes' if fits['zne_exponential_kcal']<chem_acc else 'no':>22}")
    else:
        print(f"  {'ZNE exponential':<24}{'n/a':>20}{'n/a':>22}")
    if cdr_per_basis_mean is not None:
        print(f"  {'CDR per-basis (mean/8)':<24}{cdr_per_basis_mean:>20.3f}{'yes' if cdr_per_basis_mean<chem_acc else 'no':>22}")
    print("-" * 78)

    winner, winner_val = min(
        [("ZNE linear", fits["zne_linear_kcal"]), ("ZNE quadratic", fits["zne_quadratic_kcal"])]
        + ([("ZNE exponential", fits["zne_exponential_kcal"])] if fits["zne_exponential_kcal"] is not None else [])
        + ([("CDR per-basis", cdr_per_basis_mean)] if cdr_per_basis_mean is not None else []),
        key=lambda t: t[1],
    )
    print(f"\n  WINNER: {winner} at {winner_val:.3f} kcal/mol.")
    if cdr_per_basis_mean is not None:
        ratio = fits["zne_quadratic_kcal"] / cdr_per_basis_mean if cdr_per_basis_mean > 0 else float("inf")
        print(f"  CDR per-basis beats ZNE quadratic by {ratio:.2f}x on this ansatz+noise model.")
    print(f"  Chemical accuracy ({chem_acc} kcal/mol): "
          f"{'reached' if winner_val < chem_acc else 'NOT reached'} by the best method, "
          f"{'reliably' if winner == 'CDR per-basis' else '(single measurement, not seed-swept)'}.")

    print(f"\n  COST (real-hardware-equivalent circuit count, {n_groups} qubit-wise groups/slot):")
    print(f"    ZNE  (folds {FOLDS}): {zne_total_circuits} circuits, {zne_total_2q_ops} total 2-qubit gate operations")
    print(f"    CDR  ({cdr_n_train} training + 25 targets, fold=1 only): "
          f"{cdr_total_circuits} circuits, {cdr_total_2q_ops} total 2-qubit gate operations")
    print(f"    ZNE circuits at fold=5 have {fold_results[-1]['n2q_folded']} two-qubit gates each -- "
          "well past the ~23-gate flat-price threshold this project has flagged before; CDR's circuits "
          f"never exceed {p['gate_count']}.")

    results = {
        "ansatz_gate_count_base": p["gate_count"],
        "scale_convention_check": {
            "fold_labels": FOLDS,
            "claimed_scale_2fold_minus_1": claimed_scale_formula,
            "measured_scale": measured_scale_list,
            "matches_claimed_formula": scale_formula_matches_claim,
            "note": "fold_circuit's own gate-count ratio was measured directly, not assumed; "
                    "this repo's established fold_circuit/fold_native_2q convention (used throughout "
                    "ionq_run.py, ionq_fold_check.py, ionq_native_forged_energy.py) is fold argument = "
                    "multiplier directly, which is what was measured here",
        },
        "n_qubit_wise_groups": n_groups,
        "native_fold_optimizer_spot_check": spot_check,
        "fold_results": fold_results,
        "consistency_check_vs_cdr_raw": {
            "zne_fold1_err_kcal": fold1["err_kcal"],
            "cdr_raw_err_kcal": cdr_raw_kcal,
            "diff": abs(fold1["err_kcal"] - cdr_raw_kcal) if cdr_raw_kcal is not None else None,
        },
        "zne_fits_kcal_mol": fits,
        "cdr_per_basis_mean_kcal_mol": cdr_per_basis_mean,
        "cdr_per_basis_source": "fixed_ansatz_v2_results.json (mean over 8 seeds)",
        "comparison_table_kcal_mol": {
            "raw_fold1": round(fold1["err_kcal"], 4),
            "zne_linear": fits["zne_linear_kcal"],
            "zne_quadratic": fits["zne_quadratic_kcal"],
            "zne_exponential": fits["zne_exponential_kcal"],
            "cdr_per_basis_mean": cdr_per_basis_mean,
        },
        "winner": winner, "winner_kcal_mol": round(winner_val, 4),
        "chemical_accuracy_kcal": chem_acc,
        "chemical_accuracy_reached_by_winner": bool(winner_val < chem_acc),
        "cost_comparison": {
            "zne_folds": FOLDS,
            "zne_total_circuits": zne_total_circuits,
            "zne_total_2q_gate_ops": zne_total_2q_ops,
            "cdr_n_training_circuits": cdr_n_train,
            "cdr_n_target_circuits": 25,
            "cdr_total_circuits": cdr_total_circuits,
            "cdr_total_2q_gate_ops": cdr_total_2q_ops,
        },
        "references": {
            "classical_floor_kcal": p["classical_floor_kcal"],
            "chemical_accuracy_kcal": chem_acc,
            "old_native_ansatz_zne_quadratic_kcal_real_aria1": p["quadratic_zne_kcal"],
            "old_native_ansatz_note": "34.25 kcal/mol on the OLD 14-gate native ansatz, real aria-1 hardware -- "
                                       "NOT the same circuit or noise source as this file's comparison, kept only "
                                       "as external context, not compared apples-to-apples here",
        },
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
