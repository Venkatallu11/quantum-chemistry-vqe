#!/usr/bin/env python3
"""
task48_gc_native_hw_prep.py -- iteration 48. Prepares and verifies the
GC group-3 ("the hard one" -- 8 labels, the diagonalizer with the most
extra gates: 3 extra 1q / 5 extra 2q abstract) real-hardware circuit for
3 real kept slots, BEFORE any real money is spent.

WHY A NEW VERIFICATION PASS: task46/47 verified the GC diagonalizer only
at the ABSTRACT-gate level (H/S/Sdg/CX via Statevector). Real hardware
submission requires NATIVE gates (GPi/GPi2/ZZ), and native transpilation
is a real, separate translation step (to_native()) that has never been
checked against the GC diagonalizer specifically -- it could introduce
its own bug, exactly the class of bug already found twice in this
integration (task46's qubit-indexing bug, this project's own established
"opt_level>=1 collapse risk" history). Verifying BEFORE submitting,
matching this project's hard-earned discipline of never trusting an
unverified circuit with real money.

CIRCUIT: register (optimized_native_circuit, IDENTICAL to Task 41-44's
own real-hardware circuit construction) + ancilla-parity (native) +
GC group-3 diagonalizer (native-transpiled via to_native(), composed on
qubits [0,1,2,3] only) + measure_all. Uses the STANDARD (not 1q-reduced)
register on purpose: keeps this real-money test inside already-verified
territory rather than combining with the separate, not-yet-jointly-
verified 1q-reduction from Task 44.

VERIFICATION: statevector-equivalence of the NATIVE circuit's GC
reconstruction against the exact ideal Pauli expectations, for all 3
target slots -- must match to <1e-9, matching this project's own
established bar, before anything is queued for real submission.

Run:
    PYTHONHASHSEED=0 python vqe/task48_gc_native_hw_prep.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from general_commuting_measurements import build_general_commuting_measurement_plan
from native_stateprep import to_native
from fixed_ansatz import build_ansatz
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, Pauli

K = 6
GATE_NAME = "zz"
TARGET_SLOTS = ["u_0", "u_1", "u_2"]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task48_gc_native_hw_prep_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task48_gc_native_hw_prep.py -- native-gate GC group-3 circuit, verified before real submission")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    assert n_ok == 36

    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    hard_idx = max(range(len(diagonalizers)), key=lambda i: diagonalizers[i].two_qubit_count)
    group, d = groups[hard_idx], diagonalizers[hard_idx]
    print(f"\n  'hard' GC group (index {hard_idx}): {len(group)} labels, "
          f"abstract 1q={d.one_qubit_count} 2q={d.two_qubit_count}")
    print(f"  group: {group}")

    native_diag = to_native(d.to_circuit(), GATE_NAME)
    native_diag_counts = dict(native_diag.count_ops())
    print(f"  native-transpiled diagonalizer gate counts: {native_diag_counts}")

    ancilla_native = ancilla_cnots_native(GATE_NAME)

    worst_err = 0.0
    circuits = []
    for name in TARGET_SLOTS:
        angles = fixed_solutions[name]["angles"]
        register = optimized_native_circuit(angles, GATE_NAME)
        full5 = with_ancilla_parity_native(register, ancilla_native)
        full5 = full5.compose(native_diag, qubits=[0, 1, 2, 3])
        gate_counts = dict(full5.count_ops())
        print(f"\n  slot={name}: full circuit native gate counts: {gate_counts}")

        # verification: exact ideal reconstruction vs direct exact expectation
        sv5 = Statevector.from_instruction(full5)
        probs5 = sv5.probabilities_dict()
        p_anc1 = sum(v for k, v in probs5.items() if k[0] == "1")
        filtered = {bs[1:]: v for bs, v in probs5.items() if bs[0] == "0"}
        total = sum(filtered.values())
        from general_commuting_measurements import expectations_from_counts
        counts_like = {k: v / total for k, v in filtered.items()}
        recon = expectations_from_counts(counts_like, d)

        sv4 = Statevector.from_instruction(build_ansatz(angles))
        for label in group:
            exact_val = float(np.real(sv4.expectation_value(Pauli(label))))
            err = abs(recon[label] - exact_val)
            worst_err = max(worst_err, err)
        print(f"    p(ancilla=1) ideal: {p_anc1:.3e}  (must be ~0)")
        print(f"    worst reconstruction error this slot: "
              f"{max(abs(recon[l] - float(np.real(sv4.expectation_value(Pauli(l))))) for l in group):.3e}")

        full5_meas = full5.copy()
        full5_meas.measure_all()
        circuits.append(full5_meas)

    print(f"\n  -- OVERALL: worst reconstruction error across all {len(TARGET_SLOTS)} slots: {worst_err:.3e}  "
          f"{'PASS -- SAFE TO SUBMIT' if worst_err < 1e-9 else 'FAIL -- DO NOT SUBMIT, investigate'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "hard_group_index": hard_idx, "group": group,
            "native_diag_counts": native_diag_counts,
            "target_slots": TARGET_SLOTS, "worst_recon_err": worst_err,
            "safe_to_submit": bool(worst_err < 1e-9),
        }, f, indent=2)
    print(f"\n  saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
