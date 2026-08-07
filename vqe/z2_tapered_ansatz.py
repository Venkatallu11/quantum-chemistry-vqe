#!/usr/bin/env python3
"""
z2_tapered_ansatz.py — Task A continued: re-solve a fixed ansatz on the
Z2-tapered (3-qubit) register and report the new two-qubit gate count.
============================================================================
z2_tapering.py verified the taper itself is exact (2.29e-11 kcal/mol,
machine precision) for the 6 diagonal Schmidt vectors. This file extends
the SAME transform to the full 36-target K=6 slot set (6 diagonal + 15
pairs x 2 real-gauge phase states), examines the reduced (3-qubit,
8-dim) target vectors' structure, and reports a real, fitted circuit's
two-qubit gate count -- both a generic StatePreparation baseline and,
if the reduced targets' structure allows it, a hand-derived
lower-gate-count circuit (matching this project's own established
practice, fixed_ansatz.py, of preferring a hand-derived structure once
one is found, but reporting the generic baseline honestly either way).

Run:
    python vqe/z2_tapered_ansatz.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
import ef_fragment as effrag
from z2_tapering import find_register_symmetry, taper_register, verify_round_trip
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector
from qiskit import transpile

HARTREE_TO_KCAL_MOL = 627.5094740631
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "z2_tapered_ansatz_results.json")


def main():
    print("\n" + "=" * 96)
    print("  z2_tapered_ansatz.py -- fit a circuit on the Z2-tapered (3-qubit) alpha register")
    print("=" * 96)

    qop_bare, qop_pen, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    K = int(np.sum(lambdas > 1e-9))
    print(f"  K={K}")

    z2_alpha = find_register_symmetry(alpha_labels, "alpha")
    reduced_diag, sq_qubit, U, tapered_val = taper_register(u_vecs[:K], alpha_labels, z2_alpha, K, "alpha")

    # Build ALL 36 target vectors in the FULL (untapered) 4-qubit space first (real gauge would
    # normally be used here, but this project's real_gauge machinery is defined on the FULL
    # 8-qubit state -- for THIS specific check we work directly with the (possibly complex)
    # u_vecs, consistent with z2_tapering.py's own choice to stay faithful to the untapered
    # Schmidt vectors rather than introduce a different basis convention mid-analysis)
    full_targets = {}
    for n in range(K):
        full_targets[f"u_{n}"] = u_vecs[n]
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        full_targets[f"(u{n}+u{m})"] = (u_vecs[n] + u_vecs[m]) / np.sqrt(2)
        full_targets[f"(u{n}-u{m})"] = (u_vecs[n] - u_vecs[m]) / np.sqrt(2)
    print(f"  {len(full_targets)} total target slots (6 diagonal + 15 pairs x 2)")

    # taper EVERY target the same way, verify uniform tapered-qubit value + exact round trip
    all_vecs = np.array(list(full_targets.values()))
    reduced_all, sq_qubit2, U2, tapered_val2 = taper_register(all_vecs, alpha_labels, z2_alpha, len(full_targets), "alpha (all 36 targets)")
    rt_err = verify_round_trip(all_vecs, reduced_all, sq_qubit2, U2, tapered_val2, len(full_targets))
    print(f"  round-trip error across all 36 tapered targets: {rt_err:.2e}")
    assert rt_err < 1e-9

    reduced_targets = dict(zip(full_targets.keys(), reduced_all))

    # examine sparsity: how many of the 8 reduced basis states are ever populated?
    used_indices = set()
    for vec in reduced_targets.values():
        nz = np.where(np.abs(vec) > 1e-9)[0]
        used_indices.update(int(i) for i in nz)
    print(f"  reduced (3-qubit, 8-dim) basis states used across all 36 targets: {sorted(used_indices)} "
          f"({len(used_indices)}/8)")

    # -- generic StatePreparation baseline, abstract u3/cx, opt_level=0 (project invariant) --
    BASIS_GATES = ["u3", "cx"]
    gate_counts = {}
    for name, vec in reduced_targets.items():
        qc = QuantumCircuit(3)
        qc.append(StatePreparation(vec), range(3))
        t = transpile(qc, basis_gates=BASIS_GATES, optimization_level=0)
        gate_counts[name] = t.count_ops().get("cx", 0)
        # verify correctness
        sv = np.asarray(Statevector.from_instruction(t))
        idx = int(np.argmax(np.abs(vec)))
        phase = sv[idx] / vec[idx] if abs(vec[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv / phase - vec)))
        assert err < 1e-9, f"{name}: StatePreparation circuit does not match target, err={err:.2e}"

    counts_list = list(gate_counts.values())
    print(f"\n  generic StatePreparation (3-qubit, abstract u3/cx, opt_level=0): "
          f"min={min(counts_list)} max={max(counts_list)} mean={np.mean(counts_list):.2f} CX gates")
    print(f"  (compare to the UNTAPERED 4-qubit fixed ansatz's constant 11 CX gates)")

    constant = len(set(counts_list)) == 1
    print(f"  constant across all 36 targets: {constant} "
          f"({'CDR-compatible' if constant else 'NOT CDR-compatible as-is, same caveat found for TrappedIonOptimizerPlugin'})")

    results = {
        "K": K, "n_targets": len(full_targets),
        "reduced_basis_states_used": sorted(used_indices),
        "n_reduced_basis_states_used": len(used_indices),
        "state_prep_cx_counts": gate_counts,
        "state_prep_cx_min": min(counts_list), "state_prep_cx_max": max(counts_list),
        "state_prep_cx_mean": float(np.mean(counts_list)),
        "constant_gate_count": bool(constant),
        "untapered_baseline_cx": 11,
        "round_trip_verification_err": rt_err,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
