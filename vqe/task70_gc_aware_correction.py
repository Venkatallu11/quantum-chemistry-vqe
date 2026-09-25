#!/usr/bin/env python3
"""
task70_gc_aware_correction.py -- iteration 70. Fixes a real, confirmed
methodological gap caught by external review of the 56-circuit hardware
plan: `task39e_conditioned_correction.analytic_A_and_B_conditioned`
(used by task59/task60/task69's own correction step) propagates the
assumed noise model ONLY through the state-prep gates, then conditions
on even weight and reads Pauli(l) directly off the conditioned density
matrix -- implicitly treating the MEASUREMENT circuit (ancilla-parity +
basis-change/diagonalizer) as noiseless. task47_gc_noisy_champion_
comparison.py's own docstring already flags this as exact for QWC's
basis change (0 extra gates) but WRONG for GC1/GC2/GC3 (3, 3, 5 real
extra native 2q gates respectively) -- confirmed by re-reading that
file's own derivation, not just asserted here.

task47 only ever built a FORWARD SIMULATION of this (comparing a
hypothetical true-vs-assumed noise draw), never turned it into an A/B
RATIO correction usable on real measured data. This does that: same
gate-by-gate state-prep propagation as analytic_A_and_B_conditioned,
but continues propagating the SAME assumed-noise model through the
group's own native diagonalizer circuit (applied AFTER even-weight
conditioning, matching task47's own physically-derived ordering:
conditioning encodes the ancilla-parity CNOTs' effect exactly, so the
diagonalizer -- which happens after those CNOTs in the real circuit --
must be applied to the CONDITIONED density matrix, not the raw one).

REGRESSION CHECK, must pass before trusting this for anything: on GC0
(the trivial group, 0 extra native gates), this must reproduce
analytic_A_and_B_conditioned's own result EXACTLY (a 0-gate diagonalizer
changes nothing physically) -- the same check task47 already ran and
passed for its own forward-simulation version, independently repeated
here for this A/B-ratio version.

Run:
    PYTHONHASHSEED=0 python vqe/task70_gc_aware_correction.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task39e_conditioned_correction import condition_on_even_weight, analytic_A_and_B_conditioned
from qiskit.quantum_info import DensityMatrix, Operator, Pauli


def analytic_A_and_B_conditioned_gc(angles, gate_name, p_zz, p_gpi, p_gpi2, delta_zz, delta_gpi2,
                                      diag_native_circuit, transformed_items):
    """GC-aware extension of analytic_A_and_B_conditioned: continues the
    SAME assumed-noise gate-by-gate propagation through the group's own
    native diagonalizer, applied to the even-weight-CONDITIONED density
    matrix (matching task47's physically-derived ordering), then reads
    off each original label via its diagonal Z-string representative and
    sign (qiskit handles internal qubit indexing, matching task47's own
    approach, sidestepping any manual string-indexing bug)."""
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native

    qc = to_native(build_ansatz(angles), gate_name)
    n = qc.num_qubits
    dm_A = DensityMatrix.from_label("0" * n)
    dm_B = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == gate_name:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz))
            p_here, n_here = p_zz, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))
            p_here, n_here = p_gpi, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2))
            p_here, n_here = p_gpi2, 1
        else:
            dm_A = dm_A.evolve(Operator(op.to_matrix()), qargs=qargs)
            dm_B = dm_B.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm_A = dm_A.evolve(U, qargs=qargs)
        dm_B = dm_B.evolve(U, qargs=qargs)
        dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(p_here, n_here))
        dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(p_here, n_here))
        dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(p_here, n_here))

    dm_A_cond, retained_A = condition_on_even_weight(dm_A)
    dm_B_cond, retained_B = condition_on_even_weight(dm_B)

    # NEW: continue propagation through the group's own diagonalizer,
    # applied to the CONDITIONED density matrices.
    for instr in diag_native_circuit.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [diag_native_circuit.find_bit(q).index for q in instr.qubits]
        if op.name == gate_name:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz))
            p_here, n_here = p_zz, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))
            p_here, n_here = p_gpi, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2))
            p_here, n_here = p_gpi2, 1
        else:
            dm_A_cond = dm_A_cond.evolve(Operator(op.to_matrix()), qargs=qargs)
            dm_B_cond = dm_B_cond.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm_A_cond = dm_A_cond.evolve(U, qargs=qargs)
        dm_B_cond = dm_B_cond.evolve(U, qargs=qargs)
        dm_A_cond = apply_pauli_mixture(dm_A_cond, qargs, depolarizing_weights(p_here, n_here))
        dm_B_cond = apply_pauli_mixture(dm_B_cond, qargs, depolarizing_weights(p_here, n_here))
        dm_B_cond = apply_pauli_mixture(dm_B_cond, qargs, pec_inverse_weights(p_here, n_here))

    A, B = {}, {}
    for item in transformed_items:
        Z_op = np.asarray(Pauli(item.diagonal).to_matrix())
        A[item.original] = item.sign * float(np.real(np.trace(Z_op @ dm_A_cond.data)))
        B[item.original] = item.sign * float(np.real(np.trace(Z_op @ dm_B_cond.data)))
    return A, B, retained_A, retained_B


def _regression_check():
    print("\n" + "=" * 96)
    print("  task70_gc_aware_correction.py -- regression check: GC0 (0-gate) must exactly match the old correction")
    print("=" * 96)
    from qforge import setup_fragment, fit_all_targets
    from general_commuting_measurements import build_general_commuting_measurement_plan
    from native_stateprep import to_native
    from task37b_h4_noise_model import GPI_REAL_MEAN

    K = 6
    GATE_NAME = "zz"
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)

    trivial_idx = [i for i, d in enumerate(diagonalizers) if d.two_qubit_count == 0 and d.one_qubit_count == 0][0]
    diag_native = to_native(diagonalizers[trivial_idx].to_circuit(), GATE_NAME)
    print(f"  trivial group index={trivial_idx}, native gate count={dict(diag_native.count_ops())} (should be empty)")

    angles = fixed_solutions["u_0"]["angles"]
    p_zz, p_gpi2 = 0.014593, 0.0006
    A_old, B_old, retA_old, retB_old = analytic_A_and_B_conditioned(
        angles, GATE_NAME, p_zz, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0, groups[trivial_idx])
    A_new, B_new, retA_new, retB_new = analytic_A_and_B_conditioned_gc(
        angles, GATE_NAME, p_zz, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0, diag_native, diagonalizers[trivial_idx].transformed)

    max_diff = max(max(abs(A_old[l] - A_new[l]) for l in groups[trivial_idx]),
                    max(abs(B_old[l] - B_new[l]) for l in groups[trivial_idx]))
    print(f"  max diff (old vs new, GC0 trivial group): {max_diff:.3e}   "
          f"retained: old=({retA_old:.6f},{retB_old:.6f}) new=({retA_new:.6f},{retB_new:.6f})")
    print(f"  {'PASS' if max_diff < 1e-9 else 'FAIL -- STOP, do not trust the GC-aware correction'}")
    if max_diff >= 1e-9:
        raise RuntimeError("GC-aware correction does not reduce to the old correction on GC0 -- bug, stop")

    print(f"\n  -- sanity: hard group (index 3) SHOULD now differ from the old (noiseless-diagonalizer) correction --")
    hard_idx = max(range(len(diagonalizers)), key=lambda i: diagonalizers[i].two_qubit_count)
    diag_native_hard = to_native(diagonalizers[hard_idx].to_circuit(), GATE_NAME)
    A_old_h, B_old_h, _, _ = analytic_A_and_B_conditioned(
        angles, GATE_NAME, p_zz, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0, groups[hard_idx])
    A_new_h, B_new_h, retA_h, retB_h = analytic_A_and_B_conditioned_gc(
        angles, GATE_NAME, p_zz, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0, diag_native_hard, diagonalizers[hard_idx].transformed)
    max_diff_h = max(max(abs(A_old_h[l] - A_new_h[l]) for l in groups[hard_idx]),
                       max(abs(B_old_h[l] - B_new_h[l]) for l in groups[hard_idx]))
    print(f"  hard group (index {hard_idx}, 5 extra 2q gates): max diff (old vs new) = {max_diff_h:.4f}  "
          f"(should be clearly nonzero -- old model was ignoring real diagonalizer noise here)")
    print(f"  retained fraction: new A={retA_h:.6f} B={retB_h:.6f} (should be lower than the state-prep-only "
          f"retained fraction, since the diagonalizer's own gates add more depolarizing loss)")


if __name__ == "__main__":
    _regression_check()
