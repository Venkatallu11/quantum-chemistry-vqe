#!/usr/bin/env python3
"""
task39e_conditioned_correction.py -- iteration 39, Task E. Task D's real
result: naively applying the existing UNCONDITIONAL analytic PEC
correction to ancilla-POSTSELECTED real data made things measurably
WORSE, not better (forte-1: 5.90 -> 14.39 kcal/mol). Diagnosed cause,
now fixed here: postselecting on the ancilla reading 0 is mathematically
equivalent to projecting the register's density matrix onto the
EVEN-WEIGHT subspace (weight in {0,2,4}) and renormalizing, BEFORE
computing any Pauli expectation -- the existing `analytic_A_and_B_
5param` has no such step, so it was computing the ratio for the WRONG
(unconditional) ensemble and applying it to conditioned data.

THE FIX, directly derived, not guessed: P_even = (I + Z(0)Z(1)Z(2)Z(3))/2
is exactly the projector onto weight-even computational-basis states
(eigenvalue of the all-Z parity operator is (-1)^weight, so +1 for even
weight, matching the ancilla's own "reads 0 <=> even weight" convention,
Task 39B/iteration 18's `postselect_counts`). Applying THIS projector to
the SAME dm_A/dm_B density matrices `analytic_A_and_B_5param` already
builds (gate-by-gate, before any basis rotation -- the correct point,
since the ancilla was entangled with the PRE-rotation state) and
renormalizing gives the exact CONDITIONAL ensemble a real ancilla-0
postselection produces. A/B are then computed by tracing Pauli
operators against these CONDITIONED density matrices instead of the raw
ones -- everything else about the function (gate-by-gate propagation,
depolarizing+PEC-inverse per gate) is reused UNCHANGED.

SANITY CHECK, before trusting this for real data: at theta=0 (fully
ideal, no noise at all), the state is EXACTLY in the even-weight sector
already (P(even)=1, Task 39B's own verified finding), so conditioning
must be a mathematical no-op there -- checked directly below, not
assumed.

Run:
    PYTHONHASHSEED=0 python vqe/task39e_conditioned_correction.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

GATE_NAME = "zz"


def even_weight_projector(n=4):
    ZZZZ = np.array([1.0])
    for _ in range(n):
        ZZZZ = np.kron(ZZZZ, np.array([1.0, -1.0]))
    return np.diag((1.0 + ZZZZ) / 2.0)


_P_EVEN = even_weight_projector(4)


def condition_on_even_weight(dm):
    rho = dm.data
    rho_cond = _P_EVEN @ rho @ _P_EVEN
    norm = np.real(np.trace(rho_cond))
    if norm < 1e-12:
        return dm, 0.0  # degenerate, should not happen for any physically reasonable noise level
    return DensityMatrix(rho_cond / norm), float(norm)


def analytic_A_and_B_conditioned(angles, gate_name, p_zz, p_gpi, p_gpi2, delta_zz, delta_gpi2, labels):
    """Strict extension of Task 37C's analytic_A_and_B_5param: identical
    gate-by-gate propagation, but A and B are traced against the
    EVEN-WEIGHT-CONDITIONED density matrices, matching what a real
    ancilla-0 postselection actually measures. Also returns the retained
    (even-weight) probability for A and B, for direct comparison against
    the real observed retained-shot fraction."""
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
    A = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_A_cond.data))) for l in labels}
    B = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_B_cond.data))) for l in labels}
    return A, B, retained_A, retained_B


def _self_test():
    print("\n" + "=" * 96)
    print("  task39e_conditioned_correction.py -- self-test: no-op at theta=0, retained-fraction sanity")
    print("=" * 96)
    from qforge import setup_fragment, fit_all_targets
    from task37c_extended_forward_model import analytic_A_and_B_5param

    K = 6
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)
    angles = fixed_solutions["u_0"]["angles"]
    labels = ["XYYX", "IYYI"]

    print("\n  -- check 1: at theta=0 (zero noise), conditioning must be a no-op --")
    A_unc, B_unc = analytic_A_and_B_5param(angles, GATE_NAME, 0.0, 0.0, 0.0, 0.0, 0.0, labels)
    A_c, B_c, ret_A, ret_B = analytic_A_and_B_conditioned(angles, GATE_NAME, 0.0, 0.0, 0.0, 0.0, 0.0, labels)
    max_diff = max(max(abs(A_unc[l] - A_c[l]) for l in labels), max(abs(B_unc[l] - B_c[l]) for l in labels))
    print(f"    max diff (unconditioned vs conditioned) at theta=0: {max_diff:.3e}   retained_A={ret_A:.6f}  retained_B={ret_B:.6f}")
    print(f"    {'PASS' if max_diff < 1e-9 and abs(ret_A - 1.0) < 1e-9 else 'FAIL -- STOP, do not trust this correction'}")

    print("\n  -- check 2: at real, nonzero calibration (ZZ=0.014593, GPi2=0), retained fraction should be well under 1.0 --")
    A_c2, B_c2, ret_A2, ret_B2 = analytic_A_and_B_conditioned(angles, GATE_NAME, 0.014593, 0.000119, 0.0, 0.0, 0.0, labels)
    print(f"    retained_A={ret_A2:.6f}  retained_B={ret_B2:.6f}  (real ancilla data retained 88.9%/89.8% of shots -- "
          f"should be in a broadly comparable range, not wildly different, if this model is realistic)")


if __name__ == "__main__":
    _self_test()
