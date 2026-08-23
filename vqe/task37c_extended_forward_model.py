#!/usr/bin/env python3
"""
task37c_extended_forward_model.py -- iteration 37, Task C, step 1.
Extends Task 30B's own well-tested `analytic_A_and_B` (2 parameters:
p_zz, p_gpi2, both purely INCOHERENT/depolarizing) to the full 5-parameter
model Task 37B defined: adds coherent per-gate-type angle bias
(delta_zz, delta_gpi2) and readout error (p_readout) -- built as a
strict, verified generalization, NOT a rewrite: with delta_zz=delta_gpi2=
p_readout=0 this reproduces `analytic_A_and_B` bit-for-bit (checked below,
not assumed).

COHERENT BIAS MODEL, verified before use (not hand-derived Euler math
trusted blind): IonQ's own GPi(phi)/GPi2(phi) gates (`qiskit_ionq.
ionq_gates`) are exactly single-qubit rotations by a FIXED angle (pi for
GPi, pi/2 for GPi2) about the axis at angle 2*pi*phi in the XY-plane --
confirmed by direct numerical comparison against the general rotation
formula R(theta,axis) = cos(theta/2)*I - i*sin(theta/2)*(cos(axis)*X +
sin(axis)*Y) for 5 sampled phi values (exact match for GPi2; GPi matches
up to a fixed global phase of i, which is physically irrelevant for
density-matrix evolution -- U*rho*U^dag is invariant under any global
phase on U). A coherent over/under-rotation bias delta is then simply
R(theta+delta, axis) -- the SAME verified formula, evaluated at a
perturbed angle, not a new construction. (`scipy.linalg.
fractional_matrix_power` was tried first and rejected: it hits the
degenerate-eigenvalue branch-cut exactly at theta=pi, giving a WRONG
answer for GPi specifically -- caught by the same before-trusting-it
numerical check, see `_self_test_rotation_formula`.) ZZGate's own bias is
even simpler and needs no manual formula at all: `ZZGate(theta)` IS
`exp(-i*pi*theta*Z⊗Z)` exactly (confirmed by inspection of its matrix,
diag(e^{-i pi theta}, e^{i pi theta}, e^{i pi theta}, e^{-i pi theta})),
so a biased ZZ gate is just `ZZGate(theta + delta_zz)`, calling IonQ's own
gate class directly -- no re-derivation at all.

READOUT MODEL: a symmetric per-qubit bit-flip readout channel with flip
probability p_readout attenuates any Pauli expectation value by
(1-2*p_readout)^w, where w = Hamming weight (number of non-identity
qubits) of the label -- a standard, directly-checkable identity (verified
below via direct density-matrix simulation on random states, not merely
cited). IMPORTANT MODELING CHOICE, made deliberately and disclosed: this
readout factor is applied as a SEPARATE multiplicative un-mixing step to
the REAL RAW measured value, NOT folded into the A/B ratio -- because
folding it into both A and B identically would cancel exactly in B/A and
make p_readout completely unidentifiable through this correction
pathway. Un-mixing the raw data directly (m_raw / (1-2p)^w) is also the
physically correct place for it: readout error is a purely classical,
end-of-circuit effect entirely decoupled from the gate-level PEC algebra
that produces A and B.

Run (self-test -- regression check + both physics checks):
    python vqe/task37c_extended_forward_model.py
"""
import os
import sys
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(__file__))
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from qiskit.quantum_info import DensityMatrix, Operator, Pauli
from qiskit_ionq.ionq_gates import GPIGate, GPI2Gate, ZZGate

_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_I2 = np.eye(2, dtype=complex)


def _rotation(theta, axis_angle):
    """Verified single-qubit rotation by angle theta about the XY-plane
    axis at angle axis_angle -- see module docstring for the check."""
    return np.cos(theta / 2) * _I2 - 1j * np.sin(theta / 2) * (
        np.cos(axis_angle) * _X + np.sin(axis_angle) * _Y)


def biased_gpi_matrix(phi, delta):
    return _rotation(np.pi + delta, 2 * np.pi * phi)


def biased_gpi2_matrix(phi, delta):
    return _rotation(np.pi / 2 + delta, 2 * np.pi * phi)


def biased_zz_matrix(theta, delta):
    return ZZGate(theta + delta).to_matrix()


def hamming_weight(label):
    return sum(1 for c in label if c != "I")


def readout_attenuation(label, p_readout):
    return (1.0 - 2.0 * p_readout) ** hamming_weight(label)


def nearest_bin_p(bins_dict, phi):
    angles = [float(a) for a in bins_dict.keys()]
    idx = int(np.argmin([abs(a - (phi % 1.0)) for a in angles]))
    key = list(bins_dict.keys())[idx]
    return bins_dict[key]["p"]


def analytic_A_and_B_5param(angles, gate_name, p_zz, p_gpi, p_gpi2, delta_zz, delta_gpi2, labels):
    """Strict generalization of task30b_pec_application.analytic_A_and_B:
    p_zz, p_gpi2 are scalar (angle-independent) depolarizing rates, not
    4-bin lookups -- Task 37B's own scope decision, disclosed there (a
    single p_gpi2 scalar rather than 4 separate bin parameters, per its
    'don't overfit' instruction). p_gpi is likewise a single real,
    FIXED-elsewhere scalar (Task 37 Phase 0's real GPi calibration), not
    a free parameter of this function. delta_zz/delta_gpi2 bias the
    UNITARY applied (coherent, gate-TYPE-wide, matches task31h's own
    convention) -- the depolarizing/PEC-inverse steps are computed with
    the SAME scalar p as before, unaffected by delta. Readout is NOT
    applied here -- see module docstring; call `readout_attenuation`
    separately on the raw measured value."""
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
        if op.name == gate_name:  # native 2q gate (zz)
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz))
            p_here = p_zz
            n_here = 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))  # GPi held fixed -- no delta_gpi in scope (Task 37B)
            p_here = p_gpi
            n_here = 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2))
            p_here = p_gpi2
            n_here = 1
        else:
            U = Operator(op.to_matrix())
            dm_A = dm_A.evolve(U, qargs=qargs)
            dm_B = dm_B.evolve(U, qargs=qargs)
            continue
        dm_A = dm_A.evolve(U, qargs=qargs)
        dm_B = dm_B.evolve(U, qargs=qargs)
        dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(p_here, n_here))
        dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(p_here, n_here))
        dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(p_here, n_here))
    A = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_A.data))) for l in labels}
    B = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_B.data))) for l in labels}
    return A, B


def _self_test_regression():
    """delta_zz=delta_gpi2=0 must reproduce task30b's original
    analytic_A_and_B (with its 4-bin gpi/gpi2 lookups collapsed to a
    single flat value, the only structural difference) bit-for-bit."""
    from task30b_pec_application import analytic_A_and_B
    from qforge import setup_fragment, fit_all_targets
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=6, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)
    angles = fixed_solutions["u_0"]["angles"]
    labels = ["XYYX", "IYYI"]
    p_zz_test, p_gpi2_test = 0.0146, 0.003
    flat_gpi_bins = {"0.0": {"p": 0.000119}, "0.25": {"p": 0.000119},
                      "0.5": {"p": 0.000119}, "0.75": {"p": 0.000119}}
    flat_gpi2_bins = {"0.0": {"p": p_gpi2_test}, "0.25": {"p": p_gpi2_test},
                       "0.5": {"p": p_gpi2_test}, "0.75": {"p": p_gpi2_test}}
    A_orig, B_orig = analytic_A_and_B(angles, "zz", p_zz_test, flat_gpi_bins, flat_gpi2_bins, labels)
    A_new, B_new = analytic_A_and_B_5param(angles, "zz", p_zz_test, 0.000119, p_gpi2_test, 0.0, 0.0, labels)
    max_diff = max(max(abs(A_orig[l] - A_new[l]) for l in labels),
                    max(abs(B_orig[l] - B_new[l]) for l in labels))
    print(f"  regression vs task30b analytic_A_and_B (delta=0): max diff = {max_diff:.3e}  "
          f"{'PASS' if max_diff < 1e-9 else 'FAIL'}")
    return max_diff < 1e-9


def _self_test_rotation_formula():
    """Confirms the biased-gate construction against fractional_matrix_power
    where that method is reliable (GPi2, theta=pi/2 -- no branch-cut issue)
    and documents/demonstrates why it was REJECTED for GPi (theta=pi,
    degenerate eigenvalues under the principal log branch)."""
    from scipy.linalg import fractional_matrix_power
    phi, delta = 0.31, 0.05
    U2 = GPI2Gate(phi).to_matrix()
    frac2 = fractional_matrix_power(U2, (np.pi / 2 + delta) / (np.pi / 2))
    direct2 = biased_gpi2_matrix(phi, delta)
    diff2 = np.max(np.abs(frac2 - direct2))
    print(f"  GPi2 biased-rotation cross-check (fractional_matrix_power, reliable here): "
          f"diff={diff2:.3e}  {'PASS' if diff2 < 1e-9 else 'FAIL'}")

    U1 = GPIGate(phi).to_matrix()
    frac1 = fractional_matrix_power(U1, (np.pi + delta) / np.pi)
    direct1 = biased_gpi_matrix(phi, delta)
    diff1_raw = np.max(np.abs(frac1 - direct1))
    print(f"  GPi fractional_matrix_power vs direct rotation (branch-cut case, EXPECTED to diverge): "
          f"diff={diff1_raw:.3e}  (this is WHY the direct rotation formula is used instead)")

    U0 = biased_gpi_matrix(phi, 0.0)
    diff0 = np.max(np.abs(GPIGate(phi).to_matrix() - 1j * U0))
    print(f"  GPi direct-rotation formula vs real gate matrix (up to global phase i, delta=0): "
          f"diff={diff0:.3e}  {'PASS' if diff0 < 1e-9 else 'FAIL'}")
    return diff2 < 1e-9 and diff0 < 1e-9


def _self_test_readout_formula():
    """Checks (1-2p)^w the way it actually applies to a real measurement:
    readout error acts on the CLASSICAL BITS after the Pauli-diagonalizing
    basis rotation (exactly Task 37A's own validated `method3` pipeline --
    reused here directly, not re-derived), not on the density matrix's
    original-basis populations (an earlier version of this test wrongly
    mixed the un-rotated diagonal, which is only correct for a
    Z-only/identity-only label -- caught by this very check failing on
    the X/Y-containing labels below before this file was trusted)."""
    from task2_fold_response_dataset import native_basis_change
    import ef_fragment as effrag_mod
    from qiskit.quantum_info import Statevector

    rng = np.random.default_rng(0)
    ok = True
    n = 2
    for trial in range(20):
        psi_vec = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
        psi_vec /= np.linalg.norm(psi_vec)
        from qiskit import QuantumCircuit
        base = QuantumCircuit(n)
        base.initialize(psi_vec, range(n))
        p_ro = rng.uniform(0, 0.02)
        for label in ["XZ", "ZX", "ZZ", "XX", "IZ", "ZI"]:
            ideal_exp = float(np.real(Statevector.from_instruction(
                QuantumCircuit(n).compose(base)).expectation_value(Pauli(label))))
            basis_qc = native_basis_change(effrag_mod.combined_basis_label([label]), "zz")
            full = base.compose(basis_qc)
            sv = Statevector.from_instruction(full)
            probs = sv.probabilities_dict()
            non_id_positions = [i for i, c in enumerate(reversed(label)) if c != "I"]
            w = hamming_weight(label)
            flip1 = np.array([[1 - p_ro, p_ro], [p_ro, 1 - p_ro]])
            M = flip1
            for _ in range(n - 1):
                M = np.kron(M, flip1)
            dense = np.zeros(2**n)
            for bitstring, prob in probs.items():
                dense[int(bitstring.zfill(n), 2)] = prob
            noisy = M @ dense
            total = 0.0
            for idx, prob in enumerate(noisy):
                bits = format(idx, f"0{n}b")
                parity = 1
                for pos in non_id_positions:
                    if bits[n - 1 - pos] == "1":
                        parity *= -1
                total += parity * prob
            rhs = ideal_exp * (1 - 2 * p_ro) ** w
            if abs(total - rhs) > 1e-9:
                ok = False
    print(f"  readout attenuation formula (1-2p)^w vs rotated-frame bit-flip simulation, 20x6 random cases: "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def _self_test():
    print("\n" + "=" * 96)
    print("  task37c_extended_forward_model.py -- self-test (regression + both physics checks)")
    print("=" * 96)
    ok1 = _self_test_regression()
    ok2 = _self_test_rotation_formula()
    ok3 = _self_test_readout_formula()
    print(f"\n  OVERALL: {'ALL PASS' if (ok1 and ok2 and ok3) else 'AT LEAST ONE FAILED -- DO NOT USE'}")


if __name__ == "__main__":
    _self_test()
