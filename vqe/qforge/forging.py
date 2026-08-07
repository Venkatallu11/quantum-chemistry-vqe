#!/usr/bin/env python3
"""
qforge.forging — entanglement-forging assembly: fragment Hamiltonian,
Schmidt decomposition in the real gauge, beta-register sign-reuse, and
qubit-wise-commuting ("frame='h'") measurement grouping. Extracted from
ef_fragment.py / rank6_symmetry_vd.py / ionq_run.py -- physics and
verification logic unchanged, consolidated into one importable surface.

WHY "frame='h'": every basis change this project uses to read out a
non-Z Pauli in the Z basis is a single-qubit Clifford diagonalization
built from H (for X) and Sdg+H (for Y) -- qiskit's own group_commuting
(qubit_wise=True) grouping, not a general commuting-set + arbitrary-
Clifford-diagonalization scheme (ef_fragment.py's own docstring documents
why that riskier alternative was NOT built). `frame` is exposed as an
explicit parameter on measurement_circuit() rather than hardcoded so a
future diagonalization scheme is a real extension point, not a rewrite --
"h" is the only value implemented and enforced by assertion.

INVARIANTS ENFORCED IN CODE:
  - the identity Pauli is NEVER rescaled: combine_matrices() special-
    cases the identity label to exactly 1.0 regardless of any per-basis/
    per-slot correction, and beta_signs() asserts the beta=sign*alpha
    shortcut actually holds before any caller relies on it.
  - qubit-wise-commuting-violation is asserted impossible in
    combined_basis_label(), not assumed from group_commuting()'s contract.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import ef_fragment as _effrag
from entanglement_forging_h4 import (  # noqa: F401 -- re-exported, unchanged
    precompute_exact_matrices,
    ef_energy_from_matrices,
    ef_energy_from_noisy_matrices,
)

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp

HARTREE_TO_KCAL_MOL = 627.5094740631

# -- fragment Hamiltonian / exact reference / Schmidt decomposition, re-exported --
build_fragment_qop = _effrag.build_fragment_qop
exact_ground_state = _effrag.exact_ground_state
real_gauge = _effrag.real_gauge
schmidt_decompose_real = _effrag.schmidt_decompose_real
decompose_pauli_terms = _effrag.decompose_pauli_terms
group_labels_qubit_wise = _effrag.group_labels_qubit_wise
combined_basis_label = _effrag.combined_basis_label


def beta_signs(u_vecs, v_vecs, K, tol=1e-8):
    """The beta-register sign-reuse shortcut: for a real-gauged state,
    verified (not assumed) that v_n = sign_n * u_n for each of the top K
    Schmidt vectors, so the beta register's matrices never need their own
    independent measurement -- derive_beta_matrices() below reuses the
    alpha measurement directly. Raises if the shortcut does not hold for
    the given decomposition (e.g. a fragment shape where it genuinely
    doesn't), rather than silently using a wrong sign."""
    signs = np.array([1.0 if np.dot(v_vecs[n], u_vecs[n]) >= 0 else -1.0 for n in range(K)])
    residual = max(float(np.max(np.abs(v_vecs[n] - signs[n] * u_vecs[n]))) for n in range(K))
    assert residual < tol, (
        f"beta=sign*alpha shortcut invariant violated (residual={residual:.3e}, tol={tol:.0e}) -- "
        "the beta register must be measured independently for this fragment, not derived"
    )
    return signs, residual


def derive_beta_matrices(alpha_matrices, signs, K):
    """beta_mat = S . alpha_mat . S, S = diag(signs) -- the entire reason
    this project's circuit count is halved vs measuring both registers."""
    S = np.diag(signs)
    return {label: S @ mat @ S for label, mat in alpha_matrices.items()}


def basis_change_h(qc, label):
    """The 'h' measurement frame: H for X, Sdg+H for Y, nothing for Z/I,
    per-qubit, MSQ-first label convention (label[i] acts on qubit
    n-1-i) -- rotates so a Z-basis measurement reads out `label`'s
    eigenbasis."""
    n = len(label)
    for i, ch in enumerate(label):
        qubit = n - 1 - i
        if ch == "X":
            qc.h(qubit)
        elif ch == "Y":
            qc.sdg(qubit)
            qc.h(qubit)


MEASUREMENT_FRAMES = {"h": basis_change_h}


def measurement_circuit(base_circuit, label, frame="h"):
    """base_circuit + a basis change for `label` (from qubit-wise-
    commuting group's combined_basis_label, typically) + measure_all().
    `frame` selects the diagonalization scheme -- only 'h' exists today,
    asserted explicitly rather than silently defaulting to it."""
    assert frame in MEASUREMENT_FRAMES, f"unknown measurement frame {frame!r}, only {list(MEASUREMENT_FRAMES)} implemented"
    qc = base_circuit.copy()
    MEASUREMENT_FRAMES[frame](qc, label)
    qc.measure_all()
    return qc


def pauli_expectation(probs, label):
    """probs: dict[bitstring -> probability]. Standard qiskit bit
    ordering: bitstring[-1] is qubit 0, MSQ-first label convention."""
    n = len(label)
    exp = 0.0
    for bitstring, p in probs.items():
        bits = bitstring[::-1]
        sign = 1
        for q in range(n):
            ch = label[n - 1 - q]
            if ch != "I" and q < len(bits) and bits[q] == "1":
                sign = -sign
        exp += sign * p
    return exp


def combine_matrices(raw, alpha_labels, identity_label, K, per_slot_scale=None, per_label_scale=None):
    """Assemble the K x K alpha matrix for every label from per-slot raw
    measurements, applying an optional CDR correction. INVARIANT: the
    identity label's diagonal is hardcoded to 1.0 -- it is NEVER read
    from `raw` and NEVER divided by a scale, regardless of what
    per_slot_scale/per_label_scale contain for that label."""
    def slot_names(K):
        names = [f"u_{n}" for n in range(K)]
        for n in range(K):
            for m in range(K):
                if n < m:
                    names += [f"(u{n}+u{m})", f"(u{n}-u{m})"]
        return names

    def corrected(name, label):
        v = raw[name][label]
        if per_slot_scale is not None:
            v = v / per_slot_scale.get(name, 1.0)
        if per_label_scale is not None:
            v = v / per_label_scale.get(label, 1.0)
        return v

    mats = {l: np.zeros((K, K)) for l in alpha_labels}
    for n in range(K):
        name = f"u_{n}"
        for l in alpha_labels:
            mats[l][n, n] = 1.0 if l == identity_label else corrected(name, l)
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            name0, name2 = f"(u{n}+u{m})", f"(u{n}-u{m})"
            for l in alpha_labels:
                re = 0.0 if l == identity_label else (corrected(name0, l) - corrected(name2, l)) / 2
                mats[l][n, m] = re
                mats[l][m, n] = re
    return mats


def energy_from_alpha_matrices(alpha_mats, terms, lambdas, enuc, signs, K,
                                exact_energy=None, noiseless_energy=None):
    beta_mats = derive_beta_matrices(alpha_mats, signs, K)
    E = ef_energy_from_noisy_matrices(terms, lambdas, alpha_mats, beta_mats, enuc, K)
    errs = {}
    if exact_energy is not None:
        errs["err_vs_exact_kcal"] = abs(E - exact_energy) * HARTREE_TO_KCAL_MOL
    if noiseless_energy is not None:
        errs["err_vs_noiseless_kcal"] = abs(E - noiseless_energy) * HARTREE_TO_KCAL_MOL
    return E, errs


def slot_names(K):
    names = [f"u_{n}" for n in range(K)]
    for n in range(K):
        for m in range(K):
            if n < m:
                names += [f"(u{n}+u{m})", f"(u{n}-u{m})"]
    return names


def target_states(u_vecs, K):
    """The K diagonal Schmidt vectors + K*(K-1) pair-phase states (2 per
    pair, real gauge) -- the full 36-slot (K=6) target set."""
    targets = {}
    for n in range(K):
        targets[f"u_{n}"] = u_vecs[n]
    for n in range(K):
        for m in range(K):
            if n < m:
                targets[f"(u{n}+u{m})"] = (u_vecs[n] + u_vecs[m]) / np.sqrt(2)
                targets[f"(u{n}-u{m})"] = (u_vecs[n] - u_vecs[m]) / np.sqrt(2)
    assert set(targets.keys()) == set(slot_names(K))
    return targets


def setup_fragment(atoms, nelec, d, K):
    """One-call setup: fragment Hamiltonian -> exact ground state -> real
    gauge -> Schmidt decomposition -> Pauli terms -> beta_signs -> target
    states -- everything needed before angle-fitting and measurement.
    Verifies (not assumes) the Schmidt rank is <=K before returning."""
    qop_bare, qop_pen, enuc = build_fragment_qop(atoms, nelec, d=d)
    n_qubits = qop_bare.num_qubits
    e_elec, psi = exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    psi_real, real_gauge_residual = real_gauge(psi)
    lambdas, u_vecs, v_vecs = schmidt_decompose_real(psi_real, n_qubits)

    max_tail = float(np.max(np.abs(lambdas[K:])))
    assert max_tail < 1e-9, f"Schmidt rank > K={K}: tail={max_tail:.3e} -- K is not exact here"

    terms = decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    assert set(alpha_labels) == set(beta_labels)
    identity_label = "I" * (n_qubits // 2)

    u_top, v_top = u_vecs[:K], v_vecs[:K]
    signs, sign_residual = beta_signs(u_top, v_top, K)

    alpha_cache, beta_cache = precompute_exact_matrices(terms, u_top, v_top)
    noiseless_energy = ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)

    return {
        "d": d, "n_qubits": n_qubits, "terms": terms, "alpha_labels": alpha_labels,
        "identity_label": identity_label, "signs": signs, "lambdas": lambdas, "u_vecs": u_top,
        "enuc": enuc, "exact_energy": exact_energy, "noiseless_energy": noiseless_energy,
        "max_schmidt_tail": max_tail, "real_gauge_residual": real_gauge_residual,
        "targets": target_states(u_top, K), "K": K,
    }
