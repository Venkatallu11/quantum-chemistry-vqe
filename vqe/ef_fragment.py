#!/usr/bin/env python3
"""
ef_fragment.py — Entanglement Forging core, generalized to ARBITRARY
fragments (not just the hardcoded H4 case in entanglement_forging_h4.py).
============================================================================
entanglement_forging_h4.py hardcodes: atoms [0,1,2,3], 8 qubits split into
two fixed 4-qubit registers, and a 4-phase-circuit (k=0,1,2,3) cross-term
reconstruction. This module generalizes the THREE fragment-specific pieces
so the same machinery works for molecular-tailoring fragments of any size
(the covalent_fragment.py H6 layout needs an 8-qubit fragment A [0,1,2,3],
an 8-qubit fragment B [2,3,4,5], AND a 4-qubit overlap [2,3] -- three
different fragment shapes in one run):

  1. build_fragment_qop(atoms, nelec, d) -- generalizes build_h4_qop():
     any real chain-position atom list (via covalent_fragment.hchain, so a
     fragment cut from the middle of a chain sits at its true bond
     distance, not shifted to the origin), any electron count.

  2. schmidt_decompose_real(psi, n_qubits) -- generalizes schmidt_decompose()
     to any even qubit count, AND adds a real-gauge rotation: eigsh returns
     the ground state up to an arbitrary GLOBAL PHASE. For a real physical
     Hamiltonian (real one/two-electron integrals -> real Jordan-Wigner
     matrix, confirmed for this codebase's chem.py/qiskit-nature pipeline),
     the true ground state is real up to that phase, so rotating the
     largest-magnitude amplitude onto the real axis makes the WHOLE vector
     real. This is checked, not assumed -- rotate, then assert the residual
     imaginary part is negligible, per-fragment, every time.

  3. decompose_pauli_terms(qop_bare, n_qubits) -- generalizes the hardcoded
     label[:4]/label[4:] split to label[:half]/label[half:] for any even
     n_qubits.

  Everything else (precompute_exact_matrices, ef_energy_from_matrices,
  ef_energy_from_noisy_matrices) is already generic in
  entanglement_forging_h4.py -- imported and reused unchanged.

WHY THE REAL GAUGE MATTERS FOR CIRCUIT COUNT: once u_n, u_m (Schmidt
vectors) are real, every Pauli label P is either a REAL matrix (even
number of Y's in the label) or a PURELY IMAGINARY matrix (odd number of
Y's) -- a standard fact about tensor products of Pauli operators, since
Y = i * (real antisymmetric matrix) and I/X/Z are real. For real vectors
u_n, u_m, the cross term x = u_n^T P u_m is then PROVABLY either purely
real (even-Y label) or purely imaginary (odd-Y label) -- never a general
complex number with both parts nonzero. Concretely, with
  E_k = (1/2)[D + D' + i^-k x + i^k x*],  D=u_n^T P u_n, D'=u_m^T P u_m
  (E0 - E2)/2 = Re(x)   (always, from k=0,2)
  (E3 - E1)/2 = Im(x)   (always, from k=1,3)
so an even-Y label's nonzero part (Re(x)) is recoverable from ONLY the
k=0,2 phase circuits, and an odd-Y label's nonzero part (Im(x)) from ONLY
k=1,3 -- each label needs exactly 2 of the 4 phase circuits, not all 4.
This is verified locally (numpy-exact) against the original all-4-circuit
reconstruction before it's trusted for any real IonQ submission.
"""
import os
import sys
import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(__file__))
import chem
from covalent_fragment import hchain
from molecules_real import _constrain_particle_number
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit.quantum_info import SparsePauliOp

# Re-exported so callers only need to import this module for the generic
# (already fragment-size-agnostic) machinery.
from entanglement_forging_h4 import (  # noqa: F401
    precompute_exact_matrices,
    ef_energy_from_matrices,
    ef_energy_from_noisy_matrices,
)

REAL_GAUGE_TOLERANCE = 1e-10


def build_fragment_qop(atoms, nelec, d=1.0):
    """Generalizes entanglement_forging_h4.build_h4_qop to any real
    chain-position atom list. `atoms` are REAL chain indices (e.g. [2,3,4,5]
    for a fragment cut from the middle of a longer chain) -- hchain()
    places them at their true bond distance, not re-centered at the origin,
    matching covalent_fragment.py's convention."""
    geom = hchain(atoms, d)
    S, T, V, eri, enuc = chem.integrals(geom)
    _ehf, C, Hc = chem.rhf(S, T, V, eri, enuc, nelec=nelec)
    h1 = C.T @ Hc @ C
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)
    ee = ElectronicEnergy.from_raw_integrals(h1, h2)
    qop_bare = JordanWignerMapper().map(ee.second_q_op())
    qop_penalized = _constrain_particle_number(qop_bare, nelec)
    return qop_bare, qop_penalized, enuc


def exact_ground_state(qop_penalized):
    """Full exact diagonalization (penalized only to select the correct
    electron-count sector -- penalty is exactly 0 on the true eigenstate)."""
    M = qop_penalized.to_matrix(sparse=True)
    val, vec = eigsh(M, k=1, which="SA")
    return float(val[0]), vec[:, 0]


def real_gauge(psi):
    """Rotate away eigsh's arbitrary global phase so the ground state is
    real. Raises if the residual imaginary part isn't negligible -- for
    this codebase's real-integral Hamiltonians the true ground state IS
    real up to phase, so a large residual means something is actually
    wrong (e.g. a genuinely complex Hamiltonian, or degeneracy), not
    something to silently paper over."""
    idx = int(np.argmax(np.abs(psi)))
    phase = np.angle(psi[idx])
    psi_rotated = psi * np.exp(-1j * phase)
    residual_imag = float(np.max(np.abs(psi_rotated.imag)))
    assert residual_imag < REAL_GAUGE_TOLERANCE, (
        f"State is not real after phase rotation (residual imag={residual_imag:.3e}, "
        f"tolerance={REAL_GAUGE_TOLERANCE:.0e}) -- real-gauge assumption does not hold "
        "for this Hamiltonian; do not use the 2-phase-circuit shortcut."
    )
    return psi_rotated.real, residual_imag


def schmidt_decompose_real(psi_real, n_qubits):
    """Generalizes entanglement_forging_h4.schmidt_decompose to any even
    n_qubits, for an already-real-gauged state vector. Same bipartite
    convention: qubits 0..half-1 = alpha (columns), qubits half..n-1 = beta
    (rows) -- reshape(dim_half, dim_half) gives rows=beta, cols=alpha,
    matching the original module's documented index convention."""
    assert n_qubits % 2 == 0, "fragment must have an even qubit count to split in half"
    assert np.max(np.abs(psi_real.imag if np.iscomplexobj(psi_real) else 0)) < REAL_GAUGE_TOLERANCE, \
        "schmidt_decompose_real expects an already real-gauged (real-valued) state"
    half = n_qubits // 2
    dim_half = 2 ** half
    mat = np.asarray(psi_real, dtype=float).reshape(dim_half, dim_half)
    U, S, Vh = np.linalg.svd(mat)
    lambdas = S
    u_vecs = Vh    # alpha Schmidt vectors, real
    v_vecs = U.T   # beta Schmidt vectors, real
    return lambdas, u_vecs, v_vecs


def decompose_pauli_terms(qop_bare, n_qubits):
    """Generalizes entanglement_forging_h4.decompose_pauli_terms's fixed
    label[:4]/label[4:] split to label[:half]/label[half:] for any even
    n_qubits. Same qiskit label convention: label[:half] acts on the
    higher-indexed (beta) qubits, label[half:] on the lower-indexed
    (alpha) qubits."""
    assert n_qubits % 2 == 0
    half = n_qubits // 2
    terms = []
    for label, coeff in qop_bare.to_list():
        beta_label = label[:half]
        alpha_label = label[half:]
        terms.append((alpha_label, beta_label, complex(coeff)))
    return terms


def label_y_count(label):
    return label.count("Y")


def label_is_real(label):
    """True if this Pauli label's matrix representation is real-valued
    (even number of Y's), False if purely imaginary (odd number of Y's).
    See module docstring for the derivation."""
    return label_y_count(label) % 2 == 0


# ---------------------------------------------------------------------------
# Qubit-wise commuting measurement grouping -- cuts circuit count by grouping
# labels that can be read out from ONE shared measurement circuit. Verified
# (see scratchpad verify_grouping.py during development) against exact
# Statevector math: max error 3e-14 across 222 sampled (pair, group-member,
# phase) combinations on H4. Uses qiskit's own group_commuting(qubit_wise=True)
# -- the same grouping BackendEstimatorV2 uses internally -- NOT the riskier
# general-commuting + Clifford-diagonalization grouping (qiskit has no
# built-in "find a diagonalizing Clifford for a general commuting set"
# utility; implementing that from scratch was judged too high-risk for a
# further ~2x reduction on top of the real-gauge (2x) and beta-reuse (2x)
# savings already banked).
# ---------------------------------------------------------------------------

def group_labels_qubit_wise(labels):
    """Group Pauli labels into qubit-wise commuting cliques. Returns a list
    of lists of member labels."""
    op = SparsePauliOp.from_list([(l, 1.0) for l in labels])
    groups = op.group_commuting(qubit_wise=True)
    return [list(g.paulis.to_labels()) for g in groups]


def combined_basis_label(group_labels):
    """The single per-qubit measurement basis that covers every member of a
    qubit-wise commuting group: for each qubit, whichever non-'I' character
    appears in any member (qubit-wise commuting guarantees no two members
    disagree on a qubit's axis), else 'I'. A basis_change() + measure_all()
    on this combined label yields probabilities from which EVERY member's
    own pauli_expectation() can be read directly (pauli_expectation already
    ignores 'I' positions in a label, so a member with 'I' at a qubit the
    combined basis rotated for a DIFFERENT member is unaffected)."""
    n = len(group_labels[0])
    combined = []
    for q in range(n):
        chars = {lab[q] for lab in group_labels if lab[q] != "I"}
        assert len(chars) <= 1, (
            f"qubit-wise commuting violation at position {q}: {chars} -- "
            "this should be impossible for a group_commuting(qubit_wise=True) group"
        )
        combined.append(chars.pop() if chars else "I")
    return "".join(combined)
