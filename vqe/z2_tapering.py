#!/usr/bin/env python3
"""
z2_tapering.py — Task A: Z2 symmetry tapering (Bravyi, Gambetta,
Mezzacapo, Temme, arXiv:1701.08213), applied to entanglement forging's
per-register (4-qubit alpha/beta) problem.
============================================================================
VERIFIED (not assumed) via qiskit's own Z2Symmetries.find_z2_symmetries,
applied to a proxy operator built from just the alpha-register's own 37
unique Pauli labels (symmetry-finding depends only on which Pauli
strings appear, not their coefficients, so this correctly finds
symmetries of the ALPHA REGISTER'S measurement structure specifically,
not the full 8-qubit joint Hamiltonian):

  alpha register: exactly ONE Z2 symmetry, generator = ZZZZ (product of Z
    over all 4 alpha qubits) -- this is alpha-electron-number parity.
  beta register: the SAME, ZZZZ, independently confirmed.
  full 8-qubit Hamiltonian: THREE independent generators (ZZZZIIII,
    ZIZIIZIZ, ZIZIZIZI) -- the first is exactly the alpha-register ZZZZ
    embedded with identity on beta; the other two have support spanning
    BOTH registers (a genuine point-group-type symmetry of the symmetric
    H4 chain, real but NOT block-local). Confirmed by direct decoding of
    their Pauli strings against this project's qubit ordering (0-3=alpha,
    4-7=beta), not assumed from the group's generator count alone.

THIS FILE exploits ONLY the two independent, block-local, per-register
symmetries (alpha-ZZZZ, beta-ZZZZ) -- each taperable INDEPENDENTLY within
its own register's own circuit, exactly matching entanglement forging's
existing separated-register measurement scheme. The cross-register
point-group symmetries are real (verified above) but NOT exploited here:
tapering them would correlate the removed qubit's value ACROSS the
alpha/beta split forging deliberately keeps independent, which is a
real, unresolved reformulation question, not attempted here -- see
ALTERNATIVES NOT TAKEN at the end of this run.

WHY THIS SHOULD BE "FREE" ACCURACY, stated as a testable prediction
before verifying it: every physical state in this problem (all 6 K=6
Schmidt vectors) lives entirely in the alpha register's weight-2 sector
(2 particles among 4 spin-orbitals) -- and weight-2 has HAMMING WEIGHT
2, which is EVEN, so the ZZZZ parity eigenvalue ((-1)^weight) is +1 for
EVERY such state, uniformly. If true, tapering costs nothing (no
information is thrown away, since every physical state already lives in
the SAME sector this symmetry would select) and removes one qubit for
free. Verified directly below, not assumed.

Run:
    python vqe/z2_tapering.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
import ef_fragment as effrag
from qiskit.quantum_info import SparsePauliOp, Pauli, Z2Symmetries, Statevector, Operator

HARTREE_TO_KCAL_MOL = 627.5094740631
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "z2_tapering_results.json")


def find_register_symmetry(labels, register_name):
    op = SparsePauliOp.from_list([(l, 1.0) for l in labels])
    z2 = Z2Symmetries.find_z2_symmetries(op)
    print(f"  {register_name} register: {len(z2.symmetries)} symmetry generator(s): "
          f"{[str(s) for s in z2.symmetries]}")
    return z2


def verify_full_hamiltonian_symmetries(qop_bare):
    full_op = SparsePauliOp.from_list(qop_bare.to_list())
    z2full = Z2Symmetries.find_z2_symmetries(full_op)
    print(f"  full 8-qubit Hamiltonian: {len(z2full.symmetries)} generator(s): "
          f"{[str(s) for s in z2full.symmetries]}")
    # verify each commutes with H directly (not just trust the library)
    H_matrix = qop_bare.to_matrix(sparse=False)
    worst_comm = 0.0
    for s in z2full.symmetries:
        Pm = np.asarray(Pauli(s).to_matrix())
        comm = np.max(np.abs(Pm @ H_matrix - H_matrix @ Pm))
        worst_comm = max(worst_comm, comm)
    print(f"  verified: max ||[P,H]|| over all found generators = {worst_comm:.2e} "
          f"({'commutes exactly' if worst_comm < 1e-8 else 'DOES NOT COMMUTE -- library result is wrong'})")
    assert worst_comm < 1e-8, "a found symmetry does not actually commute with H -- refusing to proceed"

    # decode block-locality of each generator against this project's 0-3=alpha,4-7=beta convention
    for s in z2full.symmetries:
        label = str(s)
        n = len(label)
        alpha_support = any(label[n - 1 - q] != "I" for q in range(4))
        beta_support = any(label[n - 1 - q] != "I" for q in range(4, 8))
        locality = "alpha-only" if alpha_support and not beta_support else \
                   "beta-only" if beta_support and not alpha_support else "CROSS-register"
        print(f"    {label}: {locality}")
    return z2full


def taper_register(u_vecs, alpha_labels, z2, K, register_name):
    """Apply the register's own Z2 symmetry Clifford to each Schmidt
    vector, verify the tapered qubit's value is UNIFORM and matches the
    predicted +1 sector (weight-2 parity), taper, and verify the
    round-trip reconstructs the original vectors exactly."""
    assert len(z2.symmetries) == 1, f"expected exactly 1 symmetry for {register_name}, got {len(z2.symmetries)}"
    clifford_op = z2.cliffords[0]
    U_clifford = np.asarray(Operator(clifford_op).data)
    sq_qubit = z2.sq_list[0]
    n_qubits = int(round(np.log2(u_vecs.shape[1])))

    # IMPORTANT, verified not assumed: qiskit's Z2Symmetries maps the symmetry generator to a
    # single-qubit Pauli that is NOT necessarily Z -- checked directly here (this project's H4
    # case gives sq_paulis=[IIIX], an X, confirmed by applying convert_clifford to the symmetry
    # operator itself and inspecting the result before trusting this). The tapered qubit is fixed
    # in THAT Pauli's eigenbasis after the Clifford, not automatically the Z (computational) basis
    # -- an extra single-qubit basis-change gate is needed before the qubit's value can be read
    # off as a simple bit, matching the standard tapering recipe (Bravyi et al. 2017) but easy to
    # get wrong by assuming Z without checking.
    sq_pauli_label = str(z2.sq_paulis[0])
    sq_char = sq_pauli_label[n_qubits - 1 - sq_qubit]
    print(f"  {register_name}: Clifford maps symmetry to qubit {sq_qubit}, fixed in the {sq_char} basis")
    H2 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    Sdg2 = np.array([[1, 0], [0, -1j]])
    Y_TO_Z = H2 @ Sdg2  # verified directly (not hand-derived from memory): maps |+i>->|0>, |-i>->|1>
    basis_gate = {"X": H2, "Y": Y_TO_Z, "Z": np.eye(2)}[sq_char]
    ops = [np.eye(2)] * n_qubits
    ops[n_qubits - 1 - sq_qubit] = basis_gate
    U_basis = ops[0]
    for m in ops[1:]:
        U_basis = np.kron(U_basis, m)
    U = U_basis @ U_clifford

    n_qubits = int(round(np.log2(u_vecs.shape[1])))
    tapered_qubit_values = []
    reduced_vecs = []
    for n in range(K):
        vec = u_vecs[n]
        rotated = U @ vec
        # extract the tapered qubit's basis-index bit and check it's uniform across all basis
        # components with nonzero amplitude
        nz = np.where(np.abs(rotated) > 1e-9)[0]
        bit_vals = set((i >> sq_qubit) & 1 for i in nz)
        assert len(bit_vals) == 1, f"{register_name} u_{n}: tapered qubit not fixed, bits={bit_vals}"
        bit = bit_vals.pop()
        tapered_qubit_values.append(bit)
        # extract the reduced (n_qubits-1)-dim vector: keep only components with this qubit's bit
        # value, remap remaining qubit indices by dropping bit `sq_qubit`
        dim_reduced = 2 ** (n_qubits - 1)
        reduced = np.zeros(dim_reduced, dtype=complex)
        for i in nz:
            low = i & ((1 << sq_qubit) - 1)
            high = i >> (sq_qubit + 1)
            new_idx = (high << sq_qubit) | low
            reduced[new_idx] = rotated[i]
        norm = np.linalg.norm(reduced)
        reduced_vecs.append(reduced / norm)
    assert len(set(tapered_qubit_values)) == 1, \
        f"{register_name}: tapered qubit value NOT uniform across Schmidt vectors: {tapered_qubit_values}"
    print(f"  {register_name}: tapered qubit value UNIFORM across all {K} Schmidt vectors: "
          f"{tapered_qubit_values[0]} (predicted: even-weight physical states -> ZZZZ=+1)")

    return np.array(reduced_vecs), sq_qubit, U, tapered_qubit_values[0]


def verify_round_trip(u_vecs, reduced_vecs, sq_qubit, U, tapered_value, K):
    n_qubits = int(round(np.log2(u_vecs.shape[1])))
    worst = 0.0
    for n in range(K):
        reduced = reduced_vecs[n]
        dim_full = 2 ** n_qubits
        rotated = np.zeros(dim_full, dtype=complex)
        for new_idx, amp in enumerate(reduced):
            low = new_idx & ((1 << sq_qubit) - 1)
            high = new_idx >> sq_qubit
            i = (high << (sq_qubit + 1)) | (tapered_value << sq_qubit) | low
            rotated[i] = amp
        reconstructed = U.conj().T @ rotated
        idx = int(np.argmax(np.abs(u_vecs[n])))
        phase = reconstructed[idx] / u_vecs[n][idx] if abs(u_vecs[n][idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(reconstructed / phase - u_vecs[n])))
        worst = max(worst, err)
    return worst


def taper_pauli_matrix(label, U, sq_qubit, tapered_value, n_qubits):
    """Taper a single Pauli label by DIRECT MATRIX conjugation with the
    SAME combined transform U (Clifford + basis-change) used for the
    state vectors -- avoids any separate symbolic Pauli-string tracking,
    which would need to independently re-derive the same U_basis
    conjugation applied to the states (a real risk of the two ending up
    inconsistent, exactly the kind of mismatch this project's own
    established practice exists to catch). Verifies directly that the
    rotated operator is block-diagonal in the tapered qubit's Z-basis
    (i.e. proportional to I or Z there, never X/Y) before trusting the
    reduction -- required for a consistent taper, checked not assumed."""
    P = np.asarray(Pauli(label).to_matrix())
    P_rot = U @ P @ U.conj().T
    dim = 2 ** n_qubits
    idx_for_bit = {0: [], 1: []}
    for i in range(dim):
        bit = (i >> sq_qubit) & 1
        idx_for_bit[bit].append(i)
    # verify block-diagonal (no coupling between the two tapered-qubit sectors)
    off_block = P_rot[np.ix_(idx_for_bit[0], idx_for_bit[1])]
    assert np.max(np.abs(off_block)) < 1e-9, \
        f"label {label}: tapered qubit is NOT block-diagonal after rotation (off-block max={np.max(np.abs(off_block)):.2e}) -- not a consistent taper"
    block0 = P_rot[np.ix_(idx_for_bit[0], idx_for_bit[0])]
    block1 = P_rot[np.ix_(idx_for_bit[1], idx_for_bit[1])]
    # the surviving (tapered_value) block is the reduced operator; verify the reduced operator
    # ANYWAY only needs to match on the surviving sector, but also check the general case
    reduced_matrix = block0 if tapered_value == 0 else block1
    return reduced_matrix


def main():
    print("\n" + "=" * 96)
    print("  z2_tapering.py -- Task A: Z2 symmetry tapering")
    print("=" * 96)

    qop_bare, qop_pen, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    K = int(np.sum(lambdas > 1e-9))
    print(f"\n  exact_energy = {exact_energy:.6f} Ha, K (Schmidt rank) = {K}")

    print(f"\n  -- finding Z2 symmetries (per-register and full) --")
    z2_alpha = find_register_symmetry(alpha_labels, "alpha")
    z2_beta = find_register_symmetry(beta_labels, "beta")
    z2_full = verify_full_hamiltonian_symmetries(qop_bare)

    print(f"\n  -- tapering the alpha register --")
    reduced_u, sq_qubit_a, U_a, tapered_val_a = taper_register(u_vecs[:K], alpha_labels, z2_alpha, K, "alpha")
    rt_err_a = verify_round_trip(u_vecs[:K], reduced_u, sq_qubit_a, U_a, tapered_val_a, K)
    print(f"  round-trip reconstruction error: {rt_err_a:.2e} ({'EXACT' if rt_err_a < 1e-9 else 'FAILED'})")
    assert rt_err_a < 1e-9

    print(f"\n  -- tapering the beta register --")
    reduced_v, sq_qubit_b, U_b, tapered_val_b = taper_register(v_vecs[:K], beta_labels, z2_beta, K, "beta")
    rt_err_b = verify_round_trip(v_vecs[:K], reduced_v, sq_qubit_b, U_b, tapered_val_b, K)
    print(f"  round-trip reconstruction error: {rt_err_b:.2e} ({'EXACT' if rt_err_b < 1e-9 else 'FAILED'})")
    assert rt_err_b < 1e-9

    print(f"\n  -- tapering the Hamiltonian's Pauli terms and re-verifying the exact energy --")
    alpha_cache_full, beta_cache_full = ef.precompute_exact_matrices(terms, u_vecs[:K], v_vecs[:K])
    E_untapered = ef.ef_energy_from_matrices(terms, lambdas, alpha_cache_full, beta_cache_full, enuc, K)
    err_untapered = abs(E_untapered - exact_energy) * HARTREE_TO_KCAL_MOL
    print(f"  untapered (3-register-qubit-equivalent, full 4-qubit reconstruction) energy check: "
          f"{err_untapered:.2e} kcal/mol")

    # build tapered alpha/beta caches directly from the REDUCED Schmidt vectors + tapered Pauli
    # operators, via DIRECT matrix conjugation with the SAME combined transform U used for the
    # state vectors (taper_pauli_matrix) -- guarantees consistency between how states and
    # operators were rotated, rather than tracking a separate symbolic Pauli-string transform.
    n_qubits_reg = 4
    alpha_cache_reduced, beta_cache_reduced = {}, {}
    for l in alpha_labels:
        Pa_reduced = taper_pauli_matrix(l, U_a, sq_qubit_a, tapered_val_a, n_qubits_reg)
        alpha_cache_reduced[l] = reduced_u.conj() @ Pa_reduced @ reduced_u.T
    for l in beta_labels:
        Pb_reduced = taper_pauli_matrix(l, U_b, sq_qubit_b, tapered_val_b, n_qubits_reg)
        beta_cache_reduced[l] = reduced_v.conj() @ Pb_reduced @ reduced_v.T

    E_tapered = ef.ef_energy_from_matrices(terms, lambdas, alpha_cache_reduced, beta_cache_reduced, enuc, K)
    err_tapered = abs(E_tapered - exact_energy) * HARTREE_TO_KCAL_MOL
    print(f"  TAPERED (3-qubit alpha, 3-qubit beta) reconstruction energy: E={E_tapered:.10f} Ha")
    print(f"  error vs exact: {err_tapered:.2e} kcal/mol "
          f"({'EXACT -- tapering is genuinely free here' if err_tapered < 1e-6 else 'FAILED -- tapering introduced error'})")

    results = {
        "exact_energy_ha": exact_energy, "K": K,
        "alpha_symmetry_generators": [str(s) for s in z2_alpha.symmetries],
        "beta_symmetry_generators": [str(s) for s in z2_beta.symmetries],
        "full_hamiltonian_symmetry_generators": [str(s) for s in z2_full.symmetries],
        "alpha_tapered_qubit": int(sq_qubit_a), "alpha_tapered_value": int(tapered_val_a),
        "beta_tapered_qubit": int(sq_qubit_b), "beta_tapered_value": int(tapered_val_b),
        "alpha_round_trip_err": rt_err_a, "beta_round_trip_err": rt_err_b,
        "energy_after_tapering_err_kcal": err_tapered,
        "tapering_is_exact": bool(err_tapered < 1e-6),
        "qubits_per_register_before": 4, "qubits_per_register_after": 3,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results, reduced_u, reduced_v, K


if __name__ == "__main__":
    main()
