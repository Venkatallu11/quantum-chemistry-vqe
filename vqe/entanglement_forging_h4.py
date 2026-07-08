#!/usr/bin/env python3
"""
entanglement_forging_h4.py — Entanglement Forging (EF) for the H4 fragment.
============================================================================
Entanglement forging is a way to simulate a big entangled molecule with two
SMALL quantum registers instead of one big one. The trick: split the qubits
into two halves (here: alpha-spin orbitals vs beta-spin orbitals), write the
true ground state as a Schmidt decomposition

    |psi> = sum_n  lambda_n  |u_n>_alpha  |v_n>_beta

and reconstruct any expectation value <psi|H|psi> from small pieces measured
SEPARATELY on the alpha register (4 qubits) and the beta register (4 qubits),
combined classically. You never need all 8 qubits alive at once, and you
never need them entangled with each other on real hardware -- the "cross"
terms below stand in for that entanglement, computed classically after the
fact from a handful of small circuits.

This script asks one question: if fragmentation gets you the qubit count
down (8 qubits -> two 4-qubit registers), does entanglement forging survive
HARDWARE NOISE on those small registers?

Stage 1 (noiseless): use exact (numpy) Schmidt vectors and exact (numpy)
matrix elements. This isolates whether the EF bipartite bookkeeping itself
is correct, with zero noise anywhere. Must reproduce known reference errors
before Stage 2 is trusted.

Stage 2 (noisy): SAME exact Schmidt vectors, but the alpha/beta matrix
elements are now measured from actual 4-qubit circuits (StatePreparation +
the standard EF superposition trick for off-diagonal terms) run through a
qiskit-aer Quantinuum-like depolarizing noise model. This isolates whether
noise on the small registers erodes the reconstruction.

HONEST LIMITATION (read before trusting the numbers): both stages use the
EXACT Schmidt vectors taken from full diagonalization of the real H4
Hamiltonian -- not a variationally trained ansatz. This isolates the noise
question from the optimization question. A full EF-VQE (variationally
optimizing the two 4-qubit circuits instead of reading them off exact
diagonalization) is the natural next step and is NOT done here.

Run:
    python vqe/entanglement_forging_h4.py
"""
import os
import sys
import json
import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(__file__))
import chem
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit.quantum_info import Pauli
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit import transpile
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error

from covalent_fragment import hchain
from molecules_real import _constrain_particle_number

HARTREE_TO_KCAL_MOL = 627.5094740631
K_VALUES = [1, 3, 5, 8]
REFERENCE_ERRORS_HA = {1: 0.0648, 3: 0.0044, 5: 0.0009, 8: 0.0}
STAGE1_TOLERANCE_HA = 0.0015  # how close Stage 1 must land to the reference errors to be trusted
BASIS_GATES = ["u3", "cx"]
QUANTINUUM_TWO_Q_ERROR = 0.002
QUANTINUUM_ONE_Q_ERROR = 0.00005


# ---------------------------------------------------------------------------
# Setup: H4 fragment Hamiltonian, exact ground state, Schmidt decomposition
# ---------------------------------------------------------------------------

def build_h4_qop(d=1.0):
    """H4 fragment (atoms 0-3, d=1.0 Ang), full space, 8 qubits (4 spatial
    orbitals x 2 spins). Same integral -> RHF -> qiskit-nature pipeline as
    the rest of this repo (molecules_real.py, hardware_fragmentation.py)."""
    geom = hchain([0, 1, 2, 3], d)
    nelec = 4
    S, T, V, eri, enuc = chem.integrals(geom)
    _ehf, C, Hc = chem.rhf(S, T, V, eri, enuc, nelec=nelec)
    h1 = C.T @ Hc @ C
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)
    ee = ElectronicEnergy.from_raw_integrals(h1, h2)
    qop_bare = JordanWignerMapper().map(ee.second_q_op())
    qop_penalized = _constrain_particle_number(qop_bare, nelec)
    return qop_bare, qop_penalized, enuc


def exact_ground_state(qop_penalized):
    """Full 8-qubit exact diagonalization (penalized only to select the
    correct-electron-count sector -- penalty is exactly 0 on the true
    eigenstate, so the eigenvalue returned IS the bare electronic energy)."""
    M = qop_penalized.to_matrix(sparse=True)
    val, vec = eigsh(M, k=1, which="SA")
    return float(val[0]), vec[:, 0]


def schmidt_decompose(psi):
    """
    Bipartite Schmidt decomposition across qubits 0-3 (alpha) / 4-7 (beta).

    Qiskit statevector index i decomposes as i = beta_index*16 + alpha_index
    (qubit 0 is the least-significant bit), so reshape(16, 16) directly
    gives rows=beta, cols=alpha -- exactly the bipartition requested.

    psi = sum_n lambda_n |u_n>_alpha |v_n>_beta  <=>  mat[b,a] = sum_n
    lambda_n * v_n[b] * u_n[a] (a plain bilinear outer-product expansion,
    no conjugates). SVD mat = U @ diag(S) @ Vh gives mat = sum_n U[:,n]
    S[n] Vh[n,:], so matching term-by-term: v_n = U[:, n], u_n = Vh[n, :]
    directly (Vh is already the conjugate-transposed matrix -- no further
    conjugation), lambda_n = S[n] (already sorted descending by SVD).
    """
    mat = psi.reshape(16, 16)  # rows = beta (qubits 4-7), cols = alpha (qubits 0-3)
    U, S, Vh = np.linalg.svd(mat)
    lambdas = S
    u_vecs = Vh   # u_vecs[n] = alpha Schmidt vector n (16-dim, 4 qubits)
    v_vecs = U.T  # v_vecs[n] = beta Schmidt vector n  (16-dim, 4 qubits)
    return lambdas, u_vecs, v_vecs


def decompose_pauli_terms(qop_bare):
    """
    Every Pauli string on 8 qubits is a pure tensor product P_beta (x) P_alpha
    (no sums inside one term), so it factorizes exactly. Qiskit's label
    convention is left-to-right = qubit(n-1)...qubit0, so for an 8-char
    label: label[:4] acts on qubits 7,6,5,4 (beta), label[4:] acts on
    qubits 3,2,1,0 (alpha) -- matching the local bit-ordering already baked
    into u_vecs/v_vecs above.
    """
    terms = []
    for label, coeff in qop_bare.to_list():
        beta_label = label[:4]
        alpha_label = label[4:]
        terms.append((alpha_label, beta_label, complex(coeff)))
    return terms


# ---------------------------------------------------------------------------
# STAGE 1 -- exact (numpy) matrix elements
# ---------------------------------------------------------------------------

def precompute_exact_matrices(terms, u_top, v_top):
    """
    For every unique alpha/beta Pauli label appearing in the Hamiltonian,
    compute the FULL k_max x k_max matrix of exact matrix elements
    <u_n|P|u_m}> (numpy), once, reused for every K.
    """
    k_max = u_top.shape[0]
    alpha_cache, beta_cache = {}, {}
    for label in set(a for a, _, _ in terms):
        Pa = Pauli(label).to_matrix()
        alpha_cache[label] = u_top.conj() @ Pa @ u_top.T
    for label in set(b for _, b, _ in terms):
        Pb = Pauli(label).to_matrix()
        beta_cache[label] = v_top.conj() @ Pb @ v_top.T
    return alpha_cache, beta_cache


def ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K):
    """
    EF energy for the top-K Schmidt terms: enuc + (1/norm2_K) * sum over
    Pauli terms of (diagonal + cross).

    Truncating the Schmidt sum to K < full rank leaves an UNNORMALIZED
    state (sum_{n<K} lambda_n^2 < 1), so <psi_trunc|H|psi_trunc> must be
    divided by <psi_trunc|psi_trunc> = norm2_K to get a valid energy
    expectation value -- enuc is a classical constant added after, not
    part of the state being renormalized.
    """
    norm2 = sum(lambdas[n] ** 2 for n in range(K))
    E_elec = 0.0
    for alpha_label, beta_label, coeff in terms:
        Amat = alpha_cache[alpha_label][:K, :K]
        Bmat = beta_cache[beta_label][:K, :K]
        diag = sum(lambdas[n] ** 2 * Amat[n, n] * Bmat[n, n] for n in range(K))
        cross = sum(
            2 * lambdas[n] * lambdas[m] * np.real(Amat[n, m] * Bmat[n, m])
            for n in range(K) for m in range(K) if n < m
        )
        E_elec += coeff.real * (diag.real + cross)
    return E_elec / norm2 + enuc


def run_stage1(terms, lambdas, u_vecs, v_vecs, enuc, exact_energy):
    print("\n" + "=" * 70)
    print("  STAGE 1 -- noiseless EF reconstruction (exact numpy matrix elements)")
    print("=" * 70)
    k_max = max(K_VALUES)
    u_top, v_top = u_vecs[:k_max], v_vecs[:k_max]
    alpha_cache, beta_cache = precompute_exact_matrices(terms, u_top, v_top)

    results = {}
    all_pass = True
    for K in K_VALUES:
        E = ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)
        err = float(abs(E - exact_energy))
        ref = REFERENCE_ERRORS_HA[K]
        ok = bool(abs(err - ref) < STAGE1_TOLERANCE_HA)
        all_pass &= ok
        tag = "OK" if ok else "MISMATCH"
        print(f"  K={K}: E = {E:.6f} Ha, error = {err:.6f} Ha "
              f"(reference {ref:.4f} Ha) [{tag}]")
        results[str(K)] = {"energy_ha": round(E, 6), "error_ha": round(err, 6),
                            "reference_error_ha": ref, "match": ok}

    if not all_pass:
        print("\n  STAGE 1 DID NOT MATCH THE REFERENCE ERRORS -- stopping before Stage 2.")
    return results, all_pass


# ---------------------------------------------------------------------------
# STAGE 2 -- noisy (circuit) matrix elements
# ---------------------------------------------------------------------------

def build_noise_model():
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(QUANTINUUM_TWO_Q_ERROR, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(QUANTINUUM_ONE_Q_ERROR, 1), "u3")
    return nm


def state_prep_circuit(vec):
    qc = QuantumCircuit(4)
    qc.append(StatePreparation(np.asarray(vec)), range(4))
    return qc


def transpiled_state_prep_circuit(vec):
    """Aer's noisy simulator needs circuits already decomposed into the
    noise model's basis gates -- StatePreparation itself isn't a basis
    instruction the noise model knows about."""
    return transpile(state_prep_circuit(vec), basis_gates=BASIS_GATES, optimization_level=1)


def noisy_expectation(circuits, observables_per_circuit, estimator):
    """Run a batch of (circuit, [observables]) pubs through the noisy
    estimator, return a list (one per circuit) of dict label -> value."""
    pubs = list(zip(circuits, observables_per_circuit))
    result = estimator.run(pubs).result()
    out = []
    for i, obs_list in enumerate(observables_per_circuit):
        evs = np.atleast_1d(result[i].data.evs)
        out.append(dict(zip(obs_list, evs)))
    return out


def measure_diagonal_noisy(vecs, unique_labels, estimator, K):
    """<u_n|P|u_n> (or v_n) for n in 0..K-1, all unique labels batched per circuit."""
    label_list = sorted(unique_labels)
    obs_list = [Pauli(l) for l in label_list]  # SparsePauliOp accepted via Pauli too
    circuits = [transpiled_state_prep_circuit(vecs[n]) for n in range(K)]
    per_circuit_obs = [obs_list for _ in range(K)]
    raw = noisy_expectation(circuits, per_circuit_obs, estimator)
    # raw[n] keyed by Pauli object; rekey by label string
    return [{label_list[i]: list(d.values())[i] for i in range(len(label_list))} for d in raw]


def measure_cross_noisy(vecs, unique_labels, estimator, K):
    """
    <u_n|P|u_m> (or v) for all n<m in 0..K-1, via the EF superposition trick:
    prepare (|u_n> + i^k|u_m>)/sqrt(2) for k=0,1,2,3, measure P, then
        Re(x) = (E0 - E2)/2,  Im(x) = (E3 - E1)/2
    (E0..E3 are <phi_k|P|phi_k>; this symmetric combination averages over
    both signs of the phase kickback, which is the "combine with phases"
    step -- only k=0,1 are independent, k=2,3 exist for that averaging).
    """
    label_list = sorted(unique_labels)
    obs_list = [Pauli(l) for l in label_list]
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]

    circuits, per_circuit_obs = [], []
    for (n, m) in pairs:
        for k in range(4):
            vec_k = (vecs[n] + (1j ** k) * vecs[m]) / np.sqrt(2)
            circuits.append(transpiled_state_prep_circuit(vec_k))
            per_circuit_obs.append(obs_list)

    raw = noisy_expectation(circuits, per_circuit_obs, estimator)

    out = {}
    for idx, (n, m) in enumerate(pairs):
        E = [np.array(list(raw[idx * 4 + k].values())) for k in range(4)]
        re = (E[0] - E[2]) / 2
        im = (E[3] - E[1]) / 2
        x = re + 1j * im
        out[(n, m)] = dict(zip(label_list, x))
    return out


def build_noisy_matrices(vecs, terms_labels, estimator, K):
    """Assemble the full K x K Hermitian matrix of noisy matrix elements
    for every unique Pauli label, from the diagonal + cross measurements."""
    unique_labels = sorted(set(terms_labels))
    diag = measure_diagonal_noisy(vecs, unique_labels, estimator, K)
    cross = measure_cross_noisy(vecs, unique_labels, estimator, K)

    matrices = {label: np.zeros((K, K), dtype=complex) for label in unique_labels}
    for n in range(K):
        for label in unique_labels:
            matrices[label][n, n] = diag[n][label]
    for (n, m), vals in cross.items():
        for label in unique_labels:
            matrices[label][n, m] = vals[label]
            matrices[label][m, n] = np.conj(vals[label])
    return matrices


def ef_energy_from_noisy_matrices(terms, lambdas, alpha_mats, beta_mats, enuc, K):
    """Same normalization as ef_energy_from_matrices (divide by norm2_K
    before adding enuc) -- kept consistent so noiseless vs noisy K=5 are
    an apples-to-apples comparison."""
    norm2 = sum(lambdas[n] ** 2 for n in range(K))
    E_elec = 0.0
    for alpha_label, beta_label, coeff in terms:
        Amat = alpha_mats[alpha_label]
        Bmat = beta_mats[beta_label]
        diag = sum(lambdas[n] ** 2 * Amat[n, n] * Bmat[n, n] for n in range(K))
        cross = sum(
            2 * lambdas[n] * lambdas[m] * np.real(Amat[n, m] * Bmat[n, m])
            for n in range(K) for m in range(K) if n < m
        )
        E_elec += coeff.real * (diag.real + cross)
    return E_elec / norm2 + enuc


def run_stage2(terms, lambdas, u_vecs, v_vecs, enuc, exact_energy, noiseless_k5_energy):
    print("\n" + "=" * 70)
    print("  STAGE 2 -- EF under Quantinuum-like noise (K=5, 4-qubit circuits)")
    print(f"  Noise: 2-qubit depolarizing {QUANTINUUM_TWO_Q_ERROR}, "
          f"1-qubit depolarizing {QUANTINUUM_ONE_Q_ERROR}")
    print("=" * 70)

    K = 5
    nm = build_noise_model()
    estimator = AerEstimatorV2(options={"backend_options": {"noise_model": nm, "method": "density_matrix"}})

    alpha_labels = [a for a, _, _ in terms]
    beta_labels = [b for _, b, _ in terms]

    print("  Measuring alpha-register (qubits 0-3) matrix elements under noise...")
    alpha_mats = build_noisy_matrices(u_vecs, alpha_labels, estimator, K)
    print("  Measuring beta-register (qubits 4-7) matrix elements under noise...")
    beta_mats = build_noisy_matrices(v_vecs, beta_labels, estimator, K)

    E_noisy = ef_energy_from_noisy_matrices(terms, lambdas, alpha_mats, beta_mats, enuc, K)
    err_noisy_ha = abs(E_noisy - exact_energy)
    err_noiseless_ha = abs(noiseless_k5_energy - exact_energy)

    # --- transpiled 2-qubit gate counts (representative circuits) ---
    sp_circuit = state_prep_circuit(u_vecs[0])
    sp_t = transpile(sp_circuit, basis_gates=BASIS_GATES, optimization_level=1)
    cx_state_prep = sp_t.count_ops().get("cx", 0)

    cross_vec = (u_vecs[0] + 1j * u_vecs[1]) / np.sqrt(2)
    cross_circuit = state_prep_circuit(cross_vec)
    cross_t = transpile(cross_circuit, basis_gates=BASIS_GATES, optimization_level=1)
    cx_cross = cross_t.count_ops().get("cx", 0)

    print(f"\n  exact ground state           = {exact_energy:.6f} Ha")
    print(f"  noiseless EF (K=5)            = {noiseless_k5_energy:.6f} Ha  "
          f"(error {err_noiseless_ha:.6f} Ha, {err_noiseless_ha * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print(f"  Quantinuum-noise EF (K=5)     = {E_noisy:.6f} Ha  "
          f"(error {err_noisy_ha:.6f} Ha, {err_noisy_ha * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print(f"\n  transpiled 2-qubit gates -- state-prep circuit : {cx_state_prep}")
    print(f"  transpiled 2-qubit gates -- cross-term circuit : {cx_cross}")

    return {
        "K": K,
        "exact_energy_ha": round(exact_energy, 6),
        "noiseless_ef_k5_ha": round(noiseless_k5_energy, 6),
        "noiseless_ef_k5_error_ha": round(err_noiseless_ha, 6),
        "noiseless_ef_k5_error_kcal_mol": round(err_noiseless_ha * HARTREE_TO_KCAL_MOL, 4),
        "quantinuum_noise_ef_k5_ha": round(E_noisy, 6),
        "quantinuum_noise_ef_k5_error_ha": round(err_noisy_ha, 6),
        "quantinuum_noise_ef_k5_error_kcal_mol": round(err_noisy_ha * HARTREE_TO_KCAL_MOL, 4),
        "noise_model": {"two_qubit_depolarizing": QUANTINUUM_TWO_Q_ERROR,
                         "one_qubit_depolarizing": QUANTINUUM_ONE_Q_ERROR,
                         "basis_gates": BASIS_GATES},
        "transpiled_2q_gate_counts": {"state_prep_circuit": cx_state_prep,
                                       "cross_term_circuit": cx_cross},
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 70)
    print("  Entanglement Forging -- H4 fragment (atoms 0-3, d=1.0 Ang)")
    print("  Tests whether fragmentation + EF survives hardware noise.")
    print("=" * 70)

    qop_bare, qop_penalized, enuc = build_h4_qop(d=1.0)
    n_qubits = qop_bare.num_qubits
    print(f"\n  qubits = {n_qubits}, nuclear repulsion enuc = {enuc:.6f} Ha")

    e_elec, psi = exact_ground_state(qop_penalized)
    exact_energy = e_elec + enuc
    print(f"  exact ground state (eigsh) = {exact_energy:.6f} Ha")

    lambdas, u_vecs, v_vecs = schmidt_decompose(psi)
    print(f"  Schmidt coefficients (top 8): {np.round(lambdas[:8], 6)}")

    terms = decompose_pauli_terms(qop_bare)
    print(f"  {len(terms)} Pauli terms in the qubit Hamiltonian")

    stage1_results, stage1_ok = run_stage1(terms, lambdas, u_vecs, v_vecs, enuc, exact_energy)

    results = {
        "molecule": "H4 fragment (atoms 0-3, d=1.0 Ang)",
        "n_qubits": n_qubits,
        "enuc_ha": round(enuc, 6),
        "exact_energy_ha": round(exact_energy, 6),
        "schmidt_coefficients": [round(float(x), 6) for x in lambdas[:8]],
        "n_pauli_terms": len(terms),
        "stage1_noiseless_ef": stage1_results,
        "stage1_passed": stage1_ok,
    }

    if not stage1_ok:
        out = os.path.join(os.path.dirname(__file__), "entanglement_forging_h4_results.json")
        with open(out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\n  Results (Stage 1 only, FAILED) saved -> {out}")
        sys.exit(1)

    print("\n  Stage 1 PASSED -- matches reference errors. Proceeding to Stage 2.")

    noiseless_k5 = ef_energy_from_matrices(
        terms, lambdas,
        *precompute_exact_matrices(terms, u_vecs[:8], v_vecs[:8]),
        enuc, K=5,
    )
    stage2_results = run_stage2(terms, lambdas, u_vecs, v_vecs, enuc, exact_energy, noiseless_k5)
    results["stage2_quantinuum_noise_ef"] = stage2_results

    print("\n" + "=" * 70)
    print("  NOTE: both stages use the EXACT Schmidt vectors from full")
    print("  diagonalization, not a variationally trained ansatz -- this")
    print("  isolates the noise question from the optimization question.")
    print("  A full EF-VQE (variationally optimizing the two 4-qubit")
    print("  circuits) is the natural next step, not done here.")
    print("=" * 70)

    out = os.path.join(os.path.dirname(__file__), "entanglement_forging_h4_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
