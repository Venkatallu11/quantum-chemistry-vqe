#!/usr/bin/env python3
"""
cs_vqe_h4.py — Task B: Contextual Subspace VQE (Kirby, Tranter, Love,
Quantum 5, 456 (2021)) for the FULL, monolithic 8-qubit H4 Jordan-Wigner
Hamiltonian, using the reference `symmer` package (same research group).
============================================================================
WHY THE FULL HAMILTONIAN, NOT THE ENTANGLEMENT-FORGING ALPHA REGISTER:
every other technique in this project's 19 completed iterations operates
on entanglement forging's bipartite (alpha/beta 4-qubit register) Schmidt
decomposition of the ground state. CS-VQE, as published, partitions a
MONOLITHIC qubit Hamiltonian into a noncontextual (classically solvable)
generating set plus a smaller quantum remainder -- it has no natural
bipartite-register analogue, and the "alpha register" used throughout
this project is not itself a standalone physical Hamiltonian (its Pauli
labels come from decomposing the FULL Hamiltonian's terms into
(alpha_label, beta_label) pairs for the forging reconstruction, not from
diagonalizing a real alpha-only system). Applying CS-VQE faithfully
means building it against the real, full, un-forged H4 problem -- a
genuinely different circuit paradigm from everything else in this
project, not a variant of the tapered/abstract circuits already tested.

VERIFIED, NOT ASSUMED, before trusting any partition: symmer's
PauliwordOp.from_qiskit() converts the SparsePauliOp directly (avoiding
any risk of a hand-rolled label-convention bug -- Z2 tapering's own
early development in this project hit exactly this class of bug, see
z2_tapering.py's X-vs-Z basis mixup). The resulting operator's ground
state (via to_sparse_matrix + scipy eigsh) is checked against this
project's own independent exact diagonalization
(ef_fragment.exact_ground_state, dense scipy.linalg.eigh) before
anything else is trusted.

Run:
    python vqe/cs_vqe_h4.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
from qforge import HARTREE_TO_KCAL_MOL

from symmer import PauliwordOp, ContextualSubspace
from scipy.sparse.linalg import eigsh

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "cs_vqe_h4_results.json")


def main():
    print("\n" + "=" * 96)
    print("  cs_vqe_h4.py -- Task B: Contextual Subspace VQE for the full 8-qubit H4 Hamiltonian")
    print("=" * 96)

    qop_bare, qop_pen, enuc = effrag.build_fragment_qop([0, 1, 2, 3], nelec=4, d=1.0)
    n_qubits = qop_bare.num_qubits
    print(f"  full Hamiltonian: {n_qubits} qubits, {len(qop_bare)} Pauli terms, enuc={enuc:.6f} Ha")

    e_elec_exact, psi_exact = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec_exact + enuc
    print(f"  this project's own exact diagonalization: electronic={e_elec_exact:.8f} Ha, "
          f"total={exact_energy:.8f} Ha")

    H = PauliwordOp.from_qiskit(qop_bare)
    print(f"  symmer PauliwordOp built: {H.n_qubits} qubits, {H.n_terms} terms")

    print("\n  -- verifying symmer's Hamiltonian matches this project's own exact diagonalization --")
    Hmat = H.to_sparse_matrix
    vals, vecs = eigsh(Hmat, k=1, which="SA")
    symmer_electronic = float(vals[0])
    symmer_total = symmer_electronic + enuc
    verify_err_kcal = abs(symmer_total - exact_energy) * HARTREE_TO_KCAL_MOL
    print(f"  symmer ground state: electronic={symmer_electronic:.8f} Ha, total={symmer_total:.8f} Ha")
    print(f"  match vs this project's own exact diagonalization: {verify_err_kcal:.6e} kcal/mol")
    assert verify_err_kcal < 1e-4, (
        f"symmer's Hamiltonian does NOT match this project's own exact diagonalization "
        f"({verify_err_kcal:.3e} kcal/mol off) -- refusing to trust any partition built on it"
    )
    print("  VERIFIED: symmer's Hamiltonian is the same physical operator, to machine precision.")

    print("\n  -- noncontextuality structure of the full Hamiltonian --")
    is_nc = H.is_noncontextual
    print(f"  is the FULL Hamiltonian itself noncontextual? {is_nc} "
          f"(expected False for a real molecular Hamiltonian -- if True, VQE isn't even needed)")

    results = {
        "n_qubits": n_qubits, "n_terms": len(qop_bare), "enuc": enuc,
        "exact_energy_ha": exact_energy, "symmer_energy_ha": symmer_total,
        "verify_err_kcal": verify_err_kcal, "full_hamiltonian_is_noncontextual": bool(is_nc),
    }

    print("\n  -- building the Contextual Subspace partition, sweeping n_qubits (quantum remainder size) --")
    cs_results = {}
    for n_q in range(0, n_qubits):
        try:
            cs = ContextualSubspace(H, noncontextual_strategy="StabilizeFirst")
            cs.update_stabilizers(n_qubits=n_q)
            H_cs = cs.project_onto_subspace()
            if H_cs.n_qubits == 0:
                cs_energy = float(np.real(H_cs.coeff_vec[0])) if H_cs.n_terms else 0.0
            else:
                Hcs_mat = H_cs.to_sparse_matrix
                if Hcs_mat.shape[0] <= 2:
                    Hcs_dense = np.asarray(Hcs_mat.todense())
                    cs_vals = np.linalg.eigvalsh(Hcs_dense)
                    cs_energy = float(np.min(cs_vals))
                else:
                    cs_vals, _ = eigsh(Hcs_mat, k=1, which="SA")
                    cs_energy = float(cs_vals[0])
            cs_total = cs_energy + enuc
            err_kcal = abs(cs_total - exact_energy) * HARTREE_TO_KCAL_MOL
            cs_results[n_q] = {"quantum_qubits": n_q, "cs_energy_ha": cs_total, "err_vs_exact_kcal": err_kcal}
            print(f"    quantum_remainder={n_q} qubits: E={cs_total:.6f} Ha, err_vs_exact={err_kcal:.3f} kcal/mol "
                  f"({'PASS chemical accuracy' if err_kcal < 1.0 else ''})")
        except Exception as e:
            cs_results[n_q] = {"quantum_qubits": n_q, "error": str(e)}
            print(f"    quantum_remainder={n_q} qubits: FAILED -- {e}")

    results["contextual_subspace_sweep"] = cs_results
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
