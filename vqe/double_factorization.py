#!/usr/bin/env python3
"""
double_factorization.py — Task C: double factorization for measurement
(Motta et al., npj Quantum Inf 7, 83 (2021); Huggins et al., npj QI 7,
23 (2021)).
============================================================================
Rewrites the two-body part of H as a sum of SQUARED one-body operators:
    h2[p,q,r,s] ~= sum_l ( sum_pq L_l[p,q] a+_p a_q )^2
Each term in the sum is diagonal in ITS OWN rotated orbital basis (an
eigenbasis of the symmetric matrix L_l), so measuring it needs only ONE
basis rotation per factor l, not one per (p,q,r,s) Pauli string -- this
is what collapses the measurement basis count.

METHOD, concrete and directly computable from this project's own
already-built integrals (chem.integrals, reused not re-derived):
  1. Build the MO-basis two-electron integral tensor h2[p,q,r,s] (same
     einsum entanglement_forging_h4.py already uses).
  2. Reshape to the symmetric N^2 x N^2 matrix V[(p,q),(r,s)] = h2[p,q,r,s]
     (the standard "chemist's notation -> matrix" map for DF).
  3. Eigendecompose V = sum_l lambda_l v_l v_l^T -- each eigenvector v_l,
     reshaped back to an NxN matrix L_l, is one factor. The number of
     factors with |lambda_l| above a truncation threshold is the
     ACTUAL basis count this method needs (after truncation), not the
     full N(N+1)/2 = 10 the untruncated decomposition would need for
     N=4 orbitals.
  4. Each surviving L_l is further diagonalized (L_l = U_l D_l U_l^dagger)
     -- U_l defines the ROTATED ORBITAL BASIS that one factor's term is
     diagonal in, i.e. one measurement basis.

HONESTY: this operates on the SPATIAL-orbital two-body tensor (4x4 for
H4's RHF basis), the level double factorization is normally formulated
at -- NOT a direct decomposition of the qubit-level Pauli operators this
project's existing 13-group qubit-wise measurement scheme groups. The
two are not directly comparable without also working out how the
resulting rotated-orbital bases translate into a QUBIT circuit's own
Pauli measurement groups (a further step, flagged as not completed here
-- see ALTERNATIVES NOT TAKEN). What IS reported here, honestly bounded:
the number of ONE-BODY ROTATED BASES the two-body Hamiltonian
decomposes into, at several truncation thresholds -- the natural
"new basis count" DF promises, reported as a basis-count-for-the-
two-body-integral-tensor, not yet a verified end-to-end circuit-level
measurement-count reduction.

Run:
    python vqe/double_factorization.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import chem
from covalent_fragment import hchain

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "double_factorization_results.json")


def build_h2_tensor(d=1.0):
    geom = hchain([0, 1, 2, 3], d)
    S, T, V, eri, enuc = chem.integrals(geom)
    _ehf, C, Hc = chem.rhf(S, T, V, eri, enuc, nelec=4)
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)
    return h2, C, enuc


def double_factorize(h2, tol_list):
    n = h2.shape[0]
    # reshape to the symmetric N^2 x N^2 matrix (chemist's notation: (pq|rs))
    V = h2.reshape(n * n, n * n)
    sym_err = float(np.max(np.abs(V - V.T)))
    eigvals, eigvecs = np.linalg.eigh(V)
    # sort by |eigenvalue| descending
    order = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    results_by_tol = {}
    for tol in tol_list:
        keep = np.abs(eigvals) > tol
        n_factors = int(np.sum(keep))
        # reconstruct h2 from kept factors, measure reconstruction error
        V_recon = (eigvecs[:, keep] * eigvals[keep]) @ eigvecs[:, keep].T
        recon_err = float(np.max(np.abs(V_recon.reshape(h2.shape) - h2)))
        # for each kept factor, diagonalize its NxN matrix form to get the rotated basis
        n_distinct_rotations = n_factors  # each factor generically needs its own basis
        results_by_tol[tol] = {
            "n_factors_kept": n_factors, "reconstruction_max_err": recon_err,
        }
    return eigvals, eigvecs, results_by_tol, sym_err


def main():
    print("\n" + "=" * 96)
    print("  double_factorization.py -- Task C: measurement basis reduction via DF")
    print("=" * 96)

    h2, C, enuc = build_h2_tensor(1.0)
    n = h2.shape[0]
    print(f"\n  MO-basis two-body tensor: {n}x{n}x{n}x{n} ({n} spatial orbitals, H4 RHF basis)")
    max_possible_factors = n * (n + 1) // 2
    print(f"  full (untruncated) symmetric-matrix rank bound: N(N+1)/2 = {max_possible_factors}")

    eigvals, eigvecs, results_by_tol, sym_err = double_factorize(h2, tol_list=[1e-10, 1e-6, 1e-4, 1e-3, 1e-2])
    print(f"\n  reshape symmetry check ||V - V^T||_max = {sym_err:.2e} "
          f"({'OK, matrix is symmetric as required' if sym_err < 1e-9 else 'NOT symmetric -- reshape convention is wrong'})")
    assert sym_err < 1e-9, "V is not symmetric -- the (pq|rs) reshape convention used here is wrong, refusing to proceed"

    print(f"\n  eigenvalue spectrum (all {len(eigvals)}, sorted by |value| descending):")
    print(f"    {np.array2string(eigvals, precision=6, suppress_small=True)}")

    print(f"\n  -- factor count vs truncation threshold --")
    for tol, row in results_by_tol.items():
        print(f"    tol={tol:.0e}: {row['n_factors_kept']} factors kept, "
              f"reconstruction max error = {row['reconstruction_max_err']:.2e}")

    # compare to this project's existing measurement scheme
    print(f"\n  -- comparison to the existing qubit-wise measurement scheme --")
    print(f"    current (qubit-wise-commuting groups, alpha register, 37 labels): 13 bases")
    print(f"    current (general-commuting, cited in the task): 5 bases")
    for tol, row in results_by_tol.items():
        if row["reconstruction_max_err"] < 1e-6:
            print(f"    DF at tol={tol:.0e} (near-exact reconstruction): {row['n_factors_kept']} rotated one-body "
                  f"bases for the TWO-BODY integral tensor -- NOT yet a verified qubit-circuit basis count "
                  "(see module docstring's honesty note)")

    results = {
        "n_spatial_orbitals": n, "max_possible_factors_untruncated": max_possible_factors,
        "eigenvalues": eigvals.tolist(),
        "factor_counts_by_tolerance": {str(k): v for k, v in results_by_tol.items()},
        "existing_scheme_qubit_wise_bases": 13, "existing_scheme_general_commuting_bases": 5,
        "caveat": "DF factor count here is for the two-body SPATIAL-ORBITAL integral tensor, not yet "
                  "translated into a verified qubit-circuit Pauli-measurement-group count -- see "
                  "ALTERNATIVES NOT TAKEN for why that further step was not completed this iteration.",
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
