#!/usr/bin/env python3
"""
orbital_rotation_study.py — Task 3: shrink the Hamiltonian coefficient
L1/L2 norm via orbital rotation, pure classical preprocessing.
============================================================================
MOTIVATION: shot noise's cost scales with the Hamiltonian's Pauli-
coefficient norm (empirically, this project's own forged-energy estimator
scales as ~1/sqrt(N) with SOME effective norm, measured in
shot_noise_study.py to be smaller than the naive L2 prediction but still
present and real). The only classical lever on that norm, independent of
which mitigation method is used downstream, is the ORBITAL BASIS the
Hamiltonian is expressed in: any orthogonal rotation among the 4 spatial
MOs leaves the exact physics (ground-state energy) unchanged but can
change the Pauli-coefficient structure of the mapped qubit Hamiltonian
substantially. This is pure classical preprocessing -- done once, before
ANY quantum circuit is built -- and reduces both gate noise (any ansatz
built from this Hamiltonian's structure) and shot noise (smaller
coefficients need fewer shots for the same absolute energy precision)
simultaneously.

APPROACH: parametrize an orthogonal 4x4 rotation U via the matrix
exponential of an antisymmetric generator (6 free parameters for a 4x4
antisymmetric matrix, i.e. all independent orbital-pair rotation angles),
apply it to the RHF MO coefficients (C -> C @ U), recompute the one- and
two-electron integrals in the rotated basis, remap to a qubit Hamiltonian
via the SAME Jordan-Wigner pipeline used everywhere else in this project,
and minimize the resulting L1 norm via `scipy.optimize.minimize` (Nelder-
Mead, since the objective is a nonsmooth function of a discrete Pauli
decomposition's coefficients that can develop degenerate/reordered terms).

VERIFIED, NOT ASSUMED: any orthogonal rotation among the occupied+virtual
MOs is a similarity transform on the physical Hilbert space, and the exact
ground-state energy is provably invariant -- checked directly (not just
argued) by recomputing the exact ground-state energy at the OPTIMIZED
rotation and confirming it matches the untouched-basis value to machine
precision.

Run:
    python vqe/orbital_rotation_study.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
import chem
from covalent_fragment import hchain
from molecules_real import _constrain_particle_number
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper
from scipy.sparse.linalg import eigsh

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "orbital_rotation_study_results.json")
N_ORB = 4  # H4, STO-3G: 4 spatial orbitals
N_PARAMS = N_ORB * (N_ORB - 1) // 2  # 6 independent Givens angles


def rotation_matrix(params):
    A = np.zeros((N_ORB, N_ORB))
    idx = 0
    for i in range(N_ORB):
        for j in range(i + 1, N_ORB):
            A[i, j] = params[idx]
            A[j, i] = -params[idx]
            idx += 1
    return expm(A)


def qop_and_norms(C, Hc, eri, enuc, nelec):
    h1 = C.T @ Hc @ C
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)
    ee = ElectronicEnergy.from_raw_integrals(h1, h2)
    qop_bare = JordanWignerMapper().map(ee.second_q_op())
    coeffs = np.array([complex(c).real for _, c in qop_bare.to_list()])
    L1 = float(np.sum(np.abs(coeffs)))
    L2 = float(np.sqrt(np.sum(coeffs ** 2)))
    qop_pen = _constrain_particle_number(qop_bare, nelec)
    return qop_bare, qop_pen, L1, L2


def exact_energy_of(qop_pen, enuc):
    M = qop_pen.to_matrix(sparse=True)
    val, _ = eigsh(M, k=1, which="SA")
    return float(val[0]) + enuc


def main():
    print("\n" + "=" * 78)
    print("  orbital_rotation_study.py -- shrink the Hamiltonian coefficient norm")
    print("=" * 78)

    ATOMS, NELEC, D = [0, 1, 2, 3], 4, 1.0
    geom = hchain(ATOMS, D)
    S, T, V, eri, enuc = chem.integrals(geom)
    _ehf, C0, Hc = chem.rhf(S, T, V, eri, enuc, nelec=NELEC)

    qop0_bare, qop0_pen, L1_0, L2_0 = qop_and_norms(C0, Hc, eri, enuc, NELEC)
    E0 = exact_energy_of(qop0_pen, enuc)
    print(f"\n  baseline (RHF orbitals, untouched): L1={L1_0:.4f} Ha, L2={L2_0:.4f} Ha")
    print(f"  baseline exact ground-state energy: {E0:.6f} Ha")

    def objective(params):
        U = rotation_matrix(params)
        C = C0 @ U
        _, _, L1, _ = qop_and_norms(C, Hc, eri, enuc, NELEC)
        return L1

    print(f"\n  optimizing {N_PARAMS} orbital-rotation angles to minimize L1 (Nelder-Mead)...")
    print(f"  (~0.45s per objective evaluation -- budget capped to keep this tractable; a modest, "
          "reported-as-partial optimum is a legitimate result even if not the global minimum)")
    best_result = None
    for attempt in range(3):
        rng = np.random.default_rng(attempt)
        x0 = rng.uniform(-0.3, 0.3, N_PARAMS) if attempt > 0 else np.zeros(N_PARAMS)
        t0 = time.time()
        res = minimize(objective, x0, method="Nelder-Mead",
                        options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 150, "maxfev": 150})
        print(f"    attempt {attempt}: L1={res.fun:.4f} (start L1={objective(x0):.4f}), "
              f"{res.nfev} evals, {time.time()-t0:.1f}s", flush=True)
        if best_result is None or res.fun < best_result.fun:
            best_result = res

    U_best = rotation_matrix(best_result.x)
    C_best = C0 @ U_best
    qop_bare_opt, qop_pen_opt, L1_opt, L2_opt = qop_and_norms(C_best, Hc, eri, enuc, NELEC)
    print(f"\n  best found: L1={L1_opt:.4f} Ha (was {L1_0:.4f}), L2={L2_opt:.4f} Ha (was {L2_0:.4f})")
    print(f"  reduction: L1 {L1_0/L1_opt:.3f}x, L2 {L2_0/L2_opt:.3f}x")

    # verify physics invariance
    E_opt = exact_energy_of(qop_pen_opt, enuc)
    energy_diff_kcal = abs(E_opt - E0) * 627.5094740631
    print(f"\n  VERIFY: exact ground-state energy in rotated basis = {E_opt:.8f} Ha "
          f"(baseline {E0:.8f} Ha, diff {energy_diff_kcal:.2e} kcal/mol)")
    physics_preserved = energy_diff_kcal < 1e-6
    print(f"  verdict: {'PASS -- physics exactly preserved, as required for ANY orbital rotation' if physics_preserved else 'FAIL -- investigate'}")

    # honest shot-budget implication: std(E) scales with the norm, so shots
    # needed for a fixed precision scale as (norm_ratio)^2
    shot_reduction_factor_L1 = (L1_0 / L1_opt) ** 2
    shot_reduction_factor_L2 = (L2_0 / L2_opt) ** 2
    print(f"\n  IMPLIED SHOT-BUDGET IMPACT (variance ~ norm^2/N, so shots for fixed precision ~ norm^2):")
    print(f"  L1-based reduction factor: {shot_reduction_factor_L1:.3f}x fewer shots needed")
    print(f"  L2-based reduction factor: {shot_reduction_factor_L2:.3f}x fewer shots needed")
    print(f"  NOTE: this project's OWN empirical variance (shot_noise_study.py) did not match the naive")
    print(f"  L2/sqrt(N) prediction (measured ~7x SMALLER already, from the bilinear beta=S.alpha.S reuse")
    print(f"  structure) -- so this reduction factor is a plausible, not verified-on-the-full-pipeline,")
    print(f"  estimate. Reported as a specification, not a re-measured shot count.")

    results = {
        "n_params": N_PARAMS,
        "baseline": {"L1": L1_0, "L2": L2_0, "exact_energy_ha": E0},
        "optimized": {"L1": L1_opt, "L2": L2_opt, "exact_energy_ha": E_opt,
                      "rotation_params": best_result.x.tolist()},
        "physics_preserved_kcal_diff": energy_diff_kcal,
        "physics_preserved": bool(physics_preserved),
        "reduction_factor": {"L1": L1_0 / L1_opt, "L2": L2_0 / L2_opt},
        "implied_shot_budget_reduction": {"via_L1": shot_reduction_factor_L1, "via_L2": shot_reduction_factor_L2},
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
