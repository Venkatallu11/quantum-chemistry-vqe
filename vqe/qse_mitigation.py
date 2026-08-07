#!/usr/bin/env python3
"""
qse_mitigation.py — quantum subspace expansion (McClean, Romero,
Babbush, Aspuru-Guzik, PRA 95, 042308 (2017)) applied to this project's
K=6 forged basis, as a noise-mitigation method that needs NO channel
model. Motivated directly by iteration 9's diagnosis: CDR and PEC both
found real IonQ noise larger/structured differently than what a small
Clifford-learned Pauli channel predicted -- QSE sidesteps that whole
problem by never modeling the channel at all. Noise resilience here is
STRUCTURAL: a classical generalized-eigenvalue solve naturally
suppresses components that don't overlap well with the low-noise
subspace, rather than correcting a measured value using an assumed noise
model.

HOW THIS MAPS ONTO ENTANGLEMENT FORGING, DERIVED HERE NOT ASSUMED: the
standard forged-energy formula (entanglement_forging_h4.py
ef_energy_from_noisy_matrices) is
    E(lambda) = enuc + [ sum_n lambda_n^2 A_nn B_nn
                        + sum_{n<m} 2 lambda_n lambda_m Re(A_nm B_nm) ] / norm2
which is EXACTLY the Rayleigh quotient lambda^T H_eff lambda / (lambda^T
lambda) + enuc for the symmetric K x K matrix
    H_eff[n,m] = sum_terms coeff * Re(A_mats[a][n,m] * B_mats[b][n,m])
(diagonal entries n=m fold in trivially, A/B already real there). Every
result in this project through iteration 9 evaluated this quadratic form
at the CLASSICALLY KNOWN Schmidt coefficients lambda -- i.e. trusted that
the exact-diagonalization-derived combining weights remain optimal even
when A_mats/B_mats are noisy. QSE removes that assumption: measure
H_eff (and, in the regularized variant, an overlap matrix S) from the
SAME noisy circuits already required for entanglement forging, then let
a classical eigensolve find the best combination -- MIN eigenvalue of
H_eff c = eps*S*c is the QSE-refined energy, no re-measurement of
lambda, no channel to learn, no proximity to any target (H_eff is built
from the SAME alpha/beta Pauli measurements iteration 9 already
collects; nothing here depends on which of the 36 targets is which).

VERIFIED BEFORE TRUSTING (see verify_exact_recovery() below): this
project's own K=6 finding is that the Schmidt rank is EXACTLY 6 for this
fragment (RESEARCH_LEDGER.md, rank6_symmetry_vd.py) -- the K=6 forged
construction with lambda_known is not a truncation, it reproduces the
true ground state exactly. That means, in the NOISELESS case, H_eff's
own lowest eigenvalue MUST equal the standard forging-formula energy to
near machine precision, since lambda_known must already be H_eff's
ground-state eigenvector (the K=6 subspace contains the true ground
state, so nothing beats it). This is checked directly before any noisy
or real-hardware use -- if it fails, the H_eff construction has a bug,
full stop.

TWO VARIANTS, BOTH IMPLEMENTED, BOTH FLOOR-TESTED SEPARATELY:
  1. ORDINARY (S=I): the ONLY free choice is which K x K matrices to
     build H_eff from -- no continuous parameter to sweep, so no floor
     test applies here any more than it did to loop_pec.py's exact PEC
     (deterministic given the noisy matrices; verified instead via a
     determinism check).
  2. REGULARIZED (S measured, not assumed I): needs
     K*(K-1)/2=15 new circuits per geometry -- a "compute-uncompute"
     fidelity measurement (prepare u_n, apply the INVERSE of u_m's
     ansatz, measure P(|0000>) = |<u_m|u_n>|^2), the standard ancilla-
     free way to get a state overlap MAGNITUDE from two circuits sharing
     the same parametrized family. HONEST LIMITATION, stated once here:
     this circuit gives |<u_m|u_n>|, not its SIGN -- S[n,m] is built as
     +sqrt(P0000) for n!=m, a magnitude-only proxy, not a claim of
     measuring the true signed overlap. (A Hadamard test would recover
     the sign but needs an ancilla and controlled state-prep, a real
     circuit-count/complexity cost this file does not spend without
     first checking whether the magnitude-only proxy already helps.) The
     regularization threshold epsilon (on S's eigenvalues) IS a genuine
     free parameter and IS floor-tested below, per the task's explicit
     requirement.

A REJECTED SHORTCUT, recorded so it is not tried again: the identity
Pauli's "cross term" on the (u_n+u_m)/sqrt(2) phase circuits does NOT
measure the overlap <u_n|u_m> -- <psi|I|psi>=1 is a normalization
tautology for ANY properly normalized measured probability distribution,
true whether or not the circuit is noisy, so it carries zero information
about state overlap. This was the first idea tried while designing this
file and was caught by direct derivation before writing any code that
depended on it, not by a failed experiment.

Run:
    python vqe/qse_mitigation.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, HARTREE_TO_KCAL_MOL, floor_test, transpile_fixed,
    build_ansatz, measurement_circuit, group_labels_qubit_wise, combined_basis_label,
    pauli_expectation, P2_PER_GATE, P1_PER_GATE, energy_from_alpha_matrices,
    derive_beta_matrices, combine_matrices,
)
import rank6_symmetry_vd as r  # local noisy density-matrix machinery, reused not re-derived

from qiskit.quantum_info import Statevector
from scipy.linalg import eigh as generalized_eigh

K = 6
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "qse_mitigation_results.json")


# ---------------------------------------------------------------------------
# H_eff construction -- the core, channel-model-free piece
# ---------------------------------------------------------------------------

def build_h_eff(terms, alpha_mats, beta_mats, K):
    """H_eff[n,m] = sum_terms coeff.real * Re(A[a][n,m]*B[b][a][n,m]) --
    the symmetric K x K matrix whose Rayleigh quotient (lambda^T H_eff
    lambda)/(lambda^T lambda) reproduces the standard forging-formula
    electronic energy exactly (derived in the module docstring)."""
    H_eff = np.zeros((K, K))
    for alpha_label, beta_label, coeff in terms:
        A = alpha_mats[alpha_label][:K, :K]
        B = beta_mats[beta_label][:K, :K]
        H_eff += coeff.real * np.real(A * B)
    return H_eff


def qse_energy_ordinary(H_eff, enuc):
    """S=I ordinary eigenproblem -- min eigenvalue is the QSE-refined
    electronic energy."""
    eigvals = np.linalg.eigvalsh(H_eff)
    return float(eigvals[0]) + enuc, eigvals


def qse_energy_regularized(H_eff, S, enuc, threshold):
    """Regularized generalized eigenproblem: project out S's eigenvectors
    with eigenvalue < threshold (the standard QSE regularization recipe --
    those directions are numerically ill-conditioned / near-linearly-
    dependent under the measured S, and including them can make the
    generalized eigenproblem solve nonsense, not just noisy, results),
    solve the ordinary eigenproblem in the surviving subspace, map back."""
    s_eigvals, s_eigvecs = np.linalg.eigh(S)
    keep = s_eigvals > threshold
    n_kept = int(np.sum(keep))
    if n_kept == 0:
        return None, n_kept, s_eigvals
    V = s_eigvecs[:, keep] / np.sqrt(s_eigvals[keep])  # whitening transform
    H_proj = V.T @ H_eff @ V
    eigvals = np.linalg.eigvalsh(H_proj)
    return float(eigvals[0]) + enuc, n_kept, s_eigvals


# ---------------------------------------------------------------------------
# Overlap magnitude circuit: |<u_m|u_n>|^2 via compute-uncompute, no
# ancilla, reusing the SAME fixed-ansatz circuit family entanglement
# forging already needs.
# ---------------------------------------------------------------------------

def overlap_magnitude_circuit(angles_n, angles_m, basis_gates):
    """build_ansatz(angles_n) then the INVERSE of build_ansatz(angles_m),
    transpiled at the project's invariant optimization_level=0. Measuring
    P(all-zero bitstring) on this circuit gives |<u_m|u_n>|^2 exactly:
    U_m^dagger U_n |0> has |0>-component <0|U_m^dagger U_n|0> =
    <u_m|u_n> (since U_n|0>=|u_n>, U_m|0>=|u_m>), so the |0000> outcome
    probability is |<u_m|u_n>|^2."""
    qc_n = transpile_fixed(build_ansatz(angles_n), basis_gates)
    qc_m_inv = transpile_fixed(build_ansatz(angles_m), basis_gates).inverse()
    qc = qc_n.compose(qc_m_inv)
    qc.measure_all()
    return qc


def overlap_magnitude_local(angles_n, angles_m):
    """Exact (noiseless) |<u_m|u_n>|^2 via direct statevector overlap --
    used to VERIFY the compute-uncompute circuit construction before
    trusting it for noisy/real use."""
    sv_n = np.asarray(Statevector.from_instruction(build_ansatz(angles_n)))
    sv_m = np.asarray(Statevector.from_instruction(build_ansatz(angles_m)))
    return float(np.abs(np.vdot(sv_m, sv_n)) ** 2)


def verify_overlap_circuit(solutions, basis_gates):
    """Verify the compute-uncompute CIRCUIT construction reproduces the
    exact statevector overlap to machine precision for a few pairs,
    before this circuit is trusted for any noisy or real measurement."""
    names = sorted(solutions.keys())[:4]
    worst = 0.0
    rows = []
    for i, n in enumerate(names):
        for m in names[i + 1:]:
            qc = overlap_magnitude_circuit(solutions[n]["angles"], solutions[m]["angles"], basis_gates)
            sv = np.asarray(Statevector.from_instruction(qc.remove_final_measurements(inplace=False)))
            p0000 = float(np.abs(sv[0]) ** 2)
            exact = overlap_magnitude_local(solutions[n]["angles"], solutions[m]["angles"])
            err = abs(p0000 - exact)
            worst = max(worst, err)
            rows.append({"pair": (n, m), "circuit_p0000": p0000, "exact_overlap_sq": exact, "err": err})
    return rows, worst


# ---------------------------------------------------------------------------
# Noisy (local density-matrix) matrices + overlap, for floor-testing
# BEFORE any real submission -- reuses rank6_symmetry_vd's noise model
# and density-matrix machinery (same one iteration-1-through-8 used).
# ---------------------------------------------------------------------------

def noisy_alpha_beta_mats(p, non_id_labels, noise_model):
    raw, _ = r.measure_raw_per_slot(p, non_id_labels, noise_model)
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    beta_mats = derive_beta_matrices(alpha_mats, p["signs"], K)
    return alpha_mats, beta_mats


def noisy_overlap_matrix(p, noise_model, basis_gates):
    """S[n,m] = +sqrt(P0000) for n!=m (magnitude-only proxy, see module
    docstring), 1.0 on the diagonal, from the LOCAL noisy density-matrix
    simulator (same noise model as every other local floor test in this
    project) -- NOT yet real hardware."""
    names = [f"u_{n}" for n in range(K)]
    S = np.eye(K)
    for i in range(K):
        for j in range(i + 1, K):
            qc = overlap_magnitude_circuit(p["solutions"][names[i]]["angles"],
                                            p["solutions"][names[j]]["angles"], basis_gates)
            qc2 = qc.copy()
            qc2.remove_final_measurements()
            qc2.save_density_matrix()
            from qiskit_aer import AerSimulator
            sim = AerSimulator(method="density_matrix", noise_model=noise_model)
            result = sim.run(qc2).result()
            dm = np.asarray(result.data(0)["density_matrix"])
            p0000 = float(np.real(dm[0, 0]))
            S[i, j] = S[j, i] = np.sqrt(max(p0000, 0.0))
    return S


# ---------------------------------------------------------------------------
# Main: verify, then floor-test, all locally, before any real submission.
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 96)
    print("  qse_mitigation.py -- quantum subspace expansion, no channel model")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    from qforge import fit_all_targets, verify_constant_gate_count
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    assert n_ok == 36
    counts = verify_constant_gate_count(solutions)
    assert counts == {11}
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    BASIS_GATES = r.BASIS_GATES
    print(f"  setup OK: 36/36 targets converged, gate count={counts}, exact_energy={p['exact_energy']:.6f} Ha")

    # -- STEP 1: verify H_eff's noiseless ground eigenvalue matches the standard forging formula --
    print(f"\n  -- VERIFY 1: H_eff's noiseless ground eigenvalue == standard forging-formula energy --")
    from entanglement_forging_h4 import precompute_exact_matrices
    alpha_cache_exact, beta_cache_exact = precompute_exact_matrices(p["terms"], p["u_vecs"], p["signs"][:, None] * p["u_vecs"])
    H_eff_exact = build_h_eff(p["terms"], alpha_cache_exact, beta_cache_exact, K)
    E_qse_exact, eigvals_exact = qse_energy_ordinary(H_eff_exact, p["enuc"])
    err_exact_kcal = abs(E_qse_exact - p["noiseless_energy"]) * HARTREE_TO_KCAL_MOL
    print(f"    QSE (noiseless matrices, S=I) E={E_qse_exact:.10f} Ha")
    print(f"    standard forging formula        E={p['noiseless_energy']:.10f} Ha")
    print(f"    diff = {err_exact_kcal:.2e} kcal/mol")
    verify1_pass = err_exact_kcal < 1e-6
    assert verify1_pass, "H_eff construction is WRONG -- does not recover the known-exact energy, refusing to proceed"
    print(f"    VERIFIED: H_eff construction is correct (K=6 is exact here, so lambda_known must already be "
          "H_eff's own ground-state eigenvector -- confirmed, not assumed)")

    # -- STEP 2: verify the overlap circuit against exact statevector overlap --
    print(f"\n  -- VERIFY 2: compute-uncompute overlap circuit vs exact statevector overlap --")
    overlap_rows, overlap_worst = verify_overlap_circuit(solutions, BASIS_GATES)
    print(f"    checked {len(overlap_rows)} pairs, worst error = {overlap_worst:.2e}")
    verify2_pass = overlap_worst < 1e-8
    assert verify2_pass, "overlap circuit does not match exact statevector overlap -- refusing to proceed"
    print(f"    VERIFIED: compute-uncompute circuit correctly measures |<u_m|u_n>|^2")

    # -- STEP 3: noisy comparison -- does QSE (ordinary, S=I) beat plugging lambda_known into the SAME noisy matrices? --
    print(f"\n  -- NOISY COMPARISON: QSE (ordinary) vs standard forging formula, SAME noisy matrices --")
    noise_model = r.build_noise_model()
    alpha_mats_noisy, beta_mats_noisy = noisy_alpha_beta_mats(p, non_id_labels, noise_model)
    E_standard_noisy, err_standard_noisy = energy_from_alpha_matrices(
        alpha_mats_noisy, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
        exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    H_eff_noisy = build_h_eff(p["terms"], alpha_mats_noisy, beta_mats_noisy, K)
    E_qse_noisy, eigvals_noisy = qse_energy_ordinary(H_eff_noisy, p["enuc"])
    err_qse_noisy_kcal = abs(E_qse_noisy - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    print(f"    standard forging (lambda_known, noisy matrices): err_vs_exact={err_standard_noisy['err_vs_exact_kcal']:.4f} kcal/mol")
    print(f"    QSE ordinary (re-diagonalized, SAME noisy matrices): err_vs_exact={err_qse_noisy_kcal:.4f} kcal/mol")
    qse_helps = err_qse_noisy_kcal < err_standard_noisy["err_vs_exact_kcal"]
    print(f"    QSE {'HELPS' if qse_helps else 'does NOT help'} vs the standard formula on the SAME noisy data "
          f"({err_standard_noisy['err_vs_exact_kcal']/err_qse_noisy_kcal if err_qse_noisy_kcal > 0 else float('inf'):.2f}x)")

    # determinism check (S=I QSE has no free parameter/randomness)
    alpha_mats_noisy2, beta_mats_noisy2 = noisy_alpha_beta_mats(p, non_id_labels, noise_model)
    H_eff_noisy2 = build_h_eff(p["terms"], alpha_mats_noisy2, beta_mats_noisy2, K)
    E_qse_noisy2, _ = qse_energy_ordinary(H_eff_noisy2, p["enuc"])
    determinism_diff = abs(E_qse_noisy2 - E_qse_noisy) * HARTREE_TO_KCAL_MOL
    print(f"    determinism check (exact local simulator, no shot noise here): diff={determinism_diff:.2e} kcal/mol")

    # -- STEP 4: regularized QSE -- measure S, floor-test the threshold --
    print(f"\n  -- REGULARIZED QSE: measure S (compute-uncompute), floor-test the eigenvalue threshold --")
    S_noisy = noisy_overlap_matrix(p, noise_model, BASIS_GATES)
    S_eigvals = np.linalg.eigvalsh(S_noisy)
    print(f"    measured S eigenvalues: {np.round(S_eigvals, 6)}")
    print(f"    S off-diagonal magnitude range: [{np.min(np.abs(S_noisy - np.eye(K))):.2e}, "
          f"{np.max(np.abs(S_noisy - np.eye(K))):.2e}] (0 = perfectly orthonormal, as noiseless Schmidt vectors are)")

    # NOTE: the range here is chosen to actually CROSS the measured S
    # eigenvalue spectrum (checked below, not assumed) -- an earlier
    # version of this sweep only went up to 0.5 and every threshold gave
    # an IDENTICAL n_kept=6/6 result, a VACUOUS floor test (the plateau
    # was real but meaningless: regularization never actually removed
    # any dimension across the whole tested range). Extending past the
    # largest measured eigenvalue is what makes this floor test meaningful.
    s_lo, s_hi = float(np.min(S_eigvals)), float(np.max(S_eigvals))
    thresholds = sorted(set([1e-6, 1e-3] +
                             list(np.linspace(max(1e-2, s_lo * 0.5), s_hi * 1.2, 12).round(4))))
    threshold_errs = []
    threshold_n_kept = []
    for thresh in thresholds:
        E_reg, n_kept, _ = qse_energy_regularized(H_eff_noisy, S_noisy, p["enuc"], thresh)
        if E_reg is None:
            threshold_errs.append(None)
            threshold_n_kept.append(0)
            print(f"    threshold={thresh:.0e}: ALL directions projected out (n_kept=0) -- undefined")
            continue
        err_kcal = abs(E_reg - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
        threshold_errs.append(err_kcal)
        threshold_n_kept.append(n_kept)
        print(f"    threshold={thresh:.0e}: n_kept={n_kept}/{K}, err_vs_exact={err_kcal:.4f} kcal/mol")

    valid = [(t, e, n) for t, e, n in zip(thresholds, threshold_errs, threshold_n_kept) if e is not None]
    floor_result = floor_test([t for t, _, _ in valid], [e for _, e, _ in valid])
    print(f"\n    FLOOR TEST on regularization threshold: {floor_result['verdict']}")

    results = {
        "idea": "quantum subspace expansion (McClean et al. PRA 95, 042308) on the K=6 forged basis, "
                "H_eff derived from the SAME alpha/beta matrices entanglement forging already measures, "
                "no channel model, no training data",
        "verify1_h_eff_matches_standard_formula": {"pass": bool(verify1_pass), "diff_kcal": err_exact_kcal},
        "verify2_overlap_circuit_matches_statevector": {"pass": bool(verify2_pass), "worst_err": overlap_worst,
                                                          "rows": [{"pair": r_["pair"], "err": r_["err"]} for r_ in overlap_rows]},
        "noisy_comparison": {
            "standard_forging_err_kcal": err_standard_noisy["err_vs_exact_kcal"],
            "qse_ordinary_err_kcal": err_qse_noisy_kcal,
            "qse_helps": bool(qse_helps),
            "determinism_diff_kcal": determinism_diff,
        },
        "regularized_qse": {
            "S_eigenvalues": S_eigvals.tolist(),
            "threshold_sweep": [{"threshold": t, "err_kcal": e, "n_kept": n} for t, e, n in zip(thresholds, threshold_errs, threshold_n_kept)],
            "floor_test": floor_result,
        },
        "chemical_accuracy_kcal": 1.0,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
