#!/usr/bin/env python3
"""
entanglement_forging_zne.py — Zero-Noise Extrapolation (ZNE) on top of
entanglement forging for the H4 fragment.
============================================================================
entanglement_forging_h4.py showed that EF under a Quantinuum-like noise
model (K=5, 4-qubit registers) lands ~20 kcal/mol off the exact answer --
nowhere near chemical accuracy (1 kcal/mol).

ZNE is a noise-mitigation trick: deliberately run the SAME circuits at
several artificially amplified noise strengths (scale = 1x, 2x, 3x the
real device error rates), then fit energy vs scale and extrapolate back
to scale = 0 (the noiseless limit) -- without ever needing a
noise-free device.

This reuses every building block from entanglement_forging_h4.py (the
H4 Hamiltonian, exact Schmidt vectors, Pauli-term decomposition, the
noisy-matrix-element machinery) and only changes one thing: the noise
model's error rates are scaled by s in {1, 2, 3} before measuring, then
a straight-line and a quadratic fit to (scale, energy) are each
extrapolated to scale=0.

Run:
    python vqe/entanglement_forging_zne.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error

HARTREE_TO_KCAL_MOL = 627.5094740631
K = 5
SCALES = [1.0, 2.0, 3.0]


def build_scaled_estimator(scale):
    """Same Quantinuum-like noise model as entanglement_forging_h4.py, but
    with every depolarizing error rate multiplied by `scale`. scale=1.0
    reproduces the original Stage 2 noise model exactly."""
    nm = NoiseModel(basis_gates=ef.BASIS_GATES)
    nm.add_all_qubit_quantum_error(
        depolarizing_error(ef.QUANTINUUM_TWO_Q_ERROR * scale, 2), "cx"
    )
    nm.add_all_qubit_quantum_error(
        depolarizing_error(ef.QUANTINUUM_ONE_Q_ERROR * scale, 1), "u3"
    )
    return AerEstimatorV2(
        options={"backend_options": {"noise_model": nm, "method": "density_matrix"}}
    )


def main():
    print("\n" + "=" * 70)
    print("  Zero-Noise Extrapolation on Entanglement Forging -- H4 fragment")
    print("=" * 70)

    qop_bare, qop_penalized, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_penalized)
    exact = e_elec + enuc
    print(f"\n  exact ground state = {exact:.6f} Ha")

    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = [a for a, _, _ in terms]
    beta_labels = [b for _, b, _ in terms]

    energies = []
    for s in SCALES:
        estimator = build_scaled_estimator(s)
        alpha_mats = ef.build_noisy_matrices(u_vecs, alpha_labels, estimator, K)
        beta_mats = ef.build_noisy_matrices(v_vecs, beta_labels, estimator, K)
        E = ef.ef_energy_from_noisy_matrices(
            terms, lambdas, alpha_mats, beta_mats, enuc, K
        )
        energies.append(E)
        err_kcal = abs(E - exact) * HARTREE_TO_KCAL_MOL
        print(f"  scale {s}: E = {E:.6f} Ha, error = {err_kcal:.2f} kcal/mol")

    lin_fit = np.polyfit(SCALES, energies, 1)
    quad_fit = np.polyfit(SCALES, energies, 2)
    E_lin0 = float(np.polyval(lin_fit, 0))
    E_quad0 = float(np.polyval(quad_fit, 0))

    no_mitigation_err = abs(energies[0] - exact) * HARTREE_TO_KCAL_MOL
    zne_linear_err = abs(E_lin0 - exact) * HARTREE_TO_KCAL_MOL
    zne_quad_err = abs(E_quad0 - exact) * HARTREE_TO_KCAL_MOL

    print(f"\n  NO mitigation (scale=1)   : {no_mitigation_err:.2f} kcal/mol")
    print(f"  ZNE linear extrapolation  : {zne_linear_err:.2f} kcal/mol")
    print(f"  ZNE quadratic extrapolation: {zne_quad_err:.2f} kcal/mol")

    results = {
        "molecule": "H4 fragment (atoms 0-3, d=1.0 Ang)",
        "K": K,
        "exact_ha": round(exact, 6),
        "scales": SCALES,
        "energies_ha": [round(float(x), 6) for x in energies],
        "zne_linear_fit_coeffs": [float(x) for x in lin_fit],
        "zne_quadratic_fit_coeffs": [float(x) for x in quad_fit],
        "no_mitigation_energy_ha": round(float(energies[0]), 6),
        "no_mitigation_err_kcal": round(no_mitigation_err, 4),
        "zne_linear_energy_ha": round(E_lin0, 6),
        "zne_linear_err_kcal": round(zne_linear_err, 4),
        "zne_quadratic_energy_ha": round(E_quad0, 6),
        "zne_quadratic_err_kcal": round(zne_quad_err, 4),
    }

    out = os.path.join(os.path.dirname(__file__), "entanglement_forging_zne_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
