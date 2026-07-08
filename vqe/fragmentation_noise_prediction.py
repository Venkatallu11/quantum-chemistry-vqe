#!/usr/bin/env python3
"""
fragmentation_noise_prediction.py — predict how covalent-bond
fragmentation (H6 chain, molecular tailoring) would perform on
Quantinuum hardware BEFORE spending any Azure credits, using a LOCAL
qiskit-aer depolarizing noise model.

Compares two noise scenarios against the noiseless ansatz target and the
true exact energy:
  - Quantinuum-like: 0.2% 2-qubit gate error, 0.005% 1-qubit gate error
    (trapped-ion hardware, typically much better fidelity than
    superconducting qubits)
  - IBM-like: 3% 2-qubit gate error, 0.1% 1-qubit gate error
    (roughly matches what this repo has actually observed submitting to
    real IBM hardware)

This is ALL LOCAL simulation -- no Azure, no IBM, no network calls, no
credits spent. It exists purely to inform whether spending real Azure
Quantum credit on Quantinuum is even likely to be worth it for this
fragment size, before committing real money to find out.

Run:
    python vqe/fragmentation_noise_prediction.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from hardware_fragmentation import _build_fragment_qop, hchain
from fragment_ansatz_test import exact_energy_no_penalty, _hf_bits

from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import ExcitationPreserving
from qiskit.primitives import StatevectorEstimator
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error

HARTREE_TO_KCAL_MOL = 627.5094740631
BASIS_GATES = ["u3", "cx"]
SHOTS = 4096

# Depolarizing error rates: (two-qubit gate error, one-qubit gate error)
NOISE_LEVELS = {
    "quantinuum": {"two_q": 0.002, "one_q": 0.00005},
    "ibm":        {"two_q": 0.03,  "one_q": 0.001},
}

FULL_EXACT_H6_HA = -3.236066


def build_noise_model(one_q_error, two_q_error):
    """A simple uniform depolarizing noise model on a fixed 1q+2q basis."""
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(two_q_error, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(one_q_error, 1), "u3")
    return nm


def optimize_fragment_locally(qop_bare, enuc, exact, n_qubits, nelec, reps=2, n_restarts=5, maxiter=500):
    """
    Hartree-Fock state prep + ExcitationPreserving(reps=2), optimized with
    L-BFGS-B on the NOISELESS statevector simulator -- identical approach
    to hardware_covalent.py / fragment_ansatz_test.py, so these parameters
    are directly comparable to the real-hardware attempts already made.
    """
    bits = _hf_bits(n_qubits, nelec)
    base_ansatz = ExcitationPreserving(n_qubits, reps=reps)
    n_params = base_ansatz.num_parameters
    estimator = StatevectorEstimator()

    def make_circuit(params):
        qc = QuantumCircuit(n_qubits)
        for i, b in enumerate(bits):
            if b:
                qc.x(i)
        qc.compose(base_ansatz.assign_parameters(params), inplace=True)
        return qc

    def objective(params):
        qc = make_circuit(params)
        result = estimator.run([(qc, qop_bare)]).result()
        return float(result[0].data.evs) + enuc

    rng = np.random.default_rng(42)
    best_params, best_energy = None, None
    for r in range(n_restarts):
        x0 = np.zeros(n_params) if r == 0 else rng.uniform(-np.pi, np.pi, size=n_params)
        res = minimize(
            objective, x0=x0, method="L-BFGS-B",
            options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-10},
        )
        if best_energy is None or res.fun < best_energy:
            best_energy, best_params = float(res.fun), res.x
        error = abs(best_energy - exact)
        print(f"      restart {r + 1}/{n_restarts}: E = {res.fun:.6f} Ha "
              f"(best {best_energy:.6f} Ha, error {error:.6f} Ha)")
        if error < 0.0005:
            break

    qc = make_circuit(best_params)
    return best_params, best_energy, abs(best_energy - exact), qc


def measure_with_noise(qc, qop_bare, enuc, one_q_error, two_q_error):
    """
    Transpile to a fixed 1q+2q basis and measure the bare electronic
    Hamiltonian's expectation value under a depolarizing noise model,
    entirely locally via qiskit-aer -- no network, no real hardware.
    """
    nm = build_noise_model(one_q_error, two_q_error)
    qc_t = transpile(qc, basis_gates=BASIS_GATES, optimization_level=1)

    estimator = AerEstimatorV2(options={
        "backend_options": {"noise_model": nm},
        "run_options": {"shots": SHOTS},
    })
    result = estimator.run([(qc_t, qop_bare)]).result()
    electronic_energy = float(result[0].data.evs)
    return electronic_energy + enuc


def process_fragment(name, indices, d):
    """Build one fragment's Hamiltonian, optimize noiselessly, then measure under both noise models."""
    geom = hchain(indices, d)
    nelec = len(indices)
    print(f"\n--- Fragment '{name}' (atoms {indices}, {nelec} electrons) ---")

    qop_bare, _qop_penalized, enuc = _build_fragment_qop(geom, nelec)
    n_qubits = qop_bare.num_qubits
    exact = exact_energy_no_penalty(qop_bare, enuc)
    print(f"  qubits = {n_qubits}, exact ground state = {exact:.6f} Ha")

    print("  Optimizing ExcitationPreserving(reps=2) on the noiseless simulator (L-BFGS-B)...")
    params, noiseless_energy, noiseless_error, qc = optimize_fragment_locally(
        qop_bare, enuc, exact, n_qubits, nelec
    )
    print(f"  -> noiseless best: E = {noiseless_energy:.6f} Ha, error = {noiseless_error:.6f} Ha")

    noise_energies = {}
    for label, rates in NOISE_LEVELS.items():
        print(f"  Measuring under {label}-like noise "
              f"(2q={rates['two_q'] * 100:.3f}%, 1q={rates['one_q'] * 100:.4f}%, shots={SHOTS})...")
        e = measure_with_noise(qc, qop_bare, enuc, rates["one_q"], rates["two_q"])
        err = abs(e - exact)
        print(f"    -> E = {e:.6f} Ha, error vs exact = {err:.6f} Ha")
        noise_energies[label] = e

    return {
        "name": name, "indices": indices, "nelec": nelec, "n_qubits": n_qubits,
        "exact_energy_ha": exact,
        "noiseless_energy_ha": noiseless_energy, "noiseless_error_ha": noiseless_error,
        "noise_energies_ha": noise_energies,
    }


def main():
    d = 1.0  # Angstrom, same bonded H6 chain as covalent_fragment.py / hardware_covalent.py

    print("\n" + "=" * 70)
    print("  Fragmentation noise prediction -- H6 chain")
    print("  (LOCAL qiskit-aer simulation only -- no Azure, no IBM, no cost)")
    print("=" * 70)

    fragments = {
        "block1":  process_fragment("block1 (atoms 0-3)",  [0, 1, 2, 3], d),
        "block2":  process_fragment("block2 (atoms 2-5)",  [2, 3, 4, 5], d),
        "overlap": process_fragment("overlap (atoms 2-3)", [2, 3],       d),
    }

    # --- Reassemble via molecular tailoring, once per noise scenario -----
    e_noiseless_h6 = (fragments["block1"]["noiseless_energy_ha"]
                       + fragments["block2"]["noiseless_energy_ha"]
                       - fragments["overlap"]["noiseless_energy_ha"])
    err_noiseless = abs(e_noiseless_h6 - FULL_EXACT_H6_HA)

    reassembled = {}
    for label in NOISE_LEVELS:
        e_h6 = (fragments["block1"]["noise_energies_ha"][label]
                + fragments["block2"]["noise_energies_ha"][label]
                - fragments["overlap"]["noise_energies_ha"][label])
        reassembled[label] = e_h6

    # --- Tables --------------------------------------------------------------
    print("\n" + "=" * 95)
    print("  FRAGMENT SUMMARY")
    print("=" * 95)
    print(f"  {'fragment':22s} | {'qubits':>6s} | {'exact (Ha)':>11s} | {'noiseless (Ha)':>14s} | "
          f"{'quantinuum (Ha)':>16s} | {'ibm-like (Ha)':>14s}")
    for frag in fragments.values():
        print(f"  {frag['name']:22s} | {frag['n_qubits']:6d} | {frag['exact_energy_ha']:11.6f} | "
              f"{frag['noiseless_energy_ha']:14.6f} | "
              f"{frag['noise_energies_ha']['quantinuum']:16.6f} | "
              f"{frag['noise_energies_ha']['ibm']:14.6f}")

    print("\n" + "=" * 95)
    print("  REASSEMBLED H6 (molecular tailoring)")
    print("=" * 95)
    print(f"  Full exact H6                    = {FULL_EXACT_H6_HA:.6f} Ha")
    print(f"  Noiseless tailoring H6           = {e_noiseless_h6:.6f} Ha  (error {err_noiseless:.6f} Ha)")

    for label, display_name in (("quantinuum", "Quantinuum-like noise"), ("ibm", "IBM-like noise")):
        e_h6 = reassembled[label]
        err_ha = abs(e_h6 - FULL_EXACT_H6_HA)
        err_kcal = err_ha * HARTREE_TO_KCAL_MOL
        print(f"  {display_name:28s} H6 = {e_h6:.6f} Ha  (error {err_ha:.6f} Ha, {err_kcal:.3f} kcal/mol)")
    print("=" * 95 + "\n")

    # --- Save results ----------------------------------------------------------
    results = {
        "full_exact_h6_ha": FULL_EXACT_H6_HA,
        "noiseless_tailoring_h6_ha": round(e_noiseless_h6, 6),
        "noiseless_tailoring_error_ha": round(err_noiseless, 6),
        "shots": SHOTS,
        "noise_levels": NOISE_LEVELS,
        "fragments": {
            key: {
                "name": frag["name"], "indices": frag["indices"], "n_qubits": frag["n_qubits"],
                "exact_energy_ha": round(frag["exact_energy_ha"], 6),
                "noiseless_energy_ha": round(frag["noiseless_energy_ha"], 6),
                "noiseless_error_ha": round(frag["noiseless_error_ha"], 6),
                "noise_energies_ha": {k: round(v, 6) for k, v in frag["noise_energies_ha"].items()},
            }
            for key, frag in fragments.items()
        },
        "reassembled_h6_ha": {k: round(v, 6) for k, v in reassembled.items()},
        "reassembled_error_ha": {k: round(abs(v - FULL_EXACT_H6_HA), 6) for k, v in reassembled.items()},
        "reassembled_error_kcal_mol": {
            k: round(abs(v - FULL_EXACT_H6_HA) * HARTREE_TO_KCAL_MOL, 4) for k, v in reassembled.items()
        },
    }
    out_path = os.path.join(os.path.dirname(__file__), "fragmentation_noise_prediction.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
