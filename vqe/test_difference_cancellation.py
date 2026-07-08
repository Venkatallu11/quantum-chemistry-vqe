#!/usr/bin/env python3
"""
test_difference_cancellation.py — does the fragmentation DIFFERENCE
structure (E(H6) = E(atoms 0-3) + E(atoms 2-5) - E(atoms 2-3)) cancel
SYSTEMATIC (coherent) noise better than RANDOM (incoherent) noise?

Hypothesis: block1, block2, and overlap share real structural overlap --
"overlap" (atoms 2-3) is literally the shared H-H bond present in both
blocks, at the same bond length (d=1.0 A). If a noise source is a fixed,
deterministic bias (e.g. a miscalibrated entangling-gate angle -- the
same physical over-rotation on structurally similar gates in each
fragment), it may partially CANCEL in the inclusion-exclusion sum, since
the same bias appears with the same sign pattern in the shared structure.
Purely random (incoherent) noise has no such correlation to exploit -- it
just adds independent noise to each fragment separately, with nothing to
cancel against.

Two noise models, calibrated to the SAME average 2-qubit gate infidelity
(1%), applied ONLY to 2-qubit (CX) gates:
  (a) COHERENT:   a fixed small RZZ over-rotation after every CX -- a
                   deterministic unitary error (systematic bias), angle
                   calibrated so its average gate infidelity equals the
                   target via F_avg = (d*F_pro + 1)/(d+1).
  (b) INCOHERENT: depolarizing_error at the same nominal rate (random,
                   uncorrelated shot-to-shot).

All local, free -- no Azure, no IBM, no network calls.

Run:
    python vqe/test_difference_cancellation.py
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
from qiskit.circuit.library import ExcitationPreserving, RZZGate
from qiskit.primitives import StatevectorEstimator
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error, coherent_unitary_error

HARTREE_TO_KCAL_MOL = 627.5094740631
BASIS_GATES = ["u3", "cx"]
SHOTS = 4096
TARGET_TWO_Q_INFIDELITY = 0.01
FULL_EXACT_H6_HA = -3.236066


def coherent_zz_angle_for_infidelity(target_infidelity, dim=4):
    """
    Angle theta for a residual RZZ(theta) over-rotation whose average gate
    infidelity equals target_infidelity, using the standard relation
    between process and average fidelity:
        F_avg = (d*F_pro + 1) / (d + 1),   F_pro = |Tr(U)/d|^2 = cos^2(theta/2)
    """
    F_avg = 1 - target_infidelity
    F_pro = (F_avg * (dim + 1) - 1) / dim
    return 2 * np.arccos(np.sqrt(F_pro))


def build_coherent_noise_model():
    """Fixed ZZ over-rotation after every CX -- a deterministic (systematic) unitary error."""
    theta = coherent_zz_angle_for_infidelity(TARGET_TWO_Q_INFIDELITY)
    over_rotation = RZZGate(theta).to_matrix()
    err = coherent_unitary_error(over_rotation)
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(err, "cx")
    return nm, theta


def build_incoherent_noise_model():
    """Depolarizing (random, uncorrelated) error on every CX, same nominal rate."""
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(TARGET_TWO_Q_INFIDELITY, 2), "cx")
    return nm


def optimize_fragment_locally(qop_bare, enuc, exact, n_qubits, nelec, reps=2, n_restarts=5, maxiter=500):
    """Same approach as hardware_covalent.py / fragment_ansatz_test.py / fragmentation_noise_prediction.py."""
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


def measure_with_noise_model(qc, qop_bare, enuc, noise_model):
    qc_t = transpile(qc, basis_gates=BASIS_GATES, optimization_level=1)
    estimator = AerEstimatorV2(options={
        "backend_options": {"noise_model": noise_model},
        "run_options": {"shots": SHOTS},
    })
    result = estimator.run([(qc_t, qop_bare)]).result()
    electronic_energy = float(result[0].data.evs)
    return electronic_energy + enuc


def process_fragment(name, indices, d, coherent_nm, incoherent_nm):
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

    print("  Measuring under COHERENT noise (fixed ZZ over-rotation)...")
    e_coherent = measure_with_noise_model(qc, qop_bare, enuc, coherent_nm)
    err_coherent = abs(e_coherent - exact)
    print(f"    -> E = {e_coherent:.6f} Ha, error vs exact = {err_coherent:.6f} Ha")

    print("  Measuring under INCOHERENT noise (depolarizing)...")
    e_incoherent = measure_with_noise_model(qc, qop_bare, enuc, incoherent_nm)
    err_incoherent = abs(e_incoherent - exact)
    print(f"    -> E = {e_incoherent:.6f} Ha, error vs exact = {err_incoherent:.6f} Ha")

    return {
        "name": name, "indices": indices, "nelec": nelec, "n_qubits": n_qubits,
        "exact_energy_ha": exact,
        "noiseless_energy_ha": noiseless_energy, "noiseless_error_ha": noiseless_error,
        "coherent_energy_ha": e_coherent, "coherent_error_ha": err_coherent,
        "incoherent_energy_ha": e_incoherent, "incoherent_error_ha": err_incoherent,
    }


def main():
    d = 1.0  # Angstrom, same bonded H6 chain as covalent_fragment.py / hardware_covalent.py

    print("\n" + "=" * 70)
    print("  Difference-cancellation test -- coherent vs incoherent noise")
    print("  (LOCAL qiskit-aer simulation only -- no Azure, no IBM, no cost)")
    print("=" * 70)

    coherent_nm, theta = build_coherent_noise_model()
    incoherent_nm = build_incoherent_noise_model()
    print(f"\n  Coherent noise:   fixed RZZ({theta:.6f} rad) over-rotation after every CX "
          f"(calibrated to {TARGET_TWO_Q_INFIDELITY * 100:.1f}% avg infidelity)")
    print(f"  Incoherent noise: depolarizing_error({TARGET_TWO_Q_INFIDELITY}, 2) after every CX")

    fragments = {
        "block1":  process_fragment("block1 (atoms 0-3)",  [0, 1, 2, 3], d, coherent_nm, incoherent_nm),
        "block2":  process_fragment("block2 (atoms 2-5)",  [2, 3, 4, 5], d, coherent_nm, incoherent_nm),
        "overlap": process_fragment("overlap (atoms 2-3)", [2, 3],       d, coherent_nm, incoherent_nm),
    }

    # --- Reassemble via molecular tailoring, once per noise type -----------
    e_coherent_h6 = (fragments["block1"]["coherent_energy_ha"]
                      + fragments["block2"]["coherent_energy_ha"]
                      - fragments["overlap"]["coherent_energy_ha"])
    e_incoherent_h6 = (fragments["block1"]["incoherent_energy_ha"]
                        + fragments["block2"]["incoherent_energy_ha"]
                        - fragments["overlap"]["incoherent_energy_ha"])

    err_coherent_h6 = abs(e_coherent_h6 - FULL_EXACT_H6_HA)
    err_incoherent_h6 = abs(e_incoherent_h6 - FULL_EXACT_H6_HA)

    avg_frag_err_coherent = float(np.mean([fragments[k]["coherent_error_ha"] for k in fragments]))
    avg_frag_err_incoherent = float(np.mean([fragments[k]["incoherent_error_ha"] for k in fragments]))

    # --- Tables --------------------------------------------------------------
    print("\n" + "=" * 100)
    print("  PER-FRAGMENT ERRORS (vs each fragment's own exact energy)")
    print("=" * 100)
    print(f"  {'fragment':22s} | {'qubits':>6s} | {'noiseless err':>14s} | {'coherent err':>13s} | {'incoherent err':>15s}")
    for frag in fragments.values():
        print(f"  {frag['name']:22s} | {frag['n_qubits']:6d} | {frag['noiseless_error_ha']:14.6f} | "
              f"{frag['coherent_error_ha']:13.6f} | {frag['incoherent_error_ha']:15.6f}")
    print(f"  {'AVERAGE':22s} | {'':>6s} | {'':>14s} | {avg_frag_err_coherent:13.6f} | {avg_frag_err_incoherent:15.6f}")

    print("\n" + "=" * 100)
    print("  REASSEMBLED H6 (molecular tailoring)")
    print("=" * 100)
    print(f"  Exact H6                            = {FULL_EXACT_H6_HA:.6f} Ha")
    print(f"  Reassembled under COHERENT noise    = {e_coherent_h6:.6f} Ha  (error {err_coherent_h6:.6f} Ha, "
          f"{err_coherent_h6 * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print(f"  Reassembled under INCOHERENT noise  = {e_incoherent_h6:.6f} Ha  (error {err_incoherent_h6:.6f} Ha, "
          f"{err_incoherent_h6 * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")

    print("\n" + "=" * 100)
    print("  HYPOTHESIS CHECK: does the difference structure cancel COHERENT noise MORE than INCOHERENT?")
    print("=" * 100)
    coherent_cancels = err_coherent_h6 < avg_frag_err_coherent
    incoherent_cancels = err_incoherent_h6 < avg_frag_err_incoherent
    print(f"  Coherent:   reassembled error ({err_coherent_h6:.6f} Ha) vs avg fragment error "
          f"({avg_frag_err_coherent:.6f} Ha) -> "
          f"{'CANCELS (smaller)' if coherent_cancels else 'DOES NOT CANCEL (not smaller)'}")
    print(f"  Incoherent: reassembled error ({err_incoherent_h6:.6f} Ha) vs avg fragment error "
          f"({avg_frag_err_incoherent:.6f} Ha) -> "
          f"{'CANCELS (smaller)' if incoherent_cancels else 'DOES NOT CANCEL (not smaller)'}")

    coherent_cancel_ratio = (1 - err_coherent_h6 / avg_frag_err_coherent) if avg_frag_err_coherent > 0 else 0.0
    incoherent_cancel_ratio = (1 - err_incoherent_h6 / avg_frag_err_incoherent) if avg_frag_err_incoherent > 0 else 0.0
    print(f"\n  Cancellation ratio (1 - reassembled_err/avg_fragment_err): "
          f"coherent = {coherent_cancel_ratio:.3f}, incoherent = {incoherent_cancel_ratio:.3f}")

    hypothesis_confirmed = coherent_cancel_ratio > incoherent_cancel_ratio
    print(f"\n  Hypothesis (coherent noise cancels MORE than incoherent in the reassembled total): "
          f"{'CONFIRMED' if hypothesis_confirmed else 'NOT CONFIRMED'}")
    print("=" * 100 + "\n")

    # --- Save results -------------------------------------------------------
    results = {
        "full_exact_h6_ha": FULL_EXACT_H6_HA,
        "target_two_qubit_infidelity": TARGET_TWO_Q_INFIDELITY,
        "coherent_rzz_angle_rad": theta,
        "shots": SHOTS,
        "fragments": {
            key: {
                "name": frag["name"], "indices": frag["indices"], "n_qubits": frag["n_qubits"],
                "exact_energy_ha": round(frag["exact_energy_ha"], 6),
                "noiseless_energy_ha": round(frag["noiseless_energy_ha"], 6),
                "noiseless_error_ha": round(frag["noiseless_error_ha"], 6),
                "coherent_energy_ha": round(frag["coherent_energy_ha"], 6),
                "coherent_error_ha": round(frag["coherent_error_ha"], 6),
                "incoherent_energy_ha": round(frag["incoherent_energy_ha"], 6),
                "incoherent_error_ha": round(frag["incoherent_error_ha"], 6),
            }
            for key, frag in fragments.items()
        },
        "average_fragment_error_ha": {
            "coherent": round(avg_frag_err_coherent, 6),
            "incoherent": round(avg_frag_err_incoherent, 6),
        },
        "reassembled_h6_ha": {
            "coherent": round(e_coherent_h6, 6),
            "incoherent": round(e_incoherent_h6, 6),
        },
        "reassembled_error_ha": {
            "coherent": round(err_coherent_h6, 6),
            "incoherent": round(err_incoherent_h6, 6),
        },
        "reassembled_error_kcal_mol": {
            "coherent": round(err_coherent_h6 * HARTREE_TO_KCAL_MOL, 4),
            "incoherent": round(err_incoherent_h6 * HARTREE_TO_KCAL_MOL, 4),
        },
        "cancellation_ratio": {
            "coherent": round(coherent_cancel_ratio, 4),
            "incoherent": round(incoherent_cancel_ratio, 4),
        },
        "hypothesis_confirmed": bool(hypothesis_confirmed),
    }
    out_path = os.path.join(os.path.dirname(__file__), "difference_cancellation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
