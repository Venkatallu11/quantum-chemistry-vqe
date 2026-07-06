#!/usr/bin/env python3
"""
hardware_mitigation_test.py — does error mitigation actually rescue an
8-qubit fragment on REAL hardware?

Tests block1 (H4, atoms 0-3, d=1.0 A) at FIXED, simulator-optimized
ExcitationPreserving(reps=2) parameters -- the same particle-conserving
ansatz that got to error 0.004749 Ha on the noiseless simulator in
fragment_ansatz_test.py (UCCSD converged closer, but its circuit is ~13x
deeper/noisier -- ExcitationPreserving is the realistic hardware
candidate). The SAME fixed circuit is measured on real hardware three
times, varying ONLY the resilience_level:

  (a) resilience_level=0  -- raw, no mitigation
  (b) resilience_level=1  -- basic mitigation (readout twirling)
  (c) resilience_level=2  -- ZNE (zero-noise extrapolation)

Honest framing: mitigation can at best pull the hardware result toward
the SIMULATOR-ANSATZ target, not the true exact ground state -- the gap
between the ansatz target and exact is baked into the ansatz's limited
expressivity, and no amount of error mitigation on real qubits can fix
that. This script reports error against BOTH references, so that
distinction is never blurred.

Do NOT modify ../quantum-hardware-mcp. This submits 3 real hardware jobs.

Run:
    python vqe/hardware_mitigation_test.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from hardware_fragmentation import _build_fragment_qop, hchain
from fragment_ansatz_test import exact_energy_no_penalty, _hf_bits

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import ExcitationPreserving
from qiskit.primitives import StatevectorEstimator
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from dotenv import load_dotenv
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "..", ".env"))

from mcp_backend import best_qubits_for  # read-only MCP connector, already built

HARTREE_TO_KCAL_MOL = 627.5094740631
CHEM_ACCURACY_TARGET_HA = 0.002


def optimize_ansatz_local(qop_bare, enuc, exact, n_qubits, nelec, reps=2, n_restarts=5, maxiter=500):
    """
    Same ExcitationPreserving + L-BFGS-B optimization used in
    fragment_ansatz_test.py -- Hartree-Fock state prep followed by the
    particle-conserving ansatz, optimized on the noiseless local simulator.
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
        print(f"    restart {r + 1}/{n_restarts}: E = {res.fun:.6f} Ha "
              f"(best {best_energy:.6f} Ha, error {error:.6f} Ha)")
        if error < CHEM_ACCURACY_TARGET_HA:
            break

    qc = make_circuit(best_params)
    return best_params, best_energy, abs(best_energy - exact), qc


def main():
    d = 1.0
    indices = [0, 1, 2, 3]  # block1
    nelec = len(indices)

    geom = hchain(indices, d)
    qop_bare, _qop_penalized, enuc = _build_fragment_qop(geom, nelec)
    n_qubits = qop_bare.num_qubits
    exact = exact_energy_no_penalty(qop_bare, enuc)

    print("\n" + "=" * 70)
    print("  Does error mitigation rescue an 8-qubit fragment on hardware?")
    print("=" * 70)
    print(f"  block1: {n_qubits} qubits, {nelec} electrons")
    print(f"  exact ground state = {exact:.6f} Ha\n")

    print("Step 1: optimizing ExcitationPreserving(reps=2) on the statevector simulator...")
    params, sim_energy, sim_error, qc = optimize_ansatz_local(qop_bare, enuc, exact, n_qubits, nelec)
    print(f"\n  Simulator-ansatz target: E = {sim_energy:.6f} Ha  (error vs exact = {sim_error:.6f} Ha)")
    if sim_error > 0.006:
        print(f"  NOTE: this run's simulator-ansatz error ({sim_error:.6f} Ha) is worse than the "
              f"~0.0047 Ha seen previously -- reporting honestly and continuing anyway.")

    print("\nStep 2: picking the cleanest physical qubits on ibm_marrakesh via the MCP connector...")
    best_qubits = best_qubits_for("ibm_marrakesh", n_qubits)
    if not best_qubits:
        print("ERROR: could not get best qubits via the MCP connector (see message above). Aborting.")
        sys.exit(1)
    initial_layout = [q["qubit"] for q in best_qubits]
    print(f"  Using physical qubits (cleanest first): {initial_layout}")

    print("\nStep 3: connecting to IBM Quantum and transpiling onto ibm_marrakesh...")
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        sys.exit(1)
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend("ibm_marrakesh")

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1, initial_layout=initial_layout)
    isa_circuit = pm.run(qc)
    isa_observable = qop_bare.apply_layout(isa_circuit.layout)

    print("\nStep 4: measuring on real hardware at THREE resilience levels (3 jobs total)...")
    hw_results = {}
    for level, label in ((0, "raw"), (1, "level1"), (2, "zne")):
        print(f"\n  Submitting resilience_level={level} ({label})...")
        estimator = EstimatorV2(mode=backend)
        estimator.options.resilience_level = level
        estimator.options.default_shots = 4096
        job = estimator.run([(isa_circuit, isa_observable)])
        job_id = job.job_id()
        print(f"    job_id={job_id} -- waiting for result...")
        result = job.result()
        electronic_energy = float(result[0].data.evs)
        hw_energy = electronic_energy + enuc
        hw_results[label] = {"resilience_level": level, "energy_ha": hw_energy, "job_id": job_id}
        err_exact = abs(hw_energy - exact)
        err_target = abs(hw_energy - sim_energy)
        print(f"    {label}: E = {hw_energy:.6f} Ha  err_vs_exact={err_exact:.6f} Ha  "
              f"err_vs_target={err_target:.6f} Ha")

    # --- final table --------------------------------------------------------
    print("\n" + "=" * 92)
    print("  RESULTS")
    print("=" * 92)
    print(f"  Exact ground state       = {exact:.6f} Ha")
    print(f"  Simulator-ansatz target  = {sim_energy:.6f} Ha  (error vs exact: {sim_error:.6f} Ha)\n")
    print(f"  {'':10s} | {'energy (Ha)':>12s} | {'err/exact (Ha)':>15s} | {'err/exact (kcal)':>17s} | "
          f"{'err/target (Ha)':>16s} | {'err/target (kcal)':>18s}")
    for label in ("raw", "level1", "zne"):
        r = hw_results[label]
        err_exact_ha = abs(r["energy_ha"] - exact)
        err_exact_kcal = err_exact_ha * HARTREE_TO_KCAL_MOL
        err_target_ha = abs(r["energy_ha"] - sim_energy)
        err_target_kcal = err_target_ha * HARTREE_TO_KCAL_MOL
        print(f"  {label:10s} | {r['energy_ha']:12.6f} | {err_exact_ha:15.6f} | {err_exact_kcal:17.3f} | "
              f"{err_target_ha:16.6f} | {err_target_kcal:18.3f}")

    err_target_raw = abs(hw_results["raw"]["energy_ha"] - sim_energy)
    err_target_zne = abs(hw_results["zne"]["energy_ha"] - sim_energy)
    zne_beats_raw = err_target_zne < err_target_raw

    print()
    print("  HONEST INTERPRETATION:")
    print(f"  Mitigation can at best approach the ansatz target ({sim_energy:.6f} Ha), NOT the exact")
    print(f"  ground state ({exact:.6f} Ha) -- that {sim_error:.6f} Ha gap is baked into the ansatz's")
    print(f"  limited expressivity, and no amount of hardware error mitigation can close it.")
    print(f"  ZNE (resilience_level=2) {'DID' if zne_beats_raw else 'did NOT'} meaningfully beat raw "
          f"(resilience_level=0) at approaching the ansatz target")
    print(f"  ({err_target_zne:.6f} Ha vs {err_target_raw:.6f} Ha error from target).")
    print("=" * 92)

    results = {
        "n_qubits": n_qubits,
        "nelec": nelec,
        "exact_energy_ha": round(exact, 6),
        "simulator_target_ha": round(sim_energy, 6),
        "simulator_target_error_vs_exact_ha": round(sim_error, 6),
        "physical_qubits_used": initial_layout,
        "backend": backend.name,
        "hardware_results": {
            label: {
                "resilience_level": r["resilience_level"],
                "job_id": r["job_id"],
                "energy_ha": round(r["energy_ha"], 6),
                "error_vs_exact_ha": round(abs(r["energy_ha"] - exact), 6),
                "error_vs_exact_kcal_mol": round(abs(r["energy_ha"] - exact) * HARTREE_TO_KCAL_MOL, 4),
                "error_vs_target_ha": round(abs(r["energy_ha"] - sim_energy), 6),
                "error_vs_target_kcal_mol": round(abs(r["energy_ha"] - sim_energy) * HARTREE_TO_KCAL_MOL, 4),
            }
            for label, r in hw_results.items()
        },
        "zne_beats_raw_vs_target": bool(zne_beats_raw),
    }
    out_path = os.path.join(os.path.dirname(__file__), "hardware_mitigation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
