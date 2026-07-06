#!/usr/bin/env python3
"""
h2_hardware_full.py  —  FULLY real-hardware VQE for H2 (no simulator, ever)
============================================================================
Every single energy number in this script comes from a real IBM quantum
computer. There is no simulator anywhere in the optimization loop — not
even to warm-start the optimizer. This is slower and noisier than the
simulator version (h2_vqe.py), but it's the real deal: actual qubits,
actual noise, actual measurement statistics.

How to run
----------
  python vqe/h2_hardware_full.py
"""

import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# Reuse the building blocks from the simulator script — same ansatz, same
# Hamiltonian, same exact-answer function. We do NOT reuse anything that
# runs on a simulator.
sys.path.insert(0, os.path.dirname(__file__))
from h2_vqe import make_ansatz, H2_HAMILTONIAN, exact_ground_state_energy

from dotenv import load_dotenv
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "..", ".env"))

CHEMICAL_ACCURACY_HA = 0.0016   # 1 kcal/mol, in Hartree
HARTREE_TO_KCAL_MOL = 627.5094740631


def main():
    print("\n" + "=" * 60)
    print("  VQE for H2 — 100% Real IBM Hardware")
    print("  (No simulator anywhere in this run)")
    print("=" * 60 + "\n")

    # --- Connect to IBM Quantum ---------------------------------------
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        sys.exit(1)

    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)

    # Pick the least-busy operational backend instead of hardcoding one --
    # a device that was quiet last time can be swamped now (ibm_kingston
    # hit 245 pending jobs and sat stuck in queue for hours).
    print("  Finding least-busy IBM quantum computer...")
    candidates = []
    for b in service.backends():
        try:
            s = b.status()
            if s.operational:
                candidates.append((b, s.pending_jobs))
        except Exception:
            pass
    backend, queue_depth = min(candidates, key=lambda x: x[1])
    print(f"  Selected: {backend.name} ({queue_depth} jobs in queue)")

    n_qubits = backend.num_qubits
    print(f"  Backend: {backend.name}  ({n_qubits} qubits)\n")

    exact = exact_ground_state_energy()
    print(f"  Exact ground state energy (classical reference): {exact:.6f} Ha\n")

    # --- Build the ansatz ONCE, using a symbolic angle -----------------
    # theta_param is a placeholder, not a number yet. This lets us
    # transpile the circuit exactly once, then just plug in different
    # numbers for theta on every hardware job (no re-transpiling).
    theta_param = Parameter("theta")
    ansatz = make_ansatz(theta_param)

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_ansatz = pm.run(ansatz)

    # --- Build the WHOLE Hamiltonian as ONE observable ------------------
    # Instead of measuring each Pauli term (II, IZ, ZI, ZZ, XX, YY) in a
    # separate job, we combine them into a single SparsePauliOp. The
    # estimator then hands back one number that IS the total energy —
    # one hardware job per optimizer step, not six.
    # Build it at the ORIGINAL 2-qubit width — apply_layout() below is
    # what pads/maps it out to the backend's full qubit count.
    hamiltonian_op = SparsePauliOp.from_list(list(H2_HAMILTONIAN.items()))

    # The transpiler may have moved our 2 logical qubits onto different
    # physical qubits on the chip. apply_layout() re-maps the observable
    # (and pads it out to all n_qubits) to match wherever the transpiled
    # circuit actually put them.
    isa_observable = hamiltonian_op.apply_layout(isa_ansatz.layout)

    # --- Run in JOB mode (not Session) ----------------------------------
    # Session requires a paid IBM plan ("dedicated" queue reservation).
    # This account is on the open (free) plan, which only allows
    # single-job submissions — each iteration re-queues independently.
    # Still 100% real hardware, zero simulator; just slower wall-clock
    # time since we lose the "back-to-back, no re-queuing" benefit.
    iteration_log = []
    job_count = 0

    estimator = EstimatorV2(mode=backend)
    estimator.options.resilience_level = 1   # built-in error mitigation
    estimator.options.default_shots = 4096    # measurements per job

    def objective(params):
        nonlocal job_count
        theta = float(params[0])

        # One real hardware job: our fixed transpiled circuit +
        # the combined observable + this iteration's theta value.
        pub = (isa_ansatz, isa_observable, [theta])
        job = estimator.run([pub])
        result = job.result()
        job_count += 1

        # Because the observable already has all 6 Hamiltonian terms
        # (with their coefficients) baked in, this single number IS
        # the measured total energy — no extra math needed.
        energy = float(result[0].data.evs)

        n = len(iteration_log) + 1
        iteration_log.append({
            "iteration": n,
            "theta": round(theta, 5),
            "energy_ha": round(energy, 8),
        })
        print(f"  Job {n:3d} | theta = {theta:+.4f} | real-hardware E = {energy:.6f} Ha")
        return energy

    # COBYLA: doesn't need gradients, tolerates noisy measurements —
    # a good fit for real hardware where every evaluation is noisy.
    result = minimize(
        objective,
        x0=[0.1],
        method="COBYLA",
        options={"maxiter": 30, "rhobeg": 0.5},
    )

    final_energy = float(result.fun)
    error_ha = abs(final_energy - exact)
    error_kcal = error_ha * HARTREE_TO_KCAL_MOL

    print("\n" + "=" * 60)
    print("  RESULTS — Real Hardware VQE")
    print("=" * 60)
    print(f"  Exact energy:              {exact:.6f} Ha")
    print(f"  Final real-hardware energy:{final_energy:.6f} Ha")
    print(f"  Difference:                {error_ha:.6f} Ha  ({error_kcal:.3f} kcal/mol)")
    print(f"  Backend:                   {backend.name}")
    print(f"  Total hardware jobs:       {job_count}")
    print("=" * 60 + "\n")

    results = {
        "method": "real_ibm_hardware_full_loop",
        "backend": backend.name,
        "num_backend_qubits": n_qubits,
        "resilience_level": 1,
        "shots_per_job": 4096,
        "optimal_theta": float(result.x[0]),
        "vqe_energy_ha": round(final_energy, 8),
        "exact_energy_ha": round(exact, 8),
        "error_ha": round(error_ha, 8),
        "error_kcal_mol": round(error_kcal, 4),
        "chemical_accuracy_threshold_ha": CHEMICAL_ACCURACY_HA,
        "reached_chemical_accuracy": bool(error_ha < CHEMICAL_ACCURACY_HA),
        "total_hardware_jobs": job_count,
        "history": iteration_log,
    }

    out_path = os.path.join(os.path.dirname(__file__), "h2_hardware_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Full results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
