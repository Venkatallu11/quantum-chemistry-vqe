#!/usr/bin/env python3
"""
h2_hardware_direct.py — measure H2 energy on real hardware DIRECTLY,
without going through Lokesh's MCP server, as a comparison point against
mcp_energy.py's MCP-routed result.

Same molecule, same verified Hamiltonian (-1.137284 Ha exact), same
ansatz, same resilience settings and shot count as the MCP path -- the
ONLY difference is how the job gets to IBM:

  mcp_energy.py:        our code -> server.estimate_expectation() ->
                         qiskit-ibm-runtime -> IBM
  h2_hardware_direct.py: our code -> qiskit-ibm-runtime -> IBM  (this file)

If both give similar energies (within noise), that's evidence the MCP
layer isn't adding its own distortion -- it's just a thin wrapper around
the same primitives.

Run (submits ONE real hardware job):
    python vqe/h2_hardware_direct.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import StatevectorEstimator
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

sys.path.insert(0, os.path.dirname(__file__))
from h2_vqe import make_ansatz, H2_HAMILTONIAN, exact_ground_state_energy, N_ANSATZ_PARAMS

from dotenv import load_dotenv
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "..", ".env"))

HARTREE_TO_KCAL_MOL = 627.5094740631


def _optimize_params_on_simulator():
    """
    Classical part of VQE: find the best EfficientSU2 parameters on a
    noiseless local simulator (free, instant) -- identical approach and
    tuning to mcp_energy.py, so both scripts start from the same optimized
    circuit and only differ in how the final measurement reaches hardware.
    """
    estimator = StatevectorEstimator()

    def objective(params):
        qc = make_ansatz(params)
        energy = 0.0
        for pauli_str, coeff in H2_HAMILTONIAN.items():
            if set(pauli_str) == {"I"}:
                energy += coeff
                continue
            observable = SparsePauliOp(pauli_str)
            result = estimator.run([(qc, observable)]).result()
            energy += coeff * float(result[0].data.evs)
        return energy

    result = minimize(
        objective, x0=np.full(N_ANSATZ_PARAMS, 0.1), method="COBYLA",
        options={"maxiter": 3000, "rhobeg": 0.5, "tol": 1e-10},
    )
    return result.x


def _pick_least_busy_backend(service):
    """
    Check every accessible IBM backend's queue depth directly via
    qiskit-ibm-runtime (no MCP server involved) and return the least-busy
    operational one -- same policy as the MCP path, applied natively.
    """
    candidates = []
    for b in service.backends():
        try:
            s = b.status()
            if s.operational:
                candidates.append((b, s.pending_jobs))
        except Exception:
            pass
    backend, queue_depth = min(candidates, key=lambda x: x[1])
    print(f"  Least-busy device right now: {backend.name} ({queue_depth} pending jobs)\n")
    return backend


def main():
    print("\n" + "=" * 60)
    print("  H2 energy, measured on real hardware DIRECTLY (no MCP)")
    print("  (verified Hamiltonian, exact ground state = -1.137284 Ha)")
    print("=" * 60 + "\n")

    # Classical step: optimize the ansatz on a free local simulator.
    params = _optimize_params_on_simulator()
    print(f"  Optimized params (statevector simulator) = "
          f"{np.round(params, 5).tolist()}\n")

    exact = exact_ground_state_energy()

    # --- Connect directly to IBM Quantum -------------------------------
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        sys.exit(1)

    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = _pick_least_busy_backend(service)
    n_qubits = backend.num_qubits

    # --- Build and transpile the optimized circuit ----------------------
    qc = make_ansatz(params)
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)

    # --- Build the WHOLE Hamiltonian as ONE observable ------------------
    # Built at the original 2-qubit width; apply_layout() pads it out to
    # the backend's full qubit count and maps it onto wherever the
    # transpiler actually placed our 2 logical qubits.
    hamiltonian_op = SparsePauliOp.from_list(list(H2_HAMILTONIAN.items()))
    isa_observable = hamiltonian_op.apply_layout(isa_circuit.layout)

    # --- ONE real hardware job, same settings as the MCP path -----------
    estimator = EstimatorV2(mode=backend)
    estimator.options.resilience_level = 1
    estimator.options.default_shots = 4096

    pub = (isa_circuit, isa_observable)
    job = estimator.run([pub])
    job_id = job.job_id()
    print(f"  Submitted job {job_id} on {backend.name} (direct, no MCP)")

    result = job.result()
    energy = float(result[0].data.evs)

    diff_ha = abs(energy - exact)
    diff_kcal = diff_ha * HARTREE_TO_KCAL_MOL

    print("\n" + "=" * 60)
    print("  RESULT")
    print("=" * 60)
    print(f"  Exact energy                  = {exact:.6f} Ha")
    print(f"  Hardware energy (direct)      = {energy:.6f} Ha")
    print(f"  Difference                    = {diff_ha:.6f} Ha  ({diff_kcal:.3f} kcal/mol)")
    print(f"  Device                        = {backend.name}")
    print(f"  Job ID                        = {job_id}")
    print("=" * 60 + "\n")

    results = {
        "method": "real_ibm_hardware_direct_no_mcp",
        "hamiltonian": "verified (chem.py engine, parity tapering)",
        "exact_energy_ha": round(exact, 6),
        "hardware_energy_ha": round(energy, 6),
        "error_ha": round(diff_ha, 6),
        "error_kcal_mol": round(diff_kcal, 4),
        "device": backend.name,
        "job_id": job_id,
        "optimal_params": [float(p) for p in params],
    }
    out_path = os.path.join(os.path.dirname(__file__), "h2_hardware_direct_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
