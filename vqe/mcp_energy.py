#!/usr/bin/env python3
"""
mcp_energy.py — measure a molecule's energy on real hardware THROUGH
Lokesh's Quantum Hardware MCP server (../quantum-hardware-mcp/server.py),
instead of calling qiskit-ibm-runtime directly the way h2_hardware_full.py
does.

Why go through the MCP server instead of qiskit-ibm-runtime directly:
  Lokesh's server already wraps the Estimator primitive as an MCP tool
  (estimate_expectation) plus job tracking (job_status). This file proves
  that a real chemistry measurement -- not just device metadata -- can be
  made to flow entirely through that MCP layer, the same layer an AI
  assistant like Claude Desktop would use.

The one gotcha (see server.py's estimate_expectation docstring, which is
slightly wrong about this): job_results does NOT work for Estimator jobs
-- it's built for Sampler (counts) jobs only. For Estimator jobs we have
to reach past job_results and pull the raw result straight from IBM via
server._get_service().job(job_id).result().

This module imports Lokesh's server.py directly (the same technique as
mcp_backend.py in this same folder) and never modifies anything in
../quantum-hardware-mcp.

Run the demo (submits ONE real hardware job, no loop):
    python vqe/mcp_energy.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from h2_vqe import make_ansatz, H2_HAMILTONIAN, exact_ground_state_energy, N_ANSATZ_PARAMS

from qiskit import qasm2
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import StatevectorEstimator
from scipy.optimize import minimize

# Lokesh's Quantum Hardware MCP server lives one folder up from this
# project, as a sibling repo -- same path convention as mcp_backend.py.
# We only ever import it (read-only); nothing in that folder is modified.
_MCP_SERVER_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "quantum-hardware-mcp")
)
if _MCP_SERVER_DIR not in sys.path:
    sys.path.insert(0, _MCP_SERVER_DIR)
import server  # Lokesh's server.py


def measure_energy_via_mcp(device_name: str, qasm_string: str,
                            hamiltonian_terms: dict, shots: int = 4096):
    """
    Measure <H> for a prepared quantum state (given as QASM), by routing
    the measurement through Lokesh's MCP server instead of talking to IBM
    directly.

    Args:
        device_name: IBM backend to run on, e.g. "ibm_kingston".
        qasm_string: circuit that prepares the state (no measurements --
                     the Estimator handles that internally).
        hamiltonian_terms: dict like H2_HAMILTONIAN, {pauli_string: coefficient}.
                           One entry must be the all-identity term (e.g. "II"),
                           which contributes a constant offset, not something
                           that needs to be measured.
        shots: shots per observable (this is a REAL hardware job -- this
               is the ONE job submitted, not one per Pauli term).

    Returns:
        (energy, job_id) -- energy in the same units as the Hamiltonian
        coefficients (Hartree, for the chemistry Hamiltonians in this repo).
    """
    # --- Step 1: split off the identity term as a constant ------------
    # The identity term ("II", "III", ...) doesn't need to be measured --
    # its expectation value is always exactly 1, so it just contributes
    # its coefficient directly to the total energy.
    identity_constant = 0.0
    non_identity_terms = []  # [(pauli_str, coefficient), ...], order preserved
    for pauli_str, coeff in hamiltonian_terms.items():
        if set(pauli_str) == {"I"}:
            identity_constant += coeff
        else:
            non_identity_terms.append((pauli_str, coeff))

    # --- Step 2: build the comma-joined observables string --------------
    # estimate_expectation pads each Pauli string with I's internally to
    # match the transpiled (ISA) circuit's full backend qubit count -- we
    # just need to match the circuit's OWN qubit count here.
    observables = ",".join(pauli_str for pauli_str, _coeff in non_identity_terms)

    # --- Step 3: submit ONE real hardware job through the MCP server ----
    submit_raw = server.estimate_expectation(device_name, qasm_string, observables, shots)
    submit_data = json.loads(submit_raw)
    if "error" in submit_data:
        raise RuntimeError(f"estimate_expectation failed: {submit_data['error']}")

    job_id = submit_data["job_id"]
    print(f"  Submitted via MCP: job_id={job_id}  device={device_name}  "
          f"observables=[{observables}]  shots={shots}")

    # --- Step 4: poll job_status (an MCP tool) until DONE ----------------
    while True:
        status_data = json.loads(server.job_status(job_id))
        status = status_data.get("status")
        print(f"  job_status -> {status}")

        if status == "DONE":
            break
        if status in ("ERROR", "CANCELLED"):
            raise RuntimeError(f"Job {job_id} ended with status {status}: {status_data}")

        time.sleep(20)  # ~20s between polls, as instructed

    # --- Step 5: job_results does NOT work for Estimator jobs -- pull the
    # raw result directly from IBM instead, using the same service object
    # server.py's own tools use internally.
    result = server._get_service().job(job_id).result()

    # --- Step 6: reassemble the energy, in the SAME order as `observables` --
    energy = identity_constant
    for i, (_pauli_str, coeff) in enumerate(non_identity_terms):
        value = float(result[i].data.evs)
        energy += coeff * value

    return energy, job_id


def _optimize_params_on_simulator():
    """
    Classical part of VQE: find the best EfficientSU2 parameter vector
    using a noiseless local simulator. This is standard practice -- only
    the FINAL energy measurement at the optimized parameters needs real
    hardware; there's no reason to burn hardware queue time on every
    optimizer step just to find where to point the circuit.

    8 free parameters (EfficientSU2(2, reps=1)) need more iterations than
    a single angle to converge tightly -- 3000 reaches ~1e-6 Ha on this
    Hamiltonian (see h2_vqe.py's run_vqe_local for the same tuning).
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


def _pick_least_busy_device() -> str:
    """
    Check every accessible IBM backend's queue depth (via the MCP server's
    compare_devices tool, read-only) and return the one with the fewest
    pending jobs. We check this fresh every run instead of hardcoding a
    device name -- a device that was quiet last week can be swamped today.
    """
    data = json.loads(server.compare_devices("queue"))
    best = data["devices"][0]
    print(f"  Least-busy device right now: {best['name']} "
          f"({best['pending_jobs']} pending jobs)\n")
    return best["name"]


HARTREE_TO_KCAL_MOL = 627.5094740631


def main():
    print("\n" + "=" * 60)
    print("  H2 energy, measured on real hardware THROUGH the MCP server")
    print("  (verified Hamiltonian, exact ground state = -1.137284 Ha)")
    print("=" * 60 + "\n")

    # Classical step: find the best ansatz parameters on a simulator
    # (free, instant). Only the FINAL measurement needs real hardware.
    params = _optimize_params_on_simulator()
    print(f"  Optimized params (statevector simulator) = "
          f"{np.round(params, 5).tolist()}\n")

    # Build the ansatz circuit at those parameters and export it to QASM2 --
    # this is the exact circuit whose energy we'll measure on hardware.
    qc = make_ansatz(params)
    qasm_string = qasm2.dumps(qc)

    exact = exact_ground_state_energy()

    # Check which device is least busy RIGHT NOW instead of hardcoding one.
    device_name = _pick_least_busy_device()

    # ONE real hardware measurement, submitted through Lokesh's MCP server.
    energy, job_id = measure_energy_via_mcp(device_name, qasm_string, H2_HAMILTONIAN)

    diff_ha = abs(energy - exact)
    diff_kcal = diff_ha * HARTREE_TO_KCAL_MOL

    print("\n" + "=" * 60)
    print("  RESULT")
    print("=" * 60)
    print(f"  Exact energy               = {exact:.6f} Ha")
    print(f"  Hardware energy (via MCP)  = {energy:.6f} Ha")
    print(f"  Difference                 = {diff_ha:.6f} Ha  ({diff_kcal:.3f} kcal/mol)")
    print(f"  Device                     = {device_name}")
    print(f"  Job ID                     = {job_id}")
    print("=" * 60 + "\n")

    results = {
        "method": "real_ibm_hardware_via_mcp_server",
        "hamiltonian": "verified (chem.py engine, parity tapering)",
        "exact_energy_ha": round(exact, 6),
        "hardware_energy_ha": round(energy, 6),
        "error_ha": round(diff_ha, 6),
        "error_kcal_mol": round(diff_kcal, 4),
        "device": device_name,
        "job_id": job_id,
        "optimal_params": [float(p) for p in params],
    }
    out_path = os.path.join(os.path.dirname(__file__), "h2_hardware_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
