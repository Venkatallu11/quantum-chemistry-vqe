#!/usr/bin/env python3
"""
hardware_clean.py — measure H2 on real IBM hardware ENTIRELY through our
own engine (qiskit-ibm-runtime directly). Lokesh's MCP server is used for
exactly ONE thing: picking the best device (best_device(), which wraps
server.compare_devices). The measurement itself never touches the MCP
server, and no specific qubits are forced -- the transpiler auto-places
using its own noise+connectivity-aware layout/routing.

Why auto-place instead of hand-picking qubits: hardware_mitigation_test.py
tried forcing the 8 "cleanest individual" qubits (via best_qubits_for) as
an initial_layout, and it backfired badly -- those qubits were barely
connected on the chip (1 edge out of 8 qubits), forcing a flood of SWAP
gates and producing energies off by >1 Ha. Letting
generate_preset_pass_manager pick its own layout avoids that entirely.

Uses the VERIFIED H2_HAMILTONIAN (exact = -1.137284 Ha) and EfficientSU2
ansatz from h2_vqe.py -- the real chem.py-derived Hamiltonian, not the
old toy one.

Do NOT modify ../quantum-hardware-mcp. Submits 3 real hardware jobs.

Run:
    python vqe/hardware_clean.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from h2_vqe import make_ansatz, H2_HAMILTONIAN, exact_ground_state_energy, N_ANSATZ_PARAMS

from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import StatevectorEstimator
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from dotenv import load_dotenv
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "..", ".env"))

from mcp_backend import best_device  # MCP used for exactly ONE thing: device selection

CHEM_ACCURACY_TARGET_HA = 0.0001
HARTREE_TO_KCAL_MOL = 627.5094740631
TWO_QUBIT_GATE_NAMES = ("cx", "cz", "ecr", "rzx", "swap", "iswap")


def optimize_locally():
    """
    Optimize the EfficientSU2 ansatz against the verified H2_HAMILTONIAN
    on the noiseless local simulator -- same approach as mcp_energy.py /
    h2_hardware_direct.py.
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
    return result.x, float(result.fun)


def main():
    print("\n" + "=" * 70)
    print("  H2 on REAL hardware -- our own engine, MCP used only for device pick")
    print("=" * 70)

    exact = exact_ground_state_energy()
    print(f"  exact ground state = {exact:.6f} Ha\n")

    # --- Step 1: local optimization + verification ------------------------
    print("Step 1: optimizing EfficientSU2 ansatz on the statevector simulator...")
    params, sim_energy = optimize_locally()
    sim_error = abs(sim_energy - exact)
    print(f"  simulator energy = {sim_energy:.6f} Ha  (error vs exact = {sim_error:.6f} Ha)")
    if sim_error >= CHEM_ACCURACY_TARGET_HA:
        print(f"  STOPPING: did not reach the {CHEM_ACCURACY_TARGET_HA} Ha verification target.")
        sys.exit(1)
    print(f"  -> verified within {CHEM_ACCURACY_TARGET_HA} Ha of exact.\n")

    # --- Step 2: device selection, MCP used for this ONE thing -------------
    print("Step 2: picking the best device via Lokesh's MCP server (compare_devices)...")
    device_name = best_device()
    if not device_name:
        print("ERROR: could not get best device via the MCP connector. Aborting.")
        sys.exit(1)
    print(f"  Best device: {device_name}\n")

    # --- Step 3: connect directly, transpile with auto-placement -----------
    print("Step 3: connecting DIRECTLY to IBM Quantum (no MCP routing of the measurement)...")
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        sys.exit(1)
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(device_name)

    qc = make_ansatz(params)
    # No initial_layout -- let the transpiler's own noise+connectivity-aware
    # SabreLayout choose physical qubits automatically.
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)
    isa_observable = SparsePauliOp.from_list(list(H2_HAMILTONIAN.items())).apply_layout(isa_circuit.layout)

    depth = isa_circuit.depth()
    ops = dict(isa_circuit.count_ops())
    two_qubit_gates = sum(count for name, count in ops.items() if name in TWO_QUBIT_GATE_NAMES)
    print(f"  Transpiled circuit: depth = {depth}, 2-qubit gates = {two_qubit_gates}\n")

    # --- Step 4: measure at three resilience levels -------------------------
    print("Step 4: measuring at THREE resilience levels (3 jobs total)...")
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
        energy = float(result[0].data.evs)
        hw_results[label] = {"resilience_level": level, "energy_ha": energy, "job_id": job_id}
        err = abs(energy - exact)
        print(f"    {label}: E = {energy:.6f} Ha  error = {err:.6f} Ha  ({err * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")

    # --- Final report -------------------------------------------------------
    print("\n" + "=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print(f"  Exact energy               = {exact:.6f} Ha")
    print(f"  Simulator-ansatz energy    = {sim_energy:.6f} Ha  (error {sim_error:.6f} Ha)")
    print(f"  Device                     = {device_name}")
    print(f"  Transpiled depth/2q gates  = {depth} / {two_qubit_gates}\n")
    for label in ("raw", "level1", "zne"):
        r = hw_results[label]
        err = abs(r["energy_ha"] - exact)
        print(f"  {label:8s}: E = {r['energy_ha']:.6f} Ha  err = {err:.6f} Ha  "
              f"({err * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)  job={r['job_id']}")
    print("=" * 70)

    results = {
        "exact_energy_ha": round(exact, 6),
        "simulator_energy_ha": round(sim_energy, 6),
        "simulator_error_ha": round(sim_error, 6),
        "device": device_name,
        "transpiled_depth": depth,
        "transpiled_two_qubit_gates": two_qubit_gates,
        "hardware_results": {
            label: {
                "resilience_level": r["resilience_level"],
                "job_id": r["job_id"],
                "energy_ha": round(r["energy_ha"], 6),
                "error_ha": round(abs(r["energy_ha"] - exact), 6),
                "error_kcal_mol": round(abs(r["energy_ha"] - exact) * HARTREE_TO_KCAL_MOL, 4),
            }
            for label, r in hw_results.items()
        },
    }
    out_path = os.path.join(os.path.dirname(__file__), "hardware_clean_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
