#!/usr/bin/env python3
"""
hardware_covalent.py — covalent-bond fragmentation (molecular tailoring)
on REAL IBM hardware, for a bonded H6 chain (d=1.0 A):

    E(H6) = E_hw(atoms 0-3) + E_hw(atoms 2-5) - E_hw(atoms 2-3)

This is the honest frontier experiment for this repo: three real
8-/8-/4-qubit measurements, reassembled via inclusion-exclusion. The
result WILL be noisy -- that is expected and reported truthfully, not
tuned or adjusted.

Fixes applied from prior failed attempts (see hardware_fragmentation.py
and hardware_mitigation_test.py for the history):
  - NO forced qubit layout. best_qubits_for() previously picked
    "individually clean" qubits that turned out to be almost entirely
    disconnected on the chip (1 edge out of 8 chosen qubits), forcing a
    flood of SWAP gates and >1 Ha errors. generate_preset_pass_manager's
    own layout stage is connectivity-aware -- this script lets it choose.
  - ExcitationPreserving(reps=2) is used AS-IS: whatever error it reaches
    on the simulator (~0.0047 Ha expected, per fragment_ansatz_test.py)
    is accepted and reported, not gated behind a stricter threshold that
    earlier attempts (with the wrong ansatz) could never reach.
  - Bounded job wait: previous runs discovered the IBM account can hit a
    usage limit where jobs sit QUEUED forever. This script polls with a
    timeout and STOPS with a clear message instead of blocking forever.

MCP usage is restricted to exactly ONE thing: best_device() (device
selection via server.compare_devices). The measurement itself never
routes through the MCP server. ../quantum-hardware-mcp is never modified.

Run (submits 3 real hardware jobs, ZNE only):
    python vqe/hardware_covalent.py
"""
import os
import sys
import json
import time
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

from mcp_backend import best_device  # MCP used for exactly ONE thing: device selection

HARTREE_TO_KCAL_MOL = 627.5094740631
TWO_QUBIT_GATE_NAMES = ("cx", "cz", "ecr", "rzx", "swap", "iswap")
JOB_MAX_WAIT_SECONDS = 480
JOB_POLL_INTERVAL_SECONDS = 15

FULL_EXACT_H6_HA = -3.236066
SIM_TAILORING_H6_HA = -3.231625


# --------------------------------------------------------------------------
# Local (simulator) optimization -- ExcitationPreserving, accepted as-is
# --------------------------------------------------------------------------

def optimize_fragment_locally(qop_bare, enuc, exact, n_qubits, nelec, reps=2, n_restarts=5, maxiter=500):
    """
    Hartree-Fock state prep + ExcitationPreserving(reps=2), optimized with
    L-BFGS-B over several restarts. No hard convergence gate -- whatever
    error this reaches (expected ~0.0047 Ha, per fragment_ansatz_test.py)
    is accepted and carried forward honestly.
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
        if error < 0.0005:  # essentially exact already -- no need for more restarts
            break

    qc = make_circuit(best_params)
    return best_params, best_energy, abs(best_energy - exact), qc


def process_fragment(name, indices, d):
    """Build one fragment's Hamiltonian and optimize its ansatz locally (no hardware yet)."""
    geom = hchain(indices, d)
    nelec = len(indices)
    print(f"\n--- Fragment '{name}' (atoms {indices}, {nelec} electrons) ---")

    qop_bare, _qop_penalized, enuc = _build_fragment_qop(geom, nelec)
    n_qubits = qop_bare.num_qubits
    exact = exact_energy_no_penalty(qop_bare, enuc)
    print(f"  qubits = {n_qubits}, exact ground state = {exact:.6f} Ha")

    print("  Optimizing ExcitationPreserving(reps=2) locally (L-BFGS-B)...")
    params, sim_energy, sim_error, qc = optimize_fragment_locally(qop_bare, enuc, exact, n_qubits, nelec)
    print(f"  -> simulator best: E = {sim_energy:.6f} Ha, error = {sim_error:.6f} Ha")

    return {
        "name": name, "indices": indices, "nelec": nelec, "n_qubits": n_qubits,
        "exact_energy_ha": exact, "sim_energy_ha": sim_energy, "sim_error_ha": sim_error,
        "params": params, "qop_bare": qop_bare, "enuc": enuc, "circuit": qc,
    }


# --------------------------------------------------------------------------
# Real hardware measurement, with a BOUNDED wait (no indefinite blocking)
# --------------------------------------------------------------------------

def submit_and_wait(estimator, pubs, max_wait_seconds=JOB_MAX_WAIT_SECONDS,
                     poll_interval=JOB_POLL_INTERVAL_SECONDS):
    """
    Submit a job and poll its status with a bounded timeout, instead of
    calling job.result() directly (which blocks forever). If the account
    has hit a usage limit, jobs sit QUEUED indefinitely -- this raises a
    clear TimeoutError instead of hanging, per explicit instructions to
    stop rather than wait indefinitely.
    """
    job = estimator.run(pubs)
    job_id = job.job_id()
    print(f"    job_id={job_id}")

    elapsed = 0
    status = str(job.status())
    while elapsed < max_wait_seconds:
        status = str(job.status())
        print(f"    status={status} (elapsed {elapsed}s)")
        if status == "DONE":
            return job.result(), job_id
        if status in ("ERROR", "CANCELLED"):
            raise RuntimeError(f"Job {job_id} ended with status {status} -- cannot continue.")
        time.sleep(poll_interval)
        elapsed += poll_interval

    raise TimeoutError(
        f"Job {job_id} is still '{status}' after {max_wait_seconds}s. This usually means "
        f"the IBM account has hit its usage limit (workloads queued but not executing). "
        f"STOPPING per instructions -- not waiting indefinitely."
    )


def measure_fragment_on_hardware(fragment, backend):
    """
    Transpile with AUTO-PLACEMENT (no forced qubits -- see module docstring
    for why) and measure the bare electronic Hamiltonian on real hardware,
    ONE job, resilience_level=2 (ZNE), 4096 shots.
    """
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(fragment["circuit"])
    isa_observable = fragment["qop_bare"].apply_layout(isa_circuit.layout)

    depth = isa_circuit.depth()
    ops = dict(isa_circuit.count_ops())
    two_qubit_gates = sum(count for name, count in ops.items() if name in TWO_QUBIT_GATE_NAMES)
    print(f"  Transpiled (auto-placed): depth = {depth}, 2-qubit gates = {two_qubit_gates}")

    estimator = EstimatorV2(mode=backend)
    estimator.options.resilience_level = 2  # ZNE
    estimator.options.default_shots = 4096

    result, job_id = submit_and_wait(estimator, [(isa_circuit, isa_observable)])
    electronic_energy = float(result[0].data.evs)
    hw_energy = electronic_energy + fragment["enuc"]
    return hw_energy, job_id, depth, two_qubit_gates


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    d = 1.0  # Angstrom, same bonded H6 chain as covalent_fragment.py / hardware_fragmentation.py

    print("\n" + "=" * 70)
    print("  Covalent-bond fragmentation on REAL hardware -- H6 chain")
    print("  (honest frontier experiment -- noisy result expected)")
    print("=" * 70)

    # --- Step 1: build + locally optimize each fragment (no hardware yet) --
    fragments = {
        "block1":  process_fragment("block1 (atoms 0-3)",  [0, 1, 2, 3], d),
        "block2":  process_fragment("block2 (atoms 2-5)",  [2, 3, 4, 5], d),
        "overlap": process_fragment("overlap (atoms 2-3)", [2, 3],       d),
    }

    # --- Step 2: device selection ONLY via MCP -----------------------------
    print("\n" + "=" * 70)
    print("  Picking the best device via Lokesh's MCP server (device selection only)")
    print("=" * 70)
    device_name = best_device()
    if not device_name:
        print("ERROR: could not get best device via the MCP connector. Aborting.")
        sys.exit(1)
    print(f"  Best device: {device_name}")

    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        sys.exit(1)
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(device_name)

    # --- Step 3: measure each fragment ONCE on real hardware (ZNE) --------
    print("\n" + "=" * 70)
    print("  Hardware measurements (3 jobs total, resilience_level=2 / ZNE)")
    print("=" * 70)
    try:
        for frag in fragments.values():
            print(f"\n  Fragment '{frag['name']}':")
            hw_energy, job_id, depth, two_q = measure_fragment_on_hardware(frag, backend)
            frag["hardware_energy_ha"] = hw_energy
            frag["job_id"] = job_id
            frag["transpiled_depth"] = depth
            frag["transpiled_two_qubit_gates"] = two_q
            err = abs(hw_energy - frag["exact_energy_ha"])
            print(f"    hardware energy = {hw_energy:.6f} Ha  (error vs exact = {err:.6f} Ha)")
    except (TimeoutError, RuntimeError) as e:
        print(f"\n{'=' * 70}\n  STOPPING\n{'=' * 70}")
        print(f"  {e}")
        print("  No reassembly performed, no results file written. Not waiting indefinitely.")
        sys.exit(1)

    # --- Step 4: reassemble via molecular tailoring -------------------------
    e_hw_h6 = (fragments["block1"]["hardware_energy_ha"]
               + fragments["block2"]["hardware_energy_ha"]
               - fragments["overlap"]["hardware_energy_ha"])
    err_vs_exact = abs(e_hw_h6 - FULL_EXACT_H6_HA)
    err_vs_sim_tailor = abs(e_hw_h6 - SIM_TAILORING_H6_HA)

    # --- Full table ---------------------------------------------------------
    print("\n" + "=" * 100)
    print("  FRAGMENT SUMMARY")
    print("=" * 100)
    print(f"  {'fragment':22s} | {'qubits':>6s} | {'exact (Ha)':>11s} | {'sim (Ha)':>11s} | "
          f"{'sim err':>8s} | {'hardware (Ha)':>13s} | {'hw err':>8s} | {'depth':>5s} | {'2q gates':>8s}")
    for frag in fragments.values():
        hw_err = abs(frag["hardware_energy_ha"] - frag["exact_energy_ha"])
        print(f"  {frag['name']:22s} | {frag['n_qubits']:6d} | {frag['exact_energy_ha']:11.6f} | "
              f"{frag['sim_energy_ha']:11.6f} | {frag['sim_error_ha']:8.6f} | "
              f"{frag['hardware_energy_ha']:13.6f} | {hw_err:8.6f} | "
              f"{frag['transpiled_depth']:5d} | {frag['transpiled_two_qubit_gates']:8d}")

    print("\n" + "=" * 100)
    print("  REASSEMBLED H6 (molecular tailoring on REAL hardware)")
    print("=" * 100)
    print("  E_hardware(H6) = E_hw(block1) + E_hw(block2) - E_hw(overlap)")
    print(f"                 = {fragments['block1']['hardware_energy_ha']:.6f} "
          f"+ {fragments['block2']['hardware_energy_ha']:.6f} "
          f"- {fragments['overlap']['hardware_energy_ha']:.6f}")
    print(f"                 = {e_hw_h6:.6f} Ha\n")
    print(f"  Full exact H6           = {FULL_EXACT_H6_HA:.6f} Ha")
    print(f"  Simulator tailoring H6  = {SIM_TAILORING_H6_HA:.6f} Ha")
    print(f"  Hardware tailoring H6   = {e_hw_h6:.6f} Ha\n")
    print(f"  Error vs full exact     = {err_vs_exact:.6f} Ha  ({err_vs_exact * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print(f"  Error vs sim tailoring  = {err_vs_sim_tailor:.6f} Ha  ({err_vs_sim_tailor * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print("=" * 100)

    # --- Save results -------------------------------------------------------
    results = {
        "backend": device_name,
        "resilience_level": 2,
        "shots": 4096,
        "fragments": {
            key: {
                "name": frag["name"],
                "indices": frag["indices"],
                "n_qubits": frag["n_qubits"],
                "exact_energy_ha": round(frag["exact_energy_ha"], 6),
                "sim_energy_ha": round(frag["sim_energy_ha"], 6),
                "sim_error_ha": round(frag["sim_error_ha"], 6),
                "hardware_energy_ha": round(frag["hardware_energy_ha"], 6),
                "hardware_error_vs_exact_ha": round(abs(frag["hardware_energy_ha"] - frag["exact_energy_ha"]), 6),
                "transpiled_depth": frag["transpiled_depth"],
                "transpiled_two_qubit_gates": frag["transpiled_two_qubit_gates"],
                "job_id": frag["job_id"],
            }
            for key, frag in fragments.items()
        },
        "reassembled_hardware_h6_ha": round(e_hw_h6, 6),
        "full_exact_h6_ha": FULL_EXACT_H6_HA,
        "sim_tailoring_h6_ha": SIM_TAILORING_H6_HA,
        "error_vs_full_exact_ha": round(err_vs_exact, 6),
        "error_vs_full_exact_kcal_mol": round(err_vs_exact * HARTREE_TO_KCAL_MOL, 4),
        "error_vs_sim_tailoring_ha": round(err_vs_sim_tailor, 6),
        "error_vs_sim_tailoring_kcal_mol": round(err_vs_sim_tailor * HARTREE_TO_KCAL_MOL, 4),
    }
    out_path = os.path.join(os.path.dirname(__file__), "hardware_covalent_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
