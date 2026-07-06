#!/usr/bin/env python3
"""
hardware_fragmentation_mcp.py — same covalent-bond fragmentation test as
hardware_fragmentation.py (molecular tailoring on a bonded H6 chain), but
each fragment's hardware measurement is routed THROUGH Lokesh's Quantum
Hardware MCP server (estimate_expectation + job_status), instead of
connecting to qiskit-ibm-runtime directly.

Reuses, unmodified:
  - process_fragment() from hardware_fragmentation.py -- identical fragment
    Hamiltonian construction and COBYLA-optimized EfficientSU2 local
    verification. Only the HARDWARE step differs between the two files.
  - measure_energy_via_mcp() from mcp_energy.py -- identical MCP submit/
    poll/retrieve logic already proven correct for H2.

Device selection is checked live via the MCP server every run (never
hardcoded), so this always targets whichever backend is least busy right
now -- including any new devices Lokesh's server exposes.

Run (submits 3 real hardware jobs, via MCP):
    python vqe/hardware_fragmentation_mcp.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from hardware_fragmentation import (
    process_fragment, HARTREE_TO_KCAL_MOL, FULL_EXACT_H6_HA, SIM_TAILORING_H6_HA,
)
from mcp_energy import measure_energy_via_mcp  # reuse the proven MCP submit/poll/retrieve logic

from qiskit.circuit.library import EfficientSU2
from qiskit import qasm2

# Lokesh's Quantum Hardware MCP server -- same sibling-path technique as
# mcp_backend.py / mcp_energy.py. Read-only import, never modifies it.
_MCP_SERVER_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "quantum-hardware-mcp")
)
if _MCP_SERVER_DIR not in sys.path:
    sys.path.insert(0, _MCP_SERVER_DIR)
import server


def _pick_least_busy_device_via_mcp():
    """
    Check every device Lokesh's server can currently reach (via
    compare_devices, read-only) and return the least-busy one. Checked
    fresh every run -- if he's added more devices, this picks them up
    automatically instead of relying on a cached/hardcoded list.
    """
    data = json.loads(server.compare_devices("queue"))
    for d in data["devices"]:
        print(f"    {d['name']}: {d['pending_jobs']} pending jobs")
    best = data["devices"][0]
    print(f"  Least-busy device right now (via MCP): {best['name']} ({best['pending_jobs']} pending jobs)\n")
    return best["name"]


def measure_fragment_via_mcp(fragment, device_name):
    """
    Submit ONE real hardware job for this fragment's already-optimized
    ansatz, THROUGH Lokesh's MCP server, instead of connecting to
    qiskit-ibm-runtime directly (that's what hardware_fragmentation.py does).
    """
    ansatz = EfficientSU2(fragment["n_qubits"], reps=fragment["reps"])
    qc = ansatz.assign_parameters(fragment["params"])
    qasm_string = qasm2.dumps(qc)

    # qop_bare is a SparsePauliOp -- convert to the {pauli_label: coeff}
    # dict shape measure_energy_via_mcp expects (same shape as H2_HAMILTONIAN).
    qop_bare = fragment["qop_bare"]
    hamiltonian_terms = {
        label: complex(coeff).real
        for label, coeff in zip(qop_bare.paulis.to_labels(), qop_bare.coeffs)
    }

    electronic_energy, job_id = measure_energy_via_mcp(device_name, qasm_string, hamiltonian_terms)
    hardware_energy = electronic_energy + fragment["enuc"]
    return hardware_energy, job_id


def main():
    d = 1.0  # Angstrom, same bonded H6 chain as hardware_fragmentation.py

    print("\n" + "=" * 70)
    print("  Covalent-bond fragmentation on REAL hardware, VIA MCP SERVER -- H6 chain")
    print("=" * 70)

    # --- Step 1: build + locally verify each fragment -------------------
    # Identical pipeline to hardware_fragmentation.py -- same Hamiltonians,
    # same COBYLA-optimized ansatz. Runs independently in this process
    # (does not reuse any other running script's in-progress work).
    fragments = {
        "block1":  process_fragment("block1 (atoms 0-3)",  [0, 1, 2, 3], d),
        "block2":  process_fragment("block2 (atoms 2-5)",  [2, 3, 4, 5], d),
        "overlap": process_fragment("overlap (atoms 2-3)", [2, 3],       d),
    }

    # --- Step 2: pick ONE device, checked fresh via the MCP server -------
    print("=" * 70)
    print("  Checking devices via MCP server")
    print("=" * 70)
    device_name = _pick_least_busy_device_via_mcp()

    # --- Step 3: measure each fragment ONCE on real hardware, via MCP ----
    print("=" * 70)
    print("  Hardware measurements via MCP server (3 jobs total)")
    print("=" * 70)
    for frag in fragments.values():
        hw_energy, job_id = measure_fragment_via_mcp(frag, device_name)
        frag["hardware_energy_ha"] = hw_energy
        frag["job_id"] = job_id
        err = abs(hw_energy - frag["exact_energy_ha"])
        print(f"  {frag['name']}: hardware = {hw_energy:.6f} Ha, "
              f"exact = {frag['exact_energy_ha']:.6f} Ha, error = {err:.6f} Ha")

    # --- Step 4: reassemble via molecular tailoring ----------------------
    e_hw_h6 = (fragments["block1"]["hardware_energy_ha"]
               + fragments["block2"]["hardware_energy_ha"]
               - fragments["overlap"]["hardware_energy_ha"])
    err_vs_exact = abs(e_hw_h6 - FULL_EXACT_H6_HA)
    err_vs_sim_tailor = abs(e_hw_h6 - SIM_TAILORING_H6_HA)

    print("\n" + "=" * 70)
    print("  FRAGMENT SUMMARY (via MCP)")
    print("=" * 70)
    print(f"  {'fragment':22s} | {'qubits':>6s} | {'exact (Ha)':>12s} | {'hardware (Ha)':>14s} | {'error (Ha)':>10s}")
    for frag in fragments.values():
        err = abs(frag["hardware_energy_ha"] - frag["exact_energy_ha"])
        print(f"  {frag['name']:22s} | {frag['n_qubits']:6d} | {frag['exact_energy_ha']:12.6f} | "
              f"{frag['hardware_energy_ha']:14.6f} | {err:10.6f}")

    print("\n" + "=" * 70)
    print("  REASSEMBLED H6 (molecular tailoring, via MCP server)")
    print("=" * 70)
    print("  E_hardware(H6) = E_hw(block1) + E_hw(block2) - E_hw(overlap)")
    print(f"                 = {fragments['block1']['hardware_energy_ha']:.6f} "
          f"+ {fragments['block2']['hardware_energy_ha']:.6f} "
          f"- {fragments['overlap']['hardware_energy_ha']:.6f}")
    print(f"                 = {e_hw_h6:.6f} Ha\n")
    print(f"  Full exact H6                    = {FULL_EXACT_H6_HA:.6f} Ha")
    print(f"  Simulator tailoring H6           = {SIM_TAILORING_H6_HA:.6f} Ha")
    print(f"  Hardware tailoring H6 (via MCP)  = {e_hw_h6:.6f} Ha\n")
    print(f"  Error vs full exact     = {err_vs_exact:.6f} Ha  ({err_vs_exact * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print(f"  Error vs sim tailoring  = {err_vs_sim_tailor:.6f} Ha  ({err_vs_sim_tailor * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print("=" * 70)

    # --- Save results -------------------------------------------------------
    results = {
        "device": device_name,
        "route": "mcp_server",
        "fragments": {
            key: {
                "name": frag["name"],
                "indices": frag["indices"],
                "n_qubits": frag["n_qubits"],
                "reps": frag["reps"],
                "exact_energy_ha": round(frag["exact_energy_ha"], 6),
                "sim_vqe_energy_ha": round(frag["sim_energy_ha"], 6),
                "sim_vqe_error_ha": round(frag["sim_error_ha"], 6),
                "hardware_energy_ha": round(frag["hardware_energy_ha"], 6),
                "hardware_error_ha": round(abs(frag["hardware_energy_ha"] - frag["exact_energy_ha"]), 6),
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
    out_path = os.path.join(os.path.dirname(__file__), "hardware_fragmentation_mcp_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
