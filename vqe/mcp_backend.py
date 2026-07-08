#!/usr/bin/env python3
"""
mcp_backend.py — READ-ONLY connector to Lokesh's Quantum Hardware MCP server
=============================================================================
This file never submits a hardware job and never spends queue time or
quota. It only calls "look, don't touch" tools that already exist in
Lokesh's MCP server, which lives in the sibling project folder
../quantum-hardware-mcp/server.py (see this repo's README credits section):

  - best_qubits(device, n)                 -> which physical qubits are
                                               cleanest right now (just reads
                                               live calibration data)
  - circuit_report(device, qasm, version)  -> dry-run fidelity estimate;
                                               transpiles the circuit but
                                               never submits it
  - debug_circuit(qasm, device, version)   -> static + hardware bug checks;
                                               also never submits anything

Why this matters for hardware VQE:
  Before spending real queue time (and, on some plans, real quota) on a
  chemistry circuit, we want to:
    1. Know which physical qubits are least noisy TODAY (best_qubits_for)
    2. Confirm the circuit is well-formed and get a fidelity estimate
       BEFORE ever submitting a job (validate_circuit)

How this connects to Lokesh's server:
  We import his server.py directly as a plain Python module -- the exact
  same technique this project's own tests/test_server_tools.py already
  uses for the local copy of server.py (`from server import ...`). We
  never copy, edit, or otherwise touch anything inside
  ../quantum-hardware-mcp -- we only read it via `import`.

Run a quick demo (read-only, no quota used):
    python vqe/mcp_backend.py
"""
import os
import sys
import json

# Lokesh's Quantum Hardware MCP server lives one folder up from THIS
# project, as a sibling repo. We reach it by relative path only -- nothing
# in that folder is ever written to or modified by this file.
_MCP_SERVER_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "quantum-hardware-mcp")
)


def _load_mcp_server():
    """
    Import Lokesh's server.py module on demand (not at file-import time),
    so this file itself can always be imported safely even when the
    sibling project isn't checked out, or its IBM token isn't set up yet.

    Returns:
        The imported server module on success, or None on any failure.
        Every failure path prints a clear, plain-English reason -- callers
        just need to check for None and fall back gracefully, never crash.
    """
    if not os.path.isdir(_MCP_SERVER_DIR):
        print(
            f"[mcp_backend] MCP server not found at {_MCP_SERVER_DIR} -- "
            "is the quantum-hardware-mcp repo checked out as a sibling folder "
            "next to this project?"
        )
        return None

    if _MCP_SERVER_DIR not in sys.path:
        sys.path.insert(0, _MCP_SERVER_DIR)

    try:
        import server as mcp_server  # Lokesh's server.py -- read-only import
    except Exception as e:
        print(f"[mcp_backend] Failed to import Lokesh's server.py: {e}")
        return None

    return mcp_server


def best_device() -> str:
    """
    Ask Lokesh's MCP server which IBM Quantum device is best RIGHT NOW,
    blending quality (2-qubit gate error) and availability (queue depth) --
    while skipping any device currently in MAINTENANCE.

    Why the extra check: compare_devices() ranks by operational/queue/
    error, but "operational" and "in maintenance" turned out to be
    different things -- a backend can rank at the top of compare_devices
    while still being in a maintenance window (a separate status field,
    only exposed by list_devices()). Submitting there gets a job stuck
    QUEUED indefinitely with no real chance of running (discovered the
    hard way in hardware_covalent.py). This cross-checks both read-only
    tools and picks the highest-ranked device that is NOT in maintenance.

    READ-ONLY: calls server.compare_devices('combined') and
    server.list_devices(), both read-only. No circuit is built, no job is
    submitted, no queue time or quota is used.

    Returns:
        The name of the best non-maintenance device (str), or None if the
        MCP server/token isn't available, or every device is currently in
        maintenance -- check the printed message for why.
    """
    mcp_server = _load_mcp_server()
    if mcp_server is None:
        return None

    try:
        ranked_raw = mcp_server.compare_devices("combined")
        status_raw = mcp_server.list_devices()
    except Exception as e:
        print(f"[mcp_backend] best_device() failed: {e}")
        return None

    ranked = json.loads(ranked_raw)
    if "error" in ranked:
        print(f"[mcp_backend] best_device() failed: {ranked['error']}")
        return None

    statuses = json.loads(status_raw)
    if isinstance(statuses, dict) and "error" in statuses:
        print(f"[mcp_backend] best_device() failed: {statuses['error']}")
        return None
    status_by_name = {d["name"]: str(d.get("status", "")).lower() for d in statuses}

    for device in ranked["devices"]:
        name = device["name"]
        if status_by_name.get(name) == "maintenance":
            print(f"[mcp_backend] best_device(): skipping {name} -- currently in maintenance")
            continue
        return name

    print("[mcp_backend] best_device(): every ranked device is currently in maintenance")
    return None


def best_qubits_for(device_name: str, n: int = 4) -> list:
    """
    Ask Lokesh's MCP server which `n` physical qubits on `device_name` are
    currently cleanest (lowest readout error + lowest 2-qubit gate error),
    so a hardware VQE run can be steered onto good qubits instead of
    whatever the compiler happens to pick by default.

    READ-ONLY: calls server.best_qubits(), which only reads today's live
    calibration data from IBM. No circuit is built, no job is submitted,
    no queue time or quota is used.

    Args:
        device_name: IBM backend name, e.g. "ibm_kingston".
        n: how many of the cleanest qubits to return.

    Returns:
        A list of {"qubit": <int>, "score": <float>} dicts, best (lowest
        score) first. Empty list if the MCP server or IBM token isn't
        available right now -- check the printed message for why.
    """
    mcp_server = _load_mcp_server()
    if mcp_server is None:
        return []

    try:
        raw = mcp_server.best_qubits(device_name, n)
    except Exception as e:
        # _get_service() inside server.py raises a plain ValueError if
        # IBM_QUANTUM_TOKEN is missing/invalid in the MCP server's OWN
        # .env, and IBM's API can also raise its own errors (bad token,
        # network issue). Catch everything here so this never crashes.
        print(f"[mcp_backend] best_qubits_for('{device_name}') failed: {e}")
        return []

    data = json.loads(raw)
    if "error" in data:
        print(f"[mcp_backend] best_qubits_for('{device_name}') failed: {data['error']}")
        return []

    # Keep only what callers actually need: qubit index + its score.
    return [{"qubit": q["qubit"], "score": q["score"]} for q in data["best_qubits"]]


def validate_circuit(device_name: str, qasm_string: str, qasm_version: int = 2) -> dict:
    """
    Check whether a circuit will survive on `device_name` BEFORE spending
    any queue time on it, by combining two of Lokesh's read-only MCP tools:

      - circuit_report: transpiles the circuit against the real backend
        and estimates fidelity (the probability the result comes back
        correct, based on today's calibration data).
      - debug_circuit: static + hardware sanity checks (missing
        measurements, too many qubits for the backend, circuit too deep
        relative to T2 coherence time, etc.)

    READ-ONLY: both tools only transpile/analyse the circuit. Neither one
    ever submits a job.

    Args:
        device_name: IBM backend to check against, e.g. "ibm_kingston".
        qasm_string: the circuit, in OpenQASM format.
        qasm_version: 2 (default) or 3.

    Returns:
        {
          "estimated_fidelity": float or None,
          "verdict":            str,   # plain-English recommendation
          "issues":             list,  # from debug_circuit, [] if none found
          "safe_to_submit":     bool,  # False if MCP unavailable, circuit_report
                                        # failed, or debug_circuit found an ERROR
        }
    """
    mcp_server = _load_mcp_server()
    if mcp_server is None:
        return {
            "estimated_fidelity": None,
            "verdict": "MCP server unavailable -- cannot validate this circuit.",
            "issues": [],
            "safe_to_submit": False,
        }

    try:
        report_raw = mcp_server.circuit_report(device_name, qasm_string, qasm_version)
    except Exception as e:
        report_raw = json.dumps({"error": str(e)})

    try:
        debug_raw = mcp_server.debug_circuit(qasm_string, device_name, qasm_version)
    except Exception as e:
        debug_raw = json.dumps({"issues": [], "summary": str(e), "safe_to_submit": False})

    report = json.loads(report_raw)
    debug = json.loads(debug_raw)

    issues = debug.get("issues", [])
    # Only safe to submit if BOTH tools are happy: debug_circuit found no
    # blocking issues, AND circuit_report actually succeeded (didn't error
    # out on a bad device name or unparseable QASM).
    safe_to_submit = debug.get("safe_to_submit", False) and "error" not in report

    if "error" in report:
        verdict = f"circuit_report failed: {report['error']}"
        estimated_fidelity = None
    else:
        verdict = report.get("verdict")
        estimated_fidelity = report.get("estimated_fidelity")

    return {
        "estimated_fidelity": estimated_fidelity,
        "verdict": verdict,
        "issues": issues,
        "safe_to_submit": safe_to_submit,
    }


if __name__ == "__main__":
    # Quick read-only demo: which 4 physical qubits on ibm_kingston are
    # cleanest right now? No job submitted, no quota used.
    result = best_qubits_for("ibm_kingston", 4)
    print(json.dumps(result, indent=2))
