#!/usr/bin/env python3
"""
test_bell_state.py
------------------
End-to-end smoke test for the three job-submission MCP tools.

What it does:
  1. Finds the least-busy IBM backend via queue_status()
  2. Submits a 2-qubit Bell state circuit via submit_job()
  3. Polls job_status() until the job finishes (or times out)
  4. Prints the measurement counts from job_results()

A Bell state (|00> + |11>) / sqrt(2) should produce roughly 50% "00"
and 50% "11" outcomes. Any other result means something went wrong.

Usage:
    .venv/bin/python test_bell_state.py

The script needs IBM_QUANTUM_TOKEN in your .env file.
Real hardware jobs typically take 1-15 minutes depending on queue depth.
"""

import json
import time
import sys

# Import the tool functions directly from server.py
from server import queue_status, submit_job, job_status, job_results

# --------------------------------------------------------------------------
# Bell state circuit — OpenQASM 2.0
#
# Step by step:
#   h q[0]          — Hadamard puts qubit 0 into superposition: |0> + |1>
#   cx q[0],q[1]    — CNOT entangles qubit 1 with qubit 0
#   measure ...     — collapse the superposition to a classical bit
#
# Result: measuring always gives 00 or 11, never 01 or 10.
# That's entanglement — the qubits are perfectly correlated.
# --------------------------------------------------------------------------
BELL_CIRCUIT = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
"""

SHOTS       = 1024
POLL_EVERY  = 20   # seconds between status checks
TIMEOUT_MIN = 30   # give up after this many minutes


def main():
    print("=" * 60)
    print("Bell State — End-to-End IBM Hardware Test")
    print("=" * 60)

    # ── Step 1: find least-busy backend ──────────────────────────────
    print("\n[1/4] Finding least-busy backend...")
    queues = json.loads(queue_status())

    # queue_status returns sorted by pending_jobs ascending — first entry is best
    best = next((d for d in queues if d["operational"]), None)
    if not best:
        print("ERROR: No operational backends found.")
        sys.exit(1)

    device = best["name"]
    print(f"      → {device}  ({best['pending_jobs']} jobs in queue)")

    # ── Step 2: submit the circuit ────────────────────────────────────
    print(f"\n[2/4] Submitting Bell state circuit ({SHOTS} shots)...")
    submission = json.loads(submit_job(device, BELL_CIRCUIT, shots=SHOTS))

    if "error" in submission:
        print(f"ERROR: {submission['error']}")
        sys.exit(1)

    job_id = submission["job_id"]
    print(f"      → job_id: {job_id}")
    print(f"      → initial status: {submission['status']}")

    # ── Step 3: poll until done ───────────────────────────────────────
    print(f"\n[3/4] Waiting for job to complete (polling every {POLL_EVERY}s)...")
    deadline = time.time() + TIMEOUT_MIN * 60
    last_status = None

    while time.time() < deadline:
        info = json.loads(job_status(job_id))
        status = info["status"]

        if status != last_status:
            line = f"      [{time.strftime('%H:%M:%S')}] {status}"
            if "queue_position" in info:
                line += f"  (position {info['queue_position']} in queue)"
            print(line)
            last_status = status

        if status == "DONE":
            break
        if status in ("ERROR", "CANCELLED"):
            print(f"\nERROR: Job ended with status {status}")
            if "error_message" in info:
                print(f"       {info['error_message']}")
            sys.exit(1)

        time.sleep(POLL_EVERY)
    else:
        print(f"\nTIMEOUT: Job did not complete within {TIMEOUT_MIN} minutes.")
        print(f"         Check manually with: job_status('{job_id}')")
        sys.exit(1)

    # ── Step 4: retrieve and interpret results ────────────────────────
    print(f"\n[4/4] Retrieving results...")
    result = json.loads(job_results(job_id))

    if "error" in result:
        print(f"ERROR: {result['error']}")
        sys.exit(1)

    counts      = result["counts"]
    total_shots = result["total_shots"]

    print(f"\n{'=' * 60}")
    print(f"Results from {device}")
    print(f"{'=' * 60}")
    print(f"Total shots : {total_shots}")
    print(f"Counts      : {json.dumps(counts, indent=14)}")

    # Sanity check: a real Bell state gives only "00" and "11"
    expected_keys = {"00", "11"}
    got_keys      = set(counts.keys())
    unexpected    = got_keys - expected_keys

    correlated_shots = counts.get("00", 0) + counts.get("11", 0)
    correlated_pct   = 100 * correlated_shots / total_shots

    print(f"\nCorrelated outcomes (00 + 11): {correlated_pct:.1f}% of shots")

    if correlated_pct >= 90:
        print("PASS — strong Bell correlations observed (>=90% in 00 + 11).")
    elif correlated_pct >= 70:
        print("MARGINAL — Bell correlations present but noisy (<90%). "
              "Hardware error rates may be high.")
    else:
        print("FAIL — expected 00/11 dominance but got unexpected outcomes. "
              "Check device error rates.")

    if unexpected:
        print(f"Note: unexpected outcomes seen: {unexpected}  "
              "(expected only 00 and 11 for a Bell state)")

    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
