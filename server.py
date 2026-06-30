"""
Quantum Hardware MCP Server
============================
Exposes live IBM Quantum device data to AI assistants via the MCP protocol.

Tools:
  - list_devices          : all machines + status
  - get_device_details    : deep info on one machine
  - compare_devices       : rank machines by error rate / queue / combined score
  - queue_status          : current queue depth for every machine
  - device_history        : historical snapshots for one machine over N days
  - best_qubits           : best n qubits on a machine right now (calibration-based)
  - device_on_date        : historical stats for a machine on a specific past date
  - submit_job            : compile + submit an OpenQASM 2 or 3 circuit to IBM hardware
  - job_status            : check the status of a submitted job
  - job_results           : retrieve measurement counts from a completed job
  - cancel_job            : cancel a queued or running job
  - list_jobs             : list recent jobs with status and backend
  - run_grover            : built-in Grover's search demo on real hardware
  - estimate_expectation  : run Estimator primitive to compute observable expectation values
  - circuit_report        : dry-run analysis — fidelity estimate, gate counts, qubit map
"""

import os
import json
import math
import sqlite3
import argparse
import anyio
from datetime import datetime, timezone
from typing import Optional

from qiskit import QuantumCircuit
from qiskit import qasm3 as qiskit_qasm3
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime import EstimatorV2 as Estimator

from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP
from qiskit_ibm_runtime import QiskitRuntimeService
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import JSONResponse
from starlette.requests import Request

# --------------------------------------------------------------------------
# Load .env from the same folder as this file, regardless of working directory.
# This matters because Claude Desktop may launch the server from a different
# working directory than the project root.
# --------------------------------------------------------------------------
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), ".env"))

# Create the MCP server instance.
# The name "quantum-hardware" is what Claude Desktop shows in its UI.
mcp = FastMCP("quantum-hardware")

# --------------------------------------------------------------------------
# SQLite history database
# --------------------------------------------------------------------------

# Store the database next to this file so it travels with the project.
DB_PATH = os.path.join(os.path.dirname(__file__), "devices.db")


def _init_db() -> None:
    """Create the snapshots table if it doesn't exist yet."""
    with sqlite3.connect(DB_PATH) as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS device_snapshots (
                id                INTEGER PRIMARY KEY AUTOINCREMENT,
                ts                TEXT    NOT NULL,  -- ISO 8601 UTC timestamp
                name              TEXT    NOT NULL,  -- e.g. "ibm_fez"
                num_qubits        INTEGER,
                operational       INTEGER,           -- 1 = True, 0 = False
                pending_jobs      INTEGER,
                avg_cx_error      REAL,              -- NULL when not measured
                avg_readout_error REAL               -- NULL when not measured
            )
        """)
        # Index on (name, ts) makes device_history queries fast.
        con.execute("""
            CREATE INDEX IF NOT EXISTS idx_name_ts
            ON device_snapshots (name, ts)
        """)


def _save_snapshots(rows: list[dict]) -> None:
    """
    Write one row per device into device_snapshots.

    Each dict in `rows` must have at least 'name'; all other fields are
    optional and default to None if absent.
    """
    ts = datetime.now(timezone.utc).isoformat()
    with sqlite3.connect(DB_PATH) as con:
        con.executemany(
            """
            INSERT INTO device_snapshots
                (ts, name, num_qubits, operational, pending_jobs,
                 avg_cx_error, avg_readout_error)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    ts,
                    r["name"],
                    r.get("num_qubits"),
                    int(r["operational"]) if r.get("operational") is not None else None,
                    r.get("pending_jobs"),
                    r.get("avg_cx_error"),
                    r.get("avg_readout_error"),
                )
                for r in rows
            ],
        )


# Initialise the DB as soon as the server module loads.
_init_db()


# --------------------------------------------------------------------------
# Helper: connect to IBM Quantum
# --------------------------------------------------------------------------

def _get_service() -> QiskitRuntimeService:
    """
    Build a QiskitRuntimeService from the token stored in .env.
    Raises a clear error if the token is missing so the user knows exactly
    what to fix.
    """
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        raise ValueError(
            "IBM_QUANTUM_TOKEN is not set. "
            "Create a .env file in the project folder with:\n"
            "  IBM_QUANTUM_TOKEN=your_token_here\n"
            "Get your token at https://quantum.ibm.com/account"
        )
    # channel="ibm_quantum_platform" → renamed in qiskit-ibm-runtime ≥ 0.40
    return QiskitRuntimeService(channel="ibm_quantum_platform", token=token)


def _cx_errors_for_backend(props) -> list[float]:
    """
    Pull 2-qubit gate error values from calibration properties.

    Older IBM devices use CX (CNOT) as their native 2-qubit gate.
    Newer devices (e.g. ibm_fez, ibm_marrakesh) use ECR (echoed
    cross-resonance) instead — CX is synthesised from ECR and won't
    appear in raw calibration data. We check for both, plus CZ, so
    this function works across the whole IBM fleet.

    Returns an empty list if the backend has no calibration data.
    """
    if props is None:
        return []
    TWO_QUBIT_GATES = {"cx", "ecr", "cz"}
    errors = []
    for gate in props.gates:
        if gate.gate in TWO_QUBIT_GATES and gate.parameters:
            errors.append(gate.parameters[0].value)
    return errors


# --------------------------------------------------------------------------
# Tool 1: list_devices
# --------------------------------------------------------------------------

@mcp.tool()
def list_devices() -> str:
    """
    List every IBM quantum computer this account can access.

    Returns a JSON array sorted by qubit count (largest first).
    Each entry includes: name, qubit count, operational status, queue depth.
    """
    service = _get_service()

    # service.backends() returns a list of IBMBackend objects.
    # By default it returns ALL backends you have access to.
    backends = service.backends()

    devices = []
    for backend in backends:
        # status() is a lightweight call — no calibration data, just up/down + queue.
        status = backend.status()

        devices.append({
            "name": backend.name,
            "num_qubits": backend.num_qubits,
            "status": status.status_msg,       # e.g. "active", "maintenance"
            "operational": status.operational,  # True / False
            "pending_jobs": status.pending_jobs, # current queue length
        })

    # Sort biggest machines first — handy for a quick overview
    devices.sort(key=lambda d: d["num_qubits"], reverse=True)

    _save_snapshots(devices)
    return json.dumps(devices, indent=2)


# --------------------------------------------------------------------------
# Tool 2: get_device_details
# --------------------------------------------------------------------------

@mcp.tool()
def get_device_details(device_name: str) -> str:
    """
    Deep-dive into one IBM quantum computer.

    Args:
        device_name: Machine name, e.g. "ibm_brisbane" or "ibm_sherbrooke".
                     Use list_devices first if you don't know the exact name.

    Returns JSON with:
      - Qubit count, status, queue depth
      - Average / best / worst CX (2-qubit) gate error rates
      - Average readout (measurement) error
      - Average T1 and T2 coherence times in microseconds
      - Timestamp of the last calibration run
    """
    service = _get_service()

    # Fetch this specific backend by name
    backend = service.backend(device_name)
    status = backend.status()

    # Start building the result with data that is always available
    result = {
        "name": backend.name,
        "num_qubits": backend.num_qubits,
        "status": status.status_msg,
        "operational": status.operational,
        "pending_jobs": status.pending_jobs,
    }

    # properties() returns calibration data from the most recent daily calibration.
    # Simulators and some devices return None here.
    props = backend.properties()

    if props:
        # ---- Readout error ----
        # Readout error = probability of measuring the wrong bit (0 vs 1).
        # Average across all qubits gives a device-wide quality signal.
        readout_errors = [
            props.readout_error(q)
            for q in range(backend.num_qubits)
            if props.readout_error(q) is not None
        ]
        if readout_errors:
            result["avg_readout_error"] = round(
                sum(readout_errors) / len(readout_errors), 5
            )

        # ---- CX gate error ----
        # CX error = probability the 2-qubit gate produces the wrong output.
        # Lower is better. Typical good values are < 0.01 (1%).
        cx_errors = _cx_errors_for_backend(props)
        if cx_errors:
            result["avg_cx_error"] = round(sum(cx_errors) / len(cx_errors), 5)
            result["best_cx_error"] = round(min(cx_errors), 5)   # best qubit pair
            result["worst_cx_error"] = round(max(cx_errors), 5)  # worst qubit pair

        # ---- T1 and T2 coherence times ----
        # T1 (relaxation time):  how long a qubit in |1⟩ stays in |1⟩ before
        #                         spontaneously falling to |0⟩.
        # T2 (dephasing time):   how long a qubit stays in a superposition before
        #                         the phase randomises and quantum info is lost.
        # Both are in seconds from the API; we convert to microseconds (µs)
        # because that's the conventional unit in quantum computing papers.
        t1_times = [
            props.t1(q) for q in range(backend.num_qubits) if props.t1(q) is not None
        ]
        t2_times = [
            props.t2(q) for q in range(backend.num_qubits) if props.t2(q) is not None
        ]
        if t1_times:
            result["avg_t1_us"] = round(
                sum(t1_times) / len(t1_times) * 1e6, 1  # s → µs
            )
        if t2_times:
            result["avg_t2_us"] = round(
                sum(t2_times) / len(t2_times) * 1e6, 1
            )

        # When was the last calibration run?
        result["last_calibration"] = str(props.last_update_date)

    else:
        result["note"] = "No calibration data available (simulator or uncalibrated device)"

    _save_snapshots([result])
    return json.dumps(result, indent=2)


# --------------------------------------------------------------------------
# Tool 6: best_qubits
# --------------------------------------------------------------------------

@mcp.tool()
def best_qubits(device_name: str, n: int = 5) -> str:
    """
    Return the best n individual qubits on a device based on live calibration.

    Useful for researchers who want to hand-pick qubits for a circuit rather
    than letting the compiler choose automatically.

    Args:
        device_name: Machine name, e.g. "ibm_fez".
        n:           How many qubits to return (default 5).

    Scoring formula (lower = better):
        score = readout_error + best_cx_error_for_this_qubit

    Both metrics are in the same [0, 1] range so they contribute equally.
    T1 / T2 coherence times are included as supplementary context.
    Missing metrics are penalised with 1.0 (worst possible) so qubits with
    incomplete calibration data sort to the bottom.
    """
    service = _get_service()
    backend = service.backend(device_name)
    props   = backend.properties()

    if not props:
        return json.dumps({
            "error": f"{device_name} has no calibration data available."
        })

    n = min(n, backend.num_qubits)  # can't ask for more qubits than exist

    # Build dict: qubit index → lowest 2-qubit gate error of any pair
    # involving this qubit.  Covers cx / ecr / cz (see _cx_errors_for_backend).
    TWO_QUBIT_GATES = {"cx", "ecr", "cz"}
    qubit_best_cx: dict[int, float] = {}
    for gate in props.gates:
        if gate.gate in TWO_QUBIT_GATES and gate.parameters:
            err = gate.parameters[0].value
            for q in gate.qubits:
                if q not in qubit_best_cx or err < qubit_best_cx[q]:
                    qubit_best_cx[q] = err

    # Score and collect every qubit
    qubit_data = []
    for q in range(backend.num_qubits):
        ro  = props.readout_error(q)
        cx  = qubit_best_cx.get(q)
        # T1/T2 are missing for some qubits on some devices — catch gracefully
        try:
            t1 = props.t1(q)
        except Exception:
            t1 = None
        try:
            t2 = props.t2(q)
        except Exception:
            t2 = None

        score = (ro if ro is not None else 1.0) + (cx if cx is not None else 1.0)

        qubit_data.append({
            "qubit":          q,
            "score":          round(score, 6),
            "readout_error":  round(ro, 5)       if ro  is not None else None,
            "best_cx_error":  round(cx, 5)       if cx  is not None else None,
            "t1_us":          round(t1 * 1e6, 1) if t1  is not None else None,
            "t2_us":          round(t2 * 1e6, 1) if t2  is not None else None,
        })

    qubit_data.sort(key=lambda q: q["score"])

    return json.dumps(
        {
            "device":   device_name,
            "n":        n,
            "scoring":  "readout_error + best_cx_error (lower = better). "
                        "T1/T2 shown for context but not in score.",
            "best_qubits": qubit_data[:n],
        },
        indent=2,
    )


# --------------------------------------------------------------------------
# Tool 3: compare_devices
# --------------------------------------------------------------------------

@mcp.tool()
def compare_devices(sort_by: str = "cx_error") -> str:
    """
    Rank all accessible IBM quantum computers by a quality metric.

    Args:
        sort_by: Ranking criterion. Choose one of:
                 "cx_error"  – lowest 2-qubit gate error (best quality) [default]
                 "queue"     – shortest queue (fastest turnaround)
                 "qubits"    – most qubits (largest machine)
                 "combined"  – blended score: 70% quality + 30% availability

    Returns a JSON object with the ranking and a note about what it means.

    Note: fetching calibration data for every device takes ~10–30 seconds
    because it makes one API call per device.
    """
    service = _get_service()
    backends = service.backends()

    devices = []
    for backend in backends:
        status = backend.status()

        entry = {
            "name": backend.name,
            "num_qubits": backend.num_qubits,
            "pending_jobs": status.pending_jobs,
            "operational": status.operational,
        }

        # Fetch calibration data for any mode that needs error rates
        if sort_by in ("cx_error", "combined"):
            props = backend.properties()
            cx_errors = _cx_errors_for_backend(props)
            if cx_errors:
                entry["avg_cx_error"] = round(
                    sum(cx_errors) / len(cx_errors), 5
                )

        devices.append(entry)

    # Apply the requested sort
    if sort_by == "cx_error":
        # Ascending: lower error = better rank
        # Devices without calibration data fall to the end (inf sentinel)
        devices.sort(key=lambda d: d.get("avg_cx_error", float("inf")))

    elif sort_by == "queue":
        # Ascending: fewer pending jobs = shorter wait
        devices.sort(key=lambda d: d.get("pending_jobs", float("inf")))

    elif sort_by == "qubits":
        # Descending: more qubits = higher rank
        devices.sort(key=lambda d: d["num_qubits"], reverse=True)

    elif sort_by == "combined":
        # Blended score: 70% gate quality + 30% queue availability.
        #
        # Why min-max normalisation?
        # cx_error lives in ~[0.001, 0.05]; pending_jobs in ~[0, 500].
        # A raw sum would let queue dominate just because its numbers are
        # larger. Min-max rescales each metric to [0, 1] relative to the
        # current set of devices, so the 70/30 weights actually mean what
        # they say: quality matters more than speed, but both count.
        #
        # Why 70/30?
        # For research you care most about getting a correct result (low
        # error), but a 200-job queue means hours of waiting — so
        # availability gets a meaningful but smaller weight.

        cx_vals = [d["avg_cx_error"] for d in devices if d.get("avg_cx_error") is not None]
        q_vals  = [d["pending_jobs"]  for d in devices if d.get("pending_jobs")  is not None]

        min_cx, max_cx = (min(cx_vals), max(cx_vals)) if cx_vals else (0, 1)
        min_q,  max_q  = (min(q_vals),  max(q_vals))  if q_vals  else (0, 1)

        # Avoid division by zero when all devices have identical values
        cx_range = max_cx - min_cx or 1
        q_range  = max_q  - min_q  or 1

        for d in devices:
            cx = d.get("avg_cx_error")
            q  = d.get("pending_jobs")
            # Missing metrics get worst-case penalty (1.0) so uncalibrated
            # devices sort below any device with real data
            norm_cx = (cx - min_cx) / cx_range if cx is not None else 1.0
            norm_q  = (q  - min_q)  / q_range  if q  is not None else 1.0
            d["combined_score"] = round(0.7 * norm_cx + 0.3 * norm_q, 4)

        # Ascending: 0.0 = perfectly best on both metrics, 1.0 = worst
        devices.sort(key=lambda d: d.get("combined_score", float("inf")))

    else:
        return json.dumps({
            "error": f"Unknown sort_by value '{sort_by}'. "
                     "Use 'cx_error', 'queue', 'qubits', or 'combined'."
        })

    # Stamp each entry with its rank number (1 = best)
    for i, device in enumerate(devices):
        device["rank"] = i + 1

    _save_snapshots(devices)
    return json.dumps(
        {
            "sorted_by": sort_by,
            "note": {
                "cx_error": "Rank 1 = lowest 2-qubit gate error (highest quality)",
                "queue":    "Rank 1 = fewest pending jobs (shortest wait)",
                "qubits":   "Rank 1 = most qubits (largest machine)",
                "combined": "Rank 1 = best blend of quality (70%) and availability (30%). "
                            "Score is min-max normalised across current devices.",
            }.get(sort_by, ""),
            "devices": devices,
        },
        indent=2,
    )


# --------------------------------------------------------------------------
# Tool 4: queue_status
# --------------------------------------------------------------------------

@mcp.tool()
def queue_status() -> str:
    """
    Snapshot of the job queue on every IBM quantum computer.

    Useful when you want to submit a job and need to pick the machine
    with the shortest wait.

    Returns a JSON array sorted by pending_jobs (shortest queue first).
    """
    service = _get_service()
    backends = service.backends()

    queues = []
    for backend in backends:
        # status() is fast — it does NOT fetch full calibration data
        status = backend.status()
        queues.append({
            "name": backend.name,
            "num_qubits": backend.num_qubits,
            "pending_jobs": status.pending_jobs,
            "status": status.status_msg,
            "operational": status.operational,
        })

    # Shortest queue first so the "best pick right now" is at the top
    queues.sort(key=lambda d: d["pending_jobs"])

    _save_snapshots(queues)
    return json.dumps(queues, indent=2)


# --------------------------------------------------------------------------
# Tool 5: device_history
# --------------------------------------------------------------------------

@mcp.tool()
def device_history(device_name: str, days: int = 7) -> str:
    """
    Return all saved snapshots for one IBM quantum computer over the last N days.

    Args:
        device_name: Machine name, e.g. "ibm_brisbane". Must match exactly.
        days:        How many days back to look (default 7).

    Returns a JSON object with the device name and a list of snapshots in
    chronological order. Each snapshot has the fields that were available
    when it was recorded (error rates are NULL when the recording tool
    didn't fetch calibration data).
    """
    with sqlite3.connect(DB_PATH) as con:
        con.row_factory = sqlite3.Row
        rows = con.execute(
            """
            SELECT ts, num_qubits, operational, pending_jobs,
                   avg_cx_error, avg_readout_error
            FROM   device_snapshots
            WHERE  name = ?
              AND  ts >= datetime('now', ? || ' days')
            ORDER  BY ts ASC
            """,
            (device_name, f"-{days}"),
        ).fetchall()

    snapshots = [
        {
            "ts": r["ts"],
            "num_qubits": r["num_qubits"],
            "operational": bool(r["operational"]) if r["operational"] is not None else None,
            "pending_jobs": r["pending_jobs"],
            "avg_cx_error": r["avg_cx_error"],
            "avg_readout_error": r["avg_readout_error"],
        }
        for r in rows
    ]

    return json.dumps(
        {"device": device_name, "days": days, "snapshots": snapshots},
        indent=2,
    )


# --------------------------------------------------------------------------
# Tool 7: device_on_date
# --------------------------------------------------------------------------

@mcp.tool()
def device_on_date(device_name: str, date: str) -> str:
    """
    Historical stats for a device on a specific past date, from our snapshot DB.

    Useful for reproducibility: if you ran an experiment on 2026-07-01, call
    this tool with that date to see exactly what the hardware looked like —
    queue depth, error rates — and include it in your methods section.

    Args:
        device_name: Machine name, e.g. "ibm_fez".
        date:        Date in YYYY-MM-DD format, e.g. "2026-06-10".

    Returns aggregated stats averaged across all snapshots taken that day
    (snapshots are recorded every 6 hours by the background agent).
    """
    with sqlite3.connect(DB_PATH) as con:
        con.row_factory = sqlite3.Row
        rows = con.execute(
            """
            SELECT ts, operational, pending_jobs, avg_cx_error, avg_readout_error
            FROM   device_snapshots
            WHERE  name   = ?
              AND  date(ts) = ?
            ORDER  BY ts ASC
            """,
            (device_name, date),
        ).fetchall()

    if not rows:
        return json.dumps({
            "device": device_name,
            "date":   date,
            "found":  False,
            "note":   "No snapshots found for this device on this date. "
                      "Snapshots are recorded every 6 hours by the local LaunchAgent.",
        })

    snapshots = [dict(r) for r in rows]

    def _avg(field: str):
        vals = [s[field] for s in snapshots if s[field] is not None]
        return round(sum(vals) / len(vals), 5) if vals else None

    return json.dumps(
        {
            "device":             device_name,
            "date":               date,
            "found":              True,
            "snapshots_that_day": len(snapshots),
            "first_snapshot":     snapshots[0]["ts"],
            "last_snapshot":      snapshots[-1]["ts"],
            "avg_pending_jobs":   _avg("pending_jobs"),
            "avg_cx_error":       _avg("avg_cx_error"),
            "avg_readout_error":  _avg("avg_readout_error"),
            "note": "Averaged across all snapshots taken that day. "
                    "Cite this date in your paper's methods section for reproducibility.",
        },
        indent=2,
    )


# --------------------------------------------------------------------------
# Tool 8: submit_job
# --------------------------------------------------------------------------

@mcp.tool()
def submit_job(device_name: str, qasm_string: str, shots: int = 1024,
               qasm_version: int = 2) -> str:
    """
    Compile and submit a quantum circuit to an IBM quantum computer.

    Args:
        device_name:  Machine name, e.g. "ibm_fez". Use compare_devices or
                      queue_status first to pick the best available machine.
        qasm_string:  OpenQASM circuit source code.
                      For QASM 2: must start with OPENQASM 2.0; include "qelib1.inc";
                      For QASM 3: must start with OPENQASM 3.0;
                      The circuit must include measurement gates.
        shots:        How many times to run the circuit (default 1024, max 20000).
                      More shots = more accurate probability estimates.
        qasm_version: 2 (default) for OpenQASM 2.0, 3 for OpenQASM 3.0.

    Returns JSON with:
      - job_id   Save this — needed for job_status and job_results.
      - status   Initial status (usually "INITIALIZING" or "QUEUED").
      - device   Machine the job was sent to.
      - shots    Number of shots requested.
    """
    # Parse the QASM string into a Qiskit QuantumCircuit object.
    # QASM 2 uses QuantumCircuit.from_qasm_str (legacy standard).
    # QASM 3 uses qiskit.qasm3.loads (modern standard with richer features).
    try:
        if qasm_version == 3:
            circuit = qiskit_qasm3.loads(qasm_string)
        else:
            circuit = QuantumCircuit.from_qasm_str(qasm_string)
    except Exception as e:
        return json.dumps({
            "error": f"Failed to parse QASM {qasm_version}: {e}",
            "hint": (
                'QASM 2 must start with: OPENQASM 2.0;\ninclude "qelib1.inc";'
                if qasm_version != 3 else
                "QASM 3 must start with: OPENQASM 3.0;"
            ),
        })

    # Clamp shots to IBM's allowed range
    shots = max(1, min(shots, 20000))

    service = _get_service()

    try:
        backend = service.backend(device_name)
    except Exception as e:
        return json.dumps({"error": f"Device '{device_name}' not found: {e}"})

    # Transpile the circuit to the backend's native gate set and qubit topology.
    # optimization_level=1 is a good default: fast compile, decent optimisation.
    # Level 3 gives the best circuit but is much slower to compile.
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(circuit)

    # SamplerV2 is the current IBM Runtime primitive for sampling circuits.
    # It replaces the deprecated execute() function and Sampler v1.
    # mode=backend tells it which machine to target.
    sampler = Sampler(mode=backend)
    job = sampler.run([isa_circuit], shots=shots)

    return json.dumps({
        "job_id": job.job_id(),
        "status": str(job.status()),
        "device": device_name,
        "shots": shots,
        "note": "Save the job_id. Use job_status to check progress, job_results to get counts.",
    }, indent=2)


# --------------------------------------------------------------------------
# Tool 9: job_status
# --------------------------------------------------------------------------

@mcp.tool()
def job_status(job_id: str) -> str:
    """
    Check the current status of a submitted quantum job.

    Args:
        job_id: The ID returned by submit_job.

    Status values:
      INITIALIZING  Just submitted, not yet in the queue.
      QUEUED        Waiting in the machine's queue.
      RUNNING       Actively executing on hardware right now.
      DONE          Finished — call job_results to get counts.
      ERROR         Failed — error_message field will explain why.
      CANCELLED     Was cancelled before it ran.
    """
    service = _get_service()

    try:
        job = service.job(job_id)
    except Exception as e:
        return json.dumps({"error": f"Job '{job_id}' not found: {e}"})

    status = str(job.status())
    result = {
        "job_id":  job_id,
        "status":  status,
        "backend": job.backend().name,
    }

    try:
        result["creation_date"] = str(job.creation_date)
    except Exception:
        pass

    if status == "QUEUED":
        # queue_position() tells you how many jobs are ahead of yours
        try:
            pos = job.queue_position()
            if pos is not None:
                result["queue_position"] = pos
        except Exception:
            pass
        result["note"] = "Still waiting in queue. Check again in a few minutes."

    elif status == "DONE":
        result["note"] = "Job complete — call job_results to retrieve counts."

    elif status == "ERROR":
        try:
            result["error_message"] = job.error_message()
        except Exception:
            pass

    return json.dumps(result, indent=2)


# --------------------------------------------------------------------------
# Tool 10: job_results
# --------------------------------------------------------------------------

@mcp.tool()
def job_results(job_id: str) -> str:
    """
    Retrieve measurement counts from a completed quantum job.

    Args:
        job_id: The ID returned by submit_job.

    Returns JSON with bit-string counts when the job is DONE.
    If the job is still running or queued, returns current status instead.

    Counts example:
      {"00": 502, "11": 522}  ← a Bell state: roughly 50/50 between 00 and 11.

    The bit-string length equals the number of measured qubits.
    Each key is a measurement outcome; the value is how many shots produced it.
    All values sum to the total number of shots.
    """
    service = _get_service()

    try:
        job = service.job(job_id)
    except Exception as e:
        return json.dumps({"error": f"Job '{job_id}' not found: {e}"})

    status = str(job.status())

    if status != "DONE":
        return json.dumps({
            "job_id": job_id,
            "status": status,
            "note":   "Job not complete yet. Use job_status to monitor progress.",
        }, indent=2)

    try:
        result = job.result()
    except Exception as e:
        return json.dumps({"error": f"Failed to retrieve results: {e}"})

    # SamplerV2 wraps results in a PrimitiveResult containing one PubResult per circuit.
    # Each PubResult has a DataBin with one BitArray per classical register.
    # We collect counts from every register (circuits with one register are the common case).
    try:
        pub_result = result[0]
        counts_by_register = {}
        for reg_name, bit_array in vars(pub_result.data).items():
            counts_by_register[reg_name] = bit_array.get_counts()
    except Exception as e:
        return json.dumps({
            "error": f"Failed to parse result data: {e}",
            "raw_result": str(result),
        })

    # Flatten to a single counts dict when there is only one register (usual case)
    counts = (
        list(counts_by_register.values())[0]
        if len(counts_by_register) == 1
        else counts_by_register
    )

    total_shots = sum(counts.values()) if isinstance(counts, dict) else None

    return json.dumps({
        "job_id":      job_id,
        "status":      "DONE",
        "backend":     job.backend().name,
        "total_shots": total_shots,
        "counts":      counts,
        "note": "Each key is a bit-string outcome; value is how many shots produced it.",
    }, indent=2)


# --------------------------------------------------------------------------
# Tool 11: cancel_job
# --------------------------------------------------------------------------

@mcp.tool()
def cancel_job(job_id: str) -> str:
    """
    Cancel a queued or running IBM quantum job.

    Only jobs in QUEUED or RUNNING state can be cancelled. Jobs that are
    already DONE, ERROR, or CANCELLED will return an error.

    Args:
        job_id: The ID returned by submit_job or list_jobs.

    Returns JSON confirming the cancellation or explaining why it failed.
    """
    service = _get_service()

    try:
        job = service.job(job_id)
    except Exception as e:
        return json.dumps({"error": f"Job '{job_id}' not found: {e}"})

    status = str(job.status())

    # IBM only allows cancellation before the job finishes
    if status in ("DONE", "ERROR", "CANCELLED"):
        return json.dumps({
            "job_id": job_id,
            "error": f"Cannot cancel a job with status '{status}'.",
            "current_status": status,
        })

    try:
        job.cancel()
    except Exception as e:
        return json.dumps({"error": f"Cancel request failed: {e}"})

    return json.dumps({
        "job_id": job_id,
        "status": "CANCELLED",
        "note": "Cancellation requested. The job may take a moment to fully stop.",
    }, indent=2)


# --------------------------------------------------------------------------
# Tool 12: list_jobs
# --------------------------------------------------------------------------

@mcp.tool()
def list_jobs(limit: int = 10) -> str:
    """
    List your most recently submitted IBM quantum jobs.

    Useful for finding job IDs you didn't save, or getting an overview of
    what's in the queue right now.

    Args:
        limit: How many jobs to return, newest first (default 10, max 50).

    Returns a JSON array of jobs, each with: job_id, status, backend, creation date.
    """
    limit = max(1, min(limit, 50))  # clamp to [1, 50]

    service = _get_service()

    try:
        jobs = service.jobs(limit=limit)
    except Exception as e:
        return json.dumps({"error": f"Failed to fetch jobs: {e}"})

    results = []
    for job in jobs:
        entry = {
            "job_id": job.job_id(),
            "status": str(job.status()),
        }
        # backend() and creation_date can raise on malformed jobs — guard each
        try:
            entry["backend"] = job.backend().name
        except Exception:
            entry["backend"] = None
        try:
            entry["created"] = str(job.creation_date)
        except Exception:
            entry["created"] = None

        results.append(entry)

    return json.dumps({
        "count": len(results),
        "jobs": results,
        "note": "Sorted newest first. Use job_status(job_id) for full details.",
    }, indent=2)


# --------------------------------------------------------------------------
# Tool 13: run_grover
# --------------------------------------------------------------------------

@mcp.tool()
def run_grover(n_qubits: int, target_state: str) -> str:
    """
    Build and run a Grover's search algorithm demo on real IBM hardware.

    Grover's algorithm searches an unsorted list of 2^n states in O(sqrt(2^n))
    steps — a quadratic speedup over classical search. This tool builds the
    full circuit (oracle + diffusion operator), transpiles it, and submits it
    to the least-busy backend.

    Args:
        n_qubits:     Number of qubits to use. Must be 2 or 3.
                      Capped at 3 — deeper circuits lose coherence on current
                      hardware and the result becomes dominated by noise.
        target_state: Binary string marking the state to find, e.g. "11" or "101".
                      Length must equal n_qubits. Grover's will amplify this state's
                      probability so it appears far more often than the others.

    Returns JSON with the job_id, backend chosen, circuit details, and what
    fraction of shots to expect on the target state.
    """
    # Validate and clamp inputs
    n_qubits = min(max(n_qubits, 2), 3)

    if len(target_state) != n_qubits:
        return json.dumps({
            "error": (
                f"target_state length ({len(target_state)}) must equal "
                f"n_qubits ({n_qubits}). Example: n_qubits=2, target_state='11'"
            )
        })

    if not all(c in "01" for c in target_state):
        return json.dumps({"error": "target_state must contain only '0' and '1'."})

    # Optimal number of Grover iterations for a single marked state:
    # floor(π/4 * sqrt(N)) where N = 2^n_qubits
    # n=2 → 1 iteration, n=3 → 2 iterations
    n_iterations = max(1, round(math.pi / 4 * math.sqrt(2 ** n_qubits)))

    qc = QuantumCircuit(n_qubits, n_qubits)

    # Step 1: put all qubits in equal superposition (uniform over all 2^n states)
    qc.h(range(n_qubits))

    for _ in range(n_iterations):
        # ── Oracle: phase-flip the target state ───────────────────────────
        # Strategy: X-gate every qubit whose target bit is '0', so the
        # target state maps to all-|1⟩, apply a multi-controlled-Z to flip
        # its phase, then undo the X gates.
        # reversed() because Qiskit is little-endian (qubit 0 = rightmost bit).
        for i, bit in enumerate(reversed(target_state)):
            if bit == "0":
                qc.x(i)

        if n_qubits == 2:
            qc.cz(0, 1)
        else:  # n_qubits == 3
            qc.ccz(0, 1, 2)

        for i, bit in enumerate(reversed(target_state)):
            if bit == "0":
                qc.x(i)

        # ── Diffusion operator: inversion about the mean ───────────────────
        # This amplifies the target state's amplitude at the cost of the others.
        # Circuit: H⊗n → X⊗n → multi-CZ → X⊗n → H⊗n
        qc.h(range(n_qubits))
        qc.x(range(n_qubits))

        if n_qubits == 2:
            qc.cz(0, 1)
        else:
            qc.ccz(0, 1, 2)

        qc.x(range(n_qubits))
        qc.h(range(n_qubits))

    # Measure all qubits
    qc.measure(range(n_qubits), range(n_qubits))

    # Find the least-busy operational backend
    service = _get_service()
    backends = service.backends()

    operational = []
    for b in backends:
        try:
            s = b.status()
            if s.operational:
                operational.append((b, s.pending_jobs))
        except Exception:
            pass

    if not operational:
        return json.dumps({"error": "No operational backends available."})

    best_backend, _ = min(operational, key=lambda x: x[1])

    # Transpile to the backend's native gate set and submit
    pm = generate_preset_pass_manager(backend=best_backend, optimization_level=1)
    isa_circuit = pm.run(qc)

    sampler = Sampler(mode=best_backend)
    job = sampler.run([isa_circuit], shots=1024)

    # Theoretical success probability after optimal iterations
    # P = sin²((2k+1) * arcsin(1/sqrt(N))) where k = n_iterations
    theta = math.asin(1 / math.sqrt(2 ** n_qubits))
    ideal_pct = round(100 * math.sin((2 * n_iterations + 1) * theta) ** 2, 1)

    return json.dumps({
        "job_id": job.job_id(),
        "status": str(job.status()),
        "device": best_backend.name,
        "n_qubits": n_qubits,
        "target_state": target_state,
        "grover_iterations": n_iterations,
        "shots": 1024,
        "ideal_success_pct": ideal_pct,
        "note": (
            f"Searching for |{target_state}⟩ across {2**n_qubits} states. "
            f"Ideal hardware would show '{target_state}' in {ideal_pct}% of shots. "
            f"Real hardware noise will reduce this — expect 60–85% on current IBM devices. "
            f"Use job_status then job_results to see the counts."
        ),
    }, indent=2)


# --------------------------------------------------------------------------
# Tool 14: estimate_expectation
# --------------------------------------------------------------------------

@mcp.tool()
def estimate_expectation(device_name: str, qasm_string: str,
                         observables: str, shots: int = 1024,
                         qasm_version: int = 2) -> str:
    """
    Run the Estimator primitive to compute the expectation value of one or
    more observables for a parameterised quantum state.

    Unlike submit_job (which counts measurement outcomes), the Estimator
    computes <ψ|O|ψ> — the average value of an observable O. This is what
    quantum chemistry and optimisation algorithms (VQE, QAOA) need.

    Args:
        device_name:  IBM backend to run on, e.g. "ibm_fez".
        qasm_string:  Circuit that prepares the quantum state (no measurements
                      needed — Estimator handles that internally).
        observables:  Comma-separated Pauli strings, e.g. "ZZ,XI,IZ".
                      Each string is a tensor product of single-qubit Paulis
                      (I, X, Y, Z) — length must equal the number of qubits.
        shots:        Shots per observable (default 1024, max 20000).
        qasm_version: 2 (default) or 3.

    Returns JSON with:
      - job_id       Use job_status to track, job_results won't work — check
                     status via job_status and retrieve via this tool's job_id.
      - observables  List of Pauli strings submitted.
      - device       Backend used.
      - note         Explanation of expectation values.
    """
    # Parse circuit (same logic as submit_job)
    try:
        if qasm_version == 3:
            circuit = qiskit_qasm3.loads(qasm_string)
        else:
            circuit = QuantumCircuit.from_qasm_str(qasm_string)
    except Exception as e:
        return json.dumps({"error": f"Failed to parse QASM {qasm_version}: {e}"})

    # Parse comma-separated Pauli strings, e.g. "ZZ,XI" → ["ZZ", "XI"]
    pauli_list = [p.strip().upper() for p in observables.split(",") if p.strip()]
    if not pauli_list:
        return json.dumps({"error": "observables must be a comma-separated list of Pauli strings, e.g. 'ZZ,XI'"})

    service = _get_service()
    try:
        backend = service.backend(device_name)
    except Exception as e:
        return json.dumps({"error": f"Device '{device_name}' not found: {e}"})

    shots = max(1, min(shots, 20000))

    # Transpile the circuit to the backend's native gate set.
    # Estimator requires an ISA circuit (Instruction Set Architecture).
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(circuit)

    # Build SparsePauliOp objects for each observable.
    # SparsePauliOp("ZZ") represents the Z⊗Z operator on 2 qubits.
    try:
        ops = [SparsePauliOp(p) for p in pauli_list]
    except Exception as e:
        return json.dumps({"error": f"Invalid Pauli string: {e}. Use I, X, Y, Z only."})

    # EstimatorV2 takes (circuit, observable) pairs called "PUBs"
    # (Primitive Unified Blocs). One PUB per observable.
    estimator = Estimator(mode=backend)
    estimator.options.default_shots = shots
    pubs = [(isa_circuit, op) for op in ops]

    try:
        job = estimator.run(pubs)
    except Exception as e:
        return json.dumps({"error": f"Estimator submission failed: {e}"})

    return json.dumps({
        "job_id": job.job_id(),
        "status": str(job.status()),
        "device": device_name,
        "observables": pauli_list,
        "shots": shots,
        "note": (
            "Use job_status to check progress. When DONE, retrieve results with "
            "job_results — expectation values will be in the 'values' field. "
            "Each value is a float in [-1, +1]: +1 means all qubits measured the "
            "operator's +1 eigenstate, -1 means the -1 eigenstate."
        ),
    }, indent=2)


# --------------------------------------------------------------------------
# Tool 15: circuit_report
# --------------------------------------------------------------------------

@mcp.tool()
def circuit_report(device_name: str, qasm_string: str,
                   qasm_version: int = 2) -> str:
    """
    Dry-run analysis of a circuit on a specific backend — no job submitted,
    no queue time, instant results.

    This is the "look before you leap" tool. Before waiting hours in a queue,
    use circuit_report to see:
      - How the compiler transforms your circuit (gate count, depth)
      - Which physical qubits get assigned to your logical qubits
      - The error rate on each assigned qubit pair
      - An estimated fidelity — the probability your result is correct

    Researchers use this to:
      - Compare backends before committing to one
      - Detect if the compiler is bloating their circuit
      - Know in advance if today's calibration is good enough

    Args:
        device_name:  Backend to analyse against, e.g. "ibm_fez".
        qasm_string:  Circuit in OpenQASM format (measurements optional).
        qasm_version: 2 (default) or 3.

    Returns JSON with:
      - original_gates     Gate counts before transpilation
      - transpiled_gates   Gate counts after IBM compiler (usually more gates)
      - original_depth     Circuit depth before compilation
      - transpiled_depth   Circuit depth after compilation
      - qubit_mapping      Logical qubit → physical qubit assignment
      - cx_error_per_pair  2-qubit gate error on each used qubit pair
      - estimated_fidelity Probability the circuit produces the correct result
      - verdict            Human-readable recommendation
    """
    # Parse the circuit
    try:
        if qasm_version == 3:
            circuit = qiskit_qasm3.loads(qasm_string)
        else:
            circuit = QuantumCircuit.from_qasm_str(qasm_string)
    except Exception as e:
        return json.dumps({"error": f"Failed to parse QASM {qasm_version}: {e}"})

    service = _get_service()
    try:
        backend = service.backend(device_name)
    except Exception as e:
        return json.dumps({"error": f"Device '{device_name}' not found: {e}"})

    # Transpile — this is the same step submit_job does, but we stop before
    # actually running anything.
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    try:
        isa_circuit = pm.run(circuit)
    except Exception as e:
        return json.dumps({"error": f"Transpilation failed: {e}"})

    # Gate counts before and after compilation
    original_gates = dict(circuit.count_ops())
    transpiled_gates = dict(isa_circuit.count_ops())
    original_depth = circuit.depth()
    transpiled_depth = isa_circuit.depth()

    # Extract qubit layout: logical index → physical qubit index
    layout = isa_circuit.layout
    qubit_mapping = {}
    if layout and layout.final_layout:
        for logical, physical in enumerate(layout.final_layout):
            qubit_mapping[f"q{logical}"] = int(str(physical).split("_")[-1]) if "_" in str(physical) else logical
    elif layout and layout.initial_layout:
        for logical_bit, physical_bit in layout.initial_layout.get_physical_bits().items():
            if hasattr(physical_bit, "index"):
                qubit_mapping[f"q{physical_bit.index}"] = logical_bit

    # Pull 2-qubit gate errors from backend calibration for the used qubits.
    # This tells you whether the assigned qubit pairs are in good shape today.
    cx_errors = {}
    try:
        props = backend.properties()
        if props:
            used_indices = list(qubit_mapping.values()) if qubit_mapping else list(range(circuit.num_qubits))
            for gate in props.gates:
                if gate.gate in ("cx", "ecr", "cz") and len(gate.qubits) == 2:
                    q0, q1 = gate.qubits
                    if q0 in used_indices or q1 in used_indices:
                        for param in gate.parameters:
                            if param.name == "gate_error":
                                cx_errors[f"q{q0}-q{q1}"] = round(param.value, 6)
    except Exception:
        pass  # calibration data unavailable — report without it

    # Estimate circuit fidelity using the product-of-gate-errors model:
    # fidelity ≈ ∏(1 - error_i) for each 2-qubit gate in the transpiled circuit.
    # This is a lower bound — real fidelity is often better due to error correlation.
    n_cx = transpiled_gates.get("cx", 0) + transpiled_gates.get("ecr", 0) + transpiled_gates.get("cz", 0)
    avg_cx_error = sum(cx_errors.values()) / len(cx_errors) if cx_errors else 0.005
    estimated_fidelity = round((1 - avg_cx_error) ** n_cx, 4) if n_cx > 0 else 1.0

    # Plain-English verdict based on fidelity
    if estimated_fidelity >= 0.90:
        verdict = "Excellent — this circuit should produce clean results on this backend today."
    elif estimated_fidelity >= 0.70:
        verdict = "Good — expect some noise but results should be meaningful."
    elif estimated_fidelity >= 0.50:
        verdict = "Fair — significant noise expected. Consider a lower-error backend or fewer gates."
    else:
        verdict = "Poor — high noise likely to obscure results. Try compare_devices to find a better backend."

    # Overhead: how much did the compiler bloat the circuit?
    overhead = round(
        (sum(transpiled_gates.values()) - sum(original_gates.values()))
        / max(sum(original_gates.values()), 1) * 100, 1
    )

    return json.dumps({
        "device": device_name,
        "original_gates": original_gates,
        "transpiled_gates": transpiled_gates,
        "original_depth": original_depth,
        "transpiled_depth": transpiled_depth,
        "compiler_overhead_pct": overhead,
        "qubit_mapping": qubit_mapping,
        "cx_error_per_pair": cx_errors,
        "estimated_fidelity": estimated_fidelity,
        "n_two_qubit_gates": n_cx,
        "verdict": verdict,
    }, indent=2)


# --------------------------------------------------------------------------
# API Key Authentication Middleware
# --------------------------------------------------------------------------

class APIKeyAuthMiddleware(BaseHTTPMiddleware):
    """
    Middleware to validate API key for HTTP requests.
    
    If MCP_API_KEY environment variable is set, all requests must include
    a matching X-API-Key header. If not set, all requests are allowed
    (development mode).
    """
    
    def __init__(self, app, api_key: Optional[str] = None):
        super().__init__(app)
        self.api_key = api_key or os.getenv("MCP_API_KEY")
    
    async def dispatch(self, request: Request, call_next):
        # If no API key is configured, allow all requests (development mode)
        if not self.api_key:
            return await call_next(request)
        
        # Check for API key in request headers (case-insensitive)
        request_key = request.headers.get("X-API-Key") or request.headers.get("x-api-key")
        
        if request_key != self.api_key:
            return JSONResponse(
                status_code=401,
                content={
                    "error": "Unauthorized",
                    "message": "Invalid or missing API key. Include X-API-Key header with your request.",
                }
            )
        
        return await call_next(request)


# --------------------------------------------------------------------------
# Command-Line Argument Parsing
# --------------------------------------------------------------------------

def parse_args():
    """Parse command-line arguments for transport configuration."""
    parser = argparse.ArgumentParser(
        description="Quantum Hardware MCP Server - Exposes IBM Quantum device data via MCP protocol",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # stdio mode (default, for Claude Desktop)
  python server.py
  
  # HTTP mode on localhost
  python server.py --transport http
  
  # HTTP mode on all interfaces with custom port
  python server.py --transport http --host 0.0.0.0 --port 8080
  
  # HTTP mode with specific CORS origins
  python server.py --transport http --cors-origins "https://myapp.com,https://api.myapp.com"

Environment Variables:
  IBM_QUANTUM_TOKEN   IBM Quantum API token (required)
  MCP_HTTP_HOST       HTTP server host (default: 127.0.0.1)
  MCP_HTTP_PORT       HTTP server port (default: 8000)
  MCP_CORS_ORIGINS    Comma-separated CORS origins (default: *)
  MCP_API_KEY         API key for authentication (optional, recommended for production)
        """
    )
    
    parser.add_argument(
        "--transport",
        choices=["stdio", "http"],
        default="stdio",
        help="Transport mode: 'stdio' for Claude Desktop (default), 'http' for remote clients"
    )
    
    parser.add_argument(
        "--host",
        default=os.getenv("MCP_HTTP_HOST", "127.0.0.1"),
        help="HTTP server host (default: 127.0.0.1, use 0.0.0.0 for all interfaces)"
    )
    
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.getenv("MCP_HTTP_PORT", "8000")),
        help="HTTP server port (default: 8000)"
    )
    
    parser.add_argument(
        "--cors-origins",
        default=os.getenv("MCP_CORS_ORIGINS", "*"),
        help="Comma-separated CORS origins (default: *, use specific domains in production)"
    )
    
    return parser.parse_args()


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    
    if args.transport == "stdio":
        # stdio transport for Claude Desktop integration
        # Claude Desktop launches this process and communicates over stdin/stdout
        # Note: Cannot use print() here as it would corrupt the JSON-RPC protocol stream
        mcp.run(transport="stdio")
    
    elif args.transport == "http":
        # HTTP/SSE transport for remote MCP clients
        # Enables web-based AI assistants and remote integrations
        
        # Configure server settings
        mcp.settings.host = args.host
        mcp.settings.port = args.port
        
        # Check if API key is configured
        api_key = os.getenv("MCP_API_KEY")
        api_key_configured = bool(api_key)
        
        print("=" * 70, flush=True)
        print("Quantum Hardware MCP Server - HTTP Mode", flush=True)
        print("=" * 70, flush=True)
        print(f"Server URL:    http://{args.host}:{args.port}", flush=True)
        print(f"CORS Origins:  {args.cors_origins}", flush=True)
        print(f"Authentication: {'Enabled (API key required)' if api_key_configured else 'Disabled (development mode)'}", flush=True)
        
        if not api_key_configured:
            print("\n⚠️  WARNING: No API key configured!", flush=True)
            print("   Set MCP_API_KEY environment variable for production use.", flush=True)
            print("   Generate a key with: python -c \"import secrets; print(secrets.token_urlsafe(32))\"", flush=True)
        
        print("=" * 70, flush=True)
        print("\nServer starting...\n", flush=True)
        
        # Add API key authentication middleware to the Starlette app
        # This must be done before calling run() to ensure middleware is applied
        async def run_http_with_auth():
            """Run HTTP server with authentication middleware."""
            starlette_app = mcp.sse_app()
            starlette_app.add_middleware(APIKeyAuthMiddleware, api_key=api_key)

            # The MCP SDK (transport_security.py) validates the Host header against
            # the pattern "localhost:*" — it accepts any "localhost:PORT" but NOT
            # bare "localhost". In Docker the agent uses "mcp-server:3020" as the
            # host name, which the SDK rejects with 421. We rewrite it to
            # "localhost:{port}" before the SDK sees it.
            _port = mcp.settings.port
            _host_value = f"localhost:{_port}".encode()

            class DockerHostFix:
                def __init__(self, app):
                    self.app = app

                async def __call__(self, scope, receive, send):
                    if scope["type"] in ("http", "websocket"):
                        scope = {**scope, "headers": [
                            (b"host", _host_value) if k == b"host" else (k, v)
                            for k, v in scope.get("headers", [])
                        ]}
                    await self.app(scope, receive, send)

            import uvicorn
            config = uvicorn.Config(
                DockerHostFix(starlette_app),
                host=mcp.settings.host,
                port=mcp.settings.port,
                log_level=mcp.settings.log_level.lower(),
            )
            server = uvicorn.Server(config)
            await server.serve()
        
        # Run the server with authentication
        anyio.run(run_http_with_auth)
