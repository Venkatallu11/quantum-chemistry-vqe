"""
Quantum Hardware MCP Server
============================
Exposes live IBM Quantum device data to AI assistants via the MCP protocol.

Tools:
  - list_devices       : all machines + status
  - get_device_details : deep info on one machine
  - compare_devices    : rank machines by error rate / queue
  - queue_status       : current queue depth for every machine
  - device_history     : historical snapshots for one machine over N days
"""

import os
import json
import sqlite3
from datetime import datetime, timezone

from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP
from qiskit_ibm_runtime import QiskitRuntimeService

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
    Pull all CX (CNOT) gate error values from calibration properties.
    CX is the most important 2-qubit gate; its error rate is the best
    single-number summary of device quality.

    Returns an empty list if the backend has no calibration data.
    """
    if props is None:
        return []
    errors = []
    for gate in props.gates:
        # Only look at CX gates
        if gate.gate == "cx" and gate.parameters:
            # parameters[0] is always gate_error for CX gates
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

        # Only fetch calibration data when we actually need error rates
        if sort_by == "cx_error":
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

    else:
        return json.dumps({
            "error": f"Unknown sort_by value '{sort_by}'. "
                     "Use 'cx_error', 'queue', or 'qubits'."
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
# Entry point
# --------------------------------------------------------------------------

if __name__ == "__main__":
    # stdio transport is required for Claude Desktop integration.
    # Claude Desktop launches this process and communicates over stdin/stdout.
    mcp.run(transport="stdio")
