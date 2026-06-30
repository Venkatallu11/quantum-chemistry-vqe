"""
snapshot.py
-----------
Standalone script that fetches live stats for every accessible IBM Quantum
device and writes a timestamped row to devices.db.

Run manually:
    .venv/bin/python snapshot.py

Or let the LaunchAgent call it automatically every 6 hours.
It does NOT start the MCP server — it only touches the database.
"""

import os
import sys
import json
import sqlite3
from datetime import datetime, timezone

from dotenv import load_dotenv
from qiskit_ibm_runtime import QiskitRuntimeService

# Load .env from the same directory as this file.
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), ".env"))

DB_PATH = os.path.join(os.path.dirname(__file__), "devices.db")


def _init_db() -> None:
    with sqlite3.connect(DB_PATH) as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS device_snapshots (
                id                INTEGER PRIMARY KEY AUTOINCREMENT,
                ts                TEXT    NOT NULL,
                name              TEXT    NOT NULL,
                num_qubits        INTEGER,
                operational       INTEGER,
                pending_jobs      INTEGER,
                avg_cx_error      REAL,
                avg_readout_error REAL
            )
        """)
        con.execute("""
            CREATE INDEX IF NOT EXISTS idx_name_ts
            ON device_snapshots (name, ts)
        """)


def _save_snapshots(rows: list[dict]) -> None:
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


def _cx_errors(props) -> list[float]:
    if props is None:
        return []
    return [
        g.parameters[0].value
        for g in props.gates
        if g.gate == "cx" and g.parameters
    ]


def collect() -> None:
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not set in .env", file=sys.stderr)
        sys.exit(1)

    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backends = service.backends()

    rows = []
    for backend in backends:
        status = backend.status()
        props = backend.properties()

        row = {
            "name": backend.name,
            "num_qubits": backend.num_qubits,
            "operational": status.operational,
            "pending_jobs": status.pending_jobs,
        }

        # Collect error rates while we're here — richer data than list_devices.
        if props:
            cx = _cx_errors(props)
            if cx:
                row["avg_cx_error"] = round(sum(cx) / len(cx), 5)

            readout = [
                props.readout_error(q)
                for q in range(backend.num_qubits)
                if props.readout_error(q) is not None
            ]
            if readout:
                row["avg_readout_error"] = round(sum(readout) / len(readout), 5)

        rows.append(row)

    _save_snapshots(rows)
    print(
        f"[{datetime.now(timezone.utc).isoformat()}] "
        f"Saved {len(rows)} device snapshots to {DB_PATH}"
    )


if __name__ == "__main__":
    _init_db()
    collect()
