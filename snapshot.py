"""
snapshot.py
-----------
Standalone script that fetches live stats for every accessible IBM Quantum
device and persists them.

Where it writes depends on the environment:
  - Locally (LaunchAgent):  SQLite devices.db  → feeds device_history + report.py
  - GitHub Actions (CI):    CSV data/snapshots.csv → committed to the repo as history

Run manually:
    .venv/bin/python snapshot.py

Or let the LaunchAgent / GitHub Actions call it automatically every 6 hours.
"""

import os
import sys
import csv
import sqlite3
from datetime import datetime, timezone

from dotenv import load_dotenv
from qiskit_ibm_runtime import QiskitRuntimeService

# Load .env from the same directory as this file.
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), ".env"))

BASE_DIR  = os.path.dirname(__file__)
DB_PATH   = os.path.join(BASE_DIR, "devices.db")
CSV_PATH  = os.path.join(BASE_DIR, "data", "snapshots.csv")

CSV_FIELDS = ["ts", "name", "num_qubits", "operational",
              "pending_jobs", "avg_cx_error", "avg_readout_error"]


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


def _write_csv(rows: list[dict]) -> None:
    """
    Append snapshot rows to data/snapshots.csv.
    Creates the file with a header row on the first call.
    Used by GitHub Actions so the history is committed as plain text.
    """
    os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)
    ts = datetime.now(timezone.utc).isoformat()
    new_file = not os.path.exists(CSV_PATH)

    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        for r in rows:
            writer.writerow({
                "ts":               ts,
                "name":             r["name"],
                "num_qubits":       r.get("num_qubits"),
                "operational":      r.get("operational"),
                "pending_jobs":     r.get("pending_jobs"),
                "avg_cx_error":     r.get("avg_cx_error"),
                "avg_readout_error": r.get("avg_readout_error"),
            })


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

    if os.getenv("GITHUB_ACTIONS"):
        # CI: write CSV — SQLite isn't persisted between Actions runs anyway
        _write_csv(rows)
        print(f"[{datetime.now(timezone.utc).isoformat()}] "
              f"Wrote {len(rows)} rows to {CSV_PATH}")
    else:
        # Local: write SQLite — feeds device_history MCP tool and report.py
        _save_snapshots(rows)
        print(f"[{datetime.now(timezone.utc).isoformat()}] "
              f"Saved {len(rows)} snapshots to {DB_PATH}")


if __name__ == "__main__":
    if not os.getenv("GITHUB_ACTIONS"):
        _init_db()
    collect()
