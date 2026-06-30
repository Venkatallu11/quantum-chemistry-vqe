"""
Quantum Hardware MCP Server
============================
Exposes live IBM Quantum device data to AI assistants via the MCP protocol.

Tools:
  - list_devices       : all machines + status
  - get_device_details : deep info on one machine
  - compare_devices    : rank machines by error rate / queue / combined score
  - queue_status       : current queue depth for every machine
  - device_history     : historical snapshots for one machine over N days
  - best_qubits        : best n qubits on a machine right now (calibration-based)
  - device_on_date     : historical stats for a machine on a specific past date
"""

import os
import json
import sqlite3
import argparse
import anyio
from datetime import datetime, timezone
from typing import Optional

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
            # Get the Starlette app
            starlette_app = mcp.sse_app()
            
            # Add authentication middleware
            starlette_app.add_middleware(APIKeyAuthMiddleware, api_key=api_key)
            
            # Start the server
            import uvicorn
            config = uvicorn.Config(
                starlette_app,
                host=mcp.settings.host,
                port=mcp.settings.port,
                log_level=mcp.settings.log_level.lower(),
            )
            server = uvicorn.Server(config)
            await server.serve()
        
        # Run the server with authentication
        anyio.run(run_http_with_auth)
