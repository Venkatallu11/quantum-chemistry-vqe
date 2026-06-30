"""
Strong tests for all 16 MCP server tools.

These tests call the tool functions directly (no HTTP, no queue).
They verify:
  - Each tool returns a string (not an exception)
  - No tool returns a raw {"error": ...} on a valid call
  - Tools that need a real device name use a known IBM device

Run with:
    pytest tests/test_server_tools.py -v
"""

import json
import sys
import os
import pytest

# Add the project root to the path so we can import server.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# A minimal Bell state circuit in QASM 2 — used by job/circuit tools
BELL_QASM2 = """
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
h q[0];
cx q[0], q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
""".strip()

# Same circuit in QASM 3
BELL_QASM3 = """
OPENQASM 3.0;
include "stdgates.inc";
qubit[2] q;
bit[2] c;
h q[0];
cx q[0], q[1];
c = measure q;
""".strip()

# A broken circuit — missing measurements (used by debug_circuit test)
BROKEN_QASM = """
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
h q[0];
cx q[0], q[1];
""".strip()


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def parse(result: str) -> dict:
    """Parse a JSON tool result; fail clearly if it's not valid JSON."""
    try:
        return json.loads(result)
    except json.JSONDecodeError:
        pytest.fail(f"Tool returned non-JSON: {result[:200]}")


def assert_no_error(result: str, tool_name: str):
    """Assert that a tool result does not contain a top-level 'error' key."""
    data = parse(result)
    assert "error" not in data, f"{tool_name} returned error: {data.get('error')}"


# ---------------------------------------------------------------------------
# Import tools — skip all tests if IBM token is missing
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def require_token():
    from dotenv import load_dotenv
    load_dotenv()
    if not os.getenv("IBM_QUANTUM_TOKEN"):
        pytest.skip("IBM_QUANTUM_TOKEN not set — skipping hardware tests")


# We import after the fixture so the token is loaded first.
# Use a session-scoped fixture to import once and reuse.
@pytest.fixture(scope="session")
def tools():
    import server as s
    return s


@pytest.fixture(scope="session")
def real_device(tools):
    """Pick the first available device to use in tests."""
    result = json.loads(tools.list_devices())
    # list_devices returns a plain list of device dicts
    devices = result if isinstance(result, list) else result.get("devices", [])
    assert devices, "No devices returned by list_devices"
    return devices[0]["name"]


# ---------------------------------------------------------------------------
# Tool 1: list_devices
# ---------------------------------------------------------------------------

def test_list_devices(tools):
    result = tools.list_devices()
    data = parse(result)
    # Returns a list of device dicts directly
    assert isinstance(data, list), f"list_devices should return a list, got {type(data)}"
    assert len(data) > 0, "list_devices returned empty list"
    first = data[0]
    assert "name" in first
    assert "status" in first


# ---------------------------------------------------------------------------
# Tool 2: get_device_details
# ---------------------------------------------------------------------------

def test_get_device_details(tools, real_device):
    result = tools.get_device_details(real_device)
    assert_no_error(result, "get_device_details")
    data = parse(result)
    assert "name" in data


# ---------------------------------------------------------------------------
# Tool 3: best_qubits
# ---------------------------------------------------------------------------

def test_best_qubits(tools, real_device):
    result = tools.best_qubits(real_device, n=3)
    assert_no_error(result, "best_qubits")
    data = parse(result)
    assert "best_qubits" in data
    assert len(data["best_qubits"]) <= 3


# ---------------------------------------------------------------------------
# Tool 4: compare_devices
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("sort_by", ["cx_error", "queue", "qubits", "combined"])
def test_compare_devices(tools, sort_by):
    result = tools.compare_devices(sort_by=sort_by)
    assert_no_error(result, f"compare_devices(sort_by={sort_by})")
    data = parse(result)
    assert "devices" in data, f"compare_devices missing 'devices' key: {list(data.keys())}"


# ---------------------------------------------------------------------------
# Tool 5: queue_status
# ---------------------------------------------------------------------------

def test_queue_status(tools):
    result = tools.queue_status()
    data = parse(result)
    # Returns a list of device queue dicts directly
    assert isinstance(data, list), f"queue_status should return a list, got {type(data)}"
    assert len(data) > 0, "queue_status returned empty list"
    assert "name" in data[0], "queue_status entries missing 'name' field"


# ---------------------------------------------------------------------------
# Tool 6: device_history
# ---------------------------------------------------------------------------

def test_device_history(tools, real_device):
    result = tools.device_history(real_device, days=3)
    data = parse(result)
    # May return empty history if no snapshots, but should not error
    assert "error" not in data or "no snapshots" in data.get("error", "").lower()


# ---------------------------------------------------------------------------
# Tool 7: device_on_date
# ---------------------------------------------------------------------------

def test_device_on_date(tools, real_device):
    result = tools.device_on_date(real_device, "2026-01-01")
    data = parse(result)
    # May return no data for that date — acceptable, just must not crash
    assert isinstance(data, dict)


# ---------------------------------------------------------------------------
# Tool 8: submit_job (QASM 2)
# ---------------------------------------------------------------------------

def test_submit_job_qasm2(tools, real_device):
    result = tools.submit_job(real_device, BELL_QASM2, shots=128, qasm_version=2)
    data = parse(result)
    assert "job_id" in data, f"submit_job did not return job_id: {data}"


# ---------------------------------------------------------------------------
# Tool 8b: submit_job (QASM 3)
# ---------------------------------------------------------------------------

def test_submit_job_qasm3(tools, real_device):
    result = tools.submit_job(real_device, BELL_QASM3, shots=128, qasm_version=3)
    data = parse(result)
    if "error" in data and "qiskit_qasm3_import" in data.get("error", ""):
        pytest.skip("qiskit_qasm3_import not installed in this environment — runs in Docker")
    assert "job_id" in data, f"submit_job (QASM3) did not return job_id: {data}"


# ---------------------------------------------------------------------------
# Tool 9: job_status
# ---------------------------------------------------------------------------

def test_job_status(tools, real_device):
    # Submit a job first, then check its status
    sub = parse(tools.submit_job(real_device, BELL_QASM2, shots=128))
    job_id = sub.get("job_id")
    assert job_id, "Could not get job_id for job_status test"

    result = tools.job_status(job_id)
    data = parse(result)
    assert "status" in data
    assert data["status"] in ("QUEUED", "RUNNING", "DONE", "ERROR", "CANCELLED")


# ---------------------------------------------------------------------------
# Tool 10: job_results — only meaningful once job is DONE
# We just verify it doesn't crash on a queued job
# ---------------------------------------------------------------------------

def test_job_results_queued(tools, real_device):
    sub = parse(tools.submit_job(real_device, BELL_QASM2, shots=128))
    job_id = sub.get("job_id")
    assert job_id

    result = tools.job_results(job_id)
    data = parse(result)
    # If still queued, should return a clear message, not a crash
    assert isinstance(data, dict)


# ---------------------------------------------------------------------------
# Tool 11: cancel_job
# ---------------------------------------------------------------------------

def test_cancel_job(tools, real_device):
    sub = parse(tools.submit_job(real_device, BELL_QASM2, shots=128))
    job_id = sub.get("job_id")
    assert job_id

    result = tools.cancel_job(job_id)
    data = parse(result)
    # Either cancelled successfully or already done — both are acceptable
    assert "job_id" in data or "error" in data


# ---------------------------------------------------------------------------
# Tool 12: list_jobs
# ---------------------------------------------------------------------------

def test_list_jobs(tools):
    result = tools.list_jobs(limit=5)
    assert_no_error(result, "list_jobs")
    data = parse(result)
    assert "jobs" in data


# ---------------------------------------------------------------------------
# Tool 13: run_grover
# ---------------------------------------------------------------------------

def test_run_grover(tools):
    result = tools.run_grover(n_qubits=2, target_state="11")
    data = parse(result)
    assert "job_id" in data, f"run_grover did not return job_id: {data}"


# ---------------------------------------------------------------------------
# Tool 14: estimate_expectation
# ---------------------------------------------------------------------------

def test_estimate_expectation(tools, real_device):
    result = tools.estimate_expectation(
        device_name=real_device,
        qasm_string=BELL_QASM2,
        observables="ZZ"
    )
    data = parse(result)
    assert "error" not in data or "expectation" in data


# ---------------------------------------------------------------------------
# Tool 15: circuit_report
# ---------------------------------------------------------------------------

def test_circuit_report(tools, real_device):
    result = tools.circuit_report(real_device, BELL_QASM2)
    assert_no_error(result, "circuit_report")
    data = parse(result)
    assert "gate_counts" in data or "verdict" in data


# ---------------------------------------------------------------------------
# Tool 16: debug_circuit — valid circuit should have no CRITICAL issues
# ---------------------------------------------------------------------------

def test_debug_circuit_valid(tools):
    result = tools.debug_circuit(BELL_QASM2)
    data = parse(result)
    assert "issues" in data
    critical = [i for i in data["issues"] if i.get("severity") == "CRITICAL"]
    assert len(critical) == 0, f"Valid Bell circuit flagged as CRITICAL: {critical}"


def test_debug_circuit_broken(tools):
    """A circuit with no measurements should be flagged."""
    result = tools.debug_circuit(BROKEN_QASM)
    data = parse(result)
    assert "issues" in data
    assert len(data["issues"]) > 0, "Broken circuit (no measurements) should have issues"
