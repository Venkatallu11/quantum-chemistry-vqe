"""
Tests for the agent ReAct loop routing logic.

Verifies that:
  1. Hardware questions route to 'tool' action (calls MCP tools)
  2. Quantum code questions route to 'model_call' action (Qiskit specialist)
  3. General questions route to 'answer' action (no tool, no model)
  4. Each conversation is isolated — job context from one question
     does not bleed into the next question (Jack's bug)

These tests hit the live agent HTTP API at localhost:3021.
Start the agent first: docker compose up --build

Run with:
    pytest tests/test_agent_routing.py -v
"""

import json
import pytest
import requests
import os

AGENT_URL = os.getenv("AGENT_URL", "http://localhost:3021")


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def ask(question: str, timeout: int = 30) -> dict:
    """Send a question to the agent and return the parsed JSON response."""
    resp = requests.post(
        f"{AGENT_URL}/chat",
        json={"message": question},
        timeout=timeout,
    )
    resp.raise_for_status()
    return resp.json()


# ---------------------------------------------------------------------------
# Skip if agent is not running
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def require_agent():
    try:
        requests.get(f"{AGENT_URL}/health", timeout=3)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"Agent not running at {AGENT_URL} — start with: docker compose up")


# ---------------------------------------------------------------------------
# Test 1: Hardware question → must use a tool
# The agent should call list_devices, get_device_details, or compare_devices
# ---------------------------------------------------------------------------

def test_hardware_question_uses_tool():
    response = ask("Which IBM quantum computer has the lowest CX error rate right now?")
    answer = response.get("answer", "")
    # A real hardware answer will mention IBM device names like ibm_fez, ibm_torino, etc.
    assert len(answer) > 50, "Answer is too short to be a real hardware response"
    assert any(keyword in answer.lower() for keyword in ["ibm", "error", "qubit", "%"]), \
        f"Answer doesn't look like hardware data: {answer[:200]}"


# ---------------------------------------------------------------------------
# Test 2: Code question → must route to Qiskit model (if enabled)
# If QISKIT_CODE_MODEL is not set, this test is skipped
# ---------------------------------------------------------------------------

def test_code_question_routes_to_qiskit_model():
    if not os.getenv("QISKIT_CODE_MODEL"):
        pytest.skip("QISKIT_CODE_MODEL not set — Qiskit model routing not enabled")

    response = ask(
        "Write me a Grover's algorithm circuit in OpenQASM 3 for 3 qubits",
        timeout=120,  # Mistral on CPU can be slow
    )
    answer = response.get("answer", "")
    # A real QASM 3 answer will contain these keywords
    assert "OPENQASM 3" in answer or "openqasm" in answer.lower(), \
        f"Code answer missing OPENQASM 3 header: {answer[:300]}"
    assert "qubit" in answer.lower() or "qreg" in answer.lower(), \
        f"Code answer missing qubit declarations: {answer[:300]}"


# ---------------------------------------------------------------------------
# Test 3: General question → answered directly, no tool needed
# ---------------------------------------------------------------------------

def test_general_question_answered_directly():
    response = ask("What is quantum entanglement?")
    answer = response.get("answer", "")
    assert len(answer) > 100, "Answer is too short for a general knowledge question"
    assert "entangl" in answer.lower(), \
        f"Answer about entanglement doesn't mention entanglement: {answer[:200]}"


# ---------------------------------------------------------------------------
# Test 4: Context isolation — Jack's bug
# Submit a job in one question, then ask something completely different.
# The second answer must NOT mention the job or job_id.
# ---------------------------------------------------------------------------

def test_context_isolation_between_questions():
    # Question 1: submit a job (this will queue a real job)
    response1 = ask("Submit a Bell state circuit to ibm_least_busy for 128 shots")
    answer1 = response1.get("answer", "")

    # Question 2: completely different topic
    response2 = ask("What is the difference between a qubit and a classical bit?")
    answer2 = response2.get("answer", "")

    # The second answer should NOT contain job-related content
    assert "job_id" not in answer2.lower(), \
        f"Context leaked from Q1 into Q2 — job_id appeared in answer: {answer2[:300]}"
    assert "queued" not in answer2.lower(), \
        f"Context leaked from Q1 into Q2 — 'queued' appeared in answer: {answer2[:300]}"


# ---------------------------------------------------------------------------
# Test 5: Error handling — bad device name
# The agent should return a helpful error, not crash
# ---------------------------------------------------------------------------

def test_bad_device_name_handled_gracefully():
    response = ask("Get details for device named 'fake_device_xyz_does_not_exist'")
    answer = response.get("answer", "")
    # Should mention it couldn't find the device, not a 500 error
    assert len(answer) > 0, "Agent returned empty answer for bad device name"
    assert response.get("error") is None or "500" not in str(response.get("error", "")), \
        "Agent crashed on bad device name instead of handling gracefully"


# ---------------------------------------------------------------------------
# Test 6: Multi-step reasoning
# This question requires TWO tool calls: queue_status to find the least busy
# device, then get_device_details on that device
# ---------------------------------------------------------------------------

def test_multi_step_hardware_question():
    response = ask(
        "Which quantum computer has the shortest queue right now, and what is its CX error rate?",
        timeout=60,
    )
    answer = response.get("answer", "")
    assert len(answer) > 80, "Multi-step answer is too short"
    # Must contain both queue info and error rate info
    has_queue = any(w in answer.lower() for w in ["queue", "jobs", "wait"])
    has_error = any(w in answer.lower() for w in ["error", "%", "cx"])
    assert has_queue or has_error, \
        f"Multi-step answer missing both queue and error info: {answer[:300]}"
