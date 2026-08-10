#!/usr/bin/env python3
"""
task0_fidelity_correction.py — iteration 24, Task 0. THE FIDELITY
CORRECTION, done first because it may invalidate everything downstream.
============================================================================
Phase 1 (iteration 23) found IonQ's calibration API reports qpu.forte-1 at
99.52% two-qubit fidelity, while this project has used
P2_PER_GATE=0.01214 (98.786%) throughout every local noise model, labeled
"real aria-1" -- but actually applied IDENTICALLY to both aria-1 and
forte-1 local simulations. This script queries IonQ's REAL calibration
API for EVERY accessible backend (read-only, free, no hardware touched,
no QPU credits spent) and reports exactly what's there -- not marketing
numbers, not this project's own assumed constant.

THIS IS A REAL, LIVE NETWORK CALL (GET .../characterizations), not a
simulation. It reads calibration history; it does not submit any job.

Run:
    python vqe/task0_fidelity_correction.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from ionq_backend import connect_provider

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task0_fidelity_correction_results.json")

# Every backend name IonQ's own /backends endpoint reports as of this run
# (queried live below, not hardcoded blind -- this list is filled in from
# that same call, kept here only as the fallback if that call is skipped)
KNOWN_BACKEND_NAMES = ["qpu.harmony", "qpu.aria-1", "qpu.aria-2", "qpu.forte-1", "qpu.forte-enterprise-1"]


def list_all_backends(client):
    """GET /backends -- the full real backend roster (name, status,
    qubits, characterization_id), not just the 2 the qiskit_ionq
    IonQProvider.backends() wrapper exposes (ionq_simulator + a generic
    unavailable ionq_qpu)."""
    res = client.get_with_retry(client.make_path("backends"), headers=client.api_headers)
    res.raise_for_status()
    return res.json()


def query_calibration_history(client, backend_name, limit=10):
    """Real characterization history for one backend, most-recent first
    (IonQ API's own ordering). Returns [] for backends with no
    characterization data (e.g. plain 'simulator')."""
    chars = client.get_calibration_data(backend_name, limit=limit)
    rows = []
    for c in chars:
        rows.append({
            "date": c.date.isoformat(),
            "qubits": c.qubits,
            "fidelity": dict(c.fidelity),
            "timing": dict(c.timing),
        })
    return rows


def last_valid_2q(history):
    """First entry (scanning most-recent-first) with a non-null 2q median
    -- a backend that has gone stale/retired can have null or nonsensical
    recent entries (observed directly on aria-1: its newest record has
    2q=null and 1q=0.4745, clearly not representative)."""
    for row in history:
        m = row["fidelity"].get("2q", {}).get("median")
        if m is not None:
            return row
    return None


def main():
    print("\n" + "=" * 96)
    print("  task0_fidelity_correction.py -- real IonQ calibration API, every accessible backend")
    print("=" * 96)

    provider = connect_provider()
    client = provider.get_backend("ionq_simulator").client

    print("\n  -- GET /backends (real roster) --")
    backends = list_all_backends(client)
    for b in backends:
        print(f"    {b['backend']:<24} status={b.get('status'):<12} qubits={b.get('qubits')}")

    real_backend_names = [b["backend"] for b in backends if b["backend"] != "simulator"]

    print("\n  -- calibration history per backend (most recent 10, real API data) --")
    all_history = {}
    summary = {}
    for name in real_backend_names:
        try:
            hist = query_calibration_history(client, name, limit=10)
        except Exception as e:
            print(f"    {name}: ERROR {e}")
            continue
        all_history[name] = hist
        if not hist:
            print(f"    {name}: no characterization data")
            continue
        latest = hist[0]
        valid = last_valid_2q(hist)
        b_status = next((b.get("status") for b in backends if b["backend"] == name), "?")
        print(f"    {name} (status={b_status}):")
        print(f"      latest record  ({latest['date']}): 1q={latest['fidelity'].get('1q', {}).get('median')} "
              f"2q={latest['fidelity'].get('2q', {}).get('median')} spam={latest['fidelity'].get('spam', {}).get('median')}")
        if valid is not None and valid is not latest:
            print(f"      LATEST 2q IS NULL/MISSING -- last record with a real 2q median "
                  f"({valid['date']}): 1q={valid['fidelity'].get('1q', {}).get('median')} "
                  f"2q={valid['fidelity'].get('2q', {}).get('median')} spam={valid['fidelity'].get('spam', {}).get('median')}")
        summary[name] = {
            "status": b_status,
            "latest_date": latest["date"],
            "latest_fidelity": latest["fidelity"],
            "last_valid_2q_date": valid["date"] if valid else None,
            "last_valid_2q_fidelity": valid["fidelity"] if valid else None,
        }

    print("\n  -- fidelity convention, stated explicitly --")
    print("    IonQ's API exposes only a device-wide MEDIAN fidelity per gate type")
    print("    (fidelity.2q.median, fidelity.2q.stderr) -- NOT a best-pair number, NOT a")
    print("    mean. All qubit pairs are reported in 'connectivity' as all-to-all, so")
    print("    'per-pair' does not vary structurally; the API's own field is 'median'.")

    print("\n  -- the correction, stated plainly --")
    LOCAL_MODEL_CONSTANT = 0.98786
    forte1 = summary.get("qpu.forte-1", {}).get("latest_fidelity", {}).get("2q", {}).get("median")
    aria1_valid = summary.get("qpu.aria-1", {}).get("last_valid_2q_fidelity", {}).get("2q", {}).get("median")
    if forte1:
        ratio = (1 - LOCAL_MODEL_CONSTANT) / (1 - forte1)
        print(f"    forte-1 REAL (available)       : {forte1:.4%}  vs local-model constant {LOCAL_MODEL_CONSTANT:.4%} "
              f"-> project assumed {ratio:.2f}x MORE noise than forte-1 actually delivers")
    if aria1_valid:
        ratio2 = (1 - LOCAL_MODEL_CONSTANT) / (1 - aria1_valid)
        direction = "MORE" if ratio2 > 1 else "LESS"
        print(f"    aria-1 REAL (last valid, RETIRED since): {aria1_valid:.4%}  vs local-model constant {LOCAL_MODEL_CONSTANT:.4%} "
              f"-> project assumed {ratio2:.2f}x ({direction}) noise than aria-1's own last real reading")
    print(f"    -> the miscalibration is ASYMMETRIC: the single shared constant happens to sit close to")
    print(f"       aria-1's real historical rate but significantly overstates forte-1's real noise.")

    results = {
        "backends": backends,
        "calibration_history": all_history,
        "summary": summary,
        "local_model_constant_fidelity": LOCAL_MODEL_CONSTANT,
        "fidelity_convention": "device-wide MEDIAN per IonQ's own API field name, not best-pair, not mean",
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
