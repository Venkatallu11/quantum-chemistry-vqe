#!/usr/bin/env python3
"""
task43c_batch_hw_submission_v2.py -- retry of Task 43's 3-circuit batch
submission after the user raised the real IonQ project budget to $300
(Task 43's first attempt failed with QuotaExhaustedError on all 3 child
jobs before any execution -- $0 spent that time).

Also fixes the real bug found in Task 43's first attempt: qiskit-ionq
1.1.1's `job.result()` / `get_counts()` crashes on multi-circuit QPU jobs
(`ionq.multi-circuit.v1`) because the parent job's own result payload
contains child_job_ids, not per-circuit histogram data, and the SDK does
not handle that split for QPU backends. Fix: bypass `job.result()`
entirely -- poll the parent job's raw status directly, then fetch each
CHILD job's own histogram via `IonQClient.get_results()` and convert it
with the SDK's own (correct, well-tested for single jobs) `_build_counts`
helper, applied per child. This is a results-retrieval fix only; it does
not change what was submitted (same 3 circuits, same 2000 shots each).

Run:
    PYTHONHASHSEED=0 python vqe/task43c_batch_hw_submission_v2.py
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from task39c_ancilla_real_submission import load_groups_from_baseline
from ionq_backend import connect_provider
from qiskit_ionq.ionq_job import _build_counts

K = 6
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOTS = 2000
SLOT = "u_0"
N_CIRCUITS = 3
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task43_batch_hw_submission_results.json")


def poll_until_terminal(client, job_id, timeout=14400, interval=5):
    terminal = {"completed", "failed", "canceled"}
    start = time.time()
    while True:
        raw = client.retrieve_job(job_id)
        status = raw.get("status")
        if status in terminal:
            return raw
        if time.time() - start > timeout:
            raise TimeoutError(f"job {job_id} did not reach a terminal state within {timeout}s (last status={status})")
        time.sleep(interval)


def fetch_child_counts(client, child_raw, num_qubits, num_clbits, shots):
    if child_raw.get("status") != "completed":
        return None, child_raw.get("failure")
    results = child_raw.get("results") or {}
    hist_entry = results.get("histogram")
    if not hist_entry or "url" not in hist_entry:
        return None, "no histogram results url present"
    data = client.get_results(hist_entry["url"])
    counts, probs = _build_counts(data, num_qubits, list(range(num_clbits)), shots)
    return counts, None


def main():
    print("\n" + "=" * 96)
    print("  task43c_batch_hw_submission_v2.py -- RETRY, budget raised to $300, results-fetch bug fixed")
    print("=" * 96)
    print(f"  backend={BACKEND_NAME}  shots={SHOTS}  -- {N_CIRCUITS} circuits, ONE batch job")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups_by_slot = load_groups_from_baseline(kept)
    all_groups = groups_by_slot[SLOT]

    preferred_first = ["XYYX", "IYYI"]
    ordered_groups = sorted(all_groups, key=lambda g: 0 if list(g) == preferred_first else 1)
    groups = ordered_groups[:N_CIRCUITS]
    print(f"  slot={SLOT}, groups selected: {groups}")

    register = optimized_native_circuit(fixed_solutions[SLOT]["angles"], GATE_NAME)
    ancilla_native = ancilla_cnots_native(GATE_NAME)

    circuits = []
    for group in groups:
        combined = effrag_mod.combined_basis_label(group)
        basis_qc = native_basis_change(combined, GATE_NAME)
        full5 = with_ancilla_parity_native(register, ancilla_native)
        full5 = full5.compose(basis_qc, qubits=[0, 1, 2, 3])
        full5.measure_all()
        circuits.append(full5)
    gate_counts = [dict(c.count_ops()) for c in circuits]
    num_qubits = circuits[0].num_qubits
    num_clbits = circuits[0].num_clbits
    print(f"  gate counts per circuit: {gate_counts}  (num_qubits={num_qubits}, num_clbits={num_clbits})")

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")
    client = backend.client
    print(f"  connected to {backend.name}")

    print(f"\n  SUBMITTING NOW -- real hardware, real cost will be incurred ({N_CIRCUITS} circuits, one batch job).")
    job = backend.run(circuits, shots=SHOTS)
    job_id = job.job_id()
    print(f"  submitted. job_id={job_id}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS, "groups": groups,
                    "gate_counts": gate_counts, "status": "submitted"}, f, indent=2)
    print(f"  job_id saved -> {RESULTS_PATH}")

    print(f"\n  -- polling parent job status directly (avoiding the buggy SDK .result() path) --")
    parent_raw = poll_until_terminal(client, job_id)
    print(f"  parent final status: {parent_raw.get('status')}")
    print(json.dumps(parent_raw, indent=2, default=str))

    child_ids = parent_raw.get("child_job_ids") or []
    children_raw = []
    counts_list = []
    for cid in child_ids:
        craw = poll_until_terminal(client, cid)
        children_raw.append(craw)
        counts, err = fetch_child_counts(client, craw, num_qubits, num_clbits, SHOTS)
        counts_list.append(counts)
        if counts is not None:
            print(f"\n  REAL counts for child {cid}: {counts}")
        else:
            print(f"\n  child {cid} FAILED or has no data: {err}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS, "groups": groups,
                    "gate_counts": gate_counts, "status": parent_raw.get("status"),
                    "counts_list": counts_list, "parent_raw_job_metadata": parent_raw,
                    "child_job_ids": child_ids, "children_raw_job_metadata": children_raw},
                   f, indent=2, default=str)
    print(f"\n  Full results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
