#!/usr/bin/env python3
"""
task43b_recover_batch_result.py -- recovers Task 43's real batch-job
results directly via the API, bypassing qiskit-ionq's buggy multi-circuit
`.result()` parsing (crashed on a child job id it mistook for count data).
Read-only against an ALREADY-COMPLETED real job -- zero additional cost.
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from ionq_backend import connect_provider

PARENT_JOB_ID = "01a03b16-8a92-71fa-bd15-03419a95074b"
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task43_batch_hw_submission_results.json")


def main():
    provider = connect_provider()
    backend = provider.get_backend("qpu.forte-enterprise-1", gateset="native")
    client = backend.client

    parent_raw = client.retrieve_job(PARENT_JOB_ID)
    print("\n-- PARENT JOB RAW METADATA --")
    print(json.dumps(parent_raw, indent=2, default=str))

    child_ids = parent_raw.get("child_job_ids")
    print(f"\nchild_job_ids: {child_ids}")

    children_raw = []
    if child_ids:
        for cid in child_ids:
            craw = client.retrieve_job(cid)
            children_raw.append(craw)
            print(f"\n-- CHILD JOB {cid} RAW METADATA --")
            print(json.dumps(craw, indent=2, default=str))

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": PARENT_JOB_ID, "parent_raw_job_metadata": parent_raw,
                    "child_job_ids": child_ids, "children_raw_job_metadata": children_raw},
                   f, indent=2, default=str)
    print(f"\nSaved -> {RESULTS_PATH}")


if __name__ == "__main__":
    main()
