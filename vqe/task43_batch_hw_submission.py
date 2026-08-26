#!/usr/bin/env python3
"""
task43_batch_hw_submission.py -- iteration 43. THIRD REAL HARDWARE
SUBMISSION, explicit user go-ahead given. Purpose: circuit-count scaling
probe (complements Task 41/42's shot-count probe).

Task 41 (100 shots) and Task 42 (500 shots) of the SAME single circuit
both cost ~$25.xx on qpu.forte-enterprise-1 -- real, dashboard-confirmed
evidence that cost here is overhead-DOMINATED per circuit, not per shot.
The open question this script answers: does submitting MULTIPLE circuits
together in ONE job scale ~linearly in circuit count (e.g. 3 circuits
~3x$25), or is there some shared/batch efficiency?

SCOPE: 3 real circuits (u_0 slot, 3 different measurement groups
including the SAME [XYYX,IYYI] group Task 41/42 already used, so one of
the three is directly comparable) submitted TOGETHER as one batch job
(backend.run(list_of_circuits, shots=...)) to qpu.forte-enterprise-1,
2000 shots each (shots shown cheap by Task 42, so raised for better
precision this time).

Run:
    PYTHONHASHSEED=0 python vqe/task43_batch_hw_submission.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from task39c_ancilla_real_submission import load_groups_from_baseline
from ionq_backend import connect_provider

K = 6
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOTS = 2000  # raised from Task 42's 500 -- shots shown cheap, circuit count is the variable under test
SLOT = "u_0"
N_CIRCUITS = 3
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task43_batch_hw_submission_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task43_batch_hw_submission.py -- THIRD REAL HARDWARE SUBMISSION (circuit-count scaling probe)")
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
    print(f"  gate counts per circuit: {gate_counts}")

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")
    print(f"  connected to {backend.name}")

    print(f"\n  SUBMITTING NOW -- real hardware, real cost will be incurred ({N_CIRCUITS} circuits, one batch job).")
    job = backend.run(circuits, shots=SHOTS)
    job_id = job.job_id()
    print(f"  submitted. job_id={job_id}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS, "groups": groups,
                    "gate_counts": gate_counts, "status": "submitted"}, f, indent=2)
    print(f"  job_id saved -> {RESULTS_PATH} (safe to check on this job later via its ID even if this process exits)")

    print(f"\n  -- waiting for completion --")
    job.wait_for_final_state(timeout=14400)
    status = str(job.status())
    print(f"  final status: {status}")

    counts_list = None
    if job.done():
        result = job.result()
        counts_list = result.get_counts()
        if not isinstance(counts_list, list):
            counts_list = [counts_list]
        counts_list = [dict(c) for c in counts_list]
        for g, c in zip(groups, counts_list):
            print(f"  REAL counts for group {g}: {c}")

    raw = backend.client.retrieve_job(job_id)
    print(f"\n  -- FULL RAW JOB METADATA (looking for any cost/price field) --")
    print(json.dumps(raw, indent=2, default=str))

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS, "groups": groups,
                    "gate_counts": gate_counts, "status": status,
                    "counts_list": counts_list, "raw_job_metadata": raw}, f, indent=2, default=str)
    print(f"\n  Full results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
