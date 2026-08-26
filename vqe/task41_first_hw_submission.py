#!/usr/bin/env python3
"""
task41_first_hw_submission.py -- iteration 41. FIRST REAL HARDWARE
SUBMISSION, explicit user go-ahead given, budget $170 (raised from an
initial $35 specifically because real per-shot cost could not be
verified in advance via the API's own cost-estimate endpoint, which
returned cost=None for this account on every scenario tested).

SCOPE, deliberately minimal for this FIRST real-hardware data point:
ONE real H4 circuit -- the SAME ancilla-augmented, native-gate,
u_0/(XYYX,IYYI) measurement circuit this project's own certified
pipeline uses (Task 39B/C's own construction, unchanged) -- submitted to
qpu.forte-enterprise-1 (the only backend confirmed accessible on this
account via the dashboard; qpu.forte-1's own API `.status()` reports
True but is NOT confirmed actually submittable on this account) at a
SMALL, conservative shot count (100) specifically to discover the REAL
billed cost for this account before committing to anything larger.

After completion, retrieves the FULL raw job metadata (bypassing the
SDK's thin wrapper, which exposes no cost field) via the underlying
IonQClient.retrieve_job call directly, and prints every field found --
if IonQ's API reports actual cost anywhere in that payload, it will show
up here.

Run:
    PYTHONHASHSEED=0 python vqe/task41_first_hw_submission.py
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from ionq_backend import connect_provider

K = 6
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOTS = 100  # deliberately small for this FIRST real-hardware, cost-discovery submission
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task41_first_hw_submission_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task41_first_hw_submission.py -- FIRST REAL HARDWARE SUBMISSION (cost discovery)")
    print("=" * 96)
    print(f"  backend={BACKEND_NAME}  shots={SHOTS}  -- ONE circuit only")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)

    register = optimized_native_circuit(fixed_solutions["u_0"]["angles"], GATE_NAME)
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    full5 = with_ancilla_parity_native(register, ancilla_native)
    group = ["XYYX", "IYYI"]
    combined = effrag_mod.combined_basis_label(group)
    basis_qc = native_basis_change(combined, GATE_NAME)
    circuit = full5.compose(basis_qc, qubits=[0, 1, 2, 3])
    circuit.measure_all()
    counts_ops = circuit.count_ops()
    print(f"  circuit: u_0 slot, group={group}, gate counts={dict(counts_ops)}, qubits={circuit.num_qubits}")

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")
    print(f"  connected to {backend.name}")

    print(f"\n  SUBMITTING NOW -- real hardware, real cost will be incurred.")
    job = backend.run(circuit, shots=SHOTS)
    job_id = job.job_id()
    print(f"  submitted. job_id={job_id}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS,
                    "gate_counts": dict(counts_ops), "status": "submitted"}, f, indent=2)
    print(f"  job_id saved -> {RESULTS_PATH} (safe to check on this job later via its ID even if this process exits)")

    print(f"\n  -- waiting for completion (expect ~2 hour queue based on the dashboard's own live estimate) --")
    job.wait_for_final_state(timeout=14400)  # 4 hour ceiling, generous given the ~2hr queue estimate
    status = str(job.status())
    print(f"  final status: {status}")

    if job.done():
        counts = job.get_counts()
        print(f"  REAL counts: {counts}")

    # bypass the thin SDK wrapper -- get the FULL raw job payload directly, looking for any cost field
    raw = backend.client.retrieve_job(job_id)
    print(f"\n  -- FULL RAW JOB METADATA (looking for any cost/price field) --")
    print(json.dumps(raw, indent=2, default=str))

    with open(RESULTS_PATH, "w") as f:
        json.dump({"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS,
                    "gate_counts": dict(counts_ops), "status": status,
                    "counts": counts if job.done() else None, "raw_job_metadata": raw}, f, indent=2, default=str)
    print(f"\n  Full results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
