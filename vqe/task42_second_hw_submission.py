#!/usr/bin/env python3
"""
task42_second_hw_submission.py -- iteration 42. SECOND REAL HARDWARE
SUBMISSION, explicit user go-ahead given. Purpose: cost-scaling probe.

Task 41's first real submission (100 shots, same circuit) cost $25.79
per the user's own IonQ dashboard -- the API itself exposes no dollar
field. It is unknown whether that cost was dominated by fixed per-circuit
overhead (gate/qubit count) or scaled with shots. User explicitly chose
a controlled 5x shot increase (500 shots, vs Task 41's 100) as the next
step, over a riskier 10x/max-shots jump, specifically to bound worst-case
risk against the $170 total budget while still getting a real scaling
signal: SAME circuit, SAME backend, ONLY shots changed.

If cost stays close to $25.79 -> overhead-dominated, larger future
circuits/shot counts are relatively safe. If cost scales roughly
linearly (~5x, ~$129) -> shots are expensive and any full-protocol-style
real hardware replication remains out of reach for this budget.

SCOPE: ONE real H4 circuit -- the SAME ancilla-augmented, native-gate,
u_0/(XYYX,IYYI) measurement circuit as Task 41, unchanged -- submitted to
qpu.forte-enterprise-1 at SHOTS=500.

Run:
    PYTHONHASHSEED=0 python vqe/task42_second_hw_submission.py
"""
import os
import sys
import json

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
SHOTS = 500  # 5x Task 41's 100 shots -- user-chosen, bounded cost-scaling probe
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task42_second_hw_submission_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task42_second_hw_submission.py -- SECOND REAL HARDWARE SUBMISSION (cost scaling probe)")
    print("=" * 96)
    print(f"  backend={BACKEND_NAME}  shots={SHOTS}  -- ONE circuit only, SAME circuit as Task 41")

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

    print(f"\n  -- waiting for completion --")
    job.wait_for_final_state(timeout=14400)  # 4 hour ceiling
    status = str(job.status())
    print(f"  final status: {status}")

    if job.done():
        counts = job.get_counts()
        print(f"  REAL counts: {counts}")

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
