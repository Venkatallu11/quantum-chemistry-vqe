#!/usr/bin/env python3
"""
task49_gc_hard_group_hw_submission.py -- iteration 49. FIFTH REAL
HARDWARE SUBMISSION, explicit user go-ahead given (2026-09-09). First
real-hardware test of the general-commuting (GC) measurement scheme:
the "hard" GC group (index 3, most diagonalizer gates: 3 extra 1q / 5
extra 2q abstract) on 3 real kept slots (u_0, u_1, u_2), batched as ONE
job, matching Task 43's own proven multicircuit-batch pattern exactly.

PRE-FLIGHT, already done and passed (task48_gc_native_hw_prep.py):
  - native-gate translation of the GC diagonalizer verified against
    exact ideal reconstruction, worst error 9.3e-15 across all 3 slots
  - ancilla leakage-free certainty verified, worst p(ancilla=1)=6.9e-27
  - real per-slot native gate counts discovered (they vary by slot --
    the joint TrappedIonOptimizerPlugin pass finds different fusion with
    the diagonalizer depending on the slot's own state-prep angles):
    u_0: 171x1q/16x2q, u_1: 180x1q/13x2q, u_2: 195x1q/20x2q
  - estimated real cost (validated cost model, cross-checked against
    Vadim's own real breakeven numbers): ~$288.85 for this 3-circuit
    batch at 2000 shots

SHOTS: 2000, matching Task 43's own precedent exactly for direct
comparability of retention fractions and expectation convergence against
the existing real baseline data.

Run:
    PYTHONHASHSEED=0 python vqe/task49_gc_hard_group_hw_submission.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from general_commuting_measurements import build_general_commuting_measurement_plan
from native_stateprep import to_native
from ionq_backend import connect_provider

K = 6
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOTS = 2000
TARGET_SLOTS = ["u_0", "u_1", "u_2"]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task49_gc_hard_group_hw_submission_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task49_gc_hard_group_hw_submission.py -- FIFTH REAL HARDWARE SUBMISSION (GC hard-group spot-check)")
    print("=" * 96)
    print(f"  backend={BACKEND_NAME}  shots={SHOTS}  -- 3 circuits (GC group 3 on u_0/u_1/u_2), ONE batch job")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    assert n_ok == 36

    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    hard_idx = max(range(len(diagonalizers)), key=lambda i: diagonalizers[i].two_qubit_count)
    group, d = groups[hard_idx], diagonalizers[hard_idx]
    native_diag = to_native(d.to_circuit(), GATE_NAME)
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    print(f"  GC group {hard_idx} (the 'hard' one): {len(group)} labels: {group}")

    circuits = []
    gate_counts_per_slot = {}
    for name in TARGET_SLOTS:
        angles = fixed_solutions[name]["angles"]
        register = optimized_native_circuit(angles, GATE_NAME)
        full5 = with_ancilla_parity_native(register, ancilla_native)
        full5 = full5.compose(native_diag, qubits=[0, 1, 2, 3])
        gate_counts_per_slot[name] = dict(full5.count_ops())
        full5.measure_all()
        circuits.append(full5)
        print(f"    slot={name}: gate counts {gate_counts_per_slot[name]}")

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")
    print(f"\n  connected to {backend.name}")

    print(f"\n  SUBMITTING NOW -- real hardware, real cost will be incurred (3 circuits, one batch job).")
    job = backend.run(circuits, shots=SHOTS)
    job_id = job.job_id()
    print(f"  submitted. job_id={job_id}")

    results = {"job_id": job_id, "backend": BACKEND_NAME, "shots": SHOTS,
               "hard_group_index": hard_idx, "group": group, "target_slots": TARGET_SLOTS,
               "gate_counts_per_slot": gate_counts_per_slot, "status": "submitted"}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
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
        for name, c in zip(TARGET_SLOTS, counts_list):
            total = sum(c.values())
            retained = sum(v for k, v in c.items() if k[0] == "0")
            print(f"  slot={name}: total_shots={total}  ancilla-retained={retained} ({retained/total*100:.1f}%)")

    raw = backend.client.retrieve_job(job_id)
    results["status"] = status
    results["counts_list"] = counts_list
    results["raw_job_metadata"] = raw
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Full results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
