#!/usr/bin/env python3
"""
task44_1q_gate_cost_split.py -- iteration 44. FOURTH REAL HARDWARE
SUBMISSION, explicit user go-ahead given (2026-09-09). Tests Vadim
Karpusenko's (IonQ Research) Sep 8 hypothesis: our real per-shot billing
is dominated by 1-qubit gate count (120 gpi/gpi2), not 2-qubit gate count
(11 zz), because IonQ's real gate-shot cost model counts every native gate,
not just the entangling ones.

DESIGN: two single-circuit real jobs on qpu.forte-enterprise-1, 2000
shots each (same shot count as Task 43's batch, where per-circuit gate
cost -- not the flat per-job minimum charge -- actually drives price).

  Circuit A (baseline): EXACTLY Task 43's own u_0/[XYYX,IYYI] circuit --
  register optimized via TrappedIonOptimizerPlugin(opt_level=3) ALONE,
  then composed with the ancilla-CNOT sub-circuit and the basis-change
  sub-circuit, each of which was ONLY ever transpiled in isolation via
  to_native()'s plain optimization_level=1 (see native_stateprep.to_native,
  task2_fold_response_dataset.native_basis_change, and
  task39b_native_ancilla_parity's own docstring, which deliberately never
  lets the transpiler see the composed circuit, for compiler-safety
  reasons unrelated to gate-count minimization).

  Circuit B (test): the IDENTICAL composed 5-qubit circuit, but with
  TrappedIonOptimizerPlugin(opt_level=3) -- the SAME official IonQ plugin
  already trusted elsewhere in this codebase (task28d_all_gate_zne.py's
  optimized_native_circuit) -- run ONCE across the FULL composed circuit
  (register + ancilla + basis-change together), instead of stitching
  three separately-optimized-in-isolation pieces. This is the genuine,
  previously-unexploited lever Vadim pointed at: nothing in this
  construction ever let IonQ's own strongest optimizer fuse/cancel 1-qubit
  gates ACROSS the seams between register, ancilla, and basis-change.

VERIFICATION (refuse to submit if this fails -- matches this project's
own established discipline, e.g. task28d's <1e-12 fold-equivalence gate
and task39b's marginal-state/ancilla-certainty checks): Circuit B must be
statevector-identical to Circuit A (up to global phase) to <1e-12 BEFORE
either circuit is queued for submission. This is the only thing that
makes "same physics, cheaper gates" a fair, honest A/B cost comparison
rather than two different circuits.

COST: real hardware, real money. At 2000 shots and Task 43's own real
$64.02/circuit rate, ESTIMATED ~$64 x 2 = ~$128 total (no explicit cost
field appears anywhere in the real API response for either job -- as
established in Tasks 41-43, actual dollar cost can only be confirmed via
the user's own IonQ Cloud dashboard after billing settles).

Run:
    PYTHONHASHSEED=0 python vqe/task44_1q_gate_cost_split.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
from native_stateprep import native_target
import ef_fragment as effrag_mod
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from task39c_ancilla_real_submission import load_groups_from_baseline
from ionq_backend import connect_provider
from qiskit.quantum_info import Statevector

K = 6
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOTS = 2000
SLOT = "u_0"
GROUP = ["XYYX", "IYYI"]  # identical to Task 41/42/43's own first circuit
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task44_1q_gate_cost_split_results.json")


def build_circuit_a(fixed_solutions):
    """EXACTLY Task 43's construction: three separately-optimized pieces
    composed, never jointly re-optimized."""
    register = optimized_native_circuit(fixed_solutions[SLOT]["angles"], GATE_NAME)
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    combined = effrag_mod.combined_basis_label(GROUP)
    basis_qc = native_basis_change(combined, GATE_NAME)
    full5 = with_ancilla_parity_native(register, ancilla_native)
    full5 = full5.compose(basis_qc, qubits=[0, 1, 2, 3])
    return full5


def build_circuit_b(circuit_a):
    """Same composed 5-qubit circuit, but jointly re-optimized by IonQ's
    own TrappedIonOptimizerPlugin at opt_level=3 -- the same plugin/level
    task28d already trusts for the register alone, now applied across the
    full seam-including circuit."""
    from qiskit.transpiler import PassManagerConfig
    from qiskit_ionq import TrappedIonOptimizerPlugin
    tgt5 = native_target(5, GATE_NAME)
    pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt5), optimization_level=3)
    return pm.run(circuit_a)


def verify_equivalent(qc_a, qc_b, tol=1e-12):
    sv_a = np.asarray(Statevector.from_instruction(qc_a))
    sv_b = np.asarray(Statevector.from_instruction(qc_b))
    idx = int(np.argmax(np.abs(sv_a)))
    phase = sv_b[idx] / sv_a[idx] if abs(sv_a[idx]) > 1e-9 else 1.0
    err = float(np.max(np.abs(sv_b / phase - sv_a)))
    return err


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                         help="build + verify + print gate counts only, no real submission / no money spent")
    args = parser.parse_args()

    print("\n" + "=" * 96)
    print("  task44_1q_gate_cost_split.py -- FOURTH REAL HARDWARE SUBMISSION (1q-gate cost-split test)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)

    qc_a = build_circuit_a(fixed_solutions)
    qc_b = build_circuit_b(qc_a)

    err = verify_equivalent(qc_a, qc_b)
    counts_a = dict(qc_a.count_ops())
    counts_b = dict(qc_b.count_ops())
    print(f"  Circuit A gate counts (Task 43 baseline): {counts_a}")
    print(f"  Circuit B gate counts (jointly re-optimized): {counts_b}")
    print(f"  statevector-equivalence check (must be <1e-12): {err:.2e}  "
          f"{'PASS' if err < 1e-12 else 'FAIL -- REFUSING TO SUBMIT'}")

    if err >= 1e-12:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"status": "ABORTED_verification_failed", "verify_err": err,
                        "counts_a": counts_a, "counts_b": counts_b}, f, indent=2)
        print(f"\n  ABORTED: circuits are not equivalent -- no submission made. See {RESULTS_PATH}\n")
        return

    n1q_a = counts_a.get("gpi", 0) + counts_a.get("gpi2", 0)
    n1q_b = counts_b.get("gpi", 0) + counts_b.get("gpi2", 0)
    n2q_a = counts_a.get(GATE_NAME, 0)
    n2q_b = counts_b.get(GATE_NAME, 0)
    print(f"\n  1q gates: A={n1q_a}  B={n1q_b}  (reduction={n1q_a - n1q_b})")
    print(f"  2q gates: A={n2q_a}  B={n2q_b}  (must match for a fair test: {'OK' if n2q_a == n2q_b else 'MISMATCH -- interpret with caution'})")

    if n1q_b >= n1q_a:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"status": "ABORTED_no_reduction_found", "verify_err": err,
                        "counts_a": counts_a, "counts_b": counts_b}, f, indent=2)
        print(f"\n  ABORTED: joint re-optimization found no 1q-gate reduction (B >= A) -- "
              f"nothing to test, not submitting real money for an identical circuit. See {RESULTS_PATH}\n")
        return

    if args.dry_run:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"status": "DRY_RUN_ok_not_submitted", "verify_err": err,
                        "counts_a": counts_a, "counts_b": counts_b,
                        "n1q_a": n1q_a, "n1q_b": n1q_b, "n2q_a": n2q_a, "n2q_b": n2q_b}, f, indent=2)
        print(f"\n  DRY RUN: verified equivalent, reduction found, NOT submitted (no money spent). See {RESULTS_PATH}\n")
        return

    qc_a_meas = qc_a.copy()
    qc_a_meas.measure_all()
    qc_b_meas = qc_b.copy()
    qc_b_meas.measure_all()

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")
    print(f"\n  connected to {backend.name}")

    results = {"verify_err": err, "counts_a": counts_a, "counts_b": counts_b,
               "n1q_a": n1q_a, "n1q_b": n1q_b, "n2q_a": n2q_a, "n2q_b": n2q_b,
               "shots": SHOTS, "backend": BACKEND_NAME, "slot": SLOT, "group": GROUP}

    for tag, qc in [("A_baseline", qc_a_meas), ("B_reduced_1q", qc_b_meas)]:
        print(f"\n  SUBMITTING NOW -- real hardware, real cost will be incurred (circuit {tag}).")
        job = backend.run([qc], shots=SHOTS)
        job_id = job.job_id()
        print(f"  submitted. job_id={job_id}")
        results[f"job_id_{tag}"] = job_id
        with open(RESULTS_PATH, "w") as f:
            json.dump(results, f, indent=2, default=str)

        job.wait_for_final_state(timeout=14400)
        status = str(job.status())
        print(f"  final status: {status}")
        counts = None
        if job.done():
            counts = dict(job.result().get_counts())
            print(f"  REAL counts: {counts}")
        raw = backend.client.retrieve_job(job_id)
        results[f"status_{tag}"] = status
        results[f"counts_{tag}"] = counts
        results[f"raw_job_metadata_{tag}"] = raw
        with open(RESULTS_PATH, "w") as f:
            json.dump(results, f, indent=2, default=str)

    print(f"\n  Full results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
