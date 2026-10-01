#!/usr/bin/env python3
"""
task74_hw_submission.py -- iteration 74. Builds, checks and (only when asked)
submits the real-hardware H4 run analysed by task73.

DESIGN (default "full"): all 4 GC groups on all 21 slots = 84 circuits. The
3 GC3 circuits on u_0/u_1/u_2 already ran on Forte Enterprise 1 (task49,
10 Sep 2026) and are NOT resubmitted, so 81 new circuits go in one bundled
job at 1,100 shots. "subset" = task71's 54-circuit plan (51 new).

CIRCUITS: built by task59's own build_gc_measurement_circuits -- the exact
construction that produced the rehearsal data (native state prep, parity
ancilla, native GC diagonalizer, measure last).

PRE-SUBMISSION CHECKS (all must pass; no network needed):
  1. parity-ancilla native circuit reproduces the abstract CX version
  2. Schmidt sign gate (task72): stored ideal data fits at the exact frame
  3. reused circuits: today's rebuild of task49's 3 circuits has identical
     gate counts
  4. every circuit, simulated noiselessly and postselected on ancilla=0,
     reproduces the exact expectation values <t|P|t> of its target state for
     every label in its group (catches wrong angles, bit order or
     diagonalizer; NOT a Schmidt sign flip, which moves circuits and targets
     together -- check 2 guards that)
  5. estimated cost <= --max-cost

MODES
  (default)                      dry run: build + check + cost; nothing sent
  --submit --confirm-cost X      submit; X must equal the printed estimate
                                 (rounded to the dollar), a deliberate
                                 second confirmation that real money is spent
  --fetch JOB_ID                 retrieve a finished job and write
                                 task73_hw_counts.json for task73

Needs IONQ_API_KEY in the environment (.env) only for --submit / --fetch.

Run:
    PYTHONHASHSEED=0 python vqe/task74_hw_submission.py
    PYTHONHASHSEED=0 python vqe/task74_hw_submission.py --submit --confirm-cost 3882
    PYTHONHASHSEED=0 python vqe/task74_hw_submission.py --fetch <job_id>
    PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py --counts vqe/task73_hw_counts.json
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qiskit.quantum_info import Statevector, Pauli
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, verify_native_ancilla
from task59_h4_gc_no_frame_fit import build_gc_measurement_circuits, CKPT_PATH
from task72_rehearsal_sign_and_chi2_check import exact_frame_gate, EXACT_FRAME_GATE
from task73_hw_analysis import Context, key, REUSED_TASK49_PATH, K
from general_commuting_measurements import expectations_from_counts

HERE = os.path.dirname(__file__)
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOTS = 1100
RATE_1Q = 0.0001645      # $ per 1q gate-shot (task71's precise rate; task68 card is 0.000164)
RATE_2Q = 0.001121       # $ per 2q gate-shot
JOB_FLOOR = 25.7899      # applies once per bundled job
BUDGET_TOTAL = 5352.16   # $2,352.16 remaining + $3,000 approved for the full run
IDEAL_TOL = 1e-6
MANIFEST_PATH = os.path.join(HERE, "task74_hw_submission_manifest.json")
COUNTS_OUT = os.path.join(HERE, "task73_hw_counts.json")


def build_all(ctx):
    """Every circuit the design needs, keyed '<slot>|g<gi>', in task59's construction."""
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    circuits = {}
    for name in ctx.kept:
        register = optimized_native_circuit(ctx.fixed_solutions[name]["angles"], GATE_NAME)
        per_group = build_gc_measurement_circuits(register, ancilla_native, ctx.diagonalizers, ctx.diag_natives)
        for gi in ctx.measured[name]:
            circuits[key(name, gi)] = per_group[gi]
    return circuits, ancilla_native


def gate_counts(qc):
    ops = dict(qc.count_ops())
    return {"n1q": int(ops.get("gpi", 0) + ops.get("gpi2", 0)), "n2q": int(ops.get(GATE_NAME, 0)), "ops": ops}


def cost(counts_list, shots):
    total = sum(shots * (c["n1q"] * RATE_1Q + c["n2q"] * RATE_2Q) for c in counts_list)
    return max(total, JOB_FLOOR)


def ideal_check(circuits, ctx):
    """Noiseless simulation of each circuit vs the exact target expectation values."""
    worst, worst_key = 0.0, None
    for k, qc in circuits.items():
        name, g = k.split("|g")
        gi = int(g)
        probs = Statevector(qc.remove_final_measurements(inplace=False)).probabilities_dict()
        kept = {bs[1:]: p for bs, p in probs.items() if bs[0] == "0"}
        exp = expectations_from_counts(kept, ctx.diagonalizers[gi])
        t = np.asarray(ctx.p["targets"][name])
        for l in ctx.groups[gi]:
            exact = float(np.real(t.conj() @ Pauli(l).to_matrix() @ t))
            dev = abs(exp[l] - exact)
            if dev > worst:
                worst, worst_key = dev, f"{k}:{l}"
    return worst, worst_key


def reused_check(circuits):
    with open(REUSED_TASK49_PATH) as f:
        d = json.load(f)
    gi = d["hard_group_index"]
    out = {}
    for slot in d["target_slots"]:
        ops = gate_counts(circuits[key(slot, gi)])["ops"]
        stored = d["gate_counts_per_slot"][slot]
        out[key(slot, gi)] = all(ops.get(g, 0) == n for g, n in stored.items())
    return out


def prepare(design, max_cost):
    ctx = Context(design)
    circuits, ancilla_native = build_all(ctx)
    print(f"  design={design}: {len(circuits)} circuits in the analysis")

    diff = verify_native_ancilla(ancilla_native)
    print(f"  check 1 (ancilla circuit): diff={diff:.2e}")
    assert diff < 1e-8, "ancilla circuit regression failed -- STOP"

    with open(CKPT_PATH) as f:
        gate = exact_frame_gate(json.load(f), ctx.kept, ctx.groups, ctx.diagonalizers, ctx.P_S, K, ctx.non_id)
    print(f"  check 2 (Schmidt signs): exact-frame residual {gate:.2e}")
    assert gate < EXACT_FRAME_GATE, "Schmidt sign convention mismatch -- STOP"

    reused = reused_check(circuits)
    print(f"  check 3 (reused task49 circuits identical today): {reused}")
    assert all(reused.values()), "today's circuits differ from task49's -- cannot reuse, STOP"

    worst, where = ideal_check(circuits, ctx)
    print(f"  check 4 (noiseless circuits reproduce exact values): worst |dev| {worst:.2e} at {where}")
    assert worst < IDEAL_TOL, "a circuit does not prepare/measure its target state -- STOP"

    submit_keys = [k for k in circuits if k not in reused]
    counts = {k: gate_counts(circuits[k]) for k in submit_keys}
    est = cost(list(counts.values()), SHOTS)
    print(f"  check 5 (cost): {len(submit_keys)} new circuits x {SHOTS} shots -> ${est:,.2f} "
          f"(max ${max_cost:,.2f}; total budget ${BUDGET_TOTAL:,.2f}, margin ${BUDGET_TOTAL - est:,.2f})")
    assert est <= max_cost, "estimated cost exceeds --max-cost -- STOP"
    return ctx, circuits, submit_keys, counts, est


def submit(circuits, submit_keys, counts, est, design):
    from ionq_backend import connect_provider
    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")
    print(f"  connected to {backend.name}")
    manifest = {"backend": BACKEND_NAME, "shots": SHOTS, "design": design, "circuit_keys": submit_keys,
                "gate_counts": counts, "estimated_cost_usd": round(est, 2), "status": "submitting"}
    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n  SUBMITTING NOW -- real hardware, real cost (~${est:,.2f}), {len(submit_keys)} circuits, one job")
    job = backend.run([circuits[k] for k in submit_keys], shots=SHOTS)
    manifest["job_id"] = job.job_id()
    manifest["status"] = "submitted"
    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  submitted. job_id={manifest['job_id']}  (saved -> {MANIFEST_PATH})")
    print(f"  fetch later with: python vqe/task74_hw_submission.py --fetch {manifest['job_id']}")


def fetch(job_id):
    from ionq_backend import connect_provider
    with open(MANIFEST_PATH) as f:
        manifest = json.load(f)
    assert manifest.get("job_id") == job_id, "job id does not match the saved manifest -- STOP"
    provider = connect_provider()
    backend = provider.get_backend(manifest["backend"], gateset="native")
    job = backend.retrieve_job(job_id)
    job.wait_for_final_state(timeout=4 * 3600)
    assert job.done(), f"job not done: {job.status()}"
    counts_list = job.result().get_counts()
    if not isinstance(counts_list, list):
        counts_list = [counts_list]
    keys = manifest["circuit_keys"]
    assert len(counts_list) == len(keys), f"{len(counts_list)} results for {len(keys)} circuits -- STOP"
    raw = backend.client.retrieve_job(job_id)
    out = {"job_id": job_id, "backend": manifest["backend"], "shots": manifest["shots"],
           "design": manifest["design"], "submitted_at": raw.get("submitted_at"),
           "completed_at": raw.get("completed_at"), "billed_cost_usd": raw.get("cost"),
           "estimated_cost_usd": manifest["estimated_cost_usd"],
           "counts": {k: dict(c) for k, c in zip(keys, counts_list)}, "raw_job_metadata": raw}
    with open(COUNTS_OUT, "w") as f:
        json.dump(out, f, indent=2, default=str)
    manifest["status"] = "fetched"
    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  {len(keys)} circuits' counts saved -> {COUNTS_OUT}")
    print(f"  next: PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py --counts {COUNTS_OUT} "
          f"--design {manifest['design']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=["full", "subset"], default="full")
    ap.add_argument("--max-cost", type=float, default=BUDGET_TOTAL)
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--confirm-cost", type=int, help="must equal the printed estimate, rounded to the dollar")
    ap.add_argument("--fetch", metavar="JOB_ID")
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0")

    if args.fetch:
        fetch(args.fetch)
        return
    ctx, circuits, submit_keys, counts, est = prepare(args.design, args.max_cost)
    if not args.submit:
        print("\n  DRY RUN -- all checks passed, nothing submitted.")
        print(f"  to submit: --submit --confirm-cost {round(est)}")
        return
    assert args.confirm_cost == round(est), f"--confirm-cost must be {round(est)} -- STOP"
    submit(circuits, submit_keys, counts, est, args.design)


if __name__ == "__main__":
    main()
