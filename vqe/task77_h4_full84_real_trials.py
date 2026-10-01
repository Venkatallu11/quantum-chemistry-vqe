#!/usr/bin/env python3
"""
task77_h4_full84_real_trials.py -- iteration 77. The other concurrent
session's task73_hw_analysis rehearsal (bootstrap-resampled from
task59's old checkpoint) predicts the FULL 84-circuit design at 1,100
shots/circuit clears chemical accuracy on all three backends with
q95 < 1.0 kcal/mol -- and it's now affordable: $5,352.16 real approved
budget (confirmed: $2,352.16 + $3,000 the other session found approved),
vs. a real, independently-verified cost of $3,875-3,882 (81 new circuits
+ 3 reused from task49), margin ~$1,470.

This submits FRESH real trials of the full 84-circuit design to IonQ's
free ionq_simulator (same honest-data discipline as task72: no reuse of
old data, no bootstrap) to cross-check that bootstrap prediction against
genuinely independent real draws -- motivated by task75's own finding
that bootstrap resampling runs optimistic relative to real submissions
(real cross-submission drift exceeds pure shot noise).

Reuses task72's exact infrastructure (ancilla regression check,
GC-aware correction, reproducibility safeguard) with ALL 4 GC groups
measured for every one of the 21 slots (the full design, not the
6-slot-panel subset).

Run:
    PYTHONHASHSEED=0 python vqe/task77_h4_full84_real_trials.py --submit --trial 1
    PYTHONHASHSEED=0 python vqe/task77_h4_full84_real_trials.py --analyze --trial 1
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task72_h4_full_experiment as t72
from task59_h4_gc_no_frame_fit import GATE_NAME, BACKENDS
from task39b_native_ancilla_parity import ancilla_cnots_native, verify_native_ancilla
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list

PROPOSED_SHOTS = 1100


def ckpt_path(trial):
    return os.path.join(os.path.dirname(__file__), f"task77_h4_full84_trial{trial}.partial.json")


def results_path(trial):
    return os.path.join(os.path.dirname(__file__), f"task77_h4_full84_trial{trial}_results.json")


def setup_full84():
    ctx = t72.setup()
    # override: EVERY slot measures ALL groups (the full 84-circuit design)
    ctx["measured_groups_map"] = {name: list(range(ctx["n_groups"])) for name in ctx["kept"]}
    n_circuits = sum(len(v) for v in ctx["measured_groups_map"].values())
    print(f"  FULL design: {n_circuits} circuits/backend (all {ctx['n_groups']} groups x {len(ctx['kept'])} slots)")
    return ctx


def submit(ctx, trial, canary=False):
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    diff = verify_native_ancilla(ancilla_native)
    print(f"  ancilla circuit regression check: diff={diff:.3e}  {'PASS' if diff < 1e-8 else 'FAIL -- STOP'}")
    assert diff < 1e-8, "ancilla circuit failed its own regression check -- stop before submitting"

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}, shots={PROPOSED_SHOTS}/circuit, trial={trial}")

    backends_to_run = ["ideal"] if canary else BACKENDS
    slots_to_run = [ctx["kept"][0]] if canary else ctx["kept"]

    path = ckpt_path(trial)
    state = {"done": {}}
    if os.path.exists(path):
        with open(path) as f:
            state = json.load(f)

    for backend_name in backends_to_run:
        for name in slots_to_run:
            key = f"{backend_name}|{name}"
            if key in state["done"]:
                print(f"    skip (already done): {key}")
                continue
            group_idxs = ctx["measured_groups_map"][name]
            from task28d_all_gate_zne import optimized_native_circuit
            register_native = optimized_native_circuit(ctx["fixed_solutions"][name]["angles"], GATE_NAME)
            circuits = t72.build_circuits_for_slot(register_native, ancilla_native, ctx["diagonalizers"],
                                                     ctx["diag_natives"], group_idxs)
            job = submit_job(circuits, backend, backend_name, shots=PROPOSED_SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"counts": counts, "group_idxs": group_idxs}
            with open(path, "w") as f:
                json.dump(state, f, indent=2)
            print(f"    done: {key} ({len(circuits)} circuits, groups={group_idxs}, {PROPOSED_SHOTS} shots each)")

    print(f"\n  {'CANARY PASSED' if canary else 'FULL SWEEP COMPLETE'} -- saved -> {path}")


def analyze(ctx, trial):
    path = ckpt_path(trial)
    with open(path) as f:
        state = json.load(f)
    missing = [bn for bn in BACKENDS if any(f"{bn}|{name}" not in state["done"] for name in ctx["kept"])]
    if missing:
        print(f"  ERROR: real checkpoint incomplete for backends {missing} -- run --submit first.")
        return None

    print(f"\n  H4 exact_energy={ctx['p']['exact_energy']:.6f} Ha")
    results = {bn: t72.analyze_backend(ctx, bn, state) for bn in BACKENDS}
    results2 = {bn: t72.analyze_backend(ctx, bn, state) for bn in BACKENDS}
    for bn in BACKENDS:
        for k in ("raw_err_kcal", "no_frame_err_kcal", "shared_frame_err_kcal"):
            if abs(results[bn][k] - results2[bn][k]) > 1e-9:
                raise RuntimeError(f"REPRODUCIBILITY CHECK FAILED for trial {trial}, {bn}.{k}")
    print("  reproducibility check: two independent passes agree exactly -- PASS")

    for bn in BACKENDS:
        e = results[bn]
        print(f"    {bn}: accept={e['mean_accept']:.4f}  raw={e['raw_err_kcal']:+.4f}  "
              f"no_frame={e['no_frame_err_kcal']:+.4f}  shared_frame={e['shared_frame_err_kcal']:+.4f} "
              f"kcal/mol  chi2/dof={e['chi2_dof']:.4f}")

    with open(results_path(trial), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {results_path(trial)}\n")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--canary", action="store_true")
    ap.add_argument("--trial", type=int, default=1)
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = setup_full84()
    if args.submit or args.canary:
        submit(ctx, args.trial, canary=args.canary)
        if not args.canary:
            analyze(ctx, args.trial)
    elif args.analyze:
        analyze(ctx, args.trial)


if __name__ == "__main__":
    main()
