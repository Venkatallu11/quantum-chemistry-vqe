#!/usr/bin/env python3
"""
task39c_ancilla_real_submission.py -- iteration 39, Task C. Real
submission of the native ancilla-parity circuit (Task 39B, verified
locally to machine precision) combined with the current production
circuit (`optimized_native_circuit`, IonQ's own TrappedIonOptimizerPlugin
included -- matching how `task28b_optimized_raw.json`'s real baseline
data was actually collected, so the two are directly comparable without
resubmitting a fresh no-ancilla control). GATE_NAME="zz" throughout,
matching this session's own established convention (works on all
backends via the free simulator regardless of a backend's own "native"
gate family).

SHOT BUDGET, disclosed: 20,000/circuit (not this project's usual
100,000) and aria-1/forte-1 ONLY (no `ideal` -- already verified exact
locally in Task 39B) -- a deliberate reduction given today's earlier
real IonQ API quota failure (`TooManyShots`, task38_readout_
calibration.py) after a large cumulative real-submission volume. This
still gives a real, usable (if less precise) measurement, and per-slot
incremental checkpointing means a quota failure partway through loses
nothing already collected.

CANARY MODE (--canary): submits just 1 slot's 13 circuits to forte-1
only, to confirm the API is currently accepting submissions before
committing to the full 21-slot x 2-backend sweep.

Run:
    PYTHONHASHSEED=0 python vqe/task39c_ancilla_real_submission.py --canary
    PYTHONHASHSEED=0 python vqe/task39c_ancilla_real_submission.py
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native, verify_native_ancilla
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list

GATE_NAME = "zz"
SHOTS = 20_000
BACKENDS = ["aria-1", "forte-1"]
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
CKPT_PATH_TEMPLATE = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission{suffix}.partial.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission_results.json")


def load_groups_from_baseline(kept):
    """Reuse the EXACT same (slot, group) structure task28b's real
    baseline used, so the ancilla submission's circuits are directly,
    label-for-label comparable to the already-collected no-ancilla data."""
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    tags = raw_ck["tags"]["forte-1"]
    per_name = {}
    for name, group in tags:
        per_name.setdefault(name, [])
        gt = tuple(group)
        if gt not in [tuple(g) for g in per_name[name]]:
            per_name[name].append(list(group))
    return {name: per_name[name] for name in kept}


def build_ancilla_measurement_circuits(register_native_4q, ancilla_native_5q, groups):
    circuits = []
    for group in groups:
        combined = effrag_mod.combined_basis_label(group)
        basis_qc = native_basis_change(combined, GATE_NAME)  # 4-qubit basis change, ancilla untouched
        full5 = with_ancilla_parity_native(register_native_4q, ancilla_native_5q)
        full5 = full5.compose(basis_qc, qubits=[0, 1, 2, 3])
        full5.measure_all()
        circuits.append(full5)
    return circuits


def load_partial(ckpt_path):
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            return json.load(f)
    return {"done": {}}


def save_partial(state, ckpt_path):
    with open(ckpt_path, "w") as f:
        json.dump(state, f, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    ap.add_argument("--draw", type=int, default=0, help="replication draw index -- 0 is the original "
                     "Task C sweep (task39c_ancilla_real_submission.partial.json, unchanged filename); "
                     "1+ writes to a separate task39c_ancilla_real_submission_drawN.partial.json so "
                     "independent replication draws never overwrite each other or the original")
    args = ap.parse_args()
    ckpt_path = CKPT_PATH_TEMPLATE.format(suffix="" if args.draw == 0 else f"_draw{args.draw}")

    print("\n" + "=" * 96)
    print(f"  task39c_ancilla_real_submission.py -- {'CANARY (1 slot, forte-1 only)' if args.canary else f'FULL SWEEP (draw {args.draw})'}")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=6, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(6)
    groups_by_slot = load_groups_from_baseline(kept)

    ancilla_native = ancilla_cnots_native(GATE_NAME)
    diff = verify_native_ancilla(ancilla_native)
    print(f"  ancilla circuit regression check: diff={diff:.3e}  {'PASS' if diff < 1e-8 else 'FAIL -- STOP'}")
    if diff >= 1e-8:
        raise RuntimeError("ancilla circuit failed its own regression check -- stop before submitting")

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    backends_to_run = ["forte-1"] if args.canary else BACKENDS
    slots_to_run = [kept[0]] if args.canary else kept

    state = load_partial(ckpt_path)
    for backend_name in backends_to_run:
        for name in slots_to_run:
            key = f"{backend_name}|{name}"
            if key in state["done"]:
                print(f"    skip (already done): {key}")
                continue
            register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
            circuits = build_ancilla_measurement_circuits(register_native, ancilla_native, groups_by_slot[name])
            job = submit_job(circuits, backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"groups": groups_by_slot[name], "counts": counts}
            save_partial(state, ckpt_path)
            print(f"    done: {key} ({len(circuits)} circuits, {SHOTS} shots each)")

    print(f"\n  {'CANARY PASSED -- API accepting submissions, safe to run full sweep' if args.canary else 'FULL SWEEP COMPLETE'}")
    print(f"  partial checkpoint saved -> {ckpt_path}")


if __name__ == "__main__":
    main()
