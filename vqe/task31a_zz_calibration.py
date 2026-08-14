#!/usr/bin/env python3
"""
task31a_zz_calibration.py -- iteration 31, Task A. BLOCKING. Resolves the
ZZ calibration conversion bug and re-measures p2(zz) IN CONTEXT, per gate
position, on the real H4 circuit -- not an isolated 2-qubit probe.
============================================================================
PART 1 -- THE CONVERSION BUG, verified mathematically and against
qiskit_aer's own documented convention before touching any code.
`fixed_ansatz.P2_PER_GATE = 1 - 0.9952 = 0.0048` is defined as raw
INFIDELITY (1-F) and passed DIRECTLY as the depolarizing parameter to
both `qiskit_aer.noise.depolarizing_error(param, n)` (confirmed via its
own docstring: E(rho)=(1-lambda)rho+lambda*I/2^n, lambda IS the
parameter) and `loop_pec.py::depolarizing_weights(p,n)` (confirmed
algebraically identical to the same convention). For this channel,
average gate fidelity F_avg = 1 - lambda*(d-1)/d, d=2^n=4 for 2 qubits --
so the CORRECT depolarizing parameter is
    lambda = (1-F_avg) * d/(d-1) = 0.0048 * 4/3 = 0.0064,
NOT 0.0048 directly. A real, confirmed 1.333x error, present in every
local-noise-model script that has ever imported P2_PER_GATE. NOT silently
rewritten into `fixed_ansatz.py` here (would retroactively alter the
documented meaning of every historical result using it without those
results being re-run) -- reported explicitly, with the corrected value
computed and used in any NEW work this task builds.

PART 2 -- re-measure p2(zz) IN CONTEXT, per gate position. Iteration 9
already showed a reduced/isolated probe under-characterizes the full
circuit. Uses slot `u_2` (11 ZZ gates -- matches this task's own "11 gate
positions" reference exactly, unlike most other slots which the
optimizer reduced to 4-9). VERIFIED (not assumed) at the OPERATOR level
before use: ZZGate(0.25)^4 = I up to global phase (matrix power check,
independent of any specific input state) -- so substituting ANY ONE gate
occurrence with N repetitions for N in {1,5,9,13} preserves the EXACT
final circuit state regardless of what real, non-trivial state precedes
that gate position, unlike the isolated-probe design this project has
used before. 11 positions x 4 N-values = 44 real circuits, forte-1 +
ideal only (aria-1 retired per this task's own instruction).

Run:
    python vqe/task31a_zz_calibration.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from qiskit_ionq.ionq_gates import ZZGate
from qiskit.quantum_info import Statevector, Pauli

K = 6
SHOTS = 300_000
N_SEEDS = 16  # doubled from Task 30B's 8, for a real confidence interval per position
N_LADDER = [1, 5, 9, 13]
MODELS = ["ideal", "forte-1"]  # aria-1 retired this iteration
SLOT = "u_2"  # the one kept slot with exactly 11 ZZ gates
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task31a_zz_calibration.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31a_zz_calibration_results.json")


def substitute_gate_with_ladder(qc, gate_idx, n_reps):
    """Returns a NEW circuit identical to qc, except the single gate at
    gate_idx is replaced by n_reps copies of that same gate."""
    new_qc = qc.copy_empty_like()
    for i, instr in enumerate(qc.data):
        if i == gate_idx:
            for _ in range(n_reps):
                new_qc.append(instr.operation, instr.qubits, instr.clbits)
        else:
            new_qc.append(instr.operation, instr.qubits, instr.clbits)
    return new_qc


def best_observable(sv, candidate_labels):
    best_l, best_v = None, 0.0
    for l in candidate_labels:
        v = float(np.real(sv.expectation_value(Pauli(l))))
        if abs(v) > abs(best_v):
            best_l, best_v = l, v
    return best_l, best_v


def main():
    print("\n" + "=" * 96)
    print("  task31a_zz_calibration.py -- BLOCKING: conversion fix + per-position ZZ calibration")
    print("=" * 96)

    # -- PART 1: the conversion --
    F = 0.9952
    d = 4
    p2_old = 1 - F
    p2_correct = (1 - F) * d / (d - 1)
    print(f"\n  PART 1 -- conversion check:")
    print(f"    P2_PER_GATE (current, treated as depolarizing param directly) = {p2_old:.4f}")
    print(f"    correct depolarizing param, F_avg = 1 - lambda*(d-1)/d, d=4    = {p2_correct:.4f}")
    print(f"    ratio (correct/current) = {p2_correct/p2_old:.4f}")
    print(f"    Task 30B's real measured p2(zz) = 0.0140")
    print(f"    ratio measured/current  = {0.014/p2_old:.3f}x   ratio measured/corrected = {0.014/p2_correct:.3f}x")
    print(f"    -- the conversion fix explains PART of the gap (2.92x -> 2.19x), not all of it")

    # -- PART 2: per-position, in-context ZZ calibration --
    print(f"\n  PART 2 -- per-position calibration on slot={SLOT}'s real H4 circuit")
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    base = optimized_native_circuit(fixed_solutions[SLOT]["angles"], "zz")
    zz_positions = [i for i, instr in enumerate(base.data) if instr.operation.name == "zz"]
    print(f"  {len(zz_positions)} ZZ gate positions found: {zz_positions}")
    assert len(zz_positions) == 11, f"expected 11, got {len(zz_positions)} -- slot choice assumption broken"

    # -- verify operator-level periodicity holds regardless of input state, then verify per-position --
    print("  verifying period-4 substitution preserves the exact final state at EVERY position...")
    calib_points = []
    for gi in zz_positions:
        sv_ref = np.asarray(Statevector.from_instruction(substitute_gate_with_ladder(base, gi, 1)))
        worst_err = 0.0
        for N in N_LADDER:
            sv = np.asarray(Statevector.from_instruction(substitute_gate_with_ladder(base, gi, N)))
            idx = int(np.argmax(np.abs(sv_ref)))
            phase = sv[idx] / sv_ref[idx] if abs(sv_ref[idx]) > 1e-9 else 1.0
            err = float(np.max(np.abs(sv / phase - sv_ref)))
            worst_err = max(worst_err, err)
        assert worst_err < 1e-10, f"position {gi}: periodicity check failed, err={worst_err:.3e}"
        sv1 = Statevector.from_instruction(substitute_gate_with_ladder(base, gi, 1))
        all_4q_labels = ["".join(t) for t in __import__("itertools").product("IXYZ", repeat=4) if t != tuple("IIII")]
        l, v = best_observable(sv1, all_4q_labels)
        calib_points.append((gi, l, v))
    print("  all 11 positions PASS periodicity (<1e-10)")
    for gi, l, v in calib_points:
        print(f"    position {gi}: observable={l}  ideal={v:.4f}")

    if os.path.exists(CKPT_PATH):
        print("\n  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"\n  connected, backend={backend.name}")

        all_jobs = {}
        for model in MODELS:
            circuits, tags = [], []
            for gi, l, v in calib_points:
                for N in N_LADDER:
                    qc = substitute_gate_with_ladder(base, gi, N)
                    qc.measure_all()
                    circuits.append(qc)
                    tags.append((gi, l, N))
            job = submit_job(circuits, backend, model, shots=SHOTS)
            all_jobs[model] = (job, tags)
            print(f"    submitted model={model}: {len(circuits)} circuits, job_id={job.job_id()}")

        all_counts = {}
        for model, (job, tags) in all_jobs.items():
            all_counts[model] = get_counts_list(job)
            print(f"    retrieved model={model}")

        ck = {"tags": {m: v[1] for m, v in all_jobs.items()}, "counts": all_counts}
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    print(f"\n  -- fitting p_ZZ per position via log-linear decay ({N_SEEDS} seeds for a real CI) --")
    results_by_model = {}
    for model in MODELS:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        by_point = {}
        for (gi, l, N), counts in zip(tags, counts_list):
            by_point.setdefault((gi, l), {})[N] = counts

        p_by_position = {}
        for (gi, l), by_N in by_point.items():
            seed_ps = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task31a", model, gi, seed))
                measured = []
                for N in N_LADDER:
                    resampled = bootstrap_counts(by_N[N], SHOTS, rng)
                    m = expectation_from_counts(resampled, l)
                    measured.append(m)
                measured = np.array(measured)
                safe = np.clip(np.abs(measured), 1e-6, None)
                slope, _ = np.polyfit(N_LADDER, np.log(safe), 1)
                p_by_position.setdefault(gi, []).append(1 - np.exp(slope))
        p_summary = {}
        for gi in sorted(p_by_position.keys()):
            vals = np.array(p_by_position[gi])
            mean_p = float(np.mean(vals))
            ci_lo, ci_hi = float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
            p_summary[gi] = {"mean": mean_p, "std": float(np.std(vals)), "ci95_lo": ci_lo, "ci95_hi": ci_hi}
            print(f"    {model:<8} position {gi:>2}: p={mean_p:.5f}  95% CI=[{ci_lo:.5f}, {ci_hi:.5f}]")
        results_by_model[model] = p_summary

    print(f"\n  -- SPREAD ACROSS POSITIONS --")
    for model in MODELS:
        vals = [v["mean"] for v in results_by_model[model].values()]
        print(f"    {model}: min={min(vals):.5f}  max={max(vals):.5f}  mean={np.mean(vals):.5f}  "
              f"spread(max/min)={max(vals)/max(min(vals),1e-9):.2f}x")
    position_dependent = None
    if "forte-1" in results_by_model:
        vals = [v["mean"] for v in results_by_model["forte-1"].values()]
        position_dependent = (max(vals) / max(min(vals), 1e-9)) > 2.0
        print(f"\n  VERDICT: position-dependent (spread > 2x)? {position_dependent}")
        if position_dependent:
            print("  -> PEC should use PER-POSITION channels, not one averaged p2(zz)")
        else:
            print("  -> a single averaged p2(zz) is a reasonable approximation across positions")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "conversion": {"p2_old": p2_old, "p2_correct": p2_correct, "measured_p2_zz_task30b": 0.014},
            "per_position": results_by_model, "position_dependent": position_dependent,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
