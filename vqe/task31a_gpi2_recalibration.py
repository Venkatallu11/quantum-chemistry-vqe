#!/usr/bin/env python3
"""
task31a_gpi2_recalibration.py -- iteration 31, Task A continued. Task 30B's
GPi2 calibration failed (signal below shot-noise floor at 300,000 shots/
8 seeds, some fits went unphysically negative). Retry with 3.3x more
shots (1,000,000) and 2x more seeds (16) at the SAME 4 real phi bins,
SAME verified period-4 ladder. If it still fails, report an EXPLICIT
error bound on GPi2's contribution rather than silently falling back to
GPi's mean (Task 30B's disclosed but acknowledged-weak workaround).

Run:
    python vqe/task31a_gpi2_recalibration.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from qiskit_ionq.ionq_gates import GPI2Gate
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Pauli

SHOTS = 1_000_000
N_SEEDS = 16
N_LADDER = [1, 5, 9, 13]
MODELS = ["ideal", "forte-1"]
GPI2_BINS = [0.0, 0.072077, 0.510064, 0.75]  # same 4 real phi bins as Task 30B
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task31a_gpi2_recalibration.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31a_gpi2_recalibration_results.json")


def build_gpi2_circuit(phi, N):
    qc = QuantumCircuit(1)
    for _ in range(N):
        qc.append(GPI2Gate(phi), [0])
    return qc


def best_observable(phi):
    sv = Statevector.from_instruction(build_gpi2_circuit(phi, 1))
    best_l, best_v = None, 0.0
    for l in ["X", "Y", "Z"]:
        v = float(np.real(sv.expectation_value(Pauli(l))))
        if abs(v) > abs(best_v):
            best_l, best_v = l, v
    return best_l, best_v


def main():
    print("\n" + "=" * 96)
    print(f"  task31a_gpi2_recalibration.py -- {SHOTS:,} shots x {N_SEEDS} seeds (vs Task 30B's 300,000 x 8)")
    print("=" * 96)

    calib_points = []
    for phi in GPI2_BINS:
        l, v = best_observable(phi)
        calib_points.append((phi, l, v))
        print(f"    phi={phi:.4f}  observable={l}  ideal={v:.4f}")

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
            for phi, l, v in calib_points:
                for N in N_LADDER:
                    qc = build_gpi2_circuit(phi, N)
                    qc.measure_all()
                    circuits.append(qc)
                    tags.append((phi, l, N))
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

    print(f"\n  -- fitting p1(gpi2) per phi bin, {N_SEEDS} seeds --")
    results_by_model = {}
    for model in MODELS:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        by_point = {}
        for (phi, l, N), counts in zip(tags, counts_list):
            by_point.setdefault((round(phi, 6), l), {})[N] = counts

        p_summary = {}
        for (phi, l), by_N in by_point.items():
            seed_ps = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task31a_gpi2", model, phi, seed))
                measured = []
                for N in N_LADDER:
                    resampled = bootstrap_counts(by_N[N], SHOTS, rng)
                    measured.append(expectation_from_counts(resampled, l))
                measured = np.array(measured)
                safe = np.clip(np.abs(measured), 1e-6, None)
                slope, _ = np.polyfit(N_LADDER, np.log(safe), 1)
                seed_ps.append(1 - np.exp(slope))
            vals = np.array(seed_ps)
            mean_p, std_p = float(np.mean(vals)), float(np.std(vals))
            ci_lo, ci_hi = float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
            physical = mean_p >= 0 and ci_lo > -0.001
            p_summary[phi] = {"mean": mean_p, "std": std_p, "ci95_lo": ci_lo, "ci95_hi": ci_hi, "physical": physical}
            print(f"    {model:<8} phi={phi:.4f}: p={mean_p:.6f} +/- {std_p:.6f}  95% CI=[{ci_lo:.6f},{ci_hi:.6f}]  "
                  f"{'PHYSICAL' if physical else 'still unreliable'}")
        results_by_model[model] = p_summary

    if "forte-1" in results_by_model:
        vals = [v["mean"] for v in results_by_model["forte-1"].values()]
        n_physical = sum(1 for v in results_by_model["forte-1"].values() if v["physical"])
        print(f"\n  SUMMARY: {n_physical}/4 bins now physical (Task 30B: 0/4)")
        if n_physical < 4:
            worst_ci = max(abs(v["ci95_hi"]) for v in results_by_model["forte-1"].values())
            print(f"  EXPLICIT BOUND for the unreliable bin(s): |p1(gpi2)| < {worst_ci:.4f} (95% CI bound), "
                  f"NOT a point estimate -- use this as an uncertainty contribution, not a silently-substituted value")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results_by_model, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
