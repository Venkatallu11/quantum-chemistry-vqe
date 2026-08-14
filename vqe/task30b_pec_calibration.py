#!/usr/bin/env python3
"""
task30b_pec_calibration.py -- iteration 30, Task B, step 1. Learn
p_2q(ZZ, theta) and p_1q(GPi/GPi2, phi) from REAL Clifford-adjacent
calibration circuits submitted to `ionq_simulator` itself -- never the
local noise model, never the H4 target circuits.
============================================================================
FIRST, A CORRECTION TO THIS TASK'S OWN PREMISE: the task assumes "Forte's
ZZ is ARBITRARY-ANGLE" and asks to calibrate p_2q(ZZ, theta) over a real
angle distribution. Checked directly against the actual optimized H4
circuit (Task 28B, forte/zz family, all 21 kept slots): EVERY ZZ gate in
this circuit uses theta=0.25 EXACTLY (196/196 instances, verified by
extracting every 'zz' instruction's parameter). This ansatz's entangling
angle is architecturally FIXED, not target-dependent -- only the
single-qubit GPi/GPi2 gates carry target-specific rotation angles (181
distinct GPi phi values / 155 distinct GPi2 phi values across the same
196-gate circuit). Angle-conditioning therefore matters for the 1-QUBIT
gates, not the 2-qubit gate -- there is exactly one ZZ angle to calibrate.
This is reported as a correction, not silently substituted.

CALIBRATION DESIGN: repeated same-gate application from a fixed input,
using VERIFIED periodicity (not assumed) so the IDEAL target is EXACTLY
CONSTANT across the calibration ladder -- any real deviation growing with
N is pure accumulated noise, the standard calibration principle, and the
SAME principle `cx_decay_circuit`'s "H then N CX" pattern already uses
locally. Checked numerically before use: H(q0) then N x ZZGate(0.25)(q0,q1),
and N x GPi(phi)(q0) or GPi2(phi)(q0) from |0>, ALL return to the exact
N=1 state (fidelity 1.0) at N in {1,5,9,13} -- a clean period-4 ladder,
verified for every angle used below, not assumed to generalize from one
spot-check.

SCOPE: the 2q gate needs one calibration (theta=0.25, exact). The 1q
gates are calibrated at 4 representative phi values each for GPi and
GPi2, chosen from the REAL circuit's own phi percentiles (10/30/70/90 %ile
region centers) -- "calibrate only the gates you actually use, over the
angle distribution the circuits use," not a generic sweep.

Run:
    python vqe/task30b_pec_calibration.py
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
from qiskit_ionq.ionq_gates import ZZGate, GPIGate, GPI2Gate
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Pauli

K = 6
SHOTS = 300_000
N_SEEDS = 8
N_LADDER = [1, 5, 9, 13]
MODELS = ["ideal", "aria-1", "forte-1"]
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task30b_pec_calibration.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")


def verify_periodicity(build_fn, angle):
    sv_ref = np.asarray(Statevector.from_instruction(build_fn(angle, 1)))
    worst = 0.0
    for N in N_LADDER:
        sv = np.asarray(Statevector.from_instruction(build_fn(angle, N)))
        idx = int(np.argmax(np.abs(sv_ref)))
        phase = sv[idx] / sv_ref[idx] if abs(sv_ref[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv / phase - sv_ref)))
        worst = max(worst, err)
    assert worst < 1e-10, f"periodicity check failed for angle={angle}: worst={worst:.3e}"
    return worst


def build_zz_circuit(theta, N):
    """NATIVE gates only -- the native simulator backend rejects 'h'.
    GPi2(0) creates the same kind of superposition on qubit 0 natively."""
    qc = QuantumCircuit(2)
    qc.append(GPI2Gate(0.0), [0])
    for _ in range(N):
        qc.append(ZZGate(theta), [0, 1])
    return qc


def build_gpi_circuit(phi, N):
    qc = QuantumCircuit(1)
    for _ in range(N):
        qc.append(GPIGate(phi), [0])
    return qc


def build_gpi2_circuit(phi, N):
    qc = QuantumCircuit(1)
    for _ in range(N):
        qc.append(GPI2Gate(phi), [0])
    return qc


def best_observable(build_fn, angle, candidate_labels):
    sv = Statevector.from_instruction(build_fn(angle, 1))
    best_l, best_v = None, 0.0
    for l in candidate_labels:
        v = float(np.real(sv.expectation_value(Pauli(l))))
        if abs(v) > abs(best_v):
            best_l, best_v = l, v
    return best_l, best_v


def main():
    print("\n" + "=" * 96)
    print("  task30b_pec_calibration.py -- real angle-conditioned p_2q/p_1q from ionq_simulator")
    print("=" * 96)

    # -- extract the REAL angle distribution from the actual optimized H4 circuit --
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    zz_angles, gpi_phis, gpi2_phis = [], [], []
    for name in kept:
        qc = optimized_native_circuit(fixed_solutions[name]["angles"], "zz")
        for instr in qc.data:
            op = instr.operation
            if op.name == "zz":
                zz_angles.append(float(op.params[0]))
            elif op.name == "gpi":
                gpi_phis.append(float(op.params[0]) % 1.0)
            elif op.name == "gpi2":
                gpi2_phis.append(float(op.params[0]) % 1.0)
    zz_angles = np.array(zz_angles)
    gpi_phis = np.array(gpi_phis)
    gpi2_phis = np.array(gpi2_phis)
    n_distinct_zz = len(set(np.round(zz_angles, 6)))
    print(f"  REAL circuit angle census: {len(zz_angles)} zz gates, {n_distinct_zz} distinct angle(s) "
          f"(theta={zz_angles[0]:.4f}) -- {'CONFIRMS single fixed angle, correcting this task premise' if n_distinct_zz == 1 else 'multiple angles found'}")
    print(f"  {len(gpi_phis)} gpi gates, {len(set(np.round(gpi_phis,3)))} distinct phi values (range "
          f"[{gpi_phis.min():.3f}, {gpi_phis.max():.3f}])")
    print(f"  {len(gpi2_phis)} gpi2 gates, {len(set(np.round(gpi2_phis,3)))} distinct phi values (range "
          f"[{gpi2_phis.min():.3f}, {gpi2_phis.max():.3f}])")

    zz_theta = float(zz_angles[0])
    gpi_bins = [float(x) for x in np.percentile(gpi_phis, [10, 30, 70, 90])]
    gpi2_bins = [float(x) for x in np.percentile(gpi2_phis, [10, 30, 70, 90])]
    print(f"\n  calibration plan: ZZ at theta={zz_theta:.4f} (the only angle used); "
          f"GPi at phi={[round(x,3) for x in gpi_bins]}; GPi2 at phi={[round(x,3) for x in gpi2_bins]}")

    # -- verify periodicity for EVERY angle to be calibrated, before submitting anything --
    print("\n  -- verifying period-4 ladder (N=1,5,9,13 all return the exact same state) for every angle --")
    verify_periodicity(build_zz_circuit, zz_theta)
    for phi in gpi_bins:
        verify_periodicity(build_gpi_circuit, phi)
    for phi in gpi2_bins:
        verify_periodicity(build_gpi2_circuit, phi)
    print("  all periodicity checks PASS (<1e-10)")

    # -- pick the best-sensitivity observable for each calibration point --
    calib_points = []
    all_2q_labels = ["".join(pair) for pair in __import__("itertools").product("IXYZ", repeat=2) if pair != ("I", "I")]
    l_zz, v_zz = best_observable(build_zz_circuit, zz_theta, all_2q_labels)
    calib_points.append(("zz", zz_theta, l_zz, v_zz, build_zz_circuit, 2))
    for phi in gpi_bins:
        l, v = best_observable(build_gpi_circuit, phi, ["X", "Y", "Z"])
        calib_points.append(("gpi", phi, l, v, build_gpi_circuit, 1))
    for phi in gpi2_bins:
        l, v = best_observable(build_gpi2_circuit, phi, ["X", "Y", "Z"])
        calib_points.append(("gpi2", phi, l, v, build_gpi2_circuit, 1))
    print("\n  calibration points (gate, angle, best observable, ideal value):")
    for gate, angle, l, v, _, nq in calib_points:
        print(f"    {gate:<5} angle={angle:.4f}  observable={l}  ideal={v:.4f}")

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
            for gate, angle, l, v, build_fn, nq in calib_points:
                for N in N_LADDER:
                    qc = build_fn(angle, N)
                    qc.measure_all()
                    circuits.append(qc)
                    tags.append((gate, angle, l, N))
            job = submit_job(circuits, backend, model, shots=SHOTS)
            all_jobs[model] = (job, tags)
            print(f"    submitted model={model}: {len(circuits)} circuits, job_id={job.job_id()}")

        all_counts = {}
        for model, (job, tags) in all_jobs.items():
            all_counts[model] = get_counts_list(job)
            print(f"    retrieved model={model}")

        ck = {
            "tags": {m: v[1] for m, v in all_jobs.items()}, "counts": all_counts,
            "zz_theta": zz_theta, "gpi_bins": gpi_bins, "gpi2_bins": gpi2_bins,
        }
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- fit p per calibration point via log-linear regression on the decay, exactly
    #    shot_noisy_calibrate's own method (slope of log|measured| vs N -> p = 1-exp(slope)) --
    print(f"\n  -- fitting p per (gate, angle) via log-linear decay regression --")
    learned = {model: {"zz": {}, "gpi": {}, "gpi2": {}} for model in MODELS}
    for model in MODELS:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        by_point = {}
        for (gate, angle, l, N), counts in zip(tags, counts_list):
            by_point.setdefault((gate, angle, l), {})[N] = counts

        for (gate, angle, l), by_N in by_point.items():
            seed_ps = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task30b_calib", model, gate, angle, seed))
                measured = []
                for N in N_LADDER:
                    resampled = bootstrap_counts(by_N[N], SHOTS, rng)
                    m = expectation_from_counts(resampled, l)
                    measured.append(m)
                measured = np.array(measured)
                safe = np.clip(np.abs(measured), 1e-6, None) * np.sign(measured + 1e-12)
                # sign-preserving log fit: fit log|value| vs N, keep the known constant sign
                slope, _ = np.polyfit(N_LADDER, np.log(np.abs(safe)), 1)
                p_seed = 1 - np.exp(slope)
                seed_ps.append(p_seed)
            p_mean, p_std = float(np.mean(seed_ps)), float(np.std(seed_ps))
            n = 2 if gate == "zz" else 1
            key = "zz" if gate == "zz" else gate
            learned[model][key][round(angle, 6)] = {"p": p_mean, "p_std": p_std, "n_qubits": n, "label": l}
            print(f"    {model:<8} {gate:<5} angle={angle:.4f}: p={p_mean:.5f} +/- {p_std:.5f}")

    print("\n" + "=" * 96)
    print("  SUMMARY: real, angle-conditioned learned channel parameters")
    for model in MODELS:
        zz_p = list(learned[model]["zz"].values())[0]["p"] if learned[model]["zz"] else None
        gpi_ps = [v["p"] for v in learned[model]["gpi"].values()]
        gpi2_ps = [v["p"] for v in learned[model]["gpi2"].values()]
        print(f"    {model}: p2(zz,{zz_theta:.3f})={zz_p:.5f}  "
              f"p1(gpi) range=[{min(gpi_ps):.5f},{max(gpi_ps):.5f}] (spread={max(gpi_ps)-min(gpi_ps):.5f})  "
              f"p1(gpi2) range=[{min(gpi2_ps):.5f},{max(gpi2_ps):.5f}] (spread={max(gpi2_ps)-min(gpi2_ps):.5f})")
    print("=" * 96 + "\n")

    with open(RESULTS_PATH, "w") as f:
        json.dump(learned, f, indent=2)
    print(f"  Results saved -> {RESULTS_PATH}\n")
    return learned


if __name__ == "__main__":
    main()
