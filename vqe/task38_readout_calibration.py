#!/usr/bin/env python3
"""
task38_readout_calibration.py -- iteration 38, complementary lever
alongside 38D. Task 38C found p_readout responsible for 24.8% of V_cal
purely because NO real calibration has ever existed for it (weak prior
std=0.01, Task 37B) -- unlike delta_zz/delta_gpi2, real readout-error
characterization is a standard, simple, well-established technique
(prepare a known basis state, measure, the misassignment rate directly
IS the readout error), not a fragile coherent-angle probe. Real
`ionq_simulator` submission, `ideal`/`aria-1`/`forte-1`, N_DRAWS=8
independent real submissions per backend from the START -- Task 37E's
own lesson (a single draw gave a false-positive "significant" reading
that did not replicate) is applied immediately here, not re-learned.

METHOD: p(1|0) = fraction of "1" outcomes measuring a real |0> (no
gates); p(0|1) = fraction of "0" outcomes measuring a real |1> (prepared
via the real native GPI(0) gate, which is exactly the X matrix --
verified below, not assumed). p_readout_hat = mean(p(1|0), p(0|1)),
matching this project's own symmetric-readout-error model
((1-2p)^w attenuation, Task 37C's `readout_attenuation`).

Run:
    PYTHONHASHSEED=0 python vqe/task38_readout_calibration.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list
from qiskit import QuantumCircuit
from qiskit_ionq.ionq_gates import GPIGate

SHOTS = 100_000
BACKENDS = ["ideal", "aria-1", "forte-1"]
N_DRAWS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task38_readout_calibration_results.json")


def frac_one(counts):
    total = sum(counts.values())
    ones = sum(v for k, v in counts.items() if k.strip().endswith("1") or k.strip() == "1")
    return ones / total


def build_zero_circuit():
    qc = QuantumCircuit(1)
    return qc  # |0>, no gates at all


def build_one_circuit():
    qc = QuantumCircuit(1)
    qc.append(GPIGate(0.0), [0])  # verified below: GPI(0) == X exactly
    return qc


def main():
    print("\n" + "=" * 96)
    print(f"  task38_readout_calibration.py -- real symmetric readout-error measurement, {N_DRAWS} draws/backend")
    print("=" * 96)

    gpi0 = GPIGate(0.0).to_matrix()
    x_check = np.max(np.abs(gpi0 - np.array([[0, 1], [1, 0]])))
    print(f"  sanity: GPI(0) vs X, max diff = {x_check:.2e}  {'PASS' if x_check < 1e-12 else 'FAIL -- STOP'}")
    if x_check >= 1e-12:
        raise RuntimeError("GPI(0) is not X -- |1> prep would be wrong, stop before submitting")

    c0 = build_zero_circuit()
    c1 = build_one_circuit()

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    p10_by_backend = {b: [] for b in BACKENDS}
    p01_by_backend = {b: [] for b in BACKENDS}

    for draw in range(N_DRAWS):
        for backend_name in BACKENDS:
            cc0 = c0.copy(); cc0.measure_all()
            cc1 = c1.copy(); cc1.measure_all()
            job = submit_job([cc0, cc1], backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            p10 = frac_one(counts[0])   # measured |0> as 1
            p01 = 1.0 - frac_one(counts[1])  # measured |1> as 0
            p10_by_backend[backend_name].append(p10)
            p01_by_backend[backend_name].append(p01)
        print(f"    draw {draw+1}/{N_DRAWS} done")

    print(f"\n  -- RESULTS: mean +/- std/SEM across {N_DRAWS} independent real submissions --")
    results = {}
    for backend_name in BACKENDS:
        p10 = np.array(p10_by_backend[backend_name])
        p01 = np.array(p01_by_backend[backend_name])
        p10_mean, p10_std = float(p10.mean()), float(p10.std(ddof=1))
        p01_mean, p01_std = float(p01.mean()), float(p01.std(ddof=1))
        p_readout_hat = float((p10_mean + p01_mean) / 2)
        pooled_draws = (p10 + p01) / 2
        p_readout_std = float(pooled_draws.std(ddof=1))
        p_readout_sem = p_readout_std / np.sqrt(N_DRAWS)
        print(f"\n  === {backend_name} ===")
        print(f"    p(1|0) = {p10_mean:.6f} +/- {p10_std:.6f}   p(0|1) = {p01_mean:.6f} +/- {p01_std:.6f}")
        print(f"    p_readout_hat = {p_readout_hat:.6f}   std_across_draws={p_readout_std:.6f}   SEM={p_readout_sem:.6f}")
        results[backend_name] = {
            "p10_mean": p10_mean, "p10_std": p10_std, "p01_mean": p01_mean, "p01_std": p01_std,
            "p_readout_hat": p_readout_hat, "p_readout_std_across_draws": p_readout_std,
            "p_readout_sem": p_readout_sem, "p10_vals": p10.tolist(), "p01_vals": p01.tolist(),
        }

    ideal_ok = results["ideal"]["p_readout_hat"] < 3 * max(results["ideal"]["p_readout_sem"], 1e-4)
    print(f"\n  -- SANITY CHECK: ideal backend should read ~0 -- {'PASS' if ideal_ok else 'FAIL -- do not trust real backends'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"ideal_ok": bool(ideal_ok), "n_draws": N_DRAWS, "shots": SHOTS, "results": results}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
