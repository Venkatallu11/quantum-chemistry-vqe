#!/usr/bin/env python3
"""
task37e_real_calibration_submission.py -- iteration 37, Task E real
follow-up. REAL submission of the 2 dedicated calibration circuits Task
37E proposed and verified analytically -- explicit user go-ahead given
("GO AHEAD DO THE REAL SUBMISSION OF THESE 2 CALIBRATION CIRCUIT").
`ionq_simulator` only (free backend, this project's standing rule --
never real IBM/IonQ hardware), noise_model in {ideal, aria-1, forte-1}.

WHY THIS IS A CLEAN, MINIMAL, DIRECT MEASUREMENT, not an indirect fit:
both circuits were chosen (Task 37E) so the ideal/noiseless prediction is
EXACTLY ZERO, and the analytic sensitivity to the bias parameter is
large and already numerically verified. That means the REAL measured
expectation value on aria-1/forte-1 IS, to good approximation, a direct
readout of the bias itself (measured_value / slope = delta_hat) -- no
model-fitting, no PEC, no frame, nothing to get wrong upstream.
  ZZ circuit:   |++> -> ZZ(theta=0.25) -> measure <XI>. Ideal = 0.
                slope d<XI>/d(delta_zz) at theta=0.25 = -2*pi (verified
                numerically in Task 37E to 4 decimal places).
                Built from REAL native gates (GPI2(0.25) x2 to prepare
                |++>, real ZZGate(0.25)) -- NOT a synthetically-biased
                gate. Any bias is the hardware's/simulator's own.
  GPi2 circuit: |0> -> GPI2(phi=0) -> measure <Z>. Ideal = 0.
                slope d<Z>/d(delta_gpi2) at phi=0 = -1 (verified).

Shot-noise floor, reported alongside every estimate (standard error of a
single Pauli expectation, N=100k shots): SE(<P>) ~= sqrt((1-<P>^2)/N) ~=
1/sqrt(1e5) ~= 3.16e-3 near <P>~0, giving SE(delta_zz) ~= 3.16e-3/(2*pi)
~= 5.0e-4 and SE(delta_gpi2) ~= 3.16e-3/1 ~= 3.16e-3 -- a measured delta
smaller than ~2x this floor is not statistically distinguishable from
zero at this shot budget, disclosed explicitly, not glossed over.

Run:
    PYTHONHASHSEED=0 python vqe/task37e_real_calibration_submission.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from task2_fold_response_dataset import native_basis_change
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, expectation_from_counts
from qiskit import QuantumCircuit
from qiskit_ionq.ionq_gates import GPI2Gate, ZZGate

GATE_NAME = "zz"
SHOTS = 100_000
BACKENDS = ["ideal", "aria-1", "forte-1"]
THETA_TEST = 0.25
SLOPE_ZZ = -2 * np.pi
SLOPE_GPI2 = -1.0
SHOT_NOISE_SE_XI = float(np.sqrt(1.0 / SHOTS))  # (1-<P>^2) ~= 1 near <P>~0
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37e_real_calibration_submission_results.json")


def build_zz_circuit():
    qc = QuantumCircuit(2)
    qc.append(GPI2Gate(0.25), [0])  # |0> -> |+>
    qc.append(GPI2Gate(0.25), [1])  # |0> -> |+>
    qc.append(ZZGate(THETA_TEST), [0, 1])
    basis_qc = native_basis_change("XI", GATE_NAME)
    return qc.compose(basis_qc)


def build_gpi2_circuit():
    qc = QuantumCircuit(1)
    qc.append(GPI2Gate(0.0), [0])
    basis_qc = native_basis_change("Z", GATE_NAME)  # trivial (Z needs no rotation), kept for consistency
    return qc.compose(basis_qc)


def main():
    print("\n" + "=" * 96)
    print("  task37e_real_calibration_submission.py -- REAL submission, 2 dedicated calibration circuits")
    print("=" * 96)

    zz_circuit = build_zz_circuit()
    gpi2_circuit = build_gpi2_circuit()
    print(f"  ZZ circuit:   {zz_circuit.count_ops()}")
    print(f"  GPi2 circuit: {gpi2_circuit.count_ops()}")

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    circuits = []
    tags = []
    for version, base in [("zz_test", zz_circuit), ("gpi2_test", gpi2_circuit)]:
        for backend_name in BACKENDS:
            c = base.copy()
            c.measure_all()
            circuits.append(c)
            tags.append((version, backend_name))

    all_counts = [None] * len(circuits)
    for backend_name in BACKENDS:
        idxs = [i for i, t in enumerate(tags) if t[1] == backend_name]
        chunk = [circuits[i] for i in idxs]
        job = submit_job(chunk, backend, backend_name, shots=SHOTS)
        counts = get_counts_list(job)
        for i, c in zip(idxs, counts):
            all_counts[i] = c
        print(f"    {backend_name}: {len(chunk)} circuits done")

    with open(RESULTS_PATH.replace(".json", "_raw.json"), "w") as f:
        json.dump({"tags": tags, "counts": all_counts}, f, indent=2)

    print(f"\n  -- RESULTS: measured <XI>/<Z>, implied delta_zz/delta_gpi2, shot-noise floor --")
    results = {}
    for backend_name in BACKENDS:
        idx_zz = tags.index(("zz_test", backend_name))
        idx_gpi2 = tags.index(("gpi2_test", backend_name))
        m_xi = expectation_from_counts(all_counts[idx_zz], "XI")
        m_z = expectation_from_counts(all_counts[idx_gpi2], "Z")
        delta_zz_hat = m_xi / SLOPE_ZZ
        delta_gpi2_hat = m_z / SLOPE_GPI2
        se_delta_zz = SHOT_NOISE_SE_XI / abs(SLOPE_ZZ)
        se_delta_gpi2 = SHOT_NOISE_SE_XI / abs(SLOPE_GPI2)
        sig_zz = abs(delta_zz_hat) > 2 * se_delta_zz
        sig_gpi2 = abs(delta_gpi2_hat) > 2 * se_delta_gpi2
        print(f"\n  === {backend_name} ===")
        print(f"    <XI> = {m_xi:+.6f}  ->  delta_zz_hat   = {delta_zz_hat:+.6f}  (SE~{se_delta_zz:.2e})  "
              f"{'SIGNIFICANT (>2 SE)' if sig_zz else 'within shot-noise floor'}")
        print(f"    <Z>  = {m_z:+.6f}  ->  delta_gpi2_hat = {delta_gpi2_hat:+.6f}  (SE~{se_delta_gpi2:.2e})  "
              f"{'SIGNIFICANT (>2 SE)' if sig_gpi2 else 'within shot-noise floor'}")
        results[backend_name] = {
            "m_xi": m_xi, "m_z": m_z, "delta_zz_hat": delta_zz_hat, "delta_gpi2_hat": delta_gpi2_hat,
            "se_delta_zz": se_delta_zz, "se_delta_gpi2": se_delta_gpi2,
            "significant_zz": bool(sig_zz), "significant_gpi2": bool(sig_gpi2),
        }

    ideal_ok = abs(results["ideal"]["m_xi"]) < 3 * SHOT_NOISE_SE_XI and abs(results["ideal"]["m_z"]) < 3 * SHOT_NOISE_SE_XI
    print(f"\n  -- SANITY CHECK: ideal backend should read ~0 for both (no real hardware noise) -- "
          f"{'PASS' if ideal_ok else 'FAIL -- something is wrong with the circuit/pipeline, do not trust aria-1/forte-1 below'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"ideal_ok": bool(ideal_ok), "results": results,
                    "slope_zz": SLOPE_ZZ, "slope_gpi2": SLOPE_GPI2, "shots": SHOTS}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
