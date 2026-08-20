#!/usr/bin/env python3
"""
task32a_circuit_worker.py -- single-circuit submit+retrieve worker for
task32a_gpi2_calibration.py, run as a SEPARATE OS PROCESS specifically so
the parent can enforce a hard wall-clock timeout (subprocess.run(...,
timeout=X) reliably kills a hung child). This exists because the parent's
own in-process retry loop got stuck for 78+ minutes on one circuit with NO
progress and NO printed retry message -- consistent with a network call
inside qiskit-ionq that has no internal read timeout, hanging indefinitely
rather than raising an exception my own try/except could ever catch.

Usage:
    python vqe/task32a_circuit_worker.py <phi> <N> <model> <shots> <output_json_path>
Writes {"counts": {...}} to output_json_path on success. Any exception
propagates as a normal nonzero exit code (visible to the parent via
subprocess.run's returncode / CalledProcessError).
"""
import sys
import os
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list
from task32a_gpi2_calibration import build_identity_circuit


def main():
    phi = float(sys.argv[1])
    N = int(sys.argv[2])
    model = sys.argv[3]
    shots = int(sys.argv[4])
    out_path = sys.argv[5]

    provider = connect_provider()
    backend = get_native_simulator(provider)
    qc = build_identity_circuit(phi, N)
    qc.measure_all()
    job = submit_job([qc], backend, model, shots=shots)
    counts = get_counts_list(job)[0]
    with open(out_path, "w") as f:
        json.dump({"counts": counts}, f)


if __name__ == "__main__":
    main()
