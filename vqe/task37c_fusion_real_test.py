#!/usr/bin/env python3
"""
task37c_fusion_real_test.py -- iteration 37, circuit-level branch, real
test. Does the verified-correct gate_fusion.py reduction (23.7% fewer
GPi2/GPi gates, exact unitary equivalence) actually show up as LOWER
real error on real noisy backends -- not just fewer gates on paper?

METHOD: for 2 representative (slot, group) cases, build BOTH the
ORIGINAL (base.compose(basis_qc), this project's standard construction)
and FUSED (gate_fusion.fuse_single_qubit_runs) circuits, submit BOTH,
RAW (no PEC -- this isolates the circuit-level effect itself, not
mixed with any correction method), to THREE backends: `ideal` (sanity
check -- both versions MUST match exactly, zero noise), `aria-1`, and
`forte-1`. Real shots (100k), matching this project's standard.

Run:
    python vqe/task37c_fusion_real_test.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit
from native_stateprep import native_target
from gate_fusion import fuse_single_qubit_runs
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, expectation_from_counts
from qiskit.quantum_info import Statevector, Pauli

K = 6
GATE_NAME = "zz"
SHOTS = 100_000
BACKENDS = ["ideal", "aria-1", "forte-1"]
TEST_CASES = [
    {"slot": "(u0+u1)", "group": ["XYYX", "IYYI"], "tag": "IYYI (established validation case)"},
    {"slot": "(u3+u5)", "group": ["XZXZ", "XZXI", "IZIZ", "IIIZ"], "tag": "(u3+u5) (Task 33D-flagged group)"},
]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37c_fusion_real_test_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task37c_fusion_real_test.py -- original vs fused circuit, real submission, ideal/aria-1/forte-1")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    circuits, tags = [], []
    exact_by_case = {}
    for case in TEST_CASES:
        slot, group = case["slot"], case["group"]
        base = optimized_native_circuit(fixed_solutions[slot]["angles"], GATE_NAME)
        basis_qc = native_basis_change(effrag_mod.combined_basis_label(group), GATE_NAME)
        full_orig = base.compose(basis_qc)
        tgt = native_target(full_orig.num_qubits, GATE_NAME)
        full_fused = fuse_single_qubit_runs(full_orig, tgt)

        exact_by_case[slot] = {}
        for l in group:
            exact_by_case[slot][l] = float(np.real(Statevector.from_instruction(full_orig)
                                                     .expectation_value(Pauli(l))))

        for version, circ in [("original", full_orig), ("fused", full_fused)]:
            for backend_name in BACKENDS:
                c = circ.copy()
                c.measure_all()
                circuits.append(c)
                tags.append((slot, group, version, backend_name))

    print(f"  built {len(circuits)} circuits (2 cases x 2 versions x {len(BACKENDS)} backends)")

    # submit each backend's circuits as its own job (ionq_simulator noise_model is a job-level option)
    all_counts = [None] * len(circuits)
    for backend_name in BACKENDS:
        idxs = [i for i, t in enumerate(tags) if t[3] == backend_name]
        chunk = [circuits[i] for i in idxs]
        job = submit_job(chunk, backend, backend_name, shots=SHOTS)
        counts = get_counts_list(job)
        for i, c in zip(idxs, counts):
            all_counts[i] = c
        print(f"    {backend_name}: {len(chunk)} circuits done")

    with open(RESULTS_PATH.replace(".json", "_raw.json"), "w") as f:
        json.dump({"tags": tags, "counts": all_counts}, f, indent=2)

    print(f"\n  -- RESULTS: raw error vs exact, original vs fused, per backend --")
    results = {}
    for case in TEST_CASES:
        slot, group = case["slot"], case["group"]
        print(f"\n  === {case['tag']} ===")
        for backend_name in BACKENDS:
            row = {}
            for version in ["original", "fused"]:
                idx = tags.index((slot, group, version, backend_name))
                counts = all_counts[idx]
                errs = []
                for l in group:
                    m = expectation_from_counts(counts, l)
                    errs.append(abs(m - exact_by_case[slot][l]))
                row[version] = float(np.mean(errs))
            delta = row["original"] - row["fused"]
            print(f"    {backend_name:<8} original_err={row['original']:.4f}  fused_err={row['fused']:.4f}  "
                  f"(fused {'BETTER' if delta>0 else 'WORSE'} by {abs(delta):.4f})")
            results.setdefault(case["tag"], {})[backend_name] = row

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
