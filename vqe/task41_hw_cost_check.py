#!/usr/bin/env python3
"""
task41_hw_cost_check.py -- iteration 41. REAL hardware cost estimation
(GET /jobs/estimate, IonQ's own official pricing endpoint) for candidate
first-hardware-test circuit sizes, on BOTH qpu.forte-1 and
qpu.forte-enterprise-1 -- READ ONLY, submits nothing, spends nothing.
Run this and inspect its output BEFORE any real submission is even
considered.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from ionq_backend import connect_provider

K = 6
GATE_NAME = "zz"


def main():
    print("\n" + "=" * 96)
    print("  task41_hw_cost_check.py -- REAL cost estimates, read-only, zero spend")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)

    register = optimized_native_circuit(fixed_solutions["u_0"]["angles"], GATE_NAME)
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    full5 = with_ancilla_parity_native(register, ancilla_native)
    group = ["XYYX", "IYYI"]
    combined = effrag_mod.combined_basis_label(group)
    basis_qc = native_basis_change(combined, GATE_NAME)
    full5_meas = full5.compose(basis_qc, qubits=[0, 1, 2, 3])
    counts5 = full5_meas.count_ops()
    n1q_5 = counts5.get("gpi", 0) + counts5.get("gpi2", 0)
    n2q_5 = counts5.get("zz", 0)
    print(f"\n  ONE ancilla-augmented H4 measurement circuit (5 qubits, u_0/(XYYX,IYYI) group):")
    print(f"    gate counts: {dict(counts5)}  -> 1q={n1q_5} 2q={n2q_5} qubits=5")

    reg_only = optimized_native_circuit(fixed_solutions["u_0"]["angles"], GATE_NAME)
    reg_meas = reg_only.compose(basis_qc)
    counts4 = reg_meas.count_ops()
    n1q_4 = counts4.get("gpi", 0) + counts4.get("gpi2", 0)
    n2q_4 = counts4.get("zz", 0)
    print(f"\n  ONE plain (no ancilla) H4 measurement circuit (4 qubits, same slot/group):")
    print(f"    gate counts: {dict(counts4)}  -> 1q={n1q_4} 2q={n2q_4} qubits=4")

    # the 2 tiny Task 37E calibration circuits
    print(f"\n  Task 37E's own 2 tiny calibration circuits (already real-verified on the free simulator):")
    print(f"    ZZ circuit:   gpi2=4 zz=1 gpi=1  -> 1q=5 2q=1 qubits=2")
    print(f"    GPi2 circuit: gpi2=1             -> 1q=1 2q=0 qubits=1")

    provider = connect_provider()
    b_forte = provider.get_backend("qpu.forte-1")
    b_forte_ent = provider.get_backend("qpu.forte-enterprise-1")

    scenarios = [
        ("Task37E ZZ-calib circuit",    5, 1, 2),
        ("Task37E GPi2-calib circuit",  1, 0, 1),
        ("ONE ancilla H4 circuit",      n1q_5, n2q_5, 5),
        ("ONE plain H4 circuit",        n1q_4, n2q_4, 4),
    ]
    shot_options = [100, 1000, 10000]

    print(f"\n  -- REAL COST ESTIMATES (GET /jobs/estimate, official IonQ pricing, per single circuit) --")
    for backend, bname in [(b_forte, "qpu.forte-1"), (b_forte_ent, "qpu.forte-enterprise-1")]:
        print(f"\n  === {bname} ===")
        for label, n1q, n2q, nq in scenarios:
            for shots in shot_options:
                try:
                    est = backend.client.estimate_job(backend=bname, oneq_gates=n1q, twoq_gates=n2q,
                                                        qubits=nq, shots=shots)
                    print(f"    {label:<28} shots={shots:>6}  ->  cost={est.cost} {est.cost_unit}   "
                          f"exec_time={est.exec_time}s  queue_time={est.queue_time}")
                except Exception as e:
                    print(f"    {label:<28} shots={shots:>6}  ->  ERROR: {repr(e)}")


if __name__ == "__main__":
    main()
