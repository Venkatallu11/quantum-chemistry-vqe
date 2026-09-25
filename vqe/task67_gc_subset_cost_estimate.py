#!/usr/bin/env python3
"""
task67_gc_subset_cost_estimate.py -- iteration 67. REAL, read-only cost
estimates (GET /jobs/estimate, official IonQ pricing endpoint, zero spend)
for candidate real-hardware subsets of the 84-circuit GC design, needed
to answer IonQ's Sep 24 request for a concrete plan that fits the real
remaining $2,352.16 balance (the original 84-circuit/1,100-shot estimate
was ~$3,830, over budget).

CANDIDATE SUBSET, chosen for a stated reason, not arbitrarily: keep ALL
21 measurement slots (dropping slots removes equations from the joint
Schmidt-frame fit's shared 15-parameter system -- the exact structure
task66 just showed rescues the GC design from its own real shot-noise
dilution problem), but only 2 of the 4 GC groups per slot: the TRIVIAL
group (0 extra gates, already the QWC baseline) and the HARDEST group
(most diagonalizer gates, already real-hardware spot-checked for
physical correctness in task48/task49) -- so the run still tests the
real, previously-untested part (does the harder native diagonalizer
behave on real hardware as the free-simulator noise model predicts?)
while cutting circuit count from 84 to 42 (21 slots x 2 groups).
Verified separately (task36's own residual code) that a joint fit with
partial per-slot label coverage stays hugely overdetermined (dof in the
hundreds vs 15 free parameters) -- not a fit-breaking cut.

Run:
    PYTHONHASHSEED=0 python vqe/task67_gc_subset_cost_estimate.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from general_commuting_measurements import build_general_commuting_measurement_plan
from native_stateprep import to_native
from ionq_backend import connect_provider

K = 6
GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
SHOT_OPTIONS = [500, 1100, 2000]


def main():
    print("\n" + "=" * 96)
    print("  task67_gc_subset_cost_estimate.py -- REAL cost estimates, read-only, zero spend")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    ancilla_native = ancilla_cnots_native(GATE_NAME)

    trivial_idx = [i for i, d in enumerate(diagonalizers) if d.two_qubit_count == 0 and d.one_qubit_count == 0][0]
    hard_idx = max(range(len(diagonalizers)), key=lambda i: diagonalizers[i].two_qubit_count)
    print(f"  trivial group index={trivial_idx} ({len(groups[trivial_idx])} labels), "
          f"hard group index={hard_idx} ({len(groups[hard_idx])} labels)")

    diag_natives = {i: to_native(diagonalizers[i].to_circuit(), GATE_NAME) for i in (trivial_idx, hard_idx)}

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")

    per_circuit_gate_counts = {}
    for name in kept:
        register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        for gi in (trivial_idx, hard_idx):
            full5 = with_ancilla_parity_native(register_native, ancilla_native)
            qc = full5.compose(diag_natives[gi], qubits=[0, 1, 2, 3])
            counts = dict(qc.count_ops())
            n1q = counts.get("gpi", 0) + counts.get("gpi2", 0)
            n2q = counts.get(GATE_NAME, 0)
            per_circuit_gate_counts[(name, gi)] = (n1q, n2q)

    print(f"\n  {len(per_circuit_gate_counts)} circuits (21 slots x 2 groups) -- querying real IonQ estimate endpoint")
    for shots in SHOT_OPTIONS:
        total = 0.0
        unit = None
        n_err = 0
        for (name, gi), (n1q, n2q) in per_circuit_gate_counts.items():
            try:
                est = backend.client.estimate_job(backend=BACKEND_NAME, oneq_gates=n1q, twoq_gates=n2q,
                                                    qubits=5, shots=shots)
                total += float(est.cost)
                unit = est.cost_unit
            except Exception as e:
                n_err += 1
                print(f"    ERROR estimating ({name}, group {gi}): {repr(e)}")
        print(f"  shots={shots:>5}: REAL total estimated cost for 42-circuit subset (21 slots x trivial+hard) "
              f"= {total:.2f} {unit}  (errors={n_err})")

    print(f"\n  -- for comparison, full 84-circuit (all 4 groups) design at the same shot options --")
    diag_natives_all = [to_native(d.to_circuit(), GATE_NAME) for d in diagonalizers]
    per_circuit_all = {}
    for name in kept:
        register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        for gi in range(len(diagonalizers)):
            full5 = with_ancilla_parity_native(register_native, ancilla_native)
            qc = full5.compose(diag_natives_all[gi], qubits=[0, 1, 2, 3])
            counts = dict(qc.count_ops())
            n1q = counts.get("gpi", 0) + counts.get("gpi2", 0)
            n2q = counts.get(GATE_NAME, 0)
            per_circuit_all[(name, gi)] = (n1q, n2q)

    for shots in SHOT_OPTIONS:
        total = 0.0
        unit = None
        n_err = 0
        for (name, gi), (n1q, n2q) in per_circuit_all.items():
            try:
                est = backend.client.estimate_job(backend=BACKEND_NAME, oneq_gates=n1q, twoq_gates=n2q,
                                                    qubits=5, shots=shots)
                total += float(est.cost)
                unit = est.cost_unit
            except Exception as e:
                n_err += 1
        print(f"  shots={shots:>5}: REAL total estimated cost for full 84-circuit design "
              f"= {total:.2f} {unit}  (errors={n_err})")


if __name__ == "__main__":
    main()
