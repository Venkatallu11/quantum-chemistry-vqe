#!/usr/bin/env python3
"""
task74_real_hardware_cost_check.py -- iteration 74. Verifies, with the
real IonQ cost-estimate endpoint (GET /jobs/estimate, free, read-only,
zero spend), whether increasing shots/circuit on the approved 54-circuit
real-hardware plan (task71/72) fits the approved $2,352.16 budget on the
real qpu.forte-enterprise-1 device -- motivated by task72's finding that
a single real 1,100-shot draw is NOT reliably inside chemical accuracy,
and task73's finding that pooling 4 trials (4,400 effective shots) is.

REAL BUG FOUND AND WORKED AROUND: the installed qiskit-ionq SDK's
JobEstimate.cost reads the response field "estimated_cost", but the
live API actually returns "estimated_total_cost" -- so .cost is always
None, which looked exactly like a rate limit on the first few attempts
(it isn't one). This reads the raw response directly instead.

METHOD: only 2 real probe calls (paced 5s apart) -- the cheapest and
priciest circuit in the panel, at the approved 1,100 shots -- solve the
resulting 2x2 linear system for the real verified per-gate rates, then
compute the exact real cost for the full 54-circuit design at several
shot levels analytically. Also surfaces the real current queue time for
forte-enterprise-1 (found incidentally: ~14.5 days at the time this ran).

Run:
    PYTHONHASHSEED=0 python vqe/task74_real_hardware_cost_check.py
"""
import os
import sys
import time
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from native_stateprep import to_native
from ionq_backend import connect_provider
import task72_h4_full_experiment as t72

GATE_NAME = "zz"
BACKEND_NAME = "qpu.forte-enterprise-1"
PROBE_SHOTS = 1100
BUDGET = 2352.16
SHOT_LEVELS = (1100, 1650, 2200, 2750, 3300, 3850, 4400)
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task74_cost_rate_check_results.json")


def build_gate_count_table(ctx):
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    diag_natives_all = {i: to_native(d.to_circuit(), GATE_NAME) for i, d in enumerate(ctx["diagonalizers"])}
    table = {}
    for name in ctx["kept"]:
        register_native = optimized_native_circuit(ctx["fixed_solutions"][name]["angles"], GATE_NAME)
        for gi in ctx["measured_groups_map"][name]:
            full5 = with_ancilla_parity_native(register_native, ancilla_native)
            qc = full5.compose(diag_natives_all[gi], qubits=[0, 1, 2, 3])
            counts = dict(qc.count_ops())
            n1q = counts.get("gpi", 0) + counts.get("gpi2", 0)
            n2q = counts.get(GATE_NAME, 0)
            table[(name, gi)] = (n1q, n2q)
    return table


def main():
    print("\n" + "=" * 96)
    print("  task74_real_hardware_cost_check.py -- real, paced, read-only cost verification")
    print("=" * 96)

    ctx = t72.setup()
    table = build_gate_count_table(ctx)

    hard_idx = max(range(len(ctx["diagonalizers"])), key=lambda i: ctx["diagonalizers"][i].two_qubit_count)
    reused = {("u_0", hard_idx), ("u_1", hard_idx), ("u_2", hard_idx)}
    print(f"  hard group index = {hard_idx}; reused-free circuits (task49's real data): {reused}")

    cheapest_key = min(table, key=lambda k: table[k])
    priciest_key = max(table, key=lambda k: table[k][0] + table[k][1])
    print(f"  cheapest circuit: {cheapest_key} -> {table[cheapest_key]}")
    print(f"  priciest circuit: {priciest_key} -> {table[priciest_key]}")

    provider = connect_provider()
    backend = provider.get_backend(BACKEND_NAME, gateset="native")

    probes = {}
    for key in (cheapest_key, priciest_key):
        n1q, n2q = table[key]
        est = backend.client.estimate_job(backend=BACKEND_NAME, oneq_gates=n1q, twoq_gates=n2q, qubits=5,
                                            shots=PROBE_SHOTS)
        raw = est.to_dict()
        cost = raw.get("estimated_total_cost")
        queue_s = raw.get("current_predicted_queue_time")
        queue_str = f"{queue_s / 86400:.2f} days" if queue_s else "n/a"
        print(f"  probe {key}: n1q={n1q} n2q={n2q} shots={PROBE_SHOTS} -> "
              f"cost=${cost} {raw.get('estimated_unit')}  queue={queue_str}")
        probes[key] = (n1q, n2q, float(cost), queue_s)
        time.sleep(5)

    (n1a, n2a, ca, _), (n1b, n2b, cb, queue_s) = probes[cheapest_key], probes[priciest_key]
    A = np.array([[n1a, n2a], [n1b, n2b]], dtype=float)
    b = np.array([ca / PROBE_SHOTS, cb / PROBE_SHOTS], dtype=float)
    r1, r2 = np.linalg.solve(A, b)
    print(f"\n  REAL verified rates: 1q=${r1:.8f}/shot  2q=${r2:.8f}/shot")

    def total_cost(shots, skip_reused_upto_shots=2000):
        total = 0.0
        for key, (n1q, n2q) in table.items():
            if key in reused and shots <= skip_reused_upto_shots:
                continue  # free reuse from task49's real, already-paid-for data
            total += shots * (n1q * r1 + n2q * r2)
        return total

    totals = {}
    print(f"\n  -- real cost for the full 54-circuit panel on {BACKEND_NAME}, budget=${BUDGET} --")
    for shots in SHOT_LEVELS:
        tc = total_cost(shots)
        totals[shots] = tc
        verdict = "FITS" if tc <= BUDGET else f"OVER by ${tc - BUDGET:.2f}"
        print(f"  shots={shots:>5} ({shots / 1100:.2f}x approved): real cost = ${tc:.2f}  [{verdict}]")

    results = {"rate_1q": float(r1), "rate_2q": float(r2), "budget": BUDGET,
               "queue_time_days": (queue_s / 86400 if queue_s else None),
               "totals": {str(s): v for s, v in totals.items()}}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
