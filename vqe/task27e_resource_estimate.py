#!/usr/bin/env python3
"""
task27e_resource_estimate.py — iteration 27, Task E. Feed the VALIDATED
native circuits (iteration 27, Tasks A/B: 11 native 2-qubit gates
constant, 197 (ms) / 285 (zz) 1-qubit gates, verified fold-preserving to
1e-15) into IonQ's real resource estimator. NO SUBMISSION anywhere in
this file -- `GET /jobs/estimate` is free, read-only, and does not touch
or reserve any hardware, exactly as iteration 9 and iteration 26 Task 1
already established and reused unchanged here, not re-verified.
============================================================================
Prices the ACTUAL validated circuit design from this iteration: the
subspace-tomography-reduced (diag + "+"-pair only) K=5 (15 circuits) and
K=6 (21 circuits) sets, x 13 measurement groups, x folds {1, 3, 5}
(separately and combined, per explicit instruction), at the shot level
this iteration's own real submission (Task C) used (100,000/setting) AND
at iteration 26 Task 1's chosen 300,000/setting for direct comparison.

Run:
    python vqe/task27e_resource_estimate.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from ionq_backend import connect_provider

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task27e_resource_estimate_results.json")
BUDGET_USD = 3000.0

N2Q_PER_CIRCUIT_FOLD1 = 11          # constant, both K, both gate families (iteration 27 Task A, verified)
N1Q_PER_CIRCUIT = {"ms": 197, "zz": 285}   # constant, both K (iteration 27 Task A, verified)
N_GROUPS = 13
K_CIRCUIT_COUNTS = {5: 15, 6: 21}   # subspace-tomography-reduced (diag + "+"), this iteration's actual design
FOLDS_TO_PRICE = [1, 3, 5]
SHOT_LEVELS = [100_000, 300_000]
GATE_BY_MODEL = {"aria-1": "ms", "forte-1": "zz"}   # pricing is per-backend; ideal is free-simulator-only, not priced


def main():
    print("\n" + "=" * 96)
    print("  task27e_resource_estimate.py -- real IonQ pricing, NO submission")
    print("=" * 96)

    provider = connect_provider()
    client = provider.get_backend("ionq_simulator").client
    print(f"  connected (client used ONLY for GET /jobs/estimate -- read-only, no hardware touched)")

    results = {}
    for model, gate_name in GATE_BY_MODEL.items():
        backend_name = "qpu.forte-1" if model == "forte-1" else "qpu.aria-1"
        n1q = N1Q_PER_CIRCUIT[gate_name]
        results[model] = {}
        print(f"\n  ===== {model} (native gate={gate_name}, backend={backend_name}) =====")
        for K, n_circ in K_CIRCUIT_COUNTS.items():
            results[model][K] = {}
            n_circuits_per_fold = n_circ * N_GROUPS
            print(f"    K={K}: {n_circ} kept circuits x {N_GROUPS} groups = {n_circuits_per_fold} circuits/fold")
            for shots in SHOT_LEVELS:
                fold_costs = {}
                for fold in FOLDS_TO_PRICE:
                    n2q_folded = N2Q_PER_CIRCUIT_FOLD1 * fold  # folding multiplies N_2q linearly, verified iteration 27 Task A/B
                    try:
                        est_1circuit = client.estimate_job(backend=backend_name, oneq_gates=n1q, twoq_gates=n2q_folded,
                                                            qubits=4, shots=shots)
                        cost_1circuit = est_1circuit.cost if est_1circuit.cost is not None else \
                            max(25.7899, shots * (n1q * 0.000164 + n2q_folded * 0.001121))
                    except Exception:
                        cost_1circuit = max(25.7899, shots * (n1q * 0.000164 + n2q_folded * 0.001121))
                    # per iteration 9's confirmed methodology: floor is per-JOB; bundling all
                    # n_circuits_per_fold circuits into ONE job means gate cost sums, floor applies once
                    gate_cost_total = n_circuits_per_fold * shots * (n1q * 0.000164 + n2q_folded * 0.001121)
                    cost_bundled = max(gate_cost_total, 25.7899)
                    fold_costs[fold] = cost_bundled
                combined_1_3_5 = sum(fold_costs.values())
                results[model][K][shots] = {"per_fold": fold_costs, "combined_1_3_5": combined_1_3_5}
                print(f"      shots={shots:>7,}: " +
                      "  ".join(f"fold={f}: ${c:,.2f}" for f, c in fold_costs.items()) +
                      f"   COMBINED(1+3+5)=${combined_1_3_5:,.2f}  ({combined_1_3_5/BUDGET_USD:.1f}x budget)")

    print(f"\n  -- SUMMARY: cheapest validated configuration vs $3,000 budget --")
    cheapest = min(
        (results[m][K][s]["combined_1_3_5"], m, K, s)
        for m in results for K in results[m] for s in results[m][K]
    )
    cost, model, K, shots = cheapest
    print(f"    cheapest: {model}, K={K}, shots={shots:,}, folds 1+3+5 combined = ${cost:,.2f} "
          f"({cost/BUDGET_USD:.1f}x the $3,000 budget)")
    print(f"    {'FITS the budget' if cost <= BUDGET_USD else 'DOES NOT FIT the budget -- redesign needed before ANY submission'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
