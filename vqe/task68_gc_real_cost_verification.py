#!/usr/bin/env python3
"""
task68_gc_real_cost_verification.py -- iteration 68. Verifies (not
assumes) everything needed to design a real, defensible hardware-budget
plan for the GC-based H4 experiment IonQ has approved spending
$2,352.16 on:

1. Reproduces the confirmed-live IonQ rate card from task1_shot_budget.py
   (1q=$0.000164, 2q=$0.001121, job_floor=$25.7899, floor applies once
   per bundled multi-circuit job) and cross-checks it against task49's
   real stored gate counts + previously-cited real cost ($288.85 for its
   3-circuit GC3 batch) -- if the formula doesn't reproduce that real
   number, STOP, don't trust anything downstream.

2. Checks whether task49's real, already-completed GC3 circuits (u_0,
   u_1, u_2) can legitimately be reused: rebuilds those exact 3 circuits
   with TODAY's code and diffs gate counts against the real stored
   values. If they match, the circuit-construction code hasn't drifted
   since Sept 10 and reuse is valid. If not, say so plainly.

3. Computes the REAL, exact, per-(slot,group) gate counts for all 84
   GC circuits and the real formula-cost of the full 84-circuit design
   at 1,100 shots -- cross-checked against Vadim's own cited ~$3,830
   estimate, as an end-to-end sanity check of this whole methodology.

4. Saves the full real per-circuit cost table so a circuit-selection
   plan can be built from exact numbers, not per-group averages.

Run:
    PYTHONHASHSEED=0 python vqe/task68_gc_real_cost_verification.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from general_commuting_measurements import build_general_commuting_measurement_plan
from native_stateprep import to_native

K = 6
GATE_NAME = "zz"
RATE_1Q = 0.000164
RATE_2Q = 0.001121
JOB_FLOOR = 25.7899
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task68_gc_real_cost_verification_results.json")


def circuit_cost(n1q, n2q, shots):
    """Single circuit's own gate-execution cost (the floor is compared
    against the SUM across a bundled job, applied separately below)."""
    return shots * (n1q * RATE_1Q + n2q * RATE_2Q)


def main():
    print("\n" + "=" * 96)
    print("  task68_gc_real_cost_verification.py -- verify, don't assume, the real cost model")
    print("=" * 96)

    # ---- 1. reproduce task49's real bill from its real stored gate counts ----
    task49_counts = {"u_0": (171, 16), "u_1": (180, 13), "u_2": (195, 20)}
    task49_total = sum(circuit_cost(n1q, n2q, 2000) for n1q, n2q in task49_counts.values())
    print(f"\n  -- STEP 1: reproduce task49's real cited cost (~$288.85) from its real stored gate counts --")
    for name, (n1q, n2q) in task49_counts.items():
        print(f"    {name}: {n1q}x1q + {n2q}x2q @ 2000 shots -> ${circuit_cost(n1q, n2q, 2000):.2f}")
    print(f"    formula total: ${task49_total:.2f}  vs previously-cited real ${288.85}  "
          f"{'MATCH (within rounding)' if abs(task49_total - 288.85) < 1.0 else 'MISMATCH -- STOP, formula not trustworthy'}")
    if abs(task49_total - 288.85) >= 1.0:
        return

    # ---- 2. rebuild task49's 3 circuits with TODAY's code, diff gate counts ----
    print(f"\n  -- STEP 2: does TODAY's code reproduce task49's real stored gate counts (reuse validity check) --")
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    ancilla_native = ancilla_cnots_native(GATE_NAME)

    hard_idx = max(range(len(diagonalizers)), key=lambda i: diagonalizers[i].two_qubit_count)
    hard_diag_native = to_native(diagonalizers[hard_idx].to_circuit(), GATE_NAME)
    print(f"    hard group index today = {hard_idx} (task49 used index 3 -- "
          f"{'CONSISTENT' if hard_idx == 3 else 'DIFFERENT -- investigate before reuse'})")

    reuse_valid = True
    for name in ["u_0", "u_1", "u_2"]:
        register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        full5 = with_ancilla_parity_native(register_native, ancilla_native)
        qc = full5.compose(hard_diag_native, qubits=[0, 1, 2, 3])
        counts = dict(qc.count_ops())
        n1q_today = counts.get("gpi", 0) + counts.get("gpi2", 0)
        n2q_today = counts.get(GATE_NAME, 0)
        n1q_real, n2q_real = task49_counts[name]
        match = (n1q_today == n1q_real) and (n2q_today == n2q_real)
        reuse_valid = reuse_valid and match
        print(f"    {name}: today={n1q_today}x1q/{n2q_today}x2q  real-stored={n1q_real}x1q/{n2q_real}x2q  "
              f"{'MATCH -- safe to reuse' if match else 'MISMATCH -- code drifted, do NOT claim reuse without re-verifying physics'}")

    # ---- 3. real gate counts + cost for all 84 (slot, group) circuits ----
    print(f"\n  -- STEP 3: real gate counts + formula cost, all 84 circuits, cross-check vs Vadim's cited ~$3,830 --")
    diag_natives = [to_native(d.to_circuit(), GATE_NAME) for d in diagonalizers]
    table = {}
    for name in kept:
        register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        for gi, dn in enumerate(diag_natives):
            full5 = with_ancilla_parity_native(register_native, ancilla_native)
            qc = full5.compose(dn, qubits=[0, 1, 2, 3])
            counts = dict(qc.count_ops())
            n1q = counts.get("gpi", 0) + counts.get("gpi2", 0)
            n2q = counts.get(GATE_NAME, 0)
            table[f"{name}|g{gi}"] = {"n1q": n1q, "n2q": n2q, "n_labels": len(groups[gi])}

    for shots in [1100, 2000]:
        total = sum(circuit_cost(v["n1q"], v["n2q"], shots) for v in table.values())
        total = max(total, JOB_FLOOR)  # floor applies once, to the whole bundled job
        print(f"    shots={shots}: real formula total for all 84 circuits (ONE bundled job) = ${total:,.2f}"
              + (f"   vs Vadim's cited ~$3,830  diff={total-3830:+.2f}" if shots == 1100 else ""))

    with open(RESULTS_PATH, "w") as f:
        json.dump({"task49_reproduction": {"formula_total": task49_total, "cited_real": 288.85},
                    "reuse_valid": reuse_valid, "hard_group_index_today": hard_idx,
                    "per_circuit_table": table, "rate_card": {"1q": RATE_1Q, "2q": RATE_2Q, "floor": JOB_FLOOR}},
                  f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
