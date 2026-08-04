#!/usr/bin/env python3
"""
ionq_resource_estimate.py — real 2-qubit gate counts for a full H4 EF+ZNE
run using GENUINE native-gate folding, vs. the old (cancelled-fold)
assumption.
============================================================================
ionq_run.py's original abstract-gate ZNE attempt submitted circuits that
LOOKED like fold=1/3/5 (more gates in the circuit object) but, per
ionq_fold_check.py's native-vs-abstract comparison, those folds were being
CANCELLED before execution -- so every "folded" abstract circuit actually
cost the SAME real 2-qubit gate work as fold=1, regardless of what fold
factor was requested. Any resource/cost estimate based on the abstract
circuits' fold factor was therefore an UNDER-estimate of what a genuine
(uncancelled) folded run actually costs.

ionq_fold_check.py --native proved gateset="native" submission is NOT
cancelled (f decays cleanly with fold on both aria-1 and forte-1). This
script computes what a full H4 entanglement-forging run would actually
cost in real 2-qubit gate operations if redone with genuine native
folding, using numbers already measured elsewhere in this repo:
  - circuits per register at K=5: 325 (ionq_tailoring.py fragment A,
    qubit-wise-grouped, real-gauge 2-phase-circuit measurement)
  - native 2-qubit gates per state-prep circuit: 14 (native_stateprep.py,
    verified on real IonQ)
Folding multiplies the state-prep circuit's native 2-qubit gate count by
the fold factor (ionq_run.fold_circuit's convention: fold every 2-qubit
gate in the transpiled circuit, not just a single "entangler").

Run:
    python vqe/ionq_resource_estimate.py
"""
import json
import os

CIRCUITS_PER_REGISTER_K5 = 325   # ionq_tailoring_results.json, fragment A, K=5
NATIVE_2Q_GATES_PER_STATEPREP = 14  # native_stateprep_results.json, verified real+local
FOLD_FACTORS_ZNE = [1, 3, 5]      # the fold set ionq_run.py's original ZNE attempt used


def estimate(circuits_per_register, gates_per_circuit, fold_factors):
    rows = []
    for fold in fold_factors:
        real_2q_gates_per_circuit = gates_per_circuit * fold
        total_2q_gates = circuits_per_register * real_2q_gates_per_circuit
        rows.append({
            "fold": fold,
            "native_2q_gates_per_circuit": real_2q_gates_per_circuit,
            "total_2q_gates_this_fold": total_2q_gates,
        })
    return rows


def main():
    print("\n" + "=" * 70)
    print("  Resource estimate: real 2-qubit gate count, genuine native folding")
    print("=" * 70)
    print(f"\n  Basis numbers (measured elsewhere in this repo, not assumed):")
    print(f"    circuits per register at K=5: {CIRCUITS_PER_REGISTER_K5}")
    print(f"    native 2-qubit gates per state-prep circuit (fold=1): "
          f"{NATIVE_2Q_GATES_PER_STATEPREP}")

    print(f"\n  {'fold':>6} {'2q gates/circuit':>18} {'total 2q gates (1 register)':>30} "
          f"{'OLD (cancelled) estimate':>26}")
    rows = estimate(CIRCUITS_PER_REGISTER_K5, NATIVE_2Q_GATES_PER_STATEPREP, FOLD_FACTORS_ZNE)
    old_total_per_fold = CIRCUITS_PER_REGISTER_K5 * NATIVE_2Q_GATES_PER_STATEPREP  # fold=1 always, cancelled
    for r in rows:
        print(f"  {r['fold']:>6} {r['native_2q_gates_per_circuit']:>18} "
              f"{r['total_2q_gates_this_fold']:>30} {old_total_per_fold:>26}")

    real_total = sum(r["total_2q_gates_this_fold"] for r in rows)
    old_total = old_total_per_fold * len(FOLD_FACTORS_ZNE)
    ratio = real_total / old_total

    print(f"\n  Summed across fold={FOLD_FACTORS_ZNE} (1 register, K=5, one noise model):")
    print(f"    REAL (genuine native folding)  : {real_total:,} two-qubit gate operations")
    print(f"    OLD assumption (folds cancelled): {old_total:,} two-qubit gate operations")
    print(f"    Real cost is {ratio:.2f}x the old estimate")

    print(f"\n  Both registers (alpha measured + beta independently, worst case "
          f"where beta-reuse doesn't hold, matching the overlap-fragment finding):")
    print(f"    REAL: {real_total * 2:,} two-qubit gate operations")

    print(f"\n  This is PER NOISE MODEL. For aria-1 + forte-1 (2 models): "
          f"{real_total * 2 * 2:,} two-qubit gate operations, alpha register only.")

    print(f"\n  HONEST CONTEXT: this is a resource estimate for ionq_simulator "
          "(free, this project's scope). It is NOT a cost estimate for ionq_qpu "
          "(real paid trapped-ion hardware) -- this project has not priced that "
          "out and does not touch ionq_qpu. The point of this estimate is that "
          "the fold count materially changes the REAL gate count once folds "
          "aren't being silently cancelled -- any future decision to run this "
          "on real hardware needs to budget for the REAL total, not the old "
          "(cancelled-fold) one.")

    results = {
        "circuits_per_register_k5": CIRCUITS_PER_REGISTER_K5,
        "native_2q_gates_per_stateprep": NATIVE_2Q_GATES_PER_STATEPREP,
        "fold_factors": FOLD_FACTORS_ZNE,
        "per_fold": rows,
        "real_total_2q_gates_one_register_one_model": real_total,
        "old_cancelled_assumption_total": old_total,
        "real_vs_old_ratio": round(ratio, 3),
        "real_total_both_registers_one_model": real_total * 2,
        "real_total_both_registers_two_models": real_total * 2 * 2,
    }
    out = os.path.join(os.path.dirname(__file__), "ionq_resource_estimate_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
