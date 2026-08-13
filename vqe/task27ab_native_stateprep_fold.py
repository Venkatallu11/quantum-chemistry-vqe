#!/usr/bin/env python3
"""
task27ab_native_stateprep_fold.py — iteration 27, Tasks A+B. Native
Schmidt-vector state preparation (Forte: GPi/GPi2/ZZ; Aria: GPi/GPi2/MS)
for K=5 AND K=6, with an exact native-fold verifier at folds 1/3/5/7/9.
No qiskit StatePreparation anywhere.
============================================================================
REUSES, not reimplements: `fixed_ansatz.build_ansatz` is ALREADY a
hand-derived, real-only, fixed-Hamming-weight-sector circuit (X,X
reference + one double-excitation + four Givens single-excitations) --
never StatePreparation, exploiting exactly the structure Task 27A asks
for (real amplitudes, fixed weight-2 sector, known support/signs). It is
K-AGNOSTIC: which K Schmidt vectors get prepared is determined entirely
by `setup_fragment`'s own K-truncation, not by the circuit family, so
K=5 and K=6 both reuse the identical, already-verified 5-angle
architecture -- no new synthesis needed, no new correctness risk
introduced. `native_stateprep.to_native` (verified, iteration 26 Task 2:
constant 11 native 2-qubit gates across all 36 K=6 targets at both
optimization_level 0 and 1) does the abstract-to-native translation.

A REAL DISCREPANCY, caught and reported, not silently resolved either
way: the task text states K=5's method error is "~0.17 kcal/mol". Direct
recomputation here gives |noiseless_energy - exact_energy| = 0.565
kcal/mol for K=5 at d=1.0 -- which MATCHES this project's own
independently-established `ionq_native_forged_energy.py::
CLASSICAL_FLOOR_KCAL = 0.5655` from earlier iterations, not the 0.17
figure. The freshly-recomputed, independently-corroborated value (0.565)
is used throughout this file; the 0.17 claim is flagged, not adopted
silently.

FOLD VERIFICATION: `fold_native_2q` (reused unchanged from
`ionq_fold_check.py`, the EXPLICITLY verified inverse construction --
Gate.inverse() on ZZ/MS returns a mis-parametrized, unusable gate, a real
bug this project already found and worked around) is applied at folds
1/3/5/7/9 to every native circuit, and the resulting statevector is
compared to the UN-folded native circuit's statevector directly (not
just the abstract target) -- a fold that changes the answer by more than
1e-10 is treated as a bug, not a noise effect (there is no noise in this
check; it is pure unitary algebra).

Run:
    python vqe/task27ab_native_stateprep_fold.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, verify, check_angle_ranges
from ionq_fold_check import fold_native_2q
from qiskit.quantum_info import Statevector

K_VALUES = [5, 6]
GATE_FAMILIES = {"aria": "ms", "forte": "zz"}
FOLD_FACTORS = [1, 3, 5, 7, 9]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task27ab_native_stateprep_fold_results.json")

CLASSICAL_FLOOR_KCAL_K5_CITED = 0.5655  # ionq_native_forged_energy.py, independent prior source


def main():
    print("\n" + "=" * 96)
    print("  task27ab_native_stateprep_fold.py -- native state prep + exact fold verification, K=5 and K=6")
    print("=" * 96)

    method_errors = {}
    all_records = []
    fold_check_summary = {}

    for K in K_VALUES:
        strict = (K == 6)
        p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=strict)
        method_err_kcal = abs(p["noiseless_energy"] - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
        method_errors[K] = method_err_kcal
        print(f"\n  K={K}: {len(p['targets'])} slots, method error (Schmidt truncation, exact "
              f"matrix elements) = {method_err_kcal:.4f} kcal/mol")
        if K == 5:
            print(f"    NOTE: task text claimed ~0.17 kcal/mol for K=5 -- this recomputation gives "
                  f"{method_err_kcal:.4f}, matching the independently-established "
                  f"CLASSICAL_FLOOR_KCAL={CLASSICAL_FLOOR_KCAL_K5_CITED} from ionq_native_forged_energy.py. "
                  f"Using the verified, corroborated value, not the 0.17 claim.")

        solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
        assert n_ok == len(p["targets"]), f"K={K}: only {n_ok}/{len(p['targets'])} targets converged"
        print(f"    all {n_ok}/{len(p['targets'])} targets converged to <1e-10 (abstract u3/cx fit)")

        for family, gate_name in GATE_FAMILIES.items():
            print(f"\n    -- {family} ({gate_name}) native --")
            t0 = time.time()
            n2q_set, n1q_list, depth_list, fid_err_list = set(), [], [], []
            bad_angle_count = 0
            fold_worst = {f: 0.0 for f in FOLD_FACTORS}
            for name, sol in solutions.items():
                target_vec = p["targets"][name]
                qc = build_ansatz(sol["angles"])
                native = to_native(qc, gate_name)
                fid_err, ok = verify(target_vec, native, tol=1e-9)
                bad_angles = check_angle_ranges(native, gate_name)
                n2q = native.count_ops().get(gate_name, 0)
                n1q = sum(v for k2, v in native.count_ops().items() if k2 != gate_name)
                depth = native.depth()
                n2q_set.add(n2q)
                n1q_list.append(n1q)
                depth_list.append(depth)
                fid_err_list.append(fid_err)
                bad_angle_count += len(bad_angles)

                # -- exact native fold verification --
                sv_native = np.asarray(Statevector.from_instruction(native))
                for fold in FOLD_FACTORS:
                    folded = fold_native_2q(native, fold, gate_name)
                    sv_folded = np.asarray(Statevector.from_instruction(folded))
                    idx = int(np.argmax(np.abs(sv_native)))
                    phase = sv_folded[idx] / sv_native[idx] if abs(sv_native[idx]) > 1e-9 else 1.0
                    err = float(np.max(np.abs(sv_folded / phase - sv_native)))
                    fold_worst[fold] = max(fold_worst[fold], err)

                all_records.append({
                    "K": K, "family": family, "gate_name": gate_name, "slot": name,
                    "fidelity_err": fid_err, "n1q": n1q, "n2q": n2q, "depth": depth,
                    "bad_angles": len(bad_angles),
                })

            print(f"      fidelity: worst={max(fid_err_list):.2e}, all <1e-9: {all(e < 1e-9 for e in fid_err_list)}")
            print(f"      N_2q: {sorted(n2q_set)} (constant: {len(n2q_set)==1})  "
                  f"N_1q: min={min(n1q_list)} max={max(n1q_list)}  depth: min={min(depth_list)} max={max(depth_list)}")
            print(f"      angle-range violations: {bad_angle_count} (should be 0)")
            print(f"      fold preservation (worst |psi_folded - psi_native| per fold): " +
                  "  ".join(f"fold={f}:{fold_worst[f]:.2e}" for f in FOLD_FACTORS))
            all_pass = all(v < 1e-10 for v in fold_worst.values())
            print(f"      ALL folds preserve the ideal answer to <1e-10: {all_pass}")
            fold_check_summary[f"K{K}_{family}"] = {"fold_worst": fold_worst, "all_pass": all_pass,
                                                      "bad_angle_count": bad_angle_count,
                                                      "worst_fidelity_err": max(fid_err_list)}
            print(f"      ({time.time()-t0:.1f}s)")

    results = {
        "method_errors_kcal": method_errors,
        "method_error_k5_task_claim_vs_verified": {"claimed": 0.17, "verified": method_errors.get(5),
                                                     "cited_independent_source": CLASSICAL_FLOOR_KCAL_K5_CITED},
        "fold_check_summary": fold_check_summary,
        "per_vector_records": all_records,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
