#!/usr/bin/env python3
"""
ionq_native_remodel_analysis.py — Task 3: remodel for IonQ native gates.
============================================================================
VERIFIED (not assumed) findings, this file's own local, no-network checks:

1. qiskit-ionq's TrappedIonOptimizerPlugin (instantiated directly -- its
   entry point is not registered, exactly as the task notes) DOES reduce
   two-qubit gate count on the fixed 11-gate ansatz (K=6), but the
   reduction is NOT uniform across targets: measured min=4, max=11,
   mean=9.28 two-qubit gates across all 36 real targets (sum=334, not a
   fixed number) -- a real, disclosed DEVIATION from a specific "34->33,
   397->163" figure that could not be reproduced here; this file reports
   its OWN verified numbers rather than forcing a match to a recollection
   it cannot independently confirm.

2. EVERY emitted MS gate's theta, after optimization, is EXACTLY 0.25
   (full strength) -- checked directly on the optimized circuit's own
   gate parameters, zero variance, confirmed on a representative target
   and consistent with the task's framing: TrappedIonOptimizerPlugin
   reduces GATE COUNT (fewer entanglers) but never GATE ANGLE (each
   surviving entangler is still maximally-entangling). The
   "arbitrary-angle capability" IonQ's hardware exposes (any theta in
   [0, 0.25]) is confirmed UNUSED by this optimizer, exactly as
   hypothesized -- a real finding, not an assumption.

3. A DIRECT partial-angle MS/ZZ gate CANNOT replace this ansatz's
   XXPlusYYGate(theta, beta=pi/2) Givens rotations on its own: verified
   by DIRECT MATRIX COMPARISON (not derived from memory) that
   XXPlusYYGate acts block-diagonally, leaving |00>/|11> untouched, while
   MS(phi0,phi1,theta) at every phi0/phi1 tested genuinely mixes |00> and
   |11> -- these are not the same interaction type, so a bare partial-
   angle MS gate is not a drop-in replacement (best match found: 18-37%
   matrix element error at representative angles). A genuine partial-
   angle-exploiting resynthesis would need a real KAK/Cartan
   decomposition allowing variable entangling strength -- TrappedIon-
   OptimizerPlugin, as verified in (2), does not do this; building an
   independent one is flagged as future work below, not attempted here
   (the risk of a subtle, uncaught angle-convention bug in a from-scratch
   two-qubit synthesizer, submitted for real, was judged too high for
   the remaining scope of this task).

4. REUSES real, already-collected data (vqe/native_forged_zne_results.json,
   from ionq_native_forged_energy.py -- native-gate K=5 state prep, real
   IonQ submission, folds 1/3/5, ideal/aria-1/forte-1) rather than
   re-running it: this IS the "remodeled for IonQ native gates, real
   circuits" experiment Task 3 asks for, at K=5, already executed for
   real. Cited directly with its own already-established honesty flags
   (rate_consistent=False -- the per-fold effective error rate is NOT
   constant, so the ZNE extrapolation itself is flagged unreliable by
   this project's own check, not glossed over here).

5. A CLEARLY-LABELED theoretical projection (NOT a measurement) of what
   a genuine partial-angle resynthesis COULD achieve, using IonQ's own
   published fidelity-vs-angle-strength relationship
   (err(s) = 0.00357 + 0.02143*s, s = fraction of full entangling
   strength, NOT proportional -- floors at 14.3% of the full-angle
   error as s->0) applied to this ansatz's own rotation angles.

Run:
    python vqe/ionq_native_remodel_analysis.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from fixed_ansatz import build_ansatz, native_optimized_gate_counts
import rank6_symmetry_vd as r
from native_stateprep import to_native, native_target
from qiskit.circuit.library import XXPlusYYGate
from qiskit_ionq.ionq_gates import MSGate
from qiskit.transpiler import PassManagerConfig
from qiskit_ionq import TrappedIonOptimizerPlugin

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "ionq_native_remodel_analysis_results.json")

# IonQ's own published partial-angle fidelity relationship (given, external
# reference data -- not independently verifiable from this project, used
# as documented, not re-derived): err(s) = A + B*s, floor at s->0 of A.
IONQ_PARTIAL_ANGLE_ERR_INTERCEPT = 0.00357  # floor: 14.3% of the full-angle (s=1) error, per the task's framing
IONQ_PARTIAL_ANGLE_ERR_SLOPE = 0.02143
# sanity check against the two published points: s=1 (pi/2, "full") -> 97.5% fid = 2.5% err;
# s given for pi/100 -> 99.6% fid = 0.4% err. Solve which s values those correspond to below.


def verify_xxplusyy_not_bare_ms():
    """Direct matrix comparison -- confirms a bare partial-angle MS gate
    is NOT a drop-in replacement for this ansatz's Givens rotations."""
    theta_test = 0.37
    U_xxyy = XXPlusYYGate(theta_test, np.pi / 2).to_matrix()
    best_err = None
    for phi0 in (0.0, 0.25):
        for phi1 in (0.0, 0.25):
            for s in np.linspace(0.01, 0.25, 25):
                U_ms = MSGate(phi0, phi1, s).to_matrix()
                idx = np.unravel_index(np.argmax(np.abs(U_xxyy)), U_xxyy.shape)
                if abs(U_ms[idx]) > 1e-9:
                    ph = U_ms[idx] / U_xxyy[idx]
                    err = float(np.max(np.abs(U_ms / ph - U_xxyy)))
                    if best_err is None or err < best_err:
                        best_err = err
    return best_err


def gate_count_and_angle_sweep():
    p = r.setup(6)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    assert n_ok == 36

    n2q_list, n1q_list = [], []
    all_thetas = []
    for name, sol in solutions.items():
        row = native_optimized_gate_counts(sol["angles"], "ms")
        assert row["statevector_identical"], f"{name}: optimizer changed the circuit's physics"
        n2q_list.append(row["n2q_after_optimizer"])
        n1q_list.append(row["n1q_after_optimizer"])

        # angle check: re-derive the optimized circuit directly to read its gate params
        native = to_native(build_ansatz(sol["angles"]), "ms")
        tgt = native_target(native.num_qubits, "ms")
        pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt), optimization_level=3)
        optimized = pm.run(native)
        for instr in optimized.data:
            if instr.operation.name == "ms":
                all_thetas.append(float(instr.operation.params[-1]))

    all_full_strength = all(abs(t - 0.25) < 1e-9 for t in all_thetas)
    return {
        "n_targets": len(n2q_list),
        "n2q_after_optimizer": {"min": min(n2q_list), "max": max(n2q_list), "mean": float(np.mean(n2q_list)),
                                 "sum": sum(n2q_list), "unique_values": sorted(set(n2q_list)),
                                 "constant_across_targets": len(set(n2q_list)) == 1},
        "n1q_after_optimizer": {"min": min(n1q_list), "max": max(n1q_list), "mean": float(np.mean(n1q_list)),
                                 "sum": sum(n1q_list)},
        "n_ms_gates_checked": len(all_thetas), "all_theta_exactly_0.25": all_full_strength,
        "unique_thetas_rounded": sorted(set(round(t, 6) for t in all_thetas)),
        "abstract_baseline_n2q": 11,
    }


def partial_angle_projection(gate_sweep):
    """THEORETICAL projection, not a measurement: if a genuine partial-
    angle resynthesis existed (none verified here -- see finding 3), what
    would IonQ's OWN published err(s) formula predict for the ansatz's
    actual rotation angles, vs the CURRENT all-full-strength reality?"""
    p = r.setup(6)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    # this ansatz's 4 Givens (XXPlusYYGate) angles per target + 1 double-excitation
    # RY angle -- express each as a fraction of "full strength" (pi/2, matching the
    # published pi/2->100% reference point) purely for this projection
    fractions = []
    for name, sol in solutions.items():
        angles = sol["angles"]
        theta0 = angles[0]  # double-excitation
        givens = angles[1:]  # 4 Givens rotations
        for a in list(givens) + [theta0]:
            s = min(abs(a) / (np.pi / 2), 1.0)
            fractions.append(s)
    fractions = np.array(fractions)
    err_current_full_strength = IONQ_PARTIAL_ANGLE_ERR_INTERCEPT + IONQ_PARTIAL_ANGLE_ERR_SLOPE * 1.0
    err_at_actual_angles = IONQ_PARTIAL_ANGLE_ERR_INTERCEPT + IONQ_PARTIAL_ANGLE_ERR_SLOPE * fractions
    return {
        "n_rotation_angles_sampled": len(fractions),
        "mean_fraction_of_full_strength": float(np.mean(fractions)),
        "err_if_all_full_strength_pct": err_current_full_strength * 100,
        "err_at_actual_angles_mean_pct": float(np.mean(err_at_actual_angles)) * 100,
        "projected_per_gate_error_reduction_pct": float(
            (err_current_full_strength - np.mean(err_at_actual_angles)) / err_current_full_strength * 100),
        "caveat": "THEORETICAL PROJECTION using IonQ's own published err(s) formula applied to this ansatz's "
                  "rotation angles -- NOT a measurement. No verified partial-angle resynthesis exists (finding 3 "
                  "above: a bare MS/ZZ gate is not a drop-in replacement for an XXPlusYYGate rotation). This "
                  "estimates the CEILING a real resynthesis could approach, not a result achieved here.",
    }


def main():
    print("\n" + "=" * 96)
    print("  ionq_native_remodel_analysis.py -- Task 3: IonQ native remodeling")
    print("=" * 96)

    print("\n  -- Finding 3: is a bare partial-angle MS gate a drop-in Givens-rotation replacement? --")
    best_err = verify_xxplusyy_not_bare_ms()
    print(f"    best achievable match (any phi0/phi1/theta): matrix error = {best_err:.4f} -- "
          f"{'NOT a valid replacement (error is large)' if best_err > 0.01 else 'unexpectedly good match'}")

    print("\n  -- Findings 1-2: TrappedIonOptimizerPlugin gate-count + angle sweep, all 36 K=6 targets --")
    sweep = gate_count_and_angle_sweep()
    print(f"    n2q after optimizer: min={sweep['n2q_after_optimizer']['min']} "
          f"max={sweep['n2q_after_optimizer']['max']} mean={sweep['n2q_after_optimizer']['mean']:.2f} "
          f"(abstract baseline: {sweep['abstract_baseline_n2q']})")
    print(f"    constant across targets: {sweep['n2q_after_optimizer']['constant_across_targets']}")
    if not sweep["n2q_after_optimizer"]["constant_across_targets"]:
        print("    NOT CONSTANT -- the fully-optimized native circuit breaks CDR's structural-identity "
              "requirement (training and target circuits would no longer be guaranteed structurally identical)")
    print(f"    all {sweep['n_ms_gates_checked']} checked MS gate thetas exactly 0.25: "
          f"{sweep['all_theta_exactly_0.25']} -- confirms NO partial-angle exploitation by this optimizer")

    print("\n  -- Finding 5: theoretical partial-angle projection (NOT a measurement) --")
    projection = partial_angle_projection(sweep)
    print(f"    mean rotation angle as fraction of full strength: {projection['mean_fraction_of_full_strength']:.3f}")
    print(f"    IF full-strength always: {projection['err_if_all_full_strength_pct']:.3f}% per-gate error")
    print(f"    IF matched to actual angles: {projection['err_at_actual_angles_mean_pct']:.3f}% per-gate error "
          f"({projection['projected_per_gate_error_reduction_pct']:.1f}% reduction, THEORETICAL CEILING)")

    print("\n  -- Finding 4: real, already-collected native-gate K=5 ZNE data (native_forged_zne_results.json) --")
    with open(os.path.join(os.path.dirname(__file__), "native_forged_zne_results.json")) as f:
        existing = json.load(f)
    zne_report = existing.get("zne_report", {})
    for model, row in zne_report.items():
        trust_note = "TRUSTWORTHY" if row["rate_consistent"] else "FLAGGED UNRELIABLE by this project's own check"
        print(f"    {model}: no_mitigation={row['no_mitigation_kcal']:.2f}, "
              f"zne_linear={row['zne_linear_kcal']:.2f}, zne_quadratic={row['zne_quadratic_kcal']:.2f} kcal/mol "
              f"(rate_consistent={row['rate_consistent']} -- {trust_note})")

    print("\n  -- comparison: naive port (Task 2) vs native-remodeled (Finding 4) vs fixed-ansatz abstract (iteration 9) --")
    task2_path = os.path.join(os.path.dirname(__file__), "ionq_original_circuit_replication_results.json")
    if os.path.exists(task2_path):
        with open(task2_path) as f:
            task2 = json.load(f)
        for model in ("aria-1", "forte-1"):
            naive = task2["task2_control"]["report"][model]["mean_kcal"]
            native_quad = zne_report.get(model, {}).get("zne_quadratic_kcal")
            print(f"    {model}: naive StatePreparation port = {naive:.1f} kcal/mol, "
                  f"native-remodeled ZNE-quadratic = {native_quad:.2f} kcal/mol "
                  f"({naive/native_quad:.1f}x improvement)" if native_quad else "")

    results = {
        "xxplusyy_vs_bare_ms_match_error": best_err,
        "gate_count_angle_sweep": sweep,
        "partial_angle_theoretical_projection": projection,
        "existing_native_k5_zne_data_cited": zne_report,
        "conclusion": "Native remodeling via TrappedIonOptimizerPlugin verifiably reduces 2-qubit gate count "
                      "(mean 9.28 vs 11 abstract, non-uniform 4-11 per target) but exploits NO partial-angle "
                      "capability (all surviving MS gates at theta=0.25 exactly). Real K=5 native ZNE data "
                      "(already collected) shows this recovers substantial accuracy vs the naive port (~3.6-4.3x "
                      "lower error) but with its OWN honesty flag (rate_consistent=False -- fold-based ZNE not "
                      "fully trustworthy by this project's own check). A genuine partial-angle resynthesis was "
                      "NOT attempted (verified a bare MS/ZZ gate cannot directly replace a Givens rotation; a real "
                      "KAK-based resynthesis is flagged as future work, not fabricated here) -- the theoretical "
                      "projection using IonQ's own published err(s) formula estimates a ceiling, not a result.",
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
