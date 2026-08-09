#!/usr/bin/env python3
"""
z2_tapered_fixed_ansatz.py — hand-derive a fixed-structure circuit for
the Z2-tapered (3-qubit) register, mirroring fixed_ansatz.py's own
derivation for the original 4-qubit problem.
============================================================================
DIAGNOSIS FIRST (gate_structure_compare.py): the tapered circuit
(generic StatePreparation) has FEWER gates of BOTH types than the
abstract 11-gate ansatz (2q: mean 3.72 vs constant 11; 1q: mean 6.72 vs
constant 51) yet performs WORSE on real IonQ hardware. Direct inspection
of the actual gate lists shows why gate COUNT alone misleads: the
abstract ansatz's 51 "1-qubit gates" are overwhelmingly FIXED,
SPECIAL-ANGLE gates (0, +-pi/2, +-pi combinations -- decomposed
Hadamard/S/Sdg/X structural gates from the Givens-rotation recipe),
with only a handful actually carrying the fitted target angles. The
tapered circuit's gates are almost ALL generic, arbitrary-angle
rotations straight out of StatePreparation's synthesis algorithm. The
hypothesis this file tests: rebuilding the tapered register with the
SAME kind of structural, mostly-special-angle Givens-rotation recipe
(not just "fewer gates") may recover real-hardware performance closer
to the abstract ansatz's, even if the raw CX count does not drop
further.

CONSTRUCTION: the reduced problem occupies exactly 6 of 8 basis states
(indices 1,2,3,4,5,6 -- weight-1: {1,2,4}, weight-2: {3,5,6}). Index 4
(|100>) and index 3 (|011>) are EXACT BIT-COMPLEMENTS (every bit flips)
-- the SAME structural relationship fixed_ansatz.py's own double-
excitation block exploits for the original 4-qubit problem, just one
qubit down (2 fan-out CX instead of 3). Reference: X on qubit 2 (index
4, weight-1). Bridge: RY(2*theta0) on qubit 2 + CX(2,1) + CX(2,0) ->
cos(theta0)|100> + sin(theta0)|011>. Then Givens (XXPlusYYGate) hops
across the 3 available qubit pairs (0,1),(0,2),(1,2) to reach the
remaining 4 target basis states.

Every candidate is verified EXACTLY (least-squares fit to <1e-10) against
the REAL, ALREADY-VERIFIED tapered target vectors from z2_tapered_zne.py
before being trusted -- not assumed to work by construction.

Run:
    python vqe/z2_tapered_fixed_ansatz.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from z2_tapered_zne import build_reduced_problem
from qforge import slot_names
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import XXPlusYYGate
from qiskit.quantum_info import Statevector
from qiskit import transpile
from scipy.optimize import least_squares

BASIS_GATES = ["u3", "cx"]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "z2_tapered_fixed_ansatz_results.json")


def build_candidate(angles):
    """5 angles: theta0 (bridge), theta1/theta2/theta3 (Givens hops).
    Only 3 distinct qubit pairs exist on 3 qubits, so this uses 3 Givens
    gates (not 4, unlike the original 4-qubit ansatz) -- 1+3=4 angles
    tested first; a 5th (repeating a pair) is added only if 4 does not
    suffice, checked empirically below, not assumed necessary.

    BUG CAUGHT AND FIXED before this worked: the discriminator qubit for
    the bridge (RY + CX fan-out) must be a qubit the reference prep
    NEVER touches (starts at |0>), not the same qubit the reference X
    gate just set to 1 -- confusing the two gave a circuit that only
    ever reached |000>/|111> (verified directly via Statevector before
    diagnosing the mismatch, not assumed from the failed fit alone).
    Reference: X on qubit 2 -> |100> (index 4). Discriminator: qubit 0
    (untouched, starts at |0>, matching index 4's own q0=0). Fan-out:
    CX(0,1), CX(0,2) flips q1 AND q2 when q0=1, landing on q1=1,q2=0 --
    together with q0=1 that is |011> (index 3), the bit-complement of
    index 4, exactly mirroring fixed_ansatz.py's own derivation."""
    theta0, theta1, theta2, theta3 = angles[:4]
    qc = QuantumCircuit(3)
    qc.x(2)  # reference q2=1, q1=0, q0=0 -> |100> = index 4
    qc.ry(2 * theta0, 0)  # discriminator: qubit 0, untouched so far, starts at |0>
    qc.cx(0, 1)
    qc.cx(0, 2)
    qc.append(XXPlusYYGate(theta1, np.pi / 2), [0, 1])
    qc.append(XXPlusYYGate(theta2, np.pi / 2), [0, 2])
    qc.append(XXPlusYYGate(theta3, np.pi / 2), [1, 2])
    return qc


def build_candidate5(angles):
    theta0, theta1, theta2, theta3, theta4 = angles
    qc = QuantumCircuit(3)
    qc.x(2)
    qc.ry(2 * theta0, 0)
    qc.cx(0, 1)
    qc.cx(0, 2)
    qc.append(XXPlusYYGate(theta1, np.pi / 2), [0, 1])
    qc.append(XXPlusYYGate(theta2, np.pi / 2), [0, 2])
    qc.append(XXPlusYYGate(theta3, np.pi / 2), [1, 2])
    qc.append(XXPlusYYGate(theta4, np.pi / 2), [0, 1])
    return qc


def statevector_of(build_fn, angles):
    return np.asarray(Statevector.from_instruction(build_fn(angles)))


def max_abs_error(build_fn, angles, target):
    sv = statevector_of(build_fn, angles)
    idx = int(np.argmax(np.abs(target)))
    phase = sv[idx] / target[idx] if abs(target[idx]) > 1e-9 else 1.0
    return float(np.max(np.abs(sv / phase - target)))


def residual_vector(angles, build_fn, target):
    sv = statevector_of(build_fn, angles)
    idx = int(np.argmax(np.abs(target)))
    phase = sv[idx] / target[idx] if abs(target[idx]) > 1e-9 else 1.0
    diff = sv / phase - target
    return np.concatenate([diff.real, diff.imag])


def fit_angles(build_fn, n_params, target, n_attempts=12, tol=1e-10):
    """A real, caught bug: the phase-alignment convention (phase =
    sv[idx]/target[idx], matching fixed_ansatz.py's own established
    pattern) divides by sv[idx], which CAN land near-zero at an unlucky
    random initial guess -- least_squares then raises ("Residuals are
    not finite") and previously crashed the whole fit instead of just
    that one attempt. Each attempt is now isolated: a bad initial guess
    is skipped and the next random seed is tried, matching the spirit of
    fixed_ansatz.py's multi-attempt retry loop, which happened to never
    trigger this edge case but was never actually protected against it
    either."""
    best_err, best_x = None, None
    for attempt in range(n_attempts):
        rng = np.random.default_rng(attempt)
        x0 = rng.uniform(-np.pi, np.pi, n_params)
        try:
            res = least_squares(residual_vector, x0, args=(build_fn, target), method="lm",
                                 xtol=1e-15, ftol=1e-15, gtol=1e-15)
        except ValueError:
            continue  # unlucky initial guess (near-zero denominator in the phase alignment) -- skip, try next seed
        err = max_abs_error(build_fn, res.x, target)
        if best_err is None or err < best_err:
            best_err, best_x = err, res.x
        if best_err < tol:
            break
    if best_x is None:
        return np.zeros(n_params), float("inf")
    return best_x, best_err


def main():
    print("\n" + "=" * 96)
    print("  z2_tapered_fixed_ansatz.py -- hand-derived circuit for the tapered register")
    print("=" * 96)

    problem = build_reduced_problem()
    reduced_targets = problem["reduced_targets"]
    K = 6
    names = slot_names(K)

    for label, build_fn, n_params in (("4-angle (theta0 + 3 Givens)", build_candidate, 4),
                                       ("5-angle (theta0 + 4 Givens, 1 pair reused)", build_candidate5, 5)):
        print(f"\n  -- trying {label} --")
        worst = 0.0
        n_ok = 0
        solutions = {}
        for name in names:
            angles, err = fit_angles(build_fn, n_params, reduced_targets[name])
            ok = err < 1e-10
            n_ok += ok
            worst = max(worst, err)
            solutions[name] = {"angles": angles.tolist(), "err": err, "ok": bool(ok)}
        print(f"    {n_ok}/{len(names)} converged to <1e-10, worst={worst:.2e}")
        if n_ok == len(names):
            print(f"    ALL TARGETS CONVERGED with {label} -- checking gate count")
            n2q_list, n1q_list = [], []
            for name in names:
                qc = build_fn(solutions[name]["angles"])
                t = transpile(qc, basis_gates=BASIS_GATES, optimization_level=0)
                ops = t.count_ops()
                n2q_list.append(ops.get("cx", 0))
                n1q_list.append(ops.get("u3", 0))
            print(f"    2q gates: {set(n2q_list)} (constant: {len(set(n2q_list))==1})")
            print(f"    1q gates: min={min(n1q_list)} max={max(n1q_list)} mean={np.mean(n1q_list):.2f} "
                  f"(constant: {len(set(n1q_list))==1})")
            results = {
                "label": label, "n_params": n_params, "n_ok": n_ok, "worst_err": worst,
                "solutions": solutions,
                "n2q_gates": sorted(set(n2q_list)), "n2q_constant": len(set(n2q_list)) == 1,
                "n1q_gates_min": min(n1q_list), "n1q_gates_max": max(n1q_list),
                "n1q_gates_mean": float(np.mean(n1q_list)), "n1q_constant": len(set(n1q_list)) == 1,
                "build_fn": label,
            }
            with open(RESULTS_PATH, "w") as f:
                json.dump(results, f, indent=2)
            print(f"\n  Results saved -> {RESULTS_PATH}\n")
            return results, build_fn, solutions

    print("\n  NEITHER candidate converged for all 36 targets -- reporting the honest failure, not forcing it")
    return None, None, None


if __name__ == "__main__":
    main()
