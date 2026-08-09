#!/usr/bin/env python3
"""
rc_zne_coherent_noise.py — does Randomized Compiling (Pauli twirling)
unlock a genuine ZNE plateau against a DELIBERATELY COHERENT local noise
model, the one noise regime this project's local synthetic tests have
never actually simulated?
============================================================================
WHY THIS FILE EXISTS: iteration 21 researched IonQ's own published
literature and found a striking match for this project's own repeated,
never-explained ZNE failures (iterations 11, 13, 14, 19: no plateau in
ANY tested (scale-range, fit-order) combination, on either the untapered,
tapered, or leakage-postselected circuit): the RC+ZNE literature (Quantum
5, 1184 (2023)) states plainly that "very small amounts of coherent
noise in VQE can cause substantially large errors that are difficult to
suppress by conventional mitigation methods," because coherent noise
"violates ZNE's assumption that errors scale predictably." EVERY local
noise model this project has used until now (zne_floor_tested.py,
z2_tapered_zne.py, leakage_zne_floor_tested.py) is a Qiskit Aer
`depolarizing_error` channel -- purely STOCHASTIC/incoherent by
construction. If real IonQ hardware's actual error budget includes a
genuine coherent component (plausible: over-rotation from laser/motional-
mode miscalibration, a well-documented real error source in trapped-ion
systems), this project's local tests would never have been able to
reproduce the mechanism breaking ZNE, because they never modeled it.
This file builds that missing local test directly, per explicit user
request, BEFORE spending anything on a real IonQ submission.

COHERENT NOISE MODEL: `qiskit_aer.noise.coherent_unitary_error` attaches
a single, FIXED (not probabilistic) unitary error after every CX gate --
a small residual ZZ-coupling rotation, exp(-i*eps*Z@Z/2), a standard,
physically-motivated coherent error form (over-rotation / always-on
crosstalk). Unlike depolarizing noise, this is the SAME error every
single time the gate is applied -- deterministic, not random -- exactly
what Pauli twirling is designed to address and depolarizing-only tests
cannot probe.

RANDOMIZED COMPILING: the standard CX Pauli-twirling group (16 elements)
was derived and VERIFIED NUMERICALLY here (not hand-derived from
literature/memory, avoiding this project's own repeated history of sign/
convention bugs in exactly this kind of derivation) -- for every
(P_control, P_target) input pair, matrix-multiplied CX@(P_c⊗P_t) and
brute-force searched all 16 candidate output pairs for the exact
(Pout)@CX@(Pin) = ±CX identity, confirmed for all 16 with |phase|=1.
Two pairs pick up a global -1 phase (harmless -- global phase never
affects any measurable probability or expectation value).

BUG CAUGHT AND FIXED: the FIRST derivation used the bare `Operator(
CXGate())` object's matrix directly, which turned out to use a DIFFERENT
internal qubit-index convention than `Statevector.from_instruction()`
uses for an actual circuit built with `qc.cx(control, target)` -- a
genuine qubit-ordering mismatch, the same class of bug this project has
hit before (documented in z2_tapering.py's own history). Caught by a
direct, minimal 2-qubit exact-identity test: with the FIRST derivation's
table, applying (pre-twirl, cx, post-twirl) on an arbitrary input state
did NOT reproduce the un-twirled circuit's output for 15 of 16 pairs
(errors of 0.5-2.9 in statevector norm, not small numerical noise) --
only the trivial (I,I) pair passed. Root cause confirmed and fixed by
re-deriving the twirl table from an ACTUAL 2-qubit circuit's Operator
(`QuantumCircuit(2); qc.cx(0,1); Operator(qc)`), matching the SAME
convention used everywhere else in this file and project. Re-verified:
the corrected table gives EXACTLY 0.0 error (not just small) for all 16
pairs on the same minimal test.

METHOD: exact (Statevector, no shot noise first) comparison of E(s) vs
noise scale s, WITH and WITHOUT RC-twirling averaged over N_TWIRLS random
instances, using the SAME 3-direction floor test (qforge.floor_test)
this project has used since iteration 13 for every ZNE claim.

Run:
    python vqe/rc_zne_coherent_noise.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, floor_test,
)
from fixed_ansatz import build_ansatz
from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import CXGate
from qiskit.quantum_info import Statevector, Operator, Pauli

K = 6
BASIS_GATES = ["u3", "cx"]
N_TWIRLS = 16
EPS_PER_GATE = 0.06  # coherent over-rotation angle (radians) per CX at scale=1 -- chosen so scale=1's
                      # induced infidelity is comparable in ORDER OF MAGNITUDE to this project's own
                      # calibrated P2_PER_GATE=0.01214 depolarizing model (sin^2(eps/2) ~ 0.0009 per
                      # gate at eps=0.06 -- smaller per-gate, deliberately: the RC+ZNE literature's own
                      # point is that even SMALL coherent noise causes large VQE errors, so this is not
                      # cherry-picked to be dramatic)
SCALE_RANGES = [
    [1, 2, 3],
    [1, 2, 3, 4],
    [1, 2, 3, 4, 5],
    [1, 2, 3, 4, 5, 6],
    [1, 2, 3, 4, 5, 6, 7],
]
MAX_SCALE = max(max(r) for r in SCALE_RANGES)
SHOTS = 100_000
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "rc_zne_coherent_noise_results.json")


# ---------------------------------------------------------------------------
# CX Pauli-twirling group: verified numerically (see module docstring),
# hardcoded here as the confirmed result, not re-derived at runtime.
# ---------------------------------------------------------------------------
PAULIS = ["I", "X", "Y", "Z"]
TWIRL_OUT = {
    ("I", "I"): ("I", "I"), ("I", "X"): ("I", "X"), ("I", "Y"): ("Z", "Y"), ("I", "Z"): ("Z", "Z"),
    ("X", "I"): ("X", "X"), ("X", "X"): ("X", "I"), ("X", "Y"): ("Y", "Z"), ("X", "Z"): ("Y", "Y"),
    ("Y", "I"): ("Y", "X"), ("Y", "X"): ("Y", "I"), ("Y", "Y"): ("X", "Z"), ("Y", "Z"): ("X", "Y"),
    ("Z", "I"): ("Z", "I"), ("Z", "X"): ("Z", "X"), ("Z", "Y"): ("I", "Y"), ("Z", "Z"): ("I", "Z"),
}
ALL_TWIRL_INPUTS = list(TWIRL_OUT.keys())


def apply_pauli(qc, pauli_char, qubit):
    if pauli_char == "X":
        qc.x(qubit)
    elif pauli_char == "Y":
        qc.y(qubit)
    elif pauli_char == "Z":
        qc.z(qubit)
    # I: no-op


def statevector_energy(p, angles_by_name, eps, rc, rng=None, n_twirl_avg=1):
    """Exact expectation values (Statevector, no shot noise) for every
    non-identity alpha_label, at the given coherent-error scale, averaged
    over n_twirl_avg random RC instances if rc=True (n_twirl_avg=1 if not)."""
    from qiskit.quantum_info import Pauli as QPauli
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    raw = {name: {l: 0.0 for l in non_id_labels} for name in angles_by_name}
    n_avg = n_twirl_avg if rc else 1
    for name, angles in angles_by_name.items():
        base = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
        for _ in range(n_avg):
            qc = build_coherent_circuit(base, eps, rc, rng)
            sv = Statevector.from_instruction(qc)
            for l in non_id_labels:
                val = float(np.real(sv.expectation_value(QPauli(l))))
                raw[name][l] += val / n_avg
    return raw


def build_coherent_circuit(base_qc, eps, rc, rng):
    """For every cx gate: [pre-twirl if rc] -> cx -> coherent ZZ-rotation
    unitary (exp(-i*eps*Z@Z/2), the deterministic coherent error this
    test isolates) -> [compensating post-twirl if rc]."""
    diag = np.array([1, -1, -1, 1])
    zz_unitary = np.diag(np.exp(-1j * eps / 2 * diag))
    qc = QuantumCircuit(base_qc.num_qubits)
    for instr in base_qc.data:
        op = instr.operation
        qubits = [base_qc.find_bit(q).index for q in instr.qubits]
        if op.name == "cx":
            c, t = qubits
            if rc:
                pin = ALL_TWIRL_INPUTS[rng.integers(0, len(ALL_TWIRL_INPUTS))]
                pout = TWIRL_OUT[pin]
                apply_pauli(qc, pin[0], c)
                apply_pauli(qc, pin[1], t)
            qc.cx(c, t)
            qc.unitary(zz_unitary, [t, c], label="coherent_zz_err")
            if rc:
                apply_pauli(qc, pout[0], c)
                apply_pauli(qc, pout[1], t)
        else:
            qc.append(op, qubits)
    return qc


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def exponential_fit(scales, energies):
    diffs = np.array(energies) - energies[-1]
    if np.all(diffs > 0) or np.all(diffs < 0):
        slope, intercept = np.polyfit(scales, np.log(np.abs(diffs)), 1)
        A = np.sign(diffs[0]) * np.exp(intercept)
        return float(energies[-1] + A), None
    return None, "non-monotonic"


def run_floor_test_suite(scheme_name, exact_raw_by_scale, p, angles_by_name):
    """Shot-noisy 8-seed sweep + the SAME mandatory 3-direction floor test
    used since iteration 13, applied to whichever scheme's exact_raw_by_scale
    dict is passed in (RC-averaged or not)."""
    from qforge import shot_sample
    range_order_results, range_exp_results = {}, {}
    for scale_range in SCALE_RANGES:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            range_order_results[(range_key, order)] = []
        range_exp_results[range_key] = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(seed * 7919 + 13)
            energies = []
            for s in scale_range:
                raw = {name: {l: shot_sample(v, SHOTS, rng) for l, v in vals.items()}
                       for name, vals in exact_raw_by_scale[s].items()}
                E, _ = energy_and_err(p, raw, K)
                energies.append(E)
            for order in range(1, max_order + 1):
                coeffs = np.polyfit(scale_range, energies, order)
                E0 = float(np.polyval(coeffs, 0))
                err = abs(E0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
                range_order_results[(range_key, order)].append(err)
            E0_exp, _ = exponential_fit(scale_range, energies)
            if E0_exp is not None:
                range_exp_results[range_key].append(abs(E0_exp - p["exact_energy"]) * HARTREE_TO_KCAL_MOL)

    summary = {}
    for scale_range in SCALE_RANGES:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            errs = range_order_results[(range_key, order)]
            summary[f"{range_key}|order{order}"] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs))}
        if range_exp_results[range_key]:
            errs = range_exp_results[range_key]
            summary[f"{range_key}|exp"] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs))}

    range_sizes = [len(r) for r in SCALE_RANGES]
    lin_errs = [summary[f"{','.join(str(s) for s in r)}|order1"]["mean_kcal"] for r in SCALE_RANGES]
    ft_range_linear = floor_test(range_sizes, lin_errs)

    widest = SCALE_RANGES[-1]
    widest_key = ",".join(str(s) for s in widest)
    orders = list(range(1, len(widest)))
    order_errs = [summary[f"{widest_key}|order{o}"]["mean_kcal"] for o in orders]
    ft_order = floor_test(orders, order_errs)

    quad_ranges = [r for r in SCALE_RANGES if len(r) > 2]
    quad_errs = [summary[f"{','.join(str(s) for s in r)}|order2"]["mean_kcal"] for r in quad_ranges]
    quad_sizes = [len(r) for r in quad_ranges]
    ft_range_quad = floor_test(quad_sizes, quad_errs)

    any_plateau = ft_range_linear["has_floor"] or ft_order["has_floor"] or ft_range_quad["has_floor"]
    print(f"\n  [{scheme_name}] range@order1: {ft_range_linear['verdict']}")
    print(f"  [{scheme_name}] order@widest-range: {ft_order['verdict']}")
    print(f"  [{scheme_name}] range@order2: {ft_range_quad['verdict']}")
    print(f"  [{scheme_name}] ANY PLATEAU: {any_plateau}")
    return {"summary": summary, "any_plateau_found": bool(any_plateau),
            "floor_test_range_at_order1": ft_range_linear, "floor_test_order_at_widest_range": ft_order,
            "floor_test_range_at_order2": ft_range_quad}


def main():
    print("\n" + "=" * 96)
    print("  rc_zne_coherent_noise.py -- does RC unlock a ZNE plateau against COHERENT noise?")
    print("=" * 96)

    t0 = time.time()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    assert n_ok == 36
    angles_by_name = {name: sol["angles"] for name, sol in solutions.items()}
    print(f"  setup OK: 36/36 converged, {time.time()-t0:.1f}s")

    print(f"\n  computing exact (no shot noise) curves, scales 1..{MAX_SCALE}, "
          f"coherent eps_per_gate={EPS_PER_GATE} rad, N_TWIRLS={N_TWIRLS}")
    t0 = time.time()
    exact_no_rc, exact_rc = {}, {}
    rng = np.random.default_rng(42)
    for s in range(1, MAX_SCALE + 1):
        eps = EPS_PER_GATE * s
        exact_no_rc[s] = statevector_energy(p, angles_by_name, eps, rc=False)
        exact_rc[s] = statevector_energy(p, angles_by_name, eps, rc=True, rng=rng, n_twirl_avg=N_TWIRLS)
        _, err_no_rc = energy_and_err(p, exact_no_rc[s], K)
        _, err_rc = energy_and_err(p, exact_rc[s], K)
        print(f"    scale={s}: NO-RC err={err_no_rc:.3f}  RC-averaged err={err_rc:.3f} kcal/mol (exact, no shot noise)")
    print(f"  exact curves done, {time.time()-t0:.1f}s")

    print(f"\n  -- 8-seed shot-noisy floor test sweep --")
    t0 = time.time()
    no_rc_results = run_floor_test_suite("NO-RC (coherent noise, raw)", exact_no_rc, p, angles_by_name)
    rc_results = run_floor_test_suite("RC-twirled (coherent noise, averaged)", exact_rc, p, angles_by_name)
    print(f"  floor test sweep done, {time.time()-t0:.1f}s")

    print(f"\n  -- OVERALL VERDICT --")
    if rc_results["any_plateau_found"] and not no_rc_results["any_plateau_found"]:
        print("  RC UNLOCKS a genuine plateau against coherent noise that raw ZNE never finds. "
              "This directly supports the RC+ZNE literature's mechanism and is worth testing for "
              "real on IonQ's free simulator.")
    elif rc_results["any_plateau_found"] and no_rc_results["any_plateau_found"]:
        print("  Both RC and non-RC find a plateau against this specific coherent noise model -- "
              "inconclusive for isolating RC's specific contribution.")
    else:
        print("  NO PLATEAU for RC either, against a genuine coherent noise model. The RC+ZNE "
              "mechanism, even though real and published, does not rescue THIS project's specific "
              "ansatz/reconstruction pipeline at this local noise level -- reported plainly.")

    results = {
        "K": K, "shots": SHOTS, "n_seeds": N_SEEDS, "scale_ranges": SCALE_RANGES,
        "eps_per_gate": EPS_PER_GATE, "n_twirls": N_TWIRLS,
        "no_rc": no_rc_results, "rc": rc_results,
        "exact_curves": {
            "no_rc": {s: energy_and_err(p, exact_no_rc[s], K)[1] for s in exact_no_rc},
            "rc": {s: energy_and_err(p, exact_rc[s], K)[1] for s in exact_rc},
        },
        "chemical_accuracy_kcal": 1.0,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
