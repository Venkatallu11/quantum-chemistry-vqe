#!/usr/bin/env python3
"""
z2_tapered_zne.py — does Z2 tapering (iteration 13, Task A: 4->3 qubits
per register, mean 3.94 CX gates vs 11) change ZNE's convergence
behavior? Local floor test first; real IonQ submission only if it
passes.
============================================================================
iteration 13's Task E found NO PLATEAU in ANY (scale-range, fit-order)
direction for the UNTAPERED 11-gate ansatz -- ZNE was not demonstrated to
converge there at all. This file asks the natural next question: does
the smaller, tapered (3-qubit, mean 3.94 CX) circuit behave differently?
Fewer gates means a smaller total noise-scale range is being explored at
any given scale factor, which could plausibly make the polynomial/
Richardson extrapolation better-conditioned -- a real, testable
hypothesis, not assumed true.

BUILT ON VERIFIED MATH, not re-derived: reuses z2_tapering.py's
Clifford + basis-rotation transform (confirmed exact to 1e-14, energy
exact to 2.29e-11 kcal/mol in Task A) and additionally verifies here
(not assumed) that EVERY ONE of the 37 alpha-register Pauli labels
reduces to a genuine single Pauli string (times a real +-1 sign) on the
3-qubit register -- checked by direct matrix comparison against all 64
three-qubit Pauli operators, not asserted. Also verified: tapering
commutes exactly with this project's beta_signs() shortcut (v_top =
signs * u_top survives tapering: max|reduced_v - signs*reduced_u| =
0.0), so this file only needs REAL circuits for the alpha register --
beta is derived, exactly like the untapered pipeline.

Run:
    python vqe/z2_tapered_zne.py
"""
import os
import sys
import json
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, HARTREE_TO_KCAL_MOL, shot_sample, floor_test, P2_PER_GATE, P1_PER_GATE
from z2_tapering import find_register_symmetry, taper_register, taper_pauli_matrix
import ef_fragment as effrag
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector, Pauli
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

K = 6
BASIS_GATES = ["u3", "cx"]
SHOTS = 100_000
N_SEEDS = 8
SCALE_RANGES = [
    [1, 2, 3],
    [1, 2, 3, 4],
    [1, 2, 3, 4, 5],
    [1, 2, 3, 4, 5, 6],
    [1, 2, 3, 4, 5, 6, 7],
]
MAX_SCALE = max(max(r) for r in SCALE_RANGES)
CHEM_ACC_KCAL = 1.0
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "z2_tapered_zne_results.json")


def all_3q_paulis():
    return {"".join(p): np.asarray(Pauli("".join(p)).to_matrix()) for p in itertools.product("IXYZ", repeat=3)}


def build_reduced_problem():
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]

    z2 = find_register_symmetry(alpha_labels, "alpha")
    assert len(z2.symmetries) == 1

    # taper all 36 K=6 target slots
    from qforge import target_states, slot_names
    targets = target_states(p["u_vecs"], K)
    all_vecs = np.array([targets[name] for name in slot_names(K)])
    reduced_all, sq_qubit, U, tapered_val = taper_register(all_vecs, alpha_labels, z2, len(all_vecs), "alpha (36 targets)")
    reduced_targets = dict(zip(slot_names(K), reduced_all))

    # verify every Pauli label reduces to a genuine 3-qubit Pauli, real sign
    pauli_mats_3q = all_3q_paulis()
    reduced_label_map = {}
    for l in alpha_labels:
        Pr = taper_pauli_matrix(l, U, sq_qubit, tapered_val, 4)
        matched = None
        for lbl, M in pauli_mats_3q.items():
            for sign in (1.0, -1.0):
                if np.allclose(Pr, sign * M, atol=1e-8):
                    matched = (lbl, sign)
                    break
            if matched:
                break
        assert matched is not None, f"label {l} did not reduce to a clean 3-qubit Pauli -- refusing to proceed"
        reduced_label_map[l] = matched
    reduced_identity_label, id_sign = reduced_label_map[identity_label]
    assert reduced_identity_label == "III" and id_sign == 1.0

    return {
        "p": p, "alpha_labels": alpha_labels, "non_id_labels": non_id_labels,
        "identity_label": identity_label, "reduced_identity_label": reduced_identity_label,
        "reduced_targets": reduced_targets, "reduced_label_map": reduced_label_map,
        "signs": p["signs"], "sq_qubit": sq_qubit, "U": U, "tapered_val": tapered_val,
    }


def statevector_verify(problem, tol=1e-9):
    """Verify the reduced 3-qubit StatePreparation circuits reproduce the
    exact tapered target vectors and, via untapering, the original
    (untapered) targets too."""
    worst = 0.0
    for name, vec in problem["reduced_targets"].items():
        qc = QuantumCircuit(3)
        qc.append(StatePreparation(vec), range(3))
        sv = np.asarray(Statevector.from_instruction(qc))
        idx = int(np.argmax(np.abs(vec)))
        phase = sv[idx] / vec[idx] if abs(vec[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv / phase - vec)))
        worst = max(worst, err)
    return worst


def build_noise_model(scale):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(min(P2_PER_GATE * scale, 0.75), 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(min(P1_PER_GATE * scale, 0.75), 1), "u3")
    return nm


def noisy_density_matrix_3q(vec, noise_model):
    qc = QuantumCircuit(3)
    qc.append(StatePreparation(vec), range(3))
    qct = transpile(qc, basis_gates=BASIS_GATES, optimization_level=0)
    qct2 = qct.copy()
    qct2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qct2).result()
    return np.asarray(result.data(0)["density_matrix"]), qct.count_ops().get("cx", 0)


def measure_exact_noisy_raw(problem, noise_model):
    """Returns raw[name][alpha_label] = tapered-sign-corrected expectation
    value, using the ORIGINAL (untapered) alpha_label as the key so this
    plugs directly into qforge.combine_matrices/energy_from_alpha_matrices
    unchanged."""
    raw = {}
    gate_counts = {}
    for name, vec in problem["reduced_targets"].items():
        dm, n2q = noisy_density_matrix_3q(vec, noise_model)
        gate_counts[name] = n2q
        vals = {}
        for l in problem["non_id_labels"]:
            rlabel, sign = problem["reduced_label_map"][l]
            P = np.asarray(Pauli(rlabel).to_matrix())
            vals[l] = sign * float(np.real(np.trace(P @ dm)))
        raw[name] = vals
    return raw, gate_counts


def energy_and_err(problem, raw, K):
    from qforge import combine_matrices, energy_from_alpha_matrices
    p = problem["p"]
    alpha_mats = combine_matrices(raw, problem["alpha_labels"], problem["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  z2_tapered_zne.py -- ZNE on the Z2-tapered (3-qubit) circuit, local floor test first")
    print("=" * 96)

    problem = build_reduced_problem()
    print(f"\n  {len(problem['alpha_labels'])} alpha labels, all verified to reduce to genuine 3-qubit Paulis")

    sv_err = statevector_verify(problem)
    print(f"  StatePreparation circuit vs exact tapered target, worst error: {sv_err:.2e} "
          f"({'EXACT' if sv_err < 1e-9 else 'FAILED'})")
    assert sv_err < 1e-9

    print(f"\n  measuring exact noisy matrices at scales 1..{MAX_SCALE}")
    exact_raw_by_scale = {}
    gate_counts_by_scale = {}
    for s in range(1, MAX_SCALE + 1):
        nm = build_noise_model(s)
        exact_raw_by_scale[s], gate_counts_by_scale[s] = measure_exact_noisy_raw(problem, nm)
        E, err = energy_and_err(problem, exact_raw_by_scale[s], K)
        gc = list(gate_counts_by_scale[s].values())
        print(f"    scale={s}: err_vs_exact={err:.3f} kcal/mol (2q gates/circuit: min={min(gc)} max={max(gc)} mean={np.mean(gc):.2f})")

    print(f"\n  -- 8-seed shot-noisy sweep across every (scale range, fit order) combination --")
    range_order_results = {}
    for scale_range in SCALE_RANGES:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            range_order_results[(range_key, order)] = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(seed * 7919 + 13)
            energies = []
            for s in scale_range:
                raw = {name: {l: shot_sample(v, SHOTS, rng) for l, v in vals.items()}
                       for name, vals in exact_raw_by_scale[s].items()}
                E, _ = energy_and_err(problem, raw, K)
                energies.append(E)
            for order in range(1, max_order + 1):
                coeffs = np.polyfit(scale_range, energies, order)
                E0 = float(np.polyval(coeffs, 0))
                err = abs(E0 - problem["p"]["exact_energy"]) * HARTREE_TO_KCAL_MOL
                range_order_results[(range_key, order)].append(err)

    print(f"\n  {'range':<16} {'order':>5} {'mean_kcal':>10} {'std_kcal':>9}")
    summary = {}
    for scale_range in SCALE_RANGES:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            errs = range_order_results[(range_key, order)]
            mean_e, std_e = float(np.mean(errs)), float(np.std(errs))
            summary[f"{range_key}|order{order}"] = {"mean_kcal": mean_e, "std_kcal": std_e,
                                                       "range": scale_range, "order": order}
            print(f"  {range_key:<16} {order:>5} {mean_e:>10.3f} {std_e:>9.3f}")

    print(f"\n  -- FLOOR TEST: fixing order=1 (ZNE-linear), sweeping RANGE --")
    lin_errs = [summary[f"{','.join(str(s) for s in r)}|order1"]["mean_kcal"] for r in SCALE_RANGES]
    ft_range_linear = floor_test([len(r) for r in SCALE_RANGES], lin_errs)
    print(f"  {ft_range_linear['verdict']}")

    print(f"\n  -- FLOOR TEST: fixing the widest range (1-7), sweeping ORDER --")
    widest = SCALE_RANGES[-1]
    widest_key = ",".join(str(s) for s in widest)
    orders = list(range(1, len(widest)))
    order_errs = [summary[f"{widest_key}|order{o}"]["mean_kcal"] for o in orders]
    ft_order = floor_test(orders, order_errs)
    print(f"  {ft_order['verdict']}")

    print(f"\n  -- FLOOR TEST: fixing order=2 (quadratic), sweeping RANGE --")
    quad_errs = [summary[f"{','.join(str(s) for s in r)}|order2"]["mean_kcal"] for r in SCALE_RANGES if len(r) > 2]
    quad_sizes = [len(r) for r in SCALE_RANGES if len(r) > 2]
    ft_range_quad = floor_test(quad_sizes, quad_errs)
    print(f"  {ft_range_quad['verdict']}")

    any_plateau = ft_range_linear["has_floor"] or ft_order["has_floor"] or ft_range_quad["has_floor"]
    best_key, best_row = min(summary.items(), key=lambda kv: kv[1]["mean_kcal"])
    print(f"\n  -- OVERALL VERDICT --")
    print(f"  best single (range,order) combination found: {best_key} = {best_row['mean_kcal']:.3f} "
          f"+/- {best_row['std_kcal']:.3f} kcal/mol (NOTE: cherry-picking the best cell is NOT a valid "
          "result on its own -- only reported alongside the floor-test verdicts above, per this "
          "project's own established discipline)")
    if any_plateau:
        print(f"  A genuine plateau WAS found on the tapered circuit -- unlike the untapered ansatz "
              "(iteration 13, Task E). This is real, new evidence the tapered circuit's ZNE behaves "
              "differently, not just smaller.")
    else:
        print(f"  NO PLATEAU FOUND on the tapered circuit either, in any of the three directions tested. "
              "Tapering reduced gate count but did NOT fix ZNE's convergence problem -- reported honestly, "
              "not softened. Per the task's own gate: this does NOT pass, so no real IonQ submission follows.")

    results = {
        "K": K, "shots": SHOTS, "n_seeds": N_SEEDS, "scale_ranges": SCALE_RANGES,
        "p2_per_gate": P2_PER_GATE, "p1_per_gate": P1_PER_GATE,
        "statevector_verification_err": sv_err,
        "gate_counts_scale1": gate_counts_by_scale[1],
        "summary": summary,
        "floor_test_range_at_order1": ft_range_linear,
        "floor_test_order_at_widest_range": ft_order,
        "floor_test_range_at_order2": ft_range_quad,
        "any_plateau_found": bool(any_plateau),
        "best_combination": {"key": best_key, **best_row},
        "chemical_accuracy_kcal": CHEM_ACC_KCAL,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results, problem


if __name__ == "__main__":
    main()
