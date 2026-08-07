#!/usr/bin/env python3
"""
zne_floor_tested.py — Task E: push ZNE-linear properly. Requires a
PLATEAU before reporting any extrapolated number, sweeping BOTH the
noise-scale range and the fit order, trying Richardson (polynomial
interpolation) and exponential fits, and reporting the extrapolation's
own uncertainty (8-seed, shot noise included) rather than a bare central
value.
============================================================================
WHY THIS FILE EXISTS: iteration 11's Task 4 found the classic 0.57
kcal/mol EF+ZNE result fails its OWN floor test -- extending the noise-
scale range from [1,2,3] to [1,2,3,4,5] changed the quadratic-fit answer
by 34x (K=5) / 5.4x (K=6) with no plateau. ZNE-linear was nonetheless the
best REAL-hardware method found (iteration 12, ~30 kcal/mol on real IonQ
noise). This file asks the mandatory follow-up question properly: is
there ANY (scale range, fit order) combination that DOES plateau for
this circuit, or is ZNE simply not converged here at any setting tested,
full stop?

METHOD: for each scale range (a set of noise-scale factors starting at
1), fit EVERY polynomial order from 1 up to len(range)-1 (order =
len(range)-1 is exact polynomial INTERPOLATION through every point --
what "Richardson extrapolation" means for ZNE, mathematically identical
to a full-degree least-squares fit through the same points) via
np.polyfit, plus a separate exponential fit E(s) = A*exp(-k*s) + E_inf
(linearized via log|E(s)-E(s_max)|, only well-posed when the shot-noisy
energies are monotonic in scale -- reported as unavailable, not faked,
when they are not). A genuine plateau requires the extrapolated error to
stop changing as EITHER the range grows (more, farther-out scale points)
OR the order grows (more fit freedom) -- both checked via
qforge.floor_test, not eyeballed.

Shot noise is included in every headline number (100,000 shots/setting,
8 real independent seeds, Binomial sampling on top of the exact noisy
density-matrix value -- shot_noise_study.py's verified methodology,
reused not re-derived).

Run:
    python vqe/zne_floor_tested.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, shot_sample, floor_test, build_ansatz,
    P2_PER_GATE, P1_PER_GATE,
)
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.quantum_info import Pauli

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
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "zne_floor_tested_results.json")


def build_noise_model(scale):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(min(P2_PER_GATE * scale, 0.75), 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(min(P1_PER_GATE * scale, 0.75), 1), "u3")
    return nm


def noisy_density_matrix(angles, noise_model):
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    qc2 = qc.copy()
    qc2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc2).result()
    return np.asarray(result.data(0)["density_matrix"])


def measure_exact_noisy_raw(p, non_id_labels, noise_model):
    raw = {}
    for name, sol in p["solutions"].items():
        dm = noisy_density_matrix(sol["angles"], noise_model)
        raw[name] = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm))) for l in non_id_labels}
    return raw


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
    return None, "energies not monotonic across scales -- exponential fit not well-posed, not faked"


def main():
    print("\n" + "=" * 96)
    print("  zne_floor_tested.py -- Task E: ZNE with a mandatory plateau requirement")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    assert n_ok == 36
    counts = verify_constant_gate_count(solutions)
    assert counts == {11}
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    print(f"  setup OK: 36/36 converged, gate count={counts}, P2_PER_GATE={P2_PER_GATE}")

    print(f"\n  measuring exact noisy matrices at scales 1..{MAX_SCALE} (shot-count independent, cached once)")
    exact_raw_by_scale = {}
    for s in range(1, MAX_SCALE + 1):
        nm = build_noise_model(s)
        exact_raw_by_scale[s] = measure_exact_noisy_raw(p, non_id_labels, nm)
        E, err = energy_and_err(p, exact_raw_by_scale[s], K)
        print(f"    scale={s}: exact (no shot noise) err_vs_exact={err:.3f} kcal/mol")

    print(f"\n  -- 8-seed shot-noisy sweep across every (scale range, fit order) combination --")
    range_order_results = {}  # (range_key, order) -> list of 8 per-seed errors
    range_exp_results = {}    # range_key -> list of 8 per-seed errors (or None)

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
            E0_exp, note = exponential_fit(scale_range, energies)
            if E0_exp is not None:
                range_exp_results[range_key].append(abs(E0_exp - p["exact_energy"]) * HARTREE_TO_KCAL_MOL)

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
        if range_exp_results[range_key]:
            errs = range_exp_results[range_key]
            mean_e, std_e = float(np.mean(errs)), float(np.std(errs))
            summary[f"{range_key}|exp"] = {"mean_kcal": mean_e, "std_kcal": std_e, "range": scale_range, "order": "exp"}
            print(f"  {range_key:<16} {'exp':>5} {mean_e:>10.3f} {std_e:>9.3f}")
        else:
            print(f"  {range_key:<16} {'exp':>5} {'n/a':>10} (non-monotonic, not faked)")

    # -- FLOOR TEST 1: fixed order=1 (linear, ZNE-linear's own regime), does extending the RANGE plateau? --
    print(f"\n  -- FLOOR TEST: fixing order=1 (ZNE-linear), sweeping RANGE --")
    range_sizes = [len(r) for r in SCALE_RANGES]
    lin_errs = [summary[f"{','.join(str(s) for s in r)}|order1"]["mean_kcal"] for r in SCALE_RANGES]
    ft_range_linear = floor_test(range_sizes, lin_errs)
    print(f"  {ft_range_linear['verdict']}")

    # -- FLOOR TEST 2: fixed range=[1..7], does increasing ORDER plateau? --
    print(f"\n  -- FLOOR TEST: fixing the widest range (1-7), sweeping ORDER --")
    widest = SCALE_RANGES[-1]
    widest_key = ",".join(str(s) for s in widest)
    orders = list(range(1, len(widest)))
    order_errs = [summary[f"{widest_key}|order{o}"]["mean_kcal"] for o in orders]
    ft_order = floor_test(orders, order_errs)
    print(f"  {ft_order['verdict']}")

    # -- FLOOR TEST 3: fixed order=2 (quadratic, iteration 11's own disqualified case), sweeping RANGE --
    print(f"\n  -- FLOOR TEST: fixing order=2 (quadratic), sweeping RANGE (reproducing iteration 11's check) --")
    quad_errs = [summary[f"{','.join(str(s) for s in r)}|order2"]["mean_kcal"] for r in SCALE_RANGES if len(r) > 2]
    quad_sizes = [len(r) for r in SCALE_RANGES if len(r) > 2]
    ft_range_quad = floor_test(quad_sizes, quad_errs)
    print(f"  {ft_range_quad['verdict']}")

    any_plateau = ft_range_linear["has_floor"] or ft_order["has_floor"] or ft_range_quad["has_floor"]
    print(f"\n  -- OVERALL VERDICT --")
    if any_plateau:
        print(f"  A genuine plateau WAS found in at least one (range, order) sweep direction -- "
              f"see which one above before trusting any specific number.")
    else:
        print(f"  NO PLATEAU FOUND in any of the three sweep directions tested (range at order=1, "
              f"order at the widest range, range at order=2). Per this project's own disqualification "
              f"rule: ZNE is NOT demonstrated to be converged for this circuit at this (local, synthetic) "
              f"noise level, full stop -- this is reported as the honest result, not softened.")

    results = {
        "K": K, "shots": SHOTS, "n_seeds": N_SEEDS, "scale_ranges": SCALE_RANGES,
        "p2_per_gate": P2_PER_GATE, "p1_per_gate": P1_PER_GATE,
        "summary": summary,
        "floor_test_range_at_order1": ft_range_linear,
        "floor_test_order_at_widest_range": ft_order,
        "floor_test_range_at_order2": ft_range_quad,
        "any_plateau_found": bool(any_plateau),
        "chemical_accuracy_kcal": 1.0,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
