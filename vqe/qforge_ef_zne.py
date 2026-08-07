#!/usr/bin/env python3
"""
qforge_ef_zne.py — Task 4: rebuild the original EF+ZNE result with this
project's current, verified machinery (vqe/qforge/), reporting what
changes and what doesn't.
============================================================================
The original 0.57 kcal/mol (entanglement_forging_zne.py) used K=5,
generic StatePreparation (no real gauge, no beta_signs() shortcut, 4-phase
cross terms, independently-measured beta register -- see Task 2's finding
that this costs real accuracy on real noise), and NO shot noise (Aer's
exact density-matrix estimator). This file rebuilds the SAME Quantinuum-
like-noise-rate ZNE experiment using: the fixed 11-gate particle-
conserving ansatz, real gauge (2 phase circuits/pair instead of 4),
beta_signs() (beta derived classically from the alpha measurement, not
independently measured), qubit-wise-commuting ("frame='h'") grouping, AND
shot noise (Binomial sampling on top of the exact noisy density-matrix
value, per shot_noise_study.py's verified methodology -- NOT re-deriving
that verification here, reusing it).

REPORTS BOTH K=5 AND K=6, honestly: K=5 sits on a classical TRUNCATION
floor (0.5655 kcal/mol -- this fragment's Schmidt rank is exactly 6, not
5, so K=5 can never reach exact agreement no matter how good the
measurement is). K=6 has NO such floor (Schmidt rank is exactly 6,
verified elsewhere in this project to ~1e-11 kcal/mol). Comparing "the
original K=5 result" to "a K=6 rebuild" without saying this would be
comparing two different ceilings, not just two different methods.

MANDATORY FLOOR TESTS: ZNE fit order (linear/quadratic/cubic -- does
error keep falling with no floor as more fit freedom is added, the same
disqualifying pattern iteration 2's training radius showed?) and the
noise-scale range (SCALES=[1,2,3] vs extending further).

Run:
    python vqe/qforge_ef_zne.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, derive_beta_matrices, shot_sample, floor_test,
    build_ansatz,
)
from qiskit import transpile
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.quantum_info import Pauli
import ef_fragment as effrag

QUANTINUUM_TWO_Q_ERROR = 0.002
QUANTINUUM_ONE_Q_ERROR = 0.00005
SCALES = [1.0, 2.0, 3.0]
BASIS_GATES = ["u3", "cx"]
N_SEEDS = 8
SHOTS = 100_000
K_VALUES = [5, 6]
CLASSICAL_FLOOR_K5_KCAL = 0.5655
CMAP4 = CouplingMap.from_full(4)
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "qforge_ef_zne_results.json")


def build_scaled_noise_model(scale):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(QUANTINUUM_TWO_Q_ERROR * scale, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(QUANTINUUM_ONE_Q_ERROR * scale, 1), "u3")
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
        vals = {}
        for l in non_id_labels:
            P = np.asarray(Pauli(l).to_matrix())
            vals[l] = float(np.real(np.trace(P @ dm)))
        raw[name] = vals
    return raw


def shot_noisy_energy(p, exact_raw, n_shots, rng, K):
    raw = {name: {l: shot_sample(v, n_shots, rng) for l, v in vals.items()} for name, vals in exact_raw.items()}
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def zne_fit(scales, energies, order):
    coeffs = np.polyfit(scales, energies, order)
    return float(np.polyval(coeffs, 0))


def run_for_K(K):
    print(f"\n  ==== K={K} ====")
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=(K == 6))
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    n_targets = len(p["targets"])
    print(f"  {n_ok}/{n_targets} targets converged, worst={worst:.2e}")
    counts = verify_constant_gate_count(solutions)
    print(f"  gate count across all targets: {counts}")
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]

    truncation_floor_kcal = abs(p["exact_energy"] - p["noiseless_energy"]) * HARTREE_TO_KCAL_MOL
    print(f"  classical truncation floor (K={K}): {truncation_floor_kcal:.4f} kcal/mol")

    # -- exact noisy matrices at each scale (shot-count-independent, computed once) --
    exact_raw_by_scale = {}
    for s in SCALES:
        nm = build_scaled_noise_model(s)
        exact_raw_by_scale[s] = measure_exact_noisy_raw(p, non_id_labels, nm)

    # -- 8-seed sweep: shot-noisy energy at each scale, ZNE fit (linear + quadratic + cubic) per seed --
    seed_results = {"raw": [], "zne_linear": [], "zne_quadratic": [], "zne_cubic": []}
    seed_results_noiseless = {"raw": [], "zne_linear": [], "zne_quadratic": [], "zne_cubic": []}
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed * 7919 + K)
        energies = []
        for s in SCALES:
            E, _ = shot_noisy_energy(p, exact_raw_by_scale[s], SHOTS, rng, K)
            energies.append(E)
        err_exact = abs(energies[0] - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
        err_noiseless = abs(energies[0] - p["noiseless_energy"]) * HARTREE_TO_KCAL_MOL
        seed_results["raw"].append(err_exact)
        seed_results_noiseless["raw"].append(err_noiseless)
        for order, key in ((1, "zne_linear"), (2, "zne_quadratic"), (3, "zne_cubic")):
            if len(SCALES) > order:
                E0 = zne_fit(SCALES, energies, order)
                seed_results[key].append(abs(E0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL)
                seed_results_noiseless[key].append(abs(E0 - p["noiseless_energy"]) * HARTREE_TO_KCAL_MOL)

    summary = {}
    for method in ("raw", "zne_linear", "zne_quadratic", "zne_cubic"):
        if seed_results[method]:
            summary[method] = {
                "err_vs_exact_mean_kcal": float(np.mean(seed_results[method])),
                "err_vs_exact_std_kcal": float(np.std(seed_results[method])),
                "err_vs_noiseless_mean_kcal": float(np.mean(seed_results_noiseless[method])),
                "err_vs_noiseless_std_kcal": float(np.std(seed_results_noiseless[method])),
            }
            print(f"    {method:>14}: err_vs_exact={summary[method]['err_vs_exact_mean_kcal']:.3f} +/- "
                  f"{summary[method]['err_vs_exact_std_kcal']:.3f} kcal/mol   "
                  f"err_vs_noiseless={summary[method]['err_vs_noiseless_mean_kcal']:.3f} +/- "
                  f"{summary[method]['err_vs_noiseless_std_kcal']:.3f} kcal/mol")

    # -- FLOOR TEST: does ZNE fit order improve without bound? --
    order_errs = [summary[k]["err_vs_exact_mean_kcal"] for k in ("zne_linear", "zne_quadratic", "zne_cubic") if k in summary]
    order_vals = [1, 2, 3][:len(order_errs)]
    fit_order_floor = floor_test(order_vals, order_errs) if len(order_errs) >= 2 else None
    if fit_order_floor:
        print(f"    FLOOR TEST (ZNE fit order): {fit_order_floor['verdict']}")

    # -- FLOOR TEST: noise-scale range -- extend past [1,2,3] to [1,2,3,4,5] --
    extended_scales = [1.0, 2.0, 3.0, 4.0, 5.0]
    exact_raw_extended = dict(exact_raw_by_scale)
    for s in extended_scales:
        if s not in exact_raw_extended:
            nm = build_scaled_noise_model(s)
            exact_raw_extended[s] = measure_exact_noisy_raw(p, non_id_labels, nm)
    range_errs = []
    for n_scales in (3, 4, 5):
        scales_here = extended_scales[:n_scales]
        rng = np.random.default_rng(9999 + K)
        energies = [shot_noisy_energy(p, exact_raw_extended[s], SHOTS, rng, K)[0] for s in scales_here]
        E0 = zne_fit(scales_here, energies, 2)
        range_errs.append(abs(E0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL)
    scale_range_floor = floor_test([3, 4, 5], range_errs)
    print(f"    FLOOR TEST (noise-scale range, quadratic fit): {scale_range_floor['verdict']}")

    return {
        "K": K, "n_targets_converged": f"{n_ok}/{n_targets}", "gate_count": list(counts),
        "exact_energy_ha": p["exact_energy"], "noiseless_energy_ha": p["noiseless_energy"],
        "truncation_floor_kcal": truncation_floor_kcal,
        "shots": SHOTS, "n_seeds": N_SEEDS, "scales": SCALES,
        "summary": summary,
        "floor_test_zne_fit_order": fit_order_floor,
        "floor_test_noise_scale_range": scale_range_floor,
    }


def main():
    print("\n" + "=" * 96)
    print("  qforge_ef_zne.py -- Task 4: EF+ZNE rebuilt with current machinery, K=5 and K=6")
    print("=" * 96)

    results = {}
    for K in K_VALUES:
        results[str(K)] = run_for_K(K)

    print(f"\n  -- SUMMARY: original (K=5, no shot noise, no real gauge) vs rebuilt --")
    print(f"    original K=5 ZNE-quadratic: 0.57 kcal/mol (NOT a real-hardware measurement -- "
          "emulator validation elsewhere in this project found the local model overestimates "
          "diagonal-term noise 0.61x and underestimates off-diagonal ~1.6x)")
    for K in K_VALUES:
        r = results[str(K)]
        zq = r["summary"].get("zne_quadratic", {})
        print(f"    rebuilt K={K} ZNE-quadratic (WITH shot noise, {SHOTS:,} shots/setting): "
              f"{zq.get('err_vs_exact_mean_kcal', float('nan')):.3f} +/- {zq.get('err_vs_exact_std_kcal', float('nan')):.3f} kcal/mol "
              f"(classical floor: {r['truncation_floor_kcal']:.4f} kcal/mol)")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
