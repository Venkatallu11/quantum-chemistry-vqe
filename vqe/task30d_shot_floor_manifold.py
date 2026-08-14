#!/usr/bin/env python3
"""
task30d_shot_floor_manifold.py -- iteration 30, Task D. REVISIT THE SHOT
FLOOR. The 1.845 kcal/mol "fundamental" shot floor was measured with the
OLD (raw per-Pauli) estimator. The manifold estimator (cleared, Task 30A)
took ideal fold-1 error from 1.456 to 0.132 kcal/mol -- re-measure the
ideal shot floor with the NEW estimator.
============================================================================
METHOD: exact statevector probabilities (already-verified equivalence to
real ionq_simulator sampling, Task 29A's own established technique) for
the optimizer-reduced (28B) circuit's ideal model, multinomial-sampled at
5 shot levels (10k/25k/50k/100k/300k), 8 seeds each, fed through the
manifold pure-state estimator (single-fold, no ZNE -- this iteration
freezes extrapolation). Fits sigma_E(N) = a/sqrt(N) + b via nonlinear
least squares and reports b explicitly -- if b < 0.25 kcal/mol, chemical
accuracy (0.5 kcal/mol) is reachable by shots alone with this estimator,
and the earlier "counting statistics make this impossible" conclusion was
estimator-specific, not fundamental. Reported either way, not assumed.

Run:
    python vqe/task30d_shot_floor_manifold.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_simulator_binding_curve import stable_seed, expectation_from_counts
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from qiskit.quantum_info import Statevector

K = 6
GATE_NAME = "zz"
SHOT_LEVELS = [10_000, 25_000, 50_000, 100_000, 300_000]
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task30d_shot_floor_manifold_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task30d_shot_floor_manifold.py -- ideal shot floor, manifold estimator, 5 shot levels")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    print("  building exact ideal probabilities for the optimizer-reduced circuit (all 21 kept slots)...")
    exact_probs = {}
    for name in kept:
        base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        exact_probs[name] = {}
        for group in groups:
            combined = effrag_mod.combined_basis_label(group)
            basis_qc = native_basis_change(combined, GATE_NAME)
            qc = base.compose(basis_qc)
            sv = Statevector.from_instruction(qc)
            exact_probs[name][tuple(group)] = sv.probabilities_dict()
    print("  done.")

    sigma_by_N = {}
    for n_shots in SHOT_LEVELS:
        errs = []
        for seed in range(N_SEEDS):
            manifold_kept = {}
            for name in kept:
                m_dict = {}
                w_dict = {}
                for group_t, probs in exact_probs[name].items():
                    bitstrings = list(probs.keys())
                    parr = np.array([probs[b] for b in bitstrings])
                    parr = parr / parr.sum()
                    rng = np.random.default_rng(stable_seed("task30d", n_shots, seed, name, group_t))
                    draws = rng.multinomial(n_shots, parr)
                    counts = {b: int(c) for b, c in zip(bitstrings, draws) if c > 0}
                    total = sum(counts.values())
                    for l in group_t:
                        m = expectation_from_counts(counts, l)
                        m_dict[l] = m
                        var = max(1 - m ** 2, 1e-4) / max(total, 1)
                        w_dict[l] = 1.0 / var
                v0 = target_coeff_vector(name, K)
                a_hat, _ = fit_pure_state(P_S, m_dict, w_dict, K, v0,
                                           seed=stable_seed("task30d_fit", n_shots, seed, name))
                manifold_kept[name] = a_hat

            full = build_full_from_a(manifold_kept, P_S, diag, K, non_id_labels)
            alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
            E, errs_d = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                     exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs.append(errs_d["err_vs_exact_kcal"])

        mean_e = float(np.mean(errs))
        std_e = float(np.std(errs))
        sigma_by_N[n_shots] = {"mean_abs_err": mean_e, "std_across_seeds": std_e, "errs": errs}
        print(f"    N={n_shots:>7,}: mean|err|={mean_e:.4f}  std={std_e:.4f} kcal/mol (n={N_SEEDS} seeds)")

    # -- fit sigma_E(N) = a/sqrt(N) + b to the STD (the shot-noise-driven spread), not the mean --
    Ns = np.array(SHOT_LEVELS, dtype=float)
    sigmas = np.array([sigma_by_N[n]["std_across_seeds"] for n in SHOT_LEVELS])

    def model(N, a, b):
        return a / np.sqrt(N) + b

    popt, pcov = curve_fit(model, Ns, sigmas, p0=[sigmas[0] * np.sqrt(Ns[0]), 0.01], bounds=([0, 0], [np.inf, np.inf]))
    a_fit, b_fit = popt
    print(f"\n  FIT: sigma_E(N) = {a_fit:.4f}/sqrt(N) + {b_fit:.4f}")
    print(f"  b (the floor as N -> infinity) = {b_fit:.4f} kcal/mol")

    reachable = b_fit < 0.25
    print(f"\n  VERDICT: b {'< ' if reachable else '>= '}0.25 kcal/mol -> "
          f"{'0.5 kcal/mol chemical accuracy IS reachable by shots alone with the manifold estimator -- the earlier counting-statistics floor was ESTIMATOR-SPECIFIC, not fundamental' if reachable else 'the floor is still too high for shots alone to close the gap -- consistent with a genuinely harder-than-shot-noise limitation'}")

    print(f"\n  cross-check vs the OLD estimator's cited 'fundamental' floor (1.845 kcal/mol, iteration 26 Task 1):")
    print(f"    manifold's own floor (b={b_fit:.4f}) is "
          f"{'far below' if b_fit < 1.845/3 else ('below' if b_fit < 1.845 else 'NOT below')} that number.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "sigma_by_N": sigma_by_N, "fit_a": float(a_fit), "fit_b": float(b_fit),
            "reachable_lt_0_25": bool(reachable), "old_estimator_floor_kcal": 1.845,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return sigma_by_N, a_fit, b_fit


if __name__ == "__main__":
    main()
