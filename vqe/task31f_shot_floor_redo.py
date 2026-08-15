#!/usr/bin/env python3
"""
task31f_shot_floor_redo.py -- iteration 31, Task F (shot-floor half).
Task 30D's shot floor (b=0.0524 kcal/mol) used only 8 seeds and showed
non-monotonic scatter (10k -> 0.312+/-0.179, 25k -> 0.357+/-0.295, 25k
WORSE than 10k) -- under-determined, not established. Re-run with 4x more
seeds (32) and report b WITH a bootstrap confidence interval, not a bare
point estimate.

Entirely local (exact ideal probabilities, multinomial-sampled -- the
established real-hardware-equivalent technique), no new real submission.

Run:
    python vqe/task31f_shot_floor_redo.py
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
N_SEEDS = 32  # 4x Task 30D's 8
N_BOOTSTRAP_CI = 2000
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31f_shot_floor_redo_results.json")


def fit_model(Ns, sigmas):
    def model(N, a, b):
        return a / np.sqrt(N) + b
    popt, _ = curve_fit(model, Ns, sigmas, p0=[sigmas[0] * np.sqrt(Ns[0]), 0.01],
                         bounds=([0, 0], [np.inf, np.inf]), maxfev=5000)
    return popt


def main():
    print("\n" + "=" * 96)
    print(f"  task31f_shot_floor_redo.py -- ideal shot floor, {N_SEEDS} seeds (4x Task 30D), with CI")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    print("  building exact ideal probabilities for all 21 slots...")
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
    print("  done.\n")

    sigma_by_N = {}
    for n_shots in SHOT_LEVELS:
        errs = []
        for seed in range(N_SEEDS):
            manifold_kept = {}
            for name in kept:
                m_dict, w_dict = {}, {}
                for group_t, probs in exact_probs[name].items():
                    bitstrings = list(probs.keys())
                    parr = np.array([probs[b] for b in bitstrings])
                    parr = parr / parr.sum()
                    rng = np.random.default_rng(stable_seed("task31f", n_shots, seed, name, group_t))
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
                                           seed=stable_seed("task31f_fit", n_shots, seed, name))
                manifold_kept[name] = a_hat

            full = build_full_from_a(manifold_kept, P_S, diag, K, non_id_labels)
            alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
            E, errs_d = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                     exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs.append(errs_d["err_vs_exact_kcal"])

        mean_e, std_e = float(np.mean(errs)), float(np.std(errs))
        sigma_by_N[n_shots] = {"mean_abs_err": mean_e, "std_across_seeds": std_e, "errs": errs}
        print(f"    N={n_shots:>7,}: mean|err|={mean_e:.4f}  std={std_e:.4f} kcal/mol (n={N_SEEDS} seeds)")

    Ns = np.array(SHOT_LEVELS, dtype=float)
    sigmas = np.array([sigma_by_N[n]["std_across_seeds"] for n in SHOT_LEVELS])
    a_fit, b_fit = fit_model(Ns, sigmas)
    print(f"\n  POINT FIT: sigma_E(N) = {a_fit:.4f}/sqrt(N) + {b_fit:.4f}")

    # -- bootstrap CI on b: resample the N_SEEDS errors at EACH shot level with replacement, refit --
    print(f"\n  -- bootstrapping a 95% CI on b ({N_BOOTSTRAP_CI} resamples) --")
    rng_boot = np.random.default_rng(2031)
    b_samples = []
    all_errs = {n: np.array(sigma_by_N[n]["errs"]) for n in SHOT_LEVELS}
    for _ in range(N_BOOTSTRAP_CI):
        sigmas_boot = []
        for n in SHOT_LEVELS:
            resampled = rng_boot.choice(all_errs[n], size=len(all_errs[n]), replace=True)
            sigmas_boot.append(float(np.std(resampled)))
        try:
            _, b_boot = fit_model(Ns, np.array(sigmas_boot))
            b_samples.append(b_boot)
        except Exception:
            continue
    b_samples = np.array(b_samples)
    ci_lo, ci_hi = float(np.percentile(b_samples, 2.5)), float(np.percentile(b_samples, 97.5))
    print(f"  b = {b_fit:.4f}  95% CI = [{ci_lo:.4f}, {ci_hi:.4f}]  (from {len(b_samples)}/{N_BOOTSTRAP_CI} successful bootstrap fits)")

    reachable = ci_hi < 0.25
    print(f"\n  VERDICT: 95% CI upper bound {'< ' if reachable else '>= '}0.25 kcal/mol -> "
          f"{'0.5 kcal/mol chemical accuracy IS reachable by shots alone, with confidence' if reachable else 'NOT established with confidence at the 0.25 threshold'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "sigma_by_N": sigma_by_N, "fit_a": float(a_fit), "fit_b": float(b_fit),
            "b_ci95_lo": ci_lo, "b_ci95_hi": ci_hi, "n_seeds": N_SEEDS,
            "reachable_with_confidence": bool(reachable),
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
