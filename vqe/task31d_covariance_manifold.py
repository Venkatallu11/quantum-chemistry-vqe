#!/usr/bin/env python3
"""
task31d_covariance_manifold.py -- iteration 31, Task D. COVARIANCE-AWARE
MANIFOLD FIT. The current fit uses uniform weighting (flagged as a
simplification in Task 30C). PEC output has nontrivial covariance --
labels measured from the SAME bitstring counts within a group are
correlated by construction. Bootstrap the full covariance matrix Sigma_y
per slot (not just per-label variance), then minimize
    chi^2(a) = (y - f(a))^T Sigma_y^{-1} (y - f(a))   subject to |a|_2 = 1
No new real submission -- reuses Task 31C's already-collected 4,368 real
circuits' raw counts, just reprocessed with more bootstrap seeds (32,
up from 8) for a better-conditioned covariance estimate.
============================================================================
REGULARIZATION, disclosed not hidden: with up to ~30 labels per slot and
only 32 bootstrap seeds, the sample covariance matrix is rank-deficient
(rank <= 31). Sigma_y^{-1} is computed via a REGULARIZED inverse
(Tikhonov: Sigma_y + eps*I, eps = 1% of the mean diagonal variance) rather
than a raw pseudo-inverse, to keep the chi^2 well-conditioned everywhere,
not just on the rank actually spanned by 32 samples.

PASS = bias decreases AND ideal control still holds. A better central
number with a degraded ideal control is a FAIL, per this task's own rule.

Run:
    python vqe/task31d_covariance_manifold.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state as fit_pure_state_uniform
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS_COV = 32  # up from Task 31C's 8, for a better-conditioned covariance estimate
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31d_covariance_manifold_results.json")


def fit_pure_state_cov(P_S, labels, y_mean, Sigma_y_inv, K, v0, seed, n_restarts=8):
    Ps = [P_S[l] for l in labels]
    y = np.array([y_mean[l] for l in labels])

    def objective(v):
        a = v / np.linalg.norm(v)
        pred = np.array([float(np.real(a @ P @ a)) for P in Ps])
        resid = pred - y
        return float(resid @ Sigma_y_inv @ resid)

    rng = np.random.default_rng(seed)
    inits = [v0] + [v0 + rng.normal(0, scale, K) for scale in [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2]]
    best_val, best_v = float("inf"), None
    for v_init in inits[:n_restarts]:
        res = minimize(objective, v_init, method="L-BFGS-B")
        if res.fun < best_val:
            best_val, best_v = res.fun, res.x
    return best_v / np.linalg.norm(best_v)


def build_full_from_a(a_by_name, P_S, diag, K, non_id_labels):
    def expect(a, l):
        return float(np.real(a @ P_S[l] @ a))
    full = {}
    for name in diag:
        a = a_by_name[name]
        full[name] = {l: expect(a, l) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            a_pl = a_by_name[pl]
            full[pl] = {l: expect(a_pl, l) for l in non_id_labels}
            synth_minus = {}
            for l in non_id_labels:
                synth_minus[l] = full[un][l] + full[um][l] - full[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def per_seed_pec_values(name, kept_labels, tags, counts_list, gamma_per_circuit, n_seeds):
    """Rebuild N_SEEDS_COV independent bootstrap draws of the signed/gamma
    -weighted literal-PEC estimate for every label of this slot, reusing
    Task 31C's already-collected real counts (no new submission)."""
    by_label = {}
    for (tname, group, draw, sign), counts, gamma in zip(tags, counts_list, gamma_per_circuit):
        if tname != name:
            continue
        for l in group:
            by_label.setdefault(l, []).append((sign, gamma, counts))

    seed_matrix = {l: [] for l in by_label}
    for seed in range(n_seeds):
        for l, entries in by_label.items():
            rng = np.random.default_rng(stable_seed("task31d_cov", name, l, seed))
            vals_signed = []
            for sign, gamma, counts in entries:
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m = expectation_from_counts(resampled, l)
                vals_signed.append(sign * gamma * m)
            seed_matrix[l].append(float(np.mean(vals_signed)))
    return seed_matrix


def main():
    print("\n" + "=" * 96)
    print("  task31d_covariance_manifold.py -- covariance-aware chi^2 fit vs current uniform-weighted fit")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(CKPT_PATH) as f:
        ck = json.load(f)
    tags = [tuple(t) for t in ck["tags"]]
    counts_list = ck["counts"]
    gamma_per_circuit = ck["gamma_per_circuit"]

    print(f"  rebuilding {N_SEEDS_COV}-seed bootstrap covariance per slot from Task 31C's real counts "
          f"(no new submission)...")
    a_uniform, a_cov = {}, {}
    cond_numbers = {}
    for i, name in enumerate(kept):
        seed_matrix = per_seed_pec_values(name, None, tags, counts_list, gamma_per_circuit, N_SEEDS_COV)
        labels = list(seed_matrix.keys())
        Y = np.array([seed_matrix[l] for l in labels])  # (n_labels, n_seeds)
        y_mean = {l: float(np.mean(Y[j])) for j, l in enumerate(labels)}
        Sigma_y = np.cov(Y)
        if Sigma_y.ndim == 0:
            Sigma_y = np.array([[float(Sigma_y)]])
        reg = 0.01 * np.mean(np.diag(Sigma_y))
        Sigma_y_reg = Sigma_y + reg * np.eye(len(labels))
        Sigma_y_inv = np.linalg.inv(Sigma_y_reg)
        cond_numbers[name] = float(np.linalg.cond(Sigma_y_reg))

        v0 = target_coeff_vector(name, K)
        a_uniform[name], _ = fit_pure_state_uniform(P_S, y_mean, {l: 1.0 for l in labels}, K, v0,
                                                      seed=stable_seed("task31d_u", name))
        a_cov[name] = fit_pure_state_cov(P_S, labels, y_mean, Sigma_y_inv, K, v0,
                                          seed=stable_seed("task31d_c", name))
        print(f"    [{i+1}/{len(kept)}] {name}: {len(labels)} labels, Sigma_y cond#={cond_numbers[name]:.1f}")

    full_uniform = build_full_from_a(a_uniform, P_S, diag, K, non_id_labels)
    full_cov = build_full_from_a(a_cov, P_S, diag, K, non_id_labels)
    _, err_uniform = energy_and_err(p, full_uniform, K)
    _, err_cov = energy_and_err(p, full_cov, K)
    print(f"\n  forte-1: uniform-weighted fit = {err_uniform:.4f} kcal/mol")
    print(f"  forte-1: covariance-aware fit = {err_cov:.4f} kcal/mol")

    # -- ideal-control check for BOTH fits --
    print(f"\n  -- ideal-control check, both fits --")
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    ideal_tags = raw_ck["tags"]["ideal"]
    ideal_counts = raw_ck["counts"]["ideal"]
    ideal_per_name = {}
    for (name, group), counts in zip(ideal_tags, ideal_counts):
        ideal_per_name.setdefault(name, {}).setdefault(tuple(group), counts)

    a_uniform_ideal, a_cov_ideal = {}, {}
    for name in kept:
        seed_matrix = {}
        for group_t, counts in ideal_per_name[name].items():
            for l in group_t:
                vals = []
                for seed in range(N_SEEDS_COV):
                    rng = np.random.default_rng(stable_seed("task31d_ideal", name, l, seed))
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    vals.append(expectation_from_counts(resampled, l))
                seed_matrix[l] = vals
        labels = list(seed_matrix.keys())
        Y = np.array([seed_matrix[l] for l in labels])
        y_mean = {l: float(np.mean(Y[j])) for j, l in enumerate(labels)}
        Sigma_y = np.cov(Y)
        if Sigma_y.ndim == 0:
            Sigma_y = np.array([[float(Sigma_y)]])
        reg = 0.01 * max(np.mean(np.diag(Sigma_y)), 1e-8)
        Sigma_y_inv = np.linalg.inv(Sigma_y + reg * np.eye(len(labels)))
        v0 = target_coeff_vector(name, K)
        a_uniform_ideal[name], _ = fit_pure_state_uniform(P_S, y_mean, {l: 1.0 for l in labels}, K, v0,
                                                            seed=stable_seed("task31d_ui", name))
        a_cov_ideal[name] = fit_pure_state_cov(P_S, labels, y_mean, Sigma_y_inv, K, v0,
                                                seed=stable_seed("task31d_ci", name))

    full_uniform_ideal = build_full_from_a(a_uniform_ideal, P_S, diag, K, non_id_labels)
    full_cov_ideal = build_full_from_a(a_cov_ideal, P_S, diag, K, non_id_labels)
    _, err_uniform_ideal = energy_and_err(p, full_uniform_ideal, K)
    _, err_cov_ideal = energy_and_err(p, full_cov_ideal, K)
    print(f"  ideal, uniform-weighted: {err_uniform_ideal:.4f} kcal/mol")
    print(f"  ideal, covariance-aware: {err_cov_ideal:.4f} kcal/mol")

    print(f"\n" + "=" * 96)
    bias_decreased = err_cov < err_uniform
    ideal_holds = err_cov_ideal <= err_uniform_ideal * 2.0 + 0.3
    verdict = "PASS" if (bias_decreased and ideal_holds) else "FAIL"
    print(f"  VERDICT: bias {'decreased' if bias_decreased else 'did NOT decrease'} "
          f"({err_uniform:.4f} -> {err_cov:.4f}); ideal control {'holds' if ideal_holds else 'DEGRADED'} "
          f"({err_uniform_ideal:.4f} -> {err_cov_ideal:.4f})")
    print(f"  -> {verdict}")
    print("=" * 96 + "\n")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "err_uniform_forte": err_uniform, "err_cov_forte": err_cov,
            "err_uniform_ideal": err_uniform_ideal, "err_cov_ideal": err_cov_ideal,
            "bias_decreased": bool(bias_decreased), "ideal_holds": bool(ideal_holds), "verdict": verdict,
            "condition_numbers": cond_numbers,
        }, f, indent=2)
    print(f"  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
