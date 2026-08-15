#!/usr/bin/env python3
"""
task31g_control_variate.py -- iteration 31, Task G. PAIRED-REFERENCE
CONTROL VARIATE. Estimate beta* = Cov(Y,X)/Var(X), report
Var(Y_cv) = Var(Y)(1-rho^2). Bringing 2.31 kcal/mol submission-to-
submission variability down to 0.25 requires rho > 0.994.
============================================================================
DESIGN, no new real submission needed: Task F's 8 independent real
submissions ALREADY measured all 21 H4 slots together, in the same job,
every time -- exactly the "interleaved reference + target within one
block" structure this task calls for, without needing a dedicated new
circuit. Uses `u_0` (a diagonal Schmidt-basis slot, exact target trivially
known: e_0) as the reference X -- its own measured deviation from its
known-exact Pauli values, real Hamiltonian-coefficient-weighted, exactly
analogous to what a dedicated reference circuit R would measure. Target Y
is Task F's own recorded full-energy error E_i for that same submission.

CAVEAT, disclosed prominently: n=8 data points is a small sample for
estimating a correlation coefficient reliably (same statistical-power
caveat already raised for Task F's autocorrelation) -- reported as such,
not overclaimed.

Run:
    python vqe/task31g_control_variate.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS = 8
N_SUBMISSIONS = 8
REFERENCE_SLOT = "u_0"
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CONV_RESULTS = os.path.join(os.path.dirname(__file__), "task31f_convergence_study_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31g_control_variate_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task31g_control_variate.py -- paired-reference control variate, reusing Task F's real data")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)

    # -- exact target values for u_0 (e_0 in the Schmidt basis: <P_l> = P_S[l][0,0]) --
    from phys_constrained_reconstruction import build_P_S
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    exact_u0 = {l: float(np.real(P_S[l][0, 0])) for l in non_id_labels}

    # -- label energy weights, for a Hamiltonian-coefficient-weighted reference discrepancy --
    label_weight = {}
    for (a_label, b_label, coeff) in p["terms"]:
        w = abs(coeff) * HARTREE_TO_KCAL_MOL
        label_weight[a_label] = label_weight.get(a_label, 0.0) + w

    with open(CONV_RESULTS) as f:
        conv = json.load(f)
    Y = np.array(conv["E_list"])  # Task F's own recorded full-energy errors, real, 8 submissions

    X = []
    for sub in range(N_SUBMISSIONS):
        ckpt_path = os.path.join(CKPT_DIR, f"task31f_convergence_submission_{sub}.json")
        with open(ckpt_path) as f:
            ck = json.load(f)
        tags = ck["tags"]
        counts_list = ck["counts"]
        per_group = {}
        for (name, group), counts in zip(tags, counts_list):
            if name == REFERENCE_SLOT:
                per_group[tuple(group)] = counts

        weighted_dev = 0.0
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task31g", sub, seed))
            seed_dev = 0.0
            for group_t, counts in per_group.items():
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    m = expectation_from_counts(resampled, l)
                    seed_dev += abs(m - exact_u0[l]) * label_weight.get(l, 0.0)
            weighted_dev += seed_dev / N_SEEDS
        X.append(weighted_dev)
        print(f"  submission {sub+1}/{N_SUBMISSIONS}: X (reference discrepancy, {REFERENCE_SLOT}) = "
              f"{weighted_dev:.4f}  Y (target energy error) = {Y[sub]:.4f}")

    X = np.array(X)
    cov_xy = float(np.cov(X, Y, ddof=1)[0, 1])
    var_x = float(np.var(X, ddof=1))
    var_y = float(np.var(Y, ddof=1))
    beta_star = cov_xy / var_x if var_x > 1e-12 else 0.0
    rho = cov_xy / np.sqrt(var_x * var_y) if var_x > 1e-12 and var_y > 1e-12 else 0.0
    var_y_cv = var_y * (1 - rho ** 2)

    print(f"\n  Cov(Y,X)={cov_xy:.6f}  Var(X)={var_x:.6f}  Var(Y)={var_y:.6f}")
    print(f"  beta* = {beta_star:.4f}")
    print(f"  rho = {rho:.4f}  (n={N_SUBMISSIONS} -- small-sample correlation estimate, real uncertainty)")
    print(f"  Var(Y_cv) = Var(Y)*(1-rho^2) = {var_y_cv:.6f}  (vs Var(Y)={var_y:.6f})")
    print(f"  std(Y_cv) = {np.sqrt(var_y_cv):.4f}  vs std(Y) = {np.sqrt(var_y):.4f}")

    rho_needed = np.sqrt(1 - (0.25 / 2.31) ** 2)
    print(f"\n  rho needed to bring 2.31 -> 0.25: {rho_needed:.4f}")
    print(f"  ACHIEVED rho = {rho:.4f}  -- {'MEETS' if abs(rho) >= rho_needed else 'DOES NOT MEET'} the bar")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "X": X.tolist(), "Y": Y.tolist(), "beta_star": beta_star, "rho": rho,
            "var_y": var_y, "var_y_cv": var_y_cv, "rho_needed_for_target": float(rho_needed),
            "meets_bar": bool(abs(rho) >= rho_needed),
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
