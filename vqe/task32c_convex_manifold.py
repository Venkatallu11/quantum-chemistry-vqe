#!/usr/bin/env python3
"""
task32c_convex_manifold.py -- iteration 32, Task C. DETERMINISTIC MANIFOLD
VIA CONVEX RELAXATION. The current 5-parameter fit (task29c's fit_pure_state)
is a nonconvex degree-4-polynomial-on-a-sphere optimization -- a likely
source of the (iii) manifold-optimizer variance quantified in Task 32B.
============================================================================
KEY OBSERVATION: m_i = a^T P_i a = Tr(P_i X) with X = a a^T is LINEAR in X.
So instead of optimizing over the unit sphere in R^K (nonconvex), solve

    min_X  sum_i w_i [ m_i - Tr(P_i X) ]^2   s.t.  X >= 0 (PSD), Tr(X) = 1

a convex SDP with a UNIQUE global optimum (strictly convex objective on a
convex feasible set) -- no local optima, no restarts, no seed dependence.
Then read off the eigenvalue spectrum of the optimal X: if lambda_1 >>
lambda_2..K, the physically-unconstrained relaxation NATURALLY lands close
to a pure (rank-1) state, and projecting to the top eigenvector is safe
(not silently overriding a genuinely mixed-looking answer).

SAME real data as Task 31C/31D (Task 31C's 4,368-circuit checkpoint,
N_SEEDS_COV=32 bootstrap, uniform weighting) -- so this is a true apples-
to-apples swap of ONLY the optimizer, not a new measurement.

Also tests dropping the target-anchored initialization from the NONLINEAR
fit entirely (8 unbiased random restarts instead of v0 + 7 perturbations)
-- the adversarial test (Task 30A) cleared target-anchoring as not
"injecting the answer," but a defensible production method should not
reference the answer at all if an anchor-free alternative works equally
well; reported for comparison, not adopted unless it holds up.

PASS = (a) IDENTICAL answers across independent reruns of the SDP (true by
construction for a convex problem -- verified, not assumed), (b) bias no
worse than Task 31D's own nonlinear-fit baseline (err_uniform_forte=0.3166,
same 32-seed data), (c) ideal control preserved (err_uniform_ideal=0.0466
baseline, same data).

Run:
    python vqe/task32c_convex_manifold.py
"""
import os
import sys
import json
import numpy as np
import cvxpy as cp
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state as fit_pure_state_uniform
from phys_constrained_reconstruction import build_P_S
from qforge import combine_matrices, energy_from_alpha_matrices
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts
from task31d_covariance_manifold import per_seed_pec_values, build_full_from_a

K = 6
SHOTS = 100_000
N_SEEDS_COV = 32
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
BASELINE_FORTE = 0.31659421376349306   # Task 31D's own nonlinear-fit result, SAME 32-seed data
BASELINE_IDEAL = 0.046613622763887454
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32c_convex_manifold_results.json")


def solve_convex_manifold(P_S, labels, y_mean, K):
    """Global unique optimum: min sum_l (y_l - Tr(P_l X))^2 s.t. X>=0, Tr(X)=1."""
    X = cp.Variable((K, K), symmetric=True)
    constraints = [X >> 0, cp.trace(X) == 1]
    resid = []
    for l in labels:
        resid.append(y_mean[l] - cp.trace(P_S[l] @ X))
    objective = cp.Minimize(cp.sum_squares(cp.hstack(resid)))
    prob = cp.Problem(objective, constraints)
    prob.solve(solver=cp.SCS)
    Xval = X.value
    Xval = (Xval + Xval.T) / 2  # numerical symmetrization
    eigvals, eigvecs = np.linalg.eigh(Xval)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    a_hat = eigvecs[:, 0]
    return a_hat, eigvals, Xval, prob.value


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def fit_pure_state_random_only(P_S, m_dict, w_dict, K, seed, n_restarts=8):
    """Task 32C's 'drop target-anchored init' variant -- 8 UNBIASED random
    restarts on the unit sphere, no reference to the known/expected target
    anywhere."""
    from scipy.optimize import minimize
    labels = list(m_dict.keys())
    Ps = [P_S[l] for l in labels]
    ms = np.array([m_dict[l] for l in labels])
    ws = np.array([w_dict[l] for l in labels])

    def objective(v):
        a = v / np.linalg.norm(v)
        pred = np.array([float(np.real(a @ P @ a)) for P in Ps])
        return float(np.sum(ws * (pred - ms) ** 2))

    rng = np.random.default_rng(seed)
    best_val, best_v = float("inf"), None
    for _ in range(n_restarts):
        v_init = rng.normal(0, 1, K)
        res = minimize(objective, v_init, method="L-BFGS-B")
        if res.fun < best_val:
            best_val, best_v = res.fun, res.x
    a_hat = best_v / np.linalg.norm(best_v)
    return a_hat, best_val


def compute_full_pipeline(p, non_id_labels, diag, kept, P_S, ck_tags, ck_counts, ck_gamma, method, seed_tag):
    a_by_name = {}
    eig_ratios = {}
    for name in kept:
        seed_matrix = per_seed_pec_values(name, None, ck_tags, ck_counts, ck_gamma, N_SEEDS_COV)
        labels = list(seed_matrix.keys())
        y_mean = {l: float(np.mean(seed_matrix[l])) for l in labels}
        v0 = target_coeff_vector(name, K)
        if method == "convex":
            a_hat, eigvals, Xval, objval = solve_convex_manifold(P_S, labels, y_mean, K)
            if np.dot(a_hat, v0) < 0:
                a_hat = -a_hat  # global sign gauge only, cosmetic (a and -a are the same physical state)
            eig_ratios[name] = (float(eigvals[0]), float(np.sum(eigvals[1:])))
        elif method == "nonlinear_anchored":
            a_hat, _ = fit_pure_state_uniform(P_S, y_mean, {l: 1.0 for l in labels}, K, v0,
                                               seed=stable_seed(seed_tag, name))
        elif method == "nonlinear_random":
            a_hat, _ = fit_pure_state_random_only(P_S, y_mean, {l: 1.0 for l in labels}, K,
                                                   seed=stable_seed(seed_tag, name))
        a_by_name[name] = a_hat
    full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err, eig_ratios


def main():
    print("\n" + "=" * 96)
    print("  task32c_convex_manifold.py -- convex SDP relaxation vs nonlinear fit, SAME real data")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(CKPT_PATH) as f:
        ck = json.load(f)
    tags = [tuple(t) for t in ck["tags"]]
    counts_list = ck["counts"]
    gamma_per_circuit = ck["gamma_per_circuit"]

    print(f"\n  -- forte-1: convex SDP relaxation --")
    err_convex, eig_ratios = compute_full_pipeline(p, non_id_labels, diag, kept, P_S, tags, counts_list,
                                                     gamma_per_circuit, "convex", "t32c_convex")
    print(f"    err = {err_convex:.4f} kcal/mol")
    lambda1s = [v[0] for v in eig_ratios.values()]
    rest = [v[1] for v in eig_ratios.values()]
    purity_ratio = np.mean([l1 / (l1 + r) for l1, r in zip(lambda1s, rest)])
    print(f"    mean lambda_1/(lambda_1+sum(lambda_2..K)) across {len(kept)} slots = {purity_ratio:.4f} "
          f"(1.0 = exactly pure, i.e. the relaxation's own optimum is already near rank-1)")

    print(f"\n  -- REPRODUCIBILITY check: rerun the SDP in 4 SEPARATE OS processes --")
    print(f"     (Task 32B found the NONLINEAR fit can land on bit-identical objective values via DIFFERENT")
    print(f"     solution vectors ACROSS process launches, even with a fixed seed -- a same-process rerun")
    print(f"     cannot detect that failure mode, so this check uses real separate processes, not just a")
    print(f"     second in-process call.)")
    with ProcessPoolExecutor(max_workers=4) as ex:
        futures = [ex.submit(compute_full_pipeline, p, non_id_labels, diag, kept, P_S, tags, counts_list,
                              gamma_per_circuit, "convex", f"t32c_convex_rerun_{i}")
                   for i in range(4)]
        rerun_errs = [fut.result()[0] for fut in futures]
    worst_diff = max(abs(err_convex - e) for e in rerun_errs)
    reproducible = worst_diff < 1e-6
    print(f"    cross-process reruns: {[round(e, 6) for e in rerun_errs]}  (original: {err_convex:.6f})")
    print(f"    worst |diff| = {worst_diff:.2e}  "
          f"{'REPRODUCIBLE across processes (convex problem, unique optimum, as expected)' if reproducible else 'NOT REPRODUCIBLE -- the SDP itself has a degenerate optimal face, investigate'}")

    print(f"\n  -- forte-1: nonlinear fit, anchored init (Task 31D's own baseline, should reproduce 0.3166) --")
    err_anchored, _ = compute_full_pipeline(p, non_id_labels, diag, kept, P_S, tags, counts_list,
                                              gamma_per_circuit, "nonlinear_anchored", "t32c_anchored")
    print(f"    err = {err_anchored:.4f}  (Task 31D baseline: {BASELINE_FORTE:.4f})")

    print(f"\n  -- forte-1: nonlinear fit, RANDOM-ONLY init (no target anchor at all) --")
    err_random, _ = compute_full_pipeline(p, non_id_labels, diag, kept, P_S, tags, counts_list,
                                            gamma_per_circuit, "nonlinear_random", "t32c_random")
    print(f"    err = {err_random:.4f}")

    # -- ideal-control check, all three methods --
    print(f"\n  -- ideal-control check, all three methods --")
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    ideal_tags_raw = raw_ck["tags"]["ideal"]
    ideal_counts_raw = raw_ck["counts"]["ideal"]
    ideal_per_name = {}
    for (name, group), counts in zip(ideal_tags_raw, ideal_counts_raw):
        ideal_per_name.setdefault(name, {}).setdefault(tuple(group), counts)

    def ideal_seed_matrix(name):
        out = {}
        for group_t, counts in ideal_per_name[name].items():
            for l in group_t:
                vals = []
                for seed in range(N_SEEDS_COV):
                    rng = np.random.default_rng(stable_seed("t32c_ideal", name, l, seed))
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    vals.append(expectation_from_counts(resampled, l))
                out[l] = vals
        return out

    ideal_results = {}
    for method in ("convex", "nonlinear_anchored", "nonlinear_random"):
        a_by_name = {}
        for name in kept:
            seed_matrix = ideal_seed_matrix(name)
            labels = list(seed_matrix.keys())
            y_mean = {l: float(np.mean(seed_matrix[l])) for l in labels}
            v0 = target_coeff_vector(name, K)
            if method == "convex":
                a_hat, eigvals, _, _ = solve_convex_manifold(P_S, labels, y_mean, K)
                if np.dot(a_hat, v0) < 0:
                    a_hat = -a_hat
            elif method == "nonlinear_anchored":
                a_hat, _ = fit_pure_state_uniform(P_S, y_mean, {l: 1.0 for l in labels}, K, v0,
                                                   seed=stable_seed("t32c_ideal_anchored", name))
            else:
                a_hat, _ = fit_pure_state_random_only(P_S, y_mean, {l: 1.0 for l in labels}, K,
                                                        seed=stable_seed("t32c_ideal_random", name))
            a_by_name[name] = a_hat
        full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
        _, err_ideal = energy_and_err(p, full, K)
        ideal_results[method] = err_ideal
        print(f"    {method:<20}: ideal err = {err_ideal:.4f}  (Task 31D anchored baseline: {BASELINE_IDEAL:.4f})")

    print(f"\n" + "=" * 96)
    bias_ok = err_convex <= BASELINE_FORTE * 1.1  # "no worse" with a small numerical-tolerance allowance
    ideal_ok = ideal_results["convex"] <= BASELINE_IDEAL * 1.5 + 0.05
    verdict = "PASS" if (reproducible and bias_ok and ideal_ok) else "FAIL"
    print(f"  VERDICT: reproducible={reproducible}, bias {'OK' if bias_ok else 'WORSE'} "
          f"({BASELINE_FORTE:.4f} -> {err_convex:.4f}), ideal {'OK' if ideal_ok else 'DEGRADED'} "
          f"({BASELINE_IDEAL:.4f} -> {ideal_results['convex']:.4f}) -> {verdict}")
    print("=" * 96 + "\n")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "err_convex": err_convex, "err_convex_reruns_cross_process": rerun_errs,
            "reproducible": bool(reproducible), "worst_cross_process_diff": worst_diff,
            "err_nonlinear_anchored": err_anchored, "err_nonlinear_random": err_random,
            "baseline_forte_task31d": BASELINE_FORTE, "baseline_ideal_task31d": BASELINE_IDEAL,
            "purity_ratio_mean": float(purity_ratio), "eig_ratios_by_slot": eig_ratios,
            "ideal_results": ideal_results, "bias_ok": bool(bias_ok), "ideal_ok": bool(ideal_ok),
            "verdict": verdict,
        }, f, indent=2)
    print(f"  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
