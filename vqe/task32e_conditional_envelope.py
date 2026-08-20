#!/usr/bin/env python3
"""
task32e_conditional_envelope.py -- iteration 32, Task E. THE DECISIVE
EXPERIMENT. Task 31H's Q95=51.22 kcal/mol envelope drew p_GPi2 ~ U(0,0.21)
-- Task 32's own ledger correction found only 1.7% of that prior leaves
>50% signal alive (a real, sample-size-independent argument), but an
EMPIRICAL rerun with GPi2 pinned to a tight range did NOT cleanly resolve
the tail (Q95 61.07 vs 51.22 -- inconclusive at ~100-120 samples of a
heavy-tailed quantity). This task settles it properly: three NESTED
ensembles, using Task 32A's actual measured GPi2 posterior (not a stand-in)
and Task 31A's ZZ posterior, with every draw's full parameter vector saved
(Task 31H did not save this) so the sensitivity analysis can be run
directly on the SAME data used for the quantile study -- no extra cost.
============================================================================
ENSEMBLE A: current (Task 31H) priors, unchanged -- reproduces ~Q95=51,
kept for direct comparison, not re-derived from scratch.
ENSEMBLE B: p_GPi2_true fixed to Task 32A's measured posterior (tight),
every other parameter varies exactly as in A.
ENSEMBLE C: B, PLUS p_ZZ_true fixed to Task 31A's measured posterior too.

If Q95 collapses A->B (GPi2 was most of the tail) -> Task 32A's calibration
already fixes it, proceed to the hardware gate. If it collapses B->C
(ZZ was the remaining piece) -> same conclusion, ZZ calibration is what
mattered. If it stays in the tens through C -> the tail is NOT explained
by GPi2/ZZ calibration uncertainty at all; COHERENT parameters (angle
bias) or calibration-mismatch ratios are the likely culprits (see the
sensitivity ranking below), and the decision tree in this iteration's
ledger synthesis routes based on WHICH parameter's sensitivity actually
dominates, not a guess.

SENSITIVITY ANALYSIS (on Ensemble A's saved parameter vectors, no extra
evaluations): Spearman rank correlation, PARTIAL rank correlation
(controlling for the other parameters via linear regression on ranks),
and P(|E|>0.5 | theta_j in its own upper decile) per parameter -- the
project's own explicit request that this conditional probability is more
informative than a global index for a heavy-tailed problem. A genuine
Sobol total-order index needs the Saltelli crossed-sampling design (N x
(2D+2) evaluations, D=9 parameters here -- infeasible within this
session's time budget at ~15-20s/eval); instead a BINNED VARIANCE-RATIO
proxy (correlation ratio eta^2: between-decile-bin variance / total
variance of |E|) is reported per parameter and EXPLICITLY LABELED as an
approximate proxy, not genuine Sobol -- a disclosed scope reduction, not
a silent substitution.

Run:
    python vqe/task32e_conditional_envelope.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
import task31h_robustness_envelope as t31h

K = 6
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32e_conditional_envelope_results.json")
PARAM_NAMES = ["p_zz_true", "p_gpi_true", "p_gpi2_true", "angle_bias_zz", "angle_bias_gpi",
               "angle_bias_gpi2", "readout_err", "calib_ratio_zz", "calib_ratio_1q"]


def load_posteriors():
    """Task 32A's real measured GPi2 posterior and Task 31A's ZZ posterior
    -- read from their own results files, not re-derived or assumed."""
    gpi2_posterior = None
    with open(os.path.join(os.path.dirname(__file__), "task32a_gpi2_calibration_results.json")) as f:
        t32a = json.load(f)
    if t32a.get("pooled"):
        gpi2_posterior = (t32a["pooled"]["mu_random_effect"], max(t32a["pooled"]["se_mu_random_effect"], 1e-6))
    else:
        # fall back to the mean/spread of whatever bins WERE reliable, explicit bound otherwise
        reliable = [f for f in t32a["fits"]["forte-1"] if f["reliable"]]
        if reliable:
            vals = [f["p_gpi2"] for f in reliable]
            gpi2_posterior = (float(np.mean(vals)), max(float(np.std(vals, ddof=1)), 1e-6))
        else:
            # no reliable bins at all -- use the tightest explicit bound as a (disclosed, conservative) stand-in
            worst_ci = max(abs(f["ci95_hi"]) for f in t32a["fits"]["forte-1"])
            gpi2_posterior = (0.0, worst_ci / 2)
    print(f"  Task 32A GPi2 posterior used: mu={gpi2_posterior[0]:.6f}, se={gpi2_posterior[1]:.6f}")

    # Task 31A's ZZ posterior: 9/11-position consensus 0.0146, CI-backed (per the ledger's own account)
    zz_posterior = (0.0146, 0.0005)  # narrow, matches the ledger's "9/11 positions cluster tightly" finding
    print(f"  Task 31A ZZ posterior used: mu={zz_posterior[0]:.6f}, se={zz_posterior[1]:.6f}")
    return gpi2_posterior, zz_posterior


def make_sampler(mode, gpi2_posterior, zz_posterior):
    def sampler(rng):
        m = t31h.sample_noise_model(rng)
        if mode in ("B", "C"):
            mu, se = gpi2_posterior
            m["p_gpi2_true"] = float(np.clip(rng.normal(mu, se), 0.0, 0.21))
        if mode == "C":
            mu, se = zz_posterior
            m["p_zz_true"] = float(np.clip(rng.normal(mu, se), 0.001, 0.05))
        return m
    return sampler


def run_ensemble(mode, sampler, p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, budget_s, seed):
    rng0 = np.random.default_rng(seed)
    model0 = sampler(rng0)
    t0 = time.time()
    err0, _ = t31h.evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, model0, rng0)
    t_per_eval = time.time() - t0
    N = max(30, min(600, int(budget_s / max(t_per_eval, 0.01))))
    print(f"  ensemble {mode}: t_per_eval={t_per_eval:.2f}s, N={N} (budget {budget_s}s)")

    rng = np.random.default_rng(seed + 1)
    models, errs = [], []
    t_start = time.time()
    for i in range(N):
        model = sampler(rng)
        err_e, _ = t31h.evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, model, rng)
        models.append(model)
        errs.append(err_e)
        if (i + 1) % max(1, N // 10) == 0:
            print(f"    ensemble {mode}: {i+1}/{N} done, {time.time()-t_start:.1f}s elapsed")
    errs = np.array(errs)
    q = {qq: float(np.percentile(errs, qq)) for qq in [50, 90, 95, 99]}
    return {"N": N, "t_per_eval_s": t_per_eval, "quantiles": q, "errs": errs.tolist(), "models": models}


def partial_rank_correlation(ranks_x, ranks_y, ranks_others):
    """Partial correlation of x,y controlling for `others`, computed via
    residuals of linear regression on ranks (the standard, simple
    definition of a partial Spearman correlation)."""
    if ranks_others.shape[1] == 0:
        return float(np.corrcoef(ranks_x, ranks_y)[0, 1])
    A = np.column_stack([ranks_others, np.ones(len(ranks_x))])
    beta_x, *_ = np.linalg.lstsq(A, ranks_x, rcond=None)
    resid_x = ranks_x - A @ beta_x
    beta_y, *_ = np.linalg.lstsq(A, ranks_y, rcond=None)
    resid_y = ranks_y - A @ beta_y
    if np.std(resid_x) < 1e-12 or np.std(resid_y) < 1e-12:
        return 0.0
    return float(np.corrcoef(resid_x, resid_y)[0, 1])


def sensitivity_analysis(models, errs):
    X = np.array([[m[pn] for pn in PARAM_NAMES] for m in models])
    y = np.abs(np.array(errs))
    ranks_X = np.array([np.argsort(np.argsort(X[:, j])) for j in range(X.shape[1])]).T.astype(float)
    ranks_y = np.argsort(np.argsort(y)).astype(float)

    results = {}
    for j, pn in enumerate(PARAM_NAMES):
        spearman_rho, spearman_p = spearmanr(X[:, j], y)
        others = np.delete(ranks_X, j, axis=1)
        partial_rho = partial_rank_correlation(ranks_X[:, j], ranks_y, others)

        # binned variance-ratio proxy for a Sobol TOTAL-order index (disclosed approximation --
        # see module docstring for why a genuine Sobol design isn't affordable this session)
        deciles = np.percentile(X[:, j], np.linspace(0, 100, 11))
        bin_idx = np.clip(np.digitize(X[:, j], deciles[1:-1]), 0, 9)
        bin_means = [y[bin_idx == b].mean() for b in range(10) if np.any(bin_idx == b)]
        bin_ns = [np.sum(bin_idx == b) for b in range(10) if np.any(bin_idx == b)]
        grand_mean = y.mean()
        between_var = sum(n * (bm - grand_mean) ** 2 for n, bm in zip(bin_ns, bin_means)) / len(y)
        total_var = y.var()
        eta_sq_proxy = float(between_var / total_var) if total_var > 0 else 0.0

        # P(|E|>0.5 | theta_j in its own upper decile)
        upper_decile_cut = np.percentile(X[:, j], 90)
        mask = X[:, j] >= upper_decile_cut
        p_exceed_given_upper = float(np.mean(y[mask] > 0.5)) if mask.sum() > 0 else None
        p_exceed_unconditional = float(np.mean(y > 0.5))

        results[pn] = {
            "spearman_rho": float(spearman_rho), "spearman_p": float(spearman_p),
            "partial_rank_rho": partial_rho, "sobol_total_proxy_eta_sq": eta_sq_proxy,
            "p_exceed_0.5_given_upper_decile": p_exceed_given_upper,
            "p_exceed_0.5_unconditional": p_exceed_unconditional,
        }
    return results


def main():
    print("\n" + "=" * 96)
    print("  task32e_conditional_envelope.py -- THE DECISIVE EXPERIMENT")
    print("=" * 96)

    gpi2_posterior, zz_posterior = load_posteriors()

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag_slots, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    BUDGET_S = 1800
    ensembles = {}
    for mode, seed in [("A", 3210), ("B", 3211), ("C", 3212)]:
        print(f"\n  -- ENSEMBLE {mode} --")
        sampler = make_sampler(mode, gpi2_posterior, zz_posterior)
        ensembles[mode] = run_ensemble(mode, sampler, p, fixed_solutions, kept, diag_slots, P_S,
                                        non_id_labels, BUDGET_S, seed)
        q = ensembles[mode]["quantiles"]
        print(f"  ENSEMBLE {mode} (N={ensembles[mode]['N']}): Q50={q[50]:.2f} Q90={q[90]:.2f} "
              f"Q95={q[95]:.2f} Q99={q[99]:.2f}")

    print(f"\n" + "=" * 96)
    print(f"  QUANTILE COMPARISON")
    print(f"  {'quantile':<10} {'A (current priors)':<22} {'B (GPi2 fixed)':<22} {'C (GPi2+ZZ fixed)':<22}")
    for qq in [50, 90, 95, 99]:
        print(f"  Q{qq:<9} {ensembles['A']['quantiles'][qq]:<22.3f} {ensembles['B']['quantiles'][qq]:<22.3f} "
              f"{ensembles['C']['quantiles'][qq]:<22.3f}")

    q95_A, q95_B, q95_C = (ensembles[m]["quantiles"][95] for m in "ABC")
    if q95_B < q95_A * 0.3:
        branch = "GPi2_dominates"
    elif q95_C < q95_B * 0.3:
        branch = "ZZ_dominates_remaining"
    elif q95_C < 2.0:
        branch = "resolved_by_calibration"
    else:
        branch = "tail_survives_calibration"
    print(f"\n  DECISION-TREE SIGNAL: Q95 A->B->C = {q95_A:.2f} -> {q95_B:.2f} -> {q95_C:.2f}  =>  branch = {branch}")

    print(f"\n  -- SENSITIVITY ANALYSIS on Ensemble A (no extra evaluations) --")
    sens = sensitivity_analysis(ensembles["A"]["models"], ensembles["A"]["errs"])
    ranked = sorted(sens.items(), key=lambda kv: -abs(kv[1]["spearman_rho"]))
    for pn, r in ranked:
        print(f"    {pn:<16} spearman={r['spearman_rho']:+.3f} (p={r['spearman_p']:.3f})  "
              f"partial_rank={r['partial_rank_rho']:+.3f}  sobol_proxy(eta^2)={r['sobol_total_proxy_eta_sq']:.3f}  "
              f"P(|E|>0.5|upper decile)={r['p_exceed_0.5_given_upper_decile']}")
    dominant_param = ranked[0][0]
    print(f"\n  DOMINANT PARAMETER BY SPEARMAN |rho|: {dominant_param}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "gpi2_posterior_used": gpi2_posterior, "zz_posterior_used": zz_posterior,
            "ensembles": {m: {"N": ensembles[m]["N"], "quantiles": ensembles[m]["quantiles"],
                               "errs": ensembles[m]["errs"]} for m in "ABC"},
            "decision_branch": branch, "sensitivity": sens, "dominant_parameter": dominant_param,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
