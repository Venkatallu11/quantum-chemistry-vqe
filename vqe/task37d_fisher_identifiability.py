#!/usr/bin/env python3
"""
task37d_fisher_identifiability.py -- iteration 37, Task D. WHY did Task
37C's PEC-corrected fit return all 5 residual noise parameters
statistically indistinguishable from their priors (|z| <= 0.30 for every
one)? Two very different explanations are consistent with that same
observation: (a) the true residual noise really is ~0 (a genuine physics
null), or (b) the data structurally CANNOT distinguish a noise-parameter
shift from some compensating adjustment of the 15-parameter Schmidt
frame, so the fit correctly falls back to the prior in a direction it
has no real information about, regardless of the true value. This file
answers that with the proposal's own prescribed tool: the Fisher
information matrix, F = J^T J (Gauss-Newton/Cramer-Rao approximation,
standard for nonlinear least squares -- the weighted residuals already
encode 1/sigma, so no separate covariance weighting is needed), evaluated
on the SAME training data and near the SAME fitted point as Task 37C's
real result, DATA ONLY (prior_lambda=0, so this measures what the H4
measurements themselves constrain, not what the prior adds).

THE KEY STEP, and the reason a naive per-parameter Fisher/variance would
be MISLEADING here: with 15 frame parameters free to move alongside the
5 noise parameters, some noise-parameter directions may be (partially)
degenerate with some frame adjustment -- a change in assumed p_gpi2, say,
rescales predictions in a way a small frame rotation could also
approximately produce. The correct way to ask "how well does the DATA
alone constrain the noise parameters, accounting for the frame being
free to compensate" is the SCHUR COMPLEMENT of the noise-noise Fisher
block after marginalizing out the frame block:
    F_eff = F_nn - F_nf @ pinv(F_ff) @ F_fn
(the same construction used throughout statistics for profiling out
nuisance parameters from a Fisher/covariance matrix). Small eigenvalues
of F_eff mean the data genuinely cannot pin down that noise direction
NO MATTER how good the frame fit is -- large eigenvalues mean the data
really does constrain it, and a near-zero fitted value is closer to
being a real physics finding. Both the naive (frame-ignoring) and
Schur-complement (frame-aware) per-parameter Fisher information are
reported side by side, specifically so the SIZE of the gap between them
is visible -- that gap IS the quantitative answer to "how much does the
frame's flexibility eat the noise identifiability."

Run:
    PYTHONHASHSEED=0 python vqe/task37d_fisher_identifiability.py
"""
import os
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task37b_h4_noise_model import PARAM_NAMES
from task37c_pec_corrected_joint_fit import (
    joint_residuals, clamp_noise, PRIORS_RESIDUAL, load_real_blended, split_train_val,
    VAL_FRACTION)
from task36_joint_schmidt_frame import N_THETA

K = 6
N_NOISE_PARAMS = len(PARAM_NAMES)


def find_fit_point(U0, P_S, K, kept, non_id_labels, fixed_solutions, train_blended, weight_unit):
    """A single, well-started LM run (not 8 restarts -- Task 37C already
    established the global optimum region; here we only need A
    representative near-optimal point for LOCAL curvature analysis, not
    to re-litigate which restart is globally best)."""
    theta0 = np.zeros(N_THETA)
    noise0 = np.array([0.0, PRIORS_RESIDUAL["p_gpi2"]["mean"], 0.0, 0.0, PRIORS_RESIDUAL["p_readout"]["mean"]])
    params0 = np.concatenate([theta0, noise0])
    sol = least_squares(joint_residuals, params0,
                         args=(U0, P_S, K, kept, non_id_labels, fixed_solutions, train_blended, weight_unit, 1.0),
                         method="lm", max_nfev=4000)
    return sol.x


def main():
    print("\n" + "=" * 96)
    print("  task37d_fisher_identifiability.py -- Fisher information, frame-marginalized, for Task 37C's noise params")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    real_blended, n_mc_total = load_real_blended(kept, non_id_labels)
    train_blended, _ = split_train_val(real_blended, kept, VAL_FRACTION, seed=37)  # SAME split as Task 37C's report
    n_train = sum(len(train_blended[n]) for n in kept)
    print(f"\n  loaded {n_mc_total} MC draws/label, using the SAME training split as Task 37C ({n_train} residuals)")

    U0 = np.eye(K)
    print("\n  -- finding a representative near-optimal fit point (single LM run) --")
    params_hat = find_fit_point(U0, P_S, K, kept, non_id_labels, fixed_solutions, train_blended, weight_unit)
    noise_hat = clamp_noise(params_hat)
    print(f"    noise_hat at this point: {dict(zip(PARAM_NAMES, noise_hat.round(6)))}")

    print("\n  -- computing DATA-ONLY Jacobian (prior_lambda=0) via scipy's verified finite-difference machinery --")
    def resid_data_only(params):
        return joint_residuals(params, U0, P_S, K, kept, non_id_labels, fixed_solutions, train_blended,
                                weight_unit, prior_lambda=0.0)
    r0 = resid_data_only(params_hat)
    # BUG FOUND AND FIXED HERE, disclosed: an earlier version passed rel_step=1e-5. Per scipy's own
    # docs, when rel_step is EXPLICITLY given (not left None), method='3-point' uses
    # h = rel_step * sign(x0) * abs(x0) -- NO max(1,|x0|) floor. For any parameter sitting at or
    # near 0 (p_zz, delta_zz, delta_gpi2, p_readout all fit near 0 here), this collapses the step
    # to ~0 or exactly 0 (sign(0)=0), silently zeroing that Jacobian column. Caught by a manual,
    # generously-sized (h=1e-3) perturbation check that showed REAL, large, nonzero sensitivity for
    # p_gpi2, p_readout, and delta_gpi2 where the auto Jacobian reported ~0 -- see the sanity-check
    # block below, kept in the script rather than deleted. Fixed by using an explicit, uniform
    # ABSOLUTE step instead (no dependence on x0's own scale at all).
    ABS_STEP = 1e-3
    J = approx_derivative(resid_data_only, params_hat, method="3-point", abs_step=ABS_STEP)
    print(f"    Jacobian shape: {J.shape}  (residuals x params, {N_THETA} frame + {N_NOISE_PARAMS} noise, abs_step={ABS_STEP})")

    print("\n  -- SANITY CHECK: manual, generously-sized (0.001) perturbation vs approx_derivative's step (now abs_step-based) --")
    print("     (catches finite-difference step-size artifacts before trusting the Jacobian below)")
    manual_step = 1e-3
    for i, pname in enumerate(PARAM_NAMES):
        col = N_THETA + i
        p_plus = params_hat.copy(); p_plus[col] += manual_step
        p_minus = params_hat.copy(); p_minus[col] -= manual_step
        r_plus = resid_data_only(p_plus)
        r_minus = resid_data_only(p_minus)
        manual_deriv_norm = float(np.linalg.norm((r_plus - r_minus) / (2 * manual_step)))
        auto_deriv_norm = float(np.linalg.norm(J[:, col]))
        flag = "" if abs(manual_deriv_norm - auto_deriv_norm) < 0.1 * max(manual_deriv_norm, 1e-12) else "  <-- MISMATCH, approx_derivative's auto step is likely wrong for this param"
        print(f"    {pname:<12} manual(h=1e-3) ||dr/dtheta||={manual_deriv_norm:.4e}   "
              f"auto(approx_derivative) ||dr/dtheta||={auto_deriv_norm:.4e}{flag}")

    print("\n  -- p_readout ONE-SIDED check (its fitted value may sit at its clamp lower bound=0.0, "
          "which would crush a naive central-difference estimate) --")
    ro_idx = N_THETA + PARAM_NAMES.index("p_readout")
    print(f"    p_readout_hat = {params_hat[ro_idx]:.6f}  (clamp lower bound = {PRIORS_RESIDUAL['p_readout']['bounds'][0]})")
    p_fwd1 = params_hat.copy(); p_fwd1[ro_idx] += 1e-3
    p_fwd2 = params_hat.copy(); p_fwd2[ro_idx] += 2e-3
    r_fwd0, r_fwd1, r_fwd2 = resid_data_only(params_hat), resid_data_only(p_fwd1), resid_data_only(p_fwd2)
    one_sided_deriv = float(np.linalg.norm((r_fwd1 - r_fwd0) / 1e-3))
    print(f"    one-sided forward derivative norm (h=1e-3): {one_sided_deriv:.4e}")
    print(f"    one-sided forward derivative norm (h=2e-3): {float(np.linalg.norm((r_fwd2 - r_fwd0) / 2e-3)):.4e}"
          f"  (should roughly agree with the h=1e-3 estimate if genuinely linear near this point)")

    F = J.T @ J  # Fisher information (Gauss-Newton approx), DATA ONLY
    F_ff = F[:N_THETA, :N_THETA]
    F_fn = F[:N_THETA, N_THETA:]
    F_nf = F[N_THETA:, :N_THETA]
    F_nn = F[N_THETA:, N_THETA:]

    # -- naive (frame-ignoring) per-parameter Fisher info: diagonal of F_nn alone --
    naive_info = np.diag(F_nn)

    # -- Schur complement: noise-parameter Fisher AFTER marginalizing out (profiling over) the
    # frame -- pinv used defensively (rcond default) in case F_ff has any genuinely flat
    # directions of its own, which would make a plain inverse blow up meaninglessly --
    F_ff_pinv = np.linalg.pinv(F_ff, rcond=1e-10)
    F_eff = F_nn - F_nf @ F_ff_pinv @ F_fn
    F_eff = (F_eff + F_eff.T) / 2  # symmetrize away finite-difference asymmetry noise

    eigvals, eigvecs = np.linalg.eigh(F_eff)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    print("\n  -- NAIVE (frame-ignoring) vs SCHUR-COMPLEMENT (frame-marginalized) per-parameter Fisher info --")
    print(f"  {'param':<12}{'naive F_ii':>14}{'naive sigma':>14}{'F_eff_ii':>14}{'eff sigma':>14}{'prior sigma':>14}{'data vs prior':>16}")
    for i, name in enumerate(PARAM_NAMES):
        naive_sigma = 1.0 / np.sqrt(naive_info[i]) if naive_info[i] > 1e-12 else float("inf")
        eff_sigma = 1.0 / np.sqrt(max(F_eff[i, i], 1e-12)) if F_eff[i, i] > 1e-12 else float("inf")
        prior_sigma = PRIORS_RESIDUAL[name]["std"]
        verdict = "DATA TIGHTER" if eff_sigma < prior_sigma else "prior still tighter"
        print(f"  {name:<12}{naive_info[i]:>14.3e}{naive_sigma:>14.3e}{F_eff[i,i]:>14.3e}{eff_sigma:>14.3e}"
              f"{prior_sigma:>14.3e}{verdict:>16}")

    print("\n  -- SCHUR-COMPLEMENT eigenvalues (marginal noise-subspace identifiability, most to least constrained) --")
    for k in range(N_NOISE_PARAMS):
        vec_str = ", ".join(f"{PARAM_NAMES[j]}={eigvecs[j,k]:+.3f}" for j in range(N_NOISE_PARAMS))
        implied_sigma = 1.0 / np.sqrt(eigvals[k]) if eigvals[k] > 1e-12 else float("inf")
        print(f"    eigval[{k}]={eigvals[k]:.3e}  (implied sigma along this direction: {implied_sigma:.3e})")
        print(f"      direction: {vec_str}")

    cond = eigvals[0] / max(eigvals[-1], 1e-300)
    print(f"\n  condition number (max/min eigenvalue) of the frame-marginalized noise Fisher block: {cond:.3e}")
    print(f"\n  -- HONEST READ --")
    n_data_tighter = sum(1 for i in range(N_NOISE_PARAMS)
                          if (1.0/np.sqrt(max(F_eff[i,i],1e-12))) < PRIORS_RESIDUAL[PARAM_NAMES[i]]["std"])
    if n_data_tighter == 0:
        print(f"    For EVERY one of the 5 noise parameters, the frame-marginalized data uncertainty is WIDER than "
              f"the prior's own uncertainty -- the fitted values landing near the prior mean (Task 37C's own "
              f"z~0 result) is consistent with the data providing close to NO additional constraint beyond the "
              f"prior in any noise direction, once the frame is allowed to compensate. This does NOT prove the "
              f"true residual noise is exactly zero -- it means this data, at this shot/MC-draw budget, with this "
              f"circuit set, cannot tell the difference between zero and anything within the prior's own width.")
    else:
        names_tighter = [PARAM_NAMES[i] for i in range(N_NOISE_PARAMS)
                          if (1.0/np.sqrt(max(F_eff[i,i],1e-12))) < PRIORS_RESIDUAL[PARAM_NAMES[i]]["std"]]
        print(f"    {n_data_tighter}/5 noise parameters ({', '.join(names_tighter)}) ARE more tightly constrained "
              f"by the data than by the prior, even after frame-marginalization -- Task 37C's near-zero fitted "
              f"values for these specific parameters are a more genuine physics finding, not just prior fallback.")


if __name__ == "__main__":
    main()
