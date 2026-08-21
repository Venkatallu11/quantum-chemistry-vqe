#!/usr/bin/env python3
"""
task37_phase1a_joint_noise_and_state.py -- iteration 37, Phase 1A. Joint
estimation of the H4 Schmidt frame AND the assumed PEC calibration
parameters, with the identifiability safeguard the proposal insisted on:
a REAL calibration prior anchors p_zz_assumed near Task 31A's actual
measurement (mean=0.014593, std=0.000124); p_gpi2_assumed gets NO
artificial prior, because none exists (Task 31A's own recalibration
attempt is flagged "physical": false) -- this is deliberately the ONE
parameter allowed to be determined by the H4 data itself, since
calibration alone has twice failed to pin it down.

SCOPE, disclosed: the proposal's own suggested first set was 5 parameters
(p_ZZ, p_GPi2, delta_ZZ, delta_GPi2, p_readout). This prototype uses 2
(p_zz_assumed, p_gpi2_assumed) -- the two parameters `analytic_A_and_B`
(Task 30B/32I's own established, already-validated correction machinery)
already supports directly, with no new noise-channel code to write and
independently verify. Coherent angle-bias and readout parameters are
real, disclosed future work, not silently dropped.

MODEL: for candidate (U, p_zz_assumed, p_gpi2_assumed), the predicted
corrected value for slot `name`, label `l` is
    m_corrected = m_raw_real * B(p_assumed)/A(p_assumed)
using `analytic_A_and_B` exactly as Task 30B/31F/32I already established
(same function, unmodified) -- NOT a new noise-forward-model. This is
compared against the joint frame's bilinear prediction v(U)@P_S[l]@v(U).
Objective = weighted sum of squared residuals + a Gaussian prior penalty
pulling p_zz_assumed toward its real calibration mean/std. p_gpi2_assumed
is bounded to [0, 0.21] (Task 31A's own established bound) but otherwise
free.

DATA: `task28b_optimized_raw.json`'s real forte-1 RAW (pre-PEC)
measurements -- the actual real hardware data this correction is meant
to act on, not the already-PEC-corrected Task 31C checkpoint Task 36
used (that data has calibration baked in already and can't be used to
re-derive it).

VALIDATION, mandatory before trusting anything (Task 36's own
discipline): ideal-data recovery, adversarial rejection on shuffled real
data.

Run:
    PYTHONHASHSEED=0 python vqe/task37_phase1a_joint_noise_and_state.py
"""
import os
import sys
import json
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state
from phys_constrained_reconstruction import build_P_S
from task30b_pec_application import analytic_A_and_B
from task36_joint_schmidt_frame import (skew_from_theta, U_from_theta, slot_vector,
                                          build_full_from_frame, procrustes_init, N_THETA)
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts, stable_seed
from qiskit.quantum_info import Statevector, Pauli

K = 6
GATE_NAME = "zz"
SHOTS = 100_000
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
ZZ_REAL_MEAN, ZZ_REAL_STD = 0.014593, 0.000124
GPI_REAL_P1 = 0.000119  # Task 30B's real GPi calibration, used as the fixed gpi_bins fallback (unchanged)
N_NOISE_PARAMS = 2  # p_zz_assumed, p_gpi2_assumed
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37_phase1a_joint_noise_and_state_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def build_full_from_frame_local(U, P_S, K, non_id_labels, kept):
    return build_full_from_frame(U, P_S, K, non_id_labels, kept)


_GLOBAL_AB_CACHE = {}  # keyed on (name, rounded p_zz, rounded p_gpi2) -- A/B do NOT depend on theta_state
# (the 15 Schmidt-frame parameters) at all, only on the 2 noise parameters and the FIXED circuit
# angles. Recomputing this per residual-function call (as an earlier version of this script did)
# meant redoing expensive 21-slot density-matrix propagation on every LM step even when only
# theta_state moved -- a real, severe performance bug, caught after Phase 1A hung with zero
# progress for 20+ minutes. Rounding to 6 decimals means near-identical noise-parameter proposals
# (common during a line search) reuse the same cached A/B instead of recomputing.


def joint_residuals(params, U0, P_S, K, kept, non_id_labels, blended_by_slot, fixed_solutions,
                     gpi_bins_fixed, weight_by_slot, prior_lambda):
    theta_state = params[:N_THETA]
    p_zz_assumed = abs(params[N_THETA])
    p_gpi2_assumed = min(0.21, max(0.0, abs(params[N_THETA + 1])))
    U = U_from_theta(theta_state, U0, K)

    res = []
    for name in kept:
        cache_key = (name, round(p_zz_assumed, 6), round(p_gpi2_assumed, 6))
        if cache_key not in _GLOBAL_AB_CACHE:
            labels_here = list(blended_by_slot[name].keys())
            A, B = analytic_A_and_B(fixed_solutions[name]["angles"], GATE_NAME, p_zz_assumed,
                                     gpi_bins_fixed, gpi_bins_fixed, labels_here)
            _GLOBAL_AB_CACHE[cache_key] = (A, B)
        A, B = _GLOBAL_AB_CACHE[cache_key]
        v = slot_vector(U, name, K)
        for l, m_raw in blended_by_slot[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            m_corrected = max(-1.0, min(1.0, m_raw * ratio))
            pred = float(np.real(v @ P_S[l] @ v))
            w = weight_by_slot[name].get(l, 1.0)
            res.append(np.sqrt(w) * (m_corrected - pred))

    # -- calibration prior: pulls p_zz_assumed toward its REAL measured value. p_gpi2_assumed
    # gets NO prior term here -- deliberately, since no real posterior exists for it. --
    res.append(np.sqrt(prior_lambda) * (p_zz_assumed - ZZ_REAL_MEAN) / ZZ_REAL_STD)
    return np.array(res)


def fit_joint(U0, P_S, K, kept, non_id_labels, blended_by_slot, fixed_solutions, gpi_bins_fixed,
              weight_by_slot, prior_lambda, rng, n_restarts=8, noise0_options=None):
    """noise0_options: list of (p_zz0, p_gpi2_0) starting points to try, one per restart (cycled
    if fewer than n_restarts). BUG FIXED here: an earlier version hard-coded the SAME noise
    starting point for every restart, so if that one start's local search couldn't reach the true
    answer, no restart ever could either -- caught when ideal-data recovery (which needs
    p_zz_assumed to reach ~0, not the calibration value ~0.0146) failed at err=0.050 with
    p_zz_hat frozen exactly at its unperturbed starting point."""
    best = None
    scales = [0.0, 0.02, 0.05, 0.08, 0.12, 0.15, 0.2, 0.3][:n_restarts]
    if noise0_options is None:
        noise0_options = [(ZZ_REAL_MEAN, GPI_REAL_P1)]
    for i, s in enumerate(scales):
        theta0 = np.zeros(N_THETA) if s == 0.0 else rng.normal(0, s, N_THETA)
        noise0 = np.array(noise0_options[i % len(noise0_options)])
        params0 = np.concatenate([theta0, noise0])
        try:
            sol = least_squares(joint_residuals, params0,
                                 args=(U0, P_S, K, kept, non_id_labels, blended_by_slot, fixed_solutions,
                                       gpi_bins_fixed, weight_by_slot, prior_lambda),
                                 method="lm", max_nfev=3000)
        except Exception:
            continue
        cost = float(2 * sol.cost)
        if best is None or cost < best[0]:
            best = (cost, sol.x)
    cost, params_hat = best
    theta_hat = params_hat[:N_THETA]
    p_zz_hat = abs(params_hat[N_THETA])
    p_gpi2_hat = min(0.21, max(0.0, abs(params_hat[N_THETA + 1])))
    U_hat = U_from_theta(theta_hat, U0, K)
    n_resid = sum(len(blended_by_slot[name]) for name in kept) + 1
    dof = max(1, n_resid - (N_THETA + N_NOISE_PARAMS))
    return U_hat, p_zz_hat, p_gpi2_hat, cost, cost / dof


def load_real_blended(kept, non_id_labels, n_seeds=8, p2_for_ideal_check=None):
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    tags = raw_ck["tags"]["forte-1"]
    counts_list = raw_ck["counts"]["forte-1"]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)
    blended = {name: {} for name in kept}
    for name in kept:
        labels_here = [l for group_t in per_name[name] for l in group_t]
        for group_t, counts in per_name[name].items():
            for seed in range(n_seeds):
                rng = np.random.default_rng(stable_seed("t37p1a", name, seed))
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    m = expectation_from_counts(resampled, l)
                    blended[name].setdefault(l, []).append(m)
    for name in kept:
        for l in blended[name]:
            blended[name][l] = float(np.mean(blended[name][l]))
    return blended


def main():
    print("\n" + "=" * 96)
    print("  task37_phase1a_joint_noise_and_state.py -- joint Schmidt frame + PEC calibration, data-role separated")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins_fixed = learned["forte-1"]["gpi"]  # GPi kept FIXED -- real, tight calibration (Task 37 Phase 0)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    # ================= VALIDATION 1: IDEAL-DATA RECOVERY =================
    print("\n  -- VALIDATION 1: ideal-data (noiseless raw, p=0) recovery --")
    ideal_blended = {}
    for name in kept:
        v_ideal = slot_vector(np.eye(K), name, K)
        ideal_blended[name] = {l: float(np.real(v_ideal @ P_S[l] @ v_ideal)) for l in non_id_labels}
    rng_val = np.random.default_rng(11)
    # at p_assumed=0, analytic_A_and_B gives A=B=exact_ideal (both the raw and PEC-inverse steps
    # are no-ops), so ratio=1 and corrected==raw==ideal target EXACTLY -- but the optimizer must
    # actually find p_assumed~0 to see this; diverse noise0_options (including near-zero) ensures
    # at least some restarts start where the true answer actually is, not just at the calibration
    # value (which is the WRONG answer for noiseless data, by construction -- PEC's correction is
    # only a no-op when the assumed p matches the true p, and true p=0 here)
    noise0_ideal = [(1e-5, 1e-5), (ZZ_REAL_MEAN, GPI_REAL_P1), (0.005, 0.05), (0.02, 0.1)]
    U_hat_i, pzz_i, pgpi2_i, cost_i, chi2dof_i = fit_joint(
        np.eye(K), P_S, K, kept, non_id_labels, ideal_blended, fixed_solutions, gpi_bins_fixed,
        weight_unit, prior_lambda=0.0, rng=rng_val, noise0_options=noise0_ideal)
    full_i = build_full_from_frame_local(U_hat_i, P_S, K, non_id_labels, kept)
    E_i, err_i = energy_and_err(p, full_i, K)
    print(f"    recovered energy err vs exact: {err_i:.6f} kcal/mol  p_zz_hat={pzz_i:.6f}  p_gpi2_hat={pgpi2_i:.6f}")
    ideal_ok = err_i < 0.05
    print(f"    IDEAL-DATA RECOVERY: {'PASS' if ideal_ok else 'FAIL -- STOP'}")
    if not ideal_ok:
        raise RuntimeError("Phase 1A fails ideal-data recovery -- a real bug, stop before trusting real data")

    # ================= VALIDATION 2: ADVERSARIAL (SHUFFLED REAL DATA) =================
    print("\n  -- VALIDATION 2: adversarial rejection (real raw data, labels shuffled within each slot) --")
    real_blended = load_real_blended(kept, non_id_labels)
    rng_adv = np.random.default_rng(2026)
    shuffled_blended = {}
    for name in kept:
        labels_here = list(real_blended[name].keys())
        vals_here = list(real_blended[name].values())
        perm = rng_adv.permutation(len(vals_here))
        shuffled_blended[name] = {l: vals_here[perm[i]] for i, l in enumerate(labels_here)}
    U0_std = np.eye(K)
    noise0_diverse = [(ZZ_REAL_MEAN, GPI_REAL_P1), (1e-5, 1e-5), (0.02, 0.1), (0.01, 0.05),
                       (ZZ_REAL_MEAN, 0.15), (0.018, 0.02), (0.012, 0.08), (ZZ_REAL_MEAN, 1e-5)]
    _, _, _, cost_real, chi2dof_real = fit_joint(U0_std, P_S, K, kept, non_id_labels, real_blended,
                                                   fixed_solutions, gpi_bins_fixed, weight_unit,
                                                   prior_lambda=1.0, rng=rng_adv, noise0_options=noise0_diverse)
    _, _, _, cost_adv, chi2dof_adv = fit_joint(U0_std, P_S, K, kept, non_id_labels, shuffled_blended,
                                                 fixed_solutions, gpi_bins_fixed, weight_unit,
                                                 prior_lambda=1.0, rng=rng_adv, noise0_options=noise0_diverse)
    print(f"    chi2/dof on REAL data (unshuffled):  {chi2dof_real:.4f}")
    print(f"    chi2/dof on SHUFFLED (adversarial):  {chi2dof_adv:.4f}")
    adv_ok = chi2dof_adv > 3 * chi2dof_real
    print(f"    ADVERSARIAL REJECTION: {'PASS' if adv_ok else 'FAIL -- STOP, not trustworthy'}")
    if not adv_ok:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"ideal_ok": bool(ideal_ok), "adv_ok": bool(adv_ok),
                       "chi2dof_real": chi2dof_real, "chi2dof_adv": chi2dof_adv,
                       "stopped_early": True}, f, indent=2)
        return

    # ================= MAIN FIT: real data, calibration prior on p_zz only =================
    print("\n  -- MAIN FIT: real H4 data + calibration prior (p_zz anchored, p_gpi2 free) --")
    U_hat, p_zz_hat, p_gpi2_hat, cost, chi2dof = fit_joint(
        U0_std, P_S, K, kept, non_id_labels, real_blended, fixed_solutions, gpi_bins_fixed,
        weight_unit, prior_lambda=1.0, rng=np.random.default_rng(37), n_restarts=8,
        noise0_options=noise0_diverse)
    full_joint = build_full_from_frame_local(U_hat, P_S, K, non_id_labels, kept)
    E_joint, err_joint = energy_and_err(p, full_joint, K)

    print(f"    p_zz_hat   = {p_zz_hat:.6f}  (calibration prior: {ZZ_REAL_MEAN:.6f} +/- {ZZ_REAL_STD:.6f}, "
          f"z-score={(p_zz_hat-ZZ_REAL_MEAN)/ZZ_REAL_STD:+.2f})")
    print(f"    p_gpi2_hat = {p_gpi2_hat:.6f}  (no calibration prior -- purely H4-data-determined; "
          f"established bound was |p|<0.21, prior fallback estimate was {GPI_REAL_P1:.6f})")
    print(f"    chi2/dof = {chi2dof:.4f}")
    print(f"    resulting H4 energy error vs exact: {err_joint:.4f} kcal/mol")

    # -- compare against: joint frame with noise FIXED at calibration values (Task 36's original) --
    print(f"\n  -- COMPARISON: does letting p_zz/p_gpi2 float (vs fixed at calibration) change the fit? --")
    U_hat_fixed, _, _, cost_fixed, chi2dof_fixed = fit_joint(
        U0_std, P_S, K, kept, non_id_labels, real_blended, fixed_solutions, gpi_bins_fixed,
        weight_unit, prior_lambda=1e8, rng=np.random.default_rng(38), n_restarts=8,
        noise0_options=noise0_diverse)  # huge prior_lambda ~ fixes p_zz; p_gpi2 still free but see note
    full_fixed = build_full_from_frame_local(U_hat_fixed, P_S, K, non_id_labels, kept)
    _, err_fixed = energy_and_err(p, full_fixed, K)
    print(f"    with p_zz strongly anchored (prior_lambda=1e8): err={err_fixed:.4f} kcal/mol, chi2/dof={chi2dof_fixed:.4f}")
    print(f"    with p_zz weakly anchored (prior_lambda=1):     err={err_joint:.4f} kcal/mol, chi2/dof={chi2dof:.4f}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "ideal_ok": ideal_ok, "adv_ok": adv_ok,
            "chi2dof_real": chi2dof_real, "chi2dof_adv": chi2dof_adv,
            "p_zz_hat": p_zz_hat, "p_gpi2_hat": p_gpi2_hat, "chi2dof_main": chi2dof,
            "err_joint_with_free_noise": err_joint, "err_joint_with_fixed_noise": err_fixed,
            "zz_real_mean": ZZ_REAL_MEAN, "zz_real_std": ZZ_REAL_STD,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
