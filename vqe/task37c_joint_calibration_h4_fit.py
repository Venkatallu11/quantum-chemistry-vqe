#!/usr/bin/env python3
"""
task37c_joint_calibration_h4_fit.py -- iteration 37, Task C. The actual
simultaneous calibration+H4 fit the proposal called for: jointly estimate
the 15-parameter Schmidt frame (Task 36) AND the 5-parameter H4-native
noise model (Task 37B: p_zz, p_gpi2, delta_zz, delta_gpi2, p_readout),
against REAL RAW (pre-PEC) H4 hardware data, using Task 37c's own
extended, self-verified forward model (`task37c_extended_forward_model.
analytic_A_and_B_5param` + `readout_attenuation`).

IDENTIFIABILITY SAFEGUARD, the proposal's central requirement, THREE
separate data roles, never mixed:
  1. CALIBRATION data: Task 31A/30B's real ZZ/GPi measurements -- already
     summarized as priors in Task 37B, used here ONLY as Gaussian penalty
     terms (mean/std), never as raw circuits refit here.
  2. H4-TRAINING data: a stratified per-slot SPLIT of the real H4 raw
     (slot, label) measurements -- ~70% of labels in every kept slot,
     fixed seed. This is what the joint optimizer actually sees.
  3. H4-VALIDATION data: the remaining ~30% of labels per slot, HELD OUT
     of the fit entirely. Used only to check chi2/dof AFTER fitting
     (generalization / overfitting check) -- this is the safeguard
     against the identifiability collapse the proposal warned about
     (Phase 1A fit ALL real data with no held-out check at all).
The known EXACT H4 energy is used for exactly ONE thing, at the very end,
after every fitting/model-selection decision above is already frozen: an
informational-only report of the resulting energy error. It is NEVER
compared, differenced against, or fed back into anything upstream of that
final print -- no target leakage.

Run:
    PYTHONHASHSEED=0 python vqe/task37c_joint_calibration_h4_fit.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import U_from_theta, slot_vector, build_full_from_frame, N_THETA
from task37b_h4_noise_model import PARAM_NAMES, PRIORS, GPI_REAL_MEAN, GPI_REAL_STD, clip_to_bounds
from task37c_extended_forward_model import analytic_A_and_B_5param, readout_attenuation, hamming_weight
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts, stable_seed

K = 6
GATE_NAME = "zz"
SHOTS = 100_000
N_NOISE_PARAMS = len(PARAM_NAMES)  # 5
VAL_FRACTION = 0.30
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37c_joint_calibration_h4_fit_results.json")

_GLOBAL_AB_CACHE = {}  # keyed on (name, rounded p_zz, p_gpi2, delta_zz, delta_gpi2) -- same performance
# fix Phase 1A needed, extended to all 4 noise-forward-model params (readout is applied afterward,
# cheaply, so it does NOT need to be part of this cache key)


def theta_from_params(params):
    return {name: float(params[N_THETA + i]) for i, name in enumerate(PARAM_NAMES)}


def clamp_noise(params):
    raw = params[N_THETA:N_THETA + N_NOISE_PARAMS]
    return clip_to_bounds(raw)


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def predict_corrected(name, l, raw_val, noise, fixed_solutions, non_id_here_cache):
    p_zz, p_gpi2, delta_zz, delta_gpi2, p_readout = noise
    cache_key = (name, round(p_zz, 6), round(p_gpi2, 6), round(delta_zz, 6), round(delta_gpi2, 6))
    if cache_key not in _GLOBAL_AB_CACHE:
        labels_here = non_id_here_cache[name]
        A, B = analytic_A_and_B_5param(fixed_solutions[name]["angles"], GATE_NAME, p_zz, GPI_REAL_MEAN,
                                        p_gpi2, delta_zz, delta_gpi2, labels_here)
        _GLOBAL_AB_CACHE[cache_key] = (A, B)
    A, B = _GLOBAL_AB_CACHE[cache_key]
    ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
    ro = readout_attenuation(l, p_readout)
    ro = ro if abs(ro) > 1e-6 else 1.0
    m_corrected = max(-1.0, min(1.0, (raw_val / ro) * ratio))
    return m_corrected


def joint_residuals(params, U0, P_S, K, kept, fixed_solutions, data_by_slot, non_id_here_cache,
                     weight_by_slot, prior_lambda):
    theta_state = params[:N_THETA]
    noise = clamp_noise(params)
    U = U_from_theta(theta_state, U0, K)

    res = []
    for name in kept:
        v = slot_vector(U, name, K)
        for l, m_raw in data_by_slot[name].items():
            m_corrected = predict_corrected(name, l, m_raw, noise, fixed_solutions, non_id_here_cache)
            pred = float(np.real(v @ P_S[l] @ v))
            w = weight_by_slot[name].get(l, 1.0)
            res.append(np.sqrt(w) * (m_corrected - pred))

    # -- calibration priors: real for p_zz, weak-but-present for the other 4 (Task 37B) --
    for i, name in enumerate(PARAM_NAMES):
        pr = PRIORS[name]
        res.append(np.sqrt(prior_lambda) * (noise[i] - pr["mean"]) / pr["std"])
    return np.array(res)


def fit_joint(U0, P_S, K, kept, fixed_solutions, data_by_slot, non_id_here_cache, weight_by_slot,
              prior_lambda, rng, n_restarts=8, noise0_options=None):
    best = None
    scales = [0.0, 0.02, 0.05, 0.08, 0.12, 0.15, 0.2, 0.3][:n_restarts]
    if noise0_options is None:
        noise0_options = [tuple(PRIORS[n]["mean"] for n in PARAM_NAMES)]
    for i, s in enumerate(scales):
        theta0 = np.zeros(N_THETA) if s == 0.0 else rng.normal(0, s, N_THETA)
        noise0 = np.array(noise0_options[i % len(noise0_options)])
        params0 = np.concatenate([theta0, noise0])
        try:
            sol = least_squares(joint_residuals, params0,
                                 args=(U0, P_S, K, kept, fixed_solutions, data_by_slot, non_id_here_cache,
                                       weight_by_slot, prior_lambda),
                                 method="lm", max_nfev=4000)
        except Exception:
            continue
        cost = float(2 * sol.cost)
        if best is None or cost < best[0]:
            best = (cost, sol.x)
    cost, params_hat = best
    theta_hat = params_hat[:N_THETA]
    noise_hat = clamp_noise(params_hat)
    U_hat = U_from_theta(theta_hat, U0, K)
    n_resid = sum(len(data_by_slot[name]) for name in kept) + N_NOISE_PARAMS
    dof = max(1, n_resid - (N_THETA + N_NOISE_PARAMS))
    return U_hat, noise_hat, cost, cost / dof


def chi2dof_on_holdout(U_hat, noise_hat, P_S, K, kept, fixed_solutions, holdout_by_slot,
                        non_id_here_cache, weight_by_slot):
    """Evaluates the ALREADY-FITTED (U_hat, noise_hat) against data it never saw during
    optimization -- no fitting happens here, pure evaluation."""
    total, n = 0.0, 0
    for name in kept:
        v = slot_vector(U_hat, name, K)
        for l, m_raw in holdout_by_slot[name].items():
            m_corrected = predict_corrected(name, l, m_raw, noise_hat, fixed_solutions, non_id_here_cache)
            pred = float(np.real(v @ P_S[l] @ v))
            w = weight_by_slot[name].get(l, 1.0)
            total += w * (m_corrected - pred) ** 2
            n += 1
    return total / max(1, n), n


def load_real_blended(kept, n_seeds=8):
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    tags = raw_ck["tags"]["forte-1"]
    counts_list = raw_ck["counts"]["forte-1"]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)
    blended = {name: {} for name in kept}
    non_id_here_cache = {}
    for name in kept:
        labels_here = [l for group_t in per_name[name] for l in group_t]
        non_id_here_cache[name] = labels_here
        for group_t, counts in per_name[name].items():
            for seed in range(n_seeds):
                rng = np.random.default_rng(stable_seed("t37c", name, seed))
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    m = expectation_from_counts(resampled, l)
                    blended[name].setdefault(l, []).append(m)
    for name in kept:
        for l in blended[name]:
            blended[name][l] = float(np.mean(blended[name][l]))
    return blended, non_id_here_cache


def split_train_val(blended, kept, val_fraction, seed):
    """Stratified PER-SLOT split -- every slot contributes to both sets,
    so the 15-parameter frame stays constrained by training data alone in
    every direction, while validation still covers every slot (a
    slot-level split would starve some frame directions entirely)."""
    rng = np.random.default_rng(seed)
    train, val = {name: {} for name in kept}, {name: {} for name in kept}
    for name in kept:
        labels = sorted(blended[name].keys())
        n_val = max(1, round(len(labels) * val_fraction)) if len(labels) > 1 else 0
        val_labels = set(rng.choice(labels, size=n_val, replace=False)) if n_val else set()
        for l in labels:
            (val if l in val_labels else train)[name][l] = blended[name][l]
    return train, val


def main():
    print("\n" + "=" * 96)
    print("  task37c_joint_calibration_h4_fit.py -- joint 15+5-param fit, train/validation SPLIT, no target leakage")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {} for name in kept}  # filled per-label below (uniform=1.0 via .get default)

    real_blended, non_id_here_cache = load_real_blended(kept)
    train_blended, val_blended = split_train_val(real_blended, kept, VAL_FRACTION, seed=37)
    n_train = sum(len(train_blended[n]) for n in kept)
    n_val = sum(len(val_blended[n]) for n in kept)
    print(f"\n  data split: {n_train} training residuals, {n_val} held-out validation residuals "
          f"({100*n_val/(n_train+n_val):.0f}% held out, stratified per slot, seed=37)")

    # ================= VALIDATION 1: IDEAL-DATA RECOVERY =================
    print("\n  -- VALIDATION 1: ideal-data (noiseless raw, p=0) recovery, full data (no split needed) --")
    ideal_blended = {}
    for name in kept:
        v_ideal = slot_vector(np.eye(K), name, K)
        ideal_blended[name] = {l: float(np.real(v_ideal @ P_S[l] @ v_ideal)) for l in non_id_here_cache[name]}
    rng_val = np.random.default_rng(11)
    noise0_ideal = [(1e-5,) * N_NOISE_PARAMS, tuple(PRIORS[n]["mean"] for n in PARAM_NAMES),
                     (0.005, 0.02, 0.0, 0.0, 0.0), (0.02, 0.04, 0.01, -0.01, 0.005)]
    U_hat_i, noise_i, cost_i, chi2dof_i = fit_joint(
        np.eye(K), P_S, K, kept, fixed_solutions, ideal_blended, non_id_here_cache, weight_unit,
        prior_lambda=0.0, rng=rng_val, noise0_options=noise0_ideal)
    full_i = build_full_from_frame(U_hat_i, P_S, K, sorted(p["alpha_labels"]), kept)
    E_i, err_i = energy_and_err(p, full_i, K)
    print(f"    recovered energy err vs exact: {err_i:.6f} kcal/mol  noise_hat={dict(zip(PARAM_NAMES, noise_i.round(6)))}")
    ideal_ok = err_i < 0.05
    print(f"    IDEAL-DATA RECOVERY: {'PASS' if ideal_ok else 'FAIL -- STOP'}")
    if not ideal_ok:
        raise RuntimeError("Task 37C fails ideal-data recovery -- a real bug, stop before trusting real data")

    # ================= VALIDATION 2: ADVERSARIAL (SHUFFLED TRAINING DATA) =================
    print("\n  -- VALIDATION 2: adversarial rejection (TRAINING data only, labels shuffled within each slot) --")
    rng_adv = np.random.default_rng(2026)
    shuffled_blended = {}
    for name in kept:
        labels_here = list(train_blended[name].keys())
        vals_here = list(train_blended[name].values())
        perm = rng_adv.permutation(len(vals_here))
        shuffled_blended[name] = {l: vals_here[perm[i]] for i, l in enumerate(labels_here)}
    U0_std = np.eye(K)
    zz_mean = PRIORS["p_zz"]["mean"]
    noise0_diverse = [tuple(PRIORS[n]["mean"] for n in PARAM_NAMES), (1e-5,) * N_NOISE_PARAMS,
                       (0.02, 0.04, 0.01, -0.01, 0.01), (0.01, 0.02, -0.005, 0.005, 0.002),
                       (zz_mean, 0.045, 0.0, 0.0, 0.015),
                       (zz_mean, 0.01, 0.02, -0.02, 0.0), (0.012, 0.03, -0.01, 0.01, 0.008),
                       (zz_mean, 1e-5, 0.0, 0.0, 0.0)]
    _, _, cost_real, chi2dof_real = fit_joint(U0_std, P_S, K, kept, fixed_solutions, train_blended,
                                                non_id_here_cache, weight_unit, prior_lambda=1.0,
                                                rng=rng_adv, noise0_options=noise0_diverse)
    _, _, cost_adv, chi2dof_adv = fit_joint(U0_std, P_S, K, kept, fixed_solutions, shuffled_blended,
                                              non_id_here_cache, weight_unit, prior_lambda=1.0,
                                              rng=rng_adv, noise0_options=noise0_diverse)
    print(f"    chi2/dof on TRAINING data (unshuffled): {chi2dof_real:.4f}")
    print(f"    chi2/dof on SHUFFLED (adversarial):     {chi2dof_adv:.4f}")
    adv_ok = chi2dof_adv > 3 * chi2dof_real
    print(f"    ADVERSARIAL REJECTION: {'PASS' if adv_ok else 'FAIL -- STOP, not trustworthy'}")
    if not adv_ok:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"ideal_ok": bool(ideal_ok), "adv_ok": bool(adv_ok),
                       "chi2dof_train": chi2dof_real, "chi2dof_adv": chi2dof_adv,
                       "stopped_early": True}, f, indent=2)
        return

    # ================= MAIN FIT: TRAINING data only =================
    print("\n  -- MAIN FIT: joint 15+5-param fit on TRAINING data only (validation held out) --")
    U_hat, noise_hat, cost, chi2dof_train = fit_joint(
        U0_std, P_S, K, kept, fixed_solutions, train_blended, non_id_here_cache, weight_unit,
        prior_lambda=1.0, rng=np.random.default_rng(37), n_restarts=8, noise0_options=noise0_diverse)
    noise_dict = dict(zip(PARAM_NAMES, noise_hat))
    print(f"    fitted noise params:")
    for name in PARAM_NAMES:
        pr = PRIORS[name]
        z = (noise_dict[name] - pr["mean"]) / pr["std"]
        print(f"      {name:<12} = {noise_dict[name]:+.6f}   (prior {pr['mean']:.5f}+/-{pr['std']:.5f}, "
              f"z={z:+.2f}, source: {pr['source']})")
    print(f"    chi2/dof (TRAINING): {chi2dof_train:.4f}")

    # ================= GENERALIZATION CHECK: held-out VALIDATION data =================
    print("\n  -- GENERALIZATION CHECK: chi2/dof on VALIDATION data the fit never saw --")
    chi2dof_val, n_val_used = chi2dof_on_holdout(U_hat, noise_hat, P_S, K, kept, fixed_solutions,
                                                   val_blended, non_id_here_cache, weight_unit)
    ratio_val_train = chi2dof_val / chi2dof_train if chi2dof_train > 1e-12 else float("inf")
    print(f"    chi2/dof (VALIDATION, n={n_val_used}): {chi2dof_val:.4f}")
    print(f"    validation/training chi2 ratio: {ratio_val_train:.2f}x  "
          f"({'no strong overfitting signal' if ratio_val_train < 3 else 'POSSIBLE OVERFITTING -- validation much worse than training'})")

    # ================= INFORMATIONAL ONLY: energy error vs EXACT (never used above) =================
    print("\n  -- INFORMATIONAL ONLY (NOT used in any decision above -- no target leakage): energy vs exact --")
    full_joint = build_full_from_frame(U_hat, P_S, K, sorted(p["alpha_labels"]), kept)
    E_joint, err_joint = energy_and_err(p, full_joint, K)
    print(f"    resulting H4 energy error vs exact: {err_joint:.4f} kcal/mol")
    print(f"    (for context -- Task 36's joint-frame-only result on the same real dataset, and Phase 1A's "
          f"2-noise-param result, are in RESEARCH_LEDGER.md; not recomputed here to avoid re-deciding anything "
          f"post-hoc based on this number)")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "ideal_ok": ideal_ok, "adv_ok": adv_ok,
            "chi2dof_train_check": chi2dof_real, "chi2dof_adv_check": chi2dof_adv,
            "chi2dof_train": chi2dof_train, "chi2dof_val": chi2dof_val,
            "val_train_ratio": ratio_val_train, "n_train": n_train, "n_val": n_val_used,
            "noise_hat": noise_dict,
            "priors": {n: {"mean": PRIORS[n]["mean"], "std": PRIORS[n]["std"], "source": PRIORS[n]["source"]}
                       for n in PARAM_NAMES},
            "err_vs_exact_informational_only": err_joint,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
