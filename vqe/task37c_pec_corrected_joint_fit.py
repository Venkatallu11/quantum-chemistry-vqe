#!/usr/bin/env python3
"""
task37c_pec_corrected_joint_fit.py -- iteration 37, Task C retry, on
PEC-CORRECTED data instead of raw. Task 37C's first attempt (`task37c_
joint_calibration_h4_fit.py`) fit RAW hardware data and failed adversarial
rejection at ~1.93x (need 3x), essentially unchanged from Phase 1A's
2-param version (~2.04x) -- adding 3 more real noise parameters barely
moved the number, suggesting the bottleneck was raw data's information
content, not noise-model richness. This retry uses Task 31C's real,
ALREADY quasi-probability-PEC-corrected checkpoint instead -- the SAME
data Task 36's frame-only fit passed adversarial rejection on by 52x.

KEY DESIGN CHANGE, disclosed precisely (this is NOT the same forward
model as the raw-data attempt, reused blindly -- that would be
conceptually wrong, see below): Task 31C's checkpoint was built by LITERAL
quasi-probability circuit twirling (each measured circuit carries its own
real `sign`/`gamma` weight from the actual PEC protocol), not the
analytic ratio-correction Phase 1A/37C used. Literal twirling PEC, by
construction, corrects exactly one thing: INCOHERENT (depolarizing-type)
gate error, using whatever p_zz/p_gpi2 calibration was assumed when the
gamma weights were computed. It does NOT, and structurally CANNOT,
correct two things this data therefore still carries in full:
  - coherent per-gate-type angle bias (delta_zz, delta_gpi2) -- twirling
    with an assumed DEPOLARIZING inverse has no effect on a genuinely
    coherent (unitary) rotation error;
  - readout error -- a purely classical, post-circuit effect, outside the
    gate-level PEC algebra entirely.
So the noise model fit here is a SECOND, SMALLER correction layered on
top of Task 31C's own real correction:
    m_final_predicted = (m_pec_corrected_real / (1-2*p_readout)^w)
                         * B(theta)/A(theta)
using the SAME verified `analytic_A_and_B_5param` (Task 37C's own,
unchanged) -- but with THETA'S PRIOR MEANING CHANGED for p_zz/p_gpi2
specifically: here they represent a RESIDUAL/leftover depolarizing
miscalibration on top of Task 31C's own real correction, not the noise
from scratch. p_zz's prior is therefore re-centered at 0 with a std equal
to its OWN real calibration uncertainty (Task 31A: std=0.000124) --
Task 31C's correction is trusted to already be close, so only a small
leftover is plausible, and that plausible scale is itself a real,
measured number, not invented. p_gpi2's prior is left UNCHANGED from
Task 37B (weak, [0,0.05]) -- because the GPi2 value baked into Task 31C's
own gamma weights came from the SAME failed recalibration attempt
(flagged "physical": false, Task 31A), so there is no more reason to
trust a "small residual" framing here than for the from-scratch case.
delta_zz, delta_gpi2, p_readout keep Task 37B's ORIGINAL priors
unchanged in both contexts -- PEC's twirling structurally cannot touch
any of the three regardless of how accurate its own depolarizing
calibration was, so the plausible residual scale for these three is the
SAME whether the data was PEC-corrected first or not.

IDENTIFIABILITY SAFEGUARD (same discipline as the raw-data attempt):
stratified per-slot 70/30 train/validation split of the real (slot,
label) data; adversarial rejection checked on TRAINING data only; a
frame-ONLY baseline (Task 36's own unchanged `fit_joint_frame`) is fit on
the SAME training split for direct comparison, to check honestly whether
adding the 5 residual noise parameters helps, hurts, or is neutral
relative to Task 36's already-working result -- not assumed either way.
Exact H4 energy used ONLY for a final informational-only report, after
every decision above is frozen. No target leakage.

Run:
    PYTHONHASHSEED=0 python vqe/task37c_pec_corrected_joint_fit.py
"""
import os
import sys
import json
import copy
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import (U_from_theta, slot_vector, build_full_from_frame, N_THETA,
                                          fit_joint_frame, CKPT_PATH)
from task37b_h4_noise_model import PARAM_NAMES, PRIORS as PRIORS_FROM_SCRATCH, GPI_REAL_MEAN
from task37c_extended_forward_model import analytic_A_and_B_5param, readout_attenuation
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts

K = 6
GATE_NAME = "zz"
SHOTS = 100_000
N_NOISE_PARAMS = len(PARAM_NAMES)  # 5
VAL_FRACTION = 0.30
ZZ_REAL_STD = 0.000124  # Task 31A -- reused as the plausible RESIDUAL scale, not invented
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37c_pec_corrected_joint_fit_results.json")

# -- residual priors: p_zz re-centered at 0 (see module docstring), everything else UNCHANGED
# from Task 37B's own from-scratch priors --
PRIORS_RESIDUAL = copy.deepcopy(PRIORS_FROM_SCRATCH)
PRIORS_RESIDUAL["p_zz"] = {"mean": 0.0, "std": ZZ_REAL_STD, "bounds": (-0.001, 0.001),
                            "source": "RESIDUAL (Task 31C's own correction already used the real "
                                      "Task 31A p_zz calibration; std reused directly from that "
                                      "calibration's own measurement uncertainty)"}

_GLOBAL_AB_CACHE = {}


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def clamp_noise(params):
    raw = params[N_THETA:N_THETA + N_NOISE_PARAMS]
    out = np.array(raw, dtype=float)
    for i, name in enumerate(PARAM_NAMES):
        lo, hi = PRIORS_RESIDUAL[name]["bounds"]
        out[i] = min(hi, max(lo, out[i]))
    return out


def predict_corrected(name, l, m_pec, noise, fixed_solutions, non_id_labels):
    p_zz, p_gpi2, delta_zz, delta_gpi2, p_readout = noise
    cache_key = (name, round(p_zz, 6), round(p_gpi2, 6), round(delta_zz, 6), round(delta_gpi2, 6))
    if cache_key not in _GLOBAL_AB_CACHE:
        # p_zz here may be slightly negative (a residual, not a from-scratch probability) --
        # analytic_A_and_B_5param's depolarizing_weights/pec_inverse_weights are evaluated
        # algebraically and remain well-defined for small |p| regardless of sign; only used to
        # compute an analytic RATIO, never as a physically-applied channel.
        A, B = analytic_A_and_B_5param(fixed_solutions[name]["angles"], GATE_NAME, p_zz, GPI_REAL_MEAN,
                                        p_gpi2, delta_zz, delta_gpi2, non_id_labels)
        _GLOBAL_AB_CACHE[cache_key] = (A, B)
    A, B = _GLOBAL_AB_CACHE[cache_key]
    ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
    ro = readout_attenuation(l, p_readout)
    ro = ro if abs(ro) > 1e-6 else 1.0
    return max(-1.0, min(1.0, (m_pec / ro) * ratio))


def joint_residuals(params, U0, P_S, K, kept, non_id_labels, fixed_solutions, data_by_slot, weight_by_slot,
                     prior_lambda):
    theta_state = params[:N_THETA]
    noise = clamp_noise(params)
    U = U_from_theta(theta_state, U0, K)
    res = []
    for name in kept:
        v = slot_vector(U, name, K)
        for l, m_pec in data_by_slot[name].items():
            pred_raw_side = predict_corrected(name, l, m_pec, noise, fixed_solutions, non_id_labels)
            pred_frame_side = float(np.real(v @ P_S[l] @ v))
            w = weight_by_slot[name].get(l, 1.0)
            res.append(np.sqrt(w) * (pred_raw_side - pred_frame_side))
    for i, name in enumerate(PARAM_NAMES):
        pr = PRIORS_RESIDUAL[name]
        res.append(np.sqrt(prior_lambda) * (noise[i] - pr["mean"]) / pr["std"])
    return np.array(res)


def fit_joint(U0, P_S, K, kept, non_id_labels, fixed_solutions, data_by_slot, weight_by_slot, prior_lambda,
              rng, n_restarts=8, noise0_options=None):
    best = None
    scales = [0.0, 0.02, 0.05, 0.08, 0.12, 0.15, 0.2, 0.3][:n_restarts]
    if noise0_options is None:
        noise0_options = [(0.0, PRIORS_RESIDUAL["p_gpi2"]["mean"], 0.0, 0.0, 0.0)]
    for i, s in enumerate(scales):
        theta0 = np.zeros(N_THETA) if s == 0.0 else rng.normal(0, s, N_THETA)
        noise0 = np.array(noise0_options[i % len(noise0_options)])
        params0 = np.concatenate([theta0, noise0])
        try:
            sol = least_squares(joint_residuals, params0,
                                 args=(U0, P_S, K, kept, non_id_labels, fixed_solutions, data_by_slot,
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


def load_real_blended(kept, non_id_labels):
    """Reuses Task 36's own real-checkpoint loading logic exactly
    (literal sign/gamma-weighted PEC-corrected values), unchanged."""
    with open(CKPT_PATH) as f:
        ck = json.load(f)
    tags = [tuple(t) for t in ck["tags"]]
    counts_list = ck["counts"]
    gamma_per_circuit = ck["gamma_per_circuit"]
    by_name_label = {}
    for (name, group, draw, sign), counts, gamma in zip(tags, counts_list, gamma_per_circuit):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma, counts))
    n_mc_total = min(len(v) for v in by_name_label.values())
    rng = np.random.default_rng(2026)

    def real_blended_for_slot(name, n_mc):
        out = {}
        for l in non_id_labels:
            key = (name, l)
            if key not in by_name_label:
                continue
            entries = by_name_label[key]
            draw_idx = rng.integers(0, len(entries), size=n_mc)
            vals = []
            for idx in draw_idx:
                sign, gamma, counts = entries[idx]
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m = expectation_from_counts(resampled, l)
                vals.append(sign * gamma * m)
            out[l] = max(-1.0, min(1.0, float(np.mean(vals))))
        return out

    return {name: real_blended_for_slot(name, n_mc_total) for name in kept}, n_mc_total


def split_train_val(blended, kept, val_fraction, seed):
    rng = np.random.default_rng(seed)
    train, val = {name: {} for name in kept}, {name: {} for name in kept}
    for name in kept:
        labels = sorted(blended[name].keys())
        n_val = max(1, round(len(labels) * val_fraction)) if len(labels) > 1 else 0
        val_labels = set(rng.choice(labels, size=n_val, replace=False)) if n_val else set()
        for l in labels:
            (val if l in val_labels else train)[name][l] = blended[name][l]
    return train, val


def chi2dof_on_holdout(U_hat, noise_hat, P_S, K, kept, non_id_labels, fixed_solutions, holdout_by_slot,
                        weight_by_slot):
    total, n = 0.0, 0
    for name in kept:
        v = slot_vector(U_hat, name, K)
        for l, m_pec in holdout_by_slot[name].items():
            pred_raw_side = predict_corrected(name, l, m_pec, noise_hat, fixed_solutions, non_id_labels)
            pred_frame_side = float(np.real(v @ P_S[l] @ v))
            w = weight_by_slot[name].get(l, 1.0)
            total += w * (pred_raw_side - pred_frame_side) ** 2
            n += 1
    return total / max(1, n), n


def main():
    print("\n" + "=" * 96)
    print("  task37c_pec_corrected_joint_fit.py -- 15+5-param fit on Task 31C's PEC-corrected data")
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
    print(f"\n  loaded Task 31C's real PEC-corrected checkpoint: {n_mc_total} MC draws/label, {len(kept)} slots")
    train_blended, val_blended = split_train_val(real_blended, kept, VAL_FRACTION, seed=37)
    n_train = sum(len(train_blended[n]) for n in kept)
    n_val = sum(len(val_blended[n]) for n in kept)
    print(f"  data split: {n_train} training residuals, {n_val} held-out validation residuals "
          f"({100*n_val/(n_train+n_val):.0f}% held out, stratified per slot, seed=37)")

    # ================= VALIDATION 1: IDEAL-DATA RECOVERY =================
    print("\n  -- VALIDATION 1: ideal-data (noiseless, no residual correction needed) recovery --")
    ideal_blended = {name: {l: float(np.real(slot_vector(np.eye(K), name, K) @ P_S[l] @ slot_vector(np.eye(K), name, K)))
                             for l in non_id_labels} for name in kept}
    rng_val = np.random.default_rng(11)
    noise0_ideal = [(0.0,) * N_NOISE_PARAMS, (0.0, 0.001, 0.0, 0.0, 0.005), (0.0002, 0.01, 0.005, -0.005, 0.002)]
    U_hat_i, noise_i, cost_i, chi2dof_i = fit_joint(
        np.eye(K), P_S, K, kept, non_id_labels, fixed_solutions, ideal_blended, weight_unit,
        prior_lambda=0.0, rng=rng_val, noise0_options=noise0_ideal)
    full_i = build_full_from_frame(U_hat_i, P_S, K, non_id_labels, kept)
    E_i, err_i = energy_and_err(p, full_i, K)
    print(f"    recovered energy err vs exact: {err_i:.6f} kcal/mol  noise_hat={dict(zip(PARAM_NAMES, noise_i.round(6)))}")
    ideal_ok = err_i < 0.05
    print(f"    IDEAL-DATA RECOVERY: {'PASS' if ideal_ok else 'FAIL -- STOP'}")
    if not ideal_ok:
        raise RuntimeError("Task 37C (PEC-corrected) fails ideal-data recovery -- a real bug, stop before trusting real data")

    # ================= VALIDATION 2: ADVERSARIAL (SHUFFLED TRAINING DATA) =================
    print("\n  -- VALIDATION 2: adversarial rejection (TRAINING data only) -- frame-only baseline AND frame+residual-noise --")
    rng_adv = np.random.default_rng(2026)
    shuffled_blended = {}
    for name in kept:
        labels_here = list(train_blended[name].keys())
        vals_here = list(train_blended[name].values())
        perm = rng_adv.permutation(len(vals_here))
        shuffled_blended[name] = {l: vals_here[perm[i]] for i, l in enumerate(labels_here)}
    U0_std = np.eye(K)

    # frame-ONLY baseline, Task 36's own unchanged function -- sanity reference
    _, _, chi2dof_frameonly_real = fit_joint_frame(U0_std, P_S, K, kept, non_id_labels, train_blended,
                                                      weight_unit, rng_adv)
    _, _, chi2dof_frameonly_adv = fit_joint_frame(U0_std, P_S, K, kept, non_id_labels, shuffled_blended,
                                                     weight_unit, rng_adv)
    print(f"    [frame-only baseline, Task 36's own fit] real={chi2dof_frameonly_real:.4f}  "
          f"shuffled={chi2dof_frameonly_adv:.4f}  ratio={chi2dof_frameonly_adv/max(chi2dof_frameonly_real,1e-12):.1f}x")

    noise0_diverse = [(0.0, PRIORS_RESIDUAL["p_gpi2"]["mean"], 0.0, 0.0, 0.0), (0.0,) * N_NOISE_PARAMS,
                       (0.0002, 0.02, 0.005, -0.005, 0.005), (-0.0002, 0.01, -0.005, 0.005, 0.002),
                       (0.0005, 0.04, 0.0, 0.0, 0.01), (-0.0005, 0.005, 0.01, -0.01, 0.0),
                       (0.0003, 0.03, -0.005, 0.005, 0.008), (0.0, 1e-5, 0.0, 0.0, 0.0)]
    _, _, cost_real, chi2dof_real = fit_joint(U0_std, P_S, K, kept, non_id_labels, fixed_solutions, train_blended,
                                                weight_unit, prior_lambda=1.0, rng=rng_adv,
                                                noise0_options=noise0_diverse)
    _, _, cost_adv, chi2dof_adv = fit_joint(U0_std, P_S, K, kept, non_id_labels, fixed_solutions, shuffled_blended,
                                              weight_unit, prior_lambda=1.0, rng=rng_adv,
                                              noise0_options=noise0_diverse)
    print(f"    [frame+residual-noise, this script]     real={chi2dof_real:.4f}  shuffled={chi2dof_adv:.4f}  "
          f"ratio={chi2dof_adv/max(chi2dof_real,1e-12):.1f}x")
    adv_ok = chi2dof_adv > 3 * chi2dof_real
    print(f"    ADVERSARIAL REJECTION: {'PASS' if adv_ok else 'FAIL -- STOP, not trustworthy'}")
    if not adv_ok:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"ideal_ok": bool(ideal_ok), "adv_ok": bool(adv_ok),
                       "chi2dof_frameonly_real": chi2dof_frameonly_real, "chi2dof_frameonly_adv": chi2dof_frameonly_adv,
                       "chi2dof_train": chi2dof_real, "chi2dof_adv": chi2dof_adv, "stopped_early": True}, f, indent=2)
        return

    # ================= MAIN FIT: TRAINING data only =================
    print("\n  -- MAIN FIT: joint 15+5-param fit on TRAINING data only (validation held out) --")
    U_hat, noise_hat, cost, chi2dof_train = fit_joint(
        U0_std, P_S, K, kept, non_id_labels, fixed_solutions, train_blended, weight_unit,
        prior_lambda=1.0, rng=np.random.default_rng(37), n_restarts=8, noise0_options=noise0_diverse)
    noise_dict = dict(zip(PARAM_NAMES, noise_hat))
    print(f"    fitted residual noise params:")
    for name in PARAM_NAMES:
        pr = PRIORS_RESIDUAL[name]
        z = (noise_dict[name] - pr["mean"]) / pr["std"]
        print(f"      {name:<12} = {noise_dict[name]:+.6f}   (prior {pr['mean']:.5f}+/-{pr['std']:.5f}, "
              f"z={z:+.2f}, source: {pr['source']})")
    print(f"    chi2/dof (TRAINING): {chi2dof_train:.4f}")

    # ================= GENERALIZATION CHECK =================
    print("\n  -- GENERALIZATION CHECK: chi2/dof on VALIDATION data the fit never saw --")
    chi2dof_val, n_val_used = chi2dof_on_holdout(U_hat, noise_hat, P_S, K, kept, non_id_labels, fixed_solutions,
                                                   val_blended, weight_unit)
    ratio_val_train = chi2dof_val / chi2dof_train if chi2dof_train > 1e-12 else float("inf")
    print(f"    chi2/dof (VALIDATION, n={n_val_used}): {chi2dof_val:.4f}")
    print(f"    validation/training chi2 ratio: {ratio_val_train:.2f}x  "
          f"({'no strong overfitting signal' if ratio_val_train < 3 else 'POSSIBLE OVERFITTING'})")

    # ================= INFORMATIONAL ONLY =================
    print("\n  -- INFORMATIONAL ONLY (NOT used in any decision above -- no target leakage): energy vs exact --")
    full_joint = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    E_joint, err_joint = energy_and_err(p, full_joint, K)
    print(f"    resulting H4 energy error vs exact: {err_joint:.4f} kcal/mol")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "ideal_ok": bool(ideal_ok), "adv_ok": bool(adv_ok),
            "chi2dof_frameonly_real": chi2dof_frameonly_real, "chi2dof_frameonly_adv": chi2dof_frameonly_adv,
            "chi2dof_train_check": chi2dof_real, "chi2dof_adv_check": chi2dof_adv,
            "chi2dof_train": chi2dof_train, "chi2dof_val": chi2dof_val,
            "val_train_ratio": ratio_val_train, "n_train": n_train, "n_val": n_val_used,
            "noise_hat": {k: float(v) for k, v in noise_dict.items()},
            "priors": {n: {"mean": PRIORS_RESIDUAL[n]["mean"], "std": PRIORS_RESIDUAL[n]["std"],
                            "source": PRIORS_RESIDUAL[n]["source"]} for n in PARAM_NAMES},
            "err_vs_exact_informational_only": err_joint,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
