#!/usr/bin/env python3
"""
task36_joint_schmidt_frame.py -- iteration 36. JOINT SCHMIDT-FRAME
ESTIMATOR: instead of fitting 21 independent 6-dim states (one per
measurement slot, each its own nonconvex `fit_pure_state` optimization),
fit ONE shared 6x6 orthonormal frame U (15 free parameters via a
skew-symmetric Lie-algebra parameterization, orthogonality exact by
construction) jointly against ALL 21 slots' real measured data at once.
Every slot's target vector is then DERIVED from U via the exact known
algebraic combination (u_n = U[:,n], (u_n+u_m)/sqrt2 = (U[:,n]+U[:,m])/
sqrt2, etc. -- the SAME formula `target_coeff_vector` already uses for
the noiseless ideal target, just applied to the FITTED frame instead of
the standard basis).

WHY THIS IS A REAL BIAS-VARIANCE TRADE, NOT A FREE LUNCH -- stated up
front, not discovered after the fact: the shared-frame algebraic identity
is EXACTLY true for the noiseless ideal targets (confirmed directly:
<P>_{(u_n-u_m)/sqrt2} = <P>_{u_n} + <P>_{u_m} - <P>_{(u_n+u_m)/sqrt2},
already exploited by this project's own `build_full_from_a` for the
minus-slots, which are never independently measured). But under REAL
noise, each of the 21 slots is a SEPARATE, INDEPENDENTLY-EXECUTED
circuit -- there is no physical law forcing slot (u_n+u_m)'s actual noisy
deviation to equal the algebraic combination of slots u_n and u_m's own
independent noise deviations. Imposing the shared-frame constraint on
noisy data is therefore a REGULARIZER (trading some potential bias, if
real noise violates the shared structure, for a large reduction in
independent per-slot variance -- 21 independent 5-dof fits = 105
directions collapsed to 15 shared dof), not a mathematically guaranteed
improvement. This is measured directly below, not assumed.

VALIDATION, mandatory before trusting anything on real data (this
project's own established discipline, violated nowhere here):
  1. IDEAL-DATA RECOVERY: fit against noiseless (p=0) exact expectation
     values. Must recover the exact energy near-perfectly with near-zero
     residual chi2 -- if the parameterization or optimizer has a bug,
     this is where it shows up cheaply, before touching real data.
  2. ADVERSARIAL REJECTION: fit against SHUFFLED real data (same real
     numbers, wrong labels). A working joint-frame fit should NOT explain
     scrambled data well -- chi2/dof should be large, confirming the
     model isn't just curve-fitting noise into looking like anything.

Then: a REAL, paired bootstrap comparison (same trial produces both the
standard 21-independent-fit energy AND the joint-frame energy, from Task
31C's real checkpoint, no new submissions) -- bias, std, MSE, chi2/dof,
and explicit outlier-rate reporting (Task 35E's cross-fitting looked fine
on paper and then produced a 6% catastrophic-outlier rate; this is
checked here with the same rigor, not assumed away).

Run:
    python vqe/task36_joint_schmidt_frame.py
"""
import os
import sys
import json
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts

K = 6
N_THETA = K * (K - 1) // 2  # 15
SHOTS = 100_000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
N_BOOT = 80
N_RESTARTS = 12
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task36_joint_schmidt_frame_results.json")


# ---------- parameterization ----------

def skew_from_theta(theta, K):
    A = np.zeros((K, K))
    idx = 0
    for i in range(K):
        for j in range(i + 1, K):
            A[i, j] = theta[idx]
            A[j, i] = -theta[idx]
            idx += 1
    return A


def U_from_theta(theta, U0, K):
    A = skew_from_theta(theta, K)
    return U0 @ expm(A)


def slot_vector(U, name, K):
    coeff = target_coeff_vector(name, K)
    return U @ coeff


def build_full_from_frame(U, P_S, K, non_id_labels, kept):
    """Every slot's reconstructed label value, direct bilinear form of
    that slot's frame-derived vector -- diag AND plus slots alike, no
    separate synthesis step needed (the minus-slot algebraic identity is
    handled automatically since it's the SAME bilinear form applied to
    the algebraically-combined vector)."""
    full = {}
    for name in kept:
        v = slot_vector(U, name, K)
        full[name] = {l: float(np.real(v @ P_S[l] @ v)) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
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


# ---------- joint fit ----------

def procrustes_init(a_by_name_diag, K):
    """Nearest-orthogonal-matrix projection (SVD) of the 6 diagonal
    slots' independently-fitted vectors -- a physically-motivated,
    non-target-peeking starting point for the joint optimizer."""
    U_tilde = np.column_stack([a_by_name_diag[f"u_{n}"] for n in range(K)])
    W, S, Vt = np.linalg.svd(U_tilde)
    return W @ Vt


def joint_residuals(theta, U0, P_S, K, kept, non_id_labels, blended_by_slot, weight_by_slot):
    U = U_from_theta(theta, U0, K)
    res = []
    for name in kept:
        v = slot_vector(U, name, K)
        blended = blended_by_slot[name]
        for l, y in blended.items():
            pred = float(np.real(v @ P_S[l] @ v))
            w = weight_by_slot[name].get(l, 1.0)
            res.append(np.sqrt(w) * (pred - y))
    return np.array(res)


def fit_joint_frame(U0, P_S, K, kept, non_id_labels, blended_by_slot, weight_by_slot, rng, n_restarts=N_RESTARTS):
    best = None
    scales = [0.0, 0.02, 0.05, 0.08, 0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9][:n_restarts]
    inits = [np.zeros(N_THETA) if s == 0.0 else rng.normal(0, s, N_THETA) for s in scales]
    for theta0 in inits:
        try:
            sol = least_squares(joint_residuals, theta0, args=(U0, P_S, K, kept, non_id_labels, blended_by_slot, weight_by_slot),
                                 method="lm", max_nfev=2000)
        except Exception:
            continue
        cost = float(2 * sol.cost)  # sum of squared residuals
        if best is None or cost < best[0]:
            best = (cost, sol.x)
    cost, theta_hat = best
    U_hat = U_from_theta(theta_hat, U0, K)
    n_resid = sum(len(blended_by_slot[name]) for name in kept)
    dof = max(1, n_resid - N_THETA)
    return U_hat, cost, cost / dof


def main():
    print("\n" + "=" * 96)
    print("  task36_joint_schmidt_frame.py -- joint 15-parameter Schmidt-frame fit vs 21 independent fits")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    # sorted(): fixes a real, diagnosed bug -- p["alpha_labels"] iteration order is not
    # guaranteed stable across process launches (Python hash randomization), which was
    # silently changing residual ORDER fed into the nonconvex joint-frame optimizer and
    # causing genuine cross-process convergence to different local optima (verified: with
    # PYTHONHASHSEED randomized, 5/10 repeated identical-seed runs landed in a basin with
    # 18-24x worse chi2/dof; with sorted(), or equivalently PYTHONHASHSEED=0, 6/6 identical).
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    # ================= VALIDATION 1: IDEAL-DATA RECOVERY =================
    print("\n  -- VALIDATION 1: ideal-data (noiseless) recovery --")
    ideal_blended = {}
    for name in kept:
        v_ideal = slot_vector(np.eye(K), name, K)  # standard basis frame = exact ideal targets
        ideal_blended[name] = {l: float(np.real(v_ideal @ P_S[l] @ v_ideal)) for l in non_id_labels}
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    rng_val = np.random.default_rng(11)
    U_hat_ideal, cost_ideal, chi2dof_ideal = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels,
                                                              ideal_blended, weight_unit, rng_val)
    full_ideal = build_full_from_frame(U_hat_ideal, P_S, K, non_id_labels, kept)
    E_ideal, err_ideal = energy_and_err(p, full_ideal, K)
    print(f"    recovered energy err vs exact: {err_ideal:.6f} kcal/mol  (should be ~0)")
    print(f"    chi2/dof: {chi2dof_ideal:.8f}  (should be ~0)")
    ideal_ok = err_ideal < 0.01 and chi2dof_ideal < 1e-6
    print(f"    IDEAL-DATA RECOVERY: {'PASS' if ideal_ok else 'FAIL -- STOP, do not trust anything below'}")
    if not ideal_ok:
        raise RuntimeError("Joint-frame fit fails to recover the noiseless ideal case -- a real bug, not a physics finding")

    # ================= VALIDATION 2: ADVERSARIAL (SHUFFLED REAL DATA) =================
    print("\n  -- VALIDATION 2: adversarial rejection (real data, labels shuffled within each slot) --")
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
    print(f"  loaded Task 31C's real checkpoint: {n_mc_total} MC draws/label, {len(kept)} slots")

    rng_adv = np.random.default_rng(2026)

    def real_blended_for_slot(name, n_mc, rng):
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

    real_blended = {name: real_blended_for_slot(name, n_mc_total, rng_adv) for name in kept}
    U0_adv = procrustes_init({n: target_coeff_vector(n, K) for n in diag}, K)  # trivial init, not target-peeking
    shuffled_blended = {}
    for name in kept:
        labels_here = list(real_blended[name].keys())
        vals_here = list(real_blended[name].values())
        perm = rng_adv.permutation(len(vals_here))
        shuffled_blended[name] = {l: vals_here[perm[i]] for i, l in enumerate(labels_here)}
    U_hat_adv, cost_adv, chi2dof_adv = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels,
                                                        shuffled_blended, weight_unit, rng_adv)
    U_hat_real, cost_real, chi2dof_real = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels,
                                                           real_blended, weight_unit, rng_adv)
    print(f"    chi2/dof on REAL data (unshuffled):  {chi2dof_real:.4f}")
    print(f"    chi2/dof on SHUFFLED (adversarial):  {chi2dof_adv:.4f}")
    adv_ok = chi2dof_adv > 3 * chi2dof_real
    print(f"    ADVERSARIAL REJECTION: {'PASS (model correctly rejects scrambled data)' if adv_ok else 'FAIL -- model fits garbage nearly as well as real data, not trustworthy'}")

    print(f"\n  {'STOPPING -- adversarial check failed, joint frame is not a trustworthy diagnostic model' if not adv_ok else 'Both validations passed -- proceeding to real bootstrap comparison'}")
    if not adv_ok:
        with open(RESULTS_PATH, "w") as f:
            json.dump({"ideal_ok": ideal_ok, "adv_ok": adv_ok, "chi2dof_real": chi2dof_real,
                        "chi2dof_adv": chi2dof_adv, "stopped_early": True}, f, indent=2)
        return

    # ================= REAL PAIRED BOOTSTRAP COMPARISON =================
    print(f"\n  -- running {N_BOOT} paired real bootstrap trials (standard 21-fit vs joint 15-param frame) --")
    v0_by_name = {name: target_coeff_vector(name, K) for name in kept}
    with ProcessPoolExecutor(max_workers=2) as ex:
        futures = [ex.submit(_trial_worker, p, K, kept, non_id_labels, P_S, by_name_label, n_mc_total,
                              v0_by_name, weight_unit, 90_000 + i)
                   for i in range(N_BOOT)]
        results = [fut.result() for fut in futures]

    errs_std = np.array([r[0] for r in results])
    errs_joint = np.array([r[1] for r in results])
    chi2dofs = np.array([r[2] for r in results])

    bias_std, std_std = float(np.median(errs_std)), float(errs_std.std(ddof=1))
    bias_joint, std_joint = float(np.median(errs_joint)), float(errs_joint.std(ddof=1))
    mse_std = bias_std ** 2 + std_std ** 2
    mse_joint = bias_joint ** 2 + std_joint ** 2

    print(f"\n  STANDARD (21 independent fits):  median={bias_std:.4f}  std={std_std:.4f}  MSE={mse_std:.4f}")
    print(f"  JOINT (15-param shared frame):   median={bias_joint:.4f}  std={std_joint:.4f}  MSE={mse_joint:.4f}")
    print(f"  chi2/dof across trials: median={np.median(chi2dofs):.4f}  max={chi2dofs.max():.4f}")

    print(f"\n  errs_std sorted:   {np.array2string(np.sort(errs_std), precision=2, max_line_width=200)}")
    print(f"  errs_joint sorted: {np.array2string(np.sort(errs_joint), precision=2, max_line_width=200)}")

    n_outliers_std = int((errs_std > 20).sum())
    n_outliers_joint = int((errs_joint > 20).sum())
    diffs = errs_std - errs_joint
    n_joint_better = int((diffs > 0).sum())

    print(f"\n  outlier rate (>20 kcal/mol): standard {n_outliers_std}/{N_BOOT}, joint {n_outliers_joint}/{N_BOOT}")
    print(f"  paired: joint wins {n_joint_better}/{N_BOOT} trials, mean(std_err-joint_err)={diffs.mean():+.4f}")

    print(f"\n  -- HONEST READ --")
    if mse_joint < mse_std and n_outliers_joint <= n_outliers_std:
        pct = 100 * (mse_std - mse_joint) / mse_std
        print(f"    JOINT SCHMIDT-FRAME WINS: {pct:.1f}% MSE reduction, no worse outlier behavior. "
              f"The shared-frame regularization is a real, net win on this real data.")
    elif mse_joint < mse_std:
        print(f"    Mixed: joint frame reduces MSE but shows a different/worse outlier profile "
              f"({n_outliers_joint} vs {n_outliers_std}) -- a real trade, not a clean win.")
    else:
        pct = 100 * (mse_joint - mse_std) / mse_std
        print(f"    JOINT SCHMIDT-FRAME DOES NOT WIN: {pct:.1f}% worse MSE than the standard 21-independent-"
              f"fit approach. The shared-frame constraint, despite fitting real data far better than "
              f"scrambled data (Validation 2), does not translate into a lower-error energy estimator here.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "ideal_ok": ideal_ok, "adv_ok": adv_ok, "chi2dof_real": chi2dof_real, "chi2dof_adv": chi2dof_adv,
            "N_BOOT": N_BOOT, "errs_std": errs_std.tolist(), "errs_joint": errs_joint.tolist(),
            "chi2dofs": chi2dofs.tolist(),
            "bias_std": bias_std, "std_std": std_std, "mse_std": mse_std,
            "bias_joint": bias_joint, "std_joint": std_joint, "mse_joint": mse_joint,
            "n_outliers_std": n_outliers_std, "n_outliers_joint": n_outliers_joint,
            "n_joint_better": n_joint_better,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


def _trial_worker(p, K, kept, non_id_labels, P_S, by_name_label, n_mc_total, v0_by_name, weight_unit, seed):
    rng = np.random.default_rng(seed)

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

    blended_by_slot = {name: real_blended_for_slot(name, n_mc_total) for name in kept}

    # -- STANDARD: 21 independent fit_pure_state calls, same discipline as every prior bootstrap worker --
    a_by_name = {}
    for name in kept:
        a_hat, _ = fit_pure_state(P_S, blended_by_slot[name], {l: 1.0 for l in blended_by_slot[name]}, K,
                                   v0_by_name[name], seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full_std = {}
    for name in kept:
        v = a_by_name[name]
        full_std[name] = {l: float(np.real(v @ P_S[l] @ v)) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            synth_minus = {}
            for l in non_id_labels:
                synth_minus[l] = full_std[un][l] + full_std[um][l] - full_std[pl][l]
            full_std[f"(u{n}-u{m})"] = synth_minus
    _, err_std = energy_and_err(p, full_std, K)

    # -- JOINT: shared 15-param frame, Procrustes-initialized from the SAME independent fits above
    # (reusing them purely as a numerically-stable starting point, not as the estimate itself) --
    U0 = procrustes_init({n: a_by_name[n] for n in [f"u_{i}" for i in range(K)]}, K)
    U_hat, cost, chi2dof = fit_joint_frame(U0, P_S, K, kept, non_id_labels, blended_by_slot, weight_unit, rng)
    full_joint = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    _, err_joint = energy_and_err(p, full_joint, K)

    return err_std, err_joint, chi2dof


if __name__ == "__main__":
    main()
