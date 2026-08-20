#!/usr/bin/env python3
"""
task35e_cross_fitted_manifold.py -- iteration 35, Task E. CROSS-FITTED
H4 manifold reconstruction, targeting the same-sample nonlinear bias
mechanism the proposal identified as a plausible cause of the project's
own repeated 0.115-vs-0.317-vs-0.438 reproducibility problem.

THE MECHANISM, made concrete for this codebase: `build_full_from_a`
(task29c_manifold_estimator.py) reconstructs every label's expectation
as a SELF bilinear form of ONE fitted state vector, `expect(a,l) =
a @ P_S[l] @ a`, where `a` is itself a NONLINEAR (L-BFGS-fit) function of
the SAME noisy data used to evaluate that bilinear form. For a nonlinear
plug-in estimator like this, E[a_hat^T P a_hat] = mu^T P mu + Tr(P*Sigma)
in general -- an extra bias term from the estimator's own sampling
covariance Sigma that does NOT vanish just because the fit is unbiased in
some other sense, and does not average away across repeated real
submissions the way ordinary sampling noise does. This applies
independently to EACH of the 21 kept slots (6 diagonal + 15 pair-plus),
since each has its own independent `fit_pure_state` call and its own self
bilinear reconstruction.

THE FIX: split each slot's data into two independent halves (A, B), fit
a_hat_A and a_hat_B SEPARATELY, then reconstruct using the SYMMETRIZED
CROSS bilinear form:
    expect_cross(l) = (a_hat_A @ P_S[l] @ a_hat_B + a_hat_B @ P_S[l] @ a_hat_A) / 2
For independent A, B: E[a_hat_A^T P a_hat_B] = mu^T P mu exactly (no
covariance cross-term, since Cov(a_hat_A, a_hat_B)=0 by construction) --
removing the Tr(P*Sigma) bias term, at the cost of each half-fit having
roughly 2x the per-half sampling variance (half the data). Net effect on
bias/variance is an empirical question, not assumed either way -- this
task measures it directly against the SAME-TRIAL standard (same-sample,
full-data) reconstruction for a fair, paired comparison.

DATA: Task 31C's real 4,368-circuit checkpoint (16 MC draws/label, all
21 kept slots), split 8+8 per slot per trial. No new submissions.

Run:
    python vqe/task35e_cross_fitted_manifold.py
"""
import os
import sys
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
N_BOOT = 32
EXACT_ENERGY_TARGET = 0.5
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task35e_cross_fitted_manifold_results.json")


def expect(a, P):
    return float(np.real(a @ P @ a))


def expect_cross(a1, a2, P):
    return float(np.real(a1 @ P @ a2 + a2 @ P @ a1) / 2)


def build_full_standard(a_by_name, P_S, diag, K, non_id_labels):
    full = {}
    for name in diag:
        a = a_by_name[name]
        full[name] = {l: expect(a, P_S[l]) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            a_pl = a_by_name[pl]
            full[pl] = {l: expect(a_pl, P_S[l]) for l in non_id_labels}
            synth_minus = {}
            for l in non_id_labels:
                synth_minus[l] = full[un][l] + full[um][l] - full[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def build_full_cross(a_by_name_A, a_by_name_B, P_S, diag, K, non_id_labels):
    full = {}
    for name in diag:
        aA, aB = a_by_name_A[name], a_by_name_B[name]
        full[name] = {l: expect_cross(aA, aB, P_S[l]) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            aA, aB = a_by_name_A[pl], a_by_name_B[pl]
            full[pl] = {l: expect_cross(aA, aB, P_S[l]) for l in non_id_labels}
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


def _trial_worker(p, non_id_labels, diag, kept, P_S, by_name_label, n_mc_total, seed):
    rng = np.random.default_rng(seed)
    a_full, a_A, a_B = {}, {}, {}
    for name in kept:
        blended_full, blended_A, blended_B = {}, {}, {}
        for l in non_id_labels:
            key = (name, l)
            if key not in by_name_label:
                continue
            entries = by_name_label[key]
            n = len(entries)
            half = n // 2
            perm = rng.permutation(n)
            idx_A_pool, idx_B_pool = perm[:half], perm[half:2 * half]

            def blended_mean(idx_pool, size):
                draw_idx = idx_pool[rng.integers(0, len(idx_pool), size=size)] if len(idx_pool) > 0 else np.array([], dtype=int)
                vals = []
                for idx in draw_idx:
                    sign, gamma, counts = entries[idx]
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    m = expectation_from_counts(resampled, l)
                    vals.append(sign * gamma * m)
                return float(np.mean(vals)) if vals else 0.0

            blended_full[l] = max(-1.0, min(1.0, blended_mean(np.arange(n), n_mc_total)))
            blended_A[l] = max(-1.0, min(1.0, blended_mean(idx_A_pool, half)))
            blended_B[l] = max(-1.0, min(1.0, blended_mean(idx_B_pool, half)))

        v0 = target_coeff_vector(name, K)
        a_full[name], _ = fit_pure_state(P_S, blended_full, {l: 1.0 for l in blended_full}, K, v0,
                                          seed=int(rng.integers(0, 2**31)))
        a_A[name], _ = fit_pure_state(P_S, blended_A, {l: 1.0 for l in blended_A}, K, v0,
                                       seed=int(rng.integers(0, 2**31)))
        a_B[name], _ = fit_pure_state(P_S, blended_B, {l: 1.0 for l in blended_B}, K, v0,
                                       seed=int(rng.integers(0, 2**31)))

    full_std = build_full_standard(a_full, P_S, diag, K, non_id_labels)
    full_crx = build_full_cross(a_A, a_B, P_S, diag, K, non_id_labels)
    _, err_std = energy_and_err(p, full_std, K)
    _, err_crx = energy_and_err(p, full_crx, K)
    return err_std, err_crx


def main():
    print("\n" + "=" * 96)
    print("  task35e_cross_fitted_manifold.py -- cross-fitted vs same-sample H4 manifold reconstruction")
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
    by_name_label = {}
    for (name, group, draw, sign), counts, gamma in zip(tags, counts_list, gamma_per_circuit):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma, counts))
    n_mc_total = min(len(v) for v in by_name_label.values())
    print(f"  loaded Task 31C's real checkpoint: {n_mc_total} MC draws/label, {len(kept)} slots "
          f"(split {n_mc_total//2}+{n_mc_total//2} per trial for cross-fitting)")

    print(f"\n  -- running {N_BOOT} paired trials (SAME trial produces both same-sample and cross-fitted "
          f"energy, for a fair paired comparison) --")
    with ProcessPoolExecutor(max_workers=2) as ex:
        futures = [ex.submit(_trial_worker, p, non_id_labels, diag, kept, P_S, by_name_label, n_mc_total,
                              80_000 + i)
                   for i in range(N_BOOT)]
        results = [fut.result() for fut in futures]
    errs_std = np.array([r[0] for r in results])
    errs_crx = np.array([r[1] for r in results])

    bias_std, std_std = float(np.median(errs_std)), float(errs_std.std(ddof=1))
    bias_crx, std_crx = float(np.median(errs_crx)), float(errs_crx.std(ddof=1))
    mse_std = bias_std ** 2 + std_std ** 2
    mse_crx = bias_crx ** 2 + std_crx ** 2

    print(f"\n  SAME-SAMPLE (standard, full 16 draws, self bilinear):  "
          f"median={bias_std:.4f}  std={std_std:.4f}  MSE={mse_std:.4f}")
    print(f"  CROSS-FITTED (8+8 split, cross bilinear):               "
          f"median={bias_crx:.4f}  std={std_crx:.4f}  MSE={mse_crx:.4f}")

    # paired difference: same trial's std-vs-cross difference, to see if cross-fitting moves each
    # trial in a consistent direction (a real paired effect) vs just adding independent noise
    diffs = errs_std - errs_crx
    print(f"\n  PAIRED per-trial difference (same_sample_err - cross_fitted_err): "
          f"mean={diffs.mean():+.4f}  std={diffs.std(ddof=1):.4f}")
    consistent_direction = int((diffs > 0).sum())
    print(f"  cross-fitting had LOWER error in {consistent_direction}/{N_BOOT} trials "
          f"({'a real, consistent effect' if abs(consistent_direction - N_BOOT/2) > N_BOOT*0.15 else 'not a consistent direction -- likely just added independent noise'})")

    print(f"\n  -- HONEST READ --")
    bias_reduced = abs(bias_crx) < abs(bias_std)
    variance_reduced = std_crx < std_std
    mse_reduced = mse_crx < mse_std
    print(f"    bias:     {'REDUCED' if bias_reduced else 'NOT reduced'} ({abs(bias_std):.4f} -> {abs(bias_crx):.4f})")
    print(f"    variance: {'REDUCED' if variance_reduced else 'INCREASED'} ({std_std:.4f} -> {std_crx:.4f})")
    print(f"    MSE:      {'REDUCED' if mse_reduced else 'INCREASED'} ({mse_std:.4f} -> {mse_crx:.4f})")
    if mse_reduced:
        pct = 100 * (mse_std - mse_crx) / mse_std
        print(f"    -> Cross-fitting is a real MSE win, {pct:.1f}% reduction. The same-sample bias "
              f"mechanism (Tr(P*Sigma)) is real and worth removing here.")
    else:
        pct = 100 * (mse_crx - mse_std) / mse_std
        print(f"    -> Cross-fitting does NOT win on MSE here ({pct:.1f}% worse) -- the bias reduction "
              f"(if any) is outweighed by each half-fit's higher variance from using half the data. "
              f"The same-sample bias mechanism may be real but small relative to this pipeline's "
              f"dominant PEC-MC variance (Task 32B), which cross-fitting does not address at all.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N_BOOT": N_BOOT, "n_mc_total": n_mc_total,
            "errs_std": errs_std.tolist(), "errs_crx": errs_crx.tolist(),
            "bias_std": bias_std, "std_std": std_std, "mse_std": mse_std,
            "bias_crx": bias_crx, "std_crx": std_crx, "mse_crx": mse_crx,
            "paired_diff_mean": float(diffs.mean()), "paired_diff_std": float(diffs.std(ddof=1)),
            "n_trials_cross_better": consistent_direction,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
