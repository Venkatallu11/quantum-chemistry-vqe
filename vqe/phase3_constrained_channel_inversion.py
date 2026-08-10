#!/usr/bin/env python3
"""
phase3_constrained_channel_inversion.py — Phase 3 of physics-constrained
reconstruction (continuing past Phase 1's decision-rule ABANDON at the
user's explicit direction, as Phase 2 already did).
============================================================================
IDEA: learn a Pauli TRANSFER matrix M (M_ij = d<P_i>_noisy / d<P_j>_ideal)
from small calibration circuits with classically-known ideal answers,
then invert it -- via CONSTRAINED least squares (this project's own
physicality constraint, not quasi-probability sampling) -- to correct
real target-circuit measurements. Explicitly contrasted with the PEC
already tried in iterations 4/5/9: PEC inverts a LEARNED CHANNEL via
QUASI-PROBABILITY sampling, with exponential sampling overhead in the
channel's own noise strength (gamma_total); this is a constrained LEAST-
SQUARES inversion instead, with no such overhead.

DATA SOURCE -- existing, no new circuits: `calibrate.json` (iteration 9's
own real IonQ calibration submission) already contains real (ideal,
noisy) Pauli-expectation pairs from 40 independent calibration circuits
per model (8 seeds x 5 random-angle draws of the SAME fixed-structure
ansatz, all 36 non-identity labels measured for each) -- exactly the raw
material Phase 3 asks for, already sitting on disk.

A REAL, HONEST FINDING BEFORE BUILDING ANYTHING FURTHER: the 40
calibration circuits' IDEAL Pauli-expectation vectors span only a
RANK-17 subspace of the full 36-dimensional label space (singular values
computed directly: 17 nonzero, 19 exactly zero to numerical precision).
This means a plain, unconstrained fit/inversion of a full 36x36 M is
fundamentally ILL-POSED -- infinite condition number, an entire 19-
dimensional null space where no calibration data constrains the answer
at all. Reported plainly, not hidden, and used to motivate the actual
design choice below rather than papered over with ad hoc ridge
regularization.

DESIGN CHOICE THIS MOTIVATES: rather than explicitly inverting M (which
would require an arbitrary regularization choice in that 19-dimensional
null space), M is folded DIRECTLY into Phase 1's own SDP as part of the
forward (not inverse) model:

    rho_S = argmin_rho  sum_i w_i ( (M @ x(rho))_i - m_noisy_i )^2
    s.t.   rho Hermitian, rho >> 0, Tr(rho) = 1,   x(rho)_j = Tr(rho P_S[j])

Phase 1 is exactly the M = Identity special case of this. The
physicality constraint itself -- a Hermitian, PSD, trace-1 K x K matrix
has only K^2-1 = 35 real degrees of freedom, and PSD further restricts
that -- supplies the regularization the rank-17 calibration data cannot
supply on its own, with no separately-chosen ridge parameter needed.

ALSO REPORTED (per Phase 3's explicit ask): what a NAIVE, UNCONSTRAINED
inversion would look like on the SAME learned M, applied to real target
data -- to show concretely whether/how much the physicality constraint
stabilises an otherwise unstable inverse problem.

A REAL, DISCLOSED CAVEAT: the calibration circuits and the real target
circuits were measured in SEPARATE real submissions (iteration 9's own
calibration phase vs targets phase) -- the same kind of batch-to-batch
concern iteration 17 flagged for a different comparison in this project.
Not resolved here; disclosed as a limitation of what these two existing
datasets can support.

Run:
    python vqe/phase3_constrained_channel_inversion.py
"""
import os
import sys
import json
import time
import numpy as np
import cvxpy as cp

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL, slot_names,
)
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed, expectation_from_counts
from phys_constrained_reconstruction import load_checkpoint, build_group_index, build_P_S, K, SHOTS, N_SEEDS

CAL_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "calibrate.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "phase3_channel_inversion_results.json")


def load_calibration_matrix(model, labels):
    with open(CAL_PATH) as f:
        cal = json.load(f)
    tbs = cal["per_model"][model]["training_by_seed"]
    by_slot = {}
    for seed, entries in tbs.items():
        for e in entries:
            by_slot.setdefault((seed, e["slot"]), {})[e["label"]] = (e["exact"], e["noisy"])
    ideal_cols, noisy_cols = [], []
    for key, vals in by_slot.items():
        if all(l in vals for l in labels):
            ideal_cols.append([vals[l][0] for l in labels])
            noisy_cols.append([vals[l][1] for l in labels])
    X = np.array(ideal_cols).T  # (n_labels, n_circuits)
    Y = np.array(noisy_cols).T
    return X, Y


def fit_M(X, Y):
    """Minimum-norm least-squares M such that M @ X ~= Y (np.linalg.pinv
    handles the rank deficiency via SVD truncation at the numerical noise
    floor -- the standard, well-defined choice when the fitting problem
    is underdetermined, not an arbitrary one)."""
    M = Y @ np.linalg.pinv(X)
    svals_X = np.linalg.svd(X, compute_uv=False)
    resid = np.linalg.norm(M @ X - Y) / np.linalg.norm(Y)
    return M, svals_X, resid


def reconstruct_rho_slot_channel(P_S, M, labels, m_noisy_dict, w_dict, K):
    """Phase 1's SDP generalized to a learned channel M: predict
    (M @ x(rho))_i instead of x(rho)_i directly. M=Identity recovers
    Phase 1 exactly.

    SPEED FIX: the first version built the objective as a Python-level
    sum of 36 separately-squared cp.Expression terms
    (`cp.sum([w_i*cp.square(...) for i in range(36)])`), which cvxpy's
    own runtime flagged directly ("Objective contains too many
    subexpressions. Consider vectorizing") -- confirmed by direct timing:
    2.3-4.0s per solve, vs Phase 1/2's ~0.15-0.3s for a structurally
    similar-sized problem. Fixed by using a single vectorized
    cp.sum_squares atom over a weighted residual vector instead of 36
    separate square atoms summed -- the SAME mathematical objective
    (sqrt(w_i) scaling before squaring is algebraically identical to
    w_i*square(.)), just compiled as one vectorized operation."""
    rho = cp.Variable((K, K), hermitian=True)
    x_terms = [cp.real(cp.trace(rho @ P_S[l])) for l in labels]  # x(rho), one entry per label
    x_vec = cp.hstack(x_terms)
    pred = M @ x_vec  # (n_labels,) predicted noisy vector
    m_vec = np.array([m_noisy_dict[l] for l in labels])
    sqrt_w = np.array([np.sqrt(w_dict[l]) for l in labels])
    resid = cp.multiply(sqrt_w, pred - m_vec)
    obj = cp.sum_squares(resid)
    prob = cp.Problem(cp.Minimize(obj), [rho >> 0, cp.trace(rho) == 1])
    prob.solve(solver=cp.SCS)
    if rho.value is None:
        raise RuntimeError(f"SDP solve failed, status={prob.status}")
    return rho.value


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def run_scheme(ck_model, p, P_S, M, non_id_labels, groups, group_idx, names, K, n_seeds, tag):
    """8-seed bootstrap sweep for one noise model's counts, comparing raw
    vs Phase1(M=I) vs Phase3(learned M). Used for BOTH the real aria-1/
    forte-1 evaluation AND the mandatory ideal-data sanity check below --
    same code path, so the check is not a weaker/different test than the
    real evaluation."""
    counts_by_slot = ck_model
    errs_raw, errs_p1, errs_p3 = [], [], []
    for seed in range(n_seeds):
        rng = np.random.default_rng(stable_seed(tag, "seed", seed))
        raw_baseline = {name: {} for name in names}
        raw_p1 = {name: {} for name in names}
        raw_p3 = {name: {} for name in names}
        for name in names:
            resampled = [bootstrap_counts(counts_by_slot[name][gi], SHOTS, rng) for gi in range(len(groups))]
            m_dict, w_dict = {}, {}
            for l in non_id_labels:
                counts = resampled[group_idx[l]]
                m = expectation_from_counts(counts, l)
                m_dict[l] = m
                total = sum(counts.values())
                w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(total, 1)
                raw_baseline[name][l] = m
            rho_p1 = reconstruct_rho_slot_channel(P_S, np.eye(len(non_id_labels)), non_id_labels, m_dict, w_dict, K)
            rho_p3 = reconstruct_rho_slot_channel(P_S, M, non_id_labels, m_dict, w_dict, K)
            for l in non_id_labels:
                raw_p1[name][l] = float(np.real(np.trace(rho_p1 @ P_S[l])))
                raw_p3[name][l] = float(np.real(np.trace(rho_p3 @ P_S[l])))
        _, err_raw = energy_and_err(p, raw_baseline, K)
        _, err_p1 = energy_and_err(p, raw_p1, K)
        _, err_p3 = energy_and_err(p, raw_p3, K)
        errs_raw.append(err_raw["err_vs_exact_kcal"])
        errs_p1.append(err_p1["err_vs_exact_kcal"])
        errs_p3.append(err_p3["err_vs_exact_kcal"])
    return {
        "raw_mean": float(np.mean(errs_raw)), "raw_std": float(np.std(errs_raw)),
        "phase1_mean": float(np.mean(errs_p1)), "phase1_std": float(np.std(errs_p1)),
        "phase3_mean": float(np.mean(errs_p3)), "phase3_std": float(np.std(errs_p3)),
    }


def main():
    print("\n" + "=" * 96)
    print("  phase3_constrained_channel_inversion.py -- continued at user's explicit direction")
    print("=" * 96)

    ck = load_checkpoint()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=ck["d"], K=K)
    U = np.asarray(p["u_vecs"]).T
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    non_id_labels = sorted([l for l in alpha_labels if l != identity_label])
    groups = ck["groups"]
    group_idx = build_group_index(groups)
    P_S = build_P_S(alpha_labels, U)
    names = ck["target_names"]

    results = {}
    for model in ["aria-1", "forte-1"]:
        print(f"\n  -- {model} --")
        X, Y = load_calibration_matrix(model, non_id_labels)
        M, svals_X, resid = fit_M(X, Y)
        eff_rank = int(np.sum(svals_X > 1e-6))
        print(f"  calibration: {X.shape[1]} circuits, {X.shape[0]} labels, "
              f"X effective rank={eff_rank}/{X.shape[0]} (rank-deficient: {eff_rank < X.shape[0]})")
        print(f"  M fit residual (||MX-Y||/||Y||): {resid:.4f}")
        svals_M = np.linalg.svd(M, compute_uv=False)
        print(f"  M singular values: min={svals_M.min():.4f} max={svals_M.max():.4f} "
              f"condition number={svals_M.max()/max(svals_M.min(),1e-12):.2e}")
        off_diag_frac = (np.sum(np.abs(M)) - np.sum(np.abs(np.diag(M)))) / np.sum(np.abs(M))
        print(f"  fraction of |M|'s mass off-diagonal (cross-Pauli mixing, 0=pure per-label scaling like CDR): "
              f"{off_diag_frac:.3f}")

        # naive UNCONSTRAINED single-target inversion, to show what physicality buys us
        counts_by_slot = ck["counts"][model]
        rng0 = np.random.default_rng(0)
        resampled0 = [bootstrap_counts(counts_by_slot["u_0"][gi], SHOTS, rng0) for gi in range(len(groups))]
        m0 = np.array([expectation_from_counts(resampled0[group_idx[l]], l) for l in non_id_labels])
        x_hat_naive = np.linalg.pinv(M) @ m0
        n_out_of_range = int(np.sum(np.abs(x_hat_naive) > 1.0))
        print(f"  naive UNCONSTRAINED inversion on a real target (slot u_0): "
              f"{n_out_of_range}/{len(non_id_labels)} entries land outside the physical [-1,1] range "
              f"(max |value|={np.max(np.abs(x_hat_naive)):.2f}) -- unconstrained inversion is unstable")

        t0 = time.time()
        real_result = run_scheme(counts_by_slot, p, P_S, M, non_id_labels, groups, group_idx, names, K,
                                  N_SEEDS, tag=f"phase3_{model}")
        wall_clock = time.time() - t0

        # MANDATORY SANITY GATE, run before trusting the real-noise-model number at all:
        # apply this SAME learned M to the near-noiseless IDEAL data. A genuine noise-
        # correcting channel should leave already-clean data close to correct. If Phase3's
        # error on IDEAL data is much WORSE than raw/Phase1 there, that proves the fit is
        # exploiting an underdetermined/overfitting freedom rather than genuinely inverting
        # noise -- and the real-model number above must be reported as DISQUALIFIED, not an
        # achievement, regardless of how good it looks.
        ideal_result = run_scheme(ck["counts"]["ideal"], p, P_S, M, non_id_labels, groups, group_idx, names, K,
                                   3, tag=f"phase3_idealcheck_{model}")
        distortion_ratio = ideal_result["phase3_mean"] / max(ideal_result["phase1_mean"], 1e-6)
        disqualified = ideal_result["phase3_mean"] > 3 * max(ideal_result["raw_mean"], ideal_result["phase1_mean"])

        results[model] = {
            "n_calibration_circuits": X.shape[1], "X_effective_rank": eff_rank, "X_n_labels": X.shape[0],
            "M_fit_residual": float(resid), "M_condition_number": float(svals_M.max()/max(svals_M.min(),1e-12)),
            "M_off_diag_mass_fraction": float(off_diag_frac),
            "naive_inversion_n_out_of_range": n_out_of_range, "naive_inversion_max_abs": float(np.max(np.abs(x_hat_naive))),
            "real_model_result": real_result,
            "ideal_sanity_check": ideal_result,
            "ideal_distortion_ratio_phase3_over_phase1": float(distortion_ratio),
            "DISQUALIFIED": bool(disqualified),
            "wall_clock_s": wall_clock,
        }
        r = results[model]
        print(f"  [{wall_clock:.0f}s] raw={r['real_model_result']['raw_mean']:.2f}+/-{r['real_model_result']['raw_std']:.2f}  "
              f"Phase1(M=I)={r['real_model_result']['phase1_mean']:.2f}+/-{r['real_model_result']['phase1_std']:.2f}  "
              f"Phase3(learned M)={r['real_model_result']['phase3_mean']:.2f}+/-{r['real_model_result']['phase3_std']:.2f} kcal/mol")
        print(f"  MANDATORY IDEAL-DATA SANITY CHECK: raw={ideal_result['raw_mean']:.2f}  "
              f"Phase1={ideal_result['phase1_mean']:.2f}  Phase3={ideal_result['phase3_mean']:.2f} kcal/mol on near-noiseless data")
        print(f"  {'*** DISQUALIFIED ***' if disqualified else 'PASSES sanity check'} -- "
              f"{'Phase3 severely distorts already-clean data (>3x worse than raw/Phase1 there); the real-model number above is NOT a genuine correction, do not report it as an achievement' if disqualified else 'Phase3 does not distort clean data; real-model number may be trustworthy'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
