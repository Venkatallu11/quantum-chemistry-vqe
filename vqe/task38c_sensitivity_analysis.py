#!/usr/bin/env python3
"""
task38c_sensitivity_analysis.py -- iteration 38, Task C (numbered to
match the user's own "38C" proposal). The gating calculation before any
robust-PEC optimization is attempted: compute g = dE/dtheta (H4 energy
sensitivity to the 5-parameter noise model) at the nominal calibration
point, and V_cal = g^T Sigma_theta g using the BEST REAL uncertainty
estimate for each parameter collected across Tasks 37B-E (not the old
independent wide priors) -- then compare V_cal against Task 36's own
real N_BOOT=80 paired-bootstrap variance (std_joint=0.4858 kcal/mol,
V=0.236 kcal/mol^2), which captures shot-noise + PEC-sampling variance
+ manifold-fit variance at FIXED (nominal) calibration. This directly
answers the decision-tree question the proposal raised: does calibration
uncertainty dominate the error budget (-> pursue robust PEC / randomized
compiling), or does something else dominate (-> a different branch is
the right next move)?

PIPELINE USED: the STANDARD (non-joint-frame) PEC-corrected energy --
raw real hardware data (`task28b_optimized_raw.json`) corrected via
Task 37C's own verified `analytic_A_and_B_5param` + readout un-mixing,
combined directly via `qforge`'s standard machinery (matching Task 30B/
31C's original pipeline, no Schmidt-frame layered on top) -- this is the
right choice for a PEC-sensitivity analysis specifically: the question
"how sensitive is the PEC step itself to calibration uncertainty" is
independent of whether a joint frame is applied afterward, and using the
simpler pipeline avoids re-deriving Task 36's own separate frame-fitting
machinery as a nuisance parameter here.

SIGMA_THETA, disclosed source per entry (a simplified, DIAGONAL version
of the "correlated constrained posterior" the proposal called for --
correlations between parameters are a real refinement not attempted
here):
  p_zz:       0.000124  (Task 31A real calibration std)
  p_gpi2:     0.00081   (Task 37D's Schur-complement DATA-driven std --
              tighter and more real than the old weak prior, since Task
              37D found the H4 data itself constrains this parameter)
  delta_zz:   0.000249  (Task 37E's real repeated-submission SEM, more
              conservative of the two backends' SEM)
  delta_gpi2: 0.00161   (same, Task 37E's real repeated-submission SEM)
  p_readout:  0.01      (Task 37B's original weak prior -- still no real
              measurement exists anywhere for this parameter)

Run:
    PYTHONHASHSEED=0 python vqe/task38c_sensitivity_analysis.py
"""
import os
import sys
import numpy as np
from scipy.optimize._numdiff import approx_derivative

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import build_full
from task37c_extended_forward_model import analytic_A_and_B_5param, readout_attenuation
from task37c_joint_calibration_h4_fit import load_real_blended
from task37b_h4_noise_model import PARAM_NAMES, GPI_REAL_MEAN

K = 6
GATE_NAME = "zz"
THETA0 = np.array([0.014593, 0.0, 0.0, 0.0, 0.0])  # nominal: real ZZ calib, everything else the
# historical "no correction applied" default this project's standard (non-Task-37) pipeline has
# always used
SIGMA_DIAG = np.array([0.000124, 0.00081, 0.000249, 0.00161, 0.01])
BOOTSTRAP_STD_JOINT_KCAL = 0.4857515793056967  # Task 36's real N_BOOT=80 result, reused not recomputed
ABS_STEP = np.array([2e-4, 2e-4, 2e-3, 2e-3, 2e-4])  # per-param, comfortably above the cache-rounding
# resolution and small relative to each parameter's own scale -- learned directly from Task 37D's
# rel_step bug: NEVER use a bare rel_step here, always an explicit, sane abs_step per parameter
_GLOBAL_AB_CACHE = {}


def energy_of_theta(theta, kept, non_id_labels, fixed_solutions, diag, real_blended, p):
    p_zz, p_gpi2, delta_zz, delta_gpi2, p_readout = theta
    raw_kept = {}
    for name in kept:
        labels_here = list(real_blended[name].keys())
        cache_key = (name, round(p_zz, 6), round(p_gpi2, 6), round(delta_zz, 6), round(delta_gpi2, 6))
        if cache_key not in _GLOBAL_AB_CACHE:
            A, B = analytic_A_and_B_5param(fixed_solutions[name]["angles"], GATE_NAME, p_zz, GPI_REAL_MEAN,
                                            p_gpi2, delta_zz, delta_gpi2, labels_here)
            _GLOBAL_AB_CACHE[cache_key] = (A, B)
        A, B = _GLOBAL_AB_CACHE[cache_key]
        raw_kept[name] = {}
        for l, m_raw in real_blended[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            ro = readout_attenuation(l, p_readout)
            ro = ro if abs(ro) > 1e-6 else 1.0
            raw_kept[name][l] = max(-1.0, min(1.0, (m_raw / ro) * ratio))
    full = build_full(raw_kept, diag, K, non_id_labels)
    alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]  # differentiating this == differentiating E (exact is a fixed offset)


def main():
    print("\n" + "=" * 96)
    print("  task38c_sensitivity_analysis.py -- dE/dtheta and calibration-uncertainty variance, standard PEC pipeline")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    real_blended, _ = load_real_blended(kept)
    print(f"\n  loaded real raw data ({sum(len(real_blended[n]) for n in kept)} (slot,label) residuals, {len(kept)} slots)")

    E0 = energy_of_theta(THETA0, kept, non_id_labels, fixed_solutions, diag, real_blended, p)
    print(f"\n  E(theta_0) err vs exact = {E0:.4f} kcal/mol  (nominal calibration point: {dict(zip(PARAM_NAMES, THETA0))})")

    def f(theta):
        return energy_of_theta(theta, kept, non_id_labels, fixed_solutions, diag, real_blended, p)

    print("\n  -- computing g = dE/dtheta via finite differences (explicit per-param abs_step, not rel_step) --")
    # BUG CAUGHT: a scalar-valued fun makes approx_derivative return shape (5,) directly (not (1,5)) --
    # an earlier version wrapped f's return in np.array([...]) and indexed J[0], which silently pulled
    # out a single SCALAR (the first param's own derivative) instead of the whole 5-vector, crashing
    # on the very next line's g[i] indexing. Fixed by returning a plain scalar and using J directly.
    g = approx_derivative(f, THETA0, method="3-point", abs_step=ABS_STEP)  # shape (5,)

    print(f"\n  {'param':<12}{'g_i (kcal/mol per unit)':>26}{'sigma_i':>14}{'g_i*sigma_i (kcal/mol)':>26}{'source':>10}")
    contributions = g * SIGMA_DIAG
    for i, name in enumerate(PARAM_NAMES):
        print(f"  {name:<12}{g[i]:>26.4f}{SIGMA_DIAG[i]:>14.3e}{contributions[i]:>26.4f}")

    V_cal = float(np.sum((g * SIGMA_DIAG) ** 2))  # diagonal Sigma_theta: g^T Sigma g = sum (g_i sigma_i)^2
    std_cal = float(np.sqrt(V_cal))
    print(f"\n  V_cal = g^T Sigma_theta g = {V_cal:.4f} kcal/mol^2   (implied calibration-uncertainty std = {std_cal:.4f} kcal/mol)")

    print(f"\n  -- per-parameter share of V_cal --")
    for i, name in enumerate(PARAM_NAMES):
        share = 100 * contributions[i] ** 2 / V_cal if V_cal > 1e-12 else 0.0
        print(f"    {name:<12} {share:5.1f}%")

    V_other = BOOTSTRAP_STD_JOINT_KCAL ** 2
    print(f"\n  -- DECISION-TREE COMPARISON (per the proposal's own point 17) --")
    print(f"    V_cal (this calculation)                         = {V_cal:.4f} kcal/mol^2  (std={std_cal:.4f})")
    print(f"    V_other (Task 36 real bootstrap: shot+PEC-MC+manifold, FIXED calibration) = {V_other:.4f} kcal/mol^2  (std={BOOTSTRAP_STD_JOINT_KCAL:.4f})")
    frac_cal = 100 * V_cal / (V_cal + V_other)
    print(f"    calibration share of (V_cal + V_other): {frac_cal:.1f}%")
    if frac_cal > 70:
        verdict = "V_cal DOMINATES -- pursue robust PEC / randomized compiling (proposal's 38D/38G branch)"
    elif frac_cal < 20:
        verdict = "V_other DOMINATES -- return to QPD stratification / control-variate machinery, NOT robust PEC"
    else:
        verdict = "MIXED -- neither dominates cleanly; both branches plausibly worth pursuing, no clean single winner"
    print(f"    VERDICT: {verdict}")


if __name__ == "__main__":
    main()
