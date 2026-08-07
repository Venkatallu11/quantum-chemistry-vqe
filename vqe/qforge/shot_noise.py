#!/usr/bin/env python3
"""
qforge.shot_noise — the shot-noise model and shots-vs-accuracy harness,
extracted from shot_noise_study.py. Physics unchanged: Binomial sample-
mean estimator for a +-1-eigenvalue Pauli observable (the exact
distribution real shot-based hardware produces, not an approximation),
plus PEC's gamma_total^2-inflated-variance estimator (exact-by-
construction result of the quasi-probability decomposition, not a
Monte-Carlo simulation of it).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

HARTREE_TO_KCAL_MOL = 627.5094740631


def shot_sample(exact_value, n_shots, rng):
    """Exact Binomial sample-mean estimator of a +-1-eigenvalue observable
    with true expectation `exact_value`, from n_shots projective
    measurements."""
    p_plus = float(np.clip((1 + exact_value) / 2, 0.0, 1.0))
    n_plus = rng.binomial(n_shots, p_plus)
    return 2 * n_plus / n_shots - 1


def pec_shot_sample(exact_value, gamma_total, n_shots, rng):
    """PEC's quasi-probability Monte Carlo estimator: unbiased, variance
    inflated by gamma_total^2 relative to a plain shot-noise measurement
    of the same observable -- the standard, exact-by-construction PEC
    sampling-overhead result (van den Berg et al.; Temme, Bravyi,
    Gambetta)."""
    var = (gamma_total ** 2) * max(1 - exact_value ** 2, 0.0) / n_shots
    return exact_value + rng.normal(0.0, np.sqrt(var))


def verify_convergence(exact_energy_fn, shot_energy_fn, n_shots_list, seed=0):
    """MANDATORY verification before trusting any shot-noise sweep: does
    the shot-noisy estimator converge to the exact value as shots grow?
    exact_energy_fn() -> float, shot_energy_fn(n_shots, rng) -> float."""
    rng = np.random.default_rng(seed)
    E_exact = exact_energy_fn()
    rows = []
    for n_shots in n_shots_list:
        E = shot_energy_fn(n_shots, rng)
        diff_kcal = abs(E - E_exact) * HARTREE_TO_KCAL_MOL
        rows.append({"n_shots": n_shots, "E": E, "diff_from_exact_kcal": diff_kcal})
    converges = bool(rows[-1]["diff_from_exact_kcal"] < rows[0]["diff_from_exact_kcal"])
    return rows, converges, E_exact


def verify_scaling_exponent(shot_energy_fn, n_shots_pair, n_trials, seed_base=5000, tol=0.3):
    """Confirm empirical std scales as ~1/sqrt(N) between two shot
    counts -- the scaling law is more fundamental (and more robust to
    verify) than matching any specific closed-form variance prediction
    exactly, which this project's own bilinear (beta=S.alpha.S) energy
    estimator does NOT do (see shot_noise_study.py's honest mismatch
    report -- reported as measured, not gated on)."""
    stds = []
    for n_shots in n_shots_pair:
        energies = [shot_energy_fn(n_shots, np.random.default_rng(seed_base + trial))
                    for trial in range(n_trials)]
        stds.append(float(np.std(energies)))
    ratio_shots = n_shots_pair[1] / n_shots_pair[0]
    ratio_std = stds[0] / stds[1]
    expected_ratio_std = np.sqrt(ratio_shots)
    scaling_ok = abs(ratio_std - expected_ratio_std) / expected_ratio_std < tol
    return {"n_shots_pair": list(n_shots_pair), "stds": stds, "ratio_std": ratio_std,
            "expected_ratio_std": expected_ratio_std, "scaling_ok": bool(scaling_ok)}


def summarize(vals, chem_acc=1.0, target=0.30):
    vals = list(vals)
    return {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
            "min": float(np.min(vals)), "max": float(np.max(vals)),
            "reached_chem_acc": f"{sum(1 for v in vals if v < chem_acc)}/{len(vals)}",
            "reached_target": f"{sum(1 for v in vals if v < target)}/{len(vals)}"}
