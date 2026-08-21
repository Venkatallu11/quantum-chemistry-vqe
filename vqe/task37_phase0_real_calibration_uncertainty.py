#!/usr/bin/env python3
"""
task37_phase0_real_calibration_uncertainty.py -- iteration 37, Phase 0.
Task 36B's robustness envelope draws p_ZZ ~ Uniform(0.010,0.020) and
p_GPi ~ Uniform(0.00005,0.0006), independently -- but this project ALREADY
has real, measured calibration uncertainty for both, sitting unused:
  Task 31A (real ZZ calibration, 9/11 positions): mean=0.014593,
    std=0.000124 across positions -- ~23x TIGHTER than the assumed
    uniform range's implied std (~0.00289).
  Task 30B (real GPi calibration, 4 phase bins): mean=0.000119,
    std=0.000011 -- ~14x TIGHTER than the assumed range's implied std
    (~0.000159).
  GPi2 (Task 31A recalibration attempt): genuinely unresolved -- the
    file's own numbers are flagged "physical": false, with nonsensical
    values (negative probabilities, std up to 0.19). No real posterior
    exists; left UNCHANGED (still Uniform(0,0.21)), not silently
    tightened.

This tests directly: how much of Task 36B's Q95=18.29 (standard's
Q95=70.74) was an artifact of using unrealistically wide, INDEPENDENT
assumed ranges for parameters we actually know much more precisely than
that -- vs. how much is genuinely irreducible given GPi2's real,
unresolved uncertainty?

Everything else UNCHANGED from task31h_robustness_envelope.py /
task36b_joint_frame_robustness_envelope.py: angle biases, readout error,
calib_ratio_zz/1q (these reflect assumed-vs-true MISCALIBRATION risk, a
different, real thing calibration measurement alone doesn't resolve, not
touched here), the joint/standard estimator pipelines themselves.

Run:
    PYTHONHASHSEED=0 python vqe/task37_phase0_real_calibration_uncertainty.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from task31h_robustness_envelope import noisy_pec_dm, apply_readout_error, GATE_NAME
from task36_joint_schmidt_frame import procrustes_init, fit_joint_frame, build_full_from_frame
from qiskit.quantum_info import Pauli

K = 6
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37_phase0_real_calibration_uncertainty_results.json")
BUDGET_S = 9000

# -- real measured values, computed directly from task31a/task30b's own result files --
ZZ_REAL_MEAN, ZZ_REAL_STD = 0.014593, 0.000124
GPI_REAL_MEAN, GPI_REAL_STD = 0.000119, 0.000011


def sample_noise_model_realistic(rng):
    """SAME as task31h.sample_noise_model, EXCEPT p_zz_true and
    p_gpi_true are drawn from their REAL measured Gaussian posteriors
    instead of the original wide, independent uniform assumptions.
    p_gpi2_true, angle biases, readout, calib ratios: UNCHANGED."""
    return {
        "p_zz_true": max(1e-6, float(rng.normal(ZZ_REAL_MEAN, ZZ_REAL_STD))),
        "p_gpi_true": max(1e-6, float(rng.normal(GPI_REAL_MEAN, GPI_REAL_STD))),
        "p_gpi2_true": float(rng.uniform(0.0, 0.21)),  # UNCHANGED -- genuinely unresolved, not tightened
        "angle_bias_zz": float(rng.normal(0, 0.002)),
        "angle_bias_gpi": float(rng.normal(0, 0.002)),
        "angle_bias_gpi2": float(rng.normal(0, 0.002)),
        "readout_err": float(rng.uniform(0.0, 0.02)),
        "calib_ratio_zz": float(rng.uniform(0.85, 1.15)),
        "calib_ratio_1q": float(rng.uniform(0.85, 1.15)),
    }


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def evaluate_both(p, fixed_solutions, kept, diag, P_S, non_id_labels, model, rng):
    m_dict_by_slot = {}
    for name in kept:
        dm = noisy_pec_dm(fixed_solutions[name]["angles"], GATE_NAME, model, rng)
        dm = apply_readout_error(dm, dm.num_qubits, model["readout_err"], rng)
        m_dict = {}
        for l in non_id_labels:
            Pmat = np.asarray(Pauli(l).to_matrix())
            m_dict[l] = float(np.real(np.trace(Pmat @ dm.data)))
        m_dict_by_slot[name] = m_dict

    a_by_name = {}
    for name in kept:
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, m_dict_by_slot[name], {l: 1.0 for l in m_dict_by_slot[name]}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full_std = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err_std = energy_and_err(p, full_std, K)

    U0 = procrustes_init({n: a_by_name[n] for n in diag}, K)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    U_hat, cost, chi2dof = fit_joint_frame(U0, P_S, K, kept, non_id_labels, m_dict_by_slot, weight_unit, rng)
    full_joint = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    _, err_joint = energy_and_err(p, full_joint, K)

    return err_std, err_joint, chi2dof


def main():
    print("\n" + "=" * 96)
    print("  task37_phase0_real_calibration_uncertainty.py -- real vs assumed calibration width")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- Task 36's known nondeterminism bug may reappear.")

    print(f"\n  ZZ:  assumed Uniform(0.010,0.020) implied_std={0.01/np.sqrt(12):.6f}  "
          f"REAL mean={ZZ_REAL_MEAN:.6f} std={ZZ_REAL_STD:.6f}  "
          f"(assumed/real width ratio = {(0.01/np.sqrt(12))/ZZ_REAL_STD:.1f}x)")
    print(f"  GPi: assumed Uniform(0.00005,0.0006) implied_std={0.00055/np.sqrt(12):.6f}  "
          f"REAL mean={GPI_REAL_MEAN:.6f} std={GPI_REAL_STD:.6f}  "
          f"(assumed/real width ratio = {(0.00055/np.sqrt(12))/GPI_REAL_STD:.1f}x)")
    print(f"  GPi2: UNCHANGED, Uniform(0,0.21) -- no real posterior exists (recalibration flagged "
          f"'physical': false, std up to 0.19)")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    rng0 = np.random.default_rng(0)
    model0 = sample_noise_model_realistic(rng0)
    t0 = time.time()
    e0_std, e0_joint, chi2dof0 = evaluate_both(p, fixed_solutions, kept, diag, P_S, non_id_labels, model0, rng0)
    t_per_eval = time.time() - t0
    print(f"\n  ONE evaluation: std={e0_std:.4f}  joint={e0_joint:.4f}  wall-clock={t_per_eval:.2f}s")
    N = max(20, min(600, int(BUDGET_S / max(t_per_eval, 0.01))))
    print(f"  setting N={N} based on measured per-eval cost (budget ~{BUDGET_S}s)")

    PARTIAL_PATH = RESULTS_PATH + ".partial.json"
    errs_std, errs_joint, chi2dofs = [], [], []
    start_i = 0
    if os.path.exists(PARTIAL_PATH):
        with open(PARTIAL_PATH) as f:
            partial = json.load(f)
        errs_std, errs_joint, chi2dofs = partial["errs_std"], partial["errs_joint"], partial["chi2dofs"]
        start_i = len(errs_std)
        print(f"  resuming from partial checkpoint: {start_i}/{N} trials already done")
    rng = np.random.default_rng(37 + start_i * 1009)

    t_start = time.time()
    for i in range(start_i, N):
        model = sample_noise_model_realistic(rng)
        e_std, e_joint, chi2dof = evaluate_both(p, fixed_solutions, kept, diag, P_S, non_id_labels, model, rng)
        errs_std.append(e_std)
        errs_joint.append(e_joint)
        chi2dofs.append(chi2dof)
        if (i + 1) % max(1, N // 20) == 0 or (i + 1) == N:
            print(f"    {i+1}/{N} done, {time.time()-t_start:.1f}s elapsed, "
                  f"last: std={e_std:.3f} joint={e_joint:.3f}", flush=True)
            with open(PARTIAL_PATH, "w") as f:
                json.dump({"errs_std": errs_std, "errs_joint": errs_joint, "chi2dofs": chi2dofs}, f)
    if os.path.exists(PARTIAL_PATH):
        os.remove(PARTIAL_PATH)

    errs_std = np.array(errs_std)
    errs_joint = np.array(errs_joint)
    chi2dofs = np.array(chi2dofs)

    q_std = {q_: float(np.percentile(errs_std, q_)) for q_ in [50, 90, 95, 99]}
    q_joint = {q_: float(np.percentile(errs_joint, q_)) for q_ in [50, 90, 95, 99]}

    print(f"\n  -- QUANTILES under {N} noise models, REAL ZZ/GPi uncertainty, GPi2 unchanged --")
    print(f"    {'':<8}{'STANDARD':>12}{'JOINT':>12}")
    for q_ in [50, 90, 95, 99]:
        print(f"    Q{q_:<7}{q_std[q_]:>12.4f}{q_joint[q_]:>12.4f}")

    n_outliers_std = int((errs_std > 20).sum())
    n_outliers_joint = int((errs_joint > 20).sum())
    n_joint_better = int(((errs_std - errs_joint) > 0).sum())
    print(f"\n  outlier rate (>20 kcal/mol): standard {n_outliers_std}/{N}, joint {n_outliers_joint}/{N}")
    print(f"  paired: joint wins {n_joint_better}/{N} trials")

    print(f"\n  -- COMPARISON vs Task 36B's original (independent wide-uniform) result --")
    print(f"    Task 36B (N=267): standard Q95=70.74  joint Q95=18.29")
    print(f"    Task 37 Phase 0 (N={N}): standard Q95={q_std[95]:.2f}  joint Q95={q_joint[95]:.2f}")
    if q_joint[95] < 18.29:
        pct = 100 * (18.29 - q_joint[95]) / 18.29
        print(f"    -> Using REAL ZZ/GPi uncertainty recovers {pct:.1f}% of joint's Q95, even leaving "
              f"GPi2 fully unresolved. The remaining gap is attributable to GPi2 and/or angle-bias/"
              f"readout/calib-ratio uncertainty, not to unrealistically wide ZZ/GPi assumptions.")
    else:
        print(f"    -> Real ZZ/GPi uncertainty did NOT meaningfully reduce Q95 -- the tail is dominated "
              f"by GPi2 and/or other parameters, not by ZZ/GPi being modeled too pessimistically.")

    print(f"\n  TARGET: Q95 < 0.25 -> joint {'MET' if q_joint[95] < 0.25 else 'NOT MET'}")
    print(f"  LOOSENED HARDWARE BAR: Q95 < 0.5 -> joint {'MET' if q_joint[95] < 0.5 else 'NOT MET'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N": N, "t_per_eval_s": t_per_eval,
            "zz_real_mean": ZZ_REAL_MEAN, "zz_real_std": ZZ_REAL_STD,
            "gpi_real_mean": GPI_REAL_MEAN, "gpi_real_std": GPI_REAL_STD,
            "errs_std": errs_std.tolist(), "errs_joint": errs_joint.tolist(), "chi2dofs": chi2dofs.tolist(),
            "quantiles_std": q_std, "quantiles_joint": q_joint,
            "n_outliers_std": n_outliers_std, "n_outliers_joint": n_outliers_joint,
            "n_joint_better": n_joint_better,
            "task36b_reference_q95_std": 70.74, "task36b_reference_q95_joint": 18.29,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
