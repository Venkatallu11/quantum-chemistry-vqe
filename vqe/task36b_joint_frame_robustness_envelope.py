#!/usr/bin/env python3
"""
task36b_joint_frame_robustness_envelope.py -- iteration 36, Task B. Does
the joint Schmidt-frame estimator's real 88.8% MSE win (Task 36, 80 real
bootstrap trials on ONE fixed hardware collection) survive under NOISE-
MODEL uncertainty, not just resampling noise from one dataset? This is
the same open question flagged explicitly in Task 36's own writeup, and
the same methodology (`task31h_robustness_envelope.py`) that has been
this project's real reliability gate since iteration 31 -- the test that
would have caught the false-positive excitement around 0.115 kcal/mol
(iteration 31) had it been applied there sooner.

METHOD: reuses Task 31h's noise-model sampling and noisy-density-matrix
simulation machinery UNCHANGED (`sample_noise_model`, `noisy_pec_dm`,
`apply_readout_error` -- randomized p_ZZ/p_GPi/p_GPi2, coherent angle
bias per gate TYPE, PEC calibration mismatch, readout error, all over
the SAME intervals Task 31h justified from real measurements). For EACH
randomized noise-model draw, the SAME noisy 21-slot data is fed through
BOTH estimators -- a fair, paired, same-noise-draw comparison:
  STANDARD: `fit_pure_state` per slot independently (exactly Task 31h's
  own `evaluate_one_model`, unchanged).
  JOINT: the same 21 slots' data fit via the shared 15-parameter Schmidt
  frame (Task 36's `fit_joint_frame`/`build_full_from_frame`, unchanged).

Zero shot noise (exact density-matrix expectation values, matching Task
31h's own convention) -- this isolates whether the joint-frame estimator
is robust to CALIBRATION/NOISE-MODEL uncertainty specifically, a
different question from Task 36's shot+PEC-MC resampling test.

IMPORTANT: run with PYTHONHASHSEED=0 -- Task 36 found and fixed a real
cross-process nondeterminism bug in the joint-frame optimizer (hash-order-
dependent residual ordering tipping the nonconvex solver into different
local optima). That fix is required here too, not optional.

Run:
    PYTHONHASHSEED=0 python vqe/task36b_joint_frame_robustness_envelope.py
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
from task31h_robustness_envelope import sample_noise_model, noisy_pec_dm, apply_readout_error, GATE_NAME
from task36_joint_schmidt_frame import procrustes_init, fit_joint_frame, build_full_from_frame
from qiskit.quantum_info import Pauli

K = 6
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task36b_joint_frame_robustness_envelope_results.json")
BUDGET_S = 2400  # ~40 min, disclosed


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def evaluate_both(p, fixed_solutions, kept, diag, P_S, non_id_labels, model, rng):
    """Generate ONE noisy realization, evaluate STANDARD and JOINT on the SAME data."""
    m_dict_by_slot = {}
    for name in kept:
        dm = noisy_pec_dm(fixed_solutions[name]["angles"], GATE_NAME, model, rng)
        dm = apply_readout_error(dm, dm.num_qubits, model["readout_err"], rng)
        m_dict = {}
        for l in non_id_labels:
            Pmat = np.asarray(Pauli(l).to_matrix())
            m_dict[l] = float(np.real(np.trace(Pmat @ dm.data)))
        m_dict_by_slot[name] = m_dict

    # -- STANDARD: 21 independent fits, exactly Task 31h's own evaluate_one_model --
    a_by_name = {}
    for name in kept:
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, m_dict_by_slot[name], {l: 1.0 for l in m_dict_by_slot[name]}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full_std = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err_std = energy_and_err(p, full_std, K)

    # -- JOINT: same 21 slots' data, shared 15-param frame --
    U0 = procrustes_init({n: a_by_name[n] for n in diag}, K)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    U_hat, cost, chi2dof = fit_joint_frame(U0, P_S, K, kept, non_id_labels, m_dict_by_slot, weight_unit, rng)
    full_joint = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    _, err_joint = energy_and_err(p, full_joint, K)

    return err_std, err_joint, chi2dof


def main():
    print("\n" + "=" * 96)
    print("  task36b_joint_frame_robustness_envelope.py -- joint frame vs standard, under noise-model uncertainty")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- Task 36's known nondeterminism bug may reappear. Proceeding, "
              "results below carry a real risk of the same instability if this warning is not heeded.")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    rng0 = np.random.default_rng(0)
    model0 = sample_noise_model(rng0)
    t0 = time.time()
    e0_std, e0_joint, chi2dof0 = evaluate_both(p, fixed_solutions, kept, diag, P_S, non_id_labels, model0, rng0)
    t_per_eval = time.time() - t0
    print(f"  ONE evaluation: std={e0_std:.4f}  joint={e0_joint:.4f}  chi2dof={chi2dof0:.4f}  "
          f"wall-clock={t_per_eval:.2f}s")
    N = max(20, min(300, int(BUDGET_S / max(t_per_eval, 0.01))))
    print(f"  setting N={N} based on measured per-eval cost (budget ~{BUDGET_S}s)")

    rng = np.random.default_rng(36)
    errs_std, errs_joint, chi2dofs = [], [], []
    t_start = time.time()
    for i in range(N):
        model = sample_noise_model(rng)
        e_std, e_joint, chi2dof = evaluate_both(p, fixed_solutions, kept, diag, P_S, non_id_labels, model, rng)
        errs_std.append(e_std)
        errs_joint.append(e_joint)
        chi2dofs.append(chi2dof)
        if (i + 1) % max(1, N // 10) == 0:
            print(f"    {i+1}/{N} done, {time.time()-t_start:.1f}s elapsed, "
                  f"last: std={e_std:.3f} joint={e_joint:.3f}")

    errs_std = np.array(errs_std)
    errs_joint = np.array(errs_joint)
    chi2dofs = np.array(chi2dofs)

    q_std = {q_: float(np.percentile(errs_std, q_)) for q_ in [50, 90, 95, 99]}
    q_joint = {q_: float(np.percentile(errs_joint, q_)) for q_ in [50, 90, 95, 99]}

    print(f"\n  -- QUANTILES under {N} randomized Forte-like noise models --")
    print(f"    {'':<8}{'STANDARD':>12}{'JOINT':>12}")
    for q_ in [50, 90, 95, 99]:
        print(f"    Q{q_:<7}{q_std[q_]:>12.4f}{q_joint[q_]:>12.4f}")
    print(f"    chi2/dof: median={np.median(chi2dofs):.4f}  max={chi2dofs.max():.4f}")

    n_outliers_std = int((errs_std > 20).sum())
    n_outliers_joint = int((errs_joint > 20).sum())
    diffs = errs_std - errs_joint
    n_joint_better = int((diffs > 0).sum())
    print(f"\n  outlier rate (>20 kcal/mol): standard {n_outliers_std}/{N}, joint {n_outliers_joint}/{N}")
    print(f"  paired: joint wins {n_joint_better}/{N} trials")

    print(f"\n  TARGET: Q95 < 0.25 -> standard {'MET' if q_std[95] < 0.25 else 'NOT MET'}, "
          f"joint {'MET' if q_joint[95] < 0.25 else 'NOT MET'}")
    print(f"  LOOSENED HARDWARE BAR: Q95 < 0.5 -> standard {'MET' if q_std[95] < 0.5 else 'NOT MET'}, "
          f"joint {'MET' if q_joint[95] < 0.5 else 'NOT MET'}")

    print(f"\n  -- HONEST READ --")
    if q_joint[95] < q_std[95] and n_outliers_joint <= n_outliers_std:
        pct = 100 * (q_std[95] - q_joint[95]) / q_std[95]
        print(f"    Joint frame IMPROVES robustness to noise-model uncertainty too: Q95 reduced {pct:.1f}% "
              f"({q_std[95]:.2f} -> {q_joint[95]:.2f}), no worse outlier behavior. The Task 36 real-data win "
              f"generalizes to calibration/noise-model uncertainty, not just shot/PEC-MC resampling noise.")
    else:
        print(f"    Joint frame does NOT show the same clean win under noise-model uncertainty as it did on "
              f"real-data resampling (Task 36). This is real, disclosed information: the shared-frame "
              f"regularization's benefit may be specific to the actual noise realization Task 31C's real "
              f"hardware produced, not a general robustness improvement across plausible calibration "
              f"uncertainty.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N": N, "t_per_eval_s": t_per_eval,
            "errs_std": errs_std.tolist(), "errs_joint": errs_joint.tolist(), "chi2dofs": chi2dofs.tolist(),
            "quantiles_std": q_std, "quantiles_joint": q_joint,
            "n_outliers_std": n_outliers_std, "n_outliers_joint": n_outliers_joint,
            "n_joint_better": n_joint_better,
            "target_met_std": bool(q_std[95] < 0.25), "target_met_joint": bool(q_joint[95] < 0.25),
            "hardware_bar_met_std": bool(q_std[95] < 0.5), "hardware_bar_met_joint": bool(q_joint[95] < 0.5),
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
