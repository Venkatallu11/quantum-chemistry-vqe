#!/usr/bin/env python3
"""
task40_robustness_envelope_ablation.py -- iteration 40, extends the
point-estimate ablation matrix (Task 40G, real data) to the ROBUSTNESS
level: for the SAME theta_true draws, is the ~6000x Q95 collapse
(18.29 -> ~0.003 kcal/mol) coming from the joint Schmidt frame fit
ALONE, from the QED+conditioned-PEC+GPi2 correction ALONE, or genuinely
from the combination -- exactly the open question flagged (not yet
answered) in the robustness-envelope entry itself.

THREE variants, SAME theta_true draws (paired comparison, not
independent samples -- removes across-variant draw-to-draw noise from
the comparison):
  A) FRAME ONLY  -- fit the joint frame directly on TRUE RAW (no PEC
     correction, no QED conditioning at all) values.
  B) CORRECTION ONLY (no frame) -- QED+conditioned-PEC+GPi2 correction,
     standard (non-frame) energy combination -- matches Task 40G's
     point-estimate ablation, now across the FULL draw ensemble.
  C) FULL (correction + frame) -- the actual candidate pipeline, already
     measured at Q95~0.003-0.0035 kcal/mol (this task's own local
     confirmation, both prior runs).

Run:
    PYTHONHASHSEED=0 python vqe/task40_robustness_envelope_ablation.py
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from task30b_pec_application import build_full
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task39e_conditioned_correction import condition_on_even_weight
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
P_GPI2_ASSUMED = 0.0005
M_DRAWS = 12  # first pass -- 2 DM propagations/slot/draw (raw + corrected) doubles cost vs the
# frame-only robustness runs, scoped down accordingly
N_RESTARTS = 3

THETA_TRUE_MEAN = {"p_zz": 0.014593, "p_gpi2": 0.0003, "delta_zz": 0.0, "delta_gpi2": 0.0}
THETA_TRUE_STD = {"p_zz": 0.000124, "p_gpi2": 0.00081, "delta_zz": 0.000249, "delta_gpi2": 0.00161}


def true_raw_values(angles, theta_true, labels):
    """TRUE noisy expectation values, NO correction, NO conditioning --
    the pure 'what would we measure with no mitigation at all' baseline."""
    p_zz_true, p_gpi2_true, delta_zz_true, delta_gpi2_true = (
        theta_true["p_zz"], theta_true["p_gpi2"], theta_true["delta_zz"], theta_true["delta_gpi2"])
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native
    qc = to_native(build_ansatz(angles), GATE_NAME)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == GATE_NAME:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz_true))
            p_here, n_here = p_zz_true, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))
            p_here, n_here = GPI_REAL_MEAN, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2_true))
            p_here, n_here = p_gpi2_true, 1
        else:
            dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm = dm.evolve(U, qargs=qargs)
        dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p_here, n_here))
    return {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm.data))) for l in labels}


def true_corrected_conditioned_values(angles, theta_true, labels):
    """SAME as task40_robustness_envelope_new_pipeline's own function --
    TRUE noise + PEC-inverse (assumed correction) + even-weight
    conditioning."""
    p_zz_true, p_gpi2_true, delta_zz_true, delta_gpi2_true = (
        theta_true["p_zz"], theta_true["p_gpi2"], theta_true["delta_zz"], theta_true["delta_gpi2"])
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native
    qc = to_native(build_ansatz(angles), GATE_NAME)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == GATE_NAME:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz_true))
            p_true, p_assumed, n_here = p_zz_true, ZZ_ASSUMED, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))
            p_true, p_assumed, n_here = GPI_REAL_MEAN, GPI_REAL_MEAN, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2_true))
            p_true, p_assumed, n_here = p_gpi2_true, P_GPI2_ASSUMED, 1
        else:
            dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm = dm.evolve(U, qargs=qargs)
        dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p_true, n_here))
        dm = apply_pauli_mixture(dm, qargs, pec_inverse_weights(p_assumed, n_here))
    dm_cond, retained = condition_on_even_weight(dm)
    return {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_cond.data))) for l in labels}


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task40_robustness_envelope_ablation.py -- FRAME ONLY vs CORRECTION ONLY vs FULL, same draws")
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

    rng = np.random.default_rng(42)  # yet another fresh seed
    theta_true_draws = []
    for _ in range(M_DRAWS):
        draw = {}
        for name in ["p_zz", "p_gpi2", "delta_zz", "delta_gpi2"]:
            v = float(rng.normal(THETA_TRUE_MEAN[name], THETA_TRUE_STD[name]))
            if name in ("p_zz", "p_gpi2") and v < 0:
                v = 0.0
            draw[name] = v
        theta_true_draws.append(draw)
    print(f"\n  {M_DRAWS} theta_true draws (seed=42, fresh)")

    errs_frame_only, errs_correction_only, errs_full = [], [], []
    t0 = time.time()
    for i, theta_true in enumerate(theta_true_draws):
        t_draw0 = time.time()
        raw_kept, corrected_kept = {}, {}
        for name in kept:
            raw_kept[name] = true_raw_values(fixed_solutions[name]["angles"], theta_true, non_id_labels)
            vals = true_corrected_conditioned_values(fixed_solutions[name]["angles"], theta_true, non_id_labels)
            corrected_kept[name] = {l: max(-1.0, min(1.0, v)) for l, v in vals.items()}

        # A) FRAME ONLY: fit frame on RAW (uncorrected) values
        rng_fit = np.random.default_rng(2000 + i)
        U_hat_raw, _, _ = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_kept, weight_unit,
                                            rng_fit, n_restarts=N_RESTARTS)
        full_raw_frame = build_full_from_frame(U_hat_raw, P_S, K, non_id_labels, kept)
        err_a = energy_and_err(p, full_raw_frame, K)

        # B) CORRECTION ONLY, no frame: standard combine
        full_corrected_std = build_full(corrected_kept, diag, K, non_id_labels)
        err_b = energy_and_err(p, full_corrected_std, K)

        # C) FULL: correction + frame
        rng_fit2 = np.random.default_rng(3000 + i)
        U_hat_c, _, _ = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, corrected_kept, weight_unit,
                                          rng_fit2, n_restarts=N_RESTARTS)
        full_c_frame = build_full_from_frame(U_hat_c, P_S, K, non_id_labels, kept)
        err_c = energy_and_err(p, full_c_frame, K)

        errs_frame_only.append(err_a)
        errs_correction_only.append(err_b)
        errs_full.append(err_c)
        print(f"    draw {i+1}/{M_DRAWS}: FRAME-ONLY={abs(err_a):.4f}  CORRECTION-ONLY={abs(err_b):.4f}  "
              f"FULL={abs(err_c):.4f} kcal/mol  ({time.time()-t_draw0:.1f}s)")

    print(f"\n  total time: {time.time()-t0:.1f}s")

    for name, errs in [("FRAME ONLY (no correction)", errs_frame_only),
                        ("CORRECTION ONLY (no frame)", errs_correction_only),
                        ("FULL (correction + frame)", errs_full)]:
        abs_errs = np.abs(errs)
        q50, q95 = float(np.percentile(abs_errs, 50)), float(np.percentile(abs_errs, 95))
        print(f"  {name:<32} Q50={q50:.4f}  Q95={q95:.4f}  max={abs_errs.max():.4f} kcal/mol")

    print(f"\n  -- HONEST READ --")
    print(f"  If FRAME ONLY's Q95 is close to FULL's, the frame fit is doing nearly all the work and "
          f"QED+conditioned-PEC is not essential for robustness (though Task 40G's point-estimate ablation "
          f"already showed it matters for the nominal value). If FRAME ONLY's Q95 is much worse than FULL's "
          f"(closer to the old pipeline's 18.29), the correction is doing real, essential work at the "
          f"robustness level too, and the combination is genuinely synergistic, not frame-fit-dominated.")


if __name__ == "__main__":
    main()
