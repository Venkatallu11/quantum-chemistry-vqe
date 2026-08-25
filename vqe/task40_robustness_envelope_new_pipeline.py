#!/usr/bin/env python3
"""
task40_robustness_envelope_new_pipeline.py -- iteration 40, the user's
own explicitly-stated "most important number in the entire project"
(section 8): does the NEW pipeline (native ancilla-parity QED detection
+ conditioned PEC + a fixed GPi2 correction + a freely-fit joint Schmidt
frame) actually reduce Q95 under REAL calibration uncertainty, or does
the spectacular 0.011-0.018 kcal/mol only hold at nominal parameters --
i.e. is this the right architecture (Q95 collapses too), or was the
combination still a narrow calibration point (Q95 stays large)?

METHOD, fully analytic (EXACT populations, no shot noise -- matching
this project's own established robustness-envelope convention,
`task31h_robustness_envelope.py`'s own `evaluate_one_model`, for
computational tractability; isolates calibration-uncertainty sensitivity
specifically, cleanly separable from the shot-noise question Task 40's
ratio-bias check already addressed): for M draws of theta_true (SAME
real, disclosed distribution Task 38C/38D used, tighter and more
evidence-based than the original wide independent priors), build the
TRUE noisy 4-qubit density matrix gate-by-gate, apply the PEC-inverse
using the FIXED, ALREADY-SELECTED assumed correction (ZZ_ASSUMED,
GPI_REAL_MEAN, P_GPI2_ASSUMED -- the frozen candidate, not re-optimized
per draw, matching the user's own "freeze the pipeline, don't tune
toward the answer" instruction, section 1), THEN condition on even
weight (Task 39E's projector, exactly what real ancilla=0 postselection
measures), THEN trace against every Pauli label needed. Fit the 15-
parameter Schmidt frame FREELY (no leakage) on the resulting (slot,
label) values for EACH draw, compute Q50/Q95/Q99 of |E-E_exact| across
draws -- directly comparable to Task 36's own real Q95=18.29 kcal/mol
number for the OLD pipeline (no QED, no GPi2 correction).

SCOPE, disclosed (matching this project's own "N set by measured cost"
discipline): a joint-frame fit per draw is the expensive step. Starts
with a SMALL M and REDUCED restarts to get a real, honest first number
and measure actual per-draw cost before committing to a larger run.

Run:
    PYTHONHASHSEED=0 python vqe/task40_robustness_envelope_new_pipeline.py
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
P_GPI2_ASSUMED = 0.0005  # frozen, representative of the two selected candidates (0.0004/0.0006) --
# ONE fixed value used for every draw, exactly matching the "freeze the pipeline" instruction
M_DRAWS = 25  # second, larger pass -- per-draw cost now measured (~11-16s) from the first N=15 run
N_RESTARTS = 3  # reduced from Task 39H's 4 (itself already reduced from Task 36's 12) for this
# first, cost-finding pass

THETA_TRUE_MEAN = {"p_zz": 0.014593, "p_gpi2": 0.0003, "delta_zz": 0.0, "delta_gpi2": 0.0}
THETA_TRUE_STD = {"p_zz": 0.000124, "p_gpi2": 0.00081, "delta_zz": 0.000249, "delta_gpi2": 0.00161}


def simulate_true_then_correct_conditioned(angles, gate_name, theta_true, labels):
    p_zz_true, p_gpi2_true, delta_zz_true, delta_gpi2_true = (
        theta_true["p_zz"], theta_true["p_gpi2"], theta_true["delta_zz"], theta_true["delta_gpi2"])
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native
    qc = to_native(build_ansatz(angles), gate_name)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == gate_name:
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
    return {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_cond.data))) for l in labels}, retained


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task40_robustness_envelope_new_pipeline.py -- Q95 of the NEW pipeline under real calibration uncertainty")
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

    rng = np.random.default_rng(41)  # DIFFERENT seed from the first run (40) -- an independent
    # confirmation draw set, not a re-run of the same random numbers
    theta_true_draws = []
    n_clipped = 0
    for _ in range(M_DRAWS):
        draw = {}
        for name in ["p_zz", "p_gpi2", "delta_zz", "delta_gpi2"]:
            v = float(rng.normal(THETA_TRUE_MEAN[name], THETA_TRUE_STD[name]))
            if name in ("p_zz", "p_gpi2") and v < 0:
                v, n_clipped = 0.0, n_clipped + 1
            draw[name] = v
        theta_true_draws.append(draw)
    print(f"\n  {M_DRAWS} theta_true draws (clipped {n_clipped} negative probability draws to 0)")
    print(f"  frozen assumed correction: ZZ={ZZ_ASSUMED}, GPi2={P_GPI2_ASSUMED} (fixed for every draw)")

    t0 = time.time()
    errs = []
    retained_fracs = []
    for i, theta_true in enumerate(theta_true_draws):
        t_draw0 = time.time()
        raw_kept = {}
        min_retained = 1.0
        for name in kept:
            vals, retained = simulate_true_then_correct_conditioned(fixed_solutions[name]["angles"], GATE_NAME,
                                                                        theta_true, non_id_labels)
            raw_kept[name] = {l: max(-1.0, min(1.0, v)) for l, v in vals.items()}
            min_retained = min(min_retained, retained)
        retained_fracs.append(min_retained)
        rng_fit = np.random.default_rng(1000 + i)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_kept, weight_unit,
                                                  rng_fit, n_restarts=N_RESTARTS)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        err = energy_and_err(p, full, K)
        errs.append(err)
        print(f"    draw {i+1}/{M_DRAWS}: |err|={abs(err):.4f} kcal/mol  chi2/dof={chi2dof:.5f}  "
              f"min_retained_frac={min_retained:.3f}  ({time.time()-t_draw0:.1f}s)")

    total_time = time.time() - t0
    print(f"\n  total time: {total_time:.1f}s  ({total_time/M_DRAWS:.1f}s/draw)")

    abs_errs = np.abs(errs)
    q50, q90, q95, q99 = [float(np.percentile(abs_errs, q)) for q in [50, 90, 95, 99]]
    print(f"\n  -- RESULT: Q50={q50:.4f}  Q90={q90:.4f}  Q95={q95:.4f}  Q99={q99:.4f} kcal/mol  (N={M_DRAWS}) --")
    print(f"  errs sorted: {np.array2string(np.sort(abs_errs), precision=3, max_line_width=200)}")
    print(f"\n  -- COMPARISON: Task 36's OLD pipeline (no QED, no GPi2 correction), real N=267 robustness envelope: Q95=18.29 kcal/mol --")
    if q95 < 0.5:
        print(f"  NEW pipeline Q95={q95:.4f} < 0.5 -- the architecture appears to genuinely close the robustness gap, "
          f"not just hit a narrow calibration point. This would be a major finding.")
    elif q95 < 5.0:
        print(f"  NEW pipeline Q95={q95:.4f} -- large real improvement over 18.29 but not yet under the 0.5 bar -- "
          f"a real, substantial, partial win.")
    else:
        print(f"  NEW pipeline Q95={q95:.4f} -- still large; the 0.011-0.018 result may be closer to a narrow "
          f"calibration point than a robustly generalizing architecture. Honest, disclosed either way.")


if __name__ == "__main__":
    main()
