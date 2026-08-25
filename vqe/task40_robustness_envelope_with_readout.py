#!/usr/bin/env python3
"""
task40_robustness_envelope_with_readout.py -- iteration 40, closes a
disclosed gap in the earlier robustness-envelope runs: p_readout was not
included as a varying noise dimension at all (the new pipeline has no
readout correction), while Task 36's own old Q95=18.29 comparison point
DID include readout uncertainty in its noise model -- a real,
acknowledged apples-to-oranges gap. This adds it honestly.

READOUT MODEL, derived not assumed: symmetric per-qubit bit-flip
readout error affects TWO things simultaneously, both handled here --
  1. The ANCILLA's own readout: a true ancilla=0 (valid) shot can be
     misread as 1 (falsely discarded, a real cost but not a correctness
     problem) and a true ancilla=1 (leaked) shot can be misread as 0
     (falsely KEPT, contaminating the postselected sample -- the
     concerning direction). Modeled exactly via a 2x2 stochastic
     bit-flip mixing applied ONLY to the ancilla's own diagonal
     populations (grouping by the register's own bitstring, mixing the
     two ancilla outcomes for each), BEFORE conditioning on the
     (now readout-noisy) ancilla=0 outcome.
  2. The REGISTER's own readout, for the actual Pauli label being
     measured: Task 37C's own already-verified `readout_attenuation`
     multiplicative factor (1-2p)^weight, applied to the final traced
     expectation value -- unchanged, reused directly, not re-derived.

p_readout itself is drawn from Task 37B's own weak, disclosed-as-
uninformative prior (mean=0.005, bounds [0,0.02]) -- no real
calibration for it exists in this project (Task 38's own attempt was
blocked by a real API quota error), so this deliberately does NOT
pretend more precision than exists; it directly tests SENSITIVITY to a
plausible, honestly-sourced range.

Run:
    PYTHONHASHSEED=0 python vqe/task40_robustness_envelope_with_readout.py
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
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix, readout_attenuation
from task37b_h4_noise_model import GPI_REAL_MEAN, PRIORS
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
P_GPI2_ASSUMED = 0.0005
M_DRAWS = 15
N_RESTARTS = 3

THETA_TRUE_MEAN = {"p_zz": 0.014593, "p_gpi2": 0.0003, "delta_zz": 0.0, "delta_gpi2": 0.0}
THETA_TRUE_STD = {"p_zz": 0.000124, "p_gpi2": 0.00081, "delta_zz": 0.000249, "delta_gpi2": 0.00161}


def ancilla_readout_channel(dm5, p_ro):
    """BUG FOUND AND FIXED HERE, disclosed: an earlier version mixed only
    the DIAGONAL populations by hand, discarding every register-qubit
    off-diagonal coherence -- confirmed catastrophic by direct test
    (p_readout~0 draws, which should reduce to the already-verified
    readout-free case, instead gave ~18 kcal/mol, matching the OLD
    uncorrected baseline almost exactly -- the coherence Task 39E's own
    exact projector conditioning depends on was being thrown away before
    conditioning even ran). Fixed with a REAL Kraus-operator bit-flip
    channel on the ancilla qubit only (X_4 rho X_4 with probability p_ro,
    identity otherwise) -- a proper quantum channel that preserves every
    register coherence exactly, applied via DensityMatrix.evolve on
    qubit index 4 alone."""
    from qiskit.quantum_info import Kraus
    K0 = np.sqrt(1 - p_ro) * np.eye(2)
    K1 = np.sqrt(p_ro) * np.array([[0, 1], [1, 0]])
    kraus = Kraus([K0, K1])
    return dm5.evolve(kraus, qargs=[4])


def project_ancilla_zero(dm5):
    """Exact projection onto the ancilla=0 subspace (|0><0| on qubit 4,
    identity on the register), renormalized -- preserves every register
    coherence within that subspace exactly, matching Task 39E's own
    `condition_on_even_weight` construction (same idea, now conditioning
    on ancilla=0 directly since readout noise has already been applied
    to the ancilla, rather than on register weight-parity)."""
    proj0 = np.array([[1, 0], [0, 0]], dtype=complex)
    # .expand() places its ARGUMENT at the HIGHEST qubit index (the SAME convention already used
    # and verified for building the 5-qubit state itself: dm4.expand(ancilla)) -- using expand()
    # here instead of tensor() avoids re-deriving qiskit's tensor-vs-expand qubit-order relation
    P = Operator(np.eye(2**4)).expand(Operator(proj0))
    rho = dm5.data
    P_mat = np.asarray(P)
    rho_proj = P_mat @ rho @ P_mat
    norm = float(np.real(np.trace(rho_proj)))
    if norm < 1e-12:
        norm = 1e-12
    return DensityMatrix(rho_proj / norm), norm


def simulate_true_then_correct_conditioned_readout(angles, theta_true, p_readout_true, labels):
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

    # append ancilla (register dm -> 5-qubit joint state via the SAME expand()+CNOT convention
    # verified in Task 40's ratio-bias check)
    from task39b_native_ancilla_parity import ancilla_cnots_abstract
    dm5 = DensityMatrix(dm.data).expand(DensityMatrix.from_label("0"))
    ancilla_qc = ancilla_cnots_abstract()
    for instr in ancilla_qc.data:
        op = instr.operation
        qargs = [ancilla_qc.find_bit(q).index for q in instr.qubits]
        dm5 = dm5.evolve(Operator(op.to_matrix()), qargs=qargs)

    # ancilla readout error as a REAL quantum channel (preserves register coherence exactly),
    # THEN exact projection onto the (now readout-noisy) ancilla=0 outcome, THEN trace out the
    # ancilla to get the conditional 4-qubit register state
    from qiskit.quantum_info import partial_trace
    dm5_noisy = ancilla_readout_channel(dm5, p_readout_true)
    dm5_cond, retained = project_ancilla_zero(dm5_noisy)
    dm4_cond = partial_trace(dm5_cond, [4])

    out = {}
    for l in labels:
        register_val = float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ np.asarray(dm4_cond))))
        out[l] = register_val * readout_attenuation(l, p_readout_true)
    return out, retained


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task40_robustness_envelope_with_readout.py -- closing the disclosed p_readout gap")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag_slots, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    rng = np.random.default_rng(43)
    theta_true_draws = []
    p_readout_pr = PRIORS["p_readout"]
    for _ in range(M_DRAWS):
        draw = {}
        for name in ["p_zz", "p_gpi2", "delta_zz", "delta_gpi2"]:
            v = float(rng.normal(THETA_TRUE_MEAN[name], THETA_TRUE_STD[name]))
            if name in ("p_zz", "p_gpi2") and v < 0:
                v = 0.0
            draw[name] = v
        p_ro = float(np.clip(rng.normal(p_readout_pr["mean"], p_readout_pr["std"]), *p_readout_pr["bounds"]))
        draw["p_readout"] = p_ro
        theta_true_draws.append(draw)
    print(f"\n  {M_DRAWS} draws (seed=43), p_readout ~ N({p_readout_pr['mean']}, {p_readout_pr['std']}) "
          f"bounded {p_readout_pr['bounds']} (Task 37B's own weak, disclosed prior)")

    print("\n  -- REGRESSION CHECK: at p_readout=0, must exactly match the already-validated "
          "readout-free function (task40_robustness_envelope_new_pipeline) --")
    from task40_robustness_envelope_new_pipeline import simulate_true_then_correct_conditioned as _old_fn
    test_theta = theta_true_draws[0]
    test_labels = non_id_labels[:4]
    vals_new, ret_new = simulate_true_then_correct_conditioned_readout(
        fixed_solutions["u_0"]["angles"], test_theta, 0.0, test_labels)
    vals_old, ret_old = _old_fn(fixed_solutions["u_0"]["angles"], GATE_NAME, test_theta, test_labels)
    max_diff = max(abs(vals_new[l] - vals_old[l]) for l in test_labels)
    ret_diff = abs(ret_new - ret_old)
    print(f"    max value diff: {max_diff:.3e}   retained diff: {ret_diff:.3e}   "
          f"{'PASS' if max_diff < 1e-6 and ret_diff < 1e-6 else 'FAIL -- STOP, do not trust this run'}")
    if max_diff >= 1e-6 or ret_diff >= 1e-6:
        raise RuntimeError("readout-inclusive function does not reduce to the validated readout-free "
                            "function at p_readout=0 -- real bug, stop before trusting any numbers below")

    t0 = time.time()
    errs = []
    for i, theta_true in enumerate(theta_true_draws):
        t_draw0 = time.time()
        raw_kept = {}
        min_retained = 1.0
        for name in kept:
            vals, retained = simulate_true_then_correct_conditioned_readout(
                fixed_solutions[name]["angles"], theta_true, theta_true["p_readout"], non_id_labels)
            raw_kept[name] = {l: max(-1.0, min(1.0, v)) for l, v in vals.items()}
            min_retained = min(min_retained, retained)
        rng_fit = np.random.default_rng(4000 + i)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_kept, weight_unit,
                                                  rng_fit, n_restarts=N_RESTARTS)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        err = energy_and_err(p, full, K)
        errs.append(err)
        print(f"    draw {i+1}/{M_DRAWS}: |err|={abs(err):.4f} kcal/mol  p_readout={theta_true['p_readout']:.5f}  "
              f"min_retained={min_retained:.3f}  ({time.time()-t_draw0:.1f}s)")

    print(f"\n  total time: {time.time()-t0:.1f}s")
    abs_errs = np.abs(errs)
    q50, q90, q95, q99 = [float(np.percentile(abs_errs, q)) for q in [50, 90, 95, 99]]
    print(f"\n  -- RESULT (WITH p_readout uncertainty): Q50={q50:.4f}  Q90={q90:.4f}  Q95={q95:.4f}  Q99={q99:.4f} kcal/mol --")
    print(f"  errs sorted: {np.array2string(np.sort(abs_errs), precision=4, max_line_width=200)}")
    print(f"\n  -- COMPARISON --")
    print(f"  WITHOUT readout uncertainty (2 prior runs, N=15/N=25): Q95~0.0031-0.0035 kcal/mol")
    print(f"  WITH readout uncertainty (this run, N={M_DRAWS}):      Q95={q95:.4f} kcal/mol")
    if q95 < 0.5:
        print(f"  Still comfortably under the bars even with readout uncertainty included -- a real, "
              f"more complete confirmation, not just a favorable simplification.")
    else:
        print(f"  Readout uncertainty meaningfully degrades the result -- a real, disclosed limitation "
              f"the earlier readout-free runs did not capture.")


if __name__ == "__main__":
    main()
