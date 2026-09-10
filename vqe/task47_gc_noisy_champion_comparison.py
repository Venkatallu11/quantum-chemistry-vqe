#!/usr/bin/env python3
"""
task47_gc_noisy_champion_comparison.py -- iteration 47. The real A/B this
project's own framework called for: QWC+champion vs GC+champion, same
H4 slots, same shots-equivalent (fully analytic, exact-population
convention -- matching task40_robustness_envelope_new_pipeline.py's own
established method), same theta_true noise draws, same frozen assumed
correction, same joint-Schmidt-frame fit.

WHY THIS NEEDED A NEW DERIVATION, not just reuse: the existing analytic
conditioned correction (task39e.analytic_A_and_B_conditioned /
task40_robustness_envelope_new_pipeline.simulate_true_then_correct_conditioned)
propagates noise ONLY through the state-prep gates, then conditions on
even weight, then reads Pauli(l) directly off the conditioned density
matrix -- implicitly treating the MEASUREMENT circuit (ancilla-parity
CNOTs + basis change) as noiseless. That is exact for QWC's own basis
change (0 extra 2q gates -- a real, not approximate, no-op). It would be
the WRONG model for GC's diagonalizer, which adds 3-5 REAL extra native
2q gates for 3 of the 4 groups (task46's own verified diagonalizer 2q
counts: [0, 3, 3, 5]).

THE FIX: continue the SAME gate-by-gate noisy propagation (identical
depolarizing + PEC-inverse-weight treatment already used for state-prep)
through the diagonalizer's OWN native-gate-transpiled circuit, applied
to the region AFTER even-weight conditioning (physically correct order:
state-prep -> ancilla-parity CNOTs [conditioning encodes this exactly,
task39e's own derivation] -> basis-change/diagonalizer -> measure).
Reconstruct each original label's expectation via Pauli(item.diagonal)
(a real, valid I/Z Pauli operator -- Qiskit handles its own qubit
indexing internally, sidestepping task46's earlier string-indexing bug
entirely) times item.sign.

REGRESSION CHECK, run first: GC group 0 (the trivial group -- 10 labels
that are ALREADY pure Z-strings, 0 extra 1q/2q diagonalizer gates) must
give IDENTICAL results to the existing, unmodified QWC-style correction,
since a 0-gate diagonalizer changes nothing physically.

Run:
    PYTHONHASHSEED=0 python vqe/task47_gc_noisy_champion_comparison.py
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
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task39e_conditioned_correction import condition_on_even_weight
from task40_robustness_envelope_new_pipeline import simulate_true_then_correct_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from fixed_ansatz import build_ansatz
from native_stateprep import to_native
from general_commuting_measurements import build_general_commuting_measurement_plan
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
P_GPI2_ASSUMED = 0.0005
M_DRAWS = 15
N_RESTARTS = 3

THETA_TRUE_MEAN = {"p_zz": 0.014593, "p_gpi2": 0.0003, "delta_zz": 0.0, "delta_gpi2": 0.0}
THETA_TRUE_STD = {"p_zz": 0.000124, "p_gpi2": 0.00081, "delta_zz": 0.000249, "delta_gpi2": 0.00161}


def _build_state_prep_dm(angles, gate_name, theta_true):
    p_zz_true, p_gpi2_true, delta_zz_true, delta_gpi2_true = (
        theta_true["p_zz"], theta_true["p_gpi2"], theta_true["delta_zz"], theta_true["delta_gpi2"])
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
    return dm


def _apply_noisy_native_circuit(dm, native_circuit, gate_name, theta_true):
    p_zz_true, p_gpi2_true, delta_zz_true, delta_gpi2_true = (
        theta_true["p_zz"], theta_true["p_gpi2"], theta_true["delta_zz"], theta_true["delta_gpi2"])
    for instr in native_circuit.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [native_circuit.find_bit(q).index for q in instr.qubits]
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
    return dm


def simulate_true_then_correct_conditioned_gc(angles, gate_name, theta_true, native_diag_circuit, transformed):
    """GC version: state-prep (noisy) -> condition on even weight (encodes
    the real ancilla-0 postselection) -> CONTINUE noisy propagation through
    the diagonalizer's own native-transpiled gates -> read off each
    original label via its I/Z diagonal representative + sign."""
    dm = _build_state_prep_dm(angles, gate_name, theta_true)
    dm_cond, retained = condition_on_even_weight(dm)
    dm_after = _apply_noisy_native_circuit(dm_cond, native_diag_circuit, gate_name, theta_true)
    out = {}
    for item in transformed:
        Z_op = np.asarray(Pauli(item.diagonal).to_matrix())
        val = float(np.real(np.trace(Z_op @ dm_after.data)))
        out[item.original] = item.sign * val
    return out, retained


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task47_gc_noisy_champion_comparison.py -- REAL A/B: QWC+champion vs GC+champion, same draws")
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

    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    native_diags = [to_native(d.to_circuit(), GATE_NAME) for d in diagonalizers]
    for i, (g, d, nd) in enumerate(zip(groups, diagonalizers, native_diags)):
        counts = dict(nd.count_ops())
        print(f"  GC group {i}: {len(g)} labels, abstract 2q={d.two_qubit_count}, "
              f"NATIVE gate counts after to_native(): {counts}")

    print("\n  -- REGRESSION CHECK: GC group with 0 extra gates must exactly match the existing QWC-style "
          "conditioned correction (a 0-gate diagonalizer changes nothing physically) --")
    trivial_idx = [i for i, d in enumerate(diagonalizers) if d.two_qubit_count == 0 and d.one_qubit_count == 0]
    assert trivial_idx, "expected one GC group with a fully trivial (empty) diagonalizer"
    ti = trivial_idx[0]
    test_theta = {"p_zz": 0.014593, "p_gpi2": 0.0003, "delta_zz": 0.0001, "delta_gpi2": 0.0002}
    vals_gc, ret_gc = simulate_true_then_correct_conditioned_gc(
        fixed_solutions["u_0"]["angles"], GATE_NAME, test_theta, native_diags[ti], diagonalizers[ti].transformed)
    vals_qwc, ret_qwc = simulate_true_then_correct_conditioned(
        fixed_solutions["u_0"]["angles"], GATE_NAME, test_theta, groups[ti])
    max_diff = max(abs(vals_gc[l] - vals_qwc[l]) for l in groups[ti])
    print(f"    max diff (trivial-group GC vs existing QWC-style correction): {max_diff:.3e}   "
          f"retained: gc={ret_gc:.6f} qwc={ret_qwc:.6f}   "
          f"{'PASS' if max_diff < 1e-9 and abs(ret_gc - ret_qwc) < 1e-9 else 'FAIL -- STOP, do not trust this derivation'}")
    if max_diff >= 1e-9:
        raise RuntimeError("GC correction does not reduce to the existing QWC correction on a 0-gate group -- bug, stop")

    print(f"\n  -- REAL A/B: {M_DRAWS} theta_true draws, N_RESTARTS={N_RESTARTS} (matching task40's own "
          f"established convention exactly) --")
    rng = np.random.default_rng(41)
    theta_true_draws = []
    for _ in range(M_DRAWS):
        draw = {}
        for name in ["p_zz", "p_gpi2", "delta_zz", "delta_gpi2"]:
            v = float(rng.normal(THETA_TRUE_MEAN[name], THETA_TRUE_STD[name]))
            if name in ("p_zz", "p_gpi2") and v < 0:
                v = 0.0
            draw[name] = v
        theta_true_draws.append(draw)

    t0 = time.time()
    errs_qwc, errs_gc = [], []
    for i, theta_true in enumerate(theta_true_draws):
        t_draw0 = time.time()

        raw_qwc = {}
        for name in kept:
            vals, retained = simulate_true_then_correct_conditioned(fixed_solutions[name]["angles"], GATE_NAME,
                                                                        theta_true, non_id_labels)
            raw_qwc[name] = {l: max(-1.0, min(1.0, v)) for l, v in vals.items()}
        rng_fit = np.random.default_rng(2000 + i)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_qwc, weight_unit,
                                                  rng_fit, n_restarts=N_RESTARTS)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        err_qwc = energy_and_err(p, full, K)
        errs_qwc.append(err_qwc)

        raw_gc = {name: {} for name in kept}
        for name in kept:
            for g, d, nd in zip(groups, diagonalizers, native_diags):
                vals, retained = simulate_true_then_correct_conditioned_gc(
                    fixed_solutions[name]["angles"], GATE_NAME, theta_true, nd, d.transformed)
                for l, v in vals.items():
                    raw_gc[name][l] = max(-1.0, min(1.0, v))
        rng_fit2 = np.random.default_rng(3000 + i)
        U_hat2, cost2, chi2dof2 = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_gc, weight_unit,
                                                     rng_fit2, n_restarts=N_RESTARTS)
        full2 = build_full_from_frame(U_hat2, P_S, K, non_id_labels, kept)
        err_gc = energy_and_err(p, full2, K)
        errs_gc.append(err_gc)

        print(f"    draw {i+1}/{M_DRAWS}: QWC |err|={abs(err_qwc):.4f}  GC |err|={abs(err_gc):.4f} kcal/mol  "
              f"({time.time()-t_draw0:.1f}s)")

    total_time = time.time() - t0
    print(f"\n  total time: {total_time:.1f}s ({total_time/M_DRAWS:.1f}s/draw)")

    for label, errs in [("QWC (existing, 13 groups, 273 circuits)", errs_qwc),
                         ("GC (new, 4 groups, 84 circuits)", errs_gc)]:
        abs_errs = np.abs(errs)
        q50, q90, q95 = [float(np.percentile(abs_errs, q)) for q in [50, 90, 95]]
        print(f"\n  {label}: Q50={q50:.4f}  Q90={q90:.4f}  Q95={q95:.4f} kcal/mol")
        print(f"    errs sorted: {np.array2string(np.sort(abs_errs), precision=4, max_line_width=200)}")

    print(f"\n  -- PAIRED per-draw comparison (same theta_true draw for both) --")
    for i, (eq, eg) in enumerate(zip(errs_qwc, errs_gc)):
        print(f"    draw {i+1}: QWC={abs(eq):.4f}  GC={abs(eg):.4f}  ratio(GC/QWC)={abs(eg)/max(abs(eq),1e-9):.2f}x")


if __name__ == "__main__":
    main()
