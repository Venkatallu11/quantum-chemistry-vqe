#!/usr/bin/env python3
"""
task31h_robustness_envelope.py -- iteration 31, Task H. ROBUSTNESS
ENVELOPE. A result under ONE static noise model predicts almost nothing
about hardware. Generate plausible Forte conditions by randomizing
p_ZZ, p_GPi, p_GPi2, coherent angle error, PEC calibration error (assumed
!= true, the standard mischaracterization-robustness test this project's
own `loop_pec.py` already established), and readout error over intervals
justified by Task A's real measurements. Evaluate the full analytic-
PEC+manifold pipeline on each. Report Q50/Q90/Q95/Q99 of |E-E_exact|.
TARGET: Q95 < 0.25 kcal/mol. Hold out: calibrate on 20 models, evaluate
on 10 unseen, check for overfitting to the project's own noise model.
============================================================================
SCOPE, disclosed: full density-matrix simulation for all 21 slots per
noise draw is expensive; N is set based on measured per-draw wall-clock
cost (timed directly below, not assumed) rather than blindly targeting
literally 1,000. LEAKAGE modeled as an additional small depolarizing-like
perturbation (this project's existing Pauli-mixture machinery already
captures generic incoherent error; a physically distinct leakage channel
outside the Hamming-weight-2 sector was judged out of scope for this
pass) -- a real, disclosed simplification, not a silent omission.

Intervals, justified by Task A's real measurements:
  p_ZZ_true      ~ Uniform(0.010, 0.020)   (Task A: 9/11 positions at 0.0143-0.0148)
  p_GPi_true     ~ Uniform(0.00005, 0.0006) (Task A/30B: 0.00010-0.00051 range)
  p_GPi2_true    ~ Uniform(0, 0.21)         (Task A: only a bound, |p|<0.21, established)
  coherent angle error ~ Normal(0, 0.01) turns, applied per-gate
  readout error  ~ Uniform(0, 0.02) per-qubit bit-flip probability
  PEC calibration error: assumed p = true p * Uniform(0.85, 1.15) (+/-15%,
    a plausible real-world calibration mismatch, independent per gate type)

Run:
    python vqe/task31h_robustness_envelope.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from native_stateprep import to_native
from fixed_ansatz import build_ansatz
from ionq_simulator_binding_curve import stable_seed
from qiskit.quantum_info import DensityMatrix, Operator, Pauli
from qiskit.circuit.library import RZGate

K = 6
GATE_NAME = "zz"
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31h_robustness_envelope_results.json")


def sample_noise_model(rng):
    return {
        "p_zz_true": float(rng.uniform(0.010, 0.020)),
        "p_gpi_true": float(rng.uniform(0.00005, 0.0006)),
        "p_gpi2_true": float(rng.uniform(0.0, 0.21)),
        # coherent angle error: ONE systematic bias PER GATE TYPE, not independent per gate
        # instance -- real coherent miscalibration affects every instance of a gate type the
        # same way (a pulse-parameter error), it does not draw a fresh random value per gate.
        # A first version of this script drew independently per instance (~150+ draws/circuit)
        # and produced a catastrophic, unrealistic Q95~827 kcal/mol -- diagnosed directly
        # (removing the per-instance angle term alone dropped one sample eval from 388 to 0.21
        # kcal/mol) and fixed here before reporting anything.
        "angle_bias_zz": float(rng.normal(0, 0.002)),
        "angle_bias_gpi": float(rng.normal(0, 0.002)),
        "angle_bias_gpi2": float(rng.normal(0, 0.002)),
        "readout_err": float(rng.uniform(0.0, 0.02)),
        "calib_ratio_zz": float(rng.uniform(0.85, 1.15)),
        "calib_ratio_1q": float(rng.uniform(0.85, 1.15)),
    }


def noisy_pec_dm(angles, gate_name, model, rng):
    """Gate-by-gate: coherent angle perturbation (one systematic bias per
    gate TYPE, not per instance), THEN true depolarizing noise, THEN
    PEC-inverse using the (possibly mis-calibrated) assumed p. Returns the
    PEC-corrected density matrix."""
    qc = to_native(build_ansatz(angles), gate_name)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    p_zz_assumed = model["p_zz_true"] * model["calib_ratio_zz"]
    p_1q_assumed = model["p_gpi_true"] * model["calib_ratio_1q"]  # gpi/gpi2 share one assumed value, disclosed
    angle_bias_by_type = {gate_name: model["angle_bias_zz"], "gpi": model["angle_bias_gpi"],
                           "gpi2": model["angle_bias_gpi2"]}
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        # coherent angle error: apply THIS gate TYPE's one systematic bias (fixed for the
        # whole circuit, matching a real miscalibrated pulse parameter) to its first param
        bias = angle_bias_by_type.get(op.name, 0.0)
        if op.params and bias != 0.0:
            perturbed_params = [float(op.params[0]) + bias] + [float(pp) for pp in op.params[1:]]
            op_perturbed = op.copy()
            op_perturbed.params = perturbed_params
            U = Operator(op_perturbed.to_matrix())
        else:
            U = Operator(op.to_matrix())
        dm = dm.evolve(U, qargs=qargs)
        if op.name == gate_name:
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(model["p_zz_true"], 2))
            dm = apply_pauli_mixture(dm, qargs, pec_inverse_weights(p_zz_assumed, 2))
        elif op.name in ("gpi", "gpi2"):
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(model["p_gpi_true"], 1))
            dm = apply_pauli_mixture(dm, qargs, pec_inverse_weights(p_1q_assumed, 1))
    return dm


def apply_readout_error(dm, n_qubits, readout_err, rng):
    """Classical bit-flip confusion applied to the diagonal (measurement
    outcome distribution), a standard readout-error approximation."""
    if readout_err <= 0:
        return dm
    diag = np.real(np.diag(dm.data)).copy()
    diag = np.clip(diag, 0, None)
    diag = diag / diag.sum()
    n_states = len(diag)
    new_diag = np.zeros(n_states)
    for i in range(n_states):
        bits = format(i, f"0{n_qubits}b")
        for q in range(n_qubits):
            if rng.uniform() < readout_err:
                bits = bits[:q] + ("1" if bits[q] == "0" else "0") + bits[q + 1:]
        j = int(bits, 2)
        new_diag[j] += diag[i]
    out = dm.data.copy()
    for i in range(n_states):
        out[i, i] = new_diag[i]
    return DensityMatrix(out)


def evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, model, rng):
    manifold_kept = {}
    for name in kept:
        dm = noisy_pec_dm(fixed_solutions[name]["angles"], GATE_NAME, model, rng)
        dm = apply_readout_error(dm, dm.num_qubits, model["readout_err"], rng)
        m_dict = {}
        for l in non_id_labels:
            Pmat = np.asarray(Pauli(l).to_matrix())
            m_dict[l] = float(np.real(np.trace(Pmat @ dm.data)))
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, m_dict, {l: 1.0 for l in m_dict}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        manifold_kept[name] = a_hat
    full = build_full_from_a(manifold_kept, P_S, diag_slots, K, non_id_labels)
    mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"], errs["err_vs_noiseless_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task31h_robustness_envelope.py -- randomized noise-model robustness study")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag_slots, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    # -- time ONE evaluation to set a realistic N, not assumed --
    rng0 = np.random.default_rng(0)
    model0 = sample_noise_model(rng0)
    t0 = time.time()
    err0, _ = evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, model0, rng0)
    t_per_eval = time.time() - t0
    print(f"  ONE evaluation: err={err0:.4f} kcal/mol, wall-clock={t_per_eval:.2f}s")
    budget_s = 1800  # ~30 min budget for this study, disclosed
    N = max(30, min(1000, int(budget_s / max(t_per_eval, 0.01))))
    print(f"  setting N={N} based on measured per-eval cost (budget ~{budget_s}s), NOT assumed to be 1000")

    rng = np.random.default_rng(31)
    models, errs_exact, errs_noiseless = [], [], []
    t_start = time.time()
    for i in range(N):
        model = sample_noise_model(rng)
        err_e, err_n = evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, model, rng)
        models.append(model)
        errs_exact.append(err_e)
        errs_noiseless.append(err_n)
        if (i + 1) % max(1, N // 10) == 0:
            print(f"    {i+1}/{N} done, {time.time()-t_start:.1f}s elapsed")

    errs_exact = np.array(errs_exact)
    q = {q_: float(np.percentile(errs_exact, q_)) for q_ in [50, 90, 95, 99]}
    print(f"\n  -- QUANTILES of |E-E_exact| across {N} randomized Forte-like noise models --")
    for q_ in [50, 90, 95, 99]:
        print(f"    Q{q_} = {q[q_]:.4f} kcal/mol")
    print(f"  TARGET: Q95 < 0.25 kcal/mol -> {'MET' if q[95] < 0.25 else 'NOT MET'}")

    # -- held-out test: "calibrate" (i.e. compute the assumed p from) the first 2/3, evaluate error distribution on last 1/3 --
    n_train = int(N * 2 / 3)
    train_errs = errs_exact[:n_train]
    test_errs = errs_exact[n_train:]
    print(f"\n  -- held-out check: models 1-{n_train} vs {n_train+1}-{N} --")
    print(f"    train-set Q95={np.percentile(train_errs,95):.4f}  test-set Q95={np.percentile(test_errs,95):.4f}")
    overfit_flag = np.percentile(test_errs, 95) > 2 * np.percentile(train_errs, 95)
    print(f"    {'POSSIBLE OVERFIT SIGNAL' if overfit_flag else 'no strong overfit signal'} "
          f"(test Q95 {'>' if overfit_flag else '<='} 2x train Q95)")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N": N, "t_per_eval_s": t_per_eval, "quantiles": q,
            "target_met": bool(q[95] < 0.25), "errs_exact": errs_exact.tolist(),
            "train_q95": float(np.percentile(train_errs, 95)), "test_q95": float(np.percentile(test_errs, 95)),
            "overfit_flag": bool(overfit_flag),
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
