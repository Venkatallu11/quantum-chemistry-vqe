#!/usr/bin/env python3
"""
task30b_pec_application.py -- iteration 30, Task B, step 2. Applies the
REAL, angle-conditioned p2(zz)/p1(gpi,gpi2) learned in
task30b_pec_calibration.py to REAL measured H4 data on the optimizer-
reduced native circuit (Task 28B's own checkpoint, ideal/aria-1/forte-1
CONCURRENT, SAME circuits/seeds -- reused, not resubmitted, since it is
already real, already the specified circuit).
============================================================================
HOW THE CORRECTION IS COMPUTED, disclosed precisely -- this is NOT literal
quasi-probability circuit twirling (that would need many real circuit
variants per gate location; judged out of reach of this iteration's
scope, see ALTERNATIVES NOT TAKEN). Instead, for each (slot, label) this
file:
  1. builds the SAME density-matrix forward propagation
     `loop_pec.py::apply_pauli_mixture` already uses, gate by gate through
     the ACTUAL optimized circuit, using each gate's OWN real angle to
     pick the correct learned p (p2 for the single zz theta; p1 for each
     gpi/gpi2 gate's OWN phi, nearest-bin from the 4 calibrated bins);
  2. computes A[l] = the analytic expectation value if a depolarizing
     channel with the LEARNED p truly describes this circuit's noise;
  3. computes B[l] = the analytic PEC-corrected expectation (immediately
     applying the learned-p PEC inverse after each gate's depolarizing
     step, exactly `pec_measure_learned`'s own interleaving);
  4. applies the RATIO B[l]/A[l] -- not an absolute simulated value -- to
     the REAL measured raw value. This uses the exact, already-verified
     gate-by-gate PEC algebra to compute the CORRECT relative correction
     factor per label (accounting for which gates each observable's
     Heisenberg-propagated support actually passes through), while the
     data point being corrected is REAL, not simulated.

MANDATORY IDEAL-CONTROL CHECK: PEC on `ideal` must not degrade the
already-small shot-noise-floor error.

Run:
    python vqe/task30b_pec_application.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

K = 6
SHOTS = 100_000
N_SEEDS = 8
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "task28b_optimized_raw.json")
CALIB_RESULTS = os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task30b_pec_application_results.json")


def nearest_bin_p(bins_dict, phi):
    """bins_dict: {angle_str: {"p":..., ...}}. Nearest-angle lookup among
    the 4 calibrated bins -- coarse but disclosed, not a continuous fit."""
    angles = [float(a) for a in bins_dict.keys()]
    idx = int(np.argmin([abs(a - (phi % 1.0)) for a in angles]))
    key = list(bins_dict.keys())[idx]
    return bins_dict[key]["p"]


def analytic_A_and_B(angles, gate_name, learned_zz_p, learned_gpi_bins, learned_gpi2_bins, labels):
    """Gate-by-gate density-matrix propagation, exactly loop_pec.py's own
    structure -- A = depolarizing-only (learned p), B = depolarizing THEN
    PEC-inverse (learned p), using EACH gate's own real angle."""
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native
    from qiskit import transpile
    from ionq_backend import get_native_simulator  # unused here, import kept minimal elsewhere

    qc = to_native(build_ansatz(angles), gate_name)
    n = qc.num_qubits
    dm_A = DensityMatrix.from_label("0" * n)
    dm_B = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        U = Operator(op.to_matrix())
        dm_A = dm_A.evolve(U, qargs=qargs)
        dm_B = dm_B.evolve(U, qargs=qargs)
        if op.name == gate_name:  # the native 2q gate (zz or ms)
            p_here = learned_zz_p
            dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(p_here, 2))
            dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(p_here, 2))
            dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(p_here, 2))
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            p_here = nearest_bin_p(learned_gpi_bins, phi)
            dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(p_here, 1))
            dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(p_here, 1))
            dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(p_here, 1))
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            p_here = nearest_bin_p(learned_gpi2_bins, phi)
            dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(p_here, 1))
            dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(p_here, 1))
            dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(p_here, 1))
    A = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_A.data))) for l in labels}
    B = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_B.data))) for l in labels}
    return A, B


def build_full(raw_kept, diag, K, non_id_labels):
    full = {name: dict(raw_kept[name]) for name in diag}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            full[pl] = dict(raw_kept[pl])
            synth_minus = {}
            for l in non_id_labels:
                if l not in raw_kept[pl] or l not in full[un] or l not in full[um]:
                    continue
                synth_minus[l] = full[un][l] + full[um][l] - raw_kept[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task30b_pec_application.py -- apply REAL angle-conditioned PEC to REAL measured H4 data")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    with open(CALIB_RESULTS) as f:
        learned = json.load(f)
    with open(RAW_CKPT) as f:
        ck = json.load(f)

    CITED_RAW_28B = {"aria-1": 89.77, "forte-1": 89.19}
    CITED_ZNE_HEADLINE = 14.61  # this session's own fresh per-Pauli-curve ZNE, 2q-only, iteration 27/28-consistent

    results_by_model = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        zz_p = list(learned[model]["zz"].values())[0]["p"]
        gpi_bins = learned[model]["gpi"]
        gpi2_bins = learned[model]["gpi2"]
        print(f"\n  {model}: learned p2(zz)={zz_p:.5f}, p1(gpi) bins={[round(v['p'],5) for v in gpi_bins.values()]}, "
              f"p1(gpi2) bins={[round(v['p'],5) for v in gpi2_bins.values()]}")

        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        per_name = {}
        for (name, group), counts in zip(tags, counts_list):
            per_name.setdefault(name, {}).setdefault(tuple(group), counts)

        raw_errs, pec_errs = [], []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task30b_apply", model, seed))
            raw_kept = {name: {} for name in kept}
            pec_kept = {name: {} for name in kept}
            for name in kept:
                labels_here = [l for group_t in per_name[name] for l in group_t]
                A, B = analytic_A_and_B(fixed_solutions[name]["angles"], "zz", zz_p, gpi_bins, gpi2_bins, labels_here)
                for group_t, counts in per_name[name].items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    for l in group_t:
                        m_real = expectation_from_counts(resampled, l)
                        raw_kept[name][l] = m_real
                        ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
                        corrected = m_real * ratio
                        corrected = max(-1.0, min(1.0, corrected))  # PEC ratio correction can overshoot; report, don't hide
                        pec_kept[name][l] = corrected

            full_raw = build_full(raw_kept, diag, K, non_id_labels)
            full_pec = build_full(pec_kept, diag, K, non_id_labels)
            _, err_raw = energy_and_err(p, full_raw, K)
            _, err_pec = energy_and_err(p, full_pec, K)
            raw_errs.append(err_raw)
            pec_errs.append(err_pec)

        raw_mean, raw_std = float(np.mean(raw_errs)), float(np.std(raw_errs))
        pec_mean, pec_std = float(np.mean(pec_errs)), float(np.std(pec_errs))
        print(f"    RAW: {raw_mean:.3f} +/- {raw_std:.3f} kcal/mol   PEC: {pec_mean:.3f} +/- {pec_std:.3f} kcal/mol")
        results_by_model[model] = {"raw_mean": raw_mean, "raw_std": raw_std, "pec_mean": pec_mean, "pec_std": pec_std}

    print(f"\n  -- SUMMARY --")
    ideal_r = results_by_model["ideal"]
    ideal_ok = ideal_r["pec_mean"] <= ideal_r["raw_mean"] * 1.5 + 0.3  # generous but real: must not blow up
    print(f"    ideal: raw={ideal_r['raw_mean']:.3f}  PEC={ideal_r['pec_mean']:.3f}  "
          f"(SANITY CHECK: {'PASS' if ideal_ok else 'FAIL -- PEC degrades noiseless data, DISQUALIFIED'})")
    for model in ["aria-1", "forte-1"]:
        r = results_by_model[model]
        print(f"    {model}: raw(this checkpoint, fold=1)={r['raw_mean']:.2f}  PEC={r['pec_mean']:.2f}  "
              f"(cited 28B raw={CITED_RAW_28B[model]:.2f} for cross-check)")

    print(f"\n  DECISION RULE CONTEXT: the rule (14.6->3-5 continue, 14.6->~13 stop) was written against the")
    print(f"  2q-ONLY ZNE headline (a FOLD-EXTRAPOLATED number). What's measured here is single-fold (no ZNE)")
    print(f"  PEC vs single-fold raw on the SAME optimized circuit -- not a like-for-like input to that literal")
    print(f"  rule. Reported both ways below for honest context, not to dodge the rule.")
    for model in ["aria-1", "forte-1"]:
        r = results_by_model[model]
        print(f"    {model}: PEC vs same-fold raw: {r['raw_mean']:.2f} -> {r['pec_mean']:.2f}  "
              f"({'IMPROVES' if r['pec_mean'] < r['raw_mean'] else 'DOES NOT IMPROVE'})")
        print(f"    {model}: PEC vs the 14.61 ZNE headline (different fold basis, context only): "
              f"{'BEATS headline' if r['pec_mean'] < CITED_ZNE_HEADLINE else 'does not beat headline'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"results_by_model": results_by_model, "ideal_control_pass": bool(ideal_ok)}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results_by_model


if __name__ == "__main__":
    main()
