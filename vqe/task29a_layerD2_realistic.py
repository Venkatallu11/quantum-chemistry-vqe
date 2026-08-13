#!/usr/bin/env python3
"""
task29a_layerD2_realistic.py -- iteration 29, Task A, Layer D2. The first
diagnostic (task29a_ideal_control_diagnosis.py) tested Layer D by
extrapolating the AGGREGATE total-energy curve, and it passed cleanly.
But that is NOT what task28d_all_gate_zne.py actually does: production
extrapolates ~756 INDIVIDUAL Pauli-expectation curves separately (one
fit + held-out selection + exclusion check PER (slot, label) pair), then
reassembles the alpha matrices from the 756 independently-extrapolated
values. That is a fundamentally different granularity, and Layer D's
exact-input test above cannot rule it out as the failure point.

This script reproduces task28d_all_gate_zne.py's REAL analysis code
(curve construction, two-stage held-out selection, exclusion, off-
diagonal synthesis) verbatim, but feeds it REALISTIC shot noise generated
by multinomial-sampling task29a's own exact statevector probabilities
(100,000 shots x 8 seeds, matching production exactly) instead of a real
IonQ submission -- statistically identical to a real run, zero IonQ
calls, zero cost, isolates the extrapolator's behavior under real shot
noise from every other explanation already ruled out.

Run:
    python vqe/task29a_layerD2_realistic.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
from task27d_held_out_zne import MODEL_CLASSES, FOLD_STAGE1_FIT, FOLD_STAGE1_HOLDOUT, FOLD_STAGE2_FIT, FOLD_STAGE2_HOLDOUT
import ef_fragment as effrag_mod
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import stable_seed, expectation_from_counts
from task28d_all_gate_zne import optimized_native_circuit, fold_all_gates
from qiskit.quantum_info import Statevector

K = 6
FOLD_FACTORS = [1, 3, 5, 7, 9]
GATE_NAME = "ms"
SHOTS = 100_000
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task29a_layerD2_realistic_results.json")


def exact_probs(qc):
    sv = Statevector.from_instruction(qc)
    return sv.probabilities_dict()


def main():
    print("\n" + "=" * 96)
    print("  task29a_layerD2_realistic.py -- realistic shot noise, PRODUCTION per-Pauli-curve ZNE")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)
    base = {name: optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME) for name in kept}

    print(f"  building exact probability tables (no shots yet, no IonQ) for all folds/slots/groups...")
    exact_probs_cache = {}
    for fold in FOLD_FACTORS:
        exact_probs_cache[fold] = {}
        for name in kept:
            folded = fold_all_gates(base[name], fold, GATE_NAME)
            exact_probs_cache[fold][name] = {}
            for group in groups:
                combined = effrag_mod.combined_basis_label(group)
                basis_qc = native_basis_change(combined, GATE_NAME)
                qc = folded.compose(basis_qc)
                exact_probs_cache[fold][name][tuple(group)] = exact_probs(qc)
    print("  done.")

    print(f"\n  sampling REALISTIC shot noise: {SHOTS} shots x {N_SEEDS} seeds per (fold, slot, group), "
          f"multinomial from the exact probabilities above (statistically identical to a real ionq_simulator run)")
    curves = {}
    for fold in FOLD_FACTORS:
        for name in kept:
            for group in groups:
                gt = tuple(group)
                probs = exact_probs_cache[fold][name][gt]
                bitstrings = list(probs.keys())
                parr = np.array([probs[b] for b in bitstrings])
                parr = parr / parr.sum()
                for seed in range(N_SEEDS):
                    rng = np.random.default_rng(stable_seed("task29a_d2", fold, name, gt, seed))
                    draws = rng.multinomial(SHOTS, parr)
                    counts = {b: int(n) for b, n in zip(bitstrings, draws) if n > 0}
                    for l in group:
                        curves.setdefault((name, l), {}).setdefault(fold, []).append(expectation_from_counts(counts, l))
    for key2 in curves:
        for fold in curves[key2]:
            curves[key2][fold] = float(np.mean(curves[key2][fold]))
    print(f"  built {len(curves)} (slot, label) curves, {N_SEEDS} seeds averaged per fold each")

    # -- EXACT reproduction of task28d_all_gate_zne.py's two-stage held-out ZNE (verbatim logic) --
    print("\n  -- two-stage held-out selection (verbatim production logic) --")
    stage1_errors, stage2_selected, stage2_errors = {}, {}, {}
    for key2, c in curves.items():
        if not all(f in c for f in FOLD_STAGE1_FIT + [FOLD_STAGE1_HOLDOUT]):
            continue
        vals_fit = [c[f] for f in FOLD_STAGE1_FIT]
        best_name, best_err = None, float("inf")
        for name, fitter in MODEL_CLASSES.items():
            fn = fitter(FOLD_STAGE1_FIT, vals_fit)
            if fn is None:
                continue
            try:
                pred = float(np.atleast_1d(fn([FOLD_STAGE1_HOLDOUT]))[0])
                if np.isfinite(pred):
                    err = abs(pred - c[FOLD_STAGE1_HOLDOUT])
                    if err < best_err:
                        best_err, best_name = err, name
            except Exception:
                continue
        if best_name is not None:
            stage1_errors[key2] = best_err
    mean_s1 = float(np.mean(list(stage1_errors.values()))) if stage1_errors else float("inf")

    for key2, c in curves.items():
        if not all(f in c for f in FOLD_STAGE2_FIT + [FOLD_STAGE2_HOLDOUT]):
            continue
        vals_fit = [c[f] for f in FOLD_STAGE2_FIT]
        best_name, best_err = None, float("inf")
        for name, fitter in MODEL_CLASSES.items():
            fn = fitter(FOLD_STAGE2_FIT, vals_fit)
            if fn is None:
                continue
            try:
                pred = float(np.atleast_1d(fn([FOLD_STAGE2_HOLDOUT]))[0])
                if np.isfinite(pred):
                    err = abs(pred - c[FOLD_STAGE2_HOLDOUT])
                    if err < best_err:
                        best_err, best_name = err, name
            except Exception:
                continue
        if best_name is not None:
            stage2_selected[key2] = best_name
            stage2_errors[key2] = best_err
    mean_s2 = float(np.mean(list(stage2_errors.values()))) if stage2_errors else float("inf")
    class_counts = {}
    for cls in stage2_selected.values():
        class_counts[cls] = class_counts.get(cls, 0) + 1
    print(f"  stage1 mean holdout err={mean_s1:.4f}  stage2 mean holdout err={mean_s2:.4f}")
    print(f"  stage2 model classes selected: {class_counts}")

    n_excluded, n_total, extrap = 0, 0, {}
    excluded_detail = []
    for key2, cls in stage2_selected.items():
        n_total += 1
        c = curves[key2]
        all_folds = sorted(c.keys())
        vals_all = [c[f] for f in all_folds]
        fn = MODEL_CLASSES[cls](all_folds, vals_all)
        if fn is None:
            n_excluded += 1
            continue
        try:
            val0 = float(np.atleast_1d(fn([0]))[0])
        except Exception:
            n_excluded += 1
            continue
        if not np.isfinite(val0) or abs(val0) > 1.0:
            n_excluded += 1
            excluded_detail.append((key2, cls, val0))
            continue
        extrap[key2] = val0
    print(f"  excluded (unphysical |val0|>1): {n_excluded}/{n_total}")
    if excluded_detail:
        worst = sorted(excluded_detail, key=lambda t: abs(t[2]) if np.isfinite(t[2]) else 1e18, reverse=True)[:8]
        print("  worst excluded extrapolations (slot, label, class, val0):")
        for key2, cls, val0 in worst:
            print(f"    {key2}  class={cls}  val0={val0}")

    full = {name: {} for name in kept}
    for name in kept:
        for l in non_id_labels:
            key2 = (name, l)
            if key2 in extrap:
                full[name][l] = extrap[key2]
            elif key2 in curves and 1 in curves[key2]:
                full[name][l] = curves[key2][1]
    full_complete = {n2: dict(full[n2]) for n2 in diag}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            full_complete[pl] = dict(full[pl])
            synth_minus = {}
            for l in non_id_labels:
                if l not in full[pl] or l not in full_complete[un] or l not in full_complete[um]:
                    continue
                cross = full[pl][l] - (full_complete[un][l] + full_complete[um][l]) / 2
                synth_minus[l] = full[pl][l] - 2 * cross
            full_complete[f"(u{n}-u{m})"] = synth_minus
    alpha_mats = combine_matrices(full_complete, p["alpha_labels"], p["identity_label"], K)
    E_zne, errs_zne = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                  exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    err_zne_kcal = errs_zne["err_vs_exact_kcal"]

    raw1 = {name: {} for name in kept}
    for name in kept:
        for l in non_id_labels:
            key2 = (name, l)
            if key2 in curves and 1 in curves[key2]:
                raw1[name][l] = curves[key2][1]
    raw1_complete = {n2: dict(raw1[n2]) for n2 in diag}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            raw1_complete[pl] = dict(raw1[pl])
            synth_minus = {}
            for l in non_id_labels:
                if l not in raw1[pl] or l not in raw1_complete[un] or l not in raw1_complete[um]:
                    continue
                cross = raw1[pl][l] - (raw1_complete[un][l] + raw1_complete[um][l]) / 2
                synth_minus[l] = raw1[pl][l] - 2 * cross
            raw1_complete[f"(u{n}-u{m})"] = synth_minus
    _, errs_raw1 = energy_from_alpha_matrices(
        combine_matrices(raw1_complete, p["alpha_labels"], p["identity_label"], K),
        p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
        exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    err_raw1_kcal = errs_raw1["err_vs_exact_kcal"]

    print(f"\n  IDEAL model, realistic shot noise ({SHOTS} shots x {N_SEEDS} seeds):")
    print(f"    raw(fold=1)      = {err_raw1_kcal:.4f} kcal/mol")
    print(f"    ALL-GATE-ZNE(0)  = {err_zne_kcal:.4f} kcal/mol")
    print(f"    ratio            = {err_zne_kcal/err_raw1_kcal:.2f}x")
    reproduces = err_zne_kcal > 5 * err_raw1_kcal
    print(f"    reproduces the 21x-style blowup seen in the real Task 28D run: {reproduces}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "n_curves": len(curves), "stage1_mean_holdout_err": mean_s1, "stage2_mean_holdout_err": mean_s2,
            "stage2_class_counts": class_counts, "n_excluded": n_excluded, "n_total": n_total,
            "err_raw_fold1_kcal": err_raw1_kcal, "err_all_gate_zne_fold0_kcal": err_zne_kcal,
            "reproduces_blowup": reproduces,
        }, f, indent=2, default=str)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
