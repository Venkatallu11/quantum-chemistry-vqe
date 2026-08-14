#!/usr/bin/env python3
"""
task30c_ablation_matrix.py -- iteration 30, Task C. THE ABLATION MATRIX,
direct corrections only (ZNE frozen this iteration). Same real data, same
seeds, ideal/aria-1/forte-1, on Task 28B's optimizer-reduced native
circuit (285->69 1q gates, 132.49->89.19 raw): raw / +leakage /
+PEC / +PEC+manifold / +PEC+manifold+leakage / +PEC+leakage+PSD(general).
Every row gets the mandatory ideal-control check.
============================================================================
LEAKAGE, SCOPED: full 13-group leakage post-selection needs an ancilla-
augmented circuit (a NEW real submission, `spin_leakage_postselect_ionq.py`'s
own design) -- out of reach of this iteration's remaining scope. What IS
free (reuses Task 28B's already-collected data, iteration 18 Part 1's own
precedent): the ONE measurement group with no basis-rotation gates
(`IIZZ`/`ZIIZ`/`ZZII`) has bitstrings that directly reflect the physical
Z-basis electron occupation -- post-selecting THOSE on Hamming weight==2
is real leakage removal for 3 of 37 labels, disclosed as partial.

PEC: the analytic-ratio shortcut validated (mostly) in Task B, using the
REAL calibrated p2(zz)=0.014, p1(gpi)~0.0001-0.0005, p1(gpi2)=GPi-mean
fallback.

MANIFOLD: the 5-parameter pure-state fit, cleared in Task A.

The LAST row uses general PSD/trace-1 reconstruction (28G's SDP, allows
mixed states) INSTEAD of the pure-state manifold -- not stacked with it,
since a pure state is automatically PSD (rank-1), making "+manifold+PSD"
a no-op; the meaningful comparison is manifold (strong, pure-state
constraint) vs general-PSD (weaker, mixed-state-allowing constraint) as
ALTERNATIVE final reconstruction steps after the same upstream leakage+PEC
corrections -- disclosed explicitly, not silently reinterpreted.

Run:
    python vqe/task30c_ablation_matrix.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import analytic_A_and_B
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state as fit_pure_state_29c
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS = 8
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "task28b_optimized_raw.json")
CALIB_RESULTS = os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")
Z_BASIS_GROUP = ("IIZZ", "ZIIZ", "ZZII")  # the one group with no basis rotation -- direct leakage signal
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task30c_ablation_matrix_results.json")


def hamming_weight(bitstring):
    return bitstring.count("1")


def leakage_postselect(counts):
    kept = {bs: c for bs, c in counts.items() if hamming_weight(bs) == 2}
    return kept if kept else dict(counts)  # never fully empty-out a slot -- fall back to raw if postselection kills everything


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
    print("  task30c_ablation_matrix.py -- direct corrections only, ZNE frozen")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(CALIB_RESULTS) as f:
        learned = json.load(f)
    with open(RAW_CKPT) as f:
        ck = json.load(f)

    rows = ["raw", "raw+leakage", "raw+PEC", "raw+PEC+manifold", "raw+PEC+manifold+leakage", "raw+PEC+leakage+generalPSD"]
    results_by_model = {}

    for model in ["ideal", "aria-1", "forte-1"]:
        zz_p = list(learned[model]["zz"].values())[0]["p"]
        gpi_bins = learned[model]["gpi"]
        gpi2_bins = learned[model]["gpi2"]

        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        per_name = {}
        for (name, group), counts in zip(tags, counts_list):
            per_name.setdefault(name, {}).setdefault(tuple(group), counts)

        row_errs = {row: [] for row in rows}
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task30c", model, seed))

            raw_kept = {name: {} for name in kept}
            leak_kept = {name: {} for name in kept}
            for name in kept:
                for group_t, counts in per_name[name].items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    for l in group_t:
                        raw_kept[name][l] = expectation_from_counts(resampled, l)
                    if group_t == Z_BASIS_GROUP:
                        resampled_leak = leakage_postselect(resampled)
                        for l in group_t:
                            leak_kept[name][l] = expectation_from_counts(resampled_leak, l)
                    else:
                        for l in group_t:
                            leak_kept[name][l] = raw_kept[name][l]

            pec_kept = {name: {} for name in kept}
            pec_leak_kept = {name: {} for name in kept}
            for name in kept:
                labels_here = list(raw_kept[name].keys())
                A, B = analytic_A_and_B(fixed_solutions[name]["angles"], "zz", zz_p, gpi_bins, gpi2_bins, labels_here)
                for l in labels_here:
                    ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
                    pec_kept[name][l] = max(-1.0, min(1.0, raw_kept[name][l] * ratio))
                    pec_leak_kept[name][l] = max(-1.0, min(1.0, leak_kept[name][l] * ratio))

            manifold_kept = {}
            for name in kept:
                v0 = target_coeff_vector(name, K)
                a_hat, _ = fit_pure_state_29c(P_S, pec_kept[name], {l: 1.0 for l in pec_kept[name]}, K, v0,
                                               seed=stable_seed("task30c_mfd", model, name, seed))
                manifold_kept[name] = {l: float(np.real(a_hat @ P_S[l] @ a_hat)) for l in pec_kept[name]}

            manifold_leak_kept = {}
            for name in kept:
                v0 = target_coeff_vector(name, K)
                a_hat, _ = fit_pure_state_29c(P_S, pec_leak_kept[name], {l: 1.0 for l in pec_leak_kept[name]}, K, v0,
                                               seed=stable_seed("task30c_mfdleak", model, name, seed))
                manifold_leak_kept[name] = {l: float(np.real(a_hat @ P_S[l] @ a_hat)) for l in pec_leak_kept[name]}

            psd_leak_kept = {}
            for name in kept:
                w_dict = {l: 1.0 for l in pec_leak_kept[name]}
                try:
                    rho = reconstruct_rho_slot(P_S, pec_leak_kept[name], w_dict, K)
                    psd_leak_kept[name] = {l: float(np.real(np.trace(rho @ P_S[l]))) for l in pec_leak_kept[name]}
                except RuntimeError:
                    psd_leak_kept[name] = dict(pec_leak_kept[name])

            for row, data in [
                ("raw", raw_kept), ("raw+leakage", leak_kept), ("raw+PEC", pec_kept),
                ("raw+PEC+manifold", manifold_kept), ("raw+PEC+manifold+leakage", manifold_leak_kept),
                ("raw+PEC+leakage+generalPSD", psd_leak_kept),
            ]:
                full = build_full(data, diag, K, non_id_labels)
                _, err = energy_and_err(p, full, K)
                row_errs[row].append(err)

        results_by_model[model] = {row: {"mean": float(np.mean(v)), "std": float(np.std(v))} for row, v in row_errs.items()}
        print(f"\n  {model}:")
        for row in rows:
            r = results_by_model[model][row]
            print(f"    {row:<28}: {r['mean']:.3f} +/- {r['std']:.3f} kcal/mol")

    print(f"\n" + "=" * 96)
    print(f"  IDEAL-CONTROL CHECK (every row must not exceed raw's own error by much)")
    ideal_raw = results_by_model["ideal"]["raw"]["mean"]
    for row in rows:
        v = results_by_model["ideal"][row]["mean"]
        ok = v <= ideal_raw * 2.0 + 0.5
        print(f"    {row:<28}: {v:.3f} kcal/mol vs raw={ideal_raw:.3f}  {'PASS' if ok else 'FAIL -- DISQUALIFIED'}")

    print(f"\n  -- DOES COMPOUNDING HELP? (aria-1 / forte-1, kcal/mol) --")
    for model in ["aria-1", "forte-1"]:
        r = results_by_model[model]
        raw_v = r["raw"]["mean"]
        leak_v = r["raw+leakage"]["mean"]
        pec_v = r["raw+PEC"]["mean"]
        pec_mfd_v = r["raw+PEC+manifold"]["mean"]
        pec_mfd_leak_v = r["raw+PEC+manifold+leakage"]["mean"]
        pec_leak_psd_v = r["raw+PEC+leakage+generalPSD"]["mean"]
        leak_gain = raw_v - leak_v
        pec_gain = raw_v - pec_v
        mfd_gain_on_top_of_pec = pec_v - pec_mfd_v
        leak_gain_on_top_of_pec_mfd = pec_mfd_v - pec_mfd_leak_v
        naive_sum = leak_gain + pec_gain
        actual_pec_mfd_leak_gain = raw_v - pec_mfd_leak_v
        print(f"\n    {model}:")
        print(f"      leakage alone: {raw_v:.2f} -> {leak_v:.2f} (gain {leak_gain:.2f})")
        print(f"      PEC alone:     {raw_v:.2f} -> {pec_v:.2f} (gain {pec_gain:.2f})")
        print(f"      manifold ON TOP OF PEC: {pec_v:.2f} -> {pec_mfd_v:.2f} (additional gain {mfd_gain_on_top_of_pec:.2f})")
        print(f"      leakage ON TOP OF PEC+manifold: {pec_mfd_v:.2f} -> {pec_mfd_leak_v:.2f} "
              f"(additional gain {leak_gain_on_top_of_pec_mfd:.2f})")
        print(f"      naive sum of independent gains (leakage+PEC): {naive_sum:.2f} vs actual "
              f"PEC+manifold+leakage combined gain: {actual_pec_mfd_leak_gain:.2f} "
              f"-- {'gains roughly ADD' if abs(naive_sum-actual_pec_mfd_leak_gain) < 0.3*naive_sum else 'gains do NOT simply add (correlated/compounding effects)'}")
        print(f"      general-PSD (leakage+PEC, no manifold): {pec_leak_psd_v:.2f}  vs  "
              f"manifold (leakage+PEC+manifold): {pec_mfd_leak_v:.2f}  "
              f"-- {'manifold wins' if pec_mfd_leak_v < pec_leak_psd_v else 'general PSD wins'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results_by_model, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results_by_model


if __name__ == "__main__":
    main()
