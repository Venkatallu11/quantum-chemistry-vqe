#!/usr/bin/env python3
"""
task29e_combined_optimizer_reconstruction.py -- iteration 29, Task E. THE
COMBINATION NOBODY HAS RUN. Task 28B (TrappedIonOptimizerPlugin: raw
132.49 -> 89.19 kcal/mol) and Task 28G (PSD/trace-1 constrained
reconstruction: 13.86 -> 12.78 kcal/mol) were measured on SEPARATE
circuits and never combined.
============================================================================
This is the genuinely free/immediate form of that combination: Task 28G's
own SDP reconstruction (`build_P_S` / `reconstruct_rho_slot`, unchanged)
applied directly to Task 28B's already-collected fold=1 optimized-circuit
counts (`task28b_optimized_raw.json`) -- no new circuits, no new real
submission, a pure downstream reprocessing exactly like 28G's own
precedent. There is no existing 2q-only FOLD SWEEP on the optimized
circuit (28D's fold data is ALL-GATE folding, a different and currently
invalid pipeline per Task A above) -- so this is a single-point (fold=1)
reconstruction, not a fold-ZNE combination. It answers the more basic
question first: does physical-constraint reconstruction help an
ALREADY gate-count-reduced circuit's raw measurement, on top of the gate
reduction itself?

MANDATORY IDEAL-DATA SANITY CHECK: run identically on the "ideal" model's
own counts from the same checkpoint.

Run:
    python vqe/task29e_combined_optimizer_reconstruction.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS = 8
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "task28b_optimized_raw.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task29e_combined_optimizer_reconstruction_results.json")


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
    print("  task29e_combined_optimizer_reconstruction.py -- 28G's SDP reconstruction on 28B's optimized circuit")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    print(f"  K={K}: {len(kept)} kept slots, P_S projections built for {len(p['alpha_labels'])} labels")

    with open(CKPT_PATH) as f:
        ck = json.load(f)

    results_by_model = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        per_name = {}
        for (name, group), counts in zip(tags, counts_list):
            per_name.setdefault(name, {}).setdefault(tuple(group), counts)

        raw_errs, phys_errs = [], []
        n_sdp_fail = 0
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task29e", model, seed))
            raw_kept = {name: {} for name in kept}
            phys_kept = {name: {} for name in kept}
            for name in kept:
                m_dict, w_dict = {}, {}
                for group_t, counts in per_name[name].items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    total = sum(resampled.values())
                    for l in group_t:
                        m = expectation_from_counts(resampled, l)
                        m_dict[l] = m
                        raw_kept[name][l] = m
                        var = max(1 - m ** 2, 1e-4) / max(total, 1)
                        w_dict[l] = 1.0 / var
                try:
                    rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                    for l in non_id_labels:
                        phys_kept[name][l] = float(np.real(np.trace(rho_slot @ P_S[l])))
                except RuntimeError:
                    n_sdp_fail += 1
                    phys_kept[name] = dict(m_dict)

            full_raw = build_full(raw_kept, diag, K, non_id_labels)
            full_phys = build_full(phys_kept, diag, K, non_id_labels)
            _, err_raw = energy_and_err(p, full_raw, K)
            _, err_phys = energy_and_err(p, full_phys, K)
            raw_errs.append(err_raw)
            phys_errs.append(err_phys)

        n_total_slots = N_SEEDS * len(kept)
        raw_mean, raw_std = float(np.mean(raw_errs)), float(np.std(raw_errs))
        phys_mean, phys_std = float(np.mean(phys_errs)), float(np.std(phys_errs))
        print(f"\n  {model} ({n_sdp_fail}/{n_total_slots} SDP solve failures):")
        print(f"    RAW  (28B optimized circuit, no reconstruction): {raw_mean:.2f} +/- {raw_std:.2f} kcal/mol")
        print(f"    PHYS (28B optimized circuit, 28G reconstruction): {phys_mean:.2f} +/- {phys_std:.2f} kcal/mol")
        better = "PHYS BETTER" if phys_mean < raw_mean else ("RAW BETTER" if raw_mean < phys_mean else "TIE")
        print(f"    -- {better}")
        results_by_model[model] = {
            "raw_mean_kcal": raw_mean, "raw_std_kcal": raw_std,
            "phys_mean_kcal": phys_mean, "phys_std_kcal": phys_std,
            "n_sdp_fail": n_sdp_fail, "n_total_slots": n_total_slots,
        }

    print(f"\n  -- SUMMARY: does 28G's reconstruction help ON TOP OF 28B's gate-count reduction? --")
    for model in ["ideal", "aria-1", "forte-1"]:
        r = results_by_model[model]
        note = "  (SANITY CHECK: must not make the ideal control worse)" if model == "ideal" else ""
        print(f"    {model}: raw={r['raw_mean_kcal']:.2f}  phys={r['phys_mean_kcal']:.2f}{note}")

    cited_28g_raw = {"aria-1": None, "forte-1": None}   # 28G's ORIGINAL (un-optimized) circuit fold=1 raw, for context
    cited_28b_raw = {"aria-1": 89.77, "forte-1": 89.19}  # 28B's own reported raw numbers (should match this file's raw_mean)
    print(f"\n  cross-check vs 28B's own reported raw numbers (should match this file's RAW row):")
    for model in ["aria-1", "forte-1"]:
        print(f"    {model}: 28B reported={cited_28b_raw[model]:.2f}  this file's raw={results_by_model[model]['raw_mean_kcal']:.2f}")

    results = {"results_by_model": results_by_model, "cited_28b_raw": cited_28b_raw}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
