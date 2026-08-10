#!/usr/bin/env python3
"""
task4_ablation_study.py — iteration 24, Task 4. THE FULL ABLATION STUDY:
one table, same data, same seeds, every mitigation combination this
project has ever built, so it's visible which pieces actually contribute
and which are decoration.
============================================================================
ROWS (per the explicit spec):
  raw                                          -- targets_d1.0.json, no ancilla, no SDP
  raw + leakage postselection                  -- spin_leakage_targets.json, ancilla postselected
  IonQ debiasing                               -- N/A (QPU-only, never run -- iteration 21 confirmed
                                                   real-QPU-only and the user declined spending real
                                                   QPU credits on it; NOT fabricated here)
  IonQ debiasing + leakage                     -- N/A (compounds on the above)
  PSD reconstruction + leakage                 -- NEW combination: Phase 1's SDP reconstruction
                                                   applied to the ancilla-postselected counts (never
                                                   tried together before -- Phase 1 only ever ran on
                                                   the non-ancilla checkpoint)
  PSD reconstruction + leakage + debiasing     -- N/A (compounds)
  PSD + leakage + debiasing + residual ZNE     -- N/A (compounds; ALSO: ZNE has shown NO plateau in
                                                   every test this project has run -- iterations 11,
                                                   13, 14, 19, 22, Phase 4 -- so even if debiasing were
                                                   available, no ZNE number would be added on top; this
                                                   row reports N/A + "no plateau", never a fabricated
                                                   number)
  EXTRA: Phase 2 selective hybrid + leakage    -- Phase 2's dominant-term hybrid (90% cutoff) applied
                                                   to the ancilla-postselected counts (also a new
                                                   combination)

All real rows use the SAME bootstrap seed convention (stable_seed,
N_SEEDS=8, SHOTS=10_000) as Phase 1/2/3 for direct comparability, and the
SAME resampling draw per seed is reused across every row so differences
reflect the METHOD, not different noise draws.

Run:
    python vqe/task4_ablation_study.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL, slot_names
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed, expectation_from_counts
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot, build_group_index, K, SHOTS, N_SEEDS
from phase2_dominant_term_mitigation import rank_terms, head_tail_split
from spin_leakage_postselect_ionq import postselect_counts

STD_CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "targets_d1.0.json")
LEAK_CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "spin_leakage_targets.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task4_ablation_study_results.json")
HYBRID_CUTOFF = 0.90  # matches Phase 2's own headline cutoff


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def load_leakage_checkpoint_blocks(ck):
    """spin_leakage_targets.json stores counts as a FLAT list of 468
    (=36 slots x 13 groups) 5-bit (4 register + 1 ancilla) count dicts,
    with idx_map giving the slot name per entry, appearing in contiguous
    13-entry blocks in group order (verified via an ideal-model energy
    sanity check before trusting this for real data, see module-level
    smoke test in development -- ideal err ~1.85 kcal/mol, matching the
    original spin_leakage_postselect_ionq_results.json's own 1.83-1.82
    kcal/mol to within bootstrap noise)."""
    idx_map = ck["idx_map"]
    out = {}
    i = 0
    while i < len(idx_map):
        name = idx_map[i]
        out[name] = list(range(i, i + 13))  # store flat-list INDICES, model-specific counts looked up later
        i += 13
    return out


def main():
    print("\n" + "=" * 96)
    print("  task4_ablation_study.py -- the full ablation table, same data, same seeds")
    print("=" * 96)

    with open(STD_CKPT_PATH) as f:
        std_ck = json.load(f)
    with open(LEAK_CKPT_PATH) as f:
        leak_ck = json.load(f)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=std_ck["d"], K=K)
    alpha_labels = std_ck["alpha_labels"]
    identity_label = std_ck["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    names = std_ck["target_names"]
    assert set(names) == set(slot_names(K))

    std_groups = std_ck["groups"]
    std_group_idx = build_group_index(std_groups)
    leak_groups = leak_ck["groups"]
    leak_group_idx = build_group_index(leak_groups)
    assert leak_groups == std_groups, "leakage checkpoint's group structure differs from the standard one -- stop"

    leak_blocks = load_leakage_checkpoint_blocks(leak_ck)
    assert set(leak_blocks.keys()) == set(names)

    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(alpha_labels, U)

    contributions, mismatch = rank_terms(p, P_S)
    assert mismatch < 1e-6
    head_terms, head_alpha_labels = head_tail_split(contributions, HYBRID_CUTOFF)
    print(f"  Phase 2 hybrid cutoff={HYBRID_CUTOFF}: {len(head_alpha_labels)}/{len(non_id_labels)} head labels")

    rows = {
        "raw": {}, "raw_plus_leakage": {}, "psd_plus_leakage": {}, "phase2hybrid_plus_leakage": {},
    }
    for model in ["ideal", "aria-1", "forte-1"]:
        std_counts = std_ck["counts"][model]
        leak_counts_flat = leak_ck["counts"][model]

        errs_raw, errs_leak_raw, errs_leak_psd, errs_leak_hybrid = [], [], [], []
        t0 = time.time()
        for seed in range(N_SEEDS):
            rng_std = np.random.default_rng(stable_seed("task4_std", model, seed))
            rng_leak = np.random.default_rng(stable_seed("task4_leak", model, seed))

            # -- row 1: raw, standard (non-ancilla) checkpoint --
            raw_std = {name: {} for name in names}
            for name in names:
                resampled = [bootstrap_counts(std_counts[name][gi], SHOTS, rng_std) for gi in range(len(std_groups))]
                for l in non_id_labels:
                    raw_std[name][l] = expectation_from_counts(resampled[std_group_idx[l]], l)
            _, err_raw = energy_and_err(p, raw_std, K)
            errs_raw.append(err_raw["err_vs_exact_kcal"])

            # -- resample the leakage (ancilla) checkpoint ONCE per seed, reused for all 3 leakage-based rows --
            raw_leak, m_dict_by_name, w_dict_by_name = {name: {} for name in names}, {}, {}
            for name in names:
                block_indices = leak_blocks[name]
                resampled5 = [bootstrap_counts(leak_counts_flat[block_indices[gi]], SHOTS, rng_leak)
                              for gi in range(len(leak_groups))]
                m_dict, w_dict = {}, {}
                for l in non_id_labels:
                    c5 = resampled5[leak_group_idx[l]]
                    c4 = postselect_counts(c5)
                    total = sum(c4.values())
                    if total == 0:
                        m = 0.0  # degenerate: bootstrap draw happened to keep zero ancilla=0 shots for this group
                    else:
                        m = expectation_from_counts(c4, l)
                    raw_leak[name][l] = m
                    m_dict[l] = m
                    w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(total, 1)
                m_dict_by_name[name] = m_dict
                w_dict_by_name[name] = w_dict

            # -- row 2: raw + leakage postselection, no SDP --
            _, err_leak_raw = energy_and_err(p, raw_leak, K)
            errs_leak_raw.append(err_leak_raw["err_vs_exact_kcal"])

            # -- row 5: PSD reconstruction + leakage (NEW combination) --
            phys_leak = {name: {} for name in names}
            for name in names:
                rho_slot = reconstruct_rho_slot(P_S, m_dict_by_name[name], w_dict_by_name[name], K)
                for l in non_id_labels:
                    phys_leak[name][l] = float(np.real(np.trace(rho_slot @ P_S[l])))
            _, err_leak_psd = energy_and_err(p, phys_leak, K)
            errs_leak_psd.append(err_leak_psd["err_vs_exact_kcal"])

            # -- extra row: Phase 2 selective hybrid + leakage (NEW combination) --
            hybrid_leak = {name: {} for name in names}
            for name in names:
                for l in non_id_labels:
                    hybrid_leak[name][l] = phys_leak[name][l] if l in head_alpha_labels else raw_leak[name][l]
            _, err_leak_hybrid = energy_and_err(p, hybrid_leak, K)
            errs_leak_hybrid.append(err_leak_hybrid["err_vs_exact_kcal"])

        rows["raw"][model] = {"mean": float(np.mean(errs_raw)), "std": float(np.std(errs_raw))}
        rows["raw_plus_leakage"][model] = {"mean": float(np.mean(errs_leak_raw)), "std": float(np.std(errs_leak_raw))}
        rows["psd_plus_leakage"][model] = {"mean": float(np.mean(errs_leak_psd)), "std": float(np.std(errs_leak_psd))}
        rows["phase2hybrid_plus_leakage"][model] = {"mean": float(np.mean(errs_leak_hybrid)), "std": float(np.std(errs_leak_hybrid))}
        print(f"  {model} ({time.time()-t0:.1f}s): raw={rows['raw'][model]['mean']:.2f}+/-{rows['raw'][model]['std']:.2f}  "
              f"raw+leak={rows['raw_plus_leakage'][model]['mean']:.2f}+/-{rows['raw_plus_leakage'][model]['std']:.2f}  "
              f"psd+leak={rows['psd_plus_leakage'][model]['mean']:.2f}+/-{rows['psd_plus_leakage'][model]['std']:.2f}  "
              f"phase2hybrid+leak={rows['phase2hybrid_plus_leakage'][model]['mean']:.2f}+/-{rows['phase2hybrid_plus_leakage'][model]['std']:.2f}")

    # -- MANDATORY ideal-data sanity check for the NEW psd+leakage combination --
    ideal = rows["psd_plus_leakage"]["ideal"]["mean"]
    ideal_raw = rows["raw_plus_leakage"]["ideal"]["mean"]
    disqualified = ideal > 3 * max(ideal_raw, 1e-6)
    print(f"\n  MANDATORY IDEAL-DATA CHECK (psd+leakage, new combination): ideal={ideal:.3f} vs "
          f"raw+leak ideal={ideal_raw:.3f} kcal/mol -- {'*** DISQUALIFIED ***' if disqualified else 'passes'}")

    print(f"\n  -- THE ABLATION TABLE (real hardware models, mean +/- std, kcal/mol) --")
    print(f"    {'row':<38} {'aria-1':>16} {'forte-1':>16}")
    for label, key in [("raw", "raw"), ("raw + leakage postselection", "raw_plus_leakage"),
                        ("IonQ debiasing", None), ("IonQ debiasing + leakage", None),
                        ("PSD reconstruction + leakage", "psd_plus_leakage"),
                        ("PSD + leakage + debiasing", None),
                        ("PSD + leakage + debiasing + residual ZNE", None),
                        ("EXTRA: Phase2 hybrid + leakage", "phase2hybrid_plus_leakage")]:
        if key is None:
            print(f"    {label:<38} {'N/A':>16} {'N/A':>16}")
        else:
            a = rows[key]["aria-1"]
            f_ = rows[key]["forte-1"]
            print(f"    {label:<38} {a['mean']:>7.2f}+/-{a['std']:<6.2f} {f_['mean']:>7.2f}+/-{f_['std']:<6.2f}")

    out = {
        "K": K, "hybrid_cutoff": HYBRID_CUTOFF, "n_head_labels": len(head_alpha_labels),
        "rows": rows, "disqualified_psd_plus_leakage": bool(disqualified),
        "na_rows": {
            "ionq_debiasing": "N/A -- QPU-only per iteration 21's own research; never run, not fabricated",
            "ionq_debiasing_plus_leakage": "N/A -- compounds on ionq_debiasing",
            "psd_leakage_debiasing": "N/A -- compounds on ionq_debiasing",
            "psd_leakage_debiasing_zne": "N/A -- compounds; ALSO ZNE has shown NO plateau in every prior test "
                                          "(iterations 11,13,14,19,22, Phase4) -- no number would be added even if debiasing existed",
        },
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return out


if __name__ == "__main__":
    main()
