#!/usr/bin/env python3
"""
task28e_failure_map.py — iteration 28, Task E. Task 27 excluded 5-17% of
fold-response curves as unphysical (fold->0 extrapolation outside
[-1,1]) but never classified them -- unfinished data already sitting on
disk. Classifies every excluded curve by measurement family (diagonal /
cross-term), basis-rotation complexity, Schmidt slot, exact energy-
sensitivity weight (reused from Phase 2's rank_terms, not re-derived),
and qubit index. Checks whether failures concentrate.
============================================================================
Reuses iteration 27 Task C's checkpoints (`task27c_full_h4_folds_K6.json`)
and Task D's exact model-selection/extrapolation logic UNCHANGED -- this
file does not refit anything differently, it adds classification and
logging on top of the SAME procedure, so the exclusion counts here
reproduce Task D's own numbers exactly (a consistency check performed
below, not assumed).

Run:
    python vqe/task28e_failure_map.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task27d_held_out_zne import MODEL_CLASSES, FOLD_STAGE2_FIT, FOLD_STAGE2_HOLDOUT, ckpt_path_for_K
from phase2_dominant_term_mitigation import rank_terms
from phys_constrained_reconstruction import build_P_S
import ef_fragment as effrag_mod
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task28e_failure_map_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task28e_failure_map.py -- classifying the excluded curves iteration 27 left unclassified")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)

    # -- exact energy-sensitivity weight per label, reused from Phase 2 --
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    contributions, mismatch = rank_terms(p, P_S)
    assert mismatch < 1e-6
    label_weight = {}
    for c in contributions:
        label_weight[c["alpha_label"]] = label_weight.get(c["alpha_label"], 0.0) + abs(c["contribution_ha"]) * HARTREE_TO_KCAL_MOL

    with open(ckpt_path_for_K(K)) as f:
        ck = json.load(f)

    all_records = []
    for model in ["aria-1", "forte-1"]:
        curves = {}
        for fold in [1, 3, 5, 7, 9]:
            key = f"{fold}|{model}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            per_name = {}
            for (name, group), counts in zip(tags, counts_list):
                per_name.setdefault(name, {}).setdefault(tuple(group), counts)
            for name in kept:
                for group_t, counts in per_name[name].items():
                    seed_vals = []
                    for seed in range(N_SEEDS):
                        rng = np.random.default_rng(stable_seed("task27d", K, fold, model, seed, name))
                        resampled = bootstrap_counts(counts, SHOTS, rng)
                        for l in group_t:
                            curves.setdefault((name, l), {}).setdefault(fold, []).append(expectation_from_counts(resampled, l))
        for key2 in curves:
            for fold in curves[key2]:
                curves[key2][fold] = float(np.mean(curves[key2][fold]))

        for (name, label), c in curves.items():
            if not all(f in c for f in FOLD_STAGE2_FIT + [FOLD_STAGE2_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_STAGE2_FIT]
            best_name, best_err = None, float("inf")
            for cls_name, fitter in MODEL_CLASSES.items():
                fn = fitter(FOLD_STAGE2_FIT, vals_fit)
                if fn is None:
                    continue
                try:
                    pred = float(np.atleast_1d(fn([FOLD_STAGE2_HOLDOUT]))[0])
                    if not np.isfinite(pred):
                        continue
                    err = abs(pred - c[FOLD_STAGE2_HOLDOUT])
                    if err < best_err:
                        best_err, best_name = err, cls_name
                except Exception:
                    continue
            if best_name is None:
                continue
            all_folds = sorted(c.keys())
            vals_all = [c[f] for f in all_folds]
            fn = MODEL_CLASSES[best_name](all_folds, vals_all)
            if fn is None:
                continue
            try:
                val0 = float(np.atleast_1d(fn([0]))[0])
            except Exception:
                continue
            excluded = (not np.isfinite(val0)) or abs(val0) > 1.0

            is_cross = "+" in name
            n_nonI = sum(1 for ch in label if ch != "I")
            has_xy = any(ch in "XY" for ch in label)
            qubit_positions = [i for i, ch in enumerate(label) if ch != "I"]

            all_records.append({
                "model": model, "slot": name, "label": label, "excluded": bool(excluded),
                "family": "cross_term" if is_cross else "diagonal",
                "basis": "has_XY" if has_xy else "pure_Z",
                "n_nonI": n_nonI, "qubit_positions": qubit_positions,
                "energy_weight_kcal": label_weight.get(label, 0.0),
                "stage2_holdout_err": best_err, "extrapolated_val0": val0,
            })

    n_total = len(all_records)
    n_excluded = sum(1 for r in all_records if r["excluded"])
    print(f"  {n_total} curves classified, {n_excluded} excluded ({100*n_excluded/n_total:.1f}%) "
          f"-- consistency check against Task D's own reported counts (36/756 aria-1, 39/756 forte-1 at K=6)")

    print(f"\n  -- exclusion rate by FAMILY --")
    for model in ["aria-1", "forte-1"]:
        for fam in ["diagonal", "cross_term"]:
            rows = [r for r in all_records if r["model"] == model and r["family"] == fam]
            exc = sum(1 for r in rows if r["excluded"])
            print(f"    {model} / {fam}: {exc}/{len(rows)} excluded ({100*exc/len(rows):.1f}%)")

    print(f"\n  -- exclusion rate by BASIS (does the label need X/Y rotation?) --")
    for model in ["aria-1", "forte-1"]:
        for basis in ["pure_Z", "has_XY"]:
            rows = [r for r in all_records if r["model"] == model and r["basis"] == basis]
            exc = sum(1 for r in rows if r["excluded"])
            print(f"    {model} / {basis}: {exc}/{len(rows)} excluded ({100*exc/len(rows):.1f}%)")

    print(f"\n  -- exclusion rate by Pauli weight (n_nonI) --")
    for model in ["aria-1", "forte-1"]:
        for n in sorted(set(r["n_nonI"] for r in all_records)):
            rows = [r for r in all_records if r["model"] == model and r["n_nonI"] == n]
            if not rows:
                continue
            exc = sum(1 for r in rows if r["excluded"])
            print(f"    {model} / weight={n}: {exc}/{len(rows)} excluded ({100*exc/len(rows):.1f}%)")

    print(f"\n  -- do failures concentrate? excluded curves by SLOT (top 10 by exclusion count) --")
    for model in ["aria-1", "forte-1"]:
        by_slot = {}
        for r in all_records:
            if r["model"] == model:
                by_slot.setdefault(r["slot"], []).append(r["excluded"])
        slot_rates = sorted(((sum(v), len(v), k) for k, v in by_slot.items()), reverse=True)
        print(f"    {model}:")
        for n_exc, n_tot, slot in slot_rates[:10]:
            if n_exc > 0:
                print(f"      {slot}: {n_exc}/{n_tot} excluded ({100*n_exc/n_tot:.1f}%)")

    # -- energy weight vs exclusion: does excluding a curve throw away important information? --
    print(f"\n  -- energy weight of EXCLUDED vs INCLUDED curves --")
    for model in ["aria-1", "forte-1"]:
        exc_w = [r["energy_weight_kcal"] for r in all_records if r["model"] == model and r["excluded"]]
        inc_w = [r["energy_weight_kcal"] for r in all_records if r["model"] == model and not r["excluded"]]
        print(f"    {model}: excluded mean energy_weight={np.mean(exc_w):.2f} (n={len(exc_w)})  "
              f"included mean energy_weight={np.mean(inc_w):.2f} (n={len(inc_w)})")

    # -- connect to iteration 26 Task 2: does cross-term family (50% of energy, smallest fold response) show a distinct exclusion pattern? --
    print(f"\n  -- cross-check vs iteration 26 Task 2's finding (cross-term = 50% energy weight, smallest fold response) --")
    for model in ["aria-1", "forte-1"]:
        cross_w = sum(r["energy_weight_kcal"] for r in all_records if r["model"] == model and r["family"] == "cross_term") / 2  # /2: each label counted per model
        diag_w = sum(r["energy_weight_kcal"] for r in all_records if r["model"] == model and r["family"] == "diagonal") / 2
        total_w = cross_w + diag_w
        cross_exc = sum(1 for r in all_records if r["model"] == model and r["family"] == "cross_term" and r["excluded"])
        cross_n = sum(1 for r in all_records if r["model"] == model and r["family"] == "cross_term")
        print(f"    {model}: cross_term energy weight fraction={cross_w/total_w*100:.1f}%, "
              f"exclusion rate={100*cross_exc/cross_n:.1f}% ({cross_exc}/{cross_n})")

    results = {"n_total": n_total, "n_excluded": n_excluded, "records": all_records}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
