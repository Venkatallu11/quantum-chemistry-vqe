#!/usr/bin/env python3
"""
task28g_constrained_reconstruction.py — iteration 28, Task G. For each
fold, reconstruct rho_lambda subject to rho >= 0, Tr(rho) = 1, and
confinement to the weight-2 Schmidt subspace (Phase 1's `build_P_S` /
`reconstruct_rho_slot`, iteration 23, extended here to the FOLD dimension
for the first time), THEN E(lambda) = Tr(H rho). ZNE runs on those
physically-valid per-fold ENERGIES, not on fragile per-Pauli curves.
============================================================================
REUSES, UNCHANGED: `phys_constrained_reconstruction.py`'s `build_P_S` and
`reconstruct_rho_slot` (Phase 1's exact PSD/trace-1 SDP, K=6, one 6x6
convex problem per slot) -- this file adds the FOLD dimension on top,
nothing about the SDP itself is modified. Also reuses iteration 27 Task
C's own checkpoint data (`task27c_full_h4_folds_K6.json`) -- NO new
circuits, NO new real submission, this is a downstream reprocessing of
already-collected counts, exactly like Task 1's constrained-reconstruction
precedent was itself a free/no-new-data intervention.

CONSTRAINTS REDUCE VARIANCE, NEVER INJECT INFORMATION: the SDP is fit
per (fold, model, slot) independently from that fold's own resampled
counts only -- no information is shared ACROSS folds, and no exact/ideal
value is ever used as a constraint or prior. This mirrors the explicit
warning in the task text ("the earlier aggressive channel inversion was
correctly rejected for exactly that reason").

MANDATORY IDEAL-DATA SANITY CHECK (Phase 3's own precedent, run here on
every new combination): the identical pipeline is run on the "ideal"
model's own fold data. If constrained reconstruction were injecting
spurious structure, it would show up as the ideal curve getting WORSE,
not better -- checked explicitly, not assumed clean.

HELD-OUT VALIDATION on the resulting per-fold ENERGY curve, matching
Task D/Task 27D's own procedure: fit on folds [1,3,5], predict fold=7;
refit on [1,3,5,7], predict fold=9; then refit the selected model class
on ALL points for the final fold->0 extrapolation.

Run:
    python vqe/task28g_constrained_reconstruction.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K, ckpt_path_for_K
from task27d_held_out_zne import MODEL_CLASSES, FOLD_STAGE1_FIT, FOLD_STAGE1_HOLDOUT, FOLD_STAGE2_FIT, FOLD_STAGE2_HOLDOUT
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS = 8
FOLD_FACTORS = [1, 3, 5, 7, 9]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task28g_constrained_reconstruction_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


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


def zne_on_curve(fold_to_val):
    """Held-out two-stage model selection + extrapolation, exactly Task
    D/27D's procedure, applied here to a per-fold SCALAR (energy error),
    not a per-Pauli expectation curve."""
    if not all(f in fold_to_val for f in FOLD_STAGE1_FIT + [FOLD_STAGE1_HOLDOUT]):
        return None
    vals_fit1 = [fold_to_val[f] for f in FOLD_STAGE1_FIT]
    best1_name, best1_err = None, float("inf")
    for name, fitter in MODEL_CLASSES.items():
        fn = fitter(FOLD_STAGE1_FIT, vals_fit1)
        if fn is None:
            continue
        try:
            pred = float(np.atleast_1d(fn([FOLD_STAGE1_HOLDOUT]))[0])
            if np.isfinite(pred):
                err = abs(pred - fold_to_val[FOLD_STAGE1_HOLDOUT])
                if err < best1_err:
                    best1_err, best1_name = err, name
        except Exception:
            continue

    vals_fit2 = [fold_to_val[f] for f in FOLD_STAGE2_FIT]
    best2_name, best2_err = None, float("inf")
    for name, fitter in MODEL_CLASSES.items():
        fn = fitter(FOLD_STAGE2_FIT, vals_fit2)
        if fn is None:
            continue
        try:
            pred = float(np.atleast_1d(fn([FOLD_STAGE2_HOLDOUT]))[0])
            if np.isfinite(pred):
                err = abs(pred - fold_to_val[FOLD_STAGE2_HOLDOUT])
                if err < best2_err:
                    best2_err, best2_name = err, name
        except Exception:
            continue

    if best2_name is None:
        return {"stage1_holdout_err": best1_err, "stage2_holdout_err": None,
                "selected_class": None, "extrapolated_val0": None, "excluded": True}

    all_folds = sorted(fold_to_val.keys())
    vals_all = [fold_to_val[f] for f in all_folds]
    fn = MODEL_CLASSES[best2_name](all_folds, vals_all)
    val0 = None
    excluded = True
    if fn is not None:
        try:
            val0 = float(np.atleast_1d(fn([0]))[0])
            excluded = not np.isfinite(val0)
        except Exception:
            val0 = None
    return {"stage1_holdout_err": best1_err, "stage2_holdout_err": best2_err,
            "selected_class": best2_name, "extrapolated_val0": val0, "excluded": excluded}


def main():
    print("\n" + "=" * 96)
    print("  task28g_constrained_reconstruction.py -- PSD/trace-1 reconstruction BEFORE fold extrapolation")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    print(f"  K={K}: {len(kept)} kept slots, P_S projections built for {len(p['alpha_labels'])} labels")

    with open(ckpt_path_for_K(K)) as f:
        ck = json.load(f)

    results_by_model = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        t0 = time.time()
        # -- per fold: bootstrap resample, reconstruct rho per slot, get RAW and PHYS per-fold energy error --
        raw_fold_errs = {f: [] for f in FOLD_FACTORS}
        phys_fold_errs = {f: [] for f in FOLD_FACTORS}
        n_sdp_fail = 0
        for fold in FOLD_FACTORS:
            key = f"{fold}|{model}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            per_name = {}
            for (name, group), counts in zip(tags, counts_list):
                per_name.setdefault(name, {}).setdefault(tuple(group), counts)

            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task28g", fold, model, seed))
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
                raw_fold_errs[fold].append(err_raw)
                phys_fold_errs[fold].append(err_phys)

        raw_fold_mean = {f: float(np.mean(v)) for f, v in raw_fold_errs.items()}
        phys_fold_mean = {f: float(np.mean(v)) for f, v in phys_fold_errs.items()}
        n_total_slots = N_SEEDS * len(FOLD_FACTORS) * len(kept)
        print(f"\n  {model} ({time.time()-t0:.1f}s, {n_sdp_fail}/{n_total_slots} SDP solve failures):")
        print(f"    RAW per-fold error:  " + "  ".join(f"f{f}={raw_fold_mean[f]:.2f}" for f in FOLD_FACTORS))
        print(f"    PHYS per-fold error: " + "  ".join(f"f{f}={phys_fold_mean[f]:.2f}" for f in FOLD_FACTORS))

        zne_raw = zne_on_curve(raw_fold_mean)
        zne_phys = zne_on_curve(phys_fold_mean)
        print(f"    ZNE(raw fold-energy curve):  fold0={zne_raw['extrapolated_val0']}  "
              f"class={zne_raw['selected_class']}  excluded={zne_raw['excluded']}")
        print(f"    ZNE(phys fold-energy curve): fold0={zne_phys['extrapolated_val0']}  "
              f"class={zne_phys['selected_class']}  excluded={zne_phys['excluded']}")

        results_by_model[model] = {
            "raw_fold_mean": raw_fold_mean, "phys_fold_mean": phys_fold_mean,
            "n_sdp_fail": n_sdp_fail, "n_total_slots": n_total_slots,
            "zne_raw": zne_raw, "zne_phys": zne_phys,
        }

    print(f"\n  -- SUMMARY: does PHYS-CONSTRAINED fold-ZNE beat RAW fold-ZNE? --")
    for model in ["ideal", "aria-1", "forte-1"]:
        r = results_by_model[model]
        raw0 = r["zne_raw"]["extrapolated_val0"]
        phys0 = r["zne_phys"]["extrapolated_val0"]
        raw1 = r["raw_fold_mean"][1]
        note = ""
        if model == "ideal":
            note = "  (SANITY CHECK: physical reconstruction must not make the ideal control worse)"
        if raw0 is not None and phys0 is not None:
            better = "PHYS BETTER" if phys0 < raw0 else ("RAW BETTER" if raw0 < phys0 else "TIE")
            print(f"    {model}: raw(fold1)={raw1:.2f}  ZNE-raw(fold0)={raw0:.2f}  "
                  f"ZNE-phys(fold0)={phys0:.2f}  -- {better}{note}")
        else:
            print(f"    {model}: raw(fold1)={raw1:.2f}  ZNE-raw(fold0)={raw0}  ZNE-phys(fold0)={phys0} "
                  f"(one or both excluded as unphysical/non-finite){note}")

    results = {"results_by_model": results_by_model}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
