#!/usr/bin/env python3
"""
task32a_robustness_outlier_diagnosis.py -- iteration 32, Task A. Follow-up
to iteration 31 Task H, which flagged a real, undiagnosed gap: 2 of 97
robustness-envelope draws blew up to 886.4 and 932.8 kcal/mol (dominating
Q99=888.29) while the other 95 stayed under 68 kcal/mol, and per-draw noise-
model parameters were NOT saved that run, so the cause could only be
guessed at ("likely driven by GPi2's wide uncertainty bound (Task A),
not chased to full resolution").

LOCAL / FREE SIMULATOR ONLY -- no real QPU submission, no hardware budget
touched, per this session's explicit instruction.

This script:
  1. Re-runs task31h's EXACT pipeline (same sample_noise_model,
     evaluate_one_model, same seed=31) with N=97 fixed to match the
     original run's scope, but additionally logs the FULL noise-model
     dict for every draw (the concrete gap Task H flagged).
  2. Reports Spearman correlation of each of the 7 noise-model parameters
     against |E-E_exact| across all N draws -- a correlational first pass.
  3. For whichever draws land in the extreme tail (>500 kcal/mol, matching
     Task H's own informal threshold), runs a CAUSAL ablation: re-evaluates
     that exact draw's noise model with ONE parameter at a time swapped
     for the population median, holding everything else fixed, to isolate
     which single parameter's value is actually responsible (correlation
     is not causation; this tests it directly).

Run:
    python vqe/task32a_robustness_outlier_diagnosis.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
from task31h_robustness_envelope import (
    sample_noise_model, evaluate_one_model, K, GATE_NAME,
)
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task29c_manifold_estimator import build_full_from_a  # noqa: F401 (re-exported by task31h, imported here for clarity only)
from phys_constrained_reconstruction import build_P_S

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32a_robustness_outlier_diagnosis_results.json")
PARAM_KEYS = [
    "p_zz_true", "p_gpi_true", "p_gpi2_true",
    "angle_bias_zz", "angle_bias_gpi", "angle_bias_gpi2",
    "readout_err", "calib_ratio_zz", "calib_ratio_1q",
]
OUTLIER_THRESHOLD_KCAL = 500.0  # matches Task H's own informal "2 draws >800" framing, set a bit below to be inclusive
N_FIXED = 97  # matches iteration 31 Task H's actual run exactly, disclosed rather than re-derived from a fresh timing measurement


def main():
    print("\n" + "=" * 96)
    print("  task32a_robustness_outlier_diagnosis.py -- iteration 32 Task A (LOCAL/FREE SIMULATOR ONLY)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag_slots, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    print(f"  N={N_FIXED} (fixed to match iteration 31 Task H's actual run, not re-derived from timing)")

    rng = np.random.default_rng(31)  # SAME seed as task31h's main-loop rng -- reproduces the same 97 draws
    models, errs_exact = [], []
    t_start = time.time()
    for i in range(N_FIXED):
        model = sample_noise_model(rng)
        err_e, _ = evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, model, rng)
        models.append(model)
        errs_exact.append(err_e)
        if (i + 1) % 10 == 0:
            print(f"    {i+1}/{N_FIXED} done, {time.time()-t_start:.1f}s elapsed, last err={err_e:.2f}")

    errs_exact = np.array(errs_exact)
    q = {q_: float(np.percentile(errs_exact, q_)) for q_ in [50, 90, 95, 99]}
    print(f"\n  -- reproduced quantiles (should match iteration 31 Task H: Q50=4.53 Q90=25.13 Q95=51.22 Q99=888.29) --")
    for q_ in [50, 90, 95, 99]:
        print(f"    Q{q_} = {q[q_]:.4f} kcal/mol")

    # -- step 2: correlation pass --
    print(f"\n  -- Spearman correlation: each noise-model parameter vs |E-E_exact| across all {N_FIXED} draws --")
    corr = {}
    for key in PARAM_KEYS:
        vals = np.array([m[key] for m in models])
        rho, pval = spearmanr(vals, errs_exact)
        corr[key] = {"rho": float(rho), "pval": float(pval)}
        print(f"    {key:16s}  rho={rho:+.3f}  p={pval:.4f}")

    # -- step 3: identify outliers, causal ablation --
    outlier_idx = [i for i, e in enumerate(errs_exact) if e > OUTLIER_THRESHOLD_KCAL]
    print(f"\n  -- {len(outlier_idx)} draw(s) exceed {OUTLIER_THRESHOLD_KCAL} kcal/mol: indices {outlier_idx} --")

    medians = {key: float(np.median([m[key] for m in models])) for key in PARAM_KEYS}
    ablation_results = {}
    for idx in outlier_idx:
        outlier_model = models[idx]
        print(f"\n  draw #{idx}: err={errs_exact[idx]:.2f} kcal/mol")
        print(f"    full model: {json.dumps(outlier_model, indent=6)}")
        per_param = {}
        for key in PARAM_KEYS:
            swapped = dict(outlier_model)
            swapped[key] = medians[key]
            # independent, freshly-seeded rng per ablation eval -- we are testing whether THIS
            # parameter's extreme value causes the blowup, not reproducing the original random stream
            abl_rng = np.random.default_rng(10_000 + idx * 100 + PARAM_KEYS.index(key))
            err_abl, _ = evaluate_one_model(p, fixed_solutions, kept, diag_slots, P_S, non_id_labels, swapped, abl_rng)
            per_param[key] = float(err_abl)
            flag = "  <-- COLLAPSES toward normal range" if err_abl < 100 else ""
            print(f"      swap {key:16s} (draw={outlier_model[key]:+.5f} -> median={medians[key]:+.5f})"
                  f"  ->  err={err_abl:8.2f} kcal/mol{flag}")
        ablation_results[str(idx)] = {
            "original_err": float(errs_exact[idx]),
            "model": outlier_model,
            "per_param_swap_err": per_param,
        }

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N": N_FIXED,
            "quantiles": q,
            "errs_exact": errs_exact.tolist(),
            "models": models,
            "correlation": corr,
            "outlier_threshold_kcal": OUTLIER_THRESHOLD_KCAL,
            "outlier_indices": outlier_idx,
            "medians": medians,
            "ablation": ablation_results,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
