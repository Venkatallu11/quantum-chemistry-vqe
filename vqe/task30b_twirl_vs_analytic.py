#!/usr/bin/env python3
"""
task30b_twirl_vs_analytic.py -- direct, per-label comparison of literal
quasi-probability circuit twirling (REAL, task30b_literal_twirling.py)
against the analytic-ratio shortcut (task30b_pec_application.py) on the
SAME 4 representative slots' REAL raw data (Task 28B's checkpoint).
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task30b_pec_application import analytic_A_and_B
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts

K = 6
SHOTS = 100_000
N_SEEDS = 8
REPRESENTATIVE_SLOTS = ["u_0", "u_1", "(u0+u1)", "(u2+u3)"]
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "task28b_optimized_raw.json")
CALIB_RESULTS = os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")
TWIRL_RESULTS = os.path.join(os.path.dirname(__file__), "task30b_literal_twirling_results.json")

p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)

with open(CALIB_RESULTS) as f:
    learned = json.load(f)
with open(RAW_CKPT) as f:
    ck = json.load(f)
with open(TWIRL_RESULTS) as f:
    twirl = json.load(f)

print("\n" + "=" * 96)
print("  literal twirling vs analytic-ratio shortcut -- same 4 slots, same real raw data")
print("=" * 96)

for model in ["aria-1", "forte-1"]:
    zz_p = list(learned[model]["zz"].values())[0]["p"]
    gpi_bins = learned[model]["gpi"]
    gpi2_bins = learned[model]["gpi2"]
    tags = ck["tags"][model]
    counts_list = ck["counts"][model]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)

    print(f"\n  {model}:")
    diffs = []
    for name in REPRESENTATIVE_SLOTS:
        labels_here = [l for group_t in per_name[name] for l in group_t]
        A, B = analytic_A_and_B(fixed_solutions[name]["angles"], "zz", zz_p, gpi_bins, gpi2_bins, labels_here)

        for seed in range(1):  # one seed-averaged raw value per label, matching production style
            rng = np.random.default_rng(stable_seed("task30b_cmp", model, name))
            raw_here = {}
            for group_t, counts in per_name[name].items():
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    raw_here[l] = expectation_from_counts(resampled, l)

        for l in labels_here[:3]:  # a few representative labels per slot to keep the table readable
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            analytic_pec = max(-1.0, min(1.0, raw_here[l] * ratio))
            key = f"{name}|{l}"
            literal_val = twirl["literal_pec"][model].get(key)
            if literal_val is None:
                continue
            literal_clamped = max(-1.0, min(1.0, literal_val))
            diff = abs(analytic_pec - literal_clamped)
            diffs.append(diff)
            print(f"    {name:<10} {l:<6}: raw={raw_here[l]:+.4f}  analytic_PEC={analytic_pec:+.4f}  "
                  f"literal_twirl_PEC={literal_clamped:+.4f}  |diff|={diff:.4f}")

    print(f"\n    {model} SUMMARY: mean|analytic - literal|={np.mean(diffs):.4f}, "
          f"median={np.median(diffs):.4f}, max={np.max(diffs):.4f}  (n={len(diffs)} labels)")
    print(f"    for reference, shot noise alone at N_MC=8 draws is expected to give per-label "
          f"scatter of similar order -- this is NOT a zero-tolerance check")

print("\n" + "=" * 96 + "\n")
