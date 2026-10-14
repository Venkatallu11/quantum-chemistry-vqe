#!/usr/bin/env python3
"""
task79_h4_full84_pooled.py -- iteration 79. Pools the 3 fresh free-simulator
trials of the full 84-circuit design (task77, ionq_simulator -- NOT hardware)
into 3 x 1,100 = 3,300 shots/circuit and reruns the same locked analysis
(task72 analyze_backend) once.

Question: is the no-frame estimator's >1 kcal/mol error on forte-1 at
1,100 shots (1.26 / 1.14 / 2.46) shot noise, or a fixed bias like the
54-circuit subset's (task73 pooled 4,400 shots: forte-1 no_frame still 2.45)?

Result (kcal/mol; no_frame / shared_frame / chi2_dof):
    ideal    0.452 / 0.023 / 0.97
    aria-1   0.860 / 0.026 / 2.24
    forte-1  0.813 / 0.072 / 2.71
On the full design the no-frame error does fall under 1 with more shots, but
it levels off near 0.8 while shared-frame goes to ~0.05, so a residual no-frame
bias remains. A single 1,100-shot hardware run should expect no_frame ~1-2.5
on forte-1-like noise. No new data, no submission; per-trial re-analysis is
asserted to reproduce the committed task77 results first.

Run:
    PYTHONHASHSEED=0 python vqe/task79_h4_full84_pooled.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
import task77_h4_full84_real_trials as t77
import task72_h4_full_experiment as t72
from task73_h4_pooled_experiment import pool_states
from task59_h4_gc_no_frame_fit import BACKENDS
KEYS = ("raw_err_kcal", "no_frame_err_kcal", "shared_frame_err_kcal")

ctx = t77.setup_full84()
states = [json.load(open(t77.ckpt_path(t))) for t in (1, 2, 3)]

# 1) sanity: each trial re-analysed must reproduce its committed results file
for t, st in enumerate(states, 1):
    stored = json.load(open(t77.results_path(t)))
    for bn in BACKENDS:
        e = t72.analyze_backend(ctx, bn, st)
        for k in KEYS:
            assert abs(e[k] - stored[bn][k]) < 1e-6, f"trial {t} {bn} {k}: {e[k]} vs {stored[bn][k]}"
print("per-trial re-analysis reproduces committed results: PASS")

# 2) pooled (3 x 1100 = 3300 shots/circuit), analysed twice for reproducibility
pooled = pool_states(states)
tot = {sum(c.values()) for v in pooled["done"].values() for c in v["counts"]}
print("pooled shots per circuit:", tot)
out = {}
for bn in BACKENDS:
    a, b = t72.analyze_backend(ctx, bn, pooled), t72.analyze_backend(ctx, bn, pooled)
    assert all(abs(a[k] - b[k]) < 1e-9 for k in KEYS)
    out[bn] = {k: a[k] for k in KEYS + ("chi2_dof",)}
    print(f"POOLED {bn:8s} raw={a['raw_err_kcal']:+.3f} no_frame={a['no_frame_err_kcal']:+.3f} "
          f"shared={a['shared_frame_err_kcal']:+.3f} chi2/dof={a['chi2_dof']:.2f}")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "task79_h4_full84_pooled_results.json"), "w") as f:
    json.dump(out, f, indent=2)
