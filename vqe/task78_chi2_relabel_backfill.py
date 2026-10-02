#!/usr/bin/env python3
"""
task78_chi2_relabel_backfill.py -- iteration 78. Corrects the "chi2_dof"
field in the committed results of the fresh-simulator-submission tasks
(task72_h4_full_experiment trials 1-4, task73_h4_pooled_experiment,
task77_h4_full84 trials 1-2).

Those files stored fit_joint_frame's unit-weight cost / dof under the name
"chi2_dof" (values ~0.0002-0.002). With unit weights that is a mean squared
residual, not a chi2. task72_h4_full_experiment.analyze_backend now reports it
as "unit_weight_msr" and computes the real shot-noise-weighted "chi2_dof"
(same definition as task72_rehearsal_sign_and_chi2_check / task73_hw_analysis).

For every committed entry this script re-runs analyze_backend on the committed
checkpoint and only rewrites the file if BOTH
  - shared_frame_err_kcal reproduces the committed value (to 1e-6), and
  - the new unit_weight_msr equals the committed "chi2_dof" (to 1e-9),
i.e. it is the same fit. The energy fields are left exactly as committed;
only the chi2 fields are relabelled/added.

Run:
    PYTHONHASHSEED=0 python vqe/task78_chi2_relabel_backfill.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
import task72_h4_full_experiment as t72
import task73_h4_pooled_experiment as t73p
import task77_h4_full84_real_trials as t77
from task59_h4_gc_no_frame_fit import BACKENDS

HERE = os.path.dirname(__file__)
SUMMARY_PATH = os.path.join(HERE, "task78_chi2_relabel_backfill_results.json")
CHI2_KEYS = ("unit_weight_msr", "chi2_dof", "n_chi2_labels", "n_forced_zero_excluded")


def relabel(committed, fresh, where):
    assert abs(fresh["shared_frame_err_kcal"] - committed["shared_frame_err_kcal"]) < 1e-6, \
        f"{where}: shared-frame error does not reproduce -- STOP"
    assert abs(fresh["unit_weight_msr"] - committed["chi2_dof"]) < 1e-9, \
        f"{where}: frame fit differs from the committed one -- STOP"
    out = dict(committed)
    for k in CHI2_KEYS:
        out[k] = fresh[k]
    return out


def load(path):
    with open(path) as f:
        return json.load(f)


def save(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)


def main():
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0")
    summary = {}
    ctx54 = t72.setup()
    ctx84 = t77.setup_full84()

    jobs = [(f"task72 trial {t}", ctx54, t72.ckpt_path(t), t72.results_path(t)) for t in range(1, 5)]
    jobs += [(f"task77 trial {t}", ctx84, t77.ckpt_path(t), t77.results_path(t)) for t in (1, 2)]
    for label, ctx, ckpt, res_path in jobs:
        state, committed = load(ckpt), load(res_path)
        new = {bn: relabel(committed[bn], t72.analyze_backend(ctx, bn, state), f"{label} {bn}") for bn in BACKENDS}
        save(res_path, new)
        summary[label] = {bn: {k: new[bn][k] for k in CHI2_KEYS} for bn in BACKENDS}
        print(f"  {label}: " + "  ".join(f"{bn} chi2/dof={new[bn]['chi2_dof']:.2f}" for bn in BACKENDS))

    pooled_path = os.path.join(HERE, "task73_h4_pooled_experiment_results.json")
    pooled = load(pooled_path)
    state = t73p.pool_states([load(t72.ckpt_path(t)) for t in range(1, t73p.N_TRIALS + 1)])
    for bn in BACKENDS:
        pooled["pooled"][bn] = relabel(pooled["pooled"][bn], t72.analyze_backend(ctx54, bn, state), f"pooled {bn}")
    save(pooled_path, pooled)
    summary["task73 pooled (4,400 shots)"] = {bn: {k: pooled["pooled"][bn][k] for k in CHI2_KEYS} for bn in BACKENDS}
    print("  task73 pooled: " + "  ".join(f"{bn} chi2/dof={pooled['pooled'][bn]['chi2_dof']:.2f}" for bn in BACKENDS))

    save(SUMMARY_PATH, summary)
    print(f"  Saved -> {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
