#!/usr/bin/env python3
"""
task84_frame_fit_convergence.py -- iteration 84.

WHY THIS EXISTS
task83 scanned a single smooth parameter alpha through the locked
pipeline's joint Schmidt-frame fit and found chi2/dof was NOT a smooth
function of it at the locked n_restarts=4: between alpha=0.90 and 1.10 it
read 7.88, 11.89, 5.93, 12.04, 7.03. A one-parameter, continuous change
to the input data cannot produce that. The same scan at n_restarts=16 was
smooth. That points at the frame fit landing in different local minima,
which matters because the pre-registered hardware number (chi2/dof 5.929,
shared-frame error 0.864 kcal/mol) was produced at n_restarts=4 and sits
at a conspicuously sharp dip in the 4-restart curve.

This script holds the data fixed at the locked alpha=1 and varies only the
optimiser: n_restarts, and the RNG seed. If 5.929 / 0.864 are converged
values they will be stable. If they are an artifact of one lucky restart
draw, they will move.

Nothing here touches the exact energy as a fit objective; the energy is
read out as a consequence, exactly as in the locked pipeline.

Run:
    PYTHONHASHSEED=0 python vqe/task84_frame_fit_convergence.py \
        --counts vqe/task73_hw_counts.json
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task73_hw_analysis as t73
import task71_hw_plan_rehearsal_v2 as t71
import task83_hw_correction_magnitude_refit as t83
from task73_hw_analysis import K, FRAME_SEED, CORRECTION_GPI2

RESTART_LADDER = [4, 8, 16, 32, 64]
SEEDS = [FRAME_SEED, FRAME_SEED + 1, FRAME_SEED + 2, FRAME_SEED + 3, FRAME_SEED + 4]


def fit_with(ctx, corr, n_restarts, seed):
    unit = {n: {l: 1.0 for l in ctx.non_id} for n in ctx.kept}
    U, _, _ = t71.fit_joint_frame(np.eye(K), ctx.P_S, K, ctx.kept, ctx.non_id, corr, unit,
                                  np.random.default_rng(seed), n_restarts=n_restarts)
    return U


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--design", choices=["full", "subset"], default="full")
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = t73.Context(args.design)
    with open(args.counts) as f:
        blob = json.load(f)
    counts = {k: dict(v) for k, v in blob["counts"].items()}
    t73.add_reused_task49(counts, ctx)
    post, kept_shots = t83.postselect(counts, ctx)
    rat = t83.ratios(ctx, CORRECTION_GPI2)
    ng = len(ctx.groups)
    corr = t83.apply_alpha(post, rat, {gi: 1.0 for gi in range(ng)})   # locked alpha=1

    print(f"\n  task84 -- frame-fit convergence at the LOCKED alpha=1, job {blob.get('job_id')}")
    print(f"  pre-registered values: chi2/dof 5.929, shared-frame err 0.864 kcal/mol")
    print(f"  (produced at n_restarts=4, seed={FRAME_SEED})\n")

    rows = []
    print(f"  {'restarts':>9} {'seed':>6} {'chi2/dof':>9} {'err_kcal':>9}")
    for nr in RESTART_LADDER:
        for sd in SEEDS:
            U = fit_with(ctx, corr, nr, sd)
            r, _ = t83.residuals(ctx, post, corr, kept_shots, U)
            c = t83.chi2_dof(r)
            _, e = t83.energies(ctx, corr, U)
            rows.append({"n_restarts": nr, "seed": int(sd), "chi2_dof": c, "err_kcal": e})
            mark = "   <-- pre-registered config" if (nr == 4 and sd == FRAME_SEED) else ""
            print(f"  {nr:9d} {sd:6d} {c:9.3f} {e:9.3f}{mark}")
        sub = [x for x in rows if x["n_restarts"] == nr]
        cc = np.array([x["chi2_dof"] for x in sub]); ee = np.array([x["err_kcal"] for x in sub])
        print(f"  {'':9} {'spread':>6} chi2 {cc.min():.3f}-{cc.max():.3f}  "
              f"err {ee.min():.3f}-{ee.max():.3f}\n")

    print("  == summary: best (lowest chi2) found at each restart count ==")
    print(f"  {'restarts':>9} {'best chi2/dof':>14} {'its err_kcal':>13}")
    for nr in RESTART_LADDER:
        sub = [x for x in rows if x["n_restarts"] == nr]
        best = min(sub, key=lambda d: d["chi2_dof"])
        print(f"  {nr:9d} {best['chi2_dof']:14.3f} {best['err_kcal']:13.3f}")

    allc = np.array([x["chi2_dof"] for x in rows])
    alle = np.array([x["err_kcal"] for x in rows])
    print(f"\n  across all {len(rows)} optimiser configurations, data identical:")
    print(f"    chi2/dof  {allc.min():.3f} - {allc.max():.3f}")
    print(f"    err_kcal  {alle.min():.3f} - {alle.max():.3f}")
    print(f"  => the energy varies by {alle.max()-alle.min():.3f} kcal/mol from the")
    print(f"     optimiser alone, with the hardware data held fixed.")

    out = os.path.join(os.path.dirname(__file__), "task84_frame_fit_convergence_results.json")
    with open(out, "w") as f:
        json.dump({"job_id": blob.get("job_id"), "locked_alpha": 1.0,
                   "preregistered": {"chi2_dof": 5.9286, "err_kcal": 0.8643,
                                     "n_restarts": 4, "seed": int(FRAME_SEED)},
                   "rows": rows,
                   "chi2_range": [float(allc.min()), float(allc.max())],
                   "err_range": [float(alle.min()), float(alle.max())]}, f, indent=2)
    print(f"\n  Saved -> {out}\n")


if __name__ == "__main__":
    main()
