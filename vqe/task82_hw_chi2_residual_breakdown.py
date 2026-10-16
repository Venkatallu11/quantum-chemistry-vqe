#!/usr/bin/env python3
"""
task82_hw_chi2_residual_breakdown.py -- iteration 82.

WHY THIS EXISTS
The real Forte Enterprise 1 run (job 01a138e7, 16 Oct 2026) reproduced the
RAW error almost exactly as the noise model predicted (hardware 23.84 vs
forte-1 model 24.06 kcal/mol, 0.99x) but the shared-frame fit's chi2/dof
jumped from ~2 on every noise model to 5.93 on hardware -- and stayed there
in all 20 bootstrap draws (6.41-7.22). A tripled chi2 with an accurate raw
prediction means the hardware carries structure the ZZ/GPi/GPi2 correction
does not describe. This script localises it.

It re-runs task73's LOCKED pipeline unchanged and re-derives the identical
residual vector, but tags every residual with its slot, GC group, Pauli
label and Pauli weight, then aggregates. No fitting choice is altered: the
point-estimate chi2/dof printed here must equal task73's 5.9286, and the
script asserts that.

Run:
    PYTHONHASHSEED=0 python vqe/task82_hw_chi2_residual_breakdown.py \
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
from task73_hw_analysis import key, K, FRAME_SEED, CORRECTION_GPI2
from task72_rehearsal_sign_and_chi2_check import FORCED_ZERO_RATIO
from task36_joint_schmidt_frame import slot_vector, N_THETA
from general_commuting_measurements import expectations_from_counts


def tagged_residuals(counts, ctx, gpi2):
    """task73.analyze()'s residual loop, with provenance kept for each term."""
    # label -> GC group index (each Pauli label is measured in exactly one group)
    label_group = {}
    for gi, labels in enumerate(ctx.groups):
        for l in labels:
            label_group.setdefault(l, gi)

    post, kept_shots = {n: {} for n in ctx.kept}, {n: {} for n in ctx.kept}
    for name in ctx.kept:
        for gi in ctx.measured[name]:
            filtered = {bs[1:]: c for bs, c in counts[key(name, gi)].items() if bs[0] == "0"}
            n_kept = sum(filtered.values())
            exp = expectations_from_counts(filtered, ctx.diagonalizers[gi]) if filtered else {}
            for l in ctx.groups[gi]:
                post[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                kept_shots[name][l] = n_kept

    corr = t71.corrected_observables_gc_aware(
        post, gpi2, ctx.kept, ctx.fixed_solutions, ctx.diag_natives,
        ctx.diagonalizers, ctx.groups, ctx.measured)

    unit = {n: {l: 1.0 for l in ctx.non_id} for n in ctx.kept}
    U_hat, _, _ = t71.fit_joint_frame(np.eye(K), ctx.P_S, K, ctx.kept, ctx.non_id,
                                      corr, unit, np.random.default_rng(FRAME_SEED),
                                      n_restarts=4)

    rows, n_forced = [], 0
    for name in ctx.kept:
        v = slot_vector(U_hat, name, K)
        for l, y in corr[name].items():
            m = post[name][l]
            ratio = y / m if abs(m) > 1e-9 else 1.0
            if abs(ratio) < FORCED_ZERO_RATIO:
                n_forced += 1
                continue
            sigma = np.sqrt(max(1.0 - m * m, 1e-4) / max(kept_shots[name][l], 1)) * abs(ratio)
            r = (float(np.real(v @ ctx.P_S[l] @ v)) - y) / max(sigma, 1e-6)
            rows.append({
                "slot": name,
                "group": label_group.get(l, -1),
                "label": l,
                "weight": sum(1 for ch in l if ch != "I"),
                "resid": float(r),
                "raw_exp": float(m),
                "corrected_exp": float(y),
                "correction_ratio": float(ratio),
                "sigma": float(sigma),
                "kept_shots": int(kept_shots[name][l]),
            })
    return rows, n_forced


def summarise(vals):
    a = np.asarray(vals, dtype=float)
    return {"n": int(a.size), "chi2_sum": float(a @ a),
            "mean_sq": float((a * a).mean()), "mean_resid": float(a.mean()),
            "rms": float(np.sqrt((a * a).mean())), "max_abs": float(np.abs(a).max())}


def group_by(rows, field):
    out = {}
    for r in rows:
        out.setdefault(r[field], []).append(r["resid"])
    return {str(k): summarise(v) for k, v in sorted(out.items(), key=lambda kv: str(kv[0]))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--design", choices=["full", "subset"], default="full")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = t73.Context(args.design)
    with open(args.counts) as f:
        blob = json.load(f)
    counts = {k: dict(v) for k, v in blob["counts"].items()}
    reused = t73.add_reused_task49(counts, ctx)   # mutates counts, returns added keys
    gpi2 = CORRECTION_GPI2

    rows, n_forced = tagged_residuals(counts, ctx, gpi2)
    r = np.asarray([x["resid"] for x in rows])
    chi2_dof = float(r @ r / (len(r) - N_THETA))

    print(f"\n  task82 -- residual breakdown, {blob.get('backend')} job {blob.get('job_id')}")
    print(f"  residual terms: {len(rows)}  forced-zero excluded: {n_forced}  reused: {reused}")
    # task73's own value for this job. Tolerance is 1e-3 absolute, not float
    # equality: the frame fit's optimizer path moves the last digits across
    # scipy/BLAS builds (observed 1.8e-4 relative between two machines on
    # this same job), while any real pipeline divergence moves it far more.
    print(f"  chi2/dof = {chi2_dof:.4f}   (task73 on this job: 5.9286)")
    assert abs(chi2_dof - 5.9286) < 1e-3, (
        f"pipeline drifted from task73: {chi2_dof} vs 5.9286")
    print("  matches task73 -- same residual vector, only tagged\n")

    by_group = group_by(rows, "group")
    by_weight = group_by(rows, "weight")
    by_slot = group_by(rows, "slot")
    total = float(r @ r)

    print("  == by GC group ==")
    print(f"  {'group':>6} {'n':>5} {'%chi2':>7} {'mean r^2':>9} {'rms':>7} {'mean r':>8} {'max|r|':>7}")
    for g, s in by_group.items():
        print(f"  {g:>6} {s['n']:5d} {100*s['chi2_sum']/total:6.1f}% {s['mean_sq']:9.2f} "
              f"{s['rms']:7.2f} {s['mean_resid']:8.2f} {s['max_abs']:7.1f}")

    print("\n  == by Pauli weight ==")
    print(f"  {'wt':>6} {'n':>5} {'%chi2':>7} {'mean r^2':>9} {'rms':>7} {'mean r':>8} {'max|r|':>7}")
    for w, s in by_weight.items():
        print(f"  {w:>6} {s['n']:5d} {100*s['chi2_sum']/total:6.1f}% {s['mean_sq']:9.2f} "
              f"{s['rms']:7.2f} {s['mean_resid']:8.2f} {s['max_abs']:7.1f}")

    print("\n  == 12 worst single residuals ==")
    print(f"  {'slot':>12} {'g':>2} {'label':>7} {'r':>8} {'raw':>8} {'corr':>8} {'ratio':>7}")
    for x in sorted(rows, key=lambda d: -abs(d["resid"]))[:12]:
        print(f"  {x['slot']:>12} {x['group']:2d} {x['label']:>7} {x['resid']:8.1f} "
              f"{x['raw_exp']:8.3f} {x['corrected_exp']:8.3f} {x['correction_ratio']:7.3f}")

    print("\n  == 6 worst slots by mean r^2 ==")
    for name, s in sorted(by_slot.items(), key=lambda kv: -kv[1]["mean_sq"])[:6]:
        print(f"  {name:>12}  n={s['n']:3d}  mean r^2={s['mean_sq']:8.2f}  rms={s['rms']:6.2f}")

    out = args.out or os.path.join(os.path.dirname(__file__),
                                   "task82_hw_chi2_residual_breakdown_results.json")
    with open(out, "w") as f:
        json.dump({"job_id": blob.get("job_id"), "backend": blob.get("backend"),
                   "chi2_dof": chi2_dof, "n_residuals": len(rows),
                   "n_forced_zero_excluded": n_forced, "chi2_total": total,
                   "by_group": by_group, "by_weight": by_weight, "by_slot": by_slot,
                   "rows": rows}, f, indent=2)
    print(f"\n  Saved -> {out}\n")


if __name__ == "__main__":
    main()
