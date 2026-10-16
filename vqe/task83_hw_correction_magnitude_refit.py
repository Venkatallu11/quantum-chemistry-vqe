#!/usr/bin/env python3
"""
task83_hw_correction_magnitude_refit.py -- iteration 83.

HYPOTHESIS UNDER TEST (from task82)
The real Forte Enterprise 1 run (job 01a138e7) showed shared-frame
chi2/dof = 5.93 against ~2 on every noise model. task82 localised two
candidate causes:

  Fault 1  the residual grows with the SIZE of the correction applied.
           Labels at correction ratio 1.00-1.05 give mean r^2 = 3.69;
           labels at 1.15-1.25 give 8.51 and carry 48% of total chi2.
  Fault 2  GC group 0 carries a systematic SIGNED offset (mean residual
           -1.32 at weight 1, -1.73 at weight 2), not scatter.

If Fault 1 is real, shrinking the correction's magnitude should drop
chi2 sharply. The locked pipeline applies

    corrected = raw * ratio,        ratio = B[label] / A[label]

so this script replaces that with

    corrected = raw * (1 + alpha * (ratio - 1))

alpha = 1 is the locked pipeline. alpha = 0 is no correction at all.
Model B gives each GC group its own alpha, which also tests Fault 2's
group-locality.

BLINDNESS
alpha is chosen by minimising the shared-frame fit's chi2/dof ONLY. The
exact FCI energy is never an input to any fit here -- same discipline as
every correction-parameter selection in this project. The energies are
printed as a CONSEQUENCE of the blind fit, never as its objective.

HELD-OUT VALIDATION
Fitting alpha on the same residuals whose chi2 we then report would drop
chi2 by construction. So alpha is fit on half the slots and chi2 is
evaluated on the other half (and vice versa). A real miscalibration
improves the held-out half; overfitting does not.

Run:
    PYTHONHASHSEED=0 python vqe/task83_hw_correction_magnitude_refit.py \
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

ALPHA_GRID = np.round(np.arange(0.0, 1.41, 0.05), 3)   # default; see --alpha-max


def postselect(counts, ctx):
    """Identical to the locked pipeline's postselection step."""
    post, kept_shots = {n: {} for n in ctx.kept}, {n: {} for n in ctx.kept}
    for name in ctx.kept:
        for gi in ctx.measured[name]:
            filtered = {bs[1:]: c for bs, c in counts[key(name, gi)].items() if bs[0] == "0"}
            n_kept = sum(filtered.values())
            exp = expectations_from_counts(filtered, ctx.diagonalizers[gi]) if filtered else {}
            for l in ctx.groups[gi]:
                post[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                kept_shots[name][l] = n_kept
    return post, kept_shots


def ratios(ctx, gpi2):
    """The locked pipeline's per-(slot, label) correction ratio B/A."""
    out = {name: {} for name in ctx.kept}
    for name in ctx.kept:
        for gi in ctx.measured[name]:
            A, B, _, _ = t71.analytic_A_and_B_conditioned_gc(
                ctx.fixed_solutions[name]["angles"], t71.GATE_NAME, t71.ZZ_ASSUMED,
                t71.GPI_REAL_MEAN, float(gpi2), 0.0, 0.0,
                ctx.diag_natives[gi], ctx.diagonalizers[gi].transformed)
            for l in ctx.groups[gi]:
                d = float(A[l])
                out[name][l] = (float(B[l] / d) if abs(d) > 1e-6 else 1.0, gi)
    return out


def apply_alpha(post, rat, alpha_by_group):
    """corrected = raw * (1 + alpha_g * (ratio - 1)), clipped as the locked code does."""
    corr = {}
    for name, labels in rat.items():
        corr[name] = {}
        for l, (r, gi) in labels.items():
            a = alpha_by_group[gi]
            corr[name][l] = float(np.clip(post[name][l] * (1.0 + a * (r - 1.0)), -1.0, 1.0))
    return corr


N_RESTARTS = 4   # the locked pipeline's value; raise with --restarts


def frame_fit(ctx, corr):
    unit = {n: {l: 1.0 for l in ctx.non_id} for n in ctx.kept}
    U_hat, _, _ = t71.fit_joint_frame(np.eye(K), ctx.P_S, K, ctx.kept, ctx.non_id, corr, unit,
                                      np.random.default_rng(FRAME_SEED), n_restarts=N_RESTARTS)
    return U_hat


def forced_zero_mask(post, rat):
    """Labels the LOCKED pipeline (alpha=1) excludes, as a fixed set.

    The locked code drops a label when |corrected/raw| < FORCED_ZERO_RATIO.
    That quantity is 1 + alpha*(r-1), which is small only when alpha is near 1
    for the labels whose r ~ 0. So recomputing the mask per alpha excludes 24
    labels at alpha=1 and none anywhere else, and chi2 at different alpha would
    be computed over different residual sets -- not comparable. Freezing the
    mask at its alpha=1 membership keeps every alpha on the same 732 terms.
    """
    mask = set()
    for name, labels in rat.items():
        for l, (r, _gi) in labels.items():
            if abs(post[name][l]) > 1e-9 and abs(r) < FORCED_ZERO_RATIO:
                mask.add((name, l))
    return mask


def residuals(ctx, post, corr, kept_shots, U_hat, slots=None, mask=None, rat=None):
    """Shot-noise-weighted residuals.

    If `rat` is given, sigma is computed from the alpha=1 correction ratio and
    held FIXED as alpha varies. This matters: the locked pipeline sets
    sigma = sigma_raw * |corrected/raw|, which is correct error propagation for
    one given correction, but makes chi2 non-comparable across a hypothesis
    scan -- sigma vanishes wherever 1 + alpha*(r-1) = 0, i.e. at
    alpha = 1/(1-r). 26 labels in this dataset have such a pole in alpha
    within (0, 3]: 24 at alpha = 1.0 and 2 at alpha = 1.2, which is exactly
    the spurious chi2 spike seen at alpha = 1.2. The shot noise of the raw
    measurement does not change because we changed our correction hypothesis,
    so freezing sigma is the right comparison.
    """
    slots = slots or ctx.kept
    res, n_forced = [], 0
    for name in slots:
        v = slot_vector(U_hat, name, K)
        for l, y in corr[name].items():
            m = post[name][l]
            r = y / m if abs(m) > 1e-9 else 1.0
            drop = ((name, l) in mask) if mask is not None else (abs(r) < FORCED_ZERO_RATIO)
            if drop:
                n_forced += 1
                continue
            r_sig = rat[name][l][0] if rat is not None else r
            sig = np.sqrt(max(1.0 - m * m, 1e-4) / max(kept_shots[name][l], 1)) * abs(r_sig)
            res.append((float(np.real(v @ ctx.P_S[l] @ v)) - y) / max(sig, 1e-6))
    a = np.asarray(res)
    return a, n_forced


def chi2_dof(a):
    return float(a @ a / (len(a) - N_THETA))


def energies(ctx, corr, U_hat):
    full = t71.build_full_from_frame(U_hat, ctx.P_S, K, ctx.non_id, ctx.kept)
    E, e = t71.energy_report(ctx.p, full, K)
    return float(E), float(e["err_vs_exact_kcal"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--design", choices=["full", "subset"], default="full")
    ap.add_argument("--alpha-max", type=float, default=1.40)
    ap.add_argument("--alpha-step", type=float, default=0.05)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    global ALPHA_GRID, N_RESTARTS
    ALPHA_GRID = np.round(np.arange(0.0, args.alpha_max + 1e-9, args.alpha_step), 3)
    N_RESTARTS = args.restarts
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = t73.Context(args.design)
    with open(args.counts) as f:
        blob = json.load(f)
    counts = {k: dict(v) for k, v in blob["counts"].items()}
    t73.add_reused_task49(counts, ctx)

    post, kept_shots = postselect(counts, ctx)
    rat = ratios(ctx, CORRECTION_GPI2)
    ng = len(ctx.groups)
    # Frozen at alpha=1 membership so every alpha is scored on the same terms.
    # Without this the comparison is invalid: see forced_zero_mask().
    mask = forced_zero_mask(post, rat)

    print(f"\n  task83 -- correction-magnitude refit, job {blob.get('job_id')}")
    print(f"  backend {blob.get('backend')}   alpha=1 is the locked pipeline")
    print(f"  forced-zero mask frozen at alpha=1: {len(mask)} labels excluded at every alpha")
    print(f"  sigma frozen at its alpha=1 value (removes the sigma poles in alpha)\n")

    # ---- Model A: one global alpha, blind to the energy -------------------
    print("  == Model A: single global alpha (objective = chi2/dof only) ==")
    print(f"  {'alpha':>6} {'chi2/dof':>9} {'shared_frame err':>18}")
    recA = []
    for a in ALPHA_GRID:
        corr = apply_alpha(post, rat, {gi: float(a) for gi in range(ng)})
        U = frame_fit(ctx, corr)
        r, _ = residuals(ctx, post, corr, kept_shots, U, mask=mask, rat=rat)
        c = chi2_dof(r)
        _, err = energies(ctx, corr, U)
        recA.append({"alpha": float(a), "chi2_dof": c, "err_kcal": err})
        flag = "   <-- locked" if abs(a - 1.0) < 1e-9 else ""
        print(f"  {a:6.2f} {c:9.3f} {err:18.3f}{flag}")
    bestA = min(recA, key=lambda d: d["chi2_dof"])
    locked = next(d for d in recA if abs(d["alpha"] - 1.0) < 1e-9)
    print(f"\n  locked  alpha=1.00 -> chi2/dof {locked['chi2_dof']:.3f}, err {locked['err_kcal']:.3f}")
    print(f"  best    alpha={bestA['alpha']:.2f} -> chi2/dof {bestA['chi2_dof']:.3f}, "
          f"err {bestA['err_kcal']:.3f}")

    # ---- Model B: per-GC-group alpha, coordinate descent on chi2 ----------
    print("\n  == Model B: per-GC-group alpha (coordinate descent on chi2/dof) ==")
    alpha_g = {gi: 1.0 for gi in range(ng)}
    for sweep in range(3):
        for gi in range(ng):
            best = (None, np.inf)
            for a in ALPHA_GRID:
                trial = dict(alpha_g); trial[gi] = float(a)
                corr = apply_alpha(post, rat, trial)
                U = frame_fit(ctx, corr)
                r, _ = residuals(ctx, post, corr, kept_shots, U, mask=mask, rat=rat)
                c = chi2_dof(r)
                if c < best[1]:
                    best = (float(a), c)
            alpha_g[gi] = best[0]
        corr = apply_alpha(post, rat, alpha_g)
        U = frame_fit(ctx, corr)
        r, _ = residuals(ctx, post, corr, kept_shots, U, mask=mask, rat=rat)
        print(f"  sweep {sweep+1}: alpha={[round(alpha_g[g],2) for g in range(ng)]} "
              f"chi2/dof={chi2_dof(r):.3f}")
    corrB = apply_alpha(post, rat, alpha_g)
    UB = frame_fit(ctx, corrB)
    rB, _ = residuals(ctx, post, corrB, kept_shots, UB, mask=mask, rat=rat)
    _, errB = energies(ctx, corrB, UB)
    print(f"  Model B: chi2/dof {chi2_dof(rB):.3f}, err {errB:.3f}")

    # ---- Held-out validation: fit alpha on half the slots, score the other -
    print("\n  == Held-out validation (global alpha) ==")
    slots = list(ctx.kept)
    halves = [slots[0::2], slots[1::2]]
    held = []
    for i, (fit_slots, test_slots) in enumerate([(halves[0], halves[1]), (halves[1], halves[0])]):
        best = (None, np.inf)
        for a in ALPHA_GRID:
            corr = apply_alpha(post, rat, {gi: float(a) for gi in range(ng)})
            U = frame_fit(ctx, corr)
            r, _ = residuals(ctx, post, corr, kept_shots, U, slots=fit_slots, mask=mask, rat=rat)
            c = float(r @ r / max(len(r) - N_THETA, 1))
            if c < best[1]:
                best = (float(a), c)
        a_fit = best[0]
        for tag, a in (("locked alpha=1.00", 1.0), (f"fitted alpha={a_fit:.2f}", a_fit)):
            corr = apply_alpha(post, rat, {gi: float(a) for gi in range(ng)})
            U = frame_fit(ctx, corr)
            r, _ = residuals(ctx, post, corr, kept_shots, U, slots=test_slots, mask=mask, rat=rat)
            c = float(r @ r / max(len(r) - N_THETA, 1))
            print(f"  fold {i+1}: {tag:22} -> held-out chi2/dof {c:.3f}")
            held.append({"fold": i + 1, "alpha": float(a), "held_out_chi2_dof": c,
                         "fitted": tag.startswith("fitted")})

    out = args.out or os.path.join(os.path.dirname(__file__),
                       "task83_hw_correction_magnitude_refit_results.json")
    with open(out, "w") as f:
        json.dump({"job_id": blob.get("job_id"), "backend": blob.get("backend"),
                   "locked_alpha_1": locked, "model_a_grid": recA, "model_a_best": bestA,
                   "model_b_alpha_by_group": alpha_g,
                   "model_b": {"chi2_dof": chi2_dof(rB), "err_kcal": errB},
                   "held_out": held,
                   "config": {"alpha_max": args.alpha_max, "alpha_step": args.alpha_step,
                              "n_restarts": N_RESTARTS,
                              "frozen_mask_size": len(mask)}}, f, indent=2)
    print(f"\n  Saved -> {out}\n")


if __name__ == "__main__":
    main()
