#!/usr/bin/env python3
"""
task73_hw_analysis.py -- iteration 73. The analysis for the real-hardware H4
benchmark, committed BEFORE any hardware data exists, so the procedure is
fixed in advance.

INPUT (real run): a JSON file of raw hardware counts,
    {"job_id": ..., "backend": ..., "counts": {"<slot>|g<group>": {bitstring: count}}}
with the same 5-qubit layout as task59/task49 (bitstring[0] = parity ancilla,
bitstring[1:] = the 4-qubit register after the group's diagonalizer).
task49's three already-run GC3 circuits (u_0, u_1, u_2 on Forte Enterprise 1,
10 Sep 2026, 2,000 shots) fill their keys if the input lacks them.

DESIGN
    full   (default): all 4 GC groups on all 21 slots = 84 circuits
                      (81 new + 3 reused), funded by the extra $3,000
    subset          : task71's 54-circuit plan (GC0+GC1 everywhere, GC2+GC3 on
                      the 6-slot panel)

LOCKED ANALYSIS (same code paths as the task71/72 rehearsal)
    raw          per-slot physical-state fit (variance-weighted) to the
                 postselected, uncorrected values
    no_frame     same fit after the GC-aware correction (task70), noise
                 parameters FIXED here, never re-tuned on hardware data:
                 ZZ 0.014593, GPi 0.000119 (task30B), GPi2 0.0004 (forte-1)
    shared_frame joint 15-parameter Schmidt-frame fit to the corrected values
                 (unit weights, 4 restarts, as rehearsed)
    Uncertainty: nonparametric bootstrap -- every circuit's counts resampled
    at its own shot total, whole pipeline rerun, N_BOOT times.
    Fit quality: shot-noise-weighted chi2/dof of the shared-frame fit, labels
    with correction ratio |B/A| < 0.01 excluded and counted (task72).

GATES (analysis stops if either fails)
    1. Schmidt sign convention: stored ideal data fits the model at the exact
       frame (task72; guards the pinned signs in qforge.forging).
    2. Data completeness: every circuit the design requires has counts.

EXACT-STATE INFORMATION, per analysis
    all three: circuits prepare the FCI Schmidt states and pair combinations;
               fits live in the exact Schmidt-vector span; the energy uses the
               exact Schmidt coefficients and beta signs.
    no_frame : + correction factors simulated from the intended circuit.
    shared   : + one shared frame, fit started at the exact frame, depends on
               the Schmidt vector signs.

REHEARSAL MODE (no hardware data needed): --rehearsal {ideal,aria-1,forte-1}
builds the same input from the stored task59 noise-model data bootstrapped to
1,100 shots and runs the identical pipeline (with that noise model's own GPi2,
task59 GPI2_SELECTED) -- the simulator prediction to
compare the hardware against. On the subset design, the rehearsal point
estimate reproduces task72 trial 0 exactly.

Run:
    PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py --counts vqe/task73_hw_counts.json
    PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py --rehearsal forte-1
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task71_hw_plan_rehearsal_v2 as t71
from task72_rehearsal_sign_and_chi2_check import exact_frame_gate, FORCED_ZERO_RATIO, EXACT_FRAME_GATE
from task36_joint_schmidt_frame import slot_vector, N_THETA
from task60_ionq_no_frame_h4 import fit_all_slots, build_full_from_independent_states, variance_weights
from general_commuting_measurements import expectations_from_counts
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

HERE = os.path.dirname(__file__)
K = t71.K
CORRECTION_GPI2 = 0.0004          # forte-1 selection (task59 GPI2_SELECTED), fixed before the run
REHEARSAL_SHOTS = t71.PROPOSED_SHOTS  # 1,100
N_BOOT_DEFAULT = 20
FRAME_SEED = 7100
MIN_RETENTION = 0.70              # a circuit keeping fewer postselected shots than this is flagged
REUSED_TASK49_PATH = os.path.join(HERE, "task49_gc_hard_group_hw_submission_results.json")
SIM_PREDICTION_PATH = os.path.join(HERE, "task71_hw_plan_rehearsal_v2_results.json")


def key(name, gi):
    return f"{name}|g{gi}"


class Context:
    def __init__(self, design):
        self.p = t71.setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
        self.non_id = sorted(l for l in self.p["alpha_labels"] if l != self.p["identity_label"])
        self.fixed_solutions, _, _ = t71.fit_all_targets(self.p["targets"], tol=1e-10)
        self.diag, _, self.kept = t71.kept_slots_for_K(K)
        self.P_S = t71.build_P_S(self.p["alpha_labels"], np.asarray(self.p["u_vecs"]).T)
        self.groups, self.diagonalizers = t71.build_general_commuting_measurement_plan(self.non_id)
        self.diag_natives = [t71.to_native(d.to_circuit(), t71.GATE_NAME) for d in self.diagonalizers]
        n_groups = len(self.groups)
        if design == "full":
            self.measured = {name: list(range(n_groups)) for name in self.kept}
        else:
            self.measured = {name: t71.measured_groups_for_slot(name, n_groups) for name in self.kept}
        self.required = [key(n, gi) for n in self.kept for gi in self.measured[n]]


def add_reused_task49(counts, ctx):
    """Fill u_0/u_1/u_2 GC3 from the 10 Sep hardware job when the input lacks them."""
    with open(REUSED_TASK49_PATH) as f:
        d = json.load(f)
    gi = d["hard_group_index"]
    assert sorted(ctx.groups[gi]) == sorted(d["group"]), "task49 group no longer matches today's GC plan"
    added = []
    for slot, c in zip(d["target_slots"], d["counts_list"]):
        k = key(slot, gi)
        if k in ctx.required and k not in counts:
            counts[k] = dict(c)
            added.append(k)
    return added


def rehearsal_counts(backend, ctx):
    with open(t71.CKPT_PATH) as f:
        state = json.load(f)
    counts = {}
    for name in ctx.kept:
        entry = state["done"][f"{backend}|{name}"]
        for gi in ctx.measured[name]:
            rng = np.random.default_rng(stable_seed("task71", backend, name, gi, 0))
            counts[key(name, gi)] = bootstrap_counts(entry["counts"][gi], REHEARSAL_SHOTS, rng)
    return counts


def data_quality(counts, ctx):
    missing = [k for k in ctx.required if k not in counts]
    retention = {}
    for k in ctx.required:
        if k in counts:
            tot = sum(counts[k].values())
            retention[k] = sum(c for bs, c in counts[k].items() if bs[0] == "0") / max(tot, 1)
    low = sorted(k for k, r in retention.items() if r < MIN_RETENTION)
    return {"n_required": len(ctx.required), "n_present": len(ctx.required) - len(missing), "missing": missing,
            "mean_retention": float(np.mean(list(retention.values()))) if retention else None,
            "min_retention": float(min(retention.values())) if retention else None,
            "low_retention_circuits": low, "shots_per_circuit": sorted({sum(counts[k].values()) for k in counts})}


def analyze(counts, ctx, gpi2):
    """One pass of the locked pipeline on one set of counts."""
    post = {n: {} for n in ctx.kept}
    kept_shots = {n: {} for n in ctx.kept}
    for name in ctx.kept:
        for gi in ctx.measured[name]:
            filtered = {bs[1:]: c for bs, c in counts[key(name, gi)].items() if bs[0] == "0"}
            n_kept = sum(filtered.values())
            exp = expectations_from_counts(filtered, ctx.diagonalizers[gi]) if filtered else {}
            for l in ctx.groups[gi]:
                post[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                kept_shots[name][l] = n_kept

    fits_raw = fit_all_slots(post, variance_weights(post, kept_shots), ctx.P_S, ctx.kept)
    full_raw = build_full_from_independent_states(fits_raw, ctx.diag, K, ctx.P_S, ctx.non_id)
    E_raw, e_raw = t71.energy_report(ctx.p, full_raw, K)

    if gpi2 is not None:
        corr = t71.corrected_observables_gc_aware(post, gpi2, ctx.kept, ctx.fixed_solutions,
                                                  ctx.diag_natives, ctx.diagonalizers, ctx.groups, ctx.measured)
    else:
        corr = post
    fits_nf = fit_all_slots(corr, variance_weights(corr, kept_shots), ctx.P_S, ctx.kept)
    full_nf = build_full_from_independent_states(fits_nf, ctx.diag, K, ctx.P_S, ctx.non_id)
    E_nf, e_nf = t71.energy_report(ctx.p, full_nf, K)

    unit = {n: {l: 1.0 for l in ctx.non_id} for n in ctx.kept}
    U_hat, _, msr = t71.fit_joint_frame(np.eye(K), ctx.P_S, K, ctx.kept, ctx.non_id, corr, unit,
                                        np.random.default_rng(FRAME_SEED), n_restarts=4)
    full_sf = t71.build_full_from_frame(U_hat, ctx.P_S, K, ctx.non_id, ctx.kept)
    E_sf, e_sf = t71.energy_report(ctx.p, full_sf, K)

    resid, n_forced = [], 0
    for name in ctx.kept:
        v = slot_vector(U_hat, name, K)
        for l, y in corr[name].items():
            m = post[name][l]
            ratio = y / m if abs(m) > 1e-9 else 1.0
            if abs(ratio) < FORCED_ZERO_RATIO:
                n_forced += 1
                continue
            sigma = np.sqrt(max(1.0 - m * m, 1e-4) / max(kept_shots[name][l], 1)) * abs(ratio)
            resid.append((float(np.real(v @ ctx.P_S[l] @ v)) - y) / max(sigma, 1e-6))
    resid = np.asarray(resid)
    return {"raw": {"E_Ha": float(E_raw), "err_kcal": float(e_raw["err_vs_exact_kcal"])},
            "no_frame": {"E_Ha": float(E_nf), "err_kcal": float(e_nf["err_vs_exact_kcal"])},
            "shared_frame": {"E_Ha": float(E_sf), "err_kcal": float(e_sf["err_vs_exact_kcal"]),
                             "chi2_dof": float(resid @ resid / (len(resid) - N_THETA)),
                             "n_chi2_labels": int(len(resid)), "n_forced_zero_excluded": int(n_forced),
                             "unit_weight_msr": float(msr)}}


def bootstrap(counts, ctx, n_boot, gpi2, seed_tag):
    draws = []
    for b in range(n_boot):
        resampled = {}
        for k in ctx.required:
            tot = sum(counts[k].values())
            rng = np.random.default_rng(stable_seed("task73", seed_tag, k, b))
            resampled[k] = bootstrap_counts(counts[k], tot, rng)
        draws.append(analyze(resampled, ctx, gpi2))
        print(f"    bootstrap {b+1}/{n_boot}: raw={draws[-1]['raw']['err_kcal']:+.3f}  "
              f"no_frame={draws[-1]['no_frame']['err_kcal']:+.3f}  "
              f"shared_frame={draws[-1]['shared_frame']['err_kcal']:+.3f}", flush=True)
    out = {}
    for a in ("raw", "no_frame", "shared_frame"):
        errs = np.array([d[a]["err_kcal"] for d in draws])
        out[a] = {"err_kcal_boot_mean": float(errs.mean()), "err_kcal_boot_std": float(errs.std(ddof=1)),
                  "err_kcal_boot_q05": float(np.quantile(errs, 0.05)),
                  "err_kcal_boot_q95": float(np.quantile(errs, 0.95))}
    chis = np.array([d["shared_frame"]["chi2_dof"] for d in draws])
    out["shared_frame"]["chi2_dof_boot_mean"] = float(chis.mean())
    return out, draws


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--counts", help="real hardware counts JSON")
    src.add_argument("--rehearsal", choices=t71.BACKENDS, help="run the same pipeline on stored noise-model data")
    ap.add_argument("--design", choices=["full", "subset"], default="full")
    ap.add_argument("--n-boot", type=int, default=N_BOOT_DEFAULT)
    ap.add_argument("--out", help="results JSON path (default derived from mode)")
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- run with PYTHONHASHSEED=0 for reproducible results")

    ctx = Context(args.design)
    mode = "hardware" if args.counts else f"rehearsal:{args.rehearsal}"
    print(f"\n  task73 -- {mode}, design={args.design} ({len(ctx.required)} circuits), n_boot={args.n_boot}")

    with open(t71.CKPT_PATH) as f:
        gate = exact_frame_gate(json.load(f), ctx.kept, ctx.groups, ctx.diagonalizers, ctx.P_S, K, ctx.non_id)
    print(f"  gate 1 (Schmidt signs): exact-frame residual {gate:.2e} (limit {EXACT_FRAME_GATE:.0e})")
    assert gate < EXACT_FRAME_GATE, "Schmidt sign convention does not match the model -- STOP"

    meta = {}
    if args.counts:
        with open(args.counts) as f:
            hw = json.load(f)
        counts = {k: dict(v) for k, v in hw["counts"].items()}
        meta = {k: hw.get(k) for k in ("job_id", "backend", "submitted_at", "completed_at", "billed_cost_usd")}
        reused = add_reused_task49(counts, ctx)
        gpi2 = CORRECTION_GPI2
        seed_tag = hw.get("job_id", "hw")
    else:
        counts = rehearsal_counts(args.rehearsal, ctx)
        reused = []
        gpi2 = None if args.rehearsal == "ideal" else t71.GPI2_SELECTED[args.rehearsal]
        seed_tag = f"rehearsal-{args.rehearsal}"

    dq = data_quality(counts, ctx)
    print(f"  gate 2 (completeness): {dq['n_present']}/{dq['n_required']} circuits present; "
          f"mean retention {dq['mean_retention']:.3f}, min {dq['min_retention']:.3f}; reused: {reused}")
    if dq["low_retention_circuits"]:
        print(f"  WARNING: retention < {MIN_RETENTION}: {dq['low_retention_circuits']}")
    assert not dq["missing"], f"missing circuits: {dq['missing']} -- STOP"

    point = analyze(counts, ctx, gpi2)
    print(f"  point estimate: raw={point['raw']['err_kcal']:+.3f}  no_frame={point['no_frame']['err_kcal']:+.3f}  "
          f"shared_frame={point['shared_frame']['err_kcal']:+.3f} kcal/mol  "
          f"chi2/dof={point['shared_frame']['chi2_dof']:.2f}")
    boot, draws = bootstrap(counts, ctx, args.n_boot, gpi2, seed_tag) if args.n_boot > 1 else ({}, [])

    with open(SIM_PREDICTION_PATH) as f:
        sim = json.load(f)
    results = {"mode": mode, "design": args.design, "n_circuits": len(ctx.required), "meta": meta,
               "reused_task49": reused, "gate_exact_frame_msr": gate, "data_quality": dq,
               "correction_params": {"zz": t71.ZZ_ASSUMED, "gpi": t71.GPI_REAL_MEAN, "gpi2": gpi2},
               "point_estimate": point, "bootstrap": boot, "bootstrap_draws": draws,
               "simulator_prediction_task71_subset": {b: {a: sim[b][a]["mean_abs"] for a in
                                                          ("raw", "no_frame", "shared_frame")} for b in sim}}
    out = args.out or os.path.join(HERE, "task73_hw_analysis_results.json" if args.counts else
                                   f"task73_rehearsal_{args.rehearsal}_{args.design}_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  == summary ({mode}, {args.design}) ==")
    for a in ("raw", "no_frame", "shared_frame"):
        b = boot.get(a, {})
        extra = f"  bootstrap {b['err_kcal_boot_mean']:+.3f} ± {b['err_kcal_boot_std']:.3f}" if b else ""
        print(f"    {a:<13}: point {point[a]['err_kcal']:+.3f} kcal/mol{extra}")
    print(f"  Saved -> {out}\n")


if __name__ == "__main__":
    main()
