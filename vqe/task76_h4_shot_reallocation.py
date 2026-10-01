#!/usr/bin/env python3
"""
task76_h4_shot_reallocation.py -- iteration 76. task74 showed the
approved $2,352.16 budget has no headroom for more shots: the current
plan already spends it uniformly (1,100 shots on all 54 circuits). This
tests whether REALLOCATING the SAME total dollar cost non-uniformly --
more shots on circuits whose real measured values are closer to 0 (high
shot-noise variance, 1-m^2 near 1), fewer on circuits already pinned
near +/-1 (low variance) -- buys reliability without spending a cent
more. Classic Neyman (stratified-sampling) allocation: N_i proportional
to sqrt(V_i / cost_i), V_i estimated from the REAL already-collected
pooled data (task73, 4,400 real shots/circuit -- plenty to resample any
N_i <= 4,400 from honestly).

Tests the reallocated plan the same way task75 tested the uniform one:
bootstrap-resample the real pooled counts down to the new per-circuit
shot counts, multiple independent draws, report the pass rate.

Run:
    PYTHONHASHSEED=0 python vqe/task76_h4_shot_reallocation.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task72_h4_full_experiment as t72
import task73_h4_pooled_experiment as t73
from task59_h4_gc_no_frame_fit import BACKENDS
from general_commuting_measurements import expectations_from_counts
from task74_real_hardware_cost_check import build_gate_count_table
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

UNIFORM_SHOTS = 1100
MIN_SHOTS = 300          # floor so no circuit becomes too noisy to constrain the fit
MAX_SHOTS = 4400         # ceiling: this is all the real data we have to resample from
N_DRAWS = int(os.environ.get("TASK76_N_DRAWS", "8"))
CHEM_ACC_KCAL = 1.0
STRICT_KCAL = 0.25
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task76_h4_shot_reallocation_results.json")


def estimate_variance_per_circuit(ctx, pooled_state):
    """V_i = mean over this circuit's labels of (1 - m^2), m from the
    REAL pooled (4,400-shot) data -- the lowest-noise real estimate we
    have of each label's true value."""
    variances = {}
    for key, entry in pooled_state["done"].items():
        backend_name, name = key.split("|", 1)
        if backend_name != "aria-1":
            continue  # allocate shots based on the noisier real backend, not ideal
        for local_i, gi in enumerate(entry["group_idxs"]):
            dg = ctx["diagonalizers"][gi]
            filtered = {bs[1:]: c for bs, c in entry["counts"][local_i].items() if bs[0] == "0"}
            exp = expectations_from_counts(filtered, dg) if filtered else {}
            vals = [max(-1.0, min(1.0, exp.get(l, 0.0))) for l in ctx["groups"][gi]]
            v = float(np.mean([1.0 - m ** 2 for m in vals])) if vals else 1.0
            variances[(name, gi)] = max(v, 0.02)  # floor to avoid a near-zero allocation
    return variances


def neyman_allocation(variances, costs, total_budget, min_shots=MIN_SHOTS, max_shots=MAX_SHOTS):
    keys = list(variances.keys())
    raw = {k: np.sqrt(variances[k] / max(costs[k], 1e-9)) for k in keys}
    # solve for proportionality constant lambda s.t. sum(lambda*raw_k * cost_k) = total_budget
    denom = sum(raw[k] * costs[k] for k in keys)
    lam = total_budget / denom
    shots = {k: int(np.clip(round(lam * raw[k]), min_shots, max_shots)) for k in keys}
    return shots


def main():
    print("\n" + "=" * 96)
    print("  task76_h4_shot_reallocation.py -- Neyman-optimal shot allocation at FIXED total cost")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = t72.setup()
    states = [t73.load_trial(t) for t in range(1, 5)]
    pooled_state = t73.pool_states(states)

    table = build_gate_count_table(ctx)  # (name, gi) -> (n1q, n2q)
    r1, r2 = t72_rates()
    costs = {k: (n1q * r1 + n2q * r2) for k, (n1q, n2q) in table.items()}  # $/shot

    uniform_total = sum(costs[k] * UNIFORM_SHOTS for k in costs)
    print(f"  uniform baseline: {UNIFORM_SHOTS} shots/circuit, total variable cost = ${uniform_total:.2f}")

    variances = estimate_variance_per_circuit(ctx, pooled_state)
    shots_alloc = neyman_allocation(variances, costs, uniform_total)
    realloc_total = sum(costs[k] * shots_alloc[k] for k in shots_alloc)
    print(f"  Neyman reallocation: total variable cost = ${realloc_total:.2f} "
          f"(matched to uniform baseline within rounding)")
    print(f"  shot range: min={min(shots_alloc.values())}  max={max(shots_alloc.values())}  "
          f"mean={np.mean(list(shots_alloc.values())):.0f}")

    def resample_with_allocation(draw_seed):
        out = {"done": {}}
        for key, entry in pooled_state["done"].items():
            backend_name, name = key.split("|", 1)
            new_counts = []
            for local_i, gi in enumerate(entry["group_idxs"]):
                target = shots_alloc.get((name, gi), UNIFORM_SHOTS)
                rng = np.random.default_rng(stable_seed("task76", key, gi, draw_seed))
                new_counts.append(bootstrap_counts(entry["counts"][local_i], target, rng))
            out["done"][key] = {"counts": new_counts, "group_idxs": entry["group_idxs"]}
        return out

    def resample_uniform(draw_seed):
        out = {"done": {}}
        for key, entry in pooled_state["done"].items():
            new_counts = []
            for local_i in range(len(entry["group_idxs"])):
                rng = np.random.default_rng(stable_seed("task76u", key, local_i, draw_seed))
                new_counts.append(bootstrap_counts(entry["counts"][local_i], UNIFORM_SHOTS, rng))
            out["done"][key] = {"counts": new_counts, "group_idxs": entry["group_idxs"]}
        return out

    def run_draws(resampler, label):
        n_pass_all = 0
        n_pass_strict = 0
        per_backend = {bn: [] for bn in BACKENDS}
        for draw in range(N_DRAWS):
            state = resampler(draw)
            vals = {}
            for bn in BACKENDS:
                entry = t72.analyze_backend(ctx, bn, state)
                vals[bn] = abs(entry["shared_frame_err_kcal"])
                per_backend[bn].append(vals[bn])
            if all(v < CHEM_ACC_KCAL for v in vals.values()):
                n_pass_all += 1
            if all(v < STRICT_KCAL for v in vals.values()):
                n_pass_strict += 1
        print(f"\n  -- {label}, {N_DRAWS} draws, same ${uniform_total:.0f} budget --")
        print(f"    pass_rate(<1kcal)={n_pass_all / N_DRAWS:.2f}  pass_rate(<0.25kcal)={n_pass_strict / N_DRAWS:.2f}")
        for bn in BACKENDS:
            v = per_backend[bn]
            print(f"    {bn:<8}: mean={np.mean(v):.3f}  std={np.std(v, ddof=1):.3f}  max={np.max(v):.3f}")
        return {"pass_rate_1kcal": n_pass_all / N_DRAWS, "pass_rate_0.25kcal": n_pass_strict / N_DRAWS,
                "per_backend": {bn: {"mean": float(np.mean(v)), "std": float(np.std(v, ddof=1)),
                                       "max": float(np.max(v))} for bn, v in per_backend.items()}}

    uniform_results = run_draws(resample_uniform, "UNIFORM (current approved plan)")
    realloc_results = run_draws(resample_with_allocation, "NEYMAN REALLOCATION (same budget)")

    out = {"uniform": uniform_results, "reallocated": realloc_results,
           "uniform_total_cost": uniform_total, "realloc_total_cost": realloc_total,
           "shot_allocation": {f"{k[0]}|g{k[1]}": v for k, v in shots_alloc.items()}}
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


def t72_rates():
    with open(os.path.join(os.path.dirname(__file__), "task74_cost_rate_check_results.json")) as f:
        d = json.load(f)
    return d["rate_1q"], d["rate_2q"]


if __name__ == "__main__":
    main()
