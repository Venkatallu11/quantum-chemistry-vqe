#!/usr/bin/env python3
"""
task75_h4_reliability_sweep.py -- iteration 75. task73 showed pooling
all 4 real trials (4,400 effective shots/circuit) reliably clears
chemical accuracy, but that may cost more than the approved $2,352.16
budget supports (task74 checks the real price). This sweep answers the
PRACTICAL question: what's the SMALLEST real shot count, within what
we've already honestly collected, that reliably clears chemical
accuracy on ALL THREE backends -- so we know the true minimum ask
before deciding whether to request more budget.

METHOD: bootstrap-resample the POOLED real counts (4 independent real
trials x 1,100 shots = 4,400 real shots/circuit, task72) down to each
candidate shot level, N_BOOT independent draws per level, run the full
locked analysis (raw/no-frame/shared-frame) each time. This is honest
use of already-collected real data -- no new submissions, no new cost
-- and directly estimates the PASS RATE (fraction of draws clearing
chemical accuracy on all three backends at once) as a function of shot
count, which a single point estimate can't show.

Run:
    PYTHONHASHSEED=0 python vqe/task75_h4_reliability_sweep.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task72_h4_full_experiment as t72
import task73_h4_pooled_experiment as t73
from task59_h4_gc_no_frame_fit import BACKENDS
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed

N_TRIALS = 4
N_BOOT = int(os.environ.get("TASK75_N_BOOT", "12"))
CHEM_ACC_KCAL = 1.0
STRICT_KCAL = 0.25
_default_levels = "1100,1650,2200,2750,3300,3850,4400"
SHOT_LEVELS = [int(s) for s in os.environ.get("TASK75_SHOT_LEVELS", _default_levels).split(",")]
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task75_h4_reliability_sweep_results.json")


def resample_state(pooled_state, target_shots, draw_seed):
    """Bootstrap each circuit's pooled (4,400-shot) counts down to
    target_shots, independently per draw."""
    out = {"done": {}}
    for key, entry in pooled_state["done"].items():
        rng = np.random.default_rng(stable_seed("task75", key, target_shots, draw_seed))
        new_counts = [bootstrap_counts(c, target_shots, rng) for c in entry["counts"]]
        out["done"][key] = {"counts": new_counts, "group_idxs": entry["group_idxs"]}
    return out


def main():
    print("\n" + "=" * 96)
    print(f"  task75_h4_reliability_sweep.py -- pass-rate vs shot count, {N_BOOT} bootstrap draws/level")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = t72.setup()
    states = [t73.load_trial(t) for t in range(1, N_TRIALS + 1)]
    pooled_state = t73.pool_states(states)

    sweep = {}
    if os.path.exists(RESULTS_PATH):
        with open(RESULTS_PATH) as f:
            sweep = {int(k): v for k, v in json.load(f).items()}
        print(f"  loaded {len(sweep)} already-computed shot level(s) from a prior (lighter) run: {sorted(sweep)}")
    for shots in SHOT_LEVELS:
        if shots in sweep:
            print(f"  skip (already computed): shots={shots}")
            continue
        per_backend_vals = {bn: [] for bn in BACKENDS}
        n_pass_all = 0
        n_pass_strict = 0
        for draw in range(N_BOOT):
            resampled = resample_state(pooled_state, shots, draw)
            draw_vals = {}
            for bn in BACKENDS:
                entry = t72.analyze_backend(ctx, bn, resampled)
                draw_vals[bn] = abs(entry["shared_frame_err_kcal"])
                per_backend_vals[bn].append(draw_vals[bn])
            if all(v < CHEM_ACC_KCAL for v in draw_vals.values()):
                n_pass_all += 1
            if all(v < STRICT_KCAL for v in draw_vals.values()):
                n_pass_strict += 1

        pass_rate = n_pass_all / N_BOOT
        pass_rate_strict = n_pass_strict / N_BOOT
        stats = {}
        for bn in BACKENDS:
            vals = per_backend_vals[bn]
            stats[bn] = {"mean": float(np.mean(vals)), "std": float(np.std(vals, ddof=1)),
                         "q90": float(np.percentile(vals, 90)), "max": float(np.max(vals))}
        sweep[shots] = {"pass_rate_1kcal": pass_rate, "pass_rate_0.25kcal": pass_rate_strict, "per_backend": stats}
        print(f"\n  shots={shots:>5} ({shots/1100:.2f}x approved): "
              f"pass_rate(<1kcal, all 3 backends)={pass_rate:.2f}  pass_rate(<0.25kcal)={pass_rate_strict:.2f}")
        for bn in BACKENDS:
            s = stats[bn]
            print(f"    {bn:<8}: mean={s['mean']:.3f}  std={s['std']:.3f}  Q90={s['q90']:.3f}  max={s['max']:.3f}")

        with open(RESULTS_PATH, "w") as f:
            json.dump(sweep, f, indent=2)
        print(f"  (checkpointed -> {RESULTS_PATH})")

    with open(RESULTS_PATH, "w") as f:
        json.dump(sweep, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")

    passing_1kcal = [s for s in SHOT_LEVELS if sweep[s]["pass_rate_1kcal"] >= 0.95]
    print("  -- recommendation --")
    if passing_1kcal:
        print(f"  smallest shot level with >=95% pass rate (<1 kcal/mol, all 3 backends): {min(passing_1kcal)}")
    else:
        print("  NO shot level tested reached a 95% pass rate -- even 4,400 shots is not fully reliable "
              "(see pass_rate_1kcal at the top level above).")


if __name__ == "__main__":
    main()
