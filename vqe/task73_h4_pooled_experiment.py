#!/usr/bin/env python3
"""
task73_h4_pooled_experiment.py -- iteration 73. Task72's single real
1,100-shot/circuit draw of the approved 54-circuit design was NOT
reliably inside chemical accuracy: of 4 independent real trials run on
IonQ's free ionq_simulator (task72 --trial 1..4, all reproducibility-
verified), forte-1's shared-frame error was 0.42/1.75/0.17/0.53 kcal/mol
-- ONE of four draws (25%) exceeds the 1 kcal/mol chemical-accuracy bar.
This is real shot noise, not a bug (same pattern this project already
documented for other designs: see RESEARCH_LEDGER's drift sections).

FIX TESTED HERE: pool the real counts from all 4 independent trials
(summing each circuit's bitstring counts across trials -- a real,
honest aggregation of genuinely-collected data, not a re-weighting or
cherry-pick) to get ~4,400 effective shots/circuit, then run the exact
same locked analysis once on the pooled data. This answers: does simply
collecting more real shots (rather than changing the circuit design)
make the approved plan reliably clear chemical accuracy?

Also reports the 4 individual-trial numbers side by side (mean/std/Q50/
Q90, this project's standard drift-characterization), so the pooled
result's improvement is visible against the per-trial spread it's
curing.

Run:
    PYTHONHASHSEED=0 python vqe/task73_h4_pooled_experiment.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import task72_h4_full_experiment as t72
from task59_h4_gc_no_frame_fit import BACKENDS

N_TRIALS = 4
CHEM_ACC_KCAL = 1.0  # standard chemical-accuracy bar


def load_trial(trial):
    with open(t72.ckpt_path(trial)) as f:
        return json.load(f)


def pool_states(states):
    pooled = {"done": {}}
    keys = states[0]["done"].keys()
    for key in keys:
        group_idxs = states[0]["done"][key]["group_idxs"]
        n_circuits = len(group_idxs)
        merged_counts = []
        for ci in range(n_circuits):
            merged = {}
            for st in states:
                for bitstring, c in st["done"][key]["counts"][ci].items():
                    merged[bitstring] = merged.get(bitstring, 0) + c
            merged_counts.append(merged)
        pooled["done"][key] = {"counts": merged_counts, "group_idxs": group_idxs}
    return pooled


def main():
    print("\n" + "=" * 96)
    print(f"  task73_h4_pooled_experiment.py -- pooling {N_TRIALS} real trials ({N_TRIALS * t72.PROPOSED_SHOTS} "
          f"effective shots/circuit)")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = t72.setup()

    states = [load_trial(t) for t in range(1, N_TRIALS + 1)]
    for t, st in enumerate(states, start=1):
        missing = [bn for bn in BACKENDS if any(f"{bn}|{name}" not in st["done"] for name in ctx["kept"])]
        assert not missing, f"trial {t} checkpoint incomplete: {missing}"

    # ---- individual-trial table (drift characterization) ----
    per_trial = {bn: {"raw": [], "no_frame": [], "shared_frame": []} for bn in BACKENDS}
    for t, st in enumerate(states, start=1):
        for bn in BACKENDS:
            entry = t72.analyze_backend(ctx, bn, st)
            per_trial[bn]["raw"].append(entry["raw_err_kcal"])
            per_trial[bn]["no_frame"].append(entry["no_frame_err_kcal"])
            per_trial[bn]["shared_frame"].append(entry["shared_frame_err_kcal"])

    print(f"\n  -- {N_TRIALS} individual real trials, shared-frame error (kcal/mol) --")
    for bn in BACKENDS:
        vals = per_trial[bn]["shared_frame"]
        n_fail = sum(1 for v in vals if abs(v) > CHEM_ACC_KCAL)
        print(f"    {bn:<8}: {['%.4f' % v for v in vals]}  mean={np.mean(vals):.4f}  std={np.std(vals, ddof=1):.4f}  "
              f"max={max(vals, key=abs):.4f}  fails_1kcal={n_fail}/{N_TRIALS}")

    # ---- pooled real data, one combined analysis ----
    pooled_state = pool_states(states)
    print(f"\n  -- pooled analysis, {N_TRIALS * t72.PROPOSED_SHOTS} effective real shots/circuit --")
    pooled_results = {}
    for bn in BACKENDS:
        entry = t72.analyze_backend(ctx, bn, pooled_state)
        pooled_results[bn] = entry
        verdict = "PASS" if abs(entry["shared_frame_err_kcal"]) < CHEM_ACC_KCAL else "FAIL"
        print(f"    {bn:<8}: accept={entry['mean_accept']:.4f}  raw={entry['raw_err_kcal']:+.4f}  "
              f"no_frame={entry['no_frame_err_kcal']:+.4f}  shared_frame={entry['shared_frame_err_kcal']:+.4f} "
              f"kcal/mol  chi2/dof={entry['chi2_dof']:.4f}  [{verdict}, <{CHEM_ACC_KCAL} kcal/mol]")

    # reproducibility re-check on the pooled analysis too -- same discipline as task72
    pooled_results2 = {bn: t72.analyze_backend(ctx, bn, pooled_state) for bn in BACKENDS}
    for bn in BACKENDS:
        for k in ("raw_err_kcal", "no_frame_err_kcal", "shared_frame_err_kcal"):
            assert abs(pooled_results[bn][k] - pooled_results2[bn][k]) < 1e-9, \
                f"pooled reproducibility check FAILED for {bn}.{k}"
    print("\n  pooled reproducibility check: two independent passes agree exactly -- PASS")

    out = {
        "per_trial": {bn: {k: per_trial[bn][k] for k in ("raw", "no_frame", "shared_frame")} for bn in BACKENDS},
        "pooled": pooled_results,
        "chem_acc_kcal_bar": CHEM_ACC_KCAL,
        "n_trials_pooled": N_TRIALS,
        "effective_shots_per_circuit": N_TRIALS * t72.PROPOSED_SHOTS,
    }
    results_path = os.path.join(os.path.dirname(__file__), "task73_h4_pooled_experiment_results.json")
    with open(results_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Saved -> {results_path}\n")


if __name__ == "__main__":
    main()
