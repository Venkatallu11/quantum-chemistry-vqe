#!/usr/bin/env python3
"""
taskA_drift_characterization.py — iteration 25, Task A (done FIRST, per
explicit instruction: everything else depends on it). Task 5 found the
SAME Z2-tapered circuit set submitted twice, independently, to IonQ's
free ionq_simulator gave aria-1 51.05 vs 47.78 and forte-1 54.43 vs
51.50 kcal/mol -- a ~3 kcal/mol gap. That is two submissions and one
difference: an anecdote, not a distribution. This file submits the SAME
circuit set N_REPS=8 times, independently, real network calls, free
simulator, to build an actual distribution.
============================================================================
METHOD: build the Z2-tapered circuit set ONCE (identical to
z2_tapered_ionq.py's phase_targets(), reused not reimplemented), then
submit ALL N_REPS x 3 models = 24 jobs NON-BLOCKING up front (this
project's own established "submit first, retrieve later" pattern is what
makes this tractable in real wall-clock time -- retrieving 24 already-
queued jobs is far cheaper than 8 sequential submit-wait-submit-wait
cycles), then retrieve every job's counts.

For each (repetition, model), computes the raw forged energy via the
SAME 8-seed shot-noise bootstrap this project always uses (so the
per-repetition number is directly comparable to every other headline
number in this ledger) -- then looks at the SPREAD ACROSS REPETITIONS
(mean/std/min/max of the 8 per-repetition means) to isolate
submission-to-submission drift from shot noise.

DIAGNOSIS: does the "ideal" (noiseless) model drift across repetitions
too? If YES, something in this project's own pipeline is non-
deterministic -- a bug to find, not a hardware fact. If NO (ideal stays
essentially constant, only aria-1/forte-1 drift), the free simulator is
resampling a fresh noise REALIZATION per job under the SAME named
profile -- which is a genuine hardware-relevant fact (would also happen
on real QPU submissions), not a bug in this project's code.

Run:
    python vqe/taskA_drift_characterization.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from z2_tapered_ionq import build_reduced_problem, transpiled_state_prep_3q, measurement_circuit_for_group
from qforge import combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
import ef_fragment as effrag_mod  # for group_labels_qubit_wise, reused not reimplemented
from ionq_backend import connect_provider, get_simulator
from ionq_run import basis_change, pauli_expectation, IONQ_QIS_STANDARD_BASIS
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, NOISE_MODELS, submit_job, get_counts_list, bootstrap_counts,
    expectation_from_counts, stable_seed,
)

N_REPS = 8
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "taskA_drift_reps.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "taskA_drift_characterization_results.json")
K = 6


def build_circuits(problem):
    reduced_targets = problem["reduced_targets"]
    reduced_label_map = problem["reduced_label_map"]
    unique_reduced_labels = sorted(set(rl for rl, _ in reduced_label_map.values()))
    groups = effrag_mod.group_labels_qubit_wise(unique_reduced_labels)
    target_names = sorted(reduced_targets.keys())

    circuits, tags = [], []
    for name in target_names:
        base = transpiled_state_prep_3q(reduced_targets[name])
        for gi, group in enumerate(groups):
            circuits.append(measurement_circuit_for_group(base, group))
            tags.append((name, gi))
    return circuits, tags, target_names, groups


def submit_all_reps(circuits, backend):
    """Submit N_REPS x 3 models = 24 jobs, ALL non-blocking, before
    retrieving ANY of them -- maximizes queue concurrency on IonQ's side."""
    jobs = {}  # (rep, model) -> job handle
    t0 = time.time()
    for rep in range(N_REPS):
        for model in NOISE_MODELS:
            jobs[(rep, model)] = submit_job(circuits, backend, model, shots=SHOTS)
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted (non-blocking), {t_submit:.1f}s")
    return jobs, t_submit


def retrieve_all(jobs):
    t0 = time.time()
    counts = {}
    for i, (key, job) in enumerate(jobs.items()):
        counts[key] = get_counts_list(job)
        print(f"    retrieved {i+1}/{len(jobs)}: rep={key[0]} model={key[1]}, {time.time()-t0:.1f}s elapsed")
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")
    return counts, t_retrieve


def energy_from_counts(problem, counts_flat, tags, target_names, groups, seeds):
    """Same bootstrap/assembly methodology as z2_tapered_ionq.assemble()."""
    p = problem["p"]
    alpha_labels = problem["alpha_labels"]
    identity_label = problem["identity_label"]
    reduced_label_map = problem["reduced_label_map"]
    idx_map = [name for name in target_names for _ in groups]

    per_name = {name: [] for name in target_names}
    for i, name in enumerate(idx_map):
        per_name[name].append(counts_flat[i])

    errs = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        raw = {name: {} for name in target_names}
        for name in target_names:
            for gi, group in enumerate(groups):
                counts = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                group_vals = {rl: expectation_from_counts(counts, rl) for rl in group}
                for orig_label in alpha_labels:
                    rl, sign = reduced_label_map[orig_label]
                    if rl in group_vals:
                        raw[name][orig_label] = sign * group_vals[rl]
        alpha_mats = combine_matrices(raw, alpha_labels, identity_label, K)
        E, err = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                             exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
        errs.append(err["err_vs_exact_kcal"])
    return float(np.mean(errs)), float(np.std(errs))


def main():
    print("\n" + "=" * 96)
    print("  taskA_drift_characterization.py -- 8 independent real submissions, same circuit set")
    print("=" * 96)

    problem = build_reduced_problem()
    circuits, tags, target_names, groups = build_circuits(problem)
    print(f"  {len(circuits)} circuits/job, {N_REPS} repetitions x {len(NOISE_MODELS)} models = "
          f"{N_REPS*len(NOISE_MODELS)} total real jobs")

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"  connected, backend={backend.name}")

    ck = None
    if os.path.exists(CKPT_PATH):
        with open(CKPT_PATH) as f:
            ck = json.load(f)
        print(f"  found existing checkpoint with {len(ck.get('counts', {}))} (rep,model) entries -- reusing, not resubmitting")

    if ck is None:
        jobs, t_submit = submit_all_reps(circuits, backend)
        counts, t_retrieve = retrieve_all(jobs)
        ck = {
            "target_names": target_names, "groups": groups, "tags": [list(t) for t in tags],
            "n_reps": N_REPS, "shots": SHOTS,
            "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
            "counts": {f"{rep}|{model}": counts[(rep, model)] for rep, model in counts},
        }
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- per-repetition energy (8-seed shot-noise bootstrap, standard convention) --
    per_rep = {model: [] for model in NOISE_MODELS}
    for rep in range(N_REPS):
        for model in NOISE_MODELS:
            counts_flat = ck["counts"][f"{rep}|{model}"]
            seeds = [stable_seed("taskA", model, rep, s) for s in range(N_SEEDS)]
            mean, std = energy_from_counts(problem, counts_flat, tags, target_names, groups, seeds)
            per_rep[model].append({"rep": rep, "mean": mean, "shot_noise_std": std})
            print(f"    rep={rep} model={model}: {mean:.2f} +/- {std:.2f} kcal/mol (shot noise only)")

    print(f"\n  -- SUBMISSION-TO-SUBMISSION DRIFT (spread of the {N_REPS} per-repetition means) --")
    drift_stats = {}
    for model in NOISE_MODELS:
        means = [r["mean"] for r in per_rep[model]]
        drift_stats[model] = {
            "rep_means": means, "mean": float(np.mean(means)), "std": float(np.std(means)),
            "min": float(np.min(means)), "max": float(np.max(means)), "range": float(np.max(means) - np.min(means)),
            "mean_shot_noise_std": float(np.mean([r["shot_noise_std"] for r in per_rep[model]])),
        }
        d = drift_stats[model]
        print(f"    {model}: mean={d['mean']:.2f}  drift_std={d['std']:.2f}  min={d['min']:.2f}  max={d['max']:.2f}  "
              f"range={d['range']:.2f}  (mean shot-noise std per rep: {d['mean_shot_noise_std']:.2f})")

    print(f"\n  -- DIAGNOSIS: does ideal drift too? --")
    ideal_drift = drift_stats["ideal"]["std"]
    noisy_drift = (drift_stats["aria-1"]["std"] + drift_stats["forte-1"]["std"]) / 2
    if ideal_drift > 0.5 * noisy_drift and ideal_drift > 0.3:
        diagnosis = ("ideal ALSO drifts substantially (std={:.2f} vs noisy-model avg std={:.2f}) -- "
                      "this points to a PIPELINE bug (non-determinism somewhere in this project's own code "
                      "or in how the free simulator handles the noiseless case), not purely a hardware fact."
                      ).format(ideal_drift, noisy_drift)
    else:
        diagnosis = ("ideal is essentially stable (std={:.2f}) while aria-1/forte-1 drift much more "
                      "(avg std={:.2f}) -- consistent with the free simulator resampling a FRESH NOISE "
                      "REALIZATION per job under the same named profile, which is a genuine, real "
                      "characteristic that would also matter on real hardware, not a bug in this pipeline."
                      ).format(ideal_drift, noisy_drift)
    print(f"    {diagnosis}")

    # -- drift-aware combined error bar: shot-noise-bootstrap-std COMBINED (quadrature) with drift-std --
    print(f"\n  -- DRIFT-AWARE ERROR BAR (combined in quadrature: sqrt(shot_noise_std^2 + drift_std^2)) --")
    combined = {}
    for model in NOISE_MODELS:
        d = drift_stats[model]
        combined_std = float(np.sqrt(d["mean_shot_noise_std"] ** 2 + d["std"] ** 2))
        combined[model] = {"mean": d["mean"], "combined_std": combined_std}
        print(f"    {model}: {d['mean']:.2f} +/- {combined_std:.2f} kcal/mol (was +/- {d['mean_shot_noise_std']:.2f} shot-noise-only)")

    # -- re-quote the ablation table gaps against this new bar --
    print(f"\n  -- ARE THE ABLATION TABLE'S GAPS STILL DISTINGUISHABLE? --")
    # gaps as stated by the user: raw+leakage vs PSD+leakage, aria-1=2.19, forte-1=3.10 kcal/mol
    stated_gaps = {"aria-1": 2.19, "forte-1": 3.10}
    for model in ["aria-1", "forte-1"]:
        gap = stated_gaps[model]
        bar = combined[model]["combined_std"]
        distinguishable = gap > bar
        print(f"    {model}: gap={gap:.2f} kcal/mol vs drift-aware bar=+/-{bar:.2f} kcal/mol -- "
              f"{'gap EXCEEDS the bar (distinguishable)' if distinguishable else 'gap is INSIDE the bar -- NOT distinguishable'}")
        stated_gaps[model] = {"gap": gap, "drift_aware_bar": bar, "distinguishable": bool(distinguishable)}

    results = {
        "n_reps": N_REPS, "per_rep": per_rep, "drift_stats": drift_stats, "diagnosis": diagnosis,
        "combined_error_bar": combined, "ablation_gap_check": stated_gaps,
        "checkpoint": CKPT_PATH,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
