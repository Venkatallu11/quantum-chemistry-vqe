#!/usr/bin/env python3
"""
task27c_full_h4_folds.py — iteration 27, Task C. The FULL H4 forged
energy (not the Bell proxy) at native folds 1/3/5/7/9, on
ideal/aria-1/forte-1, CONCURRENTLY, on IonQ's free `ionq_simulator`. NO
real QPU submission anywhere in this file.
============================================================================
CIRCUIT-COUNT MINIMIZATION, per explicit instruction ("Minimise CIRCUIT
count where the physics allows"): reuses iteration 24 Task 3's subspace-
tomography design (diagonal + "+"-pair slots ONLY, dropping the
redundant "-" circuits -- Re<u_n|P|u_m> = <P>_+ - (M_nn+M_mm)/2 is exact
algebra, verified there) for BOTH K=5 (5 diag + 10 pairs = 15 slots, down
from 25) and K=6 (6 diag + 15 pairs = 21 slots, down from 36). This is
NOT a gate-count optimization (explicitly out of scope per this task's
own instruction) -- it is a circuit-COUNT reduction, the variable
iteration 26 Task 5 found actually correlates with real mitigated
results.

ONE GATESET PER MODEL, NEVER MIXED (Forte native = GPi/GPi2/ZZ, Aria =
GPi/GPi2/MS): `ideal` and `aria-1` are submitted through MS-native
circuits, `forte-1` through ZZ-native circuits -- two separate circuit
sets, reusing `task27ab_native_stateprep_fold`'s verified native
construction and `task2_fold_response_dataset`'s established
group/basis-change/fold-circuit builder unchanged.

CONCURRENCY: every (K, fold, model) job is submitted non-blocking BEFORE
any job is retrieved (this project's established pattern since iteration
25's Task A), so total wall-clock is closer to the SLOWEST single job
than the sum of all of them.

Run:
    python vqe/task27c_full_h4_folds.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task2_fold_response_dataset import build_folded_measurement_circuits
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts

K_VALUES = [5, 6]
FOLD_FACTORS = [1, 3, 5, 7, 9]
SHOTS = 100_000   # matches this project's established "clean signal" convention (ionq_native_forged_energy.py);
                   # a full-blown 300k-shot Task-1-style submission for a fold-response CHARACTERIZATION
                   # task (not a final chemical-accuracy claim) was judged not worth 3x the circuit volume here
GATE_BY_MODEL = {"ideal": "ms", "aria-1": "ms", "forte-1": "zz"}
N_SEEDS = 8
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task27c_full_h4_folds_results.json")


def kept_slots_for_K(K):
    diag = [f"u_{n}" for n in range(K)]
    plus = [f"(u{n}+u{m})" for n in range(K) for m in range(K) if n < m]
    return diag, plus, diag + plus


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def build_problem(K):
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=(K == 6))
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == len(p["targets"])
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)
    return {"p": p, "non_id_labels": non_id_labels, "fixed_solutions": fixed_solutions,
            "groups": groups, "diag": diag, "plus": plus, "kept": kept}


def ckpt_path_for_K(K):
    return os.path.join(CKPT_DIR, f"task27c_full_h4_folds_K{K}.json")


def run_submission_for_K(K, backend, prob):
    """Submits and retrieves ALL (fold, model) jobs for ONE K value,
    saving its OWN checkpoint immediately on completion -- deliberately
    NOT batched with the other K (a prior combined run was killed by an
    external process mid-retrieval with nothing saved; per-K checkpoints
    mean a future interruption loses at most one K's data, not both)."""
    ckpt_path = ckpt_path_for_K(K)
    if os.path.exists(ckpt_path):
        print(f"  K={K}: found existing checkpoint -> reusing, not resubmitting")
        with open(ckpt_path) as f:
            return json.load(f)

    import task2_fold_response_dataset as t2mod
    t2mod.FOLD_FACTORS = FOLD_FACTORS   # see module-level bug note: task2's own default is [1,3,5,9], missing fold=7

    cache = {}  # (slot, gate_name) -> {fold: [(group, qc), ...]}
    for name in prob["kept"]:
        angles = prob["fixed_solutions"][name]["angles"]
        for gate_name in ["ms", "zz"]:
            folded_circuits, meta = build_folded_measurement_circuits(angles, gate_name, prob["groups"])
            cache[(name, gate_name)] = folded_circuits

    all_jobs = {}  # (fold, model) -> (job, tags)
    t0 = time.time()
    for fold in FOLD_FACTORS:
        circuits_by_gate = {"ms": [], "zz": []}
        tags_by_gate = {"ms": [], "zz": []}
        for name in prob["kept"]:
            for gate_name in ["ms", "zz"]:
                for group, qc in cache[(name, gate_name)][fold]:
                    circuits_by_gate[gate_name].append(qc)
                    tags_by_gate[gate_name].append((name, tuple(group)))
        for model, gate_name in GATE_BY_MODEL.items():
            job = submit_job(circuits_by_gate[gate_name], backend, model, shots=SHOTS)
            all_jobs[(fold, model)] = (job, tags_by_gate[gate_name])
    t_submit = time.time() - t0
    print(f"  K={K}: all {len(all_jobs)} (fold, model) jobs submitted (non-blocking), {t_submit:.1f}s")

    t0 = time.time()
    all_counts = {}
    for i, (key, (job, tags)) in enumerate(all_jobs.items()):
        all_counts[f"{key[0]}|{key[1]}"] = get_counts_list(job)
        print(f"    K={K}: retrieved {i+1}/{len(all_jobs)}: fold={key[0]} model={key[1]}, {time.time()-t0:.1f}s elapsed")
    t_retrieve = time.time() - t0
    print(f"  K={K}: all {len(all_jobs)} jobs retrieved, {t_retrieve:.1f}s")

    ck = {
        "K": K, "fold_factors": FOLD_FACTORS, "shots": SHOTS,
        "tags": {f"{k[0]}|{k[1]}": [[t[0], list(t[1])] for t in v[1]] for k, v in all_jobs.items()},
        "counts": all_counts,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
    }
    with open(ckpt_path, "w") as f:
        json.dump(ck, f, indent=2)
    print(f"  K={K}: checkpoint saved -> {ckpt_path}")
    return ck


def main():
    print("\n" + "=" * 96)
    print("  task27c_full_h4_folds.py -- full H4 forged energy, native folds 1/3/5/7/9, concurrent real submission")
    print("=" * 96)

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    problems = {K: build_problem(K) for K in K_VALUES}
    for K, prob in problems.items():
        print(f"  K={K}: {len(prob['kept'])} kept circuits ({len(prob['diag'])} diag + {len(prob['plus'])} "
              f"'+'-pairs, down from {len(prob['p']['targets'])}) x {len(prob['groups'])} groups")

    cks = {}
    for K in K_VALUES:
        cks[K] = run_submission_for_K(K, backend, problems[K])

    # -- assemble raw Pauli values per (K, fold, model, slot, label), then reconstruct energy --
    print(f"\n  -- assembling energies (algebraic subspace-tomography cross-term derivation) --")
    energy_table = {}
    signal_death = {}
    for K in K_VALUES:
        prob = problems[K]
        ck = cks[K]
        p, non_id_labels, groups = prob["p"], prob["non_id_labels"], prob["groups"]
        diag, plus = prob["diag"], prob["plus"]
        group_idx = {}
        for gi, g in enumerate(groups):
            for l in g:
                group_idx[l] = gi

        energy_table[K] = {}
        for model in ["ideal", "aria-1", "forte-1"]:
            energy_table[K][model] = {}
            for fold in FOLD_FACTORS:
                key = f"{fold}|{model}"
                tags = ck["tags"][key]
                counts_list = ck["counts"][key]
                per_name = {}
                for (name, group), counts in zip(tags, counts_list):
                    per_name.setdefault(name, {}).setdefault(tuple(group), []).append(counts)

                errs = []
                for seed in range(N_SEEDS):
                    rng = np.random.default_rng(stable_seed("task27c", K, fold, model, seed))
                    m = {name: {} for name in diag + plus}
                    for name in diag + plus:
                        for group_t, counts_l in per_name[name].items():
                            resampled = bootstrap_counts(counts_l[0], SHOTS, rng)
                            for l in group_t:
                                m[name][l] = expectation_from_counts(resampled, l)
                    full = {name: dict(m[name]) for name in diag}
                    for n in range(K):
                        for mm in range(K):
                            if n >= mm:
                                continue
                            un, um, pl = f"u_{n}", f"u_{mm}", f"(u{n}+u{mm})"
                            full[pl] = dict(m[pl])
                            synth_minus = {}
                            for l in non_id_labels:
                                if l not in m[pl] or l not in full[un] or l not in full[um]:
                                    continue
                                cross = m[pl][l] - (full[un][l] + full[um][l]) / 2
                                synth_minus[l] = m[pl][l] - 2 * cross
                            full[f"(u{n}-u{mm})"] = synth_minus
                    _, err = energy_and_err(p, full, K)
                    errs.append(err)
                mean, std = float(np.mean(errs)), float(np.std(errs))
                energy_table[K][model][fold] = {"mean": mean, "std": std}
                print(f"    K={K} model={model} fold={fold}: {mean:.2f}+/-{std:.2f} kcal/mol")

            # -- where does the signal die? SNR = |mean_change_from_fold1| / std --
            base = energy_table[K][model][1]["mean"]
            for fold in FOLD_FACTORS:
                row = energy_table[K][model][fold]
                snr = abs(row["mean"] - base) / row["std"] if row["std"] > 1e-9 else float("inf")
                row["snr_vs_fold1"] = snr
            signal_death[f"{K}|{model}"] = {f: energy_table[K][model][f]["snr_vs_fold1"] for f in FOLD_FACTORS}

    print(f"\n  -- signal-to-noise vs fold (|mean change from fold=1| / std) --")
    for key, snrs in signal_death.items():
        dead_folds = [f for f, s in snrs.items() if s < 1.0]
        print(f"    {key}: SNR={snrs}  {'signal indistinguishable from noise at fold(s) ' + str(dead_folds) if dead_folds else 'signal remains meaningful at all tested folds'}")

    results = {"energy_table": energy_table, "signal_death_snr": signal_death,
               "shots": SHOTS, "checkpoints": {K: ckpt_path_for_K(K) for K in K_VALUES}}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
