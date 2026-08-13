#!/usr/bin/env python3
"""
task28b_ionq_optimizer.py — iteration 28, Task B. Apply qiskit-ionq's
TrappedIonOptimizerPlugin to the UNFOLDED native circuit, measured on the
RIGHT metric: raw noisy H4 energy at unchanged fidelity, not "gate count
improved."
============================================================================
Reuses `fixed_ansatz.native_optimized_gate_counts` UNCHANGED for the
gate-count/fidelity measurement (already verifies statevector-
preservation before trusting any count). Applied to ALL 36 K=6 targets,
both gate families. The user's own cited numbers ("1q 397->163 and
669->207, 2q 34->33") do not match this project's CURRENT K=6 circuit's
own gate counts and are NOT assumed to apply — this file measures fresh,
matching iteration 27's own precedent of catching an unverified number
(K=5's "0.17 kcal/mol").

RAW ENERGY, measured the RIGHT way (a real AerSimulator local proxy
cannot execute gpi/gpi2 instructions directly -- confirmed by a real
crash, "unknown instruction: gpi2" -- so this uses the SAME real,
free `ionq_simulator` submission this whole iteration is built on,
exactly as instructed ("RUN ideal/aria-1/forte-1 CONCURRENTLY on
ionq_simulator")): the "before" (un-optimized) raw energy is REUSED
directly from iteration 27 Task C's own checkpoint (fold=1, K=6, the
SAME 21-kept-circuit subspace-tomography design) — no resubmission
needed, it already exists on disk. Only the "after" (optimized) circuits
are submitted fresh, since the optimizer changes gate structure
per-target (N_2q now VARIES 4-11 across targets, not the constant 11 the
un-optimized circuit has — a real, disclosed structural change).

CRITICAL ORDERING, per explicit instruction: optimize the UNFOLDED
circuit, freeze it, THEN fold (Task 28D). This file does NOT fold
anything.

Run:
    python vqe/task28b_ionq_optimizer.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz, native_optimized_gate_counts
from native_stateprep import to_native, native_target
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from qiskit.quantum_info import Statevector

K = 6
GATE_FAMILIES = {"aria-1": "ms", "forte-1": "zz"}
SHOTS = 100_000
N_SEEDS = 8
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
BEFORE_CKPT_PATH = os.path.join(CKPT_DIR, "task27c_full_h4_folds_K6.json")
AFTER_CKPT_PATH = os.path.join(CKPT_DIR, "task28b_optimized_raw.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task28b_ionq_optimizer_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def measurement_circuits_for_optimized(optimized_qc, gate_name, groups):
    circuits = []
    for group in groups:
        combined = effrag_mod.combined_basis_label(group)
        basis_qc = native_basis_change(combined, gate_name)
        qc = optimized_qc.compose(basis_qc)
        qc.measure_all()
        circuits.append(qc)
    return circuits


def main():
    print("\n" + "=" * 96)
    print("  task28b_ionq_optimizer.py -- TrappedIonOptimizerPlugin on the RIGHT metric (raw energy)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)
    print(f"  K={K}: 36/36 targets converged; {len(kept)} kept circuits (subspace-tomography design, "
          f"same as iteration 27 Task C, so the 'before' number is directly reusable)")

    # -- gate-count / fidelity measurement, all 36 targets, both families --
    from qiskit.transpiler import PassManagerConfig
    from qiskit_ionq import TrappedIonOptimizerPlugin

    gate_stats = {}
    optimized_circuits_by_family = {}
    for family, gate_name in GATE_FAMILIES.items():
        print(f"\n  ===== {family} ({gate_name}) =====")
        before_n1q, after_n1q, before_n2q, after_n2q, sv_errs = [], [], [], [], []
        optimized_circuits = {}
        for name, sol in fixed_solutions.items():
            res = native_optimized_gate_counts(sol["angles"], gate_name)
            before_n1q.append(res["n1q_before_optimizer"])
            after_n1q.append(res["n1q_after_optimizer"])
            before_n2q.append(res["n2q_before_optimizer"])
            after_n2q.append(res["n2q_after_optimizer"])
            sv_errs.append(res["statevector_identical_err"])
            qc = build_ansatz(sol["angles"])
            native = to_native(qc, gate_name)
            tgt = native_target(qc.num_qubits, gate_name)
            pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt), optimization_level=3)
            optimized_circuits[name] = pm.run(native)
        optimized_circuits_by_family[family] = optimized_circuits

        print(f"    N_1q before: min={min(before_n1q)} max={max(before_n1q)} mean={np.mean(before_n1q):.1f}")
        print(f"    N_1q after:  min={min(after_n1q)} max={max(after_n1q)} mean={np.mean(after_n1q):.1f}")
        print(f"    N_2q before: {sorted(set(before_n2q))}   N_2q after: {sorted(set(after_n2q))} "
              f"(mean={np.mean(after_n2q):.2f})")
        print(f"    statevector-identical: worst err={max(sv_errs):.2e}, all <1e-9: {all(e < 1e-9 for e in sv_errs)}")
        gate_stats[family] = {
            "n1q_before": {"min": min(before_n1q), "max": max(before_n1q), "mean": float(np.mean(before_n1q))},
            "n1q_after": {"min": min(after_n1q), "max": max(after_n1q), "mean": float(np.mean(after_n1q))},
            "n2q_before": sorted(set(before_n2q)),
            "n2q_after": {"values": sorted(set(after_n2q)), "mean": float(np.mean(after_n2q))},
            "worst_statevector_err": max(sv_errs),
        }

    # -- REAL raw energy: "before" reused from iteration 27 Task C's checkpoint, "after" submitted fresh --
    with open(BEFORE_CKPT_PATH) as f:
        before_ck = json.load(f)

    if os.path.exists(AFTER_CKPT_PATH):
        print(f"\n  found existing 'after' checkpoint -> reusing, not resubmitting")
        with open(AFTER_CKPT_PATH) as f:
            after_ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"\n  connected, backend={backend.name}")
        all_jobs = {}
        t0 = time.time()
        for model, gate_name in GATE_FAMILIES.items():
            family = model  # GATE_FAMILIES keys are already "aria-1"/"forte-1" model names here
            circuits, tags = [], []
            for name in kept:
                optimized_qc = optimized_circuits_by_family[model][name]
                for group, qc in zip(groups, measurement_circuits_for_optimized(optimized_qc, gate_name, groups)):
                    circuits.append(qc)
                    tags.append((name, tuple(group)))
            job = submit_job(circuits, backend, model, shots=SHOTS)
            all_jobs[model] = (job, tags)
        # ideal control: use aria's (ms) optimized circuits (matches this project's own established convention)
        ideal_circuits, ideal_tags = [], []
        for name in kept:
            optimized_qc = optimized_circuits_by_family["aria-1"][name]
            for group, qc in zip(groups, measurement_circuits_for_optimized(optimized_qc, "ms", groups)):
                ideal_circuits.append(qc)
                ideal_tags.append((name, tuple(group)))
        all_jobs["ideal"] = (submit_job(ideal_circuits, backend, "ideal", shots=SHOTS), ideal_tags)
        t_submit = time.time() - t0
        print(f"  all {len(all_jobs)} model jobs submitted (non-blocking), {t_submit:.1f}s")

        t0 = time.time()
        after_counts = {}
        for i, (model, (job, tags)) in enumerate(all_jobs.items()):
            after_counts[model] = get_counts_list(job)
            print(f"    retrieved {i+1}/{len(all_jobs)}: model={model}, {time.time()-t0:.1f}s elapsed")
        t_retrieve = time.time() - t0
        after_ck = {
            "tags": {m: [[t[0], list(t[1])] for t in v[1]] for m, v in all_jobs.items()},
            "counts": after_counts, "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        }
        with open(AFTER_CKPT_PATH, "w") as f:
            json.dump(after_ck, f, indent=2)
        print(f"  checkpoint saved -> {AFTER_CKPT_PATH}")

    # -- assemble energies: BEFORE (from iteration 27's fold=1 K=6 checkpoint) and AFTER (this file) --
    group_idx = {}
    for gi, g in enumerate(groups):
        for l in g:
            group_idx[l] = gi

    def energy_from_ck(ck, key_for_model, model):
        tags = ck["tags"][key_for_model]
        counts_list = ck["counts"][key_for_model]
        per_name = {}
        for (name, group), counts in zip(tags, counts_list):
            per_name.setdefault(name, {}).setdefault(tuple(group), counts)
        errs = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task28b", model, seed))
            m = {name: {} for name in diag + plus}
            for name in diag + plus:
                for group_t, counts in per_name[name].items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
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
        return float(np.mean(errs)), float(np.std(errs))

    print(f"\n  -- RAW H4 energy: before (iteration 27 Task C, fold=1) vs after (this file's optimized circuits) --")
    raw_energy = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        before_key = f"1|{model}"
        before_mean, before_std = energy_from_ck(before_ck, before_key, model)
        after_mean, after_std = energy_from_ck(after_ck, model, model)
        raw_energy[model] = {"before": {"mean": before_mean, "std": before_std},
                              "after": {"mean": after_mean, "std": after_std}}
        print(f"    {model}: before={before_mean:.2f}+/-{before_std:.2f}  after={after_mean:.2f}+/-{after_std:.2f} kcal/mol")

    for family in ["aria-1", "forte-1"]:
        before_e = raw_energy[family]["before"]["mean"]
        after_e = raw_energy[family]["after"]["mean"]
        worst_sv = gate_stats[family]["worst_statevector_err"]
        success = after_e < before_e and worst_sv < 1e-9
        print(f"    {family}: SUCCESS (raw energy improved at unchanged fidelity) = {success} "
              f"({before_e:.2f} -> {after_e:.2f}, fidelity_err={worst_sv:.2e})")

    results = {"gate_stats": gate_stats, "raw_energy_kcal": raw_energy}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
