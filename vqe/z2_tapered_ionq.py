#!/usr/bin/env python3
"""
z2_tapered_ionq.py — the Z2-tapered (3-qubit) H4 forged circuit, RAW (no
ZNE -- iteration 14 found ZNE does not converge on this circuit either),
measured for real, concurrently, on IonQ's free ionq_simulator
(ideal/aria-1/forte-1). LOCAL BRANCH ONLY, not pushed.
============================================================================
Reuses z2_tapered_zne.py's verified tapering pipeline unchanged (tapering
commutes exactly with beta_signs(), 0.0 diff; every one of the 37 alpha
Pauli labels reduces to a genuine 3-qubit Pauli with a real +-1 sign,
verified by direct matrix comparison; StatePreparation matches the exact
tapered target to 3.31e-14). NEW here: the 37 (with duplicates) reduced
labels collapse to 27 UNIQUE 3-qubit Paulis, grouping into just 9
qubit-wise-commuting measurement bases -- down from 13 for the untapered
register, a real efficiency gain on top of the gate-count reduction
already found.

Run:
    python vqe/z2_tapered_ionq.py --targets
    python vqe/z2_tapered_ionq.py --assemble
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from z2_tapered_zne import build_reduced_problem
from qforge import HARTREE_TO_KCAL_MOL, combine_matrices, energy_from_alpha_matrices
import ef_fragment as effrag
from ionq_backend import connect_provider, get_simulator
from ionq_run import basis_change, pauli_expectation, IONQ_QIS_STANDARD_BASIS
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, NOISE_MODELS, submit_job, get_counts_list, bootstrap_counts,
    expectation_from_counts, stable_seed, save_ckpt, load_ckpt,
)

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit.transpiler import CouplingMap
from qiskit import transpile

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "z2_tapered_ionq_results.json")
CMAP3 = CouplingMap.from_full(3)


def transpiled_state_prep_3q(vec):
    qc = QuantumCircuit(3)
    qc.append(StatePreparation(np.asarray(vec)), range(3))
    return transpile(qc, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP3, optimization_level=0)


def measurement_circuit_for_group(base, group):
    combined = effrag.combined_basis_label(group)
    qc = base.copy()
    basis_change(qc, combined)
    qc.measure_all()
    return qc


def phase_targets():
    print("\n" + "=" * 96)
    print("  z2_tapered_ionq.py --targets  (RAW, no ZNE)")
    print("=" * 96)

    problem = build_reduced_problem()
    reduced_targets = problem["reduced_targets"]
    reduced_label_map = problem["reduced_label_map"]
    unique_reduced_labels = sorted(set(rl for rl, _ in reduced_label_map.values()))
    groups = effrag.group_labels_qubit_wise(unique_reduced_labels)
    target_names = sorted(reduced_targets.keys())
    print(f"  36 targets, {len(unique_reduced_labels)} unique reduced (3-qubit) labels, {len(groups)} groups")

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"  connected, backend={backend.name}")

    circuits, tags = [], []
    for name in target_names:
        base = transpiled_state_prep_3q(reduced_targets[name])
        for gi, group in enumerate(groups):
            circuits.append(measurement_circuit_for_group(base, group))
            tags.append((name, gi))
    print(f"  {len(circuits)} circuits/model")

    t0 = time.time()
    jobs = {model: submit_job(circuits, backend, model, shots=SHOTS) for model in NOISE_MODELS}
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_model = {model: get_counts_list(job) for model, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "exact_energy": problem["p"]["exact_energy"], "noiseless_energy": problem["p"]["noiseless_energy"],
        "alpha_labels": problem["alpha_labels"], "identity_label": problem["identity_label"],
        "reduced_label_map": {k: list(v) for k, v in reduced_label_map.items()},
        "target_names": target_names, "groups": groups,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "counts": {model: counts_by_model[model] for model in NOISE_MODELS},
        "tags": [list(t) for t in tags],
    }
    save_ckpt("z2_tapered_targets", out)
    print(f"\n  --targets phase complete\n")
    return out


def assemble():
    print("\n" + "=" * 96)
    print("  z2_tapered_ionq.py --assemble")
    print("=" * 96)

    ck = load_ckpt("z2_tapered_targets")
    if ck is None:
        print("  missing checkpoint -- run --targets first")
        return None

    problem = build_reduced_problem()
    p = problem["p"]
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    reduced_label_map = {k: tuple(v) for k, v in ck["reduced_label_map"].items()}
    target_names = ck["target_names"]
    groups = ck["groups"]
    tags = [tuple(t) for t in ck["tags"]]
    K = 6

    idx_map = []
    for name in target_names:
        for _ in groups:
            idx_map.append(name)

    report = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])

        errs = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("z2_tapered", model, seed))
            raw = {name: {} for name in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    counts = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                    # measure every reduced label present in this group once, then fan out to
                    # every original alpha_label that maps onto it (several can share one group)
                    group_vals = {rl: expectation_from_counts(counts, rl) for rl in group}
                    for orig_label in alpha_labels:
                        rl, sign = reduced_label_map[orig_label]
                        if rl in group_vals:
                            raw[name][orig_label] = sign * group_vals[rl]
            alpha_mats = combine_matrices(raw, alpha_labels, identity_label, K)
            E, err = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                 exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs.append(err["err_vs_exact_kcal"])
        report[model] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs)), "errs": errs}
        print(f"    {model}: {report[model]['mean_kcal']:.3f} +/- {report[model]['std_kcal']:.3f} kcal/mol (8 seeds)")

    print(f"\n  -- comparison to prior real-hardware findings --")
    print(f"    iteration 9 (abstract 11-gate ansatz, raw): aria-1=34.98, forte-1=43.03 kcal/mol")
    print(f"    iteration 12 (native-optimized, TrappedIonOptimizerPlugin, raw): aria-1=93.73, forte-1=91.43 kcal/mol")
    print(f"    THIS run (Z2-tapered, 3-qubit, mean 3.94 CX, raw): aria-1={report['aria-1']['mean_kcal']:.2f}, "
          f"forte-1={report['forte-1']['mean_kcal']:.2f} kcal/mol")
    ideal_ok = report["ideal"]["mean_kcal"] < 5.0
    print(f"\n  ideal correctness control: {report['ideal']['mean_kcal']:.3f} kcal/mol "
          f"({'PASS' if ideal_ok else 'FAIL -- pipeline bug, not noise'})")

    results = {"K": K, "n_seeds": N_SEEDS, "shots": SHOTS, "wall_clock": ck["wall_clock"], "report": report,
               "ideal_correctness_control_pass": bool(ideal_ok),
               "comparison": {
                   "iteration9_abstract_raw": {"aria-1": 34.98, "forte-1": 43.03},
                   "iteration12_native_optimized_raw": {"aria-1": 93.73, "forte-1": 91.43},
               }}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    args = parser.parse_args()
    if args.targets:
        phase_targets()
    elif args.assemble:
        assemble()
    else:
        parser.error("pass --targets or --assemble")


if __name__ == "__main__":
    main()
