#!/usr/bin/env python3
"""
z2_tapered_fixed_ansatz_ionq.py — the hand-derived, FIXED-STRUCTURE
8-angle circuit for the Z2-tapered (3-qubit) register, measured for real,
concurrently, on IonQ's free ionq_simulator (ideal/aria-1/forte-1).
LOCAL BRANCH ONLY, not pushed.
============================================================================
WHY THIS CIRCUIT EXISTS: gate_structure_compare.py found the abstract
11-CX ansatz's real-hardware edge over the generic-StatePreparation
tapered circuit (34.98/43.03 vs 47.78/51.25 kcal/mol, aria-1/forte-1)
does NOT come from raw gate count (abstract has MORE of both 2q and 1q
gates) but from gate STRUCTURE: ~83% of the abstract ansatz's u3 gates
are fixed, special-angle (0, +-pi/2, +-pi) structural gates from its
Givens-rotation recipe, vs 0% for generic StatePreparation's arbitrary
output. z2_tapered_fixed_ansatz.py ported the SAME recipe (reference
prep + discriminator-qubit bridge + XXPlusYYGate Givens hops) to the
tapered 3-qubit register and found the naive 4-5 angle version (mirroring
the original 4-qubit ansatz's parameter count 1:1) converges for only
4/36 targets (worst error up to 0.51) -- NOT because the topology is
wrong, but because it's under-parameterized: the original 4-qubit
ansatz's 25 targets all share ONE fixed bit-complement pair (a physical
fact -- 2 electrons confined to a single Hamming-weight-2 sector), while
the tapered register's 36 targets are measured (directly, not assumed)
to spread across 1, 2, or 3 of the register's 3 available bit-complement
pairs depending on target -- tapering's Clifford transform does not
preserve Hamming-weight structure. Repeating the 3 Givens pairs to 7
hops (8 angles total) DOES converge for all 36/36 targets to machine
precision (worst 2.22e-16) with a CONSTANT gate count: 16 CX, 86 u3,
83% special-angle on u_0 -- genuinely closer to the abstract ansatz's
structure than generic StatePreparation is.

HONEST CAVEAT, stated before any real submission: 16 CX is MORE than
both baselines (abstract's 11, StatePreparation's mean 3.72) -- the
naive fidelity model f=(1-p2)^n2q*(1-p1)^n1q predicts this LOSES to
both (f_new=0.801 vs f_abstract=0.861 vs f_tapered_generic=0.954, at
aria-1's measured p2=0.01214, p1=p2/40). But that SAME naive model also
predicted tapered_generic (f=0.954) should beat abstract (f=0.861) on
real hardware, and it did NOT (47.78/51.25 lost to 34.98/43.03) -- so
the naive model is already known to mispredict this exact comparison
and cannot be trusted to rule this out in advance. Running for real
(free ionq_simulator, no real credits) is the only honest way to find
out. Reported whichever way it lands.

Run:
    python vqe/z2_tapered_fixed_ansatz_ionq.py --targets
    python vqe/z2_tapered_fixed_ansatz_ionq.py --assemble
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
from qiskit.circuit.library import XXPlusYYGate
from qiskit.quantum_info import Statevector
from qiskit.transpiler import CouplingMap
from qiskit import transpile

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "z2_tapered_fixed_ansatz_ionq_results.json")
FIT_PATH = os.path.join(os.path.dirname(__file__), "z2_tapered_fixed_ansatz_results.json")
CMAP3 = CouplingMap.from_full(3)

PAIRS = [(0, 1), (0, 2), (1, 2), (0, 1), (0, 2), (1, 2), (0, 1)]


def build_candidate8(angles):
    theta0 = angles[0]
    rest = angles[1:]
    qc = QuantumCircuit(3)
    qc.x(2)
    qc.ry(2 * theta0, 0)
    qc.cx(0, 1)
    qc.cx(0, 2)
    for th, (a, b) in zip(rest, PAIRS):
        qc.append(XXPlusYYGate(th, np.pi / 2), [a, b])
    return qc


def transpiled_fixed_ansatz_3q(angles):
    qc = build_candidate8(angles)
    return transpile(qc, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP3, optimization_level=0)


def verify_fit_exact(solutions, reduced_targets, tol=1e-9):
    worst = 0.0
    for name, sol in solutions.items():
        sv = np.asarray(Statevector.from_instruction(build_candidate8(sol["angles"])))
        target = np.asarray(reduced_targets[name])
        idx = int(np.argmax(np.abs(target)))
        phase = sv[idx] / target[idx] if abs(target[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv / phase - target)))
        worst = max(worst, err)
    assert worst < tol, f"fixed-ansatz circuit does not reproduce targets exactly: worst={worst:.2e}"
    return worst


def measurement_circuit_for_group(base, group):
    combined = effrag.combined_basis_label(group)
    qc = base.copy()
    basis_change(qc, combined)
    qc.measure_all()
    return qc


def phase_targets():
    print("\n" + "=" * 96)
    print("  z2_tapered_fixed_ansatz_ionq.py --targets  (8-angle fixed-structure circuit, RAW)")
    print("=" * 96)

    problem = build_reduced_problem()
    reduced_targets = problem["reduced_targets"]
    reduced_label_map = problem["reduced_label_map"]

    with open(FIT_PATH) as f:
        fit = json.load(f)
    solutions = fit["solutions"]
    assert fit["n_ok"] == len(reduced_targets), "fit did not converge for all targets -- refusing to submit for real"
    worst = verify_fit_exact(solutions, reduced_targets)
    print(f"  fixed-ansatz circuit verified exact against all {len(reduced_targets)} targets: worst={worst:.2e}")
    print(f"  gate structure: {fit['n2q_gates']} CX (constant={fit['n2q_constant']}), "
          f"{fit['n1q_gates_mean']:.0f} u3 (constant={fit['n1q_constant']}), "
          f"{fit['u0_n_special_angle_u3']}/{fit['u0_n_u3']} special-angle on u_0")

    unique_reduced_labels = sorted(set(rl for rl, _ in reduced_label_map.values()))
    groups = effrag.group_labels_qubit_wise(unique_reduced_labels)
    target_names = sorted(reduced_targets.keys())
    print(f"  36 targets, {len(unique_reduced_labels)} unique reduced (3-qubit) labels, {len(groups)} groups")

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"  connected, backend={backend.name}")

    circuits, tags = [], []
    for name in target_names:
        base = transpiled_fixed_ansatz_3q(solutions[name]["angles"])
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
    save_ckpt("z2_tapered_fixed_ansatz_targets", out)
    print(f"\n  --targets phase complete\n")
    return out


def assemble():
    print("\n" + "=" * 96)
    print("  z2_tapered_fixed_ansatz_ionq.py --assemble")
    print("=" * 96)

    ck = load_ckpt("z2_tapered_fixed_ansatz_targets")
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
            rng = np.random.default_rng(stable_seed("z2_tapered_fixed_ansatz", model, seed))
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
        report[model] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs)), "errs": errs}
        print(f"    {model}: {report[model]['mean_kcal']:.3f} +/- {report[model]['std_kcal']:.3f} kcal/mol (8 seeds)")

    print(f"\n  -- comparison to prior real-hardware findings --")
    print(f"    abstract 11-gate ansatz (raw): aria-1=34.98, forte-1=43.03 kcal/mol (11 CX, mostly special-angle)")
    print(f"    Z2-tapered, generic StatePreparation (raw): aria-1=47.78, forte-1=51.25 kcal/mol (mean 3.94 CX, generic-angle)")
    print(f"    native-optimized (raw): aria-1=93.73, forte-1=91.43 kcal/mol")
    print(f"    THIS run (Z2-tapered, 8-angle FIXED structure, raw): aria-1={report['aria-1']['mean_kcal']:.2f}, "
          f"forte-1={report['forte-1']['mean_kcal']:.2f} kcal/mol (16 CX constant, 86 u3 constant, 83% special-angle)")
    ideal_ok = report["ideal"]["mean_kcal"] < 5.0
    print(f"\n  ideal correctness control: {report['ideal']['mean_kcal']:.3f} kcal/mol "
          f"({'PASS' if ideal_ok else 'FAIL -- pipeline bug, not noise'})")

    results = {"K": K, "n_seeds": N_SEEDS, "shots": SHOTS, "wall_clock": ck["wall_clock"], "report": report,
               "ideal_correctness_control_pass": bool(ideal_ok),
               "comparison": {
                   "abstract_raw": {"aria-1": 34.98, "forte-1": 43.03},
                   "tapered_stateprep_raw": {"aria-1": 47.78, "forte-1": 51.25},
                   "native_optimized_raw": {"aria-1": 93.73, "forte-1": 91.43},
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
