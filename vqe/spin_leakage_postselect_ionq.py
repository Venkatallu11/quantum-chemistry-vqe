#!/usr/bin/env python3
"""
spin_leakage_postselect_ionq.py — Task D (spin/particle-number leakage
projection), completed: does discarding real-hardware shots that leak
out of the physical Hamming-weight-2 sector actually improve the
untapered abstract-ansatz's real energy accuracy, across ALL 13
measurement groups (not just the one Z-basis group retroactive analysis
of iteration 9's existing data could reach)?
============================================================================
MOTIVATION (iteration 17's leading unexplained lead): the untapered
4-qubit alpha register is a 2-electron/4-orbital system, so every
PHYSICAL state has EXACTLY Hamming weight 2 in the computational (Z)
basis. Real IonQ noise can flip this. Direct analysis of iteration 9's
ALREADY-COLLECTED real hardware counts (the one measurement group that
happens to stay in the Z basis, no rotation) found REAL, non-negligible
leakage: 0% on ideal, 4.74% on aria-1, 5.12% on forte-1 of shots have
the wrong Hamming weight. Post-selecting on weight==2 for that one
group's 3 Pauli labels (IIZZ, ZIIZ, ZZII) cut RMS error vs the exact
value by ~2.2-2.4x (aria-1: 0.051->0.023, forte-1: 0.052->0.021) --
using ONLY data already on disk, no new submission.

THE LIMITATION that one-group analysis has: Hamming weight of a RAW
measured bitstring is only physically meaningful when the measurement
happened in the Z (computational) basis. The other 12 of 13 measurement
groups apply basis-rotation gates (H for X, Sdg+H for Y) before
measuring -- a raw post-rotation bitstring's weight tells us nothing
about particle number. To extend post-selection to ALL groups, this file
adds ONE ancilla qubit and 4 CNOTs (register qubit -> ancilla) BEFORE any
basis-rotation gates, so the ancilla ends up holding the PARITY of the
register's pre-rotation Z-basis weight, entirely independent of whatever
basis rotation is subsequently applied to the register qubits themselves
(gates on disjoint qubits commute; the rotation never touches the
ancilla). Detects ODD-weight leakage (the dominant single-bit-flip error
class) on every group, not just the one Z-basis group.

VERIFIED EXACTLY before building any measurement circuits (see
verify_ancilla_scheme()): the ancilla-CNOT step changes NOTHING about
the marginal state on the original 4 register qubits (partial_trace over
the ancilla matches the un-augmented circuit's density matrix to
1.5e-36), and for a noiseless physical (weight-2) state the ancilla reads
0 with probability exactly 1 -- confirmed by direct Statevector
computation, not assumed from the CNOT's textbook behavior.

Run:
    python vqe/spin_leakage_postselect_ionq.py --targets
    python vqe/spin_leakage_postselect_ionq.py --assemble
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from fixed_ansatz import build_ansatz
import rank6_symmetry_vd as r
import ef_fragment as effrag
from qforge import combine_matrices, energy_from_alpha_matrices, setup_fragment
from ionq_backend import connect_provider, get_simulator
from ionq_run import basis_change, IONQ_QIS_STANDARD_BASIS
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, NOISE_MODELS, submit_job, get_counts_list, bootstrap_counts,
    expectation_from_counts, stable_seed, save_ckpt, load_ckpt,
)

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace
from qiskit.transpiler import CouplingMap
from qiskit import transpile

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "spin_leakage_postselect_ionq_results.json")
CMAP5 = CouplingMap.from_full(5)
K = 6


def with_ancilla_parity(base_qc):
    """4-qubit state-prep circuit -> 5-qubit circuit with a parity-check
    ancilla appended (qubit 4). CNOTs must precede any basis-rotation
    gates, so this is called on the RAW state-prep circuit, before
    basis_change()."""
    qc = QuantumCircuit(5)
    qc.compose(base_qc, qubits=[0, 1, 2, 3], inplace=True)
    qc.cx(0, 4)
    qc.cx(1, 4)
    qc.cx(2, 4)
    qc.cx(3, 4)
    return qc


def verify_ancilla_scheme(angles):
    from fixed_ansatz import build_ansatz as _build
    base = _build(angles)
    qc5 = with_ancilla_parity(base)
    sv5 = Statevector.from_instruction(qc5)
    sv4 = Statevector.from_instruction(base)
    rho4_direct = np.outer(np.asarray(sv4), np.asarray(sv4).conj())
    rho4_marginal = np.asarray(partial_trace(sv5, [4]))
    marginal_err = float(np.max(np.abs(rho4_direct - rho4_marginal)))
    probs = sv5.probabilities_dict()
    p_ancilla1 = sum(v for k, v in probs.items() if k[0] == "1")
    return marginal_err, p_ancilla1


def transpiled_ansatz_with_ancilla(angles):
    base = build_ansatz(angles)
    qc5 = with_ancilla_parity(base)
    return transpile(qc5, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP5, optimization_level=0)


def measurement_circuits_for_groups(base5_circuit, groups):
    """basis_change uses len(label)=4 internally, so it only ever
    touches qubits 0-3 -- qubit 4 (ancilla) is left untouched regardless
    of circuit width, verified by construction (see module docstring)."""
    circuits = []
    for group in groups:
        combined = effrag.combined_basis_label(group)
        qc = base5_circuit.copy()
        basis_change(qc, combined)
        qc.measure_all()
        circuits.append(qc)
    return circuits


def postselect_counts(counts, keep_ancilla="0"):
    """Ancilla is qubit 4, the highest-index qubit -> leftmost character
    in qiskit's standard measure_all() bitstring (bitstring[-1]=qubit0,
    so bitstring[0]=qubit4 for a 5-qubit string). Filters to shots where
    the pre-rotation register parity was even (weight in {0,2,4} --
    catches the dominant single-bit-flip / odd-weight leakage class)."""
    return {bs[1:]: c for bs, c in counts.items() if bs[0] == keep_ancilla}


def phase_targets():
    print("\n" + "=" * 96)
    print("  spin_leakage_postselect_ionq.py --targets  (ancilla parity-check, RAW)")
    print("=" * 96)

    p = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    assert n_ok == 36, f"fit did not converge for all 36 targets: {n_ok}/36"
    print(f"  36/36 targets converged, worst={worst:.2e}")

    marginal_err, p_anc1 = verify_ancilla_scheme(solutions["u_0"]["angles"])
    print(f"  ancilla scheme verified: marginal-state err={marginal_err:.2e}, "
          f"P(ancilla=1|noiseless,physical)={p_anc1:.2e}")
    assert marginal_err < 1e-9 and p_anc1 < 1e-9, "ancilla scheme does not verify exactly -- refusing to submit"

    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    groups = effrag.group_labels_qubit_wise(alpha_labels)
    target_names = sorted(solutions.keys())
    print(f"  {len(target_names)} targets x {len(groups)} groups = {len(target_names) * len(groups)} circuits/model")

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"  connected, backend={backend.name}")

    circuits, idx_map = [], []
    for name in target_names:
        base5 = transpiled_ansatz_with_ancilla(solutions[name]["angles"])
        for gi, group in enumerate(groups):
            qc = base5.copy()
            combined = effrag.combined_basis_label(groups[gi])
            basis_change(qc, combined)
            qc.measure_all()
            circuits.append(qc)
            idx_map.append(name)

    t0 = time.time()
    jobs = {model: submit_job(circuits, backend, model, shots=SHOTS) for model in NOISE_MODELS}
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_model = {model: get_counts_list(job) for model, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "exact_energy": p["exact_energy"], "noiseless_energy": p["noiseless_numpy"],
        "alpha_labels": alpha_labels, "identity_label": identity_label,
        "target_names": target_names, "groups": groups,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "counts": {model: counts_by_model[model] for model in NOISE_MODELS},
        "idx_map": idx_map,
    }
    save_ckpt("spin_leakage_targets", out)
    print(f"\n  --targets phase complete\n")
    return out


def energy_from_counts_scheme(counts_flat, idx_map, target_names, groups, alpha_labels, identity_label,
                                p, model, postselect, rng):
    per_name = {name: [] for name in target_names}
    for i, name in enumerate(idx_map):
        per_name[name].append(counts_flat[i])

    raw = {name: {} for name in target_names}
    for name in target_names:
        for gi, group in enumerate(groups):
            counts = bootstrap_counts(per_name[name][gi], SHOTS, rng)
            if postselect:
                counts = postselect_counts(counts)
            vals = {l: expectation_from_counts(counts, l) for l in group}
            raw[name].update(vals)

    mats = combine_matrices(raw, alpha_labels, identity_label, K)
    E, err = energy_from_alpha_matrices(mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                         exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return err["err_vs_exact_kcal"]


def assemble():
    print("\n" + "=" * 96)
    print("  spin_leakage_postselect_ionq.py --assemble")
    print("=" * 96)

    ck = load_ckpt("spin_leakage_targets")
    if ck is None:
        print("  missing checkpoint -- run --targets first")
        return None

    p = setup_fragment([0, 1, 2, 3], 4, 1.0, K)
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    target_names = ck["target_names"]
    groups = ck["groups"]
    idx_map = ck["idx_map"]

    report = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]

        # real, hardware-measured leakage fraction across all circuits this run
        total, odd = 0, 0
        for c in counts_flat:
            for bs, n in c.items():
                total += n
                if bs[0] == "1":
                    odd += n
        leakage_frac = odd / total if total else float("nan")

        errs_raw, errs_post = [], []
        for seed in range(N_SEEDS):
            rng_raw = np.random.default_rng(stable_seed("spin_leakage_raw", model, seed))
            errs_raw.append(energy_from_counts_scheme(counts_flat, idx_map, target_names, groups,
                                                        alpha_labels, identity_label, p, model,
                                                        postselect=False, rng=rng_raw))
            rng_post = np.random.default_rng(stable_seed("spin_leakage_post", model, seed))
            errs_post.append(energy_from_counts_scheme(counts_flat, idx_map, target_names, groups,
                                                         alpha_labels, identity_label, p, model,
                                                         postselect=True, rng=rng_post))
        report[model] = {
            "leakage_frac": leakage_frac,
            "raw_mean_kcal": float(np.mean(errs_raw)), "raw_std_kcal": float(np.std(errs_raw)),
            "post_mean_kcal": float(np.mean(errs_post)), "post_std_kcal": float(np.std(errs_post)),
        }
        print(f"    {model}: leakage_frac(odd-weight)={leakage_frac:.4f}  "
              f"RAW={report[model]['raw_mean_kcal']:.3f}+/-{report[model]['raw_std_kcal']:.3f}  "
              f"POST-SELECTED={report[model]['post_mean_kcal']:.3f}+/-{report[model]['post_std_kcal']:.3f} kcal/mol")

    print(f"\n  -- comparison --")
    print(f"    iteration 9 baseline (no ancilla, raw): aria-1=34.98, forte-1=43.03 kcal/mol")
    print(f"    THIS run, ancilla overhead, RAW: aria-1={report['aria-1']['raw_mean_kcal']:.2f}, "
          f"forte-1={report['forte-1']['raw_mean_kcal']:.2f} kcal/mol")
    print(f"    THIS run, POST-SELECTED (odd-weight discarded): aria-1={report['aria-1']['post_mean_kcal']:.2f}, "
          f"forte-1={report['forte-1']['post_mean_kcal']:.2f} kcal/mol")
    ideal_ok = report["ideal"]["raw_mean_kcal"] < 5.0
    print(f"\n  ideal correctness control (raw): {report['ideal']['raw_mean_kcal']:.3f} kcal/mol "
          f"({'PASS' if ideal_ok else 'FAIL -- pipeline bug, not noise'})")

    results = {"n_seeds": N_SEEDS, "shots": SHOTS, "wall_clock": ck["wall_clock"], "report": report,
               "ideal_correctness_control_pass": bool(ideal_ok),
               "comparison": {"iteration9_no_ancilla_raw": {"aria-1": 34.98, "forte-1": 43.03}}}
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
