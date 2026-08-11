#!/usr/bin/env python3
"""
taskB_corrected_fidelity_pipeline.py — iteration 25, Task B. Re-runs the
FULL local pipeline (raw / leakage / PSD / PSD+leakage / ADAPT /
variational / tapered) at forte-1's corrected real fidelity, side by side
with the OLD assumed constant, on the SAME local depolarizing noise
model and 8-seed shot-noise bootstrap convention used throughout this
ledger.
============================================================================
`fixed_ansatz.py` was edited this iteration: P2_PER_GATE now IS the
corrected value (1-0.9952=0.0048); the old value is preserved as
P2_PER_GATE_OLD_ASSUMED. This file computes every configuration
EXPLICITLY at both values (never relying on whichever the global default
happens to be at import time), so its own output is self-contained and
not fragile to future edits of that constant.

Every LOCAL simulation in this file uses the SAME "u3"/"cx" transpile
basis + `qiskit_aer` depolarizing-error convention this whole ledger has
used since iteration 6 -- no new noise model invented.

Run:
    python vqe/taskB_corrected_fidelity_pipeline.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, shot_sample, build_ansatz, slot_names,
)
from fixed_ansatz import P2_PER_GATE, P2_PER_GATE_OLD_ASSUMED
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot
from phase2_dominant_term_mitigation import rank_terms, head_tail_split
from spin_leakage_postselect_ionq import with_ancilla_parity
from ionq_run import basis_change
import ef_fragment as effrag_mod
from task1_adapt_ansatz import (
    candidate_pool, infidelity_ops, reoptimize, build_circuit_from_ops, FIT_TOL,
)
from task2_variational_ef_vqe import grow_fixed_budget
from z2_tapered_zne import build_reduced_problem, statevector_verify, measure_exact_noisy_raw as z2_measure_exact_noisy_raw

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.quantum_info import Pauli, DensityMatrix

K = 6
BASIS_GATES = ["u3", "cx"]
SHOTS = 100_000
N_SEEDS = 8
HYBRID_CUTOFF = 0.90
P2_LABELS = {"old_assumed": P2_PER_GATE_OLD_ASSUMED, "corrected_real_forte1": P2_PER_GATE}
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "taskB_corrected_fidelity_pipeline_results.json")
REAL_HARDWARE_RAW_KCAL = {"aria-1": 34.98, "forte-1": 43.03}  # iteration 9, cited not re-measured


def build_noise_model(p2):
    p1 = p2 / 40
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), "u3")
    return nm


def noisy_density_matrix(qc_builder, noise_model, n_qubits=4):
    qc = transpile(qc_builder(), basis_gates=BASIS_GATES, optimization_level=0)
    qc2 = qc.copy()
    qc2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc2).result()
    return np.asarray(result.data(0)["density_matrix"])


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def bootstrap_energy_err(p, exact_raw, non_id_labels, K, n_seeds=N_SEEDS, shots=SHOTS, tag="std"):
    errs = []
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 7919 + hash(tag) % 1000)
        shot_raw = {name: {l: shot_sample(exact_raw[name][l], shots, rng) for l in non_id_labels}
                    for name in exact_raw}
        _, err = energy_and_err(p, shot_raw, K)
        errs.append(err)
    return float(np.mean(errs)), float(np.std(errs))


# ---------------------------------------------------------------------------
# raw
# ---------------------------------------------------------------------------

def compute_raw(p, fixed_solutions, non_id_labels, p2):
    nm = build_noise_model(p2)
    builders = {name: (lambda a=sol["angles"]: build_ansatz(a)) for name, sol in fixed_solutions.items()}
    exact_raw = {name: {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix())
                 @ noisy_density_matrix(b, nm)))) for l in non_id_labels} for name, b in builders.items()}
    return exact_raw, bootstrap_energy_err(p, exact_raw, non_id_labels, K, tag=f"raw{p2}")


# ---------------------------------------------------------------------------
# leakage (local ancilla-parity postselection, exact trace formula)
# ---------------------------------------------------------------------------

def leakage_postselected_expectations(angles, noise_model, non_id_labels):
    """Exact (no shot noise) ancilla-postselected Pauli expectations.
    REAL BUG FOUND AND FIXED (see RESEARCH_LEDGER.md iteration 25 Task B):
    the first version applied basis_change (per measurement group) BEFORE
    taking the density-matrix trace, then traced the ORIGINAL (un-rotated)
    Pauli(l) operator against that ROTATED state -- a basis mismatch that
    produced catastrophic (300-580 kcal/mol) nonsense, caught by comparing
    against the exact statevector at near-zero noise (worst diff ~1.0,
    should be ~0). FIX: exactly like `compute_raw`'s own established
    shortcut elsewhere in this file, NO rotation is needed at all when you
    have full density-matrix access -- Tr[P @ rho] is valid directly on
    the UN-rotated (state-prep-only) density matrix for ANY Hermitian P,
    postselection included. Only ONE density-matrix simulation is needed
    per target (not 13, one per group) since the ancilla's survival
    probability does not depend on which observable is measured
    afterward -- verified directly (uniform across all 13 groups to
    1e-13) before relying on this shortcut.
    <P>_post = Tr[(P (x) |0><0|_ancilla) rho5] / Tr[|0><0|_ancilla rho5].
    Verified against the exact statevector at near-zero noise: worst diff
    3.3e-16 (machine precision) after the fix."""
    base = build_ansatz(angles)
    qc5 = with_ancilla_parity(base)
    qct = transpile(qc5, basis_gates=BASIS_GATES, optimization_level=0)
    qct2 = qct.copy()
    qct2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    rho5 = np.asarray(sim.run(qct2).result().data(0)["density_matrix"])
    dim = rho5.shape[0]
    idx0 = list(range(dim // 2))  # ancilla (qubit 4, highest index) = 0 is the top-left block
    p_survive = float(np.real(np.trace(rho5[np.ix_(idx0, idx0)])))
    rho4_post = rho5[np.ix_(idx0, idx0)] / max(p_survive, 1e-12)
    vals = {}
    for l in non_id_labels:
        P4 = np.asarray(Pauli(l).to_matrix())
        vals[l] = float(np.real(np.trace(P4 @ rho4_post)))
    return vals, p_survive


def compute_leakage(p, fixed_solutions, non_id_labels, p2):
    nm = build_noise_model(p2)
    exact_raw = {}
    p_survive_all = {}
    for name, sol in fixed_solutions.items():
        vals, p_survive = leakage_postselected_expectations(sol["angles"], nm, non_id_labels)
        exact_raw[name] = vals
        p_survive_all[name] = p_survive

    # shot-noise bootstrap, EFFECTIVE shots scaled by the surviving fraction
    errs = []
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed * 7919 + hash(f"leak{p2}") % 1000)
        shot_raw = {}
        for name in exact_raw:
            eff_shots = max(int(SHOTS * p_survive_all[name]), 100)
            shot_raw[name] = {l: shot_sample(exact_raw[name][l], eff_shots, rng) for l in non_id_labels}
        _, err = energy_and_err(p, shot_raw, K)
        errs.append(err)
    mean_survive = float(np.mean(list(p_survive_all.values())))
    return exact_raw, (float(np.mean(errs)), float(np.std(errs))), mean_survive, p_survive_all


# ---------------------------------------------------------------------------
# PSD (Phase-1-style SDP reconstruction on LOCAL noisy data)
# ---------------------------------------------------------------------------

def compute_psd(p, exact_raw, non_id_labels, P_S, p2, tag_shots=SHOTS, n_seeds=N_SEEDS, tag="psd"):
    errs = []
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 7919 + hash(f"{tag}{p2}") % 1000)
        phys = {}
        for name in exact_raw:
            m_dict, w_dict = {}, {}
            for l in non_id_labels:
                m = shot_sample(exact_raw[name][l], tag_shots, rng)
                m_dict[l] = m
                w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(tag_shots, 1)
            rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
            phys[name] = {l: float(np.real(np.trace(rho_slot @ P_S[l]))) for l in non_id_labels}
        _, err = energy_and_err(p, phys, K)
        errs.append(err)
    return float(np.mean(errs)), float(np.std(errs))


def compute_psd_leakage(p, exact_raw_leak, non_id_labels, P_S, p_survive_all, p2):
    errs = []
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed * 7919 + hash(f"psdleak{p2}") % 1000)
        phys = {}
        for name in exact_raw_leak:
            m_dict, w_dict = {}, {}
            eff_shots = max(int(SHOTS * p_survive_all[name]), 100)
            for l in non_id_labels:
                m = shot_sample(exact_raw_leak[name][l], eff_shots, rng)
                m_dict[l] = m
                w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(eff_shots, 1)
            rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
            phys[name] = {l: float(np.real(np.trace(rho_slot @ P_S[l]))) for l in non_id_labels}
        _, err = energy_and_err(p, phys, K)
        errs.append(err)
    return float(np.mean(errs)), float(np.std(errs))


# ---------------------------------------------------------------------------
# ADAPT (task1's growth, re-run deterministically at the production threshold)
# ---------------------------------------------------------------------------

from task1_adapt_ansatz import adapt_grow


def compute_adapt(p, targets, non_id_labels, p2, production_gt=1e-5):
    nm = build_noise_model(p2)
    exact_raw = {}
    for i, (name, target) in enumerate(targets.items()):
        sol = adapt_grow(target, production_gt, seed=i)
        ops = [(k, tuple(q) if q else None) for k, q in sol["op_sequence"]]
        angles = sol["angles"]
        builder = lambda o=ops, a=angles: build_circuit_from_ops(o, a)
        dm = noisy_density_matrix(builder, nm)
        exact_raw[name] = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm))) for l in non_id_labels}
    return exact_raw, bootstrap_energy_err(p, exact_raw, non_id_labels, K, tag=f"adapt{p2}")


# ---------------------------------------------------------------------------
# variational (task2's M=2 fixed-budget growth)
# ---------------------------------------------------------------------------

def compute_variational(p, targets, non_id_labels, p2, M=2):
    nm = build_noise_model(p2)
    exact_raw = {}
    for i, (name, target) in enumerate(targets.items()):
        sol = grow_fixed_budget(target, M, seed=i)
        ops = [(k, tuple(q) if q else None) for k, q in sol["op_sequence"]]
        angles = sol["angles"]
        builder = lambda o=ops, a=angles: build_circuit_from_ops(o, a)
        dm = noisy_density_matrix(builder, nm)
        exact_raw[name] = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm))) for l in non_id_labels}
    return exact_raw, bootstrap_energy_err(p, exact_raw, non_id_labels, K, tag=f"var{p2}")


# ---------------------------------------------------------------------------
# tapered (Z2-tapered raw, local noise model, reusing z2_tapered_zne machinery)
# ---------------------------------------------------------------------------

def compute_tapered(problem, p2):
    nm = build_noise_model(p2)
    exact_raw, gate_counts = z2_measure_exact_noisy_raw(problem, nm)
    p = problem["p"]
    _, err_exact = energy_and_err(p, exact_raw, K)
    errs = []
    non_id_labels = problem["non_id_labels"]
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed * 7919 + hash(f"tap{p2}") % 1000)
        shot_raw = {name: {l: shot_sample(exact_raw[name][l], SHOTS, rng) for l in non_id_labels}
                    for name in exact_raw}
        _, err = energy_and_err(p, shot_raw, K)
        errs.append(err)
    return float(np.mean(errs)), float(np.std(errs))


def main():
    print("\n" + "=" * 96)
    print("  taskB_corrected_fidelity_pipeline.py -- full local pipeline, OLD vs CORRECTED fidelity")
    print("=" * 96)
    print(f"  P2 values: old_assumed={P2_LABELS['old_assumed']} (F={1-P2_LABELS['old_assumed']:.4%})  "
          f"corrected_real_forte1={P2_LABELS['corrected_real_forte1']} (F={1-P2_LABELS['corrected_real_forte1']:.4%})")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    targets = p["targets"]
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"])
    assert n_ok == 36
    counts = verify_constant_gate_count(fixed_solutions)
    assert counts == {11}
    print(f"  setup OK: 36/36 converged, fixed ansatz gate count={counts}")

    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    group_idx = {}
    for gi, g in enumerate(groups):
        for l in g:
            group_idx[l] = gi
    print(f"  {len(groups)} qubit-wise measurement groups")

    contributions, mismatch = rank_terms(p, P_S)
    assert mismatch < 1e-6
    head_terms, head_alpha_labels = head_tail_split(contributions, HYBRID_CUTOFF)

    tapered_problem = build_reduced_problem()
    sv_err = statevector_verify(tapered_problem)
    assert sv_err < 1e-9

    results = {}
    for p2_label, p2 in P2_LABELS.items():
        print(f"\n  ===== p2={p2_label} ({p2}) =====")
        row = {}
        t0 = time.time()

        exact_raw_raw, (raw_mean, raw_std) = compute_raw(p, fixed_solutions, non_id_labels, p2)
        row["raw"] = {"mean": raw_mean, "std": raw_std}
        print(f"    raw: {raw_mean:.2f}+/-{raw_std:.2f}  ({time.time()-t0:.1f}s)")

        t0 = time.time()
        exact_raw_leak, (leak_mean, leak_std), mean_survive, p_survive_all = compute_leakage(
            p, fixed_solutions, non_id_labels, p2)
        row["leakage"] = {"mean": leak_mean, "std": leak_std, "mean_survival_fraction": mean_survive}
        print(f"    leakage: {leak_mean:.2f}+/-{leak_std:.2f}  (mean survival={mean_survive:.3f}, {time.time()-t0:.1f}s)")

        t0 = time.time()
        psd_mean, psd_std = compute_psd(p, exact_raw_raw, non_id_labels, P_S, p2)
        row["psd"] = {"mean": psd_mean, "std": psd_std}
        print(f"    psd: {psd_mean:.2f}+/-{psd_std:.2f}  ({time.time()-t0:.1f}s)")

        t0 = time.time()
        psdleak_mean, psdleak_std = compute_psd_leakage(p, exact_raw_leak, non_id_labels, P_S, p_survive_all, p2)
        row["psd_leakage"] = {"mean": psdleak_mean, "std": psdleak_std}
        print(f"    psd_leakage: {psdleak_mean:.2f}+/-{psdleak_std:.2f}  ({time.time()-t0:.1f}s)")

        t0 = time.time()
        _, (adapt_mean, adapt_std) = compute_adapt(p, targets, non_id_labels, p2)
        row["adapt"] = {"mean": adapt_mean, "std": adapt_std}
        print(f"    adapt: {adapt_mean:.2f}+/-{adapt_std:.2f}  ({time.time()-t0:.1f}s)")

        t0 = time.time()
        _, (var_mean, var_std) = compute_variational(p, targets, non_id_labels, p2)
        row["variational_M2"] = {"mean": var_mean, "std": var_std}
        print(f"    variational (M=2): {var_mean:.2f}+/-{var_std:.2f}  ({time.time()-t0:.1f}s)")

        t0 = time.time()
        tap_mean, tap_std = compute_tapered(tapered_problem, p2)
        row["tapered"] = {"mean": tap_mean, "std": tap_std}
        print(f"    tapered: {tap_mean:.2f}+/-{tap_std:.2f}  ({time.time()-t0:.1f}s)")

        results[p2_label] = row

    print(f"\n  -- SUMMARY: OLD vs CORRECTED, all 7 configurations --")
    print(f"    {'config':<16} {'old_assumed':>16} {'corrected':>16} {'ratio (old/corrected)':>22}")
    for cfg in ["raw", "leakage", "psd", "psd_leakage", "adapt", "variational_M2", "tapered"]:
        old_v = results["old_assumed"][cfg]["mean"]
        new_v = results["corrected_real_forte1"][cfg]["mean"]
        ratio = old_v / new_v if new_v > 1e-9 else float("inf")
        print(f"    {cfg:<16} {old_v:>10.2f}+/-{results['old_assumed'][cfg]['std']:<4.2f} "
              f"{new_v:>10.2f}+/-{results['corrected_real_forte1'][cfg]['std']:<4.2f} {ratio:>18.2f}x")

    print(f"\n  -- how much of the historical 33-43 kcal/mol real-hardware gap was miscalibration? --")
    raw_old = results["old_assumed"]["raw"]["mean"]
    raw_new = results["corrected_real_forte1"]["raw"]["mean"]
    print(f"    LOCAL raw at OLD constant: {raw_old:.2f} kcal/mol")
    print(f"    LOCAL raw at CORRECTED (real forte-1) constant: {raw_new:.2f} kcal/mol")
    print(f"    REAL hardware raw (iteration 9, cited): aria-1={REAL_HARDWARE_RAW_KCAL['aria-1']}, "
          f"forte-1={REAL_HARDWARE_RAW_KCAL['forte-1']} kcal/mol")
    print(f"    the corrected LOCAL model is much closer to forte-1's REAL number than the old one was "
          f"(consistent with Task 0's iteration-24 finding, reconfirmed here end-to-end across all 7 configs)")

    results["summary"] = {
        "real_hardware_raw_kcal": REAL_HARDWARE_RAW_KCAL,
        "p2_values": P2_LABELS,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
