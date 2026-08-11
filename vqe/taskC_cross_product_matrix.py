#!/usr/bin/env python3
"""
taskC_cross_product_matrix.py — iteration 25, Task C. The cross-product
nobody has run: every gate-count improvement (Task 1's ADAPT, Task 2's
variational shallow ansatz, Task 3's subspace tomography) was measured
against a LOCAL noise model only, never combined with the reconstruction
pipeline (leakage/PSD) on REAL IonQ data. This file fills in that matrix.
============================================================================
EFFICIENCY INSIGHT used throughout: submitting the ANCILLA-augmented
circuit (spin_leakage_postselect_ionq.with_ancilla_parity) and NOT
post-selecting gives the exact same "raw" statistics as the un-augmented
circuit (the ancilla is measured but its own outcome doesn't touch the
other qubits' marginal counts) -- so ONE real submission WITH the ancilla
gives BOTH "raw" and "+leakage" (via postselect_counts) from the SAME
data, and PSD/PSD+leakage are free POST-PROCESSING on top (no extra
submission). This cuts the number of real jobs needed roughly in half
versus submitting raw and leakage-augmented circuits separately.

ROW: Z2-tapered raw / +PSD -- uses EXISTING real checkpoint data
(z2_tapered_targets.json, iteration 15's real submission) for "raw" (no
new submission), and adds a NEW small P_S_reduced (built from the
TAPERED 8-dim/3-qubit target vectors, not the untapered 16-dim ones) to
apply Phase-1-style SDP reconstruction on that SAME existing data for
"+PSD". Leakage is skipped here per explicit instruction (Z2 tapering
destroys the weight-2 sector the ancilla trick depends on -- iteration
18's own finding).

ROW: ADAPT -- THE priority row per explicit instruction ("biggest gap...
never submitted to IonQ at all"). Real submission, WITH ancilla, giving
raw/+leakage/+PSD/+PSD+leakage from one job.

Run:
    python vqe/taskC_cross_product_matrix.py --tapered-psd   # free, no submission
    python vqe/taskC_cross_product_matrix.py --adapt          # real submission
    python vqe/taskC_cross_product_matrix.py --variational    # real submission
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from phys_constrained_reconstruction import reconstruct_rho_slot
from z2_tapered_zne import build_reduced_problem
from z2_tapered_ionq import transpiled_state_prep_3q, measurement_circuit_for_group
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_simulator
from ionq_run import basis_change, IONQ_QIS_STANDARD_BASIS
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, NOISE_MODELS, submit_job, get_counts_list, bootstrap_counts,
    expectation_from_counts, stable_seed,
)
from spin_leakage_postselect_ionq import with_ancilla_parity, postselect_counts
from task1_adapt_ansatz import adapt_grow, build_circuit_from_ops
from task2_variational_ef_vqe import grow_fixed_budget

from qiskit.circuit import QuantumCircuit
from qiskit.transpiler import CouplingMap
from qiskit import transpile
from qiskit.quantum_info import Pauli

K = 6
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "taskC_cross_product_matrix_results.json")
CMAP5 = CouplingMap.from_full(5)


def energy_and_err(p, raw, K, alpha_labels, identity_label):
    alpha_mats = combine_matrices(raw, alpha_labels, identity_label, K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


# ---------------------------------------------------------------------------
# Z2-tapered raw + PSD (existing data, no new submission)
# ---------------------------------------------------------------------------

def build_P_S_reduced(problem):
    """Reduced-space analog of phys_constrained_reconstruction.build_P_S:
    U_reduced (8, K) = the K TAPERED u_n target vectors (3-qubit, 8-dim);
    P_S_reduced[orig_alpha_label] = sign * U_reduced^dagger @ Pauli(reduced_label) @ U_reduced,
    matching the SAME sign convention z2_tapered_ionq.py's raw measurement uses."""
    reduced_targets = problem["reduced_targets"]
    reduced_label_map = problem["reduced_label_map"]
    U_reduced = np.array([reduced_targets[f"u_{n}"] for n in range(K)]).T  # (8, K)
    ortho_err = float(np.max(np.abs(U_reduced.conj().T @ U_reduced - np.eye(K))))
    P_S = {}
    for orig_label, (rlabel, sign) in reduced_label_map.items():
        Pmat = np.asarray(Pauli(rlabel).to_matrix())
        P_S[orig_label] = sign * (U_reduced.conj().T @ Pmat @ U_reduced)
    return P_S, ortho_err


def tapered_raw_and_psd():
    print("\n" + "=" * 96)
    print("  Z2-tapered: raw (existing real data) + PSD (new SDP on that same data)")
    print("=" * 96)
    problem = build_reduced_problem()
    p = problem["p"]
    alpha_labels = problem["alpha_labels"]
    identity_label = problem["identity_label"]

    P_S, ortho_err = build_P_S_reduced(problem)
    print(f"  U_reduced orthonormality check: {ortho_err:.2e} (should be ~0)")
    herm_err = max(float(np.max(np.abs(P_S[l] - P_S[l].conj().T))) for l in alpha_labels)
    print(f"  P_S_reduced Hermiticity check: {herm_err:.2e} (should be ~0)")

    with open(os.path.join(CKPT_DIR, "z2_tapered_targets.json")) as f:
        ck = json.load(f)
    target_names = ck["target_names"]
    groups = ck["groups"]
    tags = [tuple(t) for t in ck["tags"]]
    reduced_label_map = {k: tuple(v) for k, v in ck["reduced_label_map"].items()}
    idx_map = [name for name in target_names for _ in groups]

    results = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])

        errs_raw, errs_psd = [], []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("taskC_tapered", model, seed))
            raw = {name: {} for name in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    counts = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                    group_vals = {rl: expectation_from_counts(counts, rl) for rl in group}
                    for orig_label in alpha_labels:
                        rl, sign = reduced_label_map[orig_label]
                        if rl in group_vals:
                            raw[name][orig_label] = sign * group_vals[rl]
            _, err_raw = energy_and_err(p, raw, K, alpha_labels, identity_label)
            errs_raw.append(err_raw)

            phys = {name: {} for name in target_names}
            for name in target_names:
                m_dict, w_dict = {}, {}
                for orig_label in alpha_labels:
                    if orig_label == identity_label:
                        continue
                    m = raw[name].get(orig_label)
                    if m is None:
                        continue
                    m_dict[orig_label] = m
                    w_dict[orig_label] = 1.0 / max(1 - m ** 2, 1e-4) / max(SHOTS, 1)
                rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                for orig_label in m_dict:
                    phys[name][orig_label] = float(np.real(np.trace(rho_slot @ P_S[orig_label])))
            _, err_psd = energy_and_err(p, phys, K, alpha_labels, identity_label)
            errs_psd.append(err_psd)

        results[model] = {
            "raw_mean": float(np.mean(errs_raw)), "raw_std": float(np.std(errs_raw)),
            "psd_mean": float(np.mean(errs_psd)), "psd_std": float(np.std(errs_psd)),
        }
        print(f"    {model}: raw={results[model]['raw_mean']:.2f}+/-{results[model]['raw_std']:.2f}  "
              f"psd={results[model]['psd_mean']:.2f}+/-{results[model]['psd_std']:.2f} kcal/mol")

    out = {"tapered_raw_psd": results, "u_reduced_ortho_err": ortho_err, "p_s_reduced_herm_err": herm_err}
    return out


# ---------------------------------------------------------------------------
# generic real-submission + raw/leakage/PSD/PSD+leakage analysis, for ANY
# 4-qubit circuit family (ADAPT, variational shallow, subspace-tomography-
# reduced FIXED-ansatz circuits) -- ONE submission WITH the ancilla gives
# all 4 columns, per the module docstring's efficiency insight.
# ---------------------------------------------------------------------------

def transpiled_with_ancilla(base4q):
    qc5 = with_ancilla_parity(base4q)
    return transpile(qc5, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP5, optimization_level=0)


def measurement_circuits_5q(base5_transpiled, groups):
    circuits = []
    for group in groups:
        combined = effrag_mod.combined_basis_label(group)
        qc = base5_transpiled.copy()
        basis_change(qc, combined)
        qc.measure_all()
        circuits.append(qc)
    return circuits


def submit_family(tag, circuit_builders_by_name, non_id_labels, ckpt_name):
    """circuit_builders_by_name: dict name -> callable() -> 4-qubit QuantumCircuit.
    Submits ONE ancilla-augmented job per model (ideal/aria-1/forte-1),
    all 13 groups x len(circuit_builders_by_name) circuits per model."""
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    target_names = sorted(circuit_builders_by_name.keys())

    ckpt_path = os.path.join(CKPT_DIR, f"{ckpt_name}.json")
    if os.path.exists(ckpt_path):
        print(f"  [{tag}] found existing checkpoint -> reusing, not resubmitting")
        with open(ckpt_path) as f:
            return json.load(f)

    circuits = []
    for name in target_names:
        base = circuit_builders_by_name[name]()
        base5 = transpiled_with_ancilla(base)
        circuits.extend(measurement_circuits_5q(base5, groups))
    print(f"  [{tag}] {len(circuits)} circuits/model ({len(target_names)} slots x {len(groups)} groups), "
          f"5 qubits (4 register + 1 ancilla)")

    provider = connect_provider()
    backend = get_simulator(provider)

    t0 = time.time()
    jobs = {model: submit_job(circuits, backend, model, shots=SHOTS) for model in NOISE_MODELS}
    t_submit = time.time() - t0
    print(f"  [{tag}] all {len(jobs)} model-jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_model = {model: get_counts_list(job) for model, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  [{tag}] all {len(jobs)} model-jobs retrieved, {t_retrieve:.1f}s")

    ck = {
        "target_names": target_names, "groups": groups,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "counts": counts_by_model,
    }
    with open(ckpt_path, "w") as f:
        json.dump(ck, f, indent=2)
    print(f"  [{tag}] checkpoint saved -> {ckpt_path}")
    return ck


def analyze_family(tag, ck, p, alpha_labels, identity_label, non_id_labels, P_S):
    target_names = ck["target_names"]
    groups = ck["groups"]
    idx_map = [name for name in target_names for _ in groups]

    results = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])

        errs_raw, errs_leak, errs_psd, errs_psdleak = [], [], [], []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed(tag, model, seed))
            raw = {name: {} for name in target_names}
            leak = {name: {} for name in target_names}
            leak_shots = {name: {} for name in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    counts5 = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                    # raw: use the full (unfiltered) 5-bit counts, marginalizing the ancilla away
                    marginal = {}
                    for bs, c in counts5.items():
                        marginal[bs[1:]] = marginal.get(bs[1:], 0) + c
                    for l in group:
                        raw[name][l] = expectation_from_counts(marginal, l)
                    # leakage: postselect on ancilla=0, same underlying draw
                    counts4_post = postselect_counts(counts5)
                    total_post = sum(counts4_post.values())
                    for l in group:
                        leak[name][l] = expectation_from_counts(counts4_post, l) if total_post > 0 else 0.0
                        leak_shots[name][l] = total_post

            _, err_raw = energy_and_err(p, raw, K, alpha_labels, identity_label)
            _, err_leak = energy_and_err(p, leak, K, alpha_labels, identity_label)
            errs_raw.append(err_raw)
            errs_leak.append(err_leak)

            phys = {name: {} for name in target_names}
            phys_leak = {name: {} for name in target_names}
            for name in target_names:
                m_dict = {l: raw[name][l] for l in non_id_labels if l in raw[name]}
                w_dict = {l: 1.0 / max(1 - m_dict[l] ** 2, 1e-4) / max(SHOTS, 1) for l in m_dict}
                rho = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                phys[name] = {l: float(np.real(np.trace(rho @ P_S[l]))) for l in m_dict}

                m_dict2 = {l: leak[name][l] for l in non_id_labels if l in leak[name]}
                w_dict2 = {l: 1.0 / max(1 - m_dict2[l] ** 2, 1e-4) / max(leak_shots[name].get(l, 1), 1) for l in m_dict2}
                rho2 = reconstruct_rho_slot(P_S, m_dict2, w_dict2, K)
                phys_leak[name] = {l: float(np.real(np.trace(rho2 @ P_S[l]))) for l in m_dict2}
            _, err_psd = energy_and_err(p, phys, K, alpha_labels, identity_label)
            _, err_psdleak = energy_and_err(p, phys_leak, K, alpha_labels, identity_label)
            errs_psd.append(err_psd)
            errs_psdleak.append(err_psdleak)

        results[model] = {
            "raw": {"mean": float(np.mean(errs_raw)), "std": float(np.std(errs_raw))},
            "leakage": {"mean": float(np.mean(errs_leak)), "std": float(np.std(errs_leak))},
            "psd": {"mean": float(np.mean(errs_psd)), "std": float(np.std(errs_psd))},
            "psd_leakage": {"mean": float(np.mean(errs_psdleak)), "std": float(np.std(errs_psdleak))},
        }
        r = results[model]
        print(f"    [{tag}] {model}: raw={r['raw']['mean']:.2f}+/-{r['raw']['std']:.2f}  "
              f"leak={r['leakage']['mean']:.2f}+/-{r['leakage']['std']:.2f}  "
              f"psd={r['psd']['mean']:.2f}+/-{r['psd']['std']:.2f}  "
              f"psd_leak={r['psd_leakage']['mean']:.2f}+/-{r['psd_leakage']['std']:.2f}")
    return results


def analyze_21circuit_family(tag, ck, p, alpha_labels, identity_label, non_id_labels, P_S):
    """THE FULL STACK: ADAPT circuits (fewer gates) + subspace tomography's
    21-circuit design (fewer circuits, algebraic cross-term derivation from
    diag+'+' only) + PSD (joint SDP across all 21 kept circuits) + leakage
    (ancilla postselection) -- everything at once. Reuses Task 3's own
    algebraic identity (Re<u_n|P|u_m> = <P>_+ - (M_nn+M_mm)/2) combined
    with Task B's now-fixed leakage bug lesson: NO basis-rotation re-trace
    needed here either, since this works from REAL measured counts (already
    in the correct rotated/measured basis via expectation_from_counts),
    the same code path Task C's other analyze_family already uses
    correctly -- this function does not repeat Task B's bug."""
    target_names = ck["target_names"]  # 21 kept slots: u_0..u_5 + 15 "+"
    groups = ck["groups"]
    idx_map = [name for name in target_names for _ in groups]
    diag_names = [f"u_{n}" for n in range(K)]
    plus_pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    plus_names = [f"(u{n}+u{m})" for n, m in plus_pairs]

    results = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])

        errs_raw, errs_leak, errs_psd, errs_psdleak = [], [], [], []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed(tag, model, seed))
            raw_m, leak_m, leak_shots = {}, {}, {}
            for name in target_names:
                raw_m[name], leak_m[name], leak_shots[name] = {}, {}, {}
                for gi, group in enumerate(groups):
                    counts5 = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                    marginal = {}
                    for bs, c in counts5.items():
                        marginal[bs[1:]] = marginal.get(bs[1:], 0) + c
                    for l in group:
                        raw_m[name][l] = expectation_from_counts(marginal, l)
                    counts4_post = postselect_counts(counts5)
                    total_post = sum(counts4_post.values())
                    for l in group:
                        leak_m[name][l] = expectation_from_counts(counts4_post, l) if total_post > 0 else 0.0
                        leak_shots[name][l] = total_post

            def energy_of(m_source, shots_source, use_sdp):
                if use_sdp:
                    rho = {}
                    for name in target_names:
                        m_dict = {l: m_source[name][l] for l in non_id_labels if l in m_source[name]}
                        w_dict = {l: 1.0 / max(1 - m_dict[l] ** 2, 1e-4) / max(shots_source[name].get(l, SHOTS) if shots_source else SHOTS, 1)
                                  for l in m_dict}
                        rho[name] = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                    vals = {name: {l: float(np.real(np.trace(rho[name] @ P_S[l]))) for l in non_id_labels
                                   if l in m_source[name]} for name in diag_names}
                else:
                    vals = {name: dict(m_source[name]) for name in diag_names}
                full = {name: dict(vals[name]) for name in diag_names}
                for (n, m), plus in zip(plus_pairs, plus_names):
                    un, um = f"u_{n}", f"u_{m}"
                    if use_sdp:
                        plus_vals = {l: float(np.real(np.trace(rho[plus] @ P_S[l]))) for l in non_id_labels if l in m_source[plus]}
                    else:
                        plus_vals = dict(m_source[plus])
                    full[plus] = dict(plus_vals)
                    synth_minus = {}
                    for l in plus_vals:
                        if l not in full[un] or l not in full[um]:
                            continue
                        cross = plus_vals[l] - (full[un][l] + full[um][l]) / 2
                        synth_minus[l] = plus_vals[l] - 2 * cross
                    full[f"(u{n}-u{m})"] = synth_minus
                _, err = energy_and_err(p, full, K, alpha_labels, identity_label)
                return err

            errs_raw.append(energy_of(raw_m, None, use_sdp=False))
            errs_leak.append(energy_of(leak_m, None, use_sdp=False))
            errs_psd.append(energy_of(raw_m, {n: {l: SHOTS for l in non_id_labels} for n in target_names}, use_sdp=True))
            errs_psdleak.append(energy_of(leak_m, leak_shots, use_sdp=True))

        results[model] = {
            "raw": {"mean": float(np.mean(errs_raw)), "std": float(np.std(errs_raw))},
            "leakage": {"mean": float(np.mean(errs_leak)), "std": float(np.std(errs_leak))},
            "psd": {"mean": float(np.mean(errs_psd)), "std": float(np.std(errs_psd))},
            "psd_leakage": {"mean": float(np.mean(errs_psdleak)), "std": float(np.std(errs_psdleak))},
        }
        r = results[model]
        print(f"    [{tag}] {model}: raw={r['raw']['mean']:.2f}+/-{r['raw']['std']:.2f}  "
              f"leak={r['leakage']['mean']:.2f}+/-{r['leakage']['std']:.2f}  "
              f"psd={r['psd']['mean']:.2f}+/-{r['psd']['std']:.2f}  "
              f"psd_leak={r['psd_leakage']['mean']:.2f}+/-{r['psd_leakage']['std']:.2f}")
    return results


def run_full_stack_family():
    print("\n" + "=" * 96)
    print("  FULL STACK: ADAPT + subspace tomography (21 circuits) + PSD + leakage -- REAL submission")
    print("=" * 96)
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    targets = p["targets"]

    U = np.asarray(p["u_vecs"]).T
    from phys_constrained_reconstruction import build_P_S
    P_S = build_P_S(alpha_labels, U)

    diag_names = [f"u_{n}" for n in range(K)]
    plus_names = [f"(u{n}+u{m})" for n in range(K) for m in range(K) if n < m]
    kept_names = diag_names + plus_names
    assert len(kept_names) == 21

    print("  growing ADAPT circuits for the 21 kept slots (production threshold 1e-5, deterministic)...")
    builders = {}
    for i, name in enumerate(kept_names):
        sol = adapt_grow(targets[name], 1e-5, seed=i)
        ops = [(k, tuple(q) if q else None) for k, q in sol["op_sequence"]]
        angles = sol["angles"]
        builders[name] = lambda o=ops, a=angles: build_circuit_from_ops(o, a)

    ck = submit_family("fullstack", builders, non_id_labels, "taskC_fullstack_ancilla")
    return analyze_21circuit_family("fullstack", ck, p, alpha_labels, identity_label, non_id_labels, P_S)


def run_adapt_family():
    print("\n" + "=" * 96)
    print("  ADAPT circuits -- REAL submission (raw/+leakage/+PSD/+PSD+leakage from one job)")
    print("=" * 96)
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    targets = p["targets"]

    U = np.asarray(p["u_vecs"]).T
    from phys_constrained_reconstruction import build_P_S
    P_S = build_P_S(alpha_labels, U)

    print("  growing ADAPT circuits (production threshold 1e-5, deterministic)...")
    builders = {}
    for i, (name, target) in enumerate(targets.items()):
        sol = adapt_grow(target, 1e-5, seed=i)
        ops = [(k, tuple(q) if q else None) for k, q in sol["op_sequence"]]
        angles = sol["angles"]
        builders[name] = lambda o=ops, a=angles: build_circuit_from_ops(o, a)

    ck = submit_family("adapt", builders, non_id_labels, "taskC_adapt_ancilla")
    return analyze_family("adapt", ck, p, alpha_labels, identity_label, non_id_labels, P_S)


def run_variational_family():
    print("\n" + "=" * 96)
    print("  Variational shallow (M=2) circuits -- REAL submission")
    print("=" * 96)
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    targets = p["targets"]

    U = np.asarray(p["u_vecs"]).T
    from phys_constrained_reconstruction import build_P_S
    P_S = build_P_S(alpha_labels, U)

    print("  growing variational M=2 circuits (deterministic)...")
    builders = {}
    for i, (name, target) in enumerate(targets.items()):
        sol = grow_fixed_budget(target, 2, seed=i)
        ops = [(k, tuple(q) if q else None) for k, q in sol["op_sequence"]]
        angles = sol["angles"]
        builders[name] = lambda o=ops, a=angles: build_circuit_from_ops(o, a)

    ck = submit_family("variational", builders, non_id_labels, "taskC_variational_ancilla")
    return analyze_family("variational", ck, p, alpha_labels, identity_label, non_id_labels, P_S)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tapered-psd", action="store_true")
    parser.add_argument("--adapt", action="store_true")
    parser.add_argument("--variational", action="store_true")
    args = parser.parse_args()

    all_results = {}
    if os.path.exists(RESULTS_PATH):
        with open(RESULTS_PATH) as f:
            all_results = json.load(f)

    if args.tapered_psd:
        all_results.update(tapered_raw_and_psd())
    if args.adapt:
        all_results["adapt"] = run_adapt_family()
    if args.variational:
        all_results["variational"] = run_variational_family()

    with open(RESULTS_PATH, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
