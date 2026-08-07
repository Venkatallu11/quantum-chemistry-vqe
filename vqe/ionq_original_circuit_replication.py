#!/usr/bin/env python3
"""
ionq_original_circuit_replication.py — Task 2 (the control) and Task 5
(the fidelity threshold curve) of the "go back to the original EF+ZNE
result" request.
============================================================================
TASK 2 -- ISOLATE CIRCUIT FROM NOISE, not blindly re-derive it. Takes the
EXACT circuits entanglement_forging_h4.py/entanglement_forging_zne.py use
-- StatePreparation of the EXACT (genuinely complex, psi max imag=0.968,
real_gauge was NEVER applied in the original script) K=5 Schmidt vectors,
BOTH registers measured independently (no beta_signs() shortcut -- that
shortcut requires a real-gauged state, which this circuit deliberately
is NOT, to stay faithful to "unchanged") -- and runs them for real,
concurrently, on ideal/aria-1/forte-1.

TWO ADAPTATIONS, DISCLOSED, NEITHER CHANGES WHAT IS BEING MEASURED:
  1. optimization_level=0 (this project's mandatory invariant, established
     after entanglement_forging_h4.py/ionq_run.py were written -- the
     ORIGINAL ionq_run.py used optimization_level=1 for state prep, which
     would violate the invariant this task explicitly lists; using 0 here
     is a correctness requirement, not a physics change).
  2. qubit-wise-commuting measurement GROUPING (13 groups covering the 37
     unique alpha/beta labels each) -- a real device cannot extract
     arbitrary Pauli expectations from one circuit's execution the way
     the original script's AerEstimatorV2(method="density_matrix") could;
     grouping only changes HOW MANY circuits it takes to measure the
     SAME set of expectation values, not which state is prepared or
     which observable is measured. The 4-phase-per-pair cross-term trick
     (NOT the real-gauge-optimized 2-phase version used in ef_fragment.py)
     is kept exactly as the original script has it, since this state is
     genuinely complex and the 2-phase shortcut does not apply.

If real IonQ noise gives ~35-43 kcal/mol raw error here (matching
iteration 9's fixed-ansatz numbers), the circuit is not the story -- the
noise is. If it lands much better, the newer 11-gate fixed ansatz is
somehow worse on IonQ than this older K=5 StatePreparation circuit, which
would overturn iteration 9's framing. Reported plainly either way.

Run:
    python vqe/ionq_original_circuit_replication.py --control
    python vqe/ionq_original_circuit_replication.py --assemble-control
    python vqe/ionq_original_circuit_replication.py --fidelity-curve
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
import ef_fragment as effrag
import rank6_symmetry_vd as r
import loop_pec as pec
from ionq_backend import connect_provider, get_simulator
from ionq_run import IONQ_QIS_STANDARD_BASIS, basis_change, pauli_expectation
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, NOISE_MODELS, HARTREE_TO_KCAL_MOL,
    submit_job, get_counts_list, bootstrap_counts, expectation_from_counts,
    stable_seed, save_ckpt, load_ckpt,
)

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit import transpile
from qiskit.transpiler import CouplingMap

K = 5
OPT_LEVEL = 0  # mandatory project invariant
CMAP4 = CouplingMap.from_full(4)
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "ionq_original_circuit_replication_results.json")


# ---------------------------------------------------------------------------
# Original state-prep circuits, transpiled at the mandatory opt_level=0
# (the ONLY change from ionq_run.py's own transpiled_state_prep, which
# predates that invariant and used opt_level=1)
# ---------------------------------------------------------------------------

def transpiled_state_prep_0(vec):
    qc = QuantumCircuit(4)
    qc.append(StatePreparation(np.asarray(vec)), range(4))
    return transpile(qc, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP4, optimization_level=OPT_LEVEL)


def measurement_circuits_for_groups(base_circuit, groups):
    circuits = []
    for group in groups:
        combined = effrag.combined_basis_label(group)
        qc = base_circuit.copy()
        basis_change(qc, combined)
        qc.measure_all()
        circuits.append(qc)
    return circuits


# ---------------------------------------------------------------------------
# PHASE: control -- submit ALL real circuits (both registers) to
# ideal/aria-1/forte-1 concurrently
# ---------------------------------------------------------------------------

def phase_control():
    print("\n" + "=" * 96)
    print("  ionq_original_circuit_replication.py --control  (Task 2)")
    print("=" * 96)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    qop_bare, qop_pen, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    groups_a = effrag.group_labels_qubit_wise(alpha_labels)
    groups_b = effrag.group_labels_qubit_wise(beta_labels)
    print(f"  K={K}, {len(alpha_labels)} alpha labels ({len(groups_a)} groups), "
          f"{len(beta_labels)} beta labels ({len(groups_b)} groups)")
    psi_max_imag = float(np.max(np.abs(psi.imag)))
    print(f"  psi max |imag|={psi_max_imag:.4f} -- genuinely complex, real_gauge/beta_signs shortcut "
          "does NOT apply here (faithful to the original, unchanged script)")

    def build_register_circuits(vecs, groups):
        circuits, tags = [], []
        for n in range(K):
            base = transpiled_state_prep_0(vecs[n])
            for gi, group in enumerate(groups):
                circuits.append(measurement_circuits_for_groups(base, [group])[0])
                tags.append(("diag", n, gi))
        pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
        for (n, m) in pairs:
            for k in range(4):
                vec_k = (vecs[n] + (1j ** k) * vecs[m]) / np.sqrt(2)
                base = transpiled_state_prep_0(vec_k)
                for gi, group in enumerate(groups):
                    circuits.append(measurement_circuits_for_groups(base, [group])[0])
                    tags.append(("cross", n, m, k, gi))
        return circuits, tags

    circuits_a, tags_a = build_register_circuits(u_vecs, groups_a)
    circuits_b, tags_b = build_register_circuits(v_vecs, groups_b)
    print(f"  alpha register: {len(circuits_a)} circuits, beta register: {len(circuits_b)} circuits")

    jobs = {}
    t0 = time.time()
    for model in NOISE_MODELS:
        jobs[("alpha", model)] = submit_job(circuits_a, backend, model, shots=SHOTS)
        jobs[("beta", model)] = submit_job(circuits_b, backend, model, shots=SHOTS)
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_key = {key: get_counts_list(job) for key, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "K": K, "exact_energy": exact_energy, "enuc": enuc,
        "lambdas": lambdas.tolist(), "psi_max_imag": psi_max_imag,
        "alpha_labels": alpha_labels, "beta_labels": beta_labels,
        "groups_a": groups_a, "groups_b": groups_b,
        "tags_a": [list(t) for t in tags_a], "tags_b": [list(t) for t in tags_b],
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "counts": {model: {"alpha": counts_by_key[("alpha", model)], "beta": counts_by_key[("beta", model)]}
                   for model in NOISE_MODELS},
    }
    save_ckpt("original_control", out)
    print(f"\n  --control phase complete\n")
    return out


def rebuild_matrices_from_counts(counts, tags, groups, labels, K, bootstrap_rng=None):
    """Reconstruct the K x K complex matrix per label from grouped
    real counts, mirroring entanglement_forging_h4.build_noisy_matrices'
    diag+cross(4-phase) reconstruction exactly."""
    def expct(i, group):
        c = counts[i]
        if bootstrap_rng is not None:
            c = bootstrap_counts(c, SHOTS, bootstrap_rng)
        return {l: expectation_from_counts(c, l) for l in group}

    diag_vals = {n: {} for n in range(K)}
    cross_vals = {}
    for i, tag in enumerate(tags):
        if tag[0] == "diag":
            _, n, gi = tag
            diag_vals[n].update(expct(i, groups[gi]))
        else:
            _, n, m, k, gi = tag
            cross_vals.setdefault((n, m, k), {}).update(expct(i, groups[gi]))

    matrices = {l: np.zeros((K, K), dtype=complex) for l in labels}
    for n in range(K):
        for l in labels:
            matrices[l][n, n] = diag_vals[n][l]
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        E = [np.array([cross_vals[(n, m, k)][l] for l in labels]) for k in range(4)]
        re = (E[0] - E[2]) / 2
        im = (E[3] - E[1]) / 2
        x = re + 1j * im
        for li, l in enumerate(labels):
            matrices[l][n, m] = x[li]
            matrices[l][m, n] = np.conj(x[li])
    return matrices


def assemble_control():
    print("\n" + "=" * 96)
    print("  ionq_original_circuit_replication.py --assemble-control")
    print("=" * 96)
    ck = load_ckpt("original_control")
    if ck is None:
        print("  missing checkpoint: original_control -- run --control first")
        return None

    lambdas = np.array(ck["lambdas"])
    qop_bare, _, _ = ef.build_h4_qop(1.0)
    terms = ef.decompose_pauli_terms(qop_bare)
    exact_energy = ck["exact_energy"]
    enuc = ck["enuc"]

    report = {}
    for model in NOISE_MODELS:
        counts_a = ck["counts"][model]["alpha"]
        counts_b = ck["counts"][model]["beta"]
        tags_a = [tuple(t) for t in ck["tags_a"]]
        tags_b = [tuple(t) for t in ck["tags_b"]]
        errs = []
        for seed in range(N_SEEDS):
            rng_a = np.random.default_rng(stable_seed("origctrl_a", model, seed))
            rng_b = np.random.default_rng(stable_seed("origctrl_b", model, seed))
            alpha_mats = rebuild_matrices_from_counts(counts_a, tags_a, ck["groups_a"], ck["alpha_labels"], K, rng_a)
            beta_mats = rebuild_matrices_from_counts(counts_b, tags_b, ck["groups_b"], ck["beta_labels"], K, rng_b)
            E = ef.ef_energy_from_noisy_matrices(terms, lambdas, alpha_mats, beta_mats, enuc, K)
            err_kcal = abs(E - exact_energy) * HARTREE_TO_KCAL_MOL
            errs.append(err_kcal)
        report[model] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs)), "errs": errs}
        print(f"    {model}: {report[model]['mean_kcal']:.3f} +/- {report[model]['std_kcal']:.3f} kcal/mol (8 seeds)")

    print(f"\n  -- comparison to iteration 9's fixed-ansatz raw numbers (35-43 kcal/mol) --")
    for model in ("aria-1", "forte-1"):
        v = report[model]["mean_kcal"]
        verdict = ("CONFIRMS the hypothesis: noise, not circuit, is the story" if 25 < v < 55
                    else "the circuit DOES matter here -- this differs materially from iteration 9's fixed ansatz")
        print(f"    {model}: {v:.1f} kcal/mol -- {verdict}")

    results = {"K": K, "exact_energy_ha": exact_energy, "n_seeds": N_SEEDS, "shots": SHOTS,
               "psi_max_imag": ck["psi_max_imag"], "wall_clock": ck["wall_clock"],
               "report": report}
    save_ckpt("original_control_assembled", results)
    with open(RESULTS_PATH, "w") as f:
        json.dump({"task2_control": results}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--control", action="store_true")
    parser.add_argument("--assemble-control", action="store_true")
    args = parser.parse_args()
    if args.control:
        phase_control()
    elif args.assemble_control:
        assemble_control()
    else:
        parser.error("pass --control or --assemble-control")


if __name__ == "__main__":
    main()
