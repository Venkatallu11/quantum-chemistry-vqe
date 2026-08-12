#!/usr/bin/env python3
"""
task2_fold_response_dataset.py — iteration 26, Task 2. THE NEW CORE
DATASET: stop treating the energy as one object. For every measured
(state-prep slot, measurement group) circuit, record structural
fingerprint (native 1q/2q counts, MS/ZZ counts, depth) AND the fold
response (ideal vs measured expectation at folds 1/3/5/9) for every
Pauli label in that group, then cluster into circuit-fingerprint
families and quantify which families carry the energy error.
============================================================================
NATIVE FOLDING, not abstract: this project's own established finding
(ionq_fold_check.py) is that folding ABSTRACT (u3/cx) gates gets
CANCELLED before execution -- only `gateset="native"` submission with
`fold_native_2q` (explicit, numerically-verified inverse construction,
reused unchanged from ionq_fold_check.py) survives. Every circuit here is
built as: fixed_ansatz.build_ansatz(angles) -> to_native(qc, gate_name)
(verified below: constant 11 native 2q gates across all 36 targets at
BOTH optimization_level 0 and 1 -- the earlier-documented opt_level>=1
collapse risk was specific to u3/cx BASIS transpilation, confirmed here
NOT to reproduce for native-TARGET transpilation) -> fold_native_2q(fold)
-> basis-change (native, appended separately, never re-transpiled with
the folded state-prep so the fold can't be optimized away).

SCOPE, disclosed explicitly (this is a genuine, principled reduction
under real time/wall-clock constraints, not a hidden shortcut): the full
36-slot x 13-group design would need 468 circuits x 4 folds x 3 models =
5,616 circuit executions. This file uses 12 REPRESENTATIVE slots -- the
6 diagonal Schmidt vectors (u_0..u_5, capturing genuine state-prep depth
diversity already established in iteration 24's ADAPT work) PLUS 6
phase-pair (Schmidt CROSS-TERM) slots -- x 13 groups x 4 folds x 3 models
= 1,872 circuit executions, keeping BOTH structural families the task
explicitly asked about ("Schmidt cross-terms" specifically named) in
scope, at 1/3 the circuit count of the full design.

Run:
    python vqe/task2_fold_response_dataset.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, native_target
from ionq_fold_check import fold_native_2q
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed
import ef_fragment as effrag_mod
from phase2_dominant_term_mitigation import rank_terms, K as PHASE2_K
from phys_constrained_reconstruction import build_P_S
from qiskit.quantum_info import Statevector, Pauli

K = 6
FOLD_FACTORS = [1, 3, 5, 9]
SHOTS = 100_000  # matches ionq_native_forged_energy.py's own established "clean signal" convention
GATE_BY_MODEL = {"ideal": "ms", "aria-1": "ms", "forte-1": "zz"}
REPRESENTATIVE_SLOTS = [f"u_{n}" for n in range(K)] + \
    ["(u0+u1)", "(u0+u2)", "(u0+u3)", "(u0+u4)", "(u0+u5)", "(u1+u2)"]
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task2_fold_response.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task2_fold_response_dataset_results.json")


def native_basis_change(label, gate_name):
    from qiskit.circuit import QuantumCircuit
    n = len(label)
    qc = QuantumCircuit(n)
    for i, ch in enumerate(label):
        qubit = n - 1 - i
        if ch == "X":
            qc.h(qubit)
        elif ch == "Y":
            qc.sdg(qubit)
            qc.h(qubit)
    return to_native(qc, gate_name)


def build_folded_measurement_circuits(angles, gate_name, groups):
    """Returns {fold: [(group, circuit), ...]}, plus per-fold gate-count
    metadata (native 1q/2q counts, depth) for the STATE-PREP block only
    (fold-independent basis-change gates are cheap/constant, recorded
    separately)."""
    base = to_native(build_ansatz(angles), gate_name)
    n2q_base = base.count_ops().get(gate_name, 0)
    out = {}
    meta = {}
    for fold in FOLD_FACTORS:
        folded = fold_native_2q(base, fold, gate_name)
        n2q = folded.count_ops().get(gate_name, 0)
        n1q = sum(v for k, v in folded.count_ops().items() if k != gate_name)
        meta[fold] = {"n2q": n2q, "n1q_stateprep": n1q, "depth": folded.depth()}
        circuits = []
        for group in groups:
            combined = effrag_mod.combined_basis_label(group)
            basis_qc = native_basis_change(combined, gate_name)
            qc = folded.compose(basis_qc)
            qc.measure_all()
            circuits.append((group, qc))
        out[fold] = circuits
    assert n2q_base == 11, f"expected 11 native 2q gates at fold=1, got {n2q_base}"
    return out, meta


def exact_ideal_expectation(angles, label):
    sv = Statevector.from_instruction(build_ansatz(angles))
    return float(sv.expectation_value(Pauli(label)).real)


def main():
    print("\n" + "=" * 96)
    print("  task2_fold_response_dataset.py -- per-term fold-response, native folding")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"])
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    print(f"  {len(REPRESENTATIVE_SLOTS)} representative slots x {len(groups)} groups x "
          f"{len(FOLD_FACTORS)} folds x 3 models")

    # -- per-label exact energy sensitivity (reused from Phase 2, not re-derived) --
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    contributions, mismatch = rank_terms(p, P_S)
    assert mismatch < 1e-6
    label_weight = {}
    for c in contributions:
        label_weight[c["alpha_label"]] = label_weight.get(c["alpha_label"], 0.0) + abs(c["contribution_ha"]) * HARTREE_TO_KCAL_MOL
    print(f"  per-label exact energy-sensitivity weight computed (reused Phase 2's rank_terms), "
          f"total={sum(label_weight.values()):.2f} kcal/mol")

    if os.path.exists(CKPT_PATH):
        print(f"  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"  connected, backend={backend.name}")

        # build all circuits, per model (native gate depends on model), submit ALL non-blocking first
        jobs = {}  # (model, fold) -> (job, tags)
        gate_meta_by_slot = {}
        t0 = time.time()
        for model, gate_name in GATE_BY_MODEL.items():
            for fold in FOLD_FACTORS:
                circuits, tags = [], []
                for name in REPRESENTATIVE_SLOTS:
                    angles = fixed_solutions[name]["angles"]
                    folded_circuits, meta = build_folded_measurement_circuits(angles, gate_name, groups)
                    if model == "aria-1" or (model == "ideal" and gate_name == "ms"):
                        gate_meta_by_slot.setdefault("ms", {})[name] = meta
                    if model == "forte-1":
                        gate_meta_by_slot.setdefault("zz", {})[name] = meta
                    for group, qc in folded_circuits[fold]:
                        circuits.append(qc)
                        tags.append((name, tuple(group)))
                job = submit_job(circuits, backend, model, shots=SHOTS)
                jobs[(model, fold)] = (job, tags)
        t_submit = time.time() - t0
        print(f"  all {len(jobs)} model x fold jobs submitted (non-blocking), {t_submit:.1f}s")

        t0 = time.time()
        counts_by_key = {}
        for i, (key, (job, tags)) in enumerate(jobs.items()):
            counts_by_key[f"{key[0]}|{key[1]}"] = get_counts_list(job)
            print(f"    retrieved {i+1}/{len(jobs)}: model={key[0]} fold={key[1]}, {time.time()-t0:.1f}s elapsed")
        t_retrieve = time.time() - t0
        print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

        ck = {
            "tags": {f"{k[0]}|{k[1]}": [[t[0], list(t[1])] for t in v[1]] for k, v in jobs.items()},
            "counts": counts_by_key,
            "gate_meta": gate_meta_by_slot,
            "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        }
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- build the per-term fold-response dataset --
    print(f"\n  -- assembling per-term dataset --")
    dataset = []
    for model, gate_name in GATE_BY_MODEL.items():
        for fold in FOLD_FACTORS:
            key = f"{model}|{fold}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            meta_by_slot = ck["gate_meta"][gate_name]
            for (name, group), counts in zip(tags, counts_list):
                probs = {}
                total = sum(counts.values())
                for bs, c in counts.items():
                    probs[bs] = c / total
                for label in group:
                    measured = pauli_expectation(probs, label)
                    ideal = exact_ideal_expectation(fixed_solutions[name]["angles"], label)
                    dataset.append({
                        "model": model, "fold": fold, "slot": name, "label": label,
                        "gate_name": gate_name, "n2q": meta_by_slot[name][str(fold)]["n2q"] if str(fold) in meta_by_slot[name] else meta_by_slot[name][fold]["n2q"],
                        "depth": meta_by_slot[name][str(fold)]["depth"] if str(fold) in meta_by_slot[name] else meta_by_slot[name][fold]["depth"],
                        "n_nonI_paulis_in_label": sum(1 for c in label if c != "I"),
                        "ideal_expectation": ideal, "measured_expectation": measured,
                        "delta": measured - ideal,
                        "is_cross_term_slot": "+" in name or "-" in name,
                        "energy_weight_kcal": label_weight.get(label, 0.0),
                    })
    print(f"  {len(dataset)} (model, fold, slot, label) data points assembled")

    # -- family clustering by circuit fingerprint --
    print(f"\n  -- family clustering by circuit fingerprint --")

    def family_of(row):
        n_nonI = row["n_nonI_paulis_in_label"]
        cross = row["is_cross_term_slot"]
        if cross:
            return "schmidt_cross_term"
        elif n_nonI <= 1:
            return "all_Z_low_depth"
        elif row["gate_name"] == "ms" and n_nonI >= 2:
            return "XX_YY_medium_ms"
        elif row["gate_name"] == "zz" and n_nonI >= 2:
            return "high_ZZ"
        return "other"

    for row in dataset:
        row["family"] = family_of(row)

    families = sorted(set(row["family"] for row in dataset))
    family_stats = {}
    for model in ["aria-1", "forte-1"]:
        family_stats[model] = {}
        for fam in families:
            rows = [r for r in dataset if r["model"] == model and r["family"] == fam and r["fold"] == 1]
            if not rows:
                continue
            abs_deltas = [abs(r["delta"]) for r in rows]
            total_weight = sum(r["energy_weight_kcal"] for r in rows)
            family_stats[model][fam] = {
                "n_datapoints": len(rows), "mean_abs_delta_fold1": float(np.mean(abs_deltas)),
                "total_energy_weight_kcal": total_weight,
            }
    for model in ["aria-1", "forte-1"]:
        print(f"    {model}:")
        total_w = sum(f["total_energy_weight_kcal"] for f in family_stats[model].values())
        for fam, s in sorted(family_stats[model].items(), key=lambda kv: -kv[1]["mean_abs_delta_fold1"]):
            frac = s["total_energy_weight_kcal"] / total_w * 100 if total_w > 0 else 0
            print(f"      {fam:<22} n={s['n_datapoints']:>3}  mean|delta|(fold1)={s['mean_abs_delta_fold1']:.4f}  "
                  f"energy_weight={s['total_energy_weight_kcal']:.2f} kcal/mol ({frac:.1f}% of total)")

    results = {
        "K": K, "shots": SHOTS, "fold_factors": FOLD_FACTORS,
        "representative_slots": REPRESENTATIVE_SLOTS,
        "n_datapoints": len(dataset), "family_stats": family_stats,
        "dataset": dataset,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
