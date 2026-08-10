#!/usr/bin/env python3
"""
task5_reproducibility_gate.py — iteration 24, Task 5. NO CLAIM WITHOUT IT.
Every headline number must satisfy |Delta E| < 1 kcal/mol across
INDEPENDENT repetitions -- separate submissions, separate seeds, separate
calibration draws -- not one lucky run.
============================================================================
HONEST DISTINCTION, stated up front: this project's standard "8-seed
mean +/- std" already used in every prior iteration is BOOTSTRAP
resampling of ONE underlying real-hardware submission's counts -- it
estimates SHOT-NOISE spread, not independence from calibration drift,
queue-time effects, or anything else that could differ between two
actually-separate API calls. Task 5 explicitly distinguishes these
("separate submissions... not one lucky run"), so this file checks BOTH,
separately, and does not conflate them:

  A) GENUINE CROSS-SUBMISSION reproducibility: this project's history
     already contains ONE case of the exact same circuit set submitted
     TWICE, independently, to IonQ's free simulator -- the Z2-tapered
     raw circuit (z2_tapered_targets.json and z2_tapered_targets_run1.json,
     confirmed independent by their differing wall-clock timestamps, not
     the same job re-read). This is REAL "separate submission" evidence,
     not simulated -- the only such case found in this project's checkpoint
     history.

  B) SPLIT-HALF seed reproducibility for this session's two most novel
     new numbers (Task 3's 21-circuit subspace tomography, Task 4's PSD+
     leakage combination): seeds 0-3 vs seeds 4-7, treated as two
     independent halves of the same bootstrap draw. This is weaker
     evidence than (A) (same underlying counts, not a new submission) but
     still catches "one lucky seed" cherry-picking, which is explicitly
     what iteration 2's 0.0636 and the classic 0.57 kcal/mol results were.

  C) AN HONEST GAP, surfaced not hidden: Phase 1/2's headline numbers and
     iteration 18's leakage-postselection best result (31.77/33.86) have
     NEVER been checked against a genuinely independent second real
     submission anywhere in this project's history -- only single-
     submission bootstrap pseudoreplication exists for them. Task 5's own
     gate, applied honestly, means these numbers do NOT yet carry (A)-type
     evidence -- stated plainly here, not glossed over.

Run:
    python vqe/task5_reproducibility_gate.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL, slot_names
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed, expectation_from_counts
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot, build_group_index, K
from z2_tapered_ionq import build_reduced_problem
from phase2_dominant_term_mitigation import rank_terms, head_tail_split
from spin_leakage_postselect_ionq import postselect_counts

CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
STD_CKPT_PATH = os.path.join(CKPT_DIR, "targets_d1.0.json")
LEAK_CKPT_PATH = os.path.join(CKPT_DIR, "spin_leakage_targets.json")
SHOTS = 10_000
GATE_KCAL = 1.0
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task5_reproducibility_gate_results.json")


def z2_tapered_energy_from_checkpoint(ckpt_name, problem, seeds):
    with open(os.path.join(CKPT_DIR, f"{ckpt_name}.json")) as f:
        ck = json.load(f)
    p = problem["p"]
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    reduced_label_map = {k: tuple(v) for k, v in ck["reduced_label_map"].items()}
    target_names = ck["target_names"]
    groups = ck["groups"]
    tags = [tuple(t) for t in ck["tags"]]
    idx_map = [name for name in target_names for _ in groups]

    report = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])
        errs = []
        for seed in seeds:
            rng = np.random.default_rng(stable_seed("z2_tapered", model, seed))
            raw = {name: {} for name in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    counts = bootstrap_counts(per_name[name][gi], 100_000, rng)
                    group_vals = {rl: expectation_from_counts(counts, rl) for rl in group}
                    for orig_label in alpha_labels:
                        rl, sign = reduced_label_map[orig_label]
                        if rl in group_vals:
                            raw[name][orig_label] = sign * group_vals[rl]
            alpha_mats = combine_matrices(raw, alpha_labels, identity_label, K)
            E, err = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                 exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs.append(err["err_vs_exact_kcal"])
        report[model] = {"mean": float(np.mean(errs)), "std": float(np.std(errs))}
    return report


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def task3_phys21_errs(p, ck, P_S, seeds):
    """Recomputes Task 3's phys21 (21-circuit, joint SDP) error for a
    given list of seeds -- reused directly from task3_subspace_tomography's
    own method, not reimplemented differently."""
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    groups = ck["groups"]
    group_idx = build_group_index(groups)
    names = ck["target_names"]
    diag_names = [f"u_{n}" for n in range(K)]
    plus_pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    plus_names = [f"(u{n}+u{m})" for n, m in plus_pairs]
    kept_names = diag_names + plus_names

    out = {}
    for model in ["aria-1", "forte-1"]:
        counts_by_slot = ck["counts"][model]
        errs = []
        for seed in seeds:
            rng = np.random.default_rng(stable_seed("subspace_tomo", model, seed))
            resampled = {}
            for name in names:
                resampled[name] = [bootstrap_counts(counts_by_slot[name][gi], SHOTS, rng) for gi in range(len(groups))]

            def m_and_w_for(name):
                m_dict, w_dict = {}, {}
                for l in non_id_labels:
                    counts = resampled[name][group_idx[l]]
                    m = expectation_from_counts(counts, l)
                    m_dict[l] = m
                    total = sum(counts.values())
                    w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(total, 1)
                return m_dict, w_dict

            raw_m = {name: m_and_w_for(name)[0] for name in names}
            rho21 = {}
            for name in kept_names:
                m_dict, w_dict = m_and_w_for(name)
                rho21[name] = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
            phys21 = {name: {l: float(np.real(np.trace(rho21[name] @ P_S[l]))) for l in non_id_labels} for name in diag_names}
            for (n, m), plus in zip(plus_pairs, plus_names):
                un, um = f"u_{n}", f"u_{m}"
                plus_vals = {l: float(np.real(np.trace(rho21[plus] @ P_S[l]))) for l in non_id_labels}
                phys21[plus] = dict(plus_vals)
                synth_minus = {}
                for l in non_id_labels:
                    cross = plus_vals[l] - (phys21[un][l] + phys21[um][l]) / 2
                    synth_minus[l] = plus_vals[l] - 2 * cross
                phys21[f"(u{n}-u{m})"] = synth_minus
            _, err = energy_and_err(p, phys21, K)
            errs.append(err["err_vs_exact_kcal"])
        out[model] = {"mean": float(np.mean(errs)), "std": float(np.std(errs)), "n": len(seeds)}
    return out


def task4_psd_leakage_errs(p, leak_ck, P_S, seeds):
    """Recomputes Task 4's PSD+leakage error for a given list of seeds."""
    alpha_labels = leak_ck["alpha_labels"]
    identity_label = leak_ck["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    leak_groups = leak_ck["groups"]
    leak_group_idx = build_group_index(leak_groups)
    names = leak_ck["target_names"]

    idx_map = leak_ck["idx_map"]
    leak_blocks = {}
    i = 0
    while i < len(idx_map):
        name = idx_map[i]
        leak_blocks[name] = list(range(i, i + 13))
        i += 13

    out = {}
    for model in ["aria-1", "forte-1"]:
        counts_flat = leak_ck["counts"][model]
        errs = []
        for seed in seeds:
            rng = np.random.default_rng(stable_seed("task4_leak", model, seed))
            phys_leak = {name: {} for name in names}
            for name in names:
                block = leak_blocks[name]
                resampled5 = [bootstrap_counts(counts_flat[block[gi]], SHOTS, rng) for gi in range(len(leak_groups))]
                m_dict, w_dict = {}, {}
                for l in non_id_labels:
                    c5 = resampled5[leak_group_idx[l]]
                    c4 = postselect_counts(c5)
                    total = sum(c4.values())
                    m = expectation_from_counts(c4, l) if total > 0 else 0.0
                    m_dict[l] = m
                    w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(total, 1)
                rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                for l in non_id_labels:
                    phys_leak[name][l] = float(np.real(np.trace(rho_slot @ P_S[l])))
            _, err = energy_and_err(p, phys_leak, K)
            errs.append(err["err_vs_exact_kcal"])
        out[model] = {"mean": float(np.mean(errs)), "std": float(np.std(errs)), "n": len(seeds)}
    return out


def main():
    print("\n" + "=" * 96)
    print("  task5_reproducibility_gate.py -- no claim without it")
    print("=" * 96)
    print(f"  GATE: |Delta E| < {GATE_KCAL} kcal/mol between independent repetitions")

    gate_results = {}

    # -- A) genuine cross-submission reproducibility: Z2-tapered raw, two real separate submissions --
    print(f"\n  -- A) GENUINE cross-submission reproducibility (Z2-tapered raw, two real IonQ submissions) --")
    problem = build_reduced_problem()
    seeds_full = list(range(8))
    rep_run0 = z2_tapered_energy_from_checkpoint("z2_tapered_targets", problem, seeds_full)
    rep_run1 = z2_tapered_energy_from_checkpoint("z2_tapered_targets_run1", problem, seeds_full)
    for model in ["aria-1", "forte-1"]:
        d = abs(rep_run0[model]["mean"] - rep_run1[model]["mean"])
        passed = d < GATE_KCAL
        gate_results[f"z2_tapered_cross_submission_{model}"] = {
            "run0_mean": rep_run0[model]["mean"], "run1_mean": rep_run1[model]["mean"],
            "delta_kcal": d, "passed": bool(passed),
        }
        print(f"    {model}: run0={rep_run0[model]['mean']:.2f}  run1={rep_run1[model]['mean']:.2f}  "
              f"Delta={d:.3f} kcal/mol  {'PASS' if passed else '*** FAIL ***'}")

    # -- B) split-half seed reproducibility for this session's new numbers --
    print(f"\n  -- B) split-half seed reproducibility (seeds 0-3 vs 4-7, same underlying submission) --")
    with open(STD_CKPT_PATH) as f:
        std_ck = json.load(f)
    with open(LEAK_CKPT_PATH) as f:
        leak_ck = json.load(f)
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=std_ck["d"], K=K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(std_ck["alpha_labels"], U)

    half_a, half_b = [0, 1, 2, 3], [4, 5, 6, 7]

    print(f"    Task 3 (21-circuit subspace tomography, phys21):")
    t3_a = task3_phys21_errs(p, std_ck, P_S, half_a)
    t3_b = task3_phys21_errs(p, std_ck, P_S, half_b)
    for model in ["aria-1", "forte-1"]:
        d = abs(t3_a[model]["mean"] - t3_b[model]["mean"])
        passed = d < GATE_KCAL
        gate_results[f"task3_phys21_splithalf_{model}"] = {
            "half_a_mean": t3_a[model]["mean"], "half_b_mean": t3_b[model]["mean"],
            "delta_kcal": d, "passed": bool(passed),
        }
        print(f"      {model}: seeds0-3={t3_a[model]['mean']:.2f}  seeds4-7={t3_b[model]['mean']:.2f}  "
              f"Delta={d:.3f} kcal/mol  {'PASS' if passed else '*** FAIL ***'}")

    print(f"    Task 4 (PSD reconstruction + leakage postselection):")
    t4_a = task4_psd_leakage_errs(p, leak_ck, P_S, half_a)
    t4_b = task4_psd_leakage_errs(p, leak_ck, P_S, half_b)
    for model in ["aria-1", "forte-1"]:
        d = abs(t4_a[model]["mean"] - t4_b[model]["mean"])
        passed = d < GATE_KCAL
        gate_results[f"task4_psd_leakage_splithalf_{model}"] = {
            "half_a_mean": t4_a[model]["mean"], "half_b_mean": t4_b[model]["mean"],
            "delta_kcal": d, "passed": bool(passed),
        }
        print(f"      {model}: seeds0-3={t4_a[model]['mean']:.2f}  seeds4-7={t4_b[model]['mean']:.2f}  "
              f"Delta={d:.3f} kcal/mol  {'PASS' if passed else '*** FAIL ***'}")

    # -- C) the honest gap --
    print(f"\n  -- C) HONEST GAP: numbers with NO independent-submission evidence in this project's history --")
    no_cross_submission_evidence = [
        "Phase 1 physics-constrained reconstruction (27.71/31.64 kcal/mol)",
        "Phase 2 selective hybrid (28.24/32.52 kcal/mol)",
        "iteration 18 leakage postselection, this project's best real result (31.77/33.86 kcal/mol)",
    ]
    for item in no_cross_submission_evidence:
        print(f"    - {item}: only single-submission 8-seed bootstrap exists, no second real submission run")

    n_pass = sum(1 for v in gate_results.values() if v["passed"])
    n_total = len(gate_results)
    print(f"\n  -- OVERALL: {n_pass}/{n_total} reproducibility checks PASS the <{GATE_KCAL} kcal/mol gate --")

    out = {
        "gate_kcal": GATE_KCAL, "gate_results": gate_results,
        "n_pass": n_pass, "n_total": n_total,
        "no_cross_submission_evidence": no_cross_submission_evidence,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return out


if __name__ == "__main__":
    main()
