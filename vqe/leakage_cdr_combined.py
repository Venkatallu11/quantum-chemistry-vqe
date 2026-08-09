#!/usr/bin/env python3
"""
leakage_cdr_combined.py — does removing particle-number leakage (Task D,
iteration 18) BEFORE fitting CDR's scale correction rescue CDR, which
made things 2.1-2.6x WORSE on real IonQ hardware (iteration 9) when
trained/applied on raw (leakage-contaminated) data?
============================================================================
MOTIVATION: CDR fits a per-label/per-slot linear scale correction
(noisy ~= f * exact, weighted least squares through the origin) from
training circuits whose exact value is classically known. This is only
a good MODEL of the noise if the noisy-vs-exact relationship really is
approximately linear/scale-like. A discrete, outlier-injecting error
process (leakage -- a shot landing entirely outside the valid sector,
iteration 18: ~7-8% of real shots) is exactly the kind of contamination
that would corrupt a LINEAR fit's slope, especially with only
N_TRAIN_PER_SLOT=5 training draws per slot (established CDR default,
`cdr_mitigation.py`) -- a handful of leakage-corrupted outliers among 5
points can swing a least-squares slope substantially. If CDR's real-
hardware failure (iteration 9) was partly caused by this, fitting CDR on
LEAKAGE-POSTSELECTED training AND target data (using iteration 18's
verified-exact ancilla scheme, ALREADY confirmed to roughly halve the
exact error at every noise scale in leakage_zne_floor_tested.py) should
recover at least some of CDR's lost value.

Tests, all on the SAME local synthetic noise model (P2_PER_GATE-scaled
depolarizing, scale=1, `leakage_zne_floor_tested.py`'s reused
noisy_probs_5q machinery) so all four numbers are directly comparable:
  (a) RAW           -- no CDR, no leakage removal
  (b) CDR only       -- iteration 9's original recipe, no leakage removal
  (c) leakage only    -- iteration 18's recipe, no CDR
  (d) CDR + leakage  -- the new combination being tested here

Run:
    python vqe/leakage_cdr_combined.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL
from leakage_zne_floor_tested import (
    noisy_probs_5q, build_noise_model, sample_counts, counts_to_probs_dict,
    expectation_raw_from_probs_dict, expectation_post_from_probs_dict,
)
from spin_leakage_postselect_ionq import with_ancilla_parity
from fixed_ansatz import build_ansatz
import ef_fragment as effrag
from qiskit.quantum_info import Statevector, Pauli

K = 6
SHOTS = 100_000
N_SEEDS = 8
N_TRAIN_PER_SLOT = 5
LOW_SIGNAL_CUTOFF = 0.05
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "leakage_cdr_combined_results.json")


def exact_labels_from_angles(angles, labels):
    sv = Statevector.from_instruction(build_ansatz(angles))
    return {l: float(sv.expectation_value(Pauli(l)).real) for l in labels}


def noisy_labels_from_angles(angles, groups, noise_model, rng, postselect, use_shots):
    """Full pipeline: ancilla circuit -> per-group density matrix (or
    exact probs if use_shots=False) -> raw or post-selected expectation
    per label. Reuses the SAME machinery as leakage_zne_floor_tested.py
    for consistency."""
    out = {}
    for group in groups:
        combined = effrag.combined_basis_label(group)
        probs32 = noisy_probs_5q(angles, combined, noise_model)
        if use_shots:
            counts = sample_counts(probs32, SHOTS, rng)
            probs_dict = counts_to_probs_dict(counts)
        else:
            probs_dict = {format(i, "05b"): probs32[i] for i in range(32) if probs32[i] > 0}
        for l in group:
            out[l] = (expectation_post_from_probs_dict(probs_dict, l) if postselect
                       else expectation_raw_from_probs_dict(probs_dict, l))
    return out


def fit_scale(pairs):
    if len(pairs) < 3:
        return None
    exact = np.array([e for e, _ in pairs])
    noisy = np.array([n for _, n in pairs])
    denom = float(np.sum(exact * exact))
    if denom < 1e-12:
        return None
    return float(np.sum(exact * noisy) / denom)


def fit_all_scales(training, non_id_labels, target_names):
    def pairs_for(label=None, slot=None):
        out = []
        for row in training:
            if label is not None and row["label"] != label:
                continue
            if slot is not None and row["slot"] != slot:
                continue
            if abs(row["exact"]) < LOW_SIGNAL_CUTOFF:
                continue
            out.append((row["exact"], row["noisy"]))
        return out

    global_scale = fit_scale(pairs_for())
    per_basis = {l: (fit_scale(pairs_for(label=l)) or global_scale) for l in non_id_labels}
    per_slot = {n: (fit_scale(pairs_for(slot=n)) or global_scale) for n in target_names}
    return global_scale, per_basis, per_slot


def energy_and_err(p, raw, K):
    from qforge import combine_matrices, energy_from_alpha_matrices
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  leakage_cdr_combined.py -- does leakage removal rescue CDR?")
    print("=" * 96)

    t0 = time.time()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    assert n_ok == 36
    print(f"  setup OK: 36/36 converged, {time.time()-t0:.1f}s")

    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    groups = effrag.group_labels_qubit_wise(alpha_labels)
    target_names = sorted(solutions.keys())
    nm1 = build_noise_model(1)

    print(f"\n  generating {N_TRAIN_PER_SLOT} training draws x {len(target_names)} slots "
          f"(exact + RAW-noisy + POST-SELECTED-noisy per draw)")
    t0 = time.time()
    rng_train = np.random.default_rng(2024)
    training_raw, training_post = [], []
    for name in target_names:
        for _ in range(N_TRAIN_PER_SLOT):
            angles = rng_train.uniform(-np.pi, np.pi, 5)
            exact_vals = exact_labels_from_angles(angles, non_id_labels)
            raw_vals = noisy_labels_from_angles(angles, groups, nm1, rng_train, postselect=False, use_shots=False)
            post_vals = noisy_labels_from_angles(angles, groups, nm1, rng_train, postselect=True, use_shots=False)
            for l in non_id_labels:
                training_raw.append({"slot": name, "label": l, "exact": exact_vals[l], "noisy": raw_vals[l]})
                training_post.append({"slot": name, "label": l, "exact": exact_vals[l], "noisy": post_vals[l]})
    print(f"  training data generated, {time.time()-t0:.1f}s")

    gs_raw, pb_raw, ps_raw = fit_all_scales(training_raw, non_id_labels, target_names)
    gs_post, pb_post, ps_post = fit_all_scales(training_post, non_id_labels, target_names)
    print(f"  fitted global scale: RAW={gs_raw:.4f}  POST-SELECTED={gs_post:.4f}")

    def corrected(name, label, val, per_basis, per_slot, global_scale):
        f = per_basis.get(label) or global_scale
        f2 = per_slot.get(name) or global_scale
        scale = f if f is not None else 1.0
        return val / scale if abs(scale) > 1e-9 else val

    print(f"\n  -- 8-seed shot-noisy evaluation on the actual 36 targets, scale=1 --")
    schemes = ["raw", "cdr_only", "leakage_only", "cdr_plus_leakage"]
    errs_by_scheme = {s: [] for s in schemes}

    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed * 104729 + 7)
        raw_meas = {n: noisy_labels_from_angles(solutions[n]["angles"], groups, nm1, rng,
                                                  postselect=False, use_shots=True) for n in target_names}
        rng2 = np.random.default_rng(seed * 104729 + 7)  # SAME draw stream -> same shots, only post-selection differs
        post_meas = {n: noisy_labels_from_angles(solutions[n]["angles"], groups, nm1, rng2,
                                                   postselect=True, use_shots=True) for n in target_names}

        raw_dict = {n: dict(raw_meas[n]) for n in target_names}
        cdr_only_dict = {n: {l: corrected(n, l, raw_meas[n][l], pb_raw, ps_raw, gs_raw) for l in non_id_labels}
                          for n in target_names}
        leakage_only_dict = {n: dict(post_meas[n]) for n in target_names}
        cdr_leakage_dict = {n: {l: corrected(n, l, post_meas[n][l], pb_post, ps_post, gs_post) for l in non_id_labels}
                             for n in target_names}

        for scheme, d in zip(schemes, [raw_dict, cdr_only_dict, leakage_only_dict, cdr_leakage_dict]):
            _, err = energy_and_err(p, d, K)
            errs_by_scheme[scheme].append(err)

    print(f"\n  {'scheme':<20} {'mean_kcal':>10} {'std_kcal':>9}")
    summary = {}
    for scheme in schemes:
        errs = errs_by_scheme[scheme]
        mean_e, std_e = float(np.mean(errs)), float(np.std(errs))
        summary[scheme] = {"mean_kcal": mean_e, "std_kcal": std_e}
        print(f"  {scheme:<20} {mean_e:>10.3f} {std_e:>9.3f}")

    results = {"K": K, "shots": SHOTS, "n_seeds": N_SEEDS, "n_train_per_slot": N_TRAIN_PER_SLOT,
               "summary": summary, "chemical_accuracy_kcal": 1.0}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
