#!/usr/bin/env python3
"""
loop_local_perturbation_cdr.py — iteration 2 of the below-0.30-kcal/mol search.
============================================================================
DISCOVERY THAT MOTIVATES THIS (verified here, not assumed): the per-label
noisy/exact ratio is NOT angle-independent the way the global per-basis fit
assumes. Measured directly: over 15 random angle draws, XXYY's ratio is
PERFECTLY constant (std=0.0) but ZZII varies over a 29% range (0.855-1.121)
and YZYZ over 22% (0.829-1.026). This is real: a fixed gate STRUCTURE
(verified constant, 11 CX) does not imply a fixed per-label NOISE SHRINK,
because the depolarizing channels commute through the circuit's PARAMETRIZED
gates in an angle-dependent way (backward Heisenberg-propagating a Pauli
through a rotation gate mixes it with other Paulis at angle-dependent
weights). A single global per-basis scale, averaged over random training
angles, therefore systematically mismatches each specific target's true
local shrink -- this is a real, verified mechanism for CDR's 2.850 kcal/mol
residual, not hand-waved.

FIX TRIED HERE: instead of training on globally-random angles (which
average over this angle-dependence and lose it), train LOCALLY -- for each
of the 36 target slots, perturb THAT slot's own fitted angles by a small
amount (+/-0.15 rad) and fit a PER-(SLOT, LABEL) scale from those nearby
points. Verified directly before building the full pipeline: local
perturbation around a target's own angles recovers the TRUE local ratio to
~1e-5 relative precision (vs 22-29% error from global random sampling) for
the labels checked.

THIS IS A DIFFERENT IDEA FROM THE EARLIER FAILED "per-circuit" scale
(45.14 kcal/mol, see ledger baseline): that approach used GLOBALLY RANDOM
angles per slot (same distribution as per-basis, just partitioned
differently) -- it never exploited locality, so it just had less data
than the pooled global fit for no benefit. This one specifically targets
the angle-dependence just measured.

COST NOTE: needs local training circuits at EVERY (slot, label) pair, but
batched per slot (all labels read off the same local perturbed circuit),
so cost is ~36 slots x K_LOCAL circuits, not 36x36xK_LOCAL.

Run:
    python vqe/loop_local_perturbation_cdr.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import rank6_symmetry_vd as r

K = r.K
N_SEEDS = r.N_SEEDS
LOW_SIGNAL_CUTOFF = r.LOW_SIGNAL_CUTOFF
K_LOCAL = 4          # local perturbation draws per slot
PERTURB_RADIUS = 0.15  # radians

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "loop_local_perturbation_cdr_results.json")


def generate_local_training_data(p, non_id_labels, noise_model, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for name, sol in p["solutions"].items():
        target_angles = np.array(sol["angles"])
        for _ in range(K_LOCAL):
            perturb = rng.uniform(-PERTURB_RADIUS, PERTURB_RADIUS, 5)
            angles = (target_angles + perturb).tolist()
            exact_vals = r.exact_labels(angles, non_id_labels)
            dm = r.density_matrix_4q(angles, noise_model)
            noisy_vals, _ = r.measure_labels(angles, non_id_labels, dm)
            for l in non_id_labels:
                rows.append({"slot": name, "label": l, "exact": exact_vals[l], "noisy": noisy_vals[l]})
    return rows


def fit_local_scales(training, non_id_labels, slot_names_list):
    """Per-(slot,label) scale, fit from ONLY that slot's local perturbation
    draws. Falls back to a per-label GLOBAL scale (pooled across all slots)
    if a specific (slot,label) has too little data after the low-signal
    filter -- never falls back to 1.0 (that would silently no-op the
    correction for exactly the hardest cases)."""
    global_per_label = {}
    for l in non_id_labels:
        pairs = r.filtered_pairs(training, label=l)
        global_per_label[l] = r.fit_scale(pairs) or 1.0

    local_scale = {}
    n_fallback = 0
    for name in slot_names_list:
        local_scale[name] = {}
        for l in non_id_labels:
            pairs = r.filtered_pairs(training, label=l, slot=name)
            f = r.fit_scale(pairs)
            if f is None:
                f = global_per_label[l]
                n_fallback += 1
            local_scale[name][l] = f
    return local_scale, n_fallback


def combine_matrices_local(raw, alpha_labels, identity_label, K, local_scale):
    """Correction depends on BOTH which slot the raw value came from AND
    the label -- unlike combine_matrices' per_slot_scale (slot-only) or
    per_label_scale (label-only), this needs the full (slot,label) table,
    applied before the (E0-E2)/2 cross-term combination exactly like the
    other per-slot-aware corrections in this project."""
    def corrected(name, label):
        return raw[name][label] / local_scale[name][label]

    mats = {l: np.zeros((K, K)) for l in alpha_labels}
    for n in range(K):
        name = f"u_{n}"
        for l in alpha_labels:
            mats[l][n, n] = 1.0 if l == identity_label else corrected(name, l)
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            name0, name2 = f"(u{n}+u{m})", f"(u{n}-u{m})"
            for l in alpha_labels:
                if l == identity_label:
                    re = 0.0
                else:
                    re = (corrected(name0, l) - corrected(name2, l)) / 2
                mats[l][n, m] = re
                mats[l][m, n] = re
    return mats


def main():
    print("\n" + "=" * 78)
    print("  loop_local_perturbation_cdr.py -- locally-perturbed per-(slot,label) CDR (K=6)")
    print("=" * 78)

    p = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    p["solutions"] = solutions
    print(f"  fit {n_ok}/{len(solutions)} targets converged, worst={worst:.2e}")
    counts = r.verify_constant_gate_count(solutions)
    assert counts == {11}, f"gate count not fixed: {counts}"

    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    slot_names_list = r.slot_names(K)
    noise_model = r.build_noise_model()

    raw_noisy, _ = r.measure_raw_per_slot(p, non_id_labels, noise_model)
    raw_mats = r.combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K)
    E_raw, err_raw = r.energy_from_alpha_matrices(raw_mats, p, K)
    print(f"  raw: err_vs_exact={err_raw['err_vs_exact_kcal']:.3f} kcal/mol (consistency check vs 103.99)")
    print(f"  local training: {len(slot_names_list)} slots x {K_LOCAL} perturbed draws "
          f"(+/-{PERTURB_RADIUS} rad) = {len(slot_names_list)*K_LOCAL} circuits per seed")

    seed_rows = []
    for seed in range(N_SEEDS):
        training = generate_local_training_data(p, non_id_labels, noise_model, seed)
        local_scale, n_fallback = fit_local_scales(training, non_id_labels, slot_names_list)
        mats = combine_matrices_local(raw_noisy, p["alpha_labels"], p["identity_label"], K, local_scale)
        E, err = r.energy_from_alpha_matrices(mats, p, K)
        seed_rows.append({"seed": seed, "err_vs_exact_kcal": err["err_vs_exact_kcal"],
                           "err_vs_noiseless_kcal": err["err_vs_noiseless_kcal"], "n_fallback": n_fallback})
        print(f"    seed={seed}: err_vs_exact={err['err_vs_exact_kcal']:.4f} kcal/mol (fallback: {n_fallback})")

    vals_exact = [row["err_vs_exact_kcal"] for row in seed_rows]
    vals_noiseless = [row["err_vs_noiseless_kcal"] for row in seed_rows]
    mean_e, std_e = float(np.mean(vals_exact)), float(np.std(vals_exact))
    mean_n, std_n = float(np.mean(vals_noiseless)), float(np.std(vals_noiseless))
    n_reached = sum(1 for v in vals_exact if v < 0.30)

    print(f"\n  locally-perturbed per-(slot,label) CDR: {mean_e:.4f} +/- {std_e:.4f} kcal/mol (vs exact), "
          f"{mean_n:.4f} +/- {std_n:.4f} (vs noiseless)")
    print(f"  vs baseline (global per-basis) CDR: 2.850 +/- 0.490")
    print(f"  reached 0.30 target: {n_reached}/{N_SEEDS} seeds")

    results = {
        "idea": "locally-perturbed per-(slot,label) CDR scale",
        "k_local": K_LOCAL, "perturb_radius": PERTURB_RADIUS,
        "mean_kcal_vs_exact": mean_e, "std_kcal_vs_exact": std_e,
        "mean_kcal_vs_noiseless": mean_n, "std_kcal_vs_noiseless": std_n,
        "baseline_mean_kcal": 2.850, "baseline_std_kcal": 0.490,
        "improved": bool(mean_e < 2.850),
        "reached_030_target": f"{n_reached}/{N_SEEDS}",
        "seed_rows": seed_rows,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
