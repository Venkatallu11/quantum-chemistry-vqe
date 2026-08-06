#!/usr/bin/env python3
"""
loop_affine_cdr.py — iteration 1 of the below-0.30-kcal/mol search.
============================================================================
IDEA: fit exact ~= a*noisy + b per label (affine), instead of the current
noisy ~= f*exact through-the-origin scale.

EXPECTATION (stated before running, not after): the noise model here is a
pure multiplicative depolarizing channel -- verified elsewhere in this
project (cdr_mitigation.py's docstring) that qiskit-aer's depolarizing_error
shrinks any non-identity Pauli expectation by EXACTLY (1-param), no additive
offset. That means the TRUE relationship between noisy and exact IS already
linear through the origin, which is exactly what the current per-basis fit
assumes. An affine fit adds a free intercept parameter with nothing true to
fit -- I expect it to be roughly NEUTRAL (absorbs a little finite-sample
correlation, possibly reducing variance slightly) or slightly WORSE (one
more parameter estimated from the same data, more overfitting risk),
not a big win. Worth 10 minutes to confirm this rather than assume it.

Identity label is still NEVER touched (forced to 1.0, no affine fit either
-- an intercept b != 0 on the identity would corrupt it just as badly as a
scale did in the original bug).

Run:
    python vqe/loop_affine_cdr.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import rank6_symmetry_vd as r

K = r.K
N_SEEDS = r.N_SEEDS
N_TRAIN_PER_SLOT = r.N_TRAIN_PER_SLOT
LOW_SIGNAL_CUTOFF = r.LOW_SIGNAL_CUTOFF

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "loop_affine_cdr_results.json")


def fit_affine(pairs):
    """exact ~= a*noisy + b, fit as noisy = a*exact + b (same regression,
    then inverted at correction time: exact_est = (noisy - b) / a).
    Returns (a, b) or None if too little data."""
    if len(pairs) < 4:
        return None
    exact = np.array([e for e, _ in pairs])
    noisy = np.array([n for _, n in pairs])
    A = np.vstack([exact, np.ones_like(exact)]).T
    result, *_ = np.linalg.lstsq(A, noisy, rcond=None)
    a, b = result
    if abs(a) < 1e-6:
        return None
    return float(a), float(b)


def fit_all_scales_affine(training, non_id_labels):
    per_basis = {}
    fallback = []
    for l in non_id_labels:
        pairs = r.filtered_pairs(training, label=l)
        ab = fit_affine(pairs)
        if ab is None:
            ab = (1.0, 0.0)
            fallback.append(l)
        per_basis[l] = ab
    return per_basis, fallback


def combine_matrices_affine(raw, alpha_labels, identity_label, K, per_label_affine):
    def corrected(name, label):
        v = raw[name][label]
        a, b = per_label_affine[label]
        return (v - b) / a

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
    print("  loop_affine_cdr.py -- affine per-basis CDR fit (K=6)")
    print("=" * 78)

    p = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    p["solutions"] = solutions
    print(f"  fit {n_ok}/{len(solutions)} targets converged, worst={worst:.2e}")
    counts = r.verify_constant_gate_count(solutions)
    assert counts == {11}, f"gate count not fixed: {counts}"

    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    noise_model = r.build_noise_model()

    raw_noisy, _ = r.measure_raw_per_slot(p, non_id_labels, noise_model)
    raw_mats = r.combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K)
    E_raw, err_raw = r.energy_from_alpha_matrices(raw_mats, p, K)
    print(f"  raw: err_vs_exact={err_raw['err_vs_exact_kcal']:.3f} kcal/mol (consistency check vs 103.99)")

    seed_rows = []
    for seed in range(N_SEEDS):
        training = r.generate_training_data(p, non_id_labels, noise_model, N_TRAIN_PER_SLOT, seed)
        per_label_affine, fallback = fit_all_scales_affine(training, non_id_labels)
        mats = combine_matrices_affine(raw_noisy, p["alpha_labels"], p["identity_label"], K, per_label_affine)
        E, err = r.energy_from_alpha_matrices(mats, p, K)
        seed_rows.append({"seed": seed, "err_vs_exact_kcal": err["err_vs_exact_kcal"],
                           "err_vs_noiseless_kcal": err["err_vs_noiseless_kcal"],
                           "n_fallback": len(fallback)})
        print(f"    seed={seed}: err_vs_exact={err['err_vs_exact_kcal']:.3f} kcal/mol "
              f"(fallback labels: {len(fallback)})")

    vals_exact = [row["err_vs_exact_kcal"] for row in seed_rows]
    vals_noiseless = [row["err_vs_noiseless_kcal"] for row in seed_rows]
    mean_e, std_e = float(np.mean(vals_exact)), float(np.std(vals_exact))
    mean_n, std_n = float(np.mean(vals_noiseless)), float(np.std(vals_noiseless))
    n_reached = sum(1 for v in vals_exact if v < 0.30)

    print(f"\n  affine per-basis CDR: {mean_e:.3f} +/- {std_e:.3f} kcal/mol (vs exact), "
          f"{mean_n:.3f} +/- {std_n:.3f} (vs noiseless)")
    print(f"  vs baseline (pure-scale) per-basis: 2.850 +/- 0.490")
    print(f"  reached 0.30 target: {n_reached}/{N_SEEDS} seeds")

    results = {
        "idea": "affine per-basis CDR fit (exact = a*noisy + b)",
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
