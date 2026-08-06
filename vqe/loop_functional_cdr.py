#!/usr/bin/env python3
"""
loop_functional_cdr.py — iteration 3 of the below-0.30-kcal/mol search.
============================================================================
CONTEXT: iteration 2 (locally-perturbed per-(slot,label) CDR) was
DISQUALIFIED -- see RESEARCH_LEDGER.md. Its free parameter (perturbation
radius) had no floor: shrinking it toward 0 just re-simulated the target
classically and reported that as the answer, which does not scale to
register sizes CDR is actually for.

IDEA HERE IS DIFFERENT IN A WAY THAT MATTERS: iteration 2's training points
were chosen BASED ON the target (perturbations of that specific target's own
angles) -- more targets or finer precision required MORE classical
simulation, scaling with target proximity. This idea fits ONE GLOBAL
FUNCTION per label from a FIXED-SIZE pool of GLOBALLY RANDOM training
angles (same cost profile, same methodology, as the original per-basis
CDR -- nothing here is chosen to be "close" to any target), then EVALUATES
that function at each target's own (already-known) angles -- a plain
function evaluation, no new classical simulation per target. The training
cost does not grow with how many targets there are or how precisely you
want each one corrected.

MOTIVATION (measured before building this, not assumed): the per-label
noisy/exact ratio is angle-dependent (established in iteration 2's
diagnosis). A constant per-basis scale is the poorest possible model of an
angle-dependent function (a 1-term / degree-0 fit). Fit a richer model
instead: noisy ~= f(angles)*exact, with f(angles) = coeffs . features(angles),
features = [1, cos(th_i), sin(th_i) for i in 0..4] (11 features -- a
physically motivated basis, since backward-propagating a Pauli through a
rotation gate generates trig functions of that gate's angle).

QUICK CHECK BEFORE BUILDING THIS FILE: on 3 sample labels, this 11-feature
linear model reduced training-residual std by 1.05x-1.86x vs the constant
model (ZZII: 1.05x, YZYZ: 1.86x, IIIZ: 1.41x) -- real, but modest. Realistic
expectation stated up front: this alone probably does NOT reach 0.30 kcal/mol;
recorded honestly either way, and used to decide whether a richer feature
set (or a different method entirely) is the next iteration.

MANDATORY FLOOR TEST: the free parameter here is N_TRAIN (total global
training draws, not proximity to any target). Swept 100/200/400/800 below
-- if the actual 8-seed ENERGY ERROR keeps falling without bound as N_TRAIN
grows, that would itself be suspicious (more global data legitimately
converging is expected and fine UP TO a floor set by the noise model's own
non-linearity/model mismatch -- but genuinely unbounded improvement with no
plateau would need investigating before trusting it, unlike iteration 2's
proximity parameter which had an obvious mechanism for cheating).

Run:
    python vqe/loop_functional_cdr.py
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
N_TRAIN_DEFAULT = 300
MIN_POINTS_FOR_FIT = 15  # >= n_features(11) with margin, else fall back to constant scale

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "loop_functional_cdr_results.json")


def angle_features(angles):
    feats = [1.0]
    for th in angles:
        feats.append(np.cos(th))
        feats.append(np.sin(th))
    return np.array(feats)  # 11 features


def generate_global_training_data(p, non_id_labels, noise_model, n_total, seed):
    """GLOBAL random draws -- NOT organized by slot, NOT chosen near any
    target. Same methodology/cost profile as the original per-basis CDR."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_total):
        angles = rng.uniform(-np.pi, np.pi, 5).tolist()
        exact_vals = r.exact_labels(angles, non_id_labels)
        dm = r.density_matrix_4q(angles, noise_model)
        noisy_vals, _ = r.measure_labels(angles, non_id_labels, dm)
        for l in non_id_labels:
            rows.append({"label": l, "exact": exact_vals[l], "noisy": noisy_vals[l], "angles": angles})
    return rows


def fit_functional_scales(training, non_id_labels):
    per_label_coeffs = {}
    n_fallback = 0
    for l in non_id_labels:
        rows = [row for row in training if row["label"] == l and abs(row["exact"]) >= LOW_SIGNAL_CUTOFF]
        if len(rows) < MIN_POINTS_FOR_FIT:
            per_label_coeffs[l] = None  # signals constant-scale fallback
            n_fallback += 1
            continue
        X = np.array([angle_features(row["angles"]) * row["exact"] for row in rows])
        y = np.array([row["noisy"] for row in rows])
        exacts = np.array([row["exact"] for row in rows])
        coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
        # constant-scale fallback value bundled in for use when the functional
        # prediction is degenerate (near-zero) at a specific target
        const_scale = float(np.sum(exacts * y) / np.sum(exacts * exacts))
        per_label_coeffs[l] = (coeffs, const_scale)
    return per_label_coeffs, n_fallback


def predict_scale(per_label_coeffs, label, angles):
    entry = per_label_coeffs[label]
    if entry is None:
        return 1.0
    coeffs, const_scale = entry
    f = float(angle_features(angles) @ coeffs)
    if abs(f) < 0.05:  # degenerate prediction -- fall back rather than blow up
        return const_scale
    return f


def combine_matrices_functional(raw, alpha_labels, identity_label, K, per_label_coeffs, target_angles_by_slot):
    def corrected(name, label):
        v = raw[name][label]
        f = predict_scale(per_label_coeffs, label, target_angles_by_slot[name])
        return v / f

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


def run_sweep(p, non_id_labels, raw_noisy, noise_model, target_angles_by_slot, n_train, n_seeds):
    seed_rows = []
    for seed in range(n_seeds):
        training = generate_global_training_data(p, non_id_labels, noise_model, n_train, seed)
        per_label_coeffs, n_fallback = fit_functional_scales(training, non_id_labels)
        mats = combine_matrices_functional(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                            per_label_coeffs, target_angles_by_slot)
        E, err = r.energy_from_alpha_matrices(mats, p, K)
        seed_rows.append({"seed": seed, "err_vs_exact_kcal": err["err_vs_exact_kcal"],
                           "err_vs_noiseless_kcal": err["err_vs_noiseless_kcal"], "n_fallback": n_fallback})
    return seed_rows


def main():
    print("\n" + "=" * 78)
    print("  loop_functional_cdr.py -- global functional (angle-feature) CDR fit (K=6)")
    print("=" * 78)

    p = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    p["solutions"] = solutions
    print(f"  fit {n_ok}/{len(solutions)} targets converged, worst={worst:.2e}")
    counts = r.verify_constant_gate_count(solutions)
    assert counts == {11}, f"gate count not fixed: {counts}"

    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    noise_model = r.build_noise_model()
    target_angles_by_slot = {name: sol["angles"] for name, sol in solutions.items()}

    raw_noisy, _ = r.measure_raw_per_slot(p, non_id_labels, noise_model)
    raw_mats = r.combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K)
    E_raw, err_raw = r.energy_from_alpha_matrices(raw_mats, p, K)
    print(f"  raw: err_vs_exact={err_raw['err_vs_exact_kcal']:.3f} kcal/mol (consistency check vs 103.99)")

    print(f"\n  -- MANDATORY FLOOR TEST: sweeping N_TRAIN (global draws, not proximity) --")
    floor_test = {}
    for n_train in (100, 200, 400, 800):
        rows = run_sweep(p, non_id_labels, raw_noisy, noise_model, target_angles_by_slot, n_train, n_seeds=2)
        vals = [row["err_vs_exact_kcal"] for row in rows]
        mean_v = float(np.mean(vals))
        floor_test[n_train] = mean_v
        print(f"    N_TRAIN={n_train}: mean(2 seeds)={mean_v:.4f} kcal/mol")

    ordered = sorted(floor_test.items())
    diffs = [ordered[i + 1][1] - ordered[i][1] for i in range(len(ordered) - 1)]
    monotone_no_floor = all(d < -0.01 for d in diffs)  # still meaningfully falling at every step
    print(f"  floor test verdict: {'NO FLOOR DETECTED YET (needs more N_TRAIN or is suspicious)' if monotone_no_floor else 'shows diminishing returns / floor -- legitimate'}")

    print(f"\n  -- Full {N_SEEDS}-seed sweep at N_TRAIN={N_TRAIN_DEFAULT} --")
    seed_rows = run_sweep(p, non_id_labels, raw_noisy, noise_model, target_angles_by_slot,
                           N_TRAIN_DEFAULT, N_SEEDS)
    for row in seed_rows:
        print(f"    seed={row['seed']}: err_vs_exact={row['err_vs_exact_kcal']:.4f} kcal/mol "
              f"(fallback labels: {row['n_fallback']})")

    vals_exact = [row["err_vs_exact_kcal"] for row in seed_rows]
    vals_noiseless = [row["err_vs_noiseless_kcal"] for row in seed_rows]
    mean_e, std_e = float(np.mean(vals_exact)), float(np.std(vals_exact))
    mean_n, std_n = float(np.mean(vals_noiseless)), float(np.std(vals_noiseless))
    n_reached = sum(1 for v in vals_exact if v < 0.30)

    print(f"\n  functional (angle-feature) CDR: {mean_e:.4f} +/- {std_e:.4f} kcal/mol (vs exact), "
          f"{mean_n:.4f} +/- {std_n:.4f} (vs noiseless)")
    print(f"  vs baseline (constant per-basis) CDR: 2.850 +/- 0.490")
    print(f"  reached 0.30 target: {n_reached}/{N_SEEDS} seeds")

    results = {
        "idea": "global functional (angle-feature, degree-1 trig) per-label CDR scale",
        "n_features": 11, "n_train_default": N_TRAIN_DEFAULT,
        "floor_test_n_train_sweep": floor_test,
        "floor_test_verdict": "no_floor_detected" if monotone_no_floor else "diminishing_returns",
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
