#!/usr/bin/env python3
"""
task30a_manifold_residual_diagnosis.py -- iteration 30, Task A. BLOCKING.
task29c's mean_pure_state_fit_residual is huge (7,725 ideal; 153,664
aria-1; 157,576 forte-1) -- noiseless data should fit a pure state in the
manifold nearly perfectly. Two competing explanations, tested directly
rather than assumed:
  (a) the residual metric is unnormalised/mis-scaled and the fits are
      actually good
  (b) the data genuinely does not lie near the manifold, and the 5-param
      constraint FORCE-PROJECTS it -- the 11x ideal improvement would then
      be the constraint dragging the estimate toward the KNOWN ground
      state, not extracting information from measurements
============================================================================
PART 1 -- explain the residual magnitude. The fit objective is
sum_l w_l (pred-m)^2 with w_l = 1/var, var = max(1-m^2,1e-4)/total. For
near-deterministic labels (|m|->1) and total=100,000 shots, var can be as
low as 1e-4/100,000 = 1e-9, i.e. w_l up to ~1e9 -- a handful of such
labels can dominate the raw SSE even when every individual |pred-m| is
tiny. Reports the UNWEIGHTED mean/median/max |pred-m| alongside the raw
weighted SSE to settle whether this is (a).

PART 2 -- THE DECISIVE TEST FOR (b). Fits the manifold to DELIBERATELY
WRONG data two ways:
  (i)  pure random noise, uniformly drawn in each label's own plausible
       range -- zero real physical information, unrelated to H4 at all
  (ii) the REAL measured H4 data with labels SHUFFLED across the group
       (same numbers, wrong physical meaning -- generically inconsistent
       with any single quantum state)
...using BOTH of production's initialization scheme (target-anchored: the
KNOWN ideal target + small perturbations, exactly what task29c does) AND
a fully UNINFORMED initialization (uniform random points on the sphere,
no knowledge of the target at all). If the target-anchored scheme returns
an energy near the true H4 ground state EVEN ON RANDOM NOISE, but the
uninformed scheme does not, the estimator's apparent skill is coming from
its initialization, not its data -- disqualified, same failure class as
iteration 2's answer-injection.

Run:
    python vqe/task30a_manifold_residual_diagnosis.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K, ckpt_path_for_K
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import stable_seed, bootstrap_counts, expectation_from_counts
from task29c_manifold_estimator import target_coeff_vector, build_full_from_a

K = 6
SHOTS = 100_000
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task30a_manifold_residual_diagnosis_results.json")


def fit_pure_state_verbose(P_S, m_dict, w_dict, K, v0, seed, n_restarts=8):
    labels = list(m_dict.keys())
    Ps = [P_S[l] for l in labels]
    ms = np.array([m_dict[l] for l in labels])
    ws = np.array([w_dict[l] for l in labels])

    def objective(v):
        a = v / np.linalg.norm(v)
        pred = np.array([float(np.real(a @ P @ a)) for P in Ps])
        return float(np.sum(ws * (pred - ms) ** 2))

    rng = np.random.default_rng(seed)
    inits = [v0] + [v0 + rng.normal(0, scale, K) for scale in [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2]]
    best_val, best_v = float("inf"), None
    for v_init in inits[:n_restarts]:
        res = minimize(objective, v_init, method="L-BFGS-B")
        if res.fun < best_val:
            best_val, best_v = res.fun, res.x
    a_hat = best_v / np.linalg.norm(best_v)
    pred = np.array([float(np.real(a_hat @ P @ a_hat)) for P in Ps])
    unweighted = np.abs(pred - ms)
    return a_hat, best_val, unweighted, ws


def fit_pure_state_uninformed(P_S, m_dict, w_dict, K, seed, n_restarts=8):
    """Same optimization, but EVERY restart is a fully random point on the
    sphere -- no informed target used anywhere, unlike production."""
    labels = list(m_dict.keys())
    Ps = [P_S[l] for l in labels]
    ms = np.array([m_dict[l] for l in labels])
    ws = np.array([w_dict[l] for l in labels])

    def objective(v):
        a = v / np.linalg.norm(v)
        pred = np.array([float(np.real(a @ P @ a)) for P in Ps])
        return float(np.sum(ws * (pred - ms) ** 2))

    rng = np.random.default_rng(seed)
    best_val, best_v = float("inf"), None
    for _ in range(n_restarts):
        v_init = rng.normal(0, 1, K)  # uninformed: random direction on the sphere
        res = minimize(objective, v_init, method="L-BFGS-B")
        if res.fun < best_val:
            best_val, best_v = res.fun, res.x
    a_hat = best_v / np.linalg.norm(best_v)
    return a_hat, best_val


def main():
    print("\n" + "=" * 96)
    print("  task30a_manifold_residual_diagnosis.py -- BLOCKING: is the manifold estimator real?")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(ckpt_path_for_K(K)) as f:
        ck = json.load(f)

    # ==================================================================
    # PART 1 -- explain the residual magnitude on REAL ideal fold=1 data
    # ==================================================================
    print("\n  -- PART 1: is the huge residual just weighting, or a bad fit? (REAL ideal fold=1 data) --")
    model, fold = "ideal", 1
    key = f"{fold}|{model}"
    tags = ck["tags"][key]
    counts_list = ck["counts"][key]
    per_name = {}
    for (name, group), counts in zip(tags, counts_list):
        per_name.setdefault(name, {}).setdefault(tuple(group), counts)

    all_unweighted = []
    all_weights = []
    per_slot_report = {}
    for name in kept[:6]:  # a representative sample -- 6 slots is enough to characterize the pattern
        seed_vals = {l: [] for l in non_id_labels if any(l in g for g in per_name[name])}
        seed_vars = {l: [] for l in seed_vals}
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task30a", fold, model, name, seed))
            for group_t, counts in per_name[name].items():
                resampled = bootstrap_counts(counts, SHOTS, rng)
                total = sum(resampled.values())
                for l in group_t:
                    m = expectation_from_counts(resampled, l)
                    seed_vals[l].append(m)
                    var = max(1 - m ** 2, 1e-4) / max(total, 1)
                    seed_vars[l].append(var)
        m_dict = {l: float(np.mean(v)) for l, v in seed_vals.items()}
        w_dict = {l: 1.0 / max(float(np.mean(seed_vars[l])), 1e-6) for l in seed_vals}
        v0 = target_coeff_vector(name, K)
        a_hat, resid, unweighted, ws = fit_pure_state_verbose(P_S, m_dict, w_dict, K, v0,
                                                                seed=stable_seed("task30a_fit", name))
        all_unweighted.extend(unweighted.tolist())
        all_weights.extend(ws.tolist())
        per_slot_report[name] = {
            "weighted_sse": resid, "mean_unweighted_abs_err": float(np.mean(unweighted)),
            "max_unweighted_abs_err": float(np.max(unweighted)), "max_weight": float(np.max(ws)),
            "min_weight": float(np.min(ws)),
        }
        print(f"    {name}: weighted_SSE={resid:.1f}  mean|pred-m|={np.mean(unweighted):.5f}  "
              f"max|pred-m|={np.max(unweighted):.5f}  weight range=[{np.min(ws):.2e}, {np.max(ws):.2e}]")

    print(f"\n    ACROSS {len(kept[:6])} SAMPLE SLOTS: mean unweighted |pred-m|={np.mean(all_unweighted):.5f}, "
          f"median={np.median(all_unweighted):.5f}, max={np.max(all_unweighted):.5f}")
    print(f"    weight range across all labels: [{np.min(all_weights):.2e}, {np.max(all_weights):.2e}] "
          f"(ratio max/min={np.max(all_weights)/np.min(all_weights):.2e})")
    explanation_a_confirmed = np.mean(all_unweighted) < 0.02 and np.max(all_weights) > 1e6
    print(f"\n    VERDICT PART 1: {'(a) CONFIRMED' if explanation_a_confirmed else 'NOT CONFIRMED'} -- "
          f"{'the fit IS good (unweighted errors small); huge SSE is inverse-variance weight blowup for near-deterministic labels' if explanation_a_confirmed else 'unweighted errors are NOT small -- residual may reflect a genuinely poor fit'}")

    # ==================================================================
    # PART 2 -- the decisive adversarial test for (b)
    # ==================================================================
    print("\n  -- PART 2: fit the manifold to DELIBERATELY WRONG data --")
    exact_energy = p["exact_energy"]
    print(f"    true H4 exact energy: {exact_energy:.6f} Ha")

    # -- build a REAL raw fold=1 dataset for reference (all kept slots) --
    def build_real_m_dicts(model_, fold_):
        key_ = f"{fold_}|{model_}"
        tags_ = ck["tags"][key_]
        counts_list_ = ck["counts"][key_]
        per_name_ = {}
        for (name_, group_), counts_ in zip(tags_, counts_list_):
            per_name_.setdefault(name_, {}).setdefault(tuple(group_), counts_)
        out = {}
        for name_ in kept:
            seed_vals_ = {l: [] for l in non_id_labels if any(l in g for g in per_name_[name_])}
            seed_vars_ = {l: [] for l in seed_vals_}
            for seed_ in range(N_SEEDS):
                rng_ = np.random.default_rng(stable_seed("task30a_real", fold_, model_, name_, seed_))
                for group_t_, counts_ in per_name_[name_].items():
                    resampled_ = bootstrap_counts(counts_, SHOTS, rng_)
                    total_ = sum(resampled_.values())
                    for l_ in group_t_:
                        m_ = expectation_from_counts(resampled_, l_)
                        seed_vals_[l_].append(m_)
                        var_ = max(1 - m_ ** 2, 1e-4) / max(total_, 1)
                        seed_vars_[l_].append(var_)
            out[name_] = {
                "m_dict": {l: float(np.mean(v)) for l, v in seed_vals_.items()},
                "w_dict": {l: 1.0 / max(float(np.mean(seed_vars_[l])), 1e-6) for l in seed_vals_},
            }
        return out

    real_data = build_real_m_dicts("ideal", 1)

    def energy_from_a_dict(a_dict):
        full = build_full_from_a(a_dict, P_S, diag, K, non_id_labels)
        alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
        E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                              exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
        return E, errs["err_vs_exact_kcal"]

    adversarial_results = {}

    # -- test (i): PURE RANDOM NOISE, unrelated to H4 at all --
    print("\n    TEST (i): fit to PURE RANDOM NOISE (uniform, no physical information)")
    rng_noise = np.random.default_rng(999)
    for init_style in ["target_anchored (production)", "uninformed (fixed)"]:
        a_dict = {}
        for name in kept:
            labels_here = list(real_data[name]["m_dict"].keys())
            m_dict_fake = {l: float(rng_noise.uniform(-0.9, 0.9)) for l in labels_here}
            w_dict_here = real_data[name]["w_dict"]
            v0 = target_coeff_vector(name, K)
            if init_style.startswith("target"):
                a_hat, _, _, _ = fit_pure_state_verbose(P_S, m_dict_fake, w_dict_here, K, v0,
                                                         seed=stable_seed("task30a_adv1", name))
            else:
                a_hat, _ = fit_pure_state_uninformed(P_S, m_dict_fake, w_dict_here, K,
                                                      seed=stable_seed("task30a_adv1u", name))
            a_dict[name] = a_hat
        E, err_kcal = energy_from_a_dict(a_dict)
        print(f"      [{init_style}]: E={E:.4f} Ha, err_vs_exact={err_kcal:.2f} kcal/mol "
              f"(raw fold=1 err for reference: 1.46 kcal/mol)")
        adversarial_results[f"random_noise__{init_style}"] = {"E": E, "err_vs_exact_kcal": err_kcal}

    # -- test (ii): REAL data, labels SHUFFLED within each slot --
    print("\n    TEST (ii): fit to REAL data with labels SHUFFLED (same numbers, wrong meaning)")
    rng_shuf = np.random.default_rng(1234)
    for init_style in ["target_anchored (production)", "uninformed (fixed)"]:
        a_dict = {}
        for name in kept:
            m_dict_real = real_data[name]["m_dict"]
            w_dict_real = real_data[name]["w_dict"]
            labels_here = list(m_dict_real.keys())
            values_here = [m_dict_real[l] for l in labels_here]
            shuffled_values = list(values_here)
            rng_shuf.shuffle(shuffled_values)
            m_dict_shuf = dict(zip(labels_here, shuffled_values))
            v0 = target_coeff_vector(name, K)
            if init_style.startswith("target"):
                a_hat, _, _, _ = fit_pure_state_verbose(P_S, m_dict_shuf, w_dict_real, K, v0,
                                                         seed=stable_seed("task30a_adv2", name))
            else:
                a_hat, _ = fit_pure_state_uninformed(P_S, m_dict_shuf, w_dict_real, K,
                                                      seed=stable_seed("task30a_adv2u", name))
            a_dict[name] = a_hat
        E, err_kcal = energy_from_a_dict(a_dict)
        print(f"      [{init_style}]: E={E:.4f} Ha, err_vs_exact={err_kcal:.2f} kcal/mol")
        adversarial_results[f"shuffled_real__{init_style}"] = {"E": E, "err_vs_exact_kcal": err_kcal}

    # -- control: REAL data, NOT shuffled, both init styles (should recover the known 0.13 result) --
    print("\n    CONTROL: fit to REAL, UNSHUFFLED data (should reproduce task29c's ~0.13 kcal/mol)")
    for init_style in ["target_anchored (production)", "uninformed (fixed)"]:
        a_dict = {}
        for name in kept:
            m_dict_real = real_data[name]["m_dict"]
            w_dict_real = real_data[name]["w_dict"]
            v0 = target_coeff_vector(name, K)
            if init_style.startswith("target"):
                a_hat, _, _, _ = fit_pure_state_verbose(P_S, m_dict_real, w_dict_real, K, v0,
                                                         seed=stable_seed("task30a_ctrl", name))
            else:
                a_hat, _ = fit_pure_state_uninformed(P_S, m_dict_real, w_dict_real, K,
                                                      seed=stable_seed("task30a_ctrlu", name))
            a_dict[name] = a_hat
        E, err_kcal = energy_from_a_dict(a_dict)
        print(f"      [{init_style}]: E={E:.4f} Ha, err_vs_exact={err_kcal:.2f} kcal/mol")
        adversarial_results[f"real_unshuffled__{init_style}"] = {"E": E, "err_vs_exact_kcal": err_kcal}

    print("\n" + "=" * 96)
    print("  VERDICT PART 2:")
    rn_target = adversarial_results["random_noise__target_anchored (production)"]["err_vs_exact_kcal"]
    rn_uninf = adversarial_results["random_noise__uninformed (fixed)"]["err_vs_exact_kcal"]
    sh_target = adversarial_results["shuffled_real__target_anchored (production)"]["err_vs_exact_kcal"]
    sh_uninf = adversarial_results["shuffled_real__uninformed (fixed)"]["err_vs_exact_kcal"]
    injects_on_noise = rn_target < 5.0
    injects_on_shuffled = sh_target < 5.0
    print(f"    Fit to PURE RANDOM NOISE, target-anchored init: {rn_target:.2f} kcal/mol "
          f"({'SUSPICIOUSLY CLOSE TO EXACT -- INJECTION SUSPECTED' if injects_on_noise else 'far from exact, as expected of real noise'})")
    print(f"    Fit to PURE RANDOM NOISE, uninformed init:      {rn_uninf:.2f} kcal/mol")
    print(f"    Fit to SHUFFLED real data, target-anchored init: {sh_target:.2f} kcal/mol "
          f"({'SUSPICIOUSLY CLOSE -- INJECTION SUSPECTED' if injects_on_shuffled else 'far from exact, as expected'})")
    print(f"    Fit to SHUFFLED real data, uninformed init:      {sh_uninf:.2f} kcal/mol")
    disqualified = injects_on_noise or injects_on_shuffled
    print(f"\n  FINAL VERDICT: {'DISQUALIFIED -- the estimator injects the answer via initialization' if disqualified else 'CLEARED -- the estimator does not inject the answer; genuinely bad data gives genuinely bad results'}")
    print("=" * 96 + "\n")

    results = {
        "part1_unweighted_stats": {
            "mean_abs_err": float(np.mean(all_unweighted)), "median_abs_err": float(np.median(all_unweighted)),
            "max_abs_err": float(np.max(all_unweighted)), "weight_min": float(np.min(all_weights)),
            "weight_max": float(np.max(all_weights)), "explanation_a_confirmed": bool(explanation_a_confirmed),
        },
        "per_slot_part1": per_slot_report,
        "part2_adversarial": adversarial_results,
        "disqualified": bool(disqualified),
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
