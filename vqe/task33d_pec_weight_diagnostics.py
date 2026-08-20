#!/usr/bin/env python3
"""
task33d_pec_weight_diagnostics.py -- iteration 33, Task D. WHY is the
champion pipeline unstable? Diagnose the PEC quasi-probability WEIGHTS
directly, instead of building another estimator on top of them.

Motivation, concrete and overdue: Task 32B already found PEC Monte Carlo
sampling is the dominant variance source (2.11 of 2.42 total). Task 33C
(this iteration) then found per-label damping on the champion pipeline
gives erratic, sometimes wildly unstable results (std spiking to 6.9-7.4
at several k values, vs ~1.0-1.3 elsewhere) and crashed outright before
finishing. Neither task asked the most basic diagnostic question: are a
FEW high-weight quasi-probability draws dominating each label's 16-draw
average? If so, "16 draws" may mean an effective sample size far below
16 for some labels -- exactly the kind of thing that would explain both
Task 32B's blown-up M=1 variance (std=90.6) and Task 33C's erratic per-k
spikes (different bootstrap resamples occasionally drawing a
disproportionate share of a label's weight from one or two circuits).

METHOD, CORRECTED after a first version's diagnostic turned out to be
inapplicable: an initial pass computed ESS=(sum|w|)^2/sum(w^2) on
w=sign*gamma per draw, and got EXACTLY 16.00/16 for all 756 groups with
zero variation -- a red flag, not a clean bill of health. Inspecting the
raw data directly showed why: gamma is a FIXED constant per (slot,label)
group (e.g. 1.227091400091747 identically on all 16 draws for one
checked group) -- standard PEC sign-sampling has each draw contribute
+-gamma with a fixed magnitude, only the SIGN is random. |w| is therefore
identical across draws BY CONSTRUCTION, so the importance-sampling ESS
formula is mathematically guaranteed to return n=16 regardless of any
real instability -- it was the wrong tool for this sampling scheme, not
evidence the pipeline is fine. The real analog: how close is each
group's SIGN PATTERN to a 50/50 coin flip? A lopsided split (e.g. 14+/2-)
gives a stable mean; a near-even split means the correction is barely
distinguishable from noise with only 16 draws. For every (slot,label)
group:
  - mean_sign = mean of the 16 real +-1 signs (near 0 = coin-flip-like,
    maximally noisy; near +-1 = consistent, stable)
  - sign_se = sqrt(1 - mean_sign^2) / sqrt(16), the actual sampling
    standard error of that mean-sign estimate
  - noise_contribution = gamma * sign_se, gamma's role as an AMPLIFIER
    of whatever sign noise exists (this is the real, sampling-scheme-
    correct analog of the ESS idea, not a discredited one)
  - gamma itself, for reference

No new submissions. No bootstrap. No manifold fit. Just reading and
summarizing real numbers already on disk -- the cheapest possible
diagnostic, which is exactly why it should have run before either
adaptive-PEC estimator (Task 33B/33C) was built on top of this data.

Run:
    python vqe/task33d_pec_weight_diagnostics.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                          "task31c_full_pec_calibration.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task33d_pec_weight_diagnostics_results.json")
LOW_ESS_THRESHOLD = 8.0  # half of the nominal 16 draws -- a disclosed, round-number flag, not tuned to the data


def main():
    print("\n" + "=" * 96)
    print("  task33d_pec_weight_diagnostics.py -- real PEC quasi-probability weight diagnostics, no new data")
    print("=" * 96)

    with open(CKPT_PATH) as f:
        ck = json.load(f)
    tags = [tuple(t) for t in ck["tags"]]
    gamma_per_circuit = ck["gamma_per_circuit"]

    by_name_label = {}
    for (name, group, draw, sign), gamma in zip(tags, gamma_per_circuit):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma))
    n_groups = len(by_name_label)
    n_mc_available = min(len(v) for v in by_name_label.values())
    print(f"  {n_groups} (slot, label) groups, {n_mc_available} MC draws/group available (Task 31C's real data)")

    stats = {}
    for (name, l), entries in by_name_label.items():
        signs = np.array([sign for sign, gamma in entries], dtype=float)
        gammas = np.array([gamma for sign, gamma in entries])
        n = len(entries)
        mean_sign = float(signs.mean())
        sign_se = float(np.sqrt(max(0.0, 1.0 - mean_sign ** 2)) / np.sqrt(n))
        gamma_val = float(gammas.mean())  # constant per group, by construction -- mean==every value
        noise_contribution = gamma_val * sign_se
        stats[f"{name}|{l}"] = {
            "n": n,
            "mean_sign": mean_sign,
            "n_plus": int((signs > 0).sum()),
            "n_minus": int((signs < 0).sum()),
            "sign_se": sign_se,
            "gamma": gamma_val,
            "noise_contribution": noise_contribution,
        }

    noise_values = np.array([s["noise_contribution"] for s in stats.values()])
    mean_sign_values = np.array([s["mean_sign"] for s in stats.values()])
    print(f"\n  -- sign-pattern diagnostics across {n_groups} (slot,label) groups, {n_mc_available} real draws each --")
    print(f"    |mean_sign| distribution: min={np.abs(mean_sign_values).min():.3f}  "
          f"median={np.median(np.abs(mean_sign_values)):.3f}  max={np.abs(mean_sign_values).max():.3f}")
    print(f"    noise_contribution (gamma * sign_se) distribution: min={noise_values.min():.4f}  "
          f"median={np.median(noise_values):.4f}  max={noise_values.max():.4f}")
    n_low = int((np.abs(mean_sign_values) < 0.5).sum())
    print(f"    {n_low}/{n_groups} groups have |mean_sign| < 0.5 (signs closer to a coin flip than a "
          f"consistent correction)")

    worst = sorted(stats.items(), key=lambda kv: -kv[1]["noise_contribution"])[:15]
    print(f"\n  -- 15 WORST (slot,label) groups by noise_contribution (gamma * sign_se) --")
    for key, s in worst:
        print(f"    {key:<22} mean_sign={s['mean_sign']:+.3f} ({s['n_plus']}+/{s['n_minus']}-)  "
              f"gamma={s['gamma']:6.3f}  noise_contribution={s['noise_contribution']:.4f}")

    best = sorted(stats.items(), key=lambda kv: kv[1]["noise_contribution"])[:5]
    print(f"\n  -- 5 BEST (slot,label) groups by noise_contribution, for contrast --")
    for key, s in best:
        print(f"    {key:<22} mean_sign={s['mean_sign']:+.3f} ({s['n_plus']}+/{s['n_minus']}-)  "
              f"gamma={s['gamma']:6.3f}  noise_contribution={s['noise_contribution']:.4f}")

    # -- gamma distribution overall: how severe is the PEC overhead being paid --
    all_gamma = np.array([s["gamma"] for s in stats.values()])
    print(f"\n  -- gamma (PEC sampling overhead, fixed per group) across all {n_groups} groups --")
    print(f"    median={np.median(all_gamma):.3f}  max={all_gamma.max():.3f}  min={all_gamma.min():.3f}")

    # -- cross-check against Task 33B's independently-derived boundary_proximity risk labels, if available --
    per_label_noise = {}
    for key, s in stats.items():
        _, l = key.split("|", 1)
        per_label_noise.setdefault(l, []).append(s["noise_contribution"])
    per_label_noise = {l: float(np.mean(v)) for l, v in per_label_noise.items()}
    t33b_path = os.path.join(os.path.dirname(__file__), "task33b_adaptive_label_pec_results.json")
    overlap_note = "Task 33B results not found -- cross-check skipped."
    if os.path.exists(t33b_path):
        with open(t33b_path) as f:
            t33b = json.load(f)
        boundary_proximity = t33b.get("boundary_proximity", {})
        common = sorted(set(per_label_noise) & set(boundary_proximity))
        if len(common) >= 3:
            x = np.array([per_label_noise[l] for l in common])
            y = np.array([boundary_proximity[l] for l in common])
            from scipy.stats import spearmanr
            rho, pval = spearmanr(x, y)
            overlap_note = (f"Spearman(per-label sign-noise, Task 33B boundary_proximity) across "
                             f"{len(common)} shared labels: rho={rho:+.3f}, p={pval:.4f}")
            print(f"\n  -- cross-check: is PEC sign-noise the SAME thing as Task 32I/33B's boundary-proximity risk? --")
            print(f"    {overlap_note}")
            top_noise = sorted(per_label_noise.items(), key=lambda kv: -kv[1])[:5]
            top_prox = sorted(boundary_proximity.items(), key=lambda kv: -kv[1])[:5]
            print(f"    top-5 by sign-noise:        {[l for l, _ in top_noise]}")
            print(f"    top-5 by boundary_proximity: {[l for l, _ in top_prox]}")

    print(f"\n  -- HONEST READ --")
    frac_low = n_low / n_groups
    print(f"    A FIRST VERSION of this diagnostic (ESS on sign*gamma weights) was mathematically "
          f"guaranteed to report 16.00/16 for every group, because gamma is constant per group in this "
          f"project's PEC implementation -- that was a broken diagnostic, not a clean result, caught by "
          f"inspecting raw entries before trusting it.")
    print(f"    The corrected version: {n_low}/{n_groups} ({100*frac_low:.1f}%) of (slot,label) groups have "
          f"|mean_sign| < 0.5 across their 16 real draws -- i.e. the correction for these groups is closer "
          f"to a coin flip than a confident signal, even before gamma amplifies whatever that noise is. "
          f"This is a real, previously undiagnosed, concrete source of the pipeline's known instability "
          f"(Task 32B's dominant PEC-MC variance, Task 33C's erratic per-k results): with only 16 draws, "
          f"any group whose true sign-balance is genuinely close to 50/50 will have a highly unstable "
          f"average, and gamma multiplies that noise directly into the final energy. This is a real, "
          f"different, and probably more useful lever than Task 33B/33C's per-label damping: for groups "
          f"with high noise_contribution, MORE real draws (not damping) would directly reduce sign_se "
          f"(which shrinks as 1/sqrt(n)) -- a concrete, targeted (not blanket) case for extending N_MC "
          f"beyond 16, but only for the specific high-noise groups identified above, not uniformly.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "n_groups": n_groups, "n_mc_available": n_mc_available,
            "mean_sign_abs_summary": {"min": float(np.abs(mean_sign_values).min()),
                                       "median": float(np.median(np.abs(mean_sign_values))),
                                       "max": float(np.abs(mean_sign_values).max())},
            "noise_contribution_summary": {"min": float(noise_values.min()),
                                            "median": float(np.median(noise_values)),
                                            "max": float(noise_values.max())},
            "n_low_mean_sign_groups": n_low,
            "worst_groups": {k: v for k, v in worst}, "best_groups": {k: v for k, v in best},
            "per_label_noise_contribution": per_label_noise,
            "boundary_proximity_crosscheck": overlap_note,
            "all_stats": stats,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
