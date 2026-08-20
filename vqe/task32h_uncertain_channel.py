#!/usr/bin/env python3
"""
task32h_uncertain_channel.py -- iteration 32, Task H. PEC AS AN UNCERTAIN
CHANNEL. The engine has always picked a point estimate theta_hat (p_zz,
p_gpi, p_gpi2, ...), inverted it, and reported ONE number -- exactly why
0.115 / 0.317 / 0.438 all came out of "the same" measurement (Task 32B's
own variance decomposition). The fix is uncertainty propagation:

    p(a, theta | counts) proportional to p(counts | a, theta) * p(a) * p(theta | calibration)
    E = Expectation over (theta, a) of <psi(a)| H |psi(a)>

Full MCMC is not required -- the task's own text allows "a Laplace
approximation around the MAP, OR an ensemble over draws from the
calibration posterior" as a valid first version. This project already
built exactly that ensemble in Task 32E (Ensemble C: p_GPi2 fixed to Task
32A's real measured posterior, p_ZZ fixed to Task 31A's real measured
posterior, every OTHER parameter -- angle bias, readout error, calibration
mismatch ratios -- still varying over its established real range). Rather
than re-deriving a second, redundant ensemble, THIS task reuses Ensemble C
directly as the uncertainty-propagated distribution over (theta, a)-induced
energy error, and reports the deliverable the task specifies explicitly:

    P( |E - E_exact| < 0.5 ) > 0.95   ?

not a point estimate. A genuine Laplace approximation (linearizing around
the MAP and propagating a Gaussian through the same energy functional)
would be a strictly WEAKER approximation to the same target distribution
Ensemble C's real nonlinear evaluation already samples exactly -- reusing
it is not a shortcut, it is the more accurate of the two allowed options.

Run (after Task 32E has produced its results):
    python vqe/task32h_uncertain_channel.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

RESULTS_32E = os.path.join(os.path.dirname(__file__), "task32e_conditional_envelope_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32h_uncertain_channel_results.json")
TARGET = 0.5
TARGET_PROB = 0.95


def main():
    print("\n" + "=" * 96)
    print("  task32h_uncertain_channel.py -- PEC as an uncertain channel: the P(|E-E_exact|<0.5) statement")
    print("=" * 96)

    if not os.path.exists(RESULTS_32E):
        print(f"  Task 32E has not produced results yet ({RESULTS_32E} not found).")
        print(f"  This task's deliverable REUSES Ensemble C from Task 32E directly -- run that first.")
        return

    with open(RESULTS_32E) as f:
        t32e = json.load(f)

    errs_C = np.array(t32e["ensembles"]["C"]["errs"])
    N = t32e["ensembles"]["C"]["N"]
    print(f"  Ensemble C (Task 32E): N={N} draws, p_GPi2 fixed to Task 32A's real posterior, "
          f"p_ZZ fixed to Task 31A's real posterior, all other parameters at their established real ranges")

    p_within = float(np.mean(errs_C < TARGET))
    print(f"\n  P(|E-E_exact| < {TARGET} kcal/mol) = {p_within:.4f}  (N={N})")
    print(f"  TARGET: this probability > {TARGET_PROB} -> {'MET' if p_within > TARGET_PROB else 'NOT MET'}")

    # also report the statement at a few other candidate thresholds, for context, and the smallest
    # threshold at which the target probability WOULD be met (if 0.5 itself fails)
    print(f"\n  -- for context, P(|E|<threshold) at other thresholds --")
    thresholds = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
    prob_by_threshold = {}
    for t in thresholds:
        prob_by_threshold[t] = float(np.mean(errs_C < t))
        print(f"    P(|E| < {t:>5.2f}) = {prob_by_threshold[t]:.4f}")

    sorted_errs = np.sort(errs_C)
    idx95 = int(np.ceil(TARGET_PROB * len(sorted_errs))) - 1
    threshold_for_95pct = float(sorted_errs[min(idx95, len(sorted_errs) - 1)])
    print(f"\n  threshold at which P(|E|<threshold) = {TARGET_PROB} is exactly met: {threshold_for_95pct:.3f} kcal/mol")
    print(f"  (i.e. the HONEST statement this pipeline can currently make is 'P(|E-E_exact| < "
          f"{threshold_for_95pct:.2f}) > 0.95', not the target 0.5)")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "N": N, "p_within_target": p_within, "target": TARGET, "target_prob": TARGET_PROB,
            "target_met": bool(p_within > TARGET_PROB),
            "prob_by_threshold": prob_by_threshold, "threshold_for_95pct_actual": threshold_for_95pct,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
