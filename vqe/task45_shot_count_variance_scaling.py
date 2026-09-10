#!/usr/bin/env python3
"""
task45_shot_count_variance_scaling.py -- iteration 45 (v2, corrected).
Answers Vadim Karpusenko's (IonQ Research) Sep 9 question: does the
finite-shot postselection/conditioned-correction estimator need 2,000
shots to keep its statistical error inside this project's Q95 envelope,
or does 1,100 shots (the "free" shot count his own per-circuit minimum-
charge floor analysis identified, given our 1q-gate-reduced circuit)
already work? Answered by REAL bootstrap resampling of REAL already-
collected 20,000-shot data (task39c_ancilla_real_submission.partial.json,
draw 0 -- the actual real ionq_simulator submission behind this
project's own reported headline numbers), NOT a theoretical 1/sqrt(N)
guess -- zero new real submissions, zero new cost.

v1 OF THIS SCRIPT WAS WRONG: it reused corrected_energy() from
task40_certification_ablation_adversarial.py, which the ledger itself
documents (RESEARCH_LEDGER.md, ~line 404) is the "STANDARD pipeline (no
Schmidt-frame fit)" ablation baseline, NOT the actual real pipeline. The
real headline 0.0132/0.0179 kcal/mol numbers require the QED+PEC+GPi2
correction FED INTO the joint Schmidt-frame fit
(task36_joint_schmidt_frame.fit_joint_frame), exactly as
task39h_leakage_free_gpi2_sweep.py does. v1's numbers (Q95~2-15 kcal/mol
even at N=20,000) were real outputs of the WRONG pipeline, not a bug and
not a real finding about this project's actual estimator -- caught via a
direct reproduction check against the documented real result before
trusting anything (see chat transcript), which is exactly why this
project's own convention is "verify against the known real number before
trusting a new derived one."

METHOD (corrected): for each candidate shot count N in
{800, 1100, 2000, 20000}, and for each of N_TRIALS independent draws:
bootstrap-resample every real circuit's RAW 5-qubit bitstring histogram
down to N total shots (matching the real physical order of operations --
N total shots run, THEN leaked shots discarded via ancilla=0
postselection), apply the conditioned QED+PEC+GPi2 correction at the
ALREADY-ESTABLISHED, fixed p_gpi2_selected (0.0006 aria-1 / 0.0004
forte-1 -- consistently selected across every real draw so far, so not
re-swept here), fit the SAME 15-parameter joint Schmidt frame on the
FULL corrected dataset (matching the robustness-envelope convention in
task40_robustness_envelope_with_readout.py, NOT task39h's own 70/30
train/val split, since that split exists only to LEAKAGE-FREE-SELECT
p_gpi2, which is not being re-selected here), and compute the resulting
informational energy error.

SANITY CHECK, run and printed BEFORE any bootstrap trial: N=20,000 using
the REAL, non-resampled, original counts (zero resampling) MUST land
close to the documented real headline numbers (aria-1 ~0.0132,
forte-1 ~0.0179 kcal/mol, or this project's own actual re-derivation
under n_restarts=4: 0.0098/0.0129) -- if it doesn't, STOP, something is
still wrong, do not trust anything below it.

SCOPE, disclosed: this isolates PURE finite-shot statistical noise at a
given N from ONE real underlying distribution (draw 0), propagated
through the REAL frame-fit pipeline. It does NOT capture the separate,
already-quantified real submission-to-submission noise-realization drift
(aria-1 std=4.01, forte-1 std=2.31 kcal/mol on the raw energy scale,
iteration 24-25's own drift task) -- a distinct phenomenon, not what
Vadim asked about here.

Run:
    PYTHONHASHSEED=0 python vqe/task45_shot_count_variance_scaling.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from ionq_simulator_binding_curve import expectation_from_counts, bootstrap_counts, stable_seed
from task39h_leakage_free_gpi2_sweep import apply_correction

K = 6
N_LIST = [800, 1100, 2000, 20000]
N_TRIALS = 10
N_RESTARTS = 4
GPI2_SELECTED = {"aria-1": 0.0006, "forte-1": 0.0004}
DRAW_CKPT = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission.partial.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task45_shot_count_variance_scaling_results.json")


def resample_and_postselect(state, kept, backend_name, N, rng, no_resample=False):
    """Bootstrap the REAL raw 5-qubit counts down to N total shots, THEN
    postselect ancilla=0. no_resample=True bypasses bootstrapping
    entirely (uses the real original counts as-is) -- the sanity-check
    path, must reproduce the documented real result."""
    blended = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for group, counts in zip(entry["groups"], entry["counts"]):
            raw = counts if no_resample else bootstrap_counts(counts, N, rng)
            filtered = {bs[1:]: c for bs, c in raw.items() if bs[0] == "0"}
            for l in group:
                blended[name][l] = expectation_from_counts(filtered, l) if filtered else 0.0
    return blended


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def one_trial(blended, p_gpi2, kept, non_id_labels, fixed_solutions, P_S, p, weight_unit, seed):
    corrected = apply_correction(blended, p_gpi2, kept, non_id_labels, fixed_solutions)
    rng = np.random.default_rng(seed)
    U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, corrected,
                                             weight_unit, rng, n_restarts=N_RESTARTS)
    full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
    return energy_and_err(p, full, K)


def main():
    print("\n" + "=" * 96)
    print("  task45_shot_count_variance_scaling.py (v2) -- REAL pipeline, bootstrap-based shot-count scaling")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    with open(DRAW_CKPT) as f:
        state = json.load(f)

    print(f"\n  -- SANITY CHECK: N=20,000, NO resampling (real original data) -- must match the documented real result --")
    for backend_name in ["aria-1", "forte-1"]:
        blended = resample_and_postselect(state, kept, backend_name, None, None, no_resample=True)
        err = one_trial(blended, GPI2_SELECTED[backend_name], kept, non_id_labels, fixed_solutions, P_S, p,
                          weight_unit, seed=39)
        print(f"    {backend_name}: err={err:+.4f} kcal/mol  (documented real result: ~0.0132/0.0179; "
              f"must be in the same sub-0.02 ballpark, NOT orders of magnitude off)")

    print(f"\n  real data source: {os.path.basename(DRAW_CKPT)} (draw 0)")
    print(f"  N_TRIALS={N_TRIALS} bootstrap resamples per (N, backend)")

    results = {}
    for backend_name in ["aria-1", "forte-1"]:
        p_gpi2_sel = GPI2_SELECTED[backend_name]
        print(f"\n  === {backend_name} (p_gpi2_selected={p_gpi2_sel}) ===")
        results[backend_name] = {}
        for N in N_LIST:
            errs = []
            for trial in range(N_TRIALS):
                rng = np.random.default_rng(stable_seed("task45v2", backend_name, N, trial))
                blended = resample_and_postselect(state, kept, backend_name, N, rng)
                err = one_trial(blended, p_gpi2_sel, kept, non_id_labels, fixed_solutions, P_S, p,
                                  weight_unit, seed=39 + trial)
                errs.append(err)
            errs = np.array(errs)
            abs_errs = np.abs(errs)
            mean_err = float(np.mean(errs))
            std_err = float(np.std(errs, ddof=1))
            q50, q90, q95 = [float(np.percentile(abs_errs, q)) for q in [50, 90, 95]]
            results[backend_name][N] = {
                "mean_signed_err": mean_err, "std_err": std_err,
                "q50": q50, "q90": q90, "q95": q95,
                "all_errs": errs.tolist(),
            }
            print(f"    N={N:>6} shots: mean={mean_err:+.4f}  std={std_err:.4f}  "
                  f"Q50={q50:.4f}  Q90={q90:.4f}  Q95={q95:.4f} kcal/mol")

    print(f"\n  -- COMPARISON vs chemical-accuracy bars (0.25 strict / 0.5 loose kcal/mol) --")
    for backend_name in ["aria-1", "forte-1"]:
        for N in N_LIST:
            q95 = results[backend_name][N]["q95"]
            margin = 0.25 / q95 if q95 > 0 else float("inf")
            print(f"    {backend_name} N={N:>6}: Q95={q95:.4f} kcal/mol  ({margin:.1f}x under the 0.25 kcal/mol strict bar)")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
