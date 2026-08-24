#!/usr/bin/env python3
"""
task39h_leakage_free_gpi2_sweep.py -- iteration 39, Task H. Redoes Task
G's p_gpi2_assumed sweep WITHOUT the target-leakage bug Task G was
caught and disclosed for: Task G selected candidates by directly
minimizing err_vs_exact_kcal (the known exact H4 energy) with no
held-out split -- exactly the failure mode that disqualified iteration
2's CDR result. Comparing against the IDEAL (U=I) per-label target
would ALSO be leakage (U=I directly encodes the exact, classically-known
ground truth Schmidt decomposition, just at finer granularity than the
aggregate energy) -- so this task does NOT do that either.

THE CORRECT, LEAKAGE-FREE CRITERION, matching Task 36/37C's own
established methodology exactly: for each candidate p_gpi2_assumed
(held fixed, not fit), apply Task 39E's conditioned correction to the
real postselected data, then FIT THE 15-PARAMETER SCHMIDT FRAME FREELY
via least-squares (`task36_joint_schmidt_frame.fit_joint_frame`, U_hat
is a FITTED unknown, never fixed to the ideal/exact value) on a
STRATIFIED 70% TRAINING split of the (slot, label) residuals. Select
the p_gpi2_assumed that gives the BEST TRAINING chi2/dof -- internal
data-self-consistency, no reference to the true state or energy
anywhere in the selection. ONLY AFTER a candidate is chosen: evaluate
its held-out 30% validation chi2/dof (generalization check) and, purely
informationally, the resulting energy vs the known exact value -- this
last number is NEVER used to pick anything, exactly Task 38C/38D's own
"informational only" convention.

Run:
    PYTHONHASHSEED=0 python vqe/task39h_leakage_free_gpi2_sweep.py
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
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from ionq_simulator_binding_curve import expectation_from_counts

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
GPI2_GRID = [0.0, 0.0004, 0.0006, 0.0007, 0.0008]  # covers the baseline plus the region Task G's
# (unvalidated) exploratory sweep flagged as interesting for each backend
VAL_FRACTION = 0.30
NEW_CKPT = os.path.join(os.path.dirname(__file__),
                          os.environ.get("TASK39_CKPT_NAME", "task39c_ancilla_real_submission.partial.json"))

_CACHE = {}


def apply_correction(raw_by_slot_label, p_gpi2, kept, non_id_labels, fixed_solutions):
    corrected = {name: {} for name in kept}
    for name in kept:
        key = (name, round(p_gpi2, 6))
        if key not in _CACHE:
            A, B, retA, retB = analytic_A_and_B_conditioned(fixed_solutions[name]["angles"], GATE_NAME,
                                                               ZZ_ASSUMED, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0,
                                                               non_id_labels)
            _CACHE[key] = (A, B)
        A, B = _CACHE[key]
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            corrected[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    return corrected


def load_postselected(kept, backend_name):
    with open(NEW_CKPT) as f:
        state = json.load(f)
    blended = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for group, counts in zip(entry["groups"], entry["counts"]):
            filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
            for l in group:
                blended[name][l] = expectation_from_counts(filtered, l)
    return blended


def split_train_val(blended, kept, val_fraction, seed):
    rng = np.random.default_rng(seed)
    train, val = {name: {} for name in kept}, {name: {} for name in kept}
    for name in kept:
        labels = sorted(blended[name].keys())
        n_val = max(1, round(len(labels) * val_fraction)) if len(labels) > 1 else 0
        val_labels = set(rng.choice(labels, size=n_val, replace=False)) if n_val else set()
        for l in labels:
            (val if l in val_labels else train)[name][l] = blended[name][l]
    return train, val


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def chi2dof_on_holdout(U_hat, P_S, K, kept, holdout_by_slot, weight_unit):
    from task36_joint_schmidt_frame import slot_vector
    total, n = 0.0, 0
    for name in kept:
        v = slot_vector(U_hat, name, K)
        for l, y in holdout_by_slot[name].items():
            pred = float(np.real(v @ P_S[l] @ v))
            total += (pred - y) ** 2
            n += 1
    return total / max(1, n), n


def run_backend(backend_name, kept, non_id_labels, fixed_solutions, diag, P_S, p):
    print(f"\n  === {backend_name} ===")
    real_blended = load_postselected(kept, backend_name)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    train_by_candidate = {}
    for p_gpi2 in GPI2_GRID:
        corrected = apply_correction(real_blended, p_gpi2, kept, non_id_labels, fixed_solutions)
        train, val = split_train_val(corrected, kept, VAL_FRACTION, seed=39)
        rng = np.random.default_rng(39)
        U_hat, cost, chi2dof_train = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, train,
                                                        weight_unit, rng, n_restarts=4)
        chi2dof_val, n_val = chi2dof_on_holdout(U_hat, P_S, K, kept, val, weight_unit)
        train_by_candidate[p_gpi2] = (U_hat, chi2dof_train, chi2dof_val, n_val)
        print(f"    p_gpi2_assumed={p_gpi2:.4f}   TRAIN chi2/dof={chi2dof_train:.5f}   "
              f"VAL chi2/dof={chi2dof_val:.5f} (n={n_val})")

    best_p_gpi2 = min(GPI2_GRID, key=lambda c: train_by_candidate[c][1])
    U_best, chi2dof_train_best, chi2dof_val_best, n_val_best = train_by_candidate[best_p_gpi2]
    print(f"\n    TRAINING-SELECTED (min TRAIN chi2/dof, no exact-energy leakage): p_gpi2_assumed={best_p_gpi2}")
    print(f"    held-out VALIDATION chi2/dof at that candidate: {chi2dof_val_best:.5f}")

    full = build_full_from_frame(U_best, P_S, K, non_id_labels, kept)
    E, err = energy_and_err(p, full, K)
    print(f"    INFORMATIONAL ONLY (not used in selection): resulting energy err vs exact = {err:+.4f} kcal/mol")
    return best_p_gpi2, chi2dof_train_best, chi2dof_val_best, err


def main():
    print("\n" + "=" * 96)
    print("  task39h_leakage_free_gpi2_sweep.py -- LEAKAGE-FREE GPi2 selection (training chi2/dof, no exact energy)")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    results = {}
    for backend_name in ["aria-1", "forte-1"]:
        results[backend_name] = run_backend(backend_name, kept, non_id_labels, fixed_solutions, diag, P_S, p)

    print(f"\n  -- SUMMARY (leakage-free selection) --")
    for backend_name, (best_p, chi2_tr, chi2_val, err) in results.items():
        print(f"    {backend_name}: selected p_gpi2_assumed={best_p}, TRAIN chi2/dof={chi2_tr:.5f}, "
              f"VAL chi2/dof={chi2_val:.5f}, informational err={err:+.4f} kcal/mol")
    print(f"\n  Chemical accuracy bars for reference: 0.25 kcal/mol (internal target), 0.5 kcal/mol (looser hardware bar).")


if __name__ == "__main__":
    main()
