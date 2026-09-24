#!/usr/bin/env python3
"""
task66_h4_gc_frame_fit.py -- iteration 66. Answers a direct question:
can the LOW-CIRCUIT GC measurement design (84 circuits/backend, task59)
be combined with the joint Schmidt-frame fit (task36) -- the thing that
actually gives this project's real, documented chemical-accuracy result
-- instead of the no-frame-fit reconstruction task59 originally used?

WHY THIS WASN'T ALREADY ANSWERED: task47 tested exactly this pairing,
but only ANALYTICALLY (exact-population noise model, zero finite-shot
statistics) -- result was GC+frame-fit at ~2x QWC+frame-fit's error,
both far under the 0.25 kcal/mol bar (worst draw 0.0172 vs 0.0090
kcal/mol). Separately, task59 found that REAL finite-shot statistics hit
the GC design hard when reconstructing WITHOUT the frame fit (denser
measurement groups dilute per-label shot counts, confirmed via an exact/
noiseless reconstruction check that ruled out a code bug). Whether the
joint frame fit -- which pools information across all 21 slots at once,
collapsing 105 independent directions to 15 shared degrees of freedom --
is robust enough to absorb that same real shot-noise dilution is a
genuinely separate, untested question. This answers it with REAL data.

DATA SOURCE: task59's own real checkpoint
(task59_h4_gc_no_frame_fit.partial.json) -- real raw bitstring counts
from real submissions to IonQ's free ionq_simulator, all 21 slots x 3
backends (ideal/aria-1/forte-1), through the actual 84-circuit GC
measurement design, 20,000 shots/circuit. NO NEW SUBMISSION, NO NEW
COST -- this script only changes what happens AFTER measurement:
corrected per-label Pauli expectations (task59's own QED+conditioned-PEC
+GPi2 correction, byte-for-byte reused) are fed into task36's
fit_joint_frame instead of task59's direct no-frame reconstruction.

Run:
    PYTHONHASHSEED=0 python vqe/task66_h4_gc_frame_fit.py
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
from general_commuting_measurements import build_general_commuting_measurement_plan, expectations_from_counts
from task59_h4_gc_no_frame_fit import CKPT_PATH, BACKENDS, GPI2_SELECTED, get_ab

K = 6
N_RESTARTS = 4
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task66_h4_gc_frame_fit_results.json")


def corrected_raw_from_real_counts(state, backend_name, kept, diagonalizers, fixed_solutions, non_id_labels):
    """Byte-for-byte the same QED postselection + conditioned PEC+GPi2
    correction task59.analyze() uses -- only the destination differs
    (fed into the frame fit here, instead of task59's own direct
    no-frame reconstruction)."""
    blended = {name: {} for name in kept}
    accept_fracs = []
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for dg, counts in zip(diagonalizers, entry["counts"]):
            filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
            total = sum(counts.values())
            accept_fracs.append(sum(filtered.values()) / total if total else 0.0)
            exp = expectations_from_counts(filtered, dg)
            blended[name].update(exp)

    if backend_name == "ideal":
        raw_kept = {name: {l: max(-1.0, min(1.0, v)) for l, v in blended[name].items()} for name in kept}
    else:
        p_gpi2_sel = GPI2_SELECTED[backend_name]
        raw_kept = {}
        for name in kept:
            raw_kept[name] = {}
            A, B = get_ab(name, p_gpi2_sel, fixed_solutions, non_id_labels)
            for l, m_raw in blended[name].items():
                ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
                raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    return raw_kept, float(np.mean(accept_fracs))


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task66_h4_gc_frame_fit.py -- REAL GC (84-circuit) measurement + joint Schmidt-frame fit")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    if not os.path.exists(CKPT_PATH):
        print(f"  ERROR: {CKPT_PATH} not found -- run task59's submit step first.")
        return

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)

    with open(CKPT_PATH) as f:
        state = json.load(f)
    missing = [bn for bn in BACKENDS if any(f"{bn}|{name}" not in state["done"] for name in kept)]
    if missing:
        print(f"  ERROR: real checkpoint incomplete for backends {missing} -- cannot proceed on partial real data.")
        return

    print(f"  H4 exact_energy={p['exact_energy']:.6f} Ha")
    print(f"  data: REAL, already-collected task59 checkpoint, 21 slots x 3 backends, "
          f"84 circuits/backend, 20,000 shots/circuit -- zero new cost")

    results = {}
    for backend_name in BACKENDS:
        raw_kept, mean_accept = corrected_raw_from_real_counts(
            state, backend_name, kept, diagonalizers, fixed_solutions, non_id_labels)

        rng = np.random.default_rng(66)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_kept,
                                                 weight_unit, rng, n_restarts=N_RESTARTS)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        E, err = energy_and_err(p, full, K)

        print(f"    {backend_name}: mean_accept={mean_accept:.4f}  chi2/dof={chi2dof:.3f}  "
              f"E={E:.6f} Ha  err_vs_exact={err:+.4f} kcal/mol")
        results[backend_name] = {"E": E, "err_kcal": err, "mean_accept": mean_accept, "chi2_dof": float(chi2dof)}

    no_frame = json.load(open(os.path.join(os.path.dirname(__file__), "task59_h4_gc_no_frame_fit_results.json")))
    print(f"\n  -- COMPARISON: same real GC circuit data, frame-fit vs task59's own no-frame reconstruction --")
    for bn in BACKENDS:
        print(f"    {bn}: frame-fit={results[bn]['err_kcal']:+.4f}  "
              f"no-frame(task59)={no_frame[bn]['err_kcal']:+.4f} kcal/mol")

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
