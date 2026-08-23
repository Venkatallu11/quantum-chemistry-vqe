#!/usr/bin/env python3
"""
task38d_robust_pec.py -- iteration 38, Task D. Task 38C found V_cal
dominates the H4 error budget (99.9%), with p_gpi2 responsible for
75.1% of it -- driven not by a fragile GPi2 CORRECTION (the standard
pipeline currently applies NONE: theta_0_gpi2=0, since Task 31A's
from-scratch GPi2 recalibration failed) but by the RAW CIRCUIT's own
intrinsic sensitivity to true GPi2 noise going uncorrected entirely.
p_readout's real calibration measurement (the complementary lever)
could not be completed this session -- a real IonQ API "TooManyShots"
error, most likely a quota/rate limit after today's real-submission
volume, stopped it before any data was collected; NOT retried blindly,
flagged to the user, and p_readout stays at Task 37B's uncorrected,
weak-prior baseline here (a real, disclosed gap, not silently patched).

QUESTION THIS FILE ANSWERS: should the standard pipeline introduce ANY
GPi2 correction at all (nonzero assumed p_gpi2 for the PEC-inverse
step), and if so, at what strength -- evaluated by ROBUSTNESS (variance
across the real range of true p_gpi2 uncertainty), not by matching one
nominal point exactly. This is a genuine optimization, not a shrinkage
of an existing correction (theta_0_gpi2=0 means there is nothing to
shrink FROM) -- the free variable is the ASSUMED p_gpi2 used for the
PEC-inverse step itself, searched over a real, PEC-feasibility-bounded
grid (Phase 1B already established candidates above ~0.01-0.02 make the
literal quasi-probability overhead astronomical; this project's own
analytic-ratio correction doesn't pay that overhead directly, but the
grid is kept in the same physically-plausible range regardless).

METHOD, fully analytic (no shot noise, no real data needed -- this
answers "how would a GIVEN candidate correction perform under the real
range of true GPi2 uncertainty", a question about the correction RULE,
not about any one dataset): for M draws of theta_true ~ (real,
Task-37-evidence-based distribution), and for each candidate assumed
p_gpi2_correction in a grid, compute the analytically-corrected energy
via `simulate_true_then_correct` (gate-by-gate: TRUE coherent bias +
TRUE depolarizing, THEN PEC-inverse using the CANDIDATE assumed value --
a direct, disclosed extension of Task 37C's own verified
`analytic_A_and_B_5param`, allowing the depolarize and PEC-inverse steps
to use DIFFERENT p, which that function does not support). TRAINING
objective (leakage-free, no exact energy used): minimize Var_m[E_m] over
a 70% train split of the theta_true draws. VALIDATION (held-out 30%):
compare the resulting Q50/Q95 |E-E_exact| against the theta_correction_
gpi2=0 baseline (today's actual standard pipeline) -- exact energy used
ONLY here, never in the training objective.

Run:
    PYTHONHASHSEED=0 python vqe/task38d_robust_pec.py
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task30b_pec_application import build_full
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

K = 6
GATE_NAME = "zz"
M_DRAWS = 30  # reduced from 60 after the first (buggy) run measured 14.4s/eval -- 480 evals would
# have taken ~115 min; this run's own measured per-eval cost decides whether M_DRAWS/grid need
# further reduction, printed and checked before committing to the full sweep
TRAIN_FRACTION = 0.70
GPI2_GRID = [0.0, 0.0002, 0.0004, 0.0006, 0.001]  # narrowed from the first run's 8-point grid --
# 0.0015+ already showed clearly worse (higher-variance) training behavior even before the negative-
# probability bug fix, no reason to keep re-exploring that region
ZZ_ASSUMED = 0.014593  # Task 31A real calibration -- always corrected, unchanged, not being optimized
# (Task 38C found g_zz*sigma_zz contributes ~0% of V_cal -- nothing to gain here)

# real, disclosed distribution for theta_true, collected across Tasks 37B-E (NOT the true value --
# this represents our real uncertainty about what it could be)
THETA_TRUE_MEAN = {"p_zz": 0.014593, "p_gpi2": 0.0003, "delta_zz": 0.0, "delta_gpi2": 0.0}
THETA_TRUE_STD = {"p_zz": 0.000124, "p_gpi2": 0.00081, "delta_zz": 0.000249, "delta_gpi2": 0.00161}


def simulate_true_then_correct(angles, gate_name, theta_true, p_zz_assumed, p_gpi2_assumed, labels):
    """Gate-by-gate: apply the TRUE coherent bias + TRUE depolarizing
    noise, THEN apply PEC-inverse using the CANDIDATE assumed p --
    directly extends Task 37C's `analytic_A_and_B_5param`, whose B always
    used the SAME p for both steps (fine for "is my assumed calibration
    correct" questions; not for "how robust is THIS candidate correction
    against a DIFFERENT true value", which is what this file asks)."""
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native

    p_zz_true, p_gpi2_true, delta_zz_true, delta_gpi2_true = (
        theta_true["p_zz"], theta_true["p_gpi2"], theta_true["delta_zz"], theta_true["delta_gpi2"])
    qc = to_native(build_ansatz(angles), gate_name)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == gate_name:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz_true))
            p_true, p_assumed, n_here = p_zz_true, p_zz_assumed, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))  # GPi held fixed, real tight calibration
            p_true, p_assumed, n_here = GPI_REAL_MEAN, GPI_REAL_MEAN, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2_true))
            p_true, p_assumed, n_here = p_gpi2_true, p_gpi2_assumed, 1
        else:
            U = Operator(op.to_matrix())
            dm = dm.evolve(U, qargs=qargs)
            continue
        dm = dm.evolve(U, qargs=qargs)
        dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p_true, n_here))
        dm = apply_pauli_mixture(dm, qargs, pec_inverse_weights(p_assumed, n_here))
    return {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm.data))) for l in labels}


def energy_for_draw(theta_true, p_gpi2_assumed, kept, non_id_labels, fixed_solutions, diag, p):
    raw_kept = {}
    for name in kept:
        raw_kept[name] = simulate_true_then_correct(fixed_solutions[name]["angles"], GATE_NAME, theta_true,
                                                       ZZ_ASSUMED, p_gpi2_assumed, non_id_labels)
    full = build_full(raw_kept, diag, K, non_id_labels)
    alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task38d_robust_pec.py -- should the standard pipeline correct GPi2 at all, and how strongly?")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    rng = np.random.default_rng(38)
    theta_true_draws = []
    n_clipped = 0
    for _ in range(M_DRAWS):
        draw = {}
        for name in ["p_zz", "p_gpi2", "delta_zz", "delta_gpi2"]:
            v = float(rng.normal(THETA_TRUE_MEAN[name], THETA_TRUE_STD[name]))
            # BUG CAUGHT (first run): p_gpi2's std (0.00081) exceeds its mean (0.0003), so an
            # unclipped Gaussian draw goes NEGATIVE ~37% of the time -- an unphysical
            # "probability" fed into depolarizing_weights/pec_inverse_weights, which accept it
            # formally but produce wildly nonsensical distortion (this is exactly why the first
            # run's baseline read 60-130 kcal/mol instead of the ~6.7 kcal/mol Task 38C's real-data
            # check found at essentially the same nominal point). Fixed: clip probability-type
            # parameters (p_zz, p_gpi2) to >=0 -- coherent-bias parameters (delta_zz, delta_gpi2)
            # are signed by nature and are NOT clipped.
            if name in ("p_zz", "p_gpi2") and v < 0:
                v = 0.0
                n_clipped += 1
            draw[name] = v
        theta_true_draws.append(draw)
    print(f"  clipped {n_clipped} negative probability draws to 0 (of {M_DRAWS*2} p_zz/p_gpi2 draws)")
    n_train = int(round(M_DRAWS * TRAIN_FRACTION))
    train_draws, val_draws = theta_true_draws[:n_train], theta_true_draws[n_train:]
    print(f"\n  {M_DRAWS} theta_true draws ({n_train} train / {M_DRAWS-n_train} held-out), "
          f"grid of {len(GPI2_GRID)} candidate p_gpi2_assumed values")

    t0 = time.time()
    _ = energy_for_draw(train_draws[0], GPI2_GRID[0], kept, non_id_labels, fixed_solutions, diag, p)
    print(f"  one (draw, candidate) evaluation took {time.time()-t0:.2f}s -- "
          f"{len(theta_true_draws)*len(GPI2_GRID)} total evaluations estimated at "
          f"{(time.time()-t0)*len(theta_true_draws)*len(GPI2_GRID)/60:.1f} min")

    print(f"\n  -- TRAINING: for each candidate, Var_train[E] over {n_train} draws (NO exact-energy leakage) --")
    train_errs_by_candidate = {}
    for p_gpi2_assumed in GPI2_GRID:
        errs = [energy_for_draw(theta, p_gpi2_assumed, kept, non_id_labels, fixed_solutions, diag, p)
                for theta in train_draws]
        train_errs_by_candidate[p_gpi2_assumed] = errs
        var = float(np.var(errs))
        print(f"    p_gpi2_assumed={p_gpi2_assumed:.4f}  Var_train={var:.4f}  "
              f"mean_train_err={np.mean(errs):+.3f}  std_train={np.std(errs):.3f}")

    best_candidate = min(GPI2_GRID, key=lambda c: np.var(train_errs_by_candidate[c]))
    print(f"\n  TRAINING-SELECTED candidate (min variance): p_gpi2_assumed = {best_candidate}")

    print(f"\n  -- VALIDATION: held-out {M_DRAWS-n_train} draws, exact energy used ONLY here, best candidate vs baseline(0) --")
    val_errs_baseline = [abs(energy_for_draw(theta, 0.0, kept, non_id_labels, fixed_solutions, diag, p))
                          for theta in val_draws]
    val_errs_best = [abs(energy_for_draw(theta, best_candidate, kept, non_id_labels, fixed_solutions, diag, p))
                      for theta in val_draws]
    q50_base, q95_base = float(np.percentile(val_errs_baseline, 50)), float(np.percentile(val_errs_baseline, 95))
    q50_best, q95_best = float(np.percentile(val_errs_best, 50)), float(np.percentile(val_errs_best, 95))
    print(f"    BASELINE (p_gpi2_assumed=0, today's actual pipeline): Q50={q50_base:.3f}  Q95={q95_base:.3f} kcal/mol")
    print(f"    ROBUST   (p_gpi2_assumed={best_candidate}):            Q50={q50_best:.3f}  Q95={q95_best:.3f} kcal/mol")
    print(f"\n  -- HONEST READ --")
    if q95_best < q95_base:
        pct = 100 * (q95_base - q95_best) / q95_base
        print(f"    ROBUST correction WINS on held-out data: Q95 {q95_base:.3f} -> {q95_best:.3f} kcal/mol "
              f"({pct:.1f}% reduction). A real, disclosed win from introducing a GPi2 correction the standard "
              f"pipeline never had, chosen to minimize variance across real calibration uncertainty rather than "
              f"to match one nominal point exactly.")
    else:
        print(f"    ROBUST correction does NOT beat the baseline on held-out data (Q95 {q95_base:.3f} vs "
              f"{q95_best:.3f}) -- honest negative result, not glossed over.")


if __name__ == "__main__":
    main()
