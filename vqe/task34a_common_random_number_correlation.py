#!/usr/bin/env python3
"""
task34a_common_random_number_correlation.py -- iteration 34, Task A. The
GATING test for a proposed new estimator architecture: PAIRED,
covariance-aware scalar ZNE on top of PEC+manifold, using COMMON RANDOM
NUMBERS (CRN) across ZNE folds to correlate the PEC Monte Carlo noise
between fold levels, so that a residual bias extrapolation sees far less
noise than an independent-per-fold design would.

WHY THIS TEST EXISTS, and why it's local/free: Task 33E already showed
"more PEC draws everywhere" doesn't fix the champion pipeline's MSE. The
proposed alternative doesn't add draws -- it correlates the draws ALREADY
being taken across fold levels 1/3/5, so that a fold-vs-fold DIFFERENCE
(used for extrapolation) has much of its shared noise cancel. This is
worth nothing to build the full pipeline for unless the correlation is
real and large. Per the proposal's own explicit instruction: "do not run
the full expensive circuit set yet -- run a small proof of concept
first." This is that proof of concept, entirely local (no ionq_simulator
submission, no cost, no wait).

AN IMPORTANT CAVEAT this project's own history makes necessary, disclosed
up front rather than glossed over: this project has already established
(iteration 25, [project_ionq_simulator_cross_submission_drift]) that
IonQ's free simulator resamples a FRESH noise realization on every
separate job submission -- there is no API to share the underlying
PHYSICAL noise process between two circuit executions. What CAN be
shared, because we control it entirely classically, is the PEC TWIRL
CHOICE (which recovery Pauli gets inserted at each noisy gate, and its
sign) -- that is what "common random numbers" means concretely here: the
SAME twirl RNG seed drives gate-recovery sampling at fold=1, fold=3, and
fold=5. Real hardware SHOT noise cannot be shared this way and is NOT
modeled in this test at all (this local test uses exact density-matrix
expectation values, zero shots, by design -- see below). This test can
therefore only tell us whether the CLASSICAL twirl-choice component of
PEC's stochastic estimator is correlated across folds when seeded
identically. That is a real, necessary, but partial question -- it is
Task 34A's job to answer just that question honestly, not to simulate
the full real-hardware pipeline.

METHOD: reuses Task 31B's own validated local noisy-twirl-circuit
simulator UNCHANGED (`noisy_dm`, `sample_twirled_circuit_local` --
density-matrix propagation through the SAME simulated depolarizing
channel the analytic PEC model assumes, for both the raw circuit gates
AND the twirl-correction gates, exactly as real hardware execution would
experience them). NEW: `fold_all_gates` (Task 28d, real circuit-level
gate folding, G->G.G^-1.G, matching every prior ZNE task 27-29) is
applied BEFORE the basis rotation/twirl step, so folds 1/3/5 are built
the same way this project's ZNE work has always built them -- no
invented "compounded probability" shortcut (none exists in this
codebase and none is introduced here).

TWO CONDITIONS, both run, so the correlation attributable to CRN
specifically can be isolated from any correlation the circuits would
show anyway:
  MATCHED:     rng re-seeded IDENTICALLY (same stable_seed) for every
               fold at draw k -- the proposed CRN construction.
  INDEPENDENT: rng seeded differently per (fold, k) -- the status-quo
               baseline (what every prior task in this project has
               effectively done, fold-blind, no shared stream).
If MATCHED and INDEPENDENT show similar correlation, the correlation
isn't coming from the shared-randomness mechanism at all (the circuits
are just similar) and CRN buys nothing. Only a MATCHED-vs-INDEPENDENT
GAP is evidence the mechanism works.

DECISION THRESHOLDS (from the proposal, applied to MATCHED rho, with
INDEPENDENT reported alongside for attribution):
  rho > 0.9         : green, proceed to scalar residual ZNE design
  0.6 <= rho <= 0.9  : yellow, still worth testing
  rho < 0.3          : red, abandon paired-CRN ZNE

Run:
    python vqe/task34a_common_random_number_correlation.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit, fold_all_gates
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from native_stateprep import to_native
from ionq_simulator_binding_curve import stable_seed
from qiskit.quantum_info import DensityMatrix, Operator, Pauli, Statevector
from qiskit_ionq.ionq_gates import GPI2Gate

K = 6
SLOT = "(u0+u1)"
GROUP = ["XYYX", "IYYI"]
LABEL = "IYYI"  # this project's established validation case (Task 31B/32I) -- real known coherent-noise complication
GATE_NAME = "zz"
P2 = 0.0146   # Task 31A consensus
P1 = 0.000119  # Task 31B's own value, unchanged
FOLDS = [1, 3, 5]
M = 256
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task34a_common_random_number_correlation_results.json")


def noisy_dm(qc, p2, p1):
    """Task 31B's own validated function, reused unchanged."""
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
        if op.name == GATE_NAME:
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p2, 2))
        elif op.name in ("gpi", "gpi2"):
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p1, 1))
    return dm


def sample_twirled_circuit_local(base_qc, p2, p1, rng):
    """Task 31B's own validated function, reused unchanged."""
    qc = base_qc.copy_empty_like()
    total_sign = 1
    total_gamma = 1.0
    for instr in base_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name in ("gpi", "gpi2"):
            weights = pec_inverse_weights(p1, 1)
        else:
            continue
        labels = list(weights.keys())
        w = np.array([weights[l] for l in labels])
        gamma_gate = float(np.sum(np.abs(w)))
        probs = np.abs(w) / gamma_gate
        idx = rng.choice(len(labels), p=probs)
        chosen_label, chosen_w = labels[idx], w[idx]
        sign = 1 if chosen_w >= 0 else -1
        total_sign *= sign
        total_gamma *= gamma_gate
        n = len(chosen_label)
        if chosen_label != "I" * n:
            qc.append(Pauli(chosen_label).to_instruction(), qargs)
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign, total_gamma


def one_draw(full_folded, rng, Pmat, p2=P2, p1=P1):
    twirled, sign, gamma = sample_twirled_circuit_local(full_folded, p2, p1, rng)
    dm_t = noisy_dm(twirled, p2, p1)
    m = float(np.real(np.trace(Pmat @ dm_t.data)))
    return sign * gamma * m


def main():
    print("\n" + "=" * 96)
    print("  task34a_common_random_number_correlation.py -- CRN-across-folds gating test (LOCAL, zero shots)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    base = optimized_native_circuit(fixed_solutions[SLOT]["angles"], GATE_NAME)
    combined = effrag_mod.combined_basis_label(GROUP)
    basis_qc = native_basis_change(combined, GATE_NAME)

    full_folded = {}
    for fold in FOLDS:
        folded_base = fold_all_gates(base, fold, GATE_NAME)
        full_folded[fold] = folded_base.compose(basis_qc)
    Pmat = np.asarray(Pauli(LABEL).to_matrix())

    for fold in FOLDS:
        n2q = sum(1 for instr in full_folded[fold].data if instr.operation.name == GATE_NAME)
        n1q = sum(1 for instr in full_folded[fold].data if instr.operation.name in ("gpi", "gpi2"))
        print(f"  fold={fold}: {n2q} ZZ gates, {n1q} 1Q gates in full circuit")

    exact_val = float(np.real(Statevector.from_instruction(full_folded[1]).expectation_value(Pauli(LABEL))))
    print(f"  TRUE exact <{LABEL}> (zero noise, fold=1): {exact_val:.4f}")

    # -- IDEAL CONTROL: at p=0, twirling/PEC must be a no-op regardless of fold or seed --
    print(f"\n  -- IDEAL CONTROL (p2=p1=0, must reproduce exact value, both fold and seed irrelevant) --")
    rng0 = np.random.default_rng(stable_seed("task34a_idealcontrol"))
    ideal_val = one_draw(full_folded[5], rng0, Pmat, p2=0.0, p1=0.0)
    print(f"    fold=5, p=0: {ideal_val:.6f}  (exact={exact_val:.6f}, |diff|={abs(ideal_val-exact_val):.6f})")
    ideal_ok = abs(ideal_val - exact_val) < 1e-6
    print(f"    IDEAL CONTROL: {'PASS' if ideal_ok else 'FAIL -- stop, do not trust anything below'}")
    if not ideal_ok:
        raise RuntimeError("Ideal control failed -- local simulator is not behaving as a no-op at p=0")

    # -- MATCHED (CRN) and INDEPENDENT (control) draws --
    print(f"\n  -- generating {M} draws per fold, MATCHED (shared seed across folds) and INDEPENDENT (control) --")
    E_matched = {fold: np.zeros(M) for fold in FOLDS}
    E_indep = {fold: np.zeros(M) for fold in FOLDS}
    for k in range(M):
        matched_seed = stable_seed("task34a_matched", SLOT, LABEL, k)
        for fold in FOLDS:
            rng = np.random.default_rng(matched_seed)
            E_matched[fold][k] = one_draw(full_folded[fold], rng, Pmat)
        for fold in FOLDS:
            rng = np.random.default_rng(stable_seed("task34a_indep", SLOT, LABEL, fold, k))
            E_indep[fold][k] = one_draw(full_folded[fold], rng, Pmat)
        if (k + 1) % 64 == 0:
            print(f"    {k+1}/{M} draws done")

    def corr(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    print(f"\n  -- CORRELATION: MATCHED (CRN) vs INDEPENDENT (control) --")
    results = {}
    for fold in [3, 5]:
        rho_m = corr(E_matched[1], E_matched[fold])
        rho_i = corr(E_indep[1], E_indep[fold])
        var1_m, varf_m = E_matched[1].var(ddof=1), E_matched[fold].var(ddof=1)
        var_diff_m = (E_matched[1] - E_matched[fold]).var(ddof=1)
        ratio_m = var_diff_m / (var1_m + varf_m)
        var1_i, varf_i = E_indep[1].var(ddof=1), E_indep[fold].var(ddof=1)
        var_diff_i = (E_indep[1] - E_indep[fold]).var(ddof=1)
        ratio_i = var_diff_i / (var1_i + varf_i)
        print(f"    fold 1 vs {fold}:")
        print(f"      MATCHED:     rho={rho_m:+.4f}   Var(E1-E{fold})/[Var(E1)+Var(E{fold})]={ratio_m:.4f} (0.5=uncorrelated, 0=perfect+corr)")
        print(f"      INDEPENDENT: rho={rho_i:+.4f}   Var(E1-E{fold})/[Var(E1)+Var(E{fold})]={ratio_i:.4f}  (control baseline)")
        print(f"      CRN GAP (matched - independent) in rho: {rho_m - rho_i:+.4f}")
        results[f"1_vs_{fold}"] = {
            "rho_matched": rho_m, "rho_independent": rho_i, "crn_gap": rho_m - rho_i,
            "var_ratio_matched": ratio_m, "var_ratio_independent": ratio_i,
            "var_E1_matched": var1_m, "var_Efold_matched": varf_m, "var_diff_matched": var_diff_m,
        }

    # -- VERDICT, per the proposal's own thresholds, applied to the CRN-attributable signal --
    print(f"\n  -- VERDICT (thresholds: green rho>0.9, yellow 0.6-0.9, red rho<0.3) --")
    rho_13_m = results["1_vs_3"]["rho_matched"]
    rho_15_m = results["1_vs_5"]["rho_matched"]
    rho_13_gap = results["1_vs_3"]["crn_gap"]
    rho_15_gap = results["1_vs_5"]["crn_gap"]
    for label, rho in [("1 vs 3", rho_13_m), ("1 vs 5", rho_15_m)]:
        verdict = "GREEN" if rho > 0.9 else ("YELLOW" if rho >= 0.6 else ("RED" if rho < 0.3 else "borderline"))
        print(f"    fold {label}: rho_matched={rho:+.4f} -> {verdict}")
    gap_meaningful = (rho_13_gap > 0.1) or (rho_15_gap > 0.1)
    print(f"\n    CRN mechanism check: matched-vs-independent rho GAP is "
          f"{'MEANINGFUL (CRN is doing real work)' if gap_meaningful else 'SMALL/NONE (any correlation seen is NOT attributable to the shared-seed mechanism -- the circuits/twirl distribution are just inherently similar, CRN buys nothing extra)'}"
          f" (fold1v3 gap={rho_13_gap:+.4f}, fold1v5 gap={rho_15_gap:+.4f})")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "slot": SLOT, "label": LABEL, "folds": FOLDS, "M": M,
            "exact_val": exact_val, "ideal_control_pass": bool(ideal_ok),
            "correlations": results,
            "caveat": "This tests ONLY the classical PEC twirl-choice component, local density-matrix "
                      "expectation (zero shots). Real hardware shot noise cannot be shared across separate "
                      "ionq_simulator job submissions (established iteration 25) and is not modeled here.",
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
