#!/usr/bin/env python3
"""
phys_constrained_reconstruction.py — Phase 1 of "physics-constrained
reconstruction": does enforcing PHYSICALITY on the reconstructed reduced
density matrix (Hermitian, positive-semidefinite, trace=1) -- something
NO prior technique in this project has done -- reduce real-hardware
energy error, using ONLY already-collected measurement data (no new
circuits, free)?
============================================================================
THE CORE INSIGHT THIS TESTS: this project's entire pipeline (all 22 prior
iterations) reconstructs the K x K alpha-register matrices ONE MATRIX
ELEMENT AT A TIME -- `combine_matrices()` takes each raw measured
Pauli expectation value m_P (for a given (slot, label) pair) and drops
it DIRECTLY into a matrix entry, with NO constraint that the resulting
per-Pauli matrices are even mutually consistent with a single physical
quantum state, let alone individually Hermitian/PSD/bounded. Under real
shot noise and gate noise, nothing stops (for example) a diagonal
"probability-like" entry from landing outside [-1,1], or the ensemble of
measured matrix elements from being inconsistent with ANY valid density
matrix at all. That is real, currently-discarded information: every
Pauli operator measured for a given slot is a linear functional of the
SAME underlying (noisy) quantum state prepared for that slot, and a
constrained joint fit across all of them should do better than an
independent per-element readout.

METHOD (exactly what phase 1 was asked for): for each of the 36 K=6
Schmidt-basis slots (u_0..u_5, and the 30 (u_n+-u_m) phase-pair states),
the exact classical Schmidt vectors u_0..u_{K-1} give a KNOWN 16xK
isometry U. For every alpha-register Pauli label P actually measured in
this project's real IonQ runs, project it into the K-dim Schmidt
subspace: P_S = U^dagger P U (a KNOWN, exactly-computable KxK Hermitian
matrix -- no estimation needed here, only linear algebra on already-known
classical quantities). Then, per slot, solve the convex problem

    rho_S = argmin_rho  sum_P w_P (Tr(rho P_S) - m_P)^2
    s.t.   rho Hermitian, rho >> 0 (PSD), Tr(rho) = 1

using cvxpy (installed this iteration) -- a genuine, global-optimum
convex semidefinite least-squares problem, not a heuristic. K=6, so this
is a 6x6 SDP: tiny and well-conditioned, solved 36 times per (model,
seed).

WHY THIS ALREADY ENFORCES THE PARTICLE-NUMBER CONSTRAINT, WITH NO EXTRA
TERM NEEDED: `fixed_ansatz.py`'s own verified finding (leakage ~1e-29)
is that EVERY Schmidt vector and EVERY phase-pair target lives EXACTLY
in the alpha register's Hamming-weight-2 sector -- i.e. the columns of U
already span (a subspace of) the physical sector, and nothing outside
it. Any convex combination of pure states expressible in the U-basis
(which is exactly what a PSD, trace-1 K x K rho parametrizes) is
therefore AUTOMATICALLY confined to that same physical sector. Restricting
to the K-dim Schmidt subspace at all IS the particle-number constraint;
no separate penalty term is required, and none is added here -- stated
explicitly rather than silently assumed.

Once rho_S is reconstructed per slot, the "cleaned" per-label expectation
value Tr(rho_S @ P_S) is fed into the EXISTING, UNCHANGED
`qforge.combine_matrices` / `energy_from_alpha_matrices` pipeline -- this
file does not reimplement the EF energy formula, only the upstream
per-slot matrix-element reconstruction step.

DECISION RULE (written BEFORE running, per the standing instruction --
do not move it afterwards): using the REAL aria-1/forte-1 raw baseline
computed on this SAME checkpoint data as the reference point, if
physics-constrained reconstruction takes that baseline error to
approximately 1/3 of its value or lower (the "33 -> ~10 kcal/mol or
better" criterion), continue to Phase 2. If the improvement is marginal
(the "33 -> 32" case, i.e. materially less than a ~2x reduction), ABANDON
this direction and say so plainly -- do not continue to Phase 2 or 3
regardless of how interesting the technique seems in principle.

Run:
    python vqe/phys_constrained_reconstruction.py
"""
import os
import sys
import json
import time
import numpy as np
import cvxpy as cp

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL,
    slot_names, floor_test,
)
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed, expectation_from_counts
from qiskit.quantum_info import Pauli

K = 6
SHOTS = 10_000
N_SEEDS = 8
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "targets_d1.0.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "phys_constrained_reconstruction_results.json")

# DECISION RULE -- fixed before running, see module docstring. Applied to the
# real, same-data raw baseline computed below, not to any external number.
DECISION_RULE_RATIO = 1.0 / 3.0  # "33 -> ~10" is roughly this ratio or better
DECISION_RULE_MIN_IMPROVEMENT = 2.0  # "33 -> 32" (essentially no change) must NOT count as a pass


def load_checkpoint():
    with open(CKPT_PATH) as f:
        return json.load(f)


def build_group_index(groups):
    """label -> index into `groups` (which group's counts to read it from)."""
    idx = {}
    for gi, group in enumerate(groups):
        for l in group:
            idx[l] = gi
    return idx


def build_P_S(alpha_labels, U):
    """U: (16, K) real matrix, columns = Schmidt vectors. Returns dict
    label -> KxK complex Hermitian numpy array = U^dagger P U (exact
    classical linear algebra, no estimation)."""
    out = {}
    for l in alpha_labels:
        Pmat = Pauli(l).to_matrix()
        out[l] = U.conj().T @ Pmat @ U
    return out


def reconstruct_rho_slot(P_S_dict, m_dict, w_dict, K):
    rho = cp.Variable((K, K), hermitian=True)
    terms = []
    for l, m in m_dict.items():
        P = P_S_dict[l]
        pred = cp.real(cp.trace(rho @ P))
        terms.append(w_dict[l] * cp.square(pred - m))
    prob = cp.Problem(cp.Minimize(cp.sum(terms)), [rho >> 0, cp.trace(rho) == 1])
    prob.solve(solver=cp.SCS)
    if rho.value is None:
        raise RuntimeError(f"SDP solve failed, status={prob.status}")
    return rho.value


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def main():
    print("\n" + "=" * 96)
    print("  phys_constrained_reconstruction.py -- Phase 1 (free, existing data only)")
    print("=" * 96)
    print(f"  DECISION RULE (fixed before running): pass requires raw-error reduction to <= "
          f"{DECISION_RULE_RATIO:.3f}x the same-data raw baseline; a reduction of less than "
          f"{DECISION_RULE_MIN_IMPROVEMENT}x counts as 'marginal' and means ABANDON.")

    ck = load_checkpoint()
    print(f"\n  loaded checkpoint: {CKPT_PATH}")
    print(f"  d={ck['d']}, {len(ck['target_names'])} slots, {len(ck['groups'])} groups, "
          f"models={list(ck['counts'].keys())}")

    t0 = time.time()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=ck["d"], K=K)
    print(f"  setup_fragment done, {time.time()-t0:.1f}s")
    exact_mismatch = abs(p["exact_energy"] - ck["exact_energy"]) * HARTREE_TO_KCAL_MOL
    print(f"  cross-check: setup_fragment's exact_energy vs checkpoint's own: "
          f"{exact_mismatch:.2e} kcal/mol (should be ~0)")
    assert exact_mismatch < 1e-6, "checkpoint and fresh setup_fragment disagree on the exact energy -- stop"

    U = np.asarray(p["u_vecs"]).T  # (16, K): columns are the K Schmidt vectors
    assert U.shape == (16, K)
    ortho_err = np.max(np.abs(U.conj().T @ U - np.eye(K)))
    print(f"  Schmidt basis orthonormality check (U^dagger U vs I_K): {ortho_err:.2e} (should be ~0)")
    assert ortho_err < 1e-9

    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    groups = ck["groups"]
    group_idx = build_group_index(groups)
    # every non-identity label must be measured in some group -- verify, don't assume
    missing = [l for l in non_id_labels if l not in group_idx]
    assert not missing, f"labels with no measurement group: {missing}"

    t0 = time.time()
    P_S = build_P_S(alpha_labels, U)
    herm_err = max(float(np.max(np.abs(P_S[l] - P_S[l].conj().T))) for l in alpha_labels)
    print(f"  P_S projections built for {len(alpha_labels)} labels, {time.time()-t0:.1f}s "
          f"(max Hermiticity violation: {herm_err:.2e}, should be ~0)")
    I_S_err = float(np.max(np.abs(P_S[identity_label] - np.eye(K))))
    print(f"  identity-label projection vs I_K: {I_S_err:.2e} (should be ~0, confirms Tr(rho)=1 "
          f"already fixes the identity-label matrix element)")

    names = ck["target_names"]
    assert set(names) == set(slot_names(K))

    results = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_by_slot = ck["counts"][model]  # slot_name -> list of group count-dicts

        errs_raw_exact, errs_raw_noiseless = [], []
        errs_phys_exact, errs_phys_noiseless = [], []
        n_sdp_fail = 0
        t0 = time.time()
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("physcon", model, seed))

            raw_baseline = {name: {} for name in names}
            raw_phys = {name: {} for name in names}
            for name in names:
                resampled_groups = [bootstrap_counts(counts_by_slot[name][gi], SHOTS, rng)
                                     for gi in range(len(groups))]
                m_dict, w_dict = {}, {}
                for l in non_id_labels:
                    gi = group_idx[l]
                    counts = resampled_groups[gi]
                    m = expectation_from_counts(counts, l)
                    m_dict[l] = m
                    raw_baseline[name][l] = m
                    total = sum(counts.values())
                    var = max(1 - m ** 2, 1e-4) / max(total, 1)
                    w_dict[l] = 1.0 / var

                try:
                    rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                except RuntimeError:
                    n_sdp_fail += 1
                    raw_phys[name] = dict(m_dict)  # fall back to raw for this slot if solve fails
                    continue
                for l in non_id_labels:
                    raw_phys[name][l] = float(np.real(np.trace(rho_slot @ P_S[l])))

            _, err_raw = energy_and_err(p, raw_baseline, K)
            _, err_phys = energy_and_err(p, raw_phys, K)
            errs_raw_exact.append(err_raw["err_vs_exact_kcal"])
            errs_raw_noiseless.append(err_raw["err_vs_noiseless_kcal"])
            errs_phys_exact.append(err_phys["err_vs_exact_kcal"])
            errs_phys_noiseless.append(err_phys["err_vs_noiseless_kcal"])

        t_elapsed = time.time() - t0
        results[model] = {
            "raw_vs_exact_mean": float(np.mean(errs_raw_exact)), "raw_vs_exact_std": float(np.std(errs_raw_exact)),
            "raw_vs_noiseless_mean": float(np.mean(errs_raw_noiseless)), "raw_vs_noiseless_std": float(np.std(errs_raw_noiseless)),
            "phys_vs_exact_mean": float(np.mean(errs_phys_exact)), "phys_vs_exact_std": float(np.std(errs_phys_exact)),
            "phys_vs_noiseless_mean": float(np.mean(errs_phys_noiseless)), "phys_vs_noiseless_std": float(np.std(errs_phys_noiseless)),
            "n_sdp_fail": n_sdp_fail, "wall_clock_s": t_elapsed,
        }
        print(f"\n  {model} ({t_elapsed:.1f}s, {n_sdp_fail} SDP solve failures out of {N_SEEDS*len(names)}):")
        print(f"    RAW (same data):    vs_exact={results[model]['raw_vs_exact_mean']:.3f}+/-{results[model]['raw_vs_exact_std']:.3f}  "
              f"vs_noiseless={results[model]['raw_vs_noiseless_mean']:.3f}+/-{results[model]['raw_vs_noiseless_std']:.3f} kcal/mol")
        print(f"    PHYS-CONSTRAINED:   vs_exact={results[model]['phys_vs_exact_mean']:.3f}+/-{results[model]['phys_vs_exact_std']:.3f}  "
              f"vs_noiseless={results[model]['phys_vs_noiseless_mean']:.3f}+/-{results[model]['phys_vs_noiseless_std']:.3f} kcal/mol")

    print(f"\n  -- DECISION RULE APPLIED (per model, vs_exact, real noise models only) --")
    verdicts = {}
    for model in ["aria-1", "forte-1"]:
        raw = results[model]["raw_vs_exact_mean"]
        phys = results[model]["phys_vs_exact_mean"]
        ratio = phys / raw if raw > 1e-9 else float("nan")
        reduction = raw / phys if phys > 1e-9 else float("inf")
        passed = (ratio <= DECISION_RULE_RATIO) or (reduction >= DECISION_RULE_MIN_IMPROVEMENT)
        verdicts[model] = {"raw": raw, "phys": phys, "ratio": ratio, "reduction_factor": reduction, "passed": bool(passed)}
        print(f"    {model}: raw={raw:.2f} -> phys={phys:.2f} kcal/mol "
              f"(reduction factor {reduction:.2f}x)  {'PASS' if passed else 'FAIL/MARGINAL'}")

    any_pass = any(v["passed"] for v in verdicts.values())
    both_pass = all(v["passed"] for v in verdicts.values())
    print(f"\n  OVERALL: {'CONTINUE to Phase 2' if both_pass else ('MIXED -- see per-model verdicts, judgment call needed' if any_pass else 'ABANDON -- decision rule not met on either real noise model')}")

    out = {
        "decision_rule": {"pass_ratio": DECISION_RULE_RATIO, "min_reduction_factor": DECISION_RULE_MIN_IMPROVEMENT},
        "results": results, "verdicts": verdicts,
        "any_pass": bool(any_pass), "both_pass": bool(both_pass),
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return out


if __name__ == "__main__":
    main()
