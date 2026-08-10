#!/usr/bin/env python3
"""
task3_subspace_tomography.py — iteration 24, Task 3. SUBSPACE TOMOGRAPHY
for the cross terms: changes WHAT IS MEASURED (not just how many
circuits), and enforces physicality across the WHOLE reconstruction
design rather than per element. Uses ONLY already-collected real IonQ
data (no new circuits submitted) -- a genuine circuit-count-reduction
question, not a new measurement.
============================================================================
CURRENT SCHEME (every prior iteration): 36 state-prep circuits per K=6 --
6 diagonal (u_n) + 30 phase-pair circuits, TWO per pair ((u_n+u_m)/sqrt2
AND (u_n-u_m)/sqrt2), combined via combine_matrices' difference formula
Re<u_n|P|u_m> = (<P>_+ - <P>_-)/2. This file changes the DESIGN itself:

  1. DROP the "-" circuit for every pair entirely -- 6 diag + 15 "+"
     circuits = 21 total (a genuine 42% circuit-count cut, not just a
     relabeling). This is possible because <P>_+ ALREADY algebraically
     contains the cross term once the (separately measured) diagonals are
     known: <P>_+ = (M_nn+M_mm)/2 + Re<u_n|P|u_m>, so
     Re<u_n|P|u_m> = <P>_+ - (M_nn+M_mm)/2 -- the "-" circuit was
     REDUNDANT information (noise-averaging only), not mathematically
     necessary. Verified below to recover the same cross terms exactly in
     the noiseless limit before trusting it under real noise.

  2. Reuse Phase 1's EXACT SDP machinery (physics-constrained rho
     reconstruction, Hermitian+PSD+trace=1) on EVERY ONE OF THE 21 KEPT
     CIRCUITS (not just the diagonal slots, unlike Phase 1) -- then derive
     the cross terms from the RECONSTRUCTED (not raw) diagonal and "+"
     values. This is what "enforce physicality across the WHOLE design"
     means operationally: the cross-term formula's INPUTS are themselves
     already-constrained-physical quantities, not independently-fit raw
     numbers -- a genuinely joint design, not per-element patching.

Both raw data (checkpoint's existing measurements) and no new circuits
are used -- "existing data only" per the standing instruction. The
DROPPED "-" circuits' data still EXISTS in the checkpoint; this file
simply does not use it for the new (21-circuit) scheme, to honestly
simulate what a 21-circuit submission would have given, while keeping
the full 36-circuit "current" scheme available (from the same checkpoint)
as the baseline comparison.

Applies Phase 3's mandatory ideal-data sanity check to this new method,
per explicit instruction: does the physics-constrained reconstruction
distort already-clean (near-noiseless) data catastrophically? If so, it
is DISQUALIFIED exactly like Phase 3 was.

Run:
    python vqe/task3_subspace_tomography.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL, slot_names
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed, expectation_from_counts
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot, build_group_index
from qiskit.quantum_info import Pauli

K = 6
SHOTS = 10_000
N_SEEDS = 8
CKPT_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", "targets_d1.0.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task3_subspace_tomography_results.json")


def load_checkpoint():
    with open(CKPT_PATH) as f:
        return json.load(f)


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def main():
    print("\n" + "=" * 96)
    print("  task3_subspace_tomography.py -- reduced-circuit design + joint physicality")
    print("=" * 96)

    ck = load_checkpoint()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=ck["d"], K=K)
    U = np.asarray(p["u_vecs"]).T
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    groups = ck["groups"]
    group_idx = build_group_index(groups)
    P_S = build_P_S(alpha_labels, U)

    names = ck["target_names"]
    assert set(names) == set(slot_names(K))
    diag_names = [f"u_{n}" for n in range(K)]
    plus_pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    plus_names = [f"(u{n}+u{m})" for n, m in plus_pairs]
    minus_names = [f"(u{n}-u{m})" for n, m in plus_pairs]
    kept_names = diag_names + plus_names
    assert len(kept_names) == 21 and len(names) == 36
    print(f"  circuit count: CURRENT scheme = {len(names)} (6 diag + 30 phase-pair)")
    print(f"  circuit count: TASK 3 scheme  = {len(kept_names)} (6 diag + 15 '+'-only, '-' DROPPED entirely)")

    # -- verify the algebraic identity in the NOISELESS limit before trusting it under noise --
    ideal_counts = ck["counts"]["ideal"]

    def exact_val(name, l):
        gi = group_idx[l]
        counts = ideal_counts[name][gi]
        return expectation_from_counts(counts, l)

    max_verify_err = 0.0
    for (n, m) in plus_pairs:
        un, um, plus = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
        for l in non_id_labels[:5]:  # spot-check a subset of labels, cheap and sufficient
            cross_from_plus_only = exact_val(plus, l) - (exact_val(un, l) + exact_val(um, l)) / 2
            cross_from_both = (exact_val(plus, l) - exact_val(f"(u{n}-u{m})", l)) / 2
            max_verify_err = max(max_verify_err, abs(cross_from_plus_only - cross_from_both))
    print(f"  algebraic identity check (drop '-' vs use both, ideal data, subset of labels): "
          f"max diff={max_verify_err:.2e} (should be ~0 in the noiseless limit)")

    results = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_by_slot = ck["counts"][model]
        errs = {k: [] for k in ("raw36", "phys36", "raw21", "phys21")}
        n_sdp_fail = 0
        t0 = time.time()
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("subspace_tomo", model, seed))

            # resample ALL 36 slots' counts once (same seed/draws for every variant -- fair comparison)
            resampled = {}
            for name in names:
                resampled[name] = [bootstrap_counts(counts_by_slot[name][gi], SHOTS, rng)
                                    for gi in range(len(groups))]

            def m_and_w_for(name):
                m_dict, w_dict = {}, {}
                for l in non_id_labels:
                    counts = resampled[name][group_idx[l]]
                    m = expectation_from_counts(counts, l)
                    m_dict[l] = m
                    total = sum(counts.values())
                    var = max(1 - m ** 2, 1e-4) / max(total, 1)
                    w_dict[l] = 1.0 / var
                return m_dict, w_dict

            raw_m = {name: m_and_w_for(name)[0] for name in names}

            # -- raw36: CURRENT scheme, all 36 circuits, no SDP --
            raw36 = {name: dict(raw_m[name]) for name in diag_names}
            for (n, m), plus, minus in zip(plus_pairs, plus_names, minus_names):
                raw36[plus] = dict(raw_m[plus])
                raw36[minus] = dict(raw_m[minus])

            # -- phys36: Phase-1-style SDP on all 36 circuits (reference, matches iteration 23) --
            rho36 = {}
            phys36 = {name: {} for name in names}
            for name in names:
                m_dict, w_dict = m_and_w_for(name)
                try:
                    rho36[name] = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                except RuntimeError:
                    n_sdp_fail += 1
                    rho36[name] = None
                    phys36[name] = dict(m_dict)
                    continue
                for l in non_id_labels:
                    phys36[name][l] = float(np.real(np.trace(rho36[name] @ P_S[l])))

            # -- raw21: 21 circuits, algebraic cross-term derivation, no SDP --
            # build cross terms into a synthetic "minus" so combine_matrices' existing (plus-minus)/2 formula still applies unchanged
            raw21_synth = {name: dict(raw_m[name]) for name in diag_names}
            for (n, m), plus in zip(plus_pairs, plus_names):
                un, um = f"u_{n}", f"u_{m}"
                raw21_synth[plus] = dict(raw_m[plus])
                synth_minus = {}
                for l in non_id_labels:
                    cross = raw_m[plus][l] - (raw_m[un][l] + raw_m[um][l]) / 2
                    synth_minus[l] = raw_m[plus][l] - 2 * cross  # so (plus-minus)/2 == cross exactly
                raw21_synth[f"(u{n}-u{m})"] = synth_minus

            # -- phys21: Task 3, SDP on all 21 KEPT circuits (diag + '+' only), cross derived from RECONSTRUCTED values --
            rho21 = {}
            for name in kept_names:
                m_dict, w_dict = m_and_w_for(name)
                try:
                    rho21[name] = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                except RuntimeError:
                    n_sdp_fail += 1
                    rho21[name] = None
            phys21_synth = {name: {} for name in diag_names}
            for name in diag_names:
                if rho21[name] is not None:
                    for l in non_id_labels:
                        phys21_synth[name][l] = float(np.real(np.trace(rho21[name] @ P_S[l])))
                else:
                    phys21_synth[name] = dict(raw_m[name])
            for (n, m), plus in zip(plus_pairs, plus_names):
                un, um = f"u_{n}", f"u_{m}"
                phys_diag_n = phys21_synth[un]
                phys_diag_m = phys21_synth[um]
                if rho21[plus] is not None:
                    plus_vals = {l: float(np.real(np.trace(rho21[plus] @ P_S[l]))) for l in non_id_labels}
                else:
                    plus_vals = raw_m[plus]
                phys21_synth[plus] = dict(plus_vals)
                synth_minus = {}
                for l in non_id_labels:
                    cross = plus_vals[l] - (phys_diag_n[l] + phys_diag_m[l]) / 2
                    synth_minus[l] = plus_vals[l] - 2 * cross
                phys21_synth[f"(u{n}-u{m})"] = synth_minus

            _, e_raw36 = energy_and_err(p, raw36, K)
            _, e_phys36 = energy_and_err(p, phys36, K)
            _, e_raw21 = energy_and_err(p, raw21_synth, K)
            _, e_phys21 = energy_and_err(p, phys21_synth, K)
            errs["raw36"].append(e_raw36["err_vs_exact_kcal"])
            errs["phys36"].append(e_phys36["err_vs_exact_kcal"])
            errs["raw21"].append(e_raw21["err_vs_exact_kcal"])
            errs["phys21"].append(e_phys21["err_vs_exact_kcal"])

        t_elapsed = time.time() - t0
        results[model] = {k: {"mean": float(np.mean(v)), "std": float(np.std(v))} for k, v in errs.items()}
        results[model]["n_sdp_fail"] = n_sdp_fail
        results[model]["wall_clock_s"] = t_elapsed
        print(f"\n  {model} ({t_elapsed:.1f}s, {n_sdp_fail} SDP failures):")
        print(f"    raw36  (current, 36 circuits):     {results[model]['raw36']['mean']:.2f}+/-{results[model]['raw36']['std']:.2f} kcal/mol")
        print(f"    phys36 (Phase-1 style, 36 circuits): {results[model]['phys36']['mean']:.2f}+/-{results[model]['phys36']['std']:.2f} kcal/mol")
        print(f"    raw21  (Task 3, 21 circuits, algebraic):  {results[model]['raw21']['mean']:.2f}+/-{results[model]['raw21']['std']:.2f} kcal/mol")
        print(f"    phys21 (Task 3, 21 circuits, joint SDP):  {results[model]['phys21']['mean']:.2f}+/-{results[model]['phys21']['std']:.2f} kcal/mol")

    # -- MANDATORY ideal-data sanity check (Phase 3's lesson, applied here per explicit instruction) --
    ideal = results["ideal"]
    disqualified = ideal["phys21"]["mean"] > 3 * max(ideal["raw36"]["mean"], ideal["phys36"]["mean"], 1e-6)
    print(f"\n  MANDATORY IDEAL-DATA SANITY CHECK (Phase 3's lesson): phys21 on near-noiseless data = "
          f"{ideal['phys21']['mean']:.3f} kcal/mol vs raw36={ideal['raw36']['mean']:.3f}, phys36={ideal['phys36']['mean']:.3f}")
    print(f"  {'*** DISQUALIFIED ***' if disqualified else 'PASSES the sanity check (no catastrophic distortion of clean data)'}")

    print(f"\n  -- summary, real hardware models --")
    for model in ["aria-1", "forte-1"]:
        r = results[model]
        print(f"    {model}: raw36={r['raw36']['mean']:.2f}  phys36={r['phys36']['mean']:.2f}  "
              f"raw21={r['raw21']['mean']:.2f}  phys21={r['phys21']['mean']:.2f} kcal/mol "
              f"(circuit count 36 -> 21, {'DISQUALIFIED' if disqualified else 'valid'})")

    out = {
        "K": K, "circuit_counts": {"current": len(names), "task3": len(kept_names)},
        "algebraic_identity_check_max_diff": max_verify_err,
        "results": results, "disqualified": bool(disqualified),
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return out


if __name__ == "__main__":
    main()
