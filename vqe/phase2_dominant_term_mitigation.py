#!/usr/bin/env python3
"""
phase2_dominant_term_mitigation.py — Phase 2 of physics-constrained
reconstruction (continued past Phase 1's decision-rule ABANDON at the
user's explicit direction — Phase 1 found a REAL, non-noise 1.20x/1.31x
error reduction, short of the pre-committed 2x/3x bar; the user chose to
continue rather than stop there. Recorded honestly: this is a deliberate
user override of the written decision rule, not a case of quietly moving
the goalposts).
============================================================================
IDEA: the H4 fragment's full 8-qubit Hamiltonian has 185 Pauli terms
(confirmed: iteration 20's CS-VQE work already established this count
independently). Each contributes coeff * (bilinear Schmidt-matrix
combination) to the total electronic energy -- and because the ENERGY IS
A LINEAR SUM over these 185 per-term contributions (even though each
individual contribution is itself bilinear in the K x K alpha/beta
matrix elements), mitigating only the alpha-register Pauli LABELS that
feed the highest-|contribution| terms, and leaving the long tail on
cheap raw estimates, should recover most of the benefit of full
mitigation at a fraction of the cost -- IF the energy truly is dominated
by a small head.

EXACT (classical, no estimation) per-term ranking: reuses Phase 1's own
P_S = U^dagger P U projection (already verified exact and Hermitian) as
the EXACT alpha matrix, and this project's own already-verified
beta-from-alpha shortcut (`derive_beta_matrices`: beta = S . alpha . S,
S = diag(signs)) to get the exact beta matrix -- no new derivation,
reusing verified machinery. Per this project's honesty rules: ranking
terms by their EXACT (not noisy-estimated) contribution uses information
that would not be available on a genuinely blind, unknown-answer
deployment -- disclosed explicitly here, not hidden. This is a validation
study on a KNOWN molecule; a real blind deployment would need a first
noisy pass to estimate which terms are dominant, which could differ from
the true ranking used here.

"Expensive mitigation" = Phase 1's physics-constrained SDP reconstruction
(applied only to slots' entries for HEAD alpha-labels). "Symmetry-
constrained raw" for the tail = the plain raw measured expectation value
(already confined to the physical weight-2 sector by the STATE PREP
circuit itself, per fixed_ansatz.py's own verified leakage~1e-29 finding
-- "symmetry-constrained" in that sense, even without going through the
Phase 1 SDP).

Sweeps the head/tail cutoff (90%, 95%, 99%, 99.9% of total |contribution|
mass) rather than reporting one cherry-picked threshold -- per this
project's "floor-test every free parameter" rule, the cutoff choice is
shown as a full tradeoff curve, not a single favorable point.

Run:
    python vqe/phase2_dominant_term_mitigation.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, combine_matrices, energy_from_alpha_matrices, derive_beta_matrices,
    HARTREE_TO_KCAL_MOL, slot_names,
)
from ionq_simulator_binding_curve import bootstrap_counts, stable_seed, expectation_from_counts
from phys_constrained_reconstruction import (
    load_checkpoint, build_group_index, build_P_S, reconstruct_rho_slot, K, SHOTS, N_SEEDS,
)

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "phase2_dominant_term_results.json")
CUTOFFS = [0.90, 0.95, 0.99, 0.999]


def rank_terms(p, P_S):
    """Exact per-term contribution to E_elec (Hartree, unnormalized by
    norm2 division is applied for absolute-energy interpretation), using
    the EXACT classical P_S as the alpha matrix and derive_beta_matrices
    for beta -- the SAME formula ef_energy_from_noisy_matrices sums over,
    evaluated per-term instead of summed."""
    terms = p["terms"]
    lambdas = p["lambdas"]
    signs = p["signs"]
    exact_alpha_mats = {l: P_S[l] for l in P_S}
    exact_beta_mats = derive_beta_matrices(exact_alpha_mats, signs, K)
    norm2 = sum(lambdas[n] ** 2 for n in range(K))

    contributions = []
    total_check = 0.0
    for alpha_label, beta_label, coeff in terms:
        Amat = exact_alpha_mats[alpha_label]
        Bmat = exact_beta_mats[beta_label]
        diag = sum(lambdas[n] ** 2 * Amat[n, n] * Bmat[n, n] for n in range(K))
        cross = sum(2 * lambdas[n] * lambdas[m] * np.real(Amat[n, m] * Bmat[n, m])
                    for n in range(K) for m in range(K) if n < m)
        contrib = coeff.real * (diag.real + cross) / norm2
        contributions.append({"alpha_label": alpha_label, "beta_label": beta_label,
                               "coeff": coeff.real, "contribution_ha": float(contrib)})
        total_check += contrib
    reconstructed_exact_energy = total_check + p["enuc"]
    mismatch = abs(reconstructed_exact_energy - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    return contributions, mismatch


def head_tail_split(contributions, cutoff_fraction):
    ranked = sorted(contributions, key=lambda c: abs(c["contribution_ha"]), reverse=True)
    total_abs = sum(abs(c["contribution_ha"]) for c in ranked)
    cum = 0.0
    head_terms = []
    for c in ranked:
        if cum >= cutoff_fraction * total_abs:
            break
        head_terms.append(c)
        cum += abs(c["contribution_ha"])
    head_alpha_labels = set(c["alpha_label"] for c in head_terms)
    return head_terms, head_alpha_labels


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def main():
    print("\n" + "=" * 96)
    print("  phase2_dominant_term_mitigation.py -- continued past Phase 1's ABANDON at user's explicit direction")
    print("=" * 96)

    ck = load_checkpoint()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=ck["d"], K=K)
    print(f"  {len(p['terms'])} Pauli terms in the full Hamiltonian (cross-check vs iteration 20's 185: "
          f"{'MATCH' if len(p['terms']) == 185 else 'MISMATCH -- ' + str(len(p['terms']))})")

    U = np.asarray(p["u_vecs"]).T
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    non_id_labels = [l for l in alpha_labels if l != identity_label]
    groups = ck["groups"]
    group_idx = build_group_index(groups)
    P_S = build_P_S(alpha_labels, U)

    contributions, mismatch = rank_terms(p, P_S)
    print(f"  exact per-term reconstruction cross-check vs setup_fragment's own exact_energy: "
          f"{mismatch:.2e} kcal/mol (should be ~0)")
    assert mismatch < 1e-6, "per-term exact energy decomposition does not sum correctly -- stop"

    names = ck["target_names"]
    exact_targets = _target_states_exact(p, K)

    # PASS 1: compute (once per model, seed, slot) the raw measured dict, the SDP-reconstructed
    # physics-constrained dict, and the exact ground-truth dict -- these do NOT depend on the
    # head/tail cutoff, so caching them here avoids re-solving 288 SDPs per cutoff (4x waste).
    print(f"\n  -- pass 1: reconstructing raw / physics-constrained / exact values (cutoff-independent) --")
    cache = {}  # (model, seed) -> {"m": {name: {label: val}}, "phys": {...}, "exact": {...}}
    for model in ["aria-1", "forte-1"]:
        counts_by_slot = ck["counts"][model]
        t0 = time.time()
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("phase2", model, seed))
            m_all = {name: {} for name in names}
            phys_all = {name: {} for name in names}
            exact_all = {name: {} for name in names}
            for name in names:
                resampled = [bootstrap_counts(counts_by_slot[name][gi], SHOTS, rng)
                             for gi in range(len(groups))]
                m_dict, w_dict = {}, {}
                for l in non_id_labels:
                    counts = resampled[group_idx[l]]
                    m = expectation_from_counts(counts, l)
                    m_dict[l] = m
                    total = sum(counts.values())
                    w_dict[l] = 1.0 / max(1 - m ** 2, 1e-4) / max(total, 1)
                rho_slot = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
                vec = exact_targets[name]
                for l in non_id_labels:
                    m_all[name][l] = m_dict[l]
                    phys_all[name][l] = float(np.real(np.trace(rho_slot @ P_S[l])))
                    exact_all[name][l] = float(np.real(np.conj(vec) @ P_S[l] @ vec))
            cache[(model, seed)] = {"m": m_all, "phys": phys_all, "exact": exact_all}
        print(f"    {model}: {time.time()-t0:.1f}s")

    print(f"\n  -- pass 2: head/tail cutoff sweep (fraction of total |contribution| mass in the head) --")
    sweep_results = {}
    for cutoff in CUTOFFS:
        head_terms, head_alpha_labels = head_tail_split(contributions, cutoff)
        print(f"\n  cutoff={cutoff*100:.1f}%: {len(head_terms)}/{len(contributions)} terms in head "
              f"({len(head_terms)/len(contributions)*100:.1f}%), "
              f"{len(head_alpha_labels)}/{len(non_id_labels)} unique alpha-labels need expensive mitigation")

        for model in ["aria-1", "forte-1"]:
            errs_raw, errs_full_phys, errs_hybrid = [], [], []
            errs_head_only, errs_tail_only = [], []
            for seed in range(N_SEEDS):
                c = cache[(model, seed)]
                raw_baseline = {name: dict(c["m"][name]) for name in names}
                raw_full_phys = {name: dict(c["phys"][name]) for name in names}
                raw_hybrid = {name: {} for name in names}
                raw_head_only = {name: {} for name in names}   # head=raw, tail=exact
                raw_tail_only = {name: {} for name in names}   # head=exact, tail=raw
                for name in names:
                    for l in non_id_labels:
                        if l in head_alpha_labels:
                            raw_hybrid[name][l] = c["phys"][name][l]
                            raw_head_only[name][l] = c["m"][name][l]
                            raw_tail_only[name][l] = c["exact"][name][l]
                        else:
                            raw_hybrid[name][l] = c["m"][name][l]
                            raw_head_only[name][l] = c["exact"][name][l]
                            raw_tail_only[name][l] = c["m"][name][l]

                _, err_raw = energy_and_err(p, raw_baseline, K)
                _, err_full = energy_and_err(p, raw_full_phys, K)
                _, err_hybrid = energy_and_err(p, raw_hybrid, K)
                _, err_head_only = energy_and_err(p, raw_head_only, K)
                _, err_tail_only = energy_and_err(p, raw_tail_only, K)
                errs_raw.append(err_raw["err_vs_exact_kcal"])
                errs_full_phys.append(err_full["err_vs_exact_kcal"])
                errs_hybrid.append(err_hybrid["err_vs_exact_kcal"])
                errs_head_only.append(err_head_only["err_vs_exact_kcal"])
                errs_tail_only.append(err_tail_only["err_vs_exact_kcal"])

            key = f"{cutoff}|{model}"
            sweep_results[key] = {
                "cutoff": cutoff, "model": model,
                "n_head_terms": len(head_terms), "n_total_terms": len(contributions),
                "n_head_labels": len(head_alpha_labels), "n_total_labels": len(non_id_labels),
                "raw_mean": float(np.mean(errs_raw)), "raw_std": float(np.std(errs_raw)),
                "full_phys_mean": float(np.mean(errs_full_phys)), "full_phys_std": float(np.std(errs_full_phys)),
                "hybrid_mean": float(np.mean(errs_hybrid)), "hybrid_std": float(np.std(errs_hybrid)),
                "head_error_isolated_mean": float(np.mean(errs_head_only)), "head_error_isolated_std": float(np.std(errs_head_only)),
                "tail_error_isolated_mean": float(np.mean(errs_tail_only)), "tail_error_isolated_std": float(np.std(errs_tail_only)),
            }
            r = sweep_results[key]
            print(f"    {model}: raw={r['raw_mean']:.2f}+/-{r['raw_std']:.2f}  "
                  f"full_phys={r['full_phys_mean']:.2f}+/-{r['full_phys_std']:.2f}  "
                  f"hybrid={r['hybrid_mean']:.2f}+/-{r['hybrid_std']:.2f} kcal/mol")
            print(f"      error-split (isolated): head-labels-raw-only={r['head_error_isolated_mean']:.2f}+/-{r['head_error_isolated_std']:.2f}  "
                  f"tail-labels-raw-only={r['tail_error_isolated_mean']:.2f}+/-{r['tail_error_isolated_std']:.2f} kcal/mol")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"sweep": sweep_results, "n_total_terms": len(contributions),
                    "n_total_labels": len(non_id_labels)}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return sweep_results


def _target_states_exact(p, K):
    """The exact K-DIM Schmidt-basis COEFFICIENT vector for each of the
    36 slots -- ground truth, no noise, no SDP. BUG CAUGHT: an earlier
    version returned the 16-dim PHYSICAL vectors (p['u_vecs'] entries
    live in the full 16-dim alpha register), which do not match P_S's
    already-projected K x K shape -- caught immediately by a matmul
    dimension-mismatch crash on the very first real run, not silently
    wrong. Fixed by using the TRIVIAL, exactly-known K-dim coefficient
    vector directly: since every slot is BY CONSTRUCTION a linear
    combination of Schmidt basis vectors, its representation in that
    same basis is just the standard unit vector (or their normalized
    sum/difference) -- no projection needed at all, simpler and correct
    by construction rather than by an unnecessary extra matmul."""
    targets = {}
    for n in range(K):
        e = np.zeros(K); e[n] = 1.0
        targets[f"u_{n}"] = e
    for n in range(K):
        for m in range(K):
            if n < m:
                e_n = np.zeros(K); e_n[n] = 1.0
                e_m = np.zeros(K); e_m[m] = 1.0
                targets[f"(u{n}+u{m})"] = (e_n + e_m) / np.sqrt(2)
                targets[f"(u{n}-u{m})"] = (e_n - e_m) / np.sqrt(2)
    return targets


if __name__ == "__main__":
    main()
