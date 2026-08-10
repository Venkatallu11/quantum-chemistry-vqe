#!/usr/bin/env python3
"""
phase4_residual_zne.py — Phase 4 of physics-constrained reconstruction
(continuing past Phase 1's decision-rule ABANDON at the user's explicit
direction, as Phases 2 and 3 already did; Phase 3 was disqualified by
its own built-in sanity gate and contributes nothing valid here).
============================================================================
IDEA: apply ZNE only to whatever error REMAINS after the genuinely
valid corrections (Phase 1's physics-constrained SDP reconstruction,
Phase 2's dominant-term selective hybrid) -- not to the full raw
~33-42 kcal/mol baseline iterations 11/13/14/19 already exhaustively
showed never plateaus. Still requires a genuine PLATEAU under this
project's own mandatory 3-direction floor test (qforge.floor_test,
unchanged since iteration 13). If no plateau is found, report NO ZNE
number at all -- do not fit anyway, per explicit instruction.

WHY THIS NEEDS A LOCAL NOISE MODEL, NOT JUST EXISTING DATA: the real
IonQ checkpoints used in Phases 1-3 are a SINGLE noise level each (no
multiple scale factors) -- ZNE fundamentally needs multiple points along
a noise-scale axis. Reuses this project's own established local ZNE
machinery UNCHANGED (`zne_floor_tested.py`'s `build_noise_model`,
`measure_exact_noisy_raw`, `SCALE_RANGES`, `floor_test`) rather than
reinventing it, combined with Phase 1/2's own already-verified
reconstruction functions.

Run:
    python vqe/phase4_residual_zne.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL, combine_matrices, \
    energy_from_alpha_matrices, floor_test, shot_sample
from zne_floor_tested import build_noise_model, measure_exact_noisy_raw, SCALE_RANGES, MAX_SCALE
from phys_constrained_reconstruction import build_P_S, reconstruct_rho_slot, K
from phase2_dominant_term_mitigation import rank_terms, head_tail_split

SHOTS = 100_000
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "phase4_residual_zne_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs


def exponential_fit(scales, energies):
    diffs = np.array(energies) - energies[-1]
    if np.all(diffs > 0) or np.all(diffs < 0):
        slope, intercept = np.polyfit(scales, np.log(np.abs(diffs)), 1)
        A = np.sign(diffs[0]) * np.exp(intercept)
        return float(energies[-1] + A), None
    return None, "non-monotonic"


def apply_scheme(raw_exact_at_scale, P_S, non_id_labels, head_labels, scheme, K, names):
    """raw_exact_at_scale: {name: {label: exact_noisy_value}} for ONE scale
    (already possibly shot-perturbed by the caller). Returns the SAME
    shape dict, transformed according to `scheme`: 'raw' (pass through),
    'phase1' (per-slot physics-constrained SDP, all labels), 'phase2'
    (selective hybrid: physics-constrained for head labels, raw for
    tail)."""
    if scheme == "raw":
        return raw_exact_at_scale
    out = {name: {} for name in names}
    for name in names:
        m_dict = raw_exact_at_scale[name]
        w_dict = {l: 1.0 / max(1 - m_dict[l] ** 2, 1e-4) for l in non_id_labels}
        rho = reconstruct_rho_slot(P_S, m_dict, w_dict, K)
        for l in non_id_labels:
            if scheme == "phase1":
                out[name][l] = float(np.real(np.trace(rho @ P_S[l])))
            else:  # phase2: selective hybrid
                out[name][l] = float(np.real(np.trace(rho @ P_S[l]))) if l in head_labels else m_dict[l]
    return out


def run_floor_test_suite(scheme_name, exact_raw_by_scale, p, P_S, non_id_labels, head_labels, names):
    """SPEED FIX: every SCALE_RANGES entry is a strict prefix of the next
    ([1,2,3] subset [1,2,3,4] subset ... subset [1..7]), and each range's
    inner loop re-seeds its RNG identically per seed (`seed*7919+13`) and
    draws scales in order starting at 1 -- so for a fixed seed, the
    sequential shot_sample draws for scales 1..k are BYTE-IDENTICAL
    across every range that includes them as a prefix. The first version
    of this function recomputed apply_scheme (36 SDP solves, ~4-9s) for
    every (range, seed, scale) triple independently -- redundant by
    roughly 5x, since the same (seed, scale) pair appears in up to 5
    different ranges. Fixed by computing each (seed, scale) energy
    EXACTLY ONCE (iterating scale 1..MAX_SCALE in order per seed, the
    same draw sequence as before) and caching it, then building each
    range's energy list by lookup -- produces IDENTICAL numbers to the
    unoptimized version (verified: same seeding, same draw order), just
    without repeating the expensive part 5 times."""
    max_scale = max(max(r) for r in SCALE_RANGES)
    energy_cache = {}  # (seed, scale) -> E
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed * 7919 + 13)
        for s in range(1, max_scale + 1):
            shot_noisy = {name: {l: shot_sample(v, SHOTS, rng) for l, v in vals.items()}
                          for name, vals in exact_raw_by_scale[s].items()}
            corrected = apply_scheme(shot_noisy, P_S, non_id_labels, head_labels, scheme_name, K, names)
            E, _ = energy_and_err(p, corrected, K)
            energy_cache[(seed, s)] = E

    range_order_results, range_exp_results = {}, {}
    for scale_range in SCALE_RANGES:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            range_order_results[(range_key, order)] = []
        range_exp_results[range_key] = []
        for seed in range(N_SEEDS):
            energies = [energy_cache[(seed, s)] for s in scale_range]
            for order in range(1, max_order + 1):
                coeffs = np.polyfit(scale_range, energies, order)
                E0 = float(np.polyval(coeffs, 0))
                err = abs(E0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
                range_order_results[(range_key, order)].append(err)
            E0_exp, _ = exponential_fit(scale_range, energies)
            if E0_exp is not None:
                range_exp_results[range_key].append(abs(E0_exp - p["exact_energy"]) * HARTREE_TO_KCAL_MOL)

    summary = {}
    for scale_range in SCALE_RANGES:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            errs = range_order_results[(range_key, order)]
            summary[f"{range_key}|order{order}"] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs))}
        if range_exp_results[range_key]:
            errs = range_exp_results[range_key]
            summary[f"{range_key}|exp"] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs))}

    range_sizes = [len(r) for r in SCALE_RANGES]
    lin_errs = [summary[f"{','.join(str(s) for s in r)}|order1"]["mean_kcal"] for r in SCALE_RANGES]
    ft_range_linear = floor_test(range_sizes, lin_errs)
    widest = SCALE_RANGES[-1]
    widest_key = ",".join(str(s) for s in widest)
    orders = list(range(1, len(widest)))
    order_errs = [summary[f"{widest_key}|order{o}"]["mean_kcal"] for o in orders]
    ft_order = floor_test(orders, order_errs)
    quad_ranges = [r for r in SCALE_RANGES if len(r) > 2]
    quad_errs = [summary[f"{','.join(str(s) for s in r)}|order2"]["mean_kcal"] for r in quad_ranges]
    quad_sizes = [len(r) for r in quad_ranges]
    ft_range_quad = floor_test(quad_sizes, quad_errs)

    any_plateau = ft_range_linear["has_floor"] or ft_order["has_floor"] or ft_range_quad["has_floor"]
    print(f"\n  [{scheme_name}] range@order1: {ft_range_linear['verdict']}")
    print(f"  [{scheme_name}] order@widest-range: {ft_order['verdict']}")
    print(f"  [{scheme_name}] range@order2: {ft_range_quad['verdict']}")
    print(f"  [{scheme_name}] ANY PLATEAU: {any_plateau}")
    return {"summary": summary, "any_plateau_found": bool(any_plateau),
            "floor_test_range_at_order1": ft_range_linear, "floor_test_order_at_widest_range": ft_order,
            "floor_test_range_at_order2": ft_range_quad}


def main():
    print("\n" + "=" * 96)
    print("  phase4_residual_zne.py -- ZNE on the RESIDUAL error after Phase 1/2, continuing at user's direction")
    print("=" * 96)

    t0 = time.time()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    assert n_ok == 36
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    names = sorted(solutions.keys())
    print(f"  setup OK: 36/36 converged, {time.time()-t0:.1f}s")

    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    contributions, mismatch = rank_terms(p, P_S)
    assert mismatch < 1e-6
    _, head_labels = head_tail_split(contributions, 0.90)
    print(f"  reusing Phase 2's 90% cutoff: {len(head_labels)}/{len(non_id_labels)} head labels get "
          f"physics-constrained treatment at every scale")

    print(f"\n  measuring exact noisy matrices at scales 1..{MAX_SCALE} (shot-count independent, cached once)")
    exact_raw_by_scale = {}
    for s in range(1, MAX_SCALE + 1):
        nm = build_noise_model(s)
        exact_raw_by_scale[s] = measure_exact_noisy_raw(p, non_id_labels, nm)
        E_raw, err_raw = energy_and_err(p, exact_raw_by_scale[s], K)
        corrected_p1 = apply_scheme(exact_raw_by_scale[s], P_S, non_id_labels, head_labels, "phase1", K, names)
        E_p1, err_p1 = energy_and_err(p, corrected_p1, K)
        corrected_p2 = apply_scheme(exact_raw_by_scale[s], P_S, non_id_labels, head_labels, "phase2", K, names)
        E_p2, err_p2 = energy_and_err(p, corrected_p2, K)
        print(f"    scale={s}: RAW={err_raw['err_vs_exact_kcal']:.3f}  "
              f"Phase1={err_p1['err_vs_exact_kcal']:.3f}  Phase2-hybrid={err_p2['err_vs_exact_kcal']:.3f} kcal/mol "
              f"(exact, no shot noise)")

    print(f"\n  -- 8-seed shot-noisy sweep across every (scale range, fit order) combination --")
    results = {}
    for scheme in ["raw", "phase1", "phase2"]:
        t0 = time.time()
        results[scheme] = run_floor_test_suite(scheme, exact_raw_by_scale, p, P_S, non_id_labels, head_labels, names)
        print(f"  [{scheme}] sweep done, {time.time()-t0:.1f}s")

    print(f"\n  -- OVERALL VERDICT --")
    for scheme in ["raw", "phase1", "phase2"]:
        found = results[scheme]["any_plateau_found"]
        print(f"  {scheme}: {'PLATEAU FOUND -- see which direction above before trusting any specific number' if found else 'NO PLATEAU -- reporting no ZNE number for this scheme, per explicit instruction'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"scale_ranges": SCALE_RANGES, "n_seeds": N_SEEDS, "shots": SHOTS,
                    "head_labels": sorted(head_labels), "n_head_labels": len(head_labels),
                    "results": results, "chemical_accuracy_kcal": 1.0}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
