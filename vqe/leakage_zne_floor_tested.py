#!/usr/bin/env python3
"""
leakage_zne_floor_tested.py — does removing particle-number LEAKAGE
(iteration 18, Task D: a real, hardware-confirmed ~5-8% of shots landing
outside the physically valid weight-2 sector) before fitting ZNE unlock
the plateau that raw ZNE has never found on this circuit (iteration 13:
NO PLATEAU in any of 3 independent sweep directions; iteration 14: same,
even after Z2 tapering)?
============================================================================
THE MATHEMATICAL MOTIVATION (why this combination hasn't been tried
before, and why it's a genuinely different thing from either prior
result alone): gate-folding ZNE's entire premise is that the noisy
expectation value E(s) is a SMOOTH, well-behaved function of the noise
scale s (polynomial or exponential, extrapolatable to s=0). That premise
implicitly assumes the DOMINANT noise process is something like
depolarizing/dephasing error that scales the SIGNAL continuously toward
zero. Leakage is a fundamentally different kind of error: a shot either
falls outside the valid sector or it doesn't -- a discrete, not
continuous, corruption of the estimator. If leakage is a meaningful
fraction of the total noise budget (iteration 18 measured 7-8% real,
hardware-confirmed), it plausibly injects exactly the kind of
non-smooth, outlier-driven behavior that would make E(s) fail to follow
any simple polynomial/exponential form -- a concrete, checkable
candidate explanation for why iteration 13's exhaustive ZNE sweep never
found a plateau. This file tests it directly: apply iteration 18's
verified-exact ancilla-based leakage detector BEFORE fitting ZNE, at
EVERY noise scale, and see whether the resulting (now leakage-cleaned)
E(s) curve behaves better.

METHOD: reuses zne_floor_tested.py's local synthetic noise model
(P2_PER_GATE/P1_PER_GATE depolarizing, scaled by integer factors 1-7,
the SAME model and SAME scale ranges already used for the raw/tapered
ZNE sweeps, for a fair apples-to-apples comparison) and iteration 18's
verified-exact ancilla parity-check circuit (with_ancilla_parity,
imported unchanged from spin_leakage_postselect_ionq.py). At each scale,
computes the EXACT noisy probability distribution over all 32 (2^5)
basis outcomes for every (target, measurement group) via a density-
matrix simulation -- not just the Pauli expectation value directly,
because post-selection requires the JOINT distribution (ancilla outcome
correlated with the measured Pauli outcome in the SAME shot), which the
project's existing shot_sample() helper (independent per-label binomial
sampling) cannot represent. 8 independent seeds draw a real multinomial
sample from this joint distribution at each scale; RAW expectation
values marginalize the ancilla away (as pauli_expectation already does
for any label shorter than the full bitstring); POST-SELECTED values
condition on ancilla=0 first, renormalize, then compute -- from the SAME
per-seed sample, so the two schemes are properly correlated draws of the
same underlying experiment, not independent noise realizations.

Applies the SAME mandatory 3-direction floor test (qforge.floor_test)
this project has used since iteration 13, to BOTH the raw and the
post-selected E(s) curves, so the comparison is apples-to-apples and a
"plateau" claim is never eyeballed.

Run:
    python vqe/leakage_zne_floor_tested.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, floor_test,
    P2_PER_GATE, P1_PER_GATE,
)
from spin_leakage_postselect_ionq import with_ancilla_parity
import ef_fragment as effrag
from ionq_run import basis_change, pauli_expectation
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

K = 6
BASIS_GATES = ["u3", "cx"]
SHOTS = 100_000
N_SEEDS = 8
SCALE_RANGES = [
    [1, 2, 3],
    [1, 2, 3, 4],
    [1, 2, 3, 4, 5],
    [1, 2, 3, 4, 5, 6],
    [1, 2, 3, 4, 5, 6, 7],
]
MAX_SCALE = max(max(r) for r in SCALE_RANGES)
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "leakage_zne_floor_tested_results.json")


def build_noise_model(scale):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(min(P2_PER_GATE * scale, 0.75), 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(min(P1_PER_GATE * scale, 0.75), 1), "u3")
    return nm


def noisy_probs_5q(angles, group_combined_label, noise_model):
    """Full 32-outcome probability vector for state-prep + ancilla-parity
    CNOTs + this group's basis-change gates, under the given (scaled)
    noise model. index i's standard qiskit bitstring is format(i,'05b')
    (qubit0 = LSB = rightmost char), matching pauli_expectation's own
    documented convention exactly."""
    from fixed_ansatz import build_ansatz
    base = build_ansatz(angles)
    qc5 = with_ancilla_parity(base)
    qc5b = qc5.copy()
    basis_change(qc5b, group_combined_label)
    qc5b = transpile(qc5b, basis_gates=BASIS_GATES, optimization_level=0)
    qc5b.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc5b).result()
    dm = np.asarray(result.data(0)["density_matrix"])
    probs = np.clip(np.real(np.diag(dm)), 0.0, None)
    probs = probs / probs.sum()
    return probs


def expectation_raw_from_probs_dict(probs_dict, label):
    return pauli_expectation(probs_dict, label)


def expectation_post_from_probs_dict(probs_dict, label):
    kept = {bs: p for bs, p in probs_dict.items() if bs[0] == "0"}  # ancilla=qubit4=leftmost char
    total = sum(kept.values())
    if total <= 0:
        return 0.0
    kept = {bs: p / total for bs, p in kept.items()}
    return pauli_expectation(kept, label)


def sample_counts(probs32, n_shots, rng):
    draws = rng.multinomial(n_shots, probs32)
    return {format(i, "05b"): int(n) for i, n in enumerate(draws) if n > 0}


def counts_to_probs_dict(counts):
    total = sum(counts.values())
    return {bs: c / total for bs, c in counts.items()} if total else {}


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def exponential_fit(scales, energies):
    diffs = np.array(energies) - energies[-1]
    if np.all(diffs > 0) or np.all(diffs < 0):
        slope, intercept = np.polyfit(scales, np.log(np.abs(diffs)), 1)
        A = np.sign(diffs[0]) * np.exp(intercept)
        return float(energies[-1] + A), None
    return None, "energies not monotonic across scales -- exponential fit not well-posed, not faked"


def run_floor_test_suite(scheme_name, per_seed_energy_fn, p, scale_ranges):
    """per_seed_energy_fn(seed, scale) -> exact energy estimate (Hartree)
    for that scheme, scale, and seed. Mirrors zne_floor_tested.py's own
    fit/aggregate/floor-test structure exactly, applied to whichever
    scheme (raw or post-selected) is passed in."""
    range_order_results, range_exp_results = {}, {}
    for scale_range in scale_ranges:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            range_order_results[(range_key, order)] = []
        range_exp_results[range_key] = []

        for seed in range(N_SEEDS):
            energies = [per_seed_energy_fn(seed, s) for s in scale_range]
            for order in range(1, max_order + 1):
                coeffs = np.polyfit(scale_range, energies, order)
                E0 = float(np.polyval(coeffs, 0))
                err = abs(E0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
                range_order_results[(range_key, order)].append(err)
            E0_exp, _ = exponential_fit(scale_range, energies)
            if E0_exp is not None:
                range_exp_results[range_key].append(abs(E0_exp - p["exact_energy"]) * HARTREE_TO_KCAL_MOL)

    summary = {}
    for scale_range in scale_ranges:
        range_key = ",".join(str(s) for s in scale_range)
        max_order = len(scale_range) - 1
        for order in range(1, max_order + 1):
            errs = range_order_results[(range_key, order)]
            summary[f"{range_key}|order{order}"] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs))}
        if range_exp_results[range_key]:
            errs = range_exp_results[range_key]
            summary[f"{range_key}|exp"] = {"mean_kcal": float(np.mean(errs)), "std_kcal": float(np.std(errs))}

    range_sizes = [len(r) for r in scale_ranges]
    lin_errs = [summary[f"{','.join(str(s) for s in r)}|order1"]["mean_kcal"] for r in scale_ranges]
    ft_range_linear = floor_test(range_sizes, lin_errs)

    widest = scale_ranges[-1]
    widest_key = ",".join(str(s) for s in widest)
    orders = list(range(1, len(widest)))
    order_errs = [summary[f"{widest_key}|order{o}"]["mean_kcal"] for o in orders]
    ft_order = floor_test(orders, order_errs)

    quad_ranges = [r for r in scale_ranges if len(r) > 2]
    quad_errs = [summary[f"{','.join(str(s) for s in r)}|order2"]["mean_kcal"] for r in quad_ranges]
    quad_sizes = [len(r) for r in quad_ranges]
    ft_range_quad = floor_test(quad_sizes, quad_errs)

    any_plateau = ft_range_linear["has_floor"] or ft_order["has_floor"] or ft_range_quad["has_floor"]
    print(f"\n  [{scheme_name}] range@order1: {ft_range_linear['verdict']}")
    print(f"  [{scheme_name}] order@widest-range: {ft_order['verdict']}")
    print(f"  [{scheme_name}] range@order2: {ft_range_quad['verdict']}")
    print(f"  [{scheme_name}] ANY PLATEAU: {any_plateau}")

    return {
        "summary": summary,
        "floor_test_range_at_order1": ft_range_linear,
        "floor_test_order_at_widest_range": ft_order,
        "floor_test_range_at_order2": ft_range_quad,
        "any_plateau_found": bool(any_plateau),
    }


def main():
    print("\n" + "=" * 96)
    print("  leakage_zne_floor_tested.py -- does leakage post-selection unlock a ZNE plateau?")
    print("=" * 96)

    t0 = time.time()
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    assert n_ok == 36, f"fit did not converge for all 36: {n_ok}/36"
    print(f"  setup OK: 36/36 converged (worst={worst:.2e}), {time.time()-t0:.1f}s")

    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    groups = effrag.group_labels_qubit_wise(alpha_labels)
    target_names = sorted(solutions.keys())
    print(f"  {len(target_names)} targets x {len(groups)} groups, scales 1..{MAX_SCALE}")

    print(f"\n  computing exact noisy 32-outcome distributions at every (target, group, scale) "
          f"-- {len(target_names)*len(groups)*MAX_SCALE} density-matrix sims, cached once")
    t0 = time.time()
    probs_cache = {}  # (name, gi, scale) -> 32-dim prob array
    for s in range(1, MAX_SCALE + 1):
        nm = build_noise_model(s)
        for name in target_names:
            angles = solutions[name]["angles"]
            for gi, group in enumerate(groups):
                combined = effrag.combined_basis_label(group)
                probs_cache[(name, gi, s)] = noisy_probs_5q(angles, combined, nm)
        print(f"    scale={s} done, {time.time()-t0:.1f}s elapsed")

    # quick exact (infinite-shot) sanity check at each scale, raw vs post-selected
    print(f"\n  -- exact (no shot noise) sanity check per scale --")
    group_labels = [effrag.combined_basis_label(g) for g in groups]
    for s in range(1, MAX_SCALE + 1):
        raw_exact, post_exact = {n: {} for n in target_names}, {n: {} for n in target_names}
        for name in target_names:
            for gi, group in enumerate(groups):
                probs32 = probs_cache[(name, gi, s)]
                probs_dict = {format(i, "05b"): probs32[i] for i in range(32) if probs32[i] > 0}
                for l in group:
                    raw_exact[name][l] = expectation_raw_from_probs_dict(probs_dict, l)
                    post_exact[name][l] = expectation_post_from_probs_dict(probs_dict, l)
        _, err_raw = energy_and_err(p, raw_exact, K)
        _, err_post = energy_and_err(p, post_exact, K)
        print(f"    scale={s}: RAW={err_raw:.3f}  POST-SELECTED={err_post:.3f} kcal/mol (exact, no shot noise)")

    print(f"\n  -- 8-seed shot-noisy sweep (joint multinomial sampling per (target,group,scale)) --")

    def per_seed_energy(scheme, seed, scale):
        rng = np.random.default_rng(seed * 7919 + 13 + scale * 104729)
        raw = {n: {} for n in target_names}
        for name in target_names:
            for gi, group in enumerate(groups):
                probs32 = probs_cache[(name, gi, scale)]
                counts = sample_counts(probs32, SHOTS, rng)
                probs_dict = counts_to_probs_dict(counts)
                for l in group:
                    if scheme == "raw":
                        raw[name][l] = expectation_raw_from_probs_dict(probs_dict, l)
                    else:
                        raw[name][l] = expectation_post_from_probs_dict(probs_dict, l)
        E, _ = energy_and_err(p, raw, K)
        return E

    t0 = time.time()
    raw_results = run_floor_test_suite("RAW (no post-selection)",
                                        lambda seed, s: per_seed_energy("raw", seed, s), p, SCALE_RANGES)
    post_results = run_floor_test_suite("POST-SELECTED (leakage discarded)",
                                         lambda seed, s: per_seed_energy("post", seed, s), p, SCALE_RANGES)
    print(f"\n  shot-noisy sweep done, {time.time()-t0:.1f}s")

    print(f"\n  -- OVERALL VERDICT --")
    if post_results["any_plateau_found"] and not raw_results["any_plateau_found"]:
        print("  Leakage post-selection UNLOCKS a genuine plateau that raw ZNE never found. "
              "This is a real, actionable improvement -- worth testing for real on IonQ.")
    elif post_results["any_plateau_found"] and raw_results["any_plateau_found"]:
        print("  Both raw and post-selected ZNE find a plateau at the local synthetic-noise level -- "
              "leakage post-selection did not change the qualitative picture here, though the specific "
              "plateau values may still differ meaningfully.")
    elif not post_results["any_plateau_found"]:
        print("  NO PLATEAU FOUND for post-selected ZNE either, at any of the 3 sweep directions tested. "
              "The leakage hypothesis for ZNE's non-convergence does NOT hold up under this local test -- "
              "reported plainly, not forced.")

    results = {
        "K": K, "shots": SHOTS, "n_seeds": N_SEEDS, "scale_ranges": SCALE_RANGES,
        "p2_per_gate": P2_PER_GATE, "p1_per_gate": P1_PER_GATE,
        "raw": raw_results, "post_selected": post_results,
        "chemical_accuracy_kcal": 1.0,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
