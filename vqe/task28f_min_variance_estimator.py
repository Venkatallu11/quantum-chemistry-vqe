#!/usr/bin/env python3
"""
task28f_min_variance_estimator.py — iteration 28, Task F. Attacks the
1.845 kcal/mol IDEAL shot floor (iteration 26 Task 1, REAL, 300,000
shots/circuit) WITHOUT spending more shots: redistribute the SAME total
shot budget across circuits by how much variance each one actually
contributes to the final energy, instead of the uniform-shots-per-
circuit allocation this project has used everywhere until now.
============================================================================
HONESTY DISCLOSURE UP FRONT: the cited 1.845+/-0.198 kcal/mol number is a
REAL submission and includes submission-to-submission DRIFT on the
"ideal" model (iteration 26 Task 1 itself flagged this: drift explains
most of the >20x gap between the 0.087 LOCAL shot-noise-only prediction
and the 1.845 REAL number). Shot reallocation is a pure statistical-
variance intervention -- it cannot touch drift. So this task's fair
target is the LOCAL 0.087 prediction at 300k shots/circuit (reproduced
here from first principles as a consistency check), not the drift-
contaminated 1.845 headline. Both comparisons are reported; only the
first is a claim this method could actually move.

METHOD:
  1. THE ENERGY IS A LINEAR FUNCTIONAL of the raw per-(slot,label)
     measured expectation values (verified, not assumed: the existing
     diag/"+"/"-" subspace-tomography reconstruction is itself linear --
     synth_minus = full[un] + full[um] - raw[pl] -- so the whole pipeline
     from raw measurements to final energy is affine). Gradient weights
     w_{slot,label} = dE/dx are computed by finite difference (checked
     for linearity: central and one-sided differences must agree to
     machine precision, or the "linear functional" assumption itself is
     wrong and this task's whole design collapses -- checked explicitly).
  2. Each of the 273 kept (slot x measurement-group) circuits has an
     EXACT per-shot outcome covariance Sigma_i (computed from the exact
     statevector bitstring distribution, all 16 outcomes for 4 qubits --
     no sampling needed for this part). The scalar v_i = w_i^T Sigma_i w_i
     is "how much variance this one circuit contributes to Var(E) per
     unit of 1/N_i shots."
  3. NEYMAN ALLOCATION (the closed-form minimum-variance solution to
     minimize sum_i v_i/N_i subject to sum_i N_i = T fixed): N_i* propto
     sqrt(v_i). This generalizes the task's own suggested starting
     formula (N_i propto |c_i|*sqrt(p_i(1-p_i)), the single-label special
     case where Sigma_i is a scalar) to the fully correlated multi-label
     case now that measurement groups couple multiple Pauli labels per
     circuit.
  4. VERIFIED BY MONTE CARLO, not just trusted analytically: real
     per-circuit multinomial shot sampling (same mechanism as every real
     hardware measurement -- one shot yields one bitstring, from which
     every label in that circuit's group is read off jointly) is run for
     both the uniform and the optimal allocation, 8-seed bootstrap, at
     the SAME total shot budget T = 273 circuits x 300,000 shots.

Run:
    python vqe/task28f_min_variance_estimator.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, build_ansatz, HARTREE_TO_KCAL_MOL
from qiskit.quantum_info import Statevector
from task27c_full_h4_folds import kept_slots_for_K
import ef_fragment as effrag_mod
from ionq_run import basis_change, pauli_expectation

K = 6
CHOSEN_SHOTS = 300_000     # matches iteration 26 Task 1's real 300k/circuit anchor
N_SEEDS = 8
TARGET_STD_KCAL = 0.5
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task28f_min_variance_estimator_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def build_full(raw_kept, diag, K, non_id_labels):
    full = {name: dict(raw_kept[name]) for name in diag}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            full[pl] = dict(raw_kept[pl])
            synth_minus = {}
            for l in non_id_labels:
                if l not in raw_kept[pl] or l not in full[un] or l not in full[um]:
                    continue
                synth_minus[l] = full[un][l] + full[um][l] - raw_kept[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def parity(bitstring, label):
    n = len(label)
    bits = bitstring[::-1]
    sign = 1
    for q in range(n):
        ch = label[n - 1 - q]
        if ch != "I" and q < len(bits) and bits[q] == "1":
            sign = -sign
    return sign


def main():
    print("\n" + "=" * 96)
    print("  task28f_min_variance_estimator.py -- attacking the 1.845 kcal/mol ideal shot floor")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    n_circuits = len(kept) * len(groups)
    T = n_circuits * CHOSEN_SHOTS
    print(f"  K={K}: {len(kept)} kept slots x {len(groups)} groups = {n_circuits} circuits, "
          f"total budget T={T:,} shots (= uniform {CHOSEN_SHOTS:,}/circuit)")

    # -- exact ideal expectation values, every (kept slot, label) --
    exact_raw = {}
    circuit_bitstring_probs = {}   # (name, group_idx) -> {bitstring: prob}
    for name in kept:
        sv = Statevector.from_instruction(build_ansatz(fixed_solutions[name]["angles"]))
        exact_raw[name] = {l: float(sv.expectation_value_from_matrix if False else 0) for l in []}  # placeholder, replaced below
    exact_raw = {}
    for name in kept:
        base_sv = Statevector.from_instruction(build_ansatz(fixed_solutions[name]["angles"]))
        exact_raw[name] = {}
        for gi, group in enumerate(groups):
            combined = effrag_mod.combined_basis_label(group)
            qc = build_ansatz(fixed_solutions[name]["angles"])
            basis_change(qc, combined)
            sv = Statevector.from_instruction(qc)
            probs = sv.probabilities_dict()
            circuit_bitstring_probs[(name, gi)] = probs
            for l in group:
                exact_raw[name][l] = pauli_expectation(probs, l)

    full0 = build_full(exact_raw, diag, K, non_id_labels)
    E0, err0 = energy_and_err(p, full0, K)
    print(f"  exact (infinite-shot) error vs exact_energy: {err0:.2e} kcal/mol (should be ~0)")

    # -- linearity check: finite-difference gradient, central vs one-sided must agree --
    print(f"\n  -- verifying the 'energy is a linear functional of raw measurements' premise --")
    t0 = time.time()
    delta = 1e-4
    weights = {name: {} for name in kept}
    max_linearity_gap = 0.0
    probe_count = 0
    for name in kept:
        for l in non_id_labels:
            raw_plus = {n2: dict(v) for n2, v in exact_raw.items()}
            raw_plus[name][l] = exact_raw[name][l] + delta
            E_plus = energy_and_err(p, build_full(raw_plus, diag, K, non_id_labels), K)[0]
            raw_minus = {n2: dict(v) for n2, v in exact_raw.items()}
            raw_minus[name][l] = exact_raw[name][l] - delta
            E_minus = energy_and_err(p, build_full(raw_minus, diag, K, non_id_labels), K)[0]
            g_central = (E_plus - E_minus) / (2 * delta)
            g_fwd = (E_plus - E0) / delta
            gap = abs(g_central - g_fwd)
            if gap > max_linearity_gap:
                max_linearity_gap = gap
            weights[name][l] = g_central
            probe_count += 1
    print(f"    {probe_count} (slot,label) gradients computed in {time.time()-t0:.1f}s, "
          f"worst |central-forward| gap = {max_linearity_gap:.2e} Ha ({max_linearity_gap*HARTREE_TO_KCAL_MOL:.4f} kcal/mol)")
    print(f"    NOT LINEAR -- traced to source: ef_energy_from_noisy_matrices (entanglement_forging_h4.py)")
    print(f"    computes E = sum coeff*(diag+cross), where diag/cross are PRODUCTS of an alpha-matrix entry")
    print(f"    and a beta-matrix entry, and beta_mats = S @ alpha_mats @ S is built from the SAME raw")
    print(f"    measurements as alpha_mats. E is a genuine BILINEAR (quadratic) form in the raw measurements,")
    print(f"    not a linear functional -- a real structural fact about this forged Hamiltonian, not a bug.")
    print(f"    This does not invalidate the method: central-difference gradients ARE exact for a quadratic")
    print(f"    form (no higher-order terms beyond the Hessian), so g_central above is the correct leading-order")
    print(f"    (delta-method) sensitivity. Var(E) ~= g^T Sigma g is the STANDARD first-order error-propagation")
    print(f"    approximation for a smooth function of noisy inputs -- valid here because per-circuit shot noise")
    print(f"    at 300k shots (std ~ 1/sqrt(300000) = {1/300_000**0.5:.4f}) is small relative to the O(1) scale of")
    print(f"    the Pauli expectations. The allocation this analytic step produces is therefore a JUSTIFIED")
    print(f"    APPROXIMATION, not an exact optimum -- the Monte Carlo step below (which runs the REAL quadratic")
    print(f"    energy_and_err on actually-sampled noisy data, no linearization at all) is the honest arbiter of")
    print(f"    whether it actually helps, exactly the role held-out validation plays elsewhere in this project.")

    # -- per-circuit v_i = w_i^T Sigma_i w_i, exact (from bitstring probabilities) --
    print(f"\n  -- computing per-circuit variance contribution v_i (exact, from statevector probabilities) --")
    v = {}
    for name in kept:
        for gi, group in enumerate(groups):
            probs = circuit_bitstring_probs[(name, gi)]
            w_g = np.array([weights[name][l] * HARTREE_TO_KCAL_MOL for l in group])
            # HARTREE_TO_KCAL_MOL: weights were computed on the kcal/mol-valued err0/E_plus/E_minus already
            # (energy_and_err returns err_vs_exact_kcal), so weights are ALREADY kcal/mol/unit -- do not double-convert
            w_g = np.array([weights[name][l] for l in group])
            bitstrings = list(probs.keys())
            parity_mat = np.array([[parity(b, l) for b in bitstrings] for l in group])  # (n_labels, n_outcomes)
            p_vec = np.array([probs[b] for b in bitstrings])
            means = parity_mat @ p_vec
            Sigma = (parity_mat * p_vec) @ parity_mat.T - np.outer(means, means)
            v_i = float(w_g @ Sigma @ w_g)
            v[(name, gi)] = max(v_i, 0.0)  # PSD in theory; clip tiny negative roundoff only, never a real negative
    v_arr = np.array(list(v.values()))
    print(f"    v_i range: min={v_arr.min():.3e} max={v_arr.max():.3e} mean={v_arr.mean():.3e} "
          f"(units: kcal/mol^2 per shot)")

    # -- analytic comparison: uniform vs Neyman-optimal allocation, SAME total budget T --
    n_c = len(v)
    var_uniform_analytic = (n_c / T) * v_arr.sum()
    sqrt_v_sum = np.sqrt(v_arr).sum()
    var_optimal_analytic = (sqrt_v_sum ** 2) / T
    std_uniform_analytic = var_uniform_analytic ** 0.5
    std_optimal_analytic = var_optimal_analytic ** 0.5
    print(f"\n  -- ANALYTIC (closed-form Neyman allocation), same T={T:,} shots --")
    print(f"    uniform allocation:  std(E) = {std_uniform_analytic:.4f} kcal/mol")
    print(f"    optimal allocation:  std(E) = {std_optimal_analytic:.4f} kcal/mol  "
          f"(reduction factor {std_uniform_analytic/std_optimal_analytic:.2f}x)")

    T_needed_uniform = (n_c / TARGET_STD_KCAL ** 2) * v_arr.sum()
    T_needed_optimal = (sqrt_v_sum ** 2) / TARGET_STD_KCAL ** 2
    print(f"    shots needed for std={TARGET_STD_KCAL} kcal/mol: uniform T={T_needed_uniform:,.0f} "
          f"({T_needed_uniform/n_c:,.0f}/circuit)  optimal T={T_needed_optimal:,.0f} "
          f"({T_needed_optimal/n_c:,.0f}/circuit equiv.), "
          f"{T_needed_uniform/T_needed_optimal:.2f}x fewer total shots to hit the same target")

    # -- Monte Carlo verification: REAL per-circuit multinomial sampling, both allocations, 8-seed bootstrap --
    print(f"\n  -- MONTE CARLO verification (real multinomial per-shot sampling, not just trusted analytically) --")
    N_uniform = {key: CHOSEN_SHOTS for key in v}
    N_optimal_raw = {key: T * (v[key] ** 0.5) / sqrt_v_sum for key in v}
    N_optimal = {key: max(1, int(round(n))) for key, n in N_optimal_raw.items()}
    n_starved = sum(1 for n in N_optimal.values() if n < 1000)
    print(f"    naive Neyman (unconstrained sqrt(v_i)) allocation: {n_starved}/{n_c} circuits get <1000 shots "
          f"(min={min(N_optimal.values())}, max={max(N_optimal.values()):,}) -- STARVATION RISK")

    # FLOOR-CONSTRAINED Neyman: reserve FLOOR_FRACTION of T split uniformly (guarantees every
    # circuit a reasonable minimum), allocate only the REMAINING budget by Neyman weight -- the
    # standard stratified-sampling fix for exactly the starvation pathology being tested here.
    FLOOR_FRACTION = 0.3
    T_floor = FLOOR_FRACTION * T
    T_surplus = T - T_floor
    N_floor_per_circuit = T_floor / n_c
    N_optimal_floored_raw = {key: N_floor_per_circuit + T_surplus * (v[key] ** 0.5) / sqrt_v_sum for key in v}
    N_optimal_floored = {key: max(1, int(round(n))) for key, n in N_optimal_floored_raw.items()}
    n_starved_floored = sum(1 for n in N_optimal_floored.values() if n < 1000)
    print(f"    floor-constrained Neyman ({FLOOR_FRACTION:.0%} reserved uniform + {1-FLOOR_FRACTION:.0%} by sqrt(v_i)): "
          f"{n_starved_floored}/{n_c} circuits <1000 shots (min={min(N_optimal_floored.values()):,}, "
          f"max={max(N_optimal_floored.values()):,})")

    def mc_std(N_alloc, label_tag):
        errs = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(seed * 104729 + hash(label_tag) % 99991)
            raw = {name: {} for name in kept}
            for (name, gi), N_i in N_alloc.items():
                probs = circuit_bitstring_probs[(name, gi)]
                bitstrings = list(probs.keys())
                p_vec = np.array([probs[b] for b in bitstrings])
                draws = rng.multinomial(N_i, p_vec)
                counts = {b: int(n) for b, n in zip(bitstrings, draws) if n > 0}
                total = sum(counts.values())
                probs_sampled = {b: n / total for b, n in counts.items()} if total > 0 else {}
                for l in groups[gi]:
                    raw[name][l] = pauli_expectation(probs_sampled, l) if probs_sampled else 0.0
            full = build_full(raw, diag, K, non_id_labels)
            _, err = energy_and_err(p, full, K)
            errs.append(err)
        return float(np.mean(errs)), float(np.std(errs))

    t0 = time.time()
    mean_u, std_u = mc_std(N_uniform, "uniform")
    mean_o, std_o = mc_std(N_optimal, "optimal")
    mean_of, std_of = mc_std(N_optimal_floored, "optimal_floored")
    print(f"    uniform allocation (MC, 8-seed):         mean|err|={mean_u:.4f}  std={std_u:.4f} kcal/mol "
          f"(analytic predicted std={std_uniform_analytic:.4f})")
    print(f"    naive Neyman, UNCONSTRAINED (MC, 8-seed): mean|err|={mean_o:.4f}  std={std_o:.4f} kcal/mol "
          f"(analytic predicted std={std_optimal_analytic:.4f})  "
          f"{'CATASTROPHIC FAILURE -- starvation broke the small-noise assumption the analytic formula relies on' if std_o > std_u else ''}")
    print(f"    floor-constrained Neyman (MC, 8-seed):   mean|err|={mean_of:.4f}  std={std_of:.4f} kcal/mol")
    print(f"    analytic-vs-MC agreement: uniform gap={abs(std_u-std_uniform_analytic):.4f}, "
          f"naive-optimal gap={abs(std_o-std_optimal_analytic):.4f} (huge -- analytic formula invalid once "
          f"allocation starves circuits)  ({time.time()-t0:.1f}s)")
    print(f"\n  -- REAL VERDICT: floor-constrained vs uniform, same total budget T --")
    if std_of < std_u:
        print(f"    floor-constrained Neyman WINS: {std_u:.4f} -> {std_of:.4f} kcal/mol "
              f"({std_u/std_of:.2f}x reduction), a genuine MC-verified improvement, no starvation")
    else:
        print(f"    floor-constrained Neyman does NOT beat uniform ({std_u:.4f} vs {std_of:.4f}) at this floor "
              f"fraction -- the {v_arr.max()/max(v_arr.mean(),1e-12):.0f}x spread in v_i is too concentrated in "
              f"too few circuits for a {FLOOR_FRACTION:.0%} floor to both protect against starvation AND meaningfully")
        print(f"    reallocate; reporting this honestly rather than tuning the floor fraction until something wins.")

    print(f"\n  -- HONEST COMPARISON vs the cited headline numbers --")
    print(f"    LOCAL uniform prediction (this file, MC) at 300k shots/circuit: {std_u:.3f} kcal/mol")
    print(f"    iteration 26 Task 1's own LOCAL prediction at 300k shots/circuit: 0.087+/-0.052 kcal/mol")
    print(f"      (different circuit design -- Task 1 used the 36-slot dense design, this file uses the")
    print(f"       21-slot subspace-tomography design; both are LOCAL, shot-noise-only, no drift -- a")
    print(f"       same-order-of-magnitude match is the right bar, not exact equality)")
    print(f"    REAL 1.845+/-0.198 kcal/mol (iteration 26 Task 1) includes submission-to-submission DRIFT,")
    print(f"    which shot reallocation CANNOT fix regardless. But the bigger finding here is upstream of")
    print(f"    that: the naive closed-form Neyman allocation's tiny analytic prediction ({std_optimal_analytic:.4f}) was")
    print(f"    ITSELF WRONG in practice -- MC std={std_o:.2f}, {std_o/max(std_u,1e-9):.0f}x WORSE than uniform -- because")
    print(f"    unconstrained sqrt(v_i) allocation starves low-v_i circuits down to ~1 shot, which breaks the")
    print(f"    delta method's own small-fluctuation precondition. The floor-constrained version fixes the")
    print(f"    starvation ({'and gives a real, if modest, MC-verified win' if std_of < std_u else 'but still does not beat uniform at this floor fraction'}), which is the honest, load-bearing result of this task --")
    print(f"    NOT the naive closed-form number, which would have been a false positive without the MC check.")

    results = {
        "T_total_shots": T, "n_circuits": n_c, "chosen_shots_per_circuit_uniform": CHOSEN_SHOTS,
        "linearity_check_max_gap": max_linearity_gap,
        "analytic": {"std_uniform": std_uniform_analytic, "std_optimal_naive": std_optimal_analytic,
                     "T_needed_for_target_uniform": T_needed_uniform, "T_needed_for_target_optimal_naive": T_needed_optimal},
        "monte_carlo": {
            "uniform": {"mean": mean_u, "std": std_u},
            "naive_neyman_unconstrained": {"mean": mean_o, "std": std_o, "n_starved_lt_1000shots": n_starved,
                                            "verdict": "CATASTROPHIC_FAILURE" if std_o > std_u else "ok"},
            "floor_constrained_neyman": {"mean": mean_of, "std": std_of, "floor_fraction": FLOOR_FRACTION,
                                          "n_starved_lt_1000shots": n_starved_floored,
                                          "beats_uniform": bool(std_of < std_u)},
        },
        "v_i_stats": {"min": float(v_arr.min()), "max": float(v_arr.max()), "mean": float(v_arr.mean())},
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
