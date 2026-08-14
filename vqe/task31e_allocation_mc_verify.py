#!/usr/bin/env python3
"""
task31e_allocation_mc_verify.py -- iteration 31, Task E. The 30/70
floor-constrained allocation was DESIGNED (iteration 30 Task D) but never
Monte-Carlo-verified. Verify it first; then improve using the FULL
Jacobian (sigma_E^2 = J^T Sigma_a J) instead of per-slot max-component
sensitivity, keeping a substantial uniform floor. NEVER unconstrained
Neyman (282x worse under MC in this project's own iteration 28 history).
============================================================================
DESIGN: fixed total shot budget (21 slots x 100,000 = 2,100,000, matching
the current uniform default), redistributed 3 ways:
  (a) UNIFORM       -- the current default, 100,000/slot
  (b) UNCONSTRAINED  -- N_i proportional to full-Jacobian sensitivity
                        ||J_i|| * sqrt(v_i), NO floor -- included
                        specifically to reproduce the known catastrophic-
                        starvation failure mode as a validation check on
                        this script's own MC methodology, not as a
                        candidate scheme
  (c) FLOOR-CONSTRAINED -- 30% of budget reserved uniform, 70% allocated
                        by the SAME full-Jacobian sensitivity, per
                        iteration 28 Task F's precedent
Local, exact ideal probabilities (the established real-hardware-
equivalent multinomial-sampling technique), N_MC_TRIALS=12 -- fewer than
this project's usual 8-32 seed range, disclosed as a scope reduction
given each trial requires a full 21-slot manifold refit (8 restarts
each) x 3 schemes.

Run:
    python vqe/task31e_allocation_mc_verify.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_simulator_binding_curve import stable_seed, expectation_from_counts
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from qiskit.quantum_info import Statevector

K = 6
GATE_NAME = "zz"
TOTAL_SHOTS_PER_SLOT_DEFAULT = 100_000
N_MC_TRIALS = 12
FLOOR_FRACTION = 0.30
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31e_allocation_mc_verify_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task31e_allocation_mc_verify.py -- MC-verify shot allocation, uniform vs Neyman vs floor-constrained")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    print("  building exact ideal probabilities for all 21 slots...")
    exact_probs = {}
    for name in kept:
        base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        exact_probs[name] = {}
        for group in groups:
            combined = effrag_mod.combined_basis_label(group)
            basis_qc = native_basis_change(combined, GATE_NAME)
            qc = base.compose(basis_qc)
            sv = Statevector.from_instruction(qc)
            exact_probs[name][tuple(group)] = sv.probabilities_dict()
    print("  done.")

    # -- full-Jacobian sensitivity per slot: ||J_slot||_2, finite difference on the assembled energy --
    print("\n  -- computing FULL Jacobian sensitivity per slot (not just max component) --")
    a_by_name_ref = {name: target_coeff_vector(name, K) for name in kept}

    def energy_from_a(a_by_name):
        full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
        mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
        E, _ = energy_from_alpha_matrices(mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K)
        return E

    eps = 1e-5
    slot_sensitivity = {}
    for name in kept:
        grad = np.zeros(K)
        a0 = a_by_name_ref[name]
        for j in range(K):
            a_plus = a0.copy(); a_plus[j] += eps; a_plus /= np.linalg.norm(a_plus)
            a_minus = a0.copy(); a_minus[j] -= eps; a_minus /= np.linalg.norm(a_minus)
            tmp = dict(a_by_name_ref)
            tmp[name] = a_plus
            Ep = energy_from_a(tmp)
            tmp[name] = a_minus
            Em = energy_from_a(tmp)
            grad[j] = (Ep - Em) / (2 * eps)
        slot_sensitivity[name] = float(np.linalg.norm(grad))  # FULL Jacobian norm, not max component
    total_sens = sum(slot_sensitivity.values())
    for name in sorted(slot_sensitivity, key=lambda n: -slot_sensitivity[n])[:5]:
        print(f"    {name}: ||J_slot||={slot_sensitivity[name]:.4f}")

    # -- 3 allocation schemes, fixed total budget --
    n_slots = len(kept)
    total_budget = n_slots * TOTAL_SHOTS_PER_SLOT_DEFAULT
    schemes = {}
    schemes["uniform"] = {name: TOTAL_SHOTS_PER_SLOT_DEFAULT for name in kept}
    unconstrained = {name: max(1, int(total_budget * slot_sensitivity[name] / total_sens)) for name in kept}
    schemes["unconstrained_neyman"] = unconstrained
    floor_share = FLOOR_FRACTION * total_budget / n_slots
    remaining = (1 - FLOOR_FRACTION) * total_budget
    floor_constrained = {name: int(floor_share + remaining * slot_sensitivity[name] / total_sens) for name in kept}
    schemes["floor_constrained_30_70"] = floor_constrained

    print(f"\n  allocation schemes (total budget={total_budget:,} shots across {n_slots} slots):")
    for scheme_name, alloc in schemes.items():
        vals = list(alloc.values())
        print(f"    {scheme_name}: min={min(vals):,}  max={max(vals):,}  "
              f"n_below_1000={sum(1 for v in vals if v < 1000)}")

    # -- MC verification: for each scheme, resample at the ALLOCATED shot count, refit, measure final-energy variance --
    print(f"\n  -- MC verification, {N_MC_TRIALS} trials per scheme (ideal model, exact probabilities) --")
    results = {}
    for scheme_name, alloc in schemes.items():
        errs = []
        for trial in range(N_MC_TRIALS):
            manifold_kept = {}
            for name in kept:
                n_shots = alloc[name]
                m_dict = {}
                for group_t, probs in exact_probs[name].items():
                    bitstrings = list(probs.keys())
                    parr = np.array([probs[b] for b in bitstrings])
                    parr = parr / parr.sum()
                    rng = np.random.default_rng(stable_seed("task31e", scheme_name, trial, name, group_t))
                    draws = rng.multinomial(n_shots, parr)
                    counts = {b: int(c) for b, c in zip(bitstrings, draws) if c > 0}
                    for l in group_t:
                        m_dict[l] = expectation_from_counts(counts, l) if counts else 0.0
                v0 = target_coeff_vector(name, K)
                w_dict = {l: 1.0 for l in m_dict}  # uniform weighting within the fit itself -- allocation is the
                                                    # variable under test here, not the fit's internal weighting
                a_hat, _ = fit_pure_state(P_S, m_dict, w_dict, K, v0,
                                           seed=stable_seed("task31e_fit", scheme_name, trial, name))
                manifold_kept[name] = a_hat
            full = build_full_from_a(manifold_kept, P_S, diag, K, non_id_labels)
            mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
            E, errs_d = energy_from_alpha_matrices(mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                     exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs.append(errs_d["err_vs_exact_kcal"])
        mean_err = float(np.mean(errs))
        std_err = float(np.std(errs))
        results[scheme_name] = {"mean": mean_err, "std": std_err, "errs": errs}
        print(f"    {scheme_name}: mean|err|={mean_err:.4f}  std={std_err:.4f} kcal/mol (n={N_MC_TRIALS} trials)")

    print("\n" + "=" * 96)
    print("  VERDICT:")
    uniform_std = results["uniform"]["std"]
    for scheme_name in ["unconstrained_neyman", "floor_constrained_30_70"]:
        ratio = results[scheme_name]["std"] / uniform_std if uniform_std > 1e-9 else float("inf")
        print(f"    {scheme_name}: std/uniform_std = {ratio:.2f}x "
              f"({'WORSE (as expected/validation check)' if ratio > 1 else 'BETTER'})")
    print("=" * 96 + "\n")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"slot_sensitivity": slot_sensitivity, "schemes": {k: {n: v for n, v in alloc.items()}
                   for k, alloc in schemes.items()}, "results": results}, f, indent=2)
    print(f"  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
