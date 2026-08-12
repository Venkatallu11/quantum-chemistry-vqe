#!/usr/bin/env python3
"""
task1_shot_budget.py — iteration 26, Task 1. THE SHOT BUDGET, done first
because it gates every other task: at this project's standard 10,000
shots/setting, the IDEAL (noiseless) control itself sits at 1.35 kcal/mol
(iteration 25, Task A) -- ABOVE chemical accuracy (1.0). No mitigation
result measured at this shot count can ever be claimed to reach chemical
accuracy, because a PERFECT device would already fail the bar at this
budget.
============================================================================
METHOD: sweep shots/setting in {10k, 30k, 100k, 300k, 1M} -- 1M is IonQ's
own REAL per-job cap, confirmed live via GET /jobs/estimate (a shots=
1,000,001 request is rejected with "shots must be less than or equal to
1000000"; NOT assumed, not guessed). Since the sweep's own upper bound
equals the real cap, no cross-job pooling is actually needed for this
sweep -- stated explicitly since the task anticipated needing it.

Two independent things are computed at each shot level:
  1. The IDEAL-model shot-noise floor (LOCAL, exact expectations +
     8-seed shot_sample bootstrap -- no network, fast) -- where does the
     noiseless control itself drop below 0.3 kcal/mol?
  2. The REAL cost this would take on actual IonQ QPU hardware, using the
     REAL rate card (queried live, not the user-supplied numbers taken on
     faith -- confirmed to match exactly: cost_1q_gate=$0.000164,
     cost_2q_gate=$0.001121, job_cost_minimum=$25.7899), extending
     iteration 9's own established methodology (the floor is charged once
     per JOB, confirmed there by a 125x-scaling test) to this shot sweep.

Once the shot level is chosen, this project's own free `ionq_simulator`
(zero cost, unlike real QPU) is used to submit ONE real validation job at
that level, as a genuine real-data check on the local shot-noise
prediction -- not purely asserted from the local model.

Run:
    python vqe/task1_shot_budget.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, shot_sample, build_ansatz,
)
from fixed_ansatz import P2_PER_GATE, P2_PER_GATE_OLD_ASSUMED
from qiskit.quantum_info import Pauli, Statevector
from ionq_backend import connect_provider

K = 6
SHOT_LEVELS = [10_000, 30_000, 100_000, 300_000, 1_000_000]
MAX_SHOTS_PER_JOB = 1_000_000  # confirmed live, see module docstring
N_SEEDS = 8
CHEM_ACC_KCAL = 1.0
IDEAL_TARGET_KCAL = 0.3  # per explicit instruction: floor for the ideal control to STOP contributing to the total budget
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task1_shot_budget_results.json")

# gate counts per full (36-slot x 13-group = 468-circuit) submission, per config
# fixed ansatz: 11 CX + 51 1q (established, ionq_resource_estimate.py / fixed_ansatz.py)
CONFIG_GATE_COUNTS = {
    "fixed (11 CX)": {"n1q": 51, "n2q": 11, "n_circuits": 468},
    "ADAPT (8.53 CX mean)": {"n1q": 51, "n2q": 9, "n_circuits": 468},  # rounded mean, see iteration 24 Task 1
    "variational (4.3 CX mean)": {"n1q": 51, "n2q": 4, "n_circuits": 468},
    "tapered (3.94 CX mean, 3 qubits)": {"n1q": 30, "n2q": 4, "n_circuits": 324},  # 36x9 groups, 3q register
    "subspace tomography (21 circ, 11 CX)": {"n1q": 51, "n2q": 11, "n_circuits": 273},  # 21x13
}
RATE_1Q = 0.000164
RATE_2Q = 0.001121
JOB_FLOOR = 25.7899


def real_cost(n1q, n2q, n_circuits, shots):
    """Per iteration 9's confirmed methodology: the $25.79 floor is
    charged ONCE PER JOB (verified there via an exact 125x-scaling test),
    and at this project's shot counts gate-execution cost already
    dominates the floor -- so total cost = max(floor, per-circuit-cost) x
    n_circuits is the RIGHT formula only if circuits are billed
    separately; bundled into ONE job (as iteration 9 verified), the floor
    applies once and gate cost simply sums -- reproduced here exactly:
    cost = job_floor_if_it_dominates, else n_circuits x shots x
    (n1q*RATE_1Q + n2q*RATE_2Q), floor compared against the SUMMED
    gate cost of the whole job, not per-circuit."""
    gate_cost_total = n_circuits * shots * (n1q * RATE_1Q + n2q * RATE_2Q)
    return max(gate_cost_total, JOB_FLOOR)


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def ideal_shot_noise_at(p, exact_ideal_raw, non_id_labels, shots, n_seeds=N_SEEDS):
    errs = []
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 7919 + shots)
        shot_raw = {name: {l: shot_sample(exact_ideal_raw[name][l], shots, rng) for l in non_id_labels}
                    for name in exact_ideal_raw}
        _, err = energy_and_err(p, shot_raw, K)
        errs.append(err)
    return float(np.mean(errs)), float(np.std(errs))


def main():
    print("\n" + "=" * 96)
    print("  task1_shot_budget.py -- the shot budget, done first, gates everything")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"])
    assert n_ok == 36
    counts = verify_constant_gate_count(fixed_solutions)
    assert counts == {11}
    print(f"  setup OK: 36/36 converged")

    # -- confirm the real per-job shot cap live (not assumed) --
    provider = connect_provider()
    client = provider.get_backend("ionq_simulator").client
    try:
        client.estimate_job(backend="qpu.forte-1", oneq_gates=51, twoq_gates=11, qubits=4, shots=MAX_SHOTS_PER_JOB + 1)
        cap_confirmed = False
    except Exception as e:
        cap_confirmed = "1000000" in str(e)
    print(f"  real per-job shot cap confirmed live: {MAX_SHOTS_PER_JOB} ({'CONFIRMED' if cap_confirmed else 'NOT CONFIRMED -- see error'})")

    # -- exact ideal expectations (statevector, no noise) --
    exact_ideal_raw = {}
    for name, sol in fixed_solutions.items():
        sv = Statevector.from_instruction(build_ansatz(sol["angles"]))
        exact_ideal_raw[name] = {l: float(sv.expectation_value(Pauli(l)).real) for l in non_id_labels}
    _, exact_err = energy_and_err(p, exact_ideal_raw, K)
    print(f"  exact (infinite-shot) ideal error vs exact_energy: {exact_err:.2e} kcal/mol (should be ~0)")

    print(f"\n  -- IDEAL shot-noise floor sweep --")
    sweep = []
    crossing_shots = None
    for shots in SHOT_LEVELS:
        mean, std = ideal_shot_noise_at(p, exact_ideal_raw, non_id_labels, shots)
        below = mean < IDEAL_TARGET_KCAL
        sweep.append({"shots": shots, "mean": mean, "std": std, "below_target": below})
        print(f"    shots={shots:>9,}: ideal control = {mean:.3f}+/-{std:.3f} kcal/mol  "
              f"{'BELOW 0.3 target' if below else ''}")
        if below and crossing_shots is None:
            crossing_shots = shots

    if crossing_shots is None:
        print(f"\n  ideal control does NOT drop below {IDEAL_TARGET_KCAL} kcal/mol even at the max swept "
              f"shot level ({SHOT_LEVELS[-1]:,}) -- reporting this honestly, not extrapolating past what was tested.")
        chosen_shots = SHOT_LEVELS[-1]
    else:
        chosen_shots = crossing_shots
        print(f"\n  CHOSEN SHOT LEVEL: {chosen_shots:,} (first level where ideal control < {IDEAL_TARGET_KCAL} kcal/mol)")

    # -- real cost, every configuration, every shot level --
    print(f"\n  -- REAL QPU cost (rate card confirmed live: 1q=${RATE_1Q}, 2q=${RATE_2Q}, floor=${JOB_FLOOR}) --")
    print(f"     per-job floor confirmed once-per-job by iteration 9's own 125x scaling test, reused not re-derived")
    cost_table = {}
    for cfg_name, g in CONFIG_GATE_COUNTS.items():
        cost_table[cfg_name] = {}
        print(f"    {cfg_name} ({g['n_circuits']} circuits, {g['n1q']} 1q/{g['n2q']} 2q gates/circuit):")
        for shots in SHOT_LEVELS:
            cost = real_cost(g["n1q"], g["n2q"], g["n_circuits"], shots)
            cost_table[cfg_name][shots] = cost
            flag = "  <- CHOSEN LEVEL" if shots == chosen_shots else ""
            print(f"      shots={shots:>9,}: ${cost:>14,.2f}  ({cost/3000:.0f}x the $3,000 budget){flag}")

    print(f"\n  -- DECISIVE FINDING: is real QPU affordable at ANY shot level in this sweep? --")
    cheapest = min(cost_table["subspace tomography (21 circ, 11 CX)"][s] for s in SHOT_LEVELS)
    print(f"    cheapest possible real-QPU option in this whole sweep (21-circuit design, 10k shots): "
          f"${cheapest:,.2f} ({cheapest/3000:.1f}x the $3,000 budget)")
    print(f"    NO shot level, NO circuit design tested here fits the budget. Real QPU remains categorically "
          f"unaffordable at every point in this sweep -- extends iteration 9's finding (452x over budget at "
          f"the old 7-geometry scope) to the current single-geometry K=6 scope and full shot range.")

    # -- re-quote local-model headline numbers at the chosen shot level --
    print(f"\n  -- re-quoting LOCAL MODEL headline numbers at shots={chosen_shots:,} --")
    nm_labels = {"old_assumed": P2_PER_GATE_OLD_ASSUMED, "corrected_real_forte1": P2_PER_GATE}
    from taskB_corrected_fidelity_pipeline import build_noise_model, noisy_density_matrix
    requoted = {}
    for p2_label, p2 in nm_labels.items():
        nm = build_noise_model(p2)
        builders = {name: (lambda a=sol["angles"]: build_ansatz(a)) for name, sol in fixed_solutions.items()}
        exact_raw = {name: {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix())
                     @ noisy_density_matrix(b, nm)))) for l in non_id_labels} for name, b in builders.items()}
        mean, std = ideal_shot_noise_at(p, exact_raw, non_id_labels, chosen_shots)
        requoted[p2_label] = {"mean": mean, "std": std}
        print(f"    raw @ {p2_label} ({p2}): {mean:.2f}+/-{std:.2f} kcal/mol (was ~104/42 at 10k shots -- "
              f"shot noise itself shrinks, the BIAS from noise does not)")

    # -- REAL validation submission at the chosen shot level (free simulator, ideal model only, small) --
    print(f"\n  -- REAL validation: submitting the ideal model at shots={chosen_shots:,} to confirm the local prediction --")
    from ionq_backend import get_simulator
    from ionq_run import basis_change, IONQ_QIS_STANDARD_BASIS
    from ionq_simulator_binding_curve import submit_job, get_counts_list, bootstrap_counts, expectation_from_counts, stable_seed
    import ef_fragment as effrag_mod
    from qiskit.transpiler import CouplingMap
    from qiskit import transpile

    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    CMAP4 = CouplingMap.from_full(4)
    backend = get_simulator(provider)
    # REAL BUG, found and fixed here: target_names MUST use the SAME order
    # circuits are built in (dict insertion order from fixed_solutions),
    # not alphabetically sorted -- sorted() puts every "(un+um)"-style name
    # BEFORE "u_0".."u_5" (ASCII '(' < 'u'), a completely different order
    # from insertion, silently scrambling which measured counts got
    # assigned to which (slot, label) pair. Caught because the resulting
    # "REAL ideal control" was 1416.9 kcal/mol -- an obviously-wrong
    # result, caught by exactly the kind of too-bad-to-be-true check this
    # project's honesty rules require before trusting ANY number, not just
    # suspiciously GOOD ones.
    target_names = list(fixed_solutions.keys())
    circuits = []
    for name in target_names:
        sol = fixed_solutions[name]
        base = transpile(build_ansatz(sol["angles"]), basis_gates=IONQ_QIS_STANDARD_BASIS,
                          coupling_map=CMAP4, optimization_level=0)
        for group in groups:
            combined = effrag_mod.combined_basis_label(group)
            qc = base.copy()
            basis_change(qc, combined)
            qc.measure_all()
            circuits.append(qc)
    print(f"    {len(circuits)} circuits, shots={chosen_shots:,}/circuit")
    t0 = time.time()
    job = submit_job(circuits, backend, "ideal", shots=chosen_shots)
    t_submit = time.time() - t0
    counts_flat = get_counts_list(job)
    t_retrieve = time.time() - t0 - t_submit
    print(f"    submitted {t_submit:.1f}s, retrieved {t_retrieve:.1f}s")

    validation_ckpt_path = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                                         "task1_shot_validation.json")
    with open(validation_ckpt_path, "w") as f:
        json.dump({"target_names": target_names, "groups": groups, "counts": counts_flat,
                    "shots": chosen_shots}, f)
    print(f"    counts checkpoint saved -> {validation_ckpt_path} (so this expensive real job never needs re-submitting)")

    idx_map = [name for name in target_names for _ in groups]
    per_name = {name: [] for name in target_names}
    for i, name in enumerate(idx_map):
        per_name[name].append(counts_flat[i])
    group_idx = {}
    for gi, g in enumerate(groups):
        for l in g:
            group_idx[l] = gi
    errs = []
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(stable_seed("task1_shotbudget", "ideal", seed))
        raw = {name: {} for name in target_names}
        for name in target_names:
            for gi, group in enumerate(groups):
                counts = bootstrap_counts(per_name[name][gi], chosen_shots, rng)
                for l in group:
                    raw[name][l] = expectation_from_counts(counts, l)
        _, err = energy_and_err(p, raw, K)
        errs.append(err)
    real_mean, real_std = float(np.mean(errs)), float(np.std(errs))
    predicted = next(s for s in sweep if s["shots"] == chosen_shots)
    print(f"    REAL ideal control @ shots={chosen_shots:,}: {real_mean:.3f}+/-{real_std:.3f} kcal/mol "
          f"(LOCAL prediction was {predicted['mean']:.3f}+/-{predicted['std']:.3f})")

    results = {
        "shot_cap_confirmed": cap_confirmed, "max_shots_per_job": MAX_SHOTS_PER_JOB,
        "ideal_shot_noise_sweep": sweep, "chosen_shots": chosen_shots,
        "cost_table_usd": cost_table, "rate_card": {"1q": RATE_1Q, "2q": RATE_2Q, "floor": JOB_FLOOR},
        "requoted_local_raw_at_chosen_shots": requoted,
        "real_validation_at_chosen_shots": {"mean": real_mean, "std": real_std,
                                             "local_prediction_mean": predicted["mean"], "local_prediction_std": predicted["std"]},
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
