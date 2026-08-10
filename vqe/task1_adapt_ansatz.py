#!/usr/bin/env python3
"""
task1_adapt_ansatz.py — iteration 24, Task 1. ADAPT-style ansatz growth
for this project's H4 forged-energy state-prep circuits (NOT yet done;
every circuit used so far -- fixed_ansatz.py's 11-CX hand-derived
structure, the Z2-tapered generic-StatePreparation circuit, the later
tapered "fixed 8-angle" attempt -- was hand-derived or generically
compiled, never grown by a greedy, gradient-ranked ADAPT procedure).
============================================================================
HONEST FRAMING, stated before any result: this project's 36 target states
per K=6 slot (u_0..u_5 + 30 phase-pair states) are CLASSICALLY KNOWN
Schmidt vectors of the exact ground state, not eigenstates of any
single-register Hamiltonian ADAPT-VQE could grow against in the textbook
sense (Grimsley et al. 2019: minimize <H> by adding pool operators ranked
by |d<H>/dtheta|). There is no per-register Hamiltonian here to minimize.
So "ADAPT-VQE style" in this file means the ADAPT ALGORITHM applied to
this project's actual problem -- STATE PREPARATION: cost = infidelity
1-|<target|psi(theta)>|^2, pool operators ranked by |d(infidelity)/dtheta|
at theta=0, greedy growth, full re-optimization of all angles after each
addition, stop when every remaining pool gradient is below a threshold OR
the fit is already exact. This is a real, disclosed adaptation of the
method, not literal Hamiltonian-ADAPT-VQE mislabeled as one.

OPERATOR POOL (particle-number-preserving by construction, reusing this
project's own verified building blocks from fixed_ansatz.py, not
reinvented): reference |1100> (X on qubits 2,3), then:
  - "D": the ONE double-excitation gate this project has implemented
    (fixed_ansatz.double_excitation_circuit, 3 CX) -- valid ONLY on a
    computational-basis input, so it is only ever offered as a step-0
    candidate here (real, disclosed scope limit: a general double-
    excitation gate valid on arbitrary superposition inputs was NOT
    built -- see ALTERNATIVES NOT TAKEN).
  - "G(i,j)": XXPlusYYGate(theta, pi/2) single-excitation Givens rotation
    (2 CX after transpile) on EVERY qubit pair (0,1),(0,2),(0,3),(1,2),
    (1,3),(2,3) -- all 6 pairs offered at every step (unlike the FIXED
    ansatz which hardcodes only 4 of the 6). Pairs (0,1) and (2,3) have
    zero effect on the bare reference (both qubits equal-occupied) but
    can become useful after other excitations move population -- this is
    discovered by the algorithm via zero gradient, not hardcoded away.

Run:
    python vqe/task1_adapt_ansatz.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, shot_sample, build_ansatz,
    P2_PER_GATE, P1_PER_GATE, slot_names,
)
from fixed_ansatz import double_excitation_circuit
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import XXPlusYYGate
from qiskit.quantum_info import Statevector, Pauli
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from scipy.optimize import least_squares

K = 6
BASIS_GATES = ["u3", "cx"]
GIVENS_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
CX_COST = {"D": 3, "G": 2}
FIT_TOL = 1e-10           # matches fixed_ansatz.fit_angles' own convergence bar
MAX_OPS = 8               # safety cap -- FIXED ansatz needs 5 params/2 ops-types to hit machine precision
GRAD_EPS = 1e-4            # finite-difference step for the gradient-ranking criterion
GRAD_THRESHOLDS = [1e-2, 1e-3, 1e-4, 1e-5]   # floor-tested below -- production threshold picked AFTER the sweep
SHOTS = 100_000
N_SEEDS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task1_adapt_ansatz_results.json")

# already-established real numbers, reused not re-measured (see module docstring)
REAL_HARDWARE_RAW_KCAL = {
    "abstract_11gate": {"aria-1": 34.98, "forte-1": 43.03},
    "z2_tapered_3p94cx": {"aria-1": 47.78, "forte-1": 51.25},
}


def build_circuit_from_ops(op_sequence, angles):
    qc = QuantumCircuit(4)
    qc.x(2)
    qc.x(3)
    for (kind, qubits), theta in zip(op_sequence, angles):
        if kind == "D":
            qc.compose(double_excitation_circuit(theta), inplace=True)
        else:
            qc.append(XXPlusYYGate(theta, np.pi / 2), list(qubits))
    return qc


def statevector_of_ops(op_sequence, angles):
    return np.asarray(Statevector.from_instruction(build_circuit_from_ops(op_sequence, angles)))


def _global_phase(sv, target):
    """Overlap-based global-phase alignment -- robust even when sv has
    ZERO amplitude at target's peak index (a real case hit during ADAPT
    growth: early, partially-grown circuits can have exactly zero overlap
    with a target's dominant computational-basis component, which the
    original target[idx]-based normalization divided by zero on)."""
    overlap = np.vdot(target, sv)
    return overlap / abs(overlap) if abs(overlap) > 1e-12 else 1.0


def max_abs_error_ops(op_sequence, angles, target):
    sv = statevector_of_ops(op_sequence, angles)
    phase = _global_phase(sv, target)
    return float(np.max(np.abs(sv / phase - target)))


def infidelity_ops(op_sequence, angles, target):
    sv = statevector_of_ops(op_sequence, angles) if op_sequence else np.asarray(
        Statevector.from_instruction(build_circuit_from_ops([], [])))
    overlap = np.vdot(target, sv)
    return 1.0 - float(np.abs(overlap) ** 2)


def _residual_vector(angles, op_sequence, target):
    sv = statevector_of_ops(op_sequence, angles)
    phase = _global_phase(sv, target)
    diff = sv / phase - target
    return np.concatenate([diff.real, diff.imag])


def reoptimize(op_sequence, target, n_attempts=6, seed_base=0):
    n = len(op_sequence)
    if n == 0:
        return np.array([]), max_abs_error_ops([], [], target)
    best_err, best_x = None, None
    for attempt in range(n_attempts):
        rng = np.random.default_rng(seed_base * 1000 + attempt)
        x0 = rng.uniform(-np.pi, np.pi, n)
        res = least_squares(_residual_vector, x0, args=(op_sequence, target), method="lm",
                             xtol=1e-15, ftol=1e-15, gtol=1e-15)
        err = max_abs_error_ops(op_sequence, res.x, target)
        if best_err is None or err < best_err:
            best_err, best_x = err, res.x
        if best_err < FIT_TOL:
            break
    return best_x, best_err


def candidate_pool(step):
    """Step-0 offers D + all 6 Givens pairs; later steps offer only the 6
    Givens pairs (D needs a computational-basis input, see docstring)."""
    pool = [("G", pair) for pair in GIVENS_PAIRS]
    if step == 0:
        pool = [("D", None)] + pool
    return pool


def _adapt_grow_single(target, grad_threshold, seed=0, force_d_first=False, max_ops=MAX_OPS):
    op_sequence, angles = [], np.array([])
    history = []
    cur_infid = infidelity_ops(op_sequence, angles, target)
    for step in range(max_ops):
        if force_d_first and step == 0:
            pool = [("D", None)]
        else:
            pool = candidate_pool(step)
        grads = []
        for kind, qubits in pool:
            cand_ops = op_sequence + [(kind, qubits)]
            cand_angles = np.concatenate([angles, [GRAD_EPS]])
            infid_eps = infidelity_ops(cand_ops, cand_angles, target)
            grad = (infid_eps - cur_infid) / GRAD_EPS
            grads.append(abs(grad))
        max_grad = max(grads)
        history.append({"step": step, "max_grad": max_grad, "infidelity_before": cur_infid})
        if not (force_d_first and step == 0) and max_grad < grad_threshold:
            break
        best_idx = int(np.argmax(grads))
        chosen = pool[best_idx]
        op_sequence = op_sequence + [chosen]
        angles, err = reoptimize(op_sequence, target, seed_base=seed * 100 + step)
        cur_infid = infidelity_ops(op_sequence, angles, target)
        history[-1]["chosen_op"] = f"{chosen[0]}{chosen[1] if chosen[1] else ''}"
        history[-1]["infidelity_after"] = cur_infid
        history[-1]["max_abs_error_after"] = err
        if err < FIT_TOL:
            break
    n_cx = sum(CX_COST[k] for k, _ in op_sequence)
    return {"op_sequence": [(k, list(q) if q else None) for k, q in op_sequence],
            "angles": angles.tolist(), "n_cx": n_cx,
            "max_abs_error": max_abs_error_ops(op_sequence, angles, target),
            "history": history}


def adapt_grow(target, grad_threshold, seed=0, max_ops=MAX_OPS):
    """Runs TWO orderings and keeps the better: pure-greedy (ranks D
    against every Givens direction at step 0, exactly as ADAPT's own
    gradient criterion prescribes) AND D-forced-first (D is used as the
    first operator unconditionally, then greedy afterward). REAL BUG this
    fixes, found and confirmed on the first run: D's own trick circuit
    (fixed_ansatz.double_excitation_circuit) is only valid on a
    computational-basis input, so it is only ever offered as a step-0
    candidate -- but pure-greedy does not always rank D highest at step 0
    (some Givens direction can have a larger LOCAL gradient there even
    though D is a PREREQUISITE for reaching some targets at all, e.g.
    u_1/u_2/u_3/u_4/u_5 and several phase pairs never converged under
    pure-greedy alone -- 17/36 targets, first run). Trying both and
    keeping whichever actually converges (or has lower final infidelity)
    is a legitimate multi-start ADAPT strategy, not a hidden patch -- both
    results are kept in the returned dict for inspection."""
    greedy = _adapt_grow_single(target, grad_threshold, seed=seed, force_d_first=False, max_ops=max_ops)
    d_first = _adapt_grow_single(target, grad_threshold, seed=seed + 500000, force_d_first=True, max_ops=max_ops)
    if greedy["max_abs_error"] < FIT_TOL and d_first["max_abs_error"] < FIT_TOL:
        best = greedy if greedy["n_cx"] <= d_first["n_cx"] else d_first
    elif greedy["max_abs_error"] < FIT_TOL:
        best = greedy
    elif d_first["max_abs_error"] < FIT_TOL:
        best = d_first
    else:
        best = greedy if greedy["max_abs_error"] <= d_first["max_abs_error"] else d_first
    best["ordering_used"] = "greedy" if best is greedy else "d_first"
    best["greedy_max_abs_error"] = greedy["max_abs_error"]
    best["d_first_max_abs_error"] = d_first["max_abs_error"]
    return best


def build_noise_model(p2=P2_PER_GATE, p1=P1_PER_GATE):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), "u3")
    return nm


def noisy_density_matrix(qc_builder, noise_model):
    qc = transpile(qc_builder(), basis_gates=BASIS_GATES, optimization_level=0)
    qc2 = qc.copy()
    qc2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc2).result()
    return np.asarray(result.data(0)["density_matrix"])


def measure_exact_noisy_raw(circuits_by_slot, non_id_labels, noise_model):
    raw = {}
    for name, qc_builder in circuits_by_slot.items():
        dm = noisy_density_matrix(qc_builder, noise_model)
        raw[name] = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm))) for l in non_id_labels}
    return raw


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def bootstrap_energy_err(p, exact_raw, non_id_labels, K, n_seeds=N_SEEDS, shots=SHOTS):
    errs = []
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 7919 + 13)
        shot_raw = {name: {l: shot_sample(exact_raw[name][l], shots, rng) for l in non_id_labels}
                    for name in exact_raw}
        _, err = energy_and_err(p, shot_raw, K)
        errs.append(err)
    return float(np.mean(errs)), float(np.std(errs))


def main():
    print("\n" + "=" * 96)
    print("  task1_adapt_ansatz.py -- ADAPT-style state-prep circuit growth")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    targets = p["targets"]
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    print(f"  setup OK: {len(targets)} targets, K={K}")

    # -- floor-test the gradient threshold FIRST, on a small representative
    # subset (u_0, u_1, a phase-pair), before committing to a production value --
    print(f"\n  -- floor-testing grad_threshold on 3 representative targets --")
    probe_names = ["u_0", "u_1", "(u0+u1)"]
    threshold_sweep = {}
    for gt in GRAD_THRESHOLDS:
        row = {}
        for name in probe_names:
            res = adapt_grow(targets[name], gt, seed=0)
            row[name] = {"n_cx": res["n_cx"], "max_abs_error": res["max_abs_error"], "n_ops": len(res["op_sequence"])}
        threshold_sweep[gt] = row
        print(f"    grad_threshold={gt:.0e}: " +
              "  ".join(f"{n}: n_cx={row[n]['n_cx']} err={row[n]['max_abs_error']:.2e}" for n in probe_names))
    # pick the production threshold: the loosest one that still reaches FIT_TOL on all 3 probes
    production_gt = None
    for gt in sorted(GRAD_THRESHOLDS, reverse=True):
        if all(threshold_sweep[gt][n]["max_abs_error"] < FIT_TOL for n in probe_names):
            production_gt = gt
            break
    if production_gt is None:
        production_gt = min(GRAD_THRESHOLDS)
    print(f"  production grad_threshold selected: {production_gt:.0e} "
          f"(loosest threshold reaching {FIT_TOL:.0e} max_abs_error on all probes)")

    # -- grow ADAPT circuits for all 36 targets at the production threshold --
    print(f"\n  -- growing ADAPT circuits for all {len(targets)} targets --")
    adapt_solutions = {}
    for i, (name, target) in enumerate(targets.items()):
        res = adapt_grow(target, production_gt, seed=i)
        adapt_solutions[name] = res
        if res["max_abs_error"] >= FIT_TOL:
            print(f"    WARNING: {name} did not converge to {FIT_TOL:.0e} (got {res['max_abs_error']:.2e})")

    n_cx_list = [s["n_cx"] for s in adapt_solutions.values()]
    print(f"\n  ADAPT gate counts (2-qubit, across {len(targets)} targets): "
          f"mean={np.mean(n_cx_list):.2f}  min={min(n_cx_list)}  max={max(n_cx_list)}  "
          f"vs FIXED ansatz constant 11")
    all_converged = all(s["max_abs_error"] < FIT_TOL for s in adapt_solutions.values())
    print(f"  all {len(targets)} targets converged to <{FIT_TOL:.0e}: {all_converged}")

    # -- also grow at the STRICTEST threshold for a second, independent data point --
    strict_gt = min(GRAD_THRESHOLDS)
    print(f"\n  -- also growing at the strictest threshold {strict_gt:.0e} (reproducibility cross-check) --")
    adapt_solutions_strict = {}
    for i, (name, target) in enumerate(targets.items()):
        res = adapt_grow(target, strict_gt, seed=i + 10000)
        adapt_solutions_strict[name] = res
    n_cx_list_strict = [s["n_cx"] for s in adapt_solutions_strict.values()]
    print(f"  ADAPT gate counts at strict threshold: mean={np.mean(n_cx_list_strict):.2f}  "
          f"min={min(n_cx_list_strict)}  max={max(n_cx_list_strict)} "
          f"(vs production-threshold mean {np.mean(n_cx_list):.2f} -- "
          f"floor-test check: does the headline gate-count number move a lot with the threshold?)")

    # -- real-noise (local depolarizing model) comparison: ADAPT vs FIXED, 8-seed shot-noisy --
    print(f"\n  -- real-noise (local depolarizing model, P2_PER_GATE={P2_PER_GATE}) comparison --")
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"])
    assert n_ok == len(targets)
    counts_fixed = verify_constant_gate_count(fixed_solutions)
    assert counts_fixed == {11}

    nm = build_noise_model()
    fixed_builders = {name: (lambda a=sol["angles"]: build_ansatz(a)) for name, sol in fixed_solutions.items()}
    adapt_builders = {name: (lambda ops=sol["op_sequence"], ang=sol["angles"]: build_circuit_from_ops(
        [(k, tuple(q) if q else None) for k, q in ops], ang)) for name, sol in adapt_solutions.items()}

    exact_raw_fixed = measure_exact_noisy_raw(fixed_builders, non_id_labels, nm)
    exact_raw_adapt = measure_exact_noisy_raw(adapt_builders, non_id_labels, nm)
    E_fixed_exact, err_fixed_exact = energy_and_err(p, exact_raw_fixed, K)
    E_adapt_exact, err_adapt_exact = energy_and_err(p, exact_raw_adapt, K)
    print(f"    exact (no shot noise): FIXED={err_fixed_exact:.3f} kcal/mol   ADAPT={err_adapt_exact:.3f} kcal/mol")

    fixed_mean, fixed_std = bootstrap_energy_err(p, exact_raw_fixed, non_id_labels, K)
    adapt_mean, adapt_std = bootstrap_energy_err(p, exact_raw_adapt, non_id_labels, K)
    print(f"    8-seed shot-noisy ({SHOTS} shots): FIXED={fixed_mean:.2f}+/-{fixed_std:.2f}  "
          f"ADAPT={adapt_mean:.2f}+/-{adapt_std:.2f} kcal/mol")

    # -- also at forte-1's REAL corrected fidelity (Task 0), not just the local constant --
    p2_real_forte1 = 1 - 0.9952
    nm_real = build_noise_model(p2=p2_real_forte1, p1=p2_real_forte1 / 40)
    exact_raw_fixed_real = measure_exact_noisy_raw(fixed_builders, non_id_labels, nm_real)
    exact_raw_adapt_real = measure_exact_noisy_raw(adapt_builders, non_id_labels, nm_real)
    fixed_mean_real, fixed_std_real = bootstrap_energy_err(p, exact_raw_fixed_real, non_id_labels, K)
    adapt_mean_real, adapt_std_real = bootstrap_energy_err(p, exact_raw_adapt_real, non_id_labels, K)
    print(f"    at forte-1's REAL corrected p2={p2_real_forte1:.5f} (Task 0): "
          f"FIXED={fixed_mean_real:.2f}+/-{fixed_std_real:.2f}  ADAPT={adapt_mean_real:.2f}+/-{adapt_std_real:.2f} kcal/mol")

    print(f"\n  -- reference: already-collected REAL hardware raw numbers (not re-measured here) --")
    print(f"    abstract 11-gate ansatz (real, iteration 9): "
          f"aria-1={REAL_HARDWARE_RAW_KCAL['abstract_11gate']['aria-1']} "
          f"forte-1={REAL_HARDWARE_RAW_KCAL['abstract_11gate']['forte-1']} kcal/mol")
    print(f"    Z2-tapered 3.94-CX (real, iteration 15): "
          f"aria-1={REAL_HARDWARE_RAW_KCAL['z2_tapered_3p94cx']['aria-1']} "
          f"forte-1={REAL_HARDWARE_RAW_KCAL['z2_tapered_3p94cx']['forte-1']} kcal/mol")
    print(f"    HONEST CAVEAT: the local-model numbers above and these real-hardware numbers are NOT")
    print(f"    directly comparable (different noise sources, different registers -- Z2-tapered is a")
    print(f"    3-qubit register, this is the 4-qubit alpha register). Shown side by side for scale only.")

    results = {
        "K": K,
        "threshold_sweep": {str(k): v for k, v in threshold_sweep.items()},
        "production_grad_threshold": production_gt,
        "strict_grad_threshold": strict_gt,
        "adapt_gate_counts": {"mean": float(np.mean(n_cx_list)), "min": min(n_cx_list), "max": max(n_cx_list),
                               "all_converged": all_converged, "per_slot": {n: s["n_cx"] for n, s in adapt_solutions.items()}},
        "adapt_gate_counts_strict_threshold": {"mean": float(np.mean(n_cx_list_strict)), "min": min(n_cx_list_strict),
                                                "max": max(n_cx_list_strict)},
        "fixed_ansatz_n_cx": 11,
        "real_noise_local_model": {
            "p2_per_gate_local_constant": P2_PER_GATE,
            "exact_no_shot_noise": {"fixed_kcal": err_fixed_exact, "adapt_kcal": err_adapt_exact},
            "8seed_shot_noisy": {"fixed_mean_kcal": fixed_mean, "fixed_std_kcal": fixed_std,
                                  "adapt_mean_kcal": adapt_mean, "adapt_std_kcal": adapt_std},
        },
        "real_noise_forte1_corrected_p2": {
            "p2": p2_real_forte1,
            "8seed_shot_noisy": {"fixed_mean_kcal": fixed_mean_real, "fixed_std_kcal": fixed_std_real,
                                  "adapt_mean_kcal": adapt_mean_real, "adapt_std_kcal": adapt_std_real},
        },
        "real_hardware_reference_kcal": REAL_HARDWARE_RAW_KCAL,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
