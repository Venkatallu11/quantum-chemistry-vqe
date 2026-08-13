#!/usr/bin/env python3
"""
task28a_1q_spectroscopy.py — iteration 28, Task A. Converts iteration
27's INFERENCE (native circuits carry far more 1-qubit gates than
abstract ones, a "strong candidate explanation" for the 9.3x ZNE floor)
into a MEASUREMENT: build circuits that are EXACTLY unitary-equivalent
(same state, same observable) but with N_1q artificially varied via
cancelling GPi pairs, and measure dE/dN_1q directly on real hardware
noise models.
============================================================================
BASELINE CHOICE, disclosed: the task's own framing ("same 11 ZZ gates...
padded... to give N_1q = 50, 100, 150, 200, 250, 300") assumes a natural
baseline near or below 50 -- but the UN-optimized native circuit's own
N_1q (197 ms / 285 zz, iteration 27 Task A) is already ABOVE the entire
50-300 target range, so "padding" from there cannot reach 50-250 at all
(padding only ADDS gates). This file instead uses Task 28B's own
TrappedIonOptimizerPlugin-reduced circuit (verified statevector-identical
to <1e-12, N_1q~47-78 depending on target/family) as the padding
baseline -- the lowest REAL, verified-correct circuit available this
session for this logical state -- and pads UPWARD from there. Only
target levels AT OR ABOVE the baseline are reachable and tested; this is
stated explicitly, not silently substituted.

PADDING MECHANISM, verified before use (not assumed): GPi(phi)^2 = I
EXACTLY for any phi (checked directly: 5 phi values, max deviation 0
to 1.1e-16, i.e. exact to machine precision) -- so inserting
GPi(0), GPi(0) pairs at the end of a circuit is a mathematically exact
identity, regardless of phase choice. Every padded circuit's statevector
is verified against the unpadded circuit's own statevector to <1e-12
BEFORE being queued for submission (a fold that fails this check is not
submitted, matching Task B's own fold-preservation discipline).

Uses 3 representative diagonal slots (u_0, u_1, u_2) x all 13 measurement
groups, at every REACHABLE padding level, submitted concurrently to
ideal/aria-1/forte-1 -- a genuine, real measurement, not a proxy.

Run:
    python vqe/task28a_1q_spectroscopy.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, native_target
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from qiskit_ionq.ionq_gates import GPIGate
from qiskit.quantum_info import Statevector

K = 6
REPRESENTATIVE_SLOTS = ["u_0", "u_1", "u_2"]
TARGET_N1Q = [50, 100, 150, 200, 250, 300]
SHOTS = 100_000
N_SEEDS = 8
GATE_BY_MODEL = {"ideal": "ms", "aria-1": "ms", "forte-1": "zz"}
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task28a_1q_spectroscopy.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task28a_1q_spectroscopy_results.json")


def verify_gpi_involution():
    worst = 0.0
    for phi in [0.0, 0.1, 0.25, 0.37, 0.5]:
        M = np.asarray(GPIGate(phi).to_matrix())
        worst = max(worst, float(np.max(np.abs(M @ M - np.eye(2)))))
    return worst


def optimized_native_circuit(angles, gate_name):
    from qiskit.transpiler import PassManagerConfig
    from qiskit_ionq import TrappedIonOptimizerPlugin
    qc = build_ansatz(angles)
    native = to_native(qc, gate_name)
    tgt = native_target(qc.num_qubits, gate_name)
    pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt), optimization_level=3)
    return pm.run(native)


def pad_to_n1q(base_circuit, target_n1q, gate_name):
    """Appends GPi(0),GPi(0) cancelling pairs (2 gates each) on qubit 0
    until N_1q reaches (or would first meet/exceed) target_n1q. Returns
    None if target_n1q is below the base circuit's own count (unreachable
    by padding)."""
    counts = base_circuit.count_ops()
    base_n1q = sum(v for k, v in counts.items() if k in ("gpi", "gpi2"))
    if target_n1q < base_n1q:
        return None, base_n1q
    n_pairs = (target_n1q - base_n1q + 1) // 2   # each pair adds 2 gates
    padded = base_circuit.copy()
    for _ in range(n_pairs):
        padded.append(GPIGate(0.0), [0])
        padded.append(GPIGate(0.0), [0])
    actual_n1q = base_n1q + 2 * n_pairs
    return padded, actual_n1q


def main():
    print("\n" + "=" * 96)
    print("  task28a_1q_spectroscopy.py -- one-qubit noise spectroscopy, converts inference to measurement")
    print("=" * 96)

    gpi_err = verify_gpi_involution()
    print(f"  GPi(phi)^2 = I verified: worst deviation {gpi_err:.2e} (machine precision)")
    assert gpi_err < 1e-12

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)

    if os.path.exists(CKPT_PATH):
        print("  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"  connected, backend={backend.name}")

        base_circuits = {}
        reachable_levels = {}
        for gate_name in ["ms", "zz"]:
            for name in REPRESENTATIVE_SLOTS:
                angles = fixed_solutions[name]["angles"]
                base = optimized_native_circuit(angles, gate_name)
                base_circuits[(name, gate_name)] = base
                base_n1q = sum(v for k, v in base.count_ops().items() if k in ("gpi", "gpi2"))
                reachable_levels[(name, gate_name)] = [t for t in TARGET_N1Q if t >= base_n1q]
                print(f"    {name}/{gate_name}: baseline N_1q={base_n1q}, reachable targets={reachable_levels[(name, gate_name)]}")

        # -- build + verify every padded circuit, per (slot, gate_name, target_n1q) --
        padded_circuits = {}
        max_verify_err = 0.0
        for (name, gate_name), base in base_circuits.items():
            sv_base = np.asarray(Statevector.from_instruction(base))
            for target in reachable_levels[(name, gate_name)]:
                padded, actual_n1q = pad_to_n1q(base, target, gate_name)
                sv_padded = np.asarray(Statevector.from_instruction(padded))
                idx = int(np.argmax(np.abs(sv_base)))
                phase = sv_padded[idx] / sv_base[idx] if abs(sv_base[idx]) > 1e-9 else 1.0
                err = float(np.max(np.abs(sv_padded / phase - sv_base)))
                max_verify_err = max(max_verify_err, err)
                padded_circuits[(name, gate_name, target)] = (padded, actual_n1q)
        print(f"  unitary-equivalence check, ALL padded circuits vs their own base: worst={max_verify_err:.2e} "
              f"(must be <1e-12): {'PASS' if max_verify_err < 1e-12 else 'FAIL -- refusing to submit'}")
        assert max_verify_err < 1e-12, "a padded circuit is NOT unitary-identical to its base -- refusing to submit"

        # -- build measurement circuits + submit, concurrent, non-blocking --
        all_jobs = {}
        t0 = time.time()
        for model, gate_name in GATE_BY_MODEL.items():
            circuits, tags = [], []
            for name in REPRESENTATIVE_SLOTS:
                for target in reachable_levels[(name, gate_name)]:
                    padded, actual_n1q = padded_circuits[(name, gate_name, target)]
                    for group in groups:
                        combined = effrag_mod.combined_basis_label(group)
                        basis_qc = native_basis_change(combined, gate_name)
                        qc = padded.compose(basis_qc)
                        qc.measure_all()
                        circuits.append(qc)
                        tags.append((name, target, actual_n1q, tuple(group)))
            job = submit_job(circuits, backend, model, shots=SHOTS)
            all_jobs[model] = (job, tags)
        t_submit = time.time() - t0
        print(f"\n  all {len(all_jobs)} model jobs submitted (non-blocking), {t_submit:.1f}s")

        t0 = time.time()
        all_counts = {}
        for i, (model, (job, tags)) in enumerate(all_jobs.items()):
            all_counts[model] = get_counts_list(job)
            print(f"    retrieved {i+1}/{len(all_jobs)}: model={model}, {time.time()-t0:.1f}s elapsed")
        t_retrieve = time.time() - t0

        ck = {
            "tags": {m: [[t[0], t[1], t[2], list(t[3])] for t in v[1]] for m, v in all_jobs.items()},
            "counts": all_counts, "reachable_levels": {f"{k[0]}|{k[1]}": v for k, v in reachable_levels.items()},
            "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        }
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- exact ideal expectations per (slot, label) for the mean-abs-delta metric --
    exact_ideal = {}
    for name in REPRESENTATIVE_SLOTS:
        sv = Statevector.from_instruction(build_ansatz(fixed_solutions[name]["angles"]))
        from qiskit.quantum_info import Pauli
        exact_ideal[name] = {l: float(sv.expectation_value(Pauli(l)).real) for l in non_id_labels}

    # -- assemble mean|delta| (measured - ideal) vs N_1q, per model --
    print(f"\n  -- mean|delta| (measured - exact ideal) vs actual N_1q, per model --")
    dose_response = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        by_n1q = {}
        per_tag = {}
        for (name, target, actual_n1q, group), counts in zip(tags, counts_list):
            per_tag.setdefault((name, actual_n1q), {}).setdefault(tuple(group), counts)
        for (name, actual_n1q), group_counts in per_tag.items():
            deltas = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task28a", model, name, actual_n1q, seed))
                for group_t, counts in group_counts.items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    for l in group_t:
                        m = expectation_from_counts(resampled, l)
                        deltas.append(abs(m - exact_ideal[name][l]))
            by_n1q.setdefault(actual_n1q, []).extend(deltas)
        curve = {n1q: {"mean": float(np.mean(v)), "std": float(np.std(v)) / np.sqrt(len(v))} for n1q, v in sorted(by_n1q.items())}
        dose_response[model] = curve
        print(f"    {model}: " + "  ".join(f"N_1q={n}: {c['mean']:.4f}+/-{c['std']:.4f}" for n, c in curve.items()))

    # -- linear fit: dE/dN_1q per model --
    print(f"\n  -- linear regression: d(mean|delta|)/dN_1q --")
    slopes = {}
    for model in ["aria-1", "forte-1"]:
        curve = dose_response[model]
        x = np.array(sorted(curve.keys()))
        y = np.array([curve[n]["mean"] for n in x])
        if len(x) >= 2:
            slope, intercept = np.polyfit(x, y, 1)
        else:
            slope, intercept = float("nan"), float("nan")
        slopes[model] = {"slope": float(slope), "intercept": float(intercept), "n_points": len(x)}
        print(f"    {model}: slope={slope:.6f} per 1q gate, intercept={intercept:.4f}, n_points={len(x)}")
        small_slope = abs(slope) < 1e-4
        print(f"    {'SMALL SLOPE -- the 1q hypothesis may be WRONG, pivot needed' if small_slope else 'non-trivial slope -- 1q hypothesis supported, proceed to 28C/28D'}")

    results = {"gpi_involution_check": gpi_err, "dose_response": dose_response, "slopes": slopes,
               "reachable_levels": ck.get("reachable_levels", {})}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
