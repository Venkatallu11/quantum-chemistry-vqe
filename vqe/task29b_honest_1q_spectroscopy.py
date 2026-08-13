#!/usr/bin/env python3
"""
task29b_honest_1q_spectroscopy.py -- iteration 29, Task B. REDO 1q
SPECTROSCOPY HONESTLY. Task 28A (see RESEARCH_LEDGER.md correction) only
tested the OPTIMIZER-REDUCED circuit, padded with cancelling IDENTICAL-
ANGLE GPi(0)-GPi(0) pairs -- the most benign possible 1q gate, nothing
like the ~234 different working gates a real synthesis compiler emits on
the ORIGINAL 285-1q-gate topology. Its slope was also reported in the
wrong units (Pauli-expectation error, not kcal/mol).
============================================================================
THREE FIXES:
  1. TWO circuit arms: the UN-OPTIMIZED native circuit (285 1q / 11 2q,
     forte's zz family -- the circuit that actually produced the 89
     kcal/mol native-vs-abstract gap) alongside the OPTIMIZED one (Task
     28B's circuit, ~61-78 1q mean).
  2. TWO padding styles, compared directly: (i) CANCELLING IDENTICAL-ANGLE
     GPi(0)-GPi(0) pairs (28A's original method, angle=0 always -- the
     most benign possible 1q gate). (ii) VARIED-ANGLE synthesis-style
     pairs: GPi2(phi_k), GPi2(phi_k+0.5) for a DIFFERENT phi_k per
     inserted pair (verified exact identity to <1e-16 for several angles
     before use, same standing/inverse fact `fold_all_gates` already
     relies on) -- a materially closer proxy for what a real synthesis
     pass emits (many distinct rotation angles) than repeating one gate.
  3. Slopes reported in kcal/mol of ENERGY, with the Hamiltonian-
     coefficient-weighted conversion shown explicitly, not left in raw
     Pauli-expectation units and mislabeled.

If the two padding styles give DIFFERENT slopes, that difference IS a
real finding (Task 29's own framing) -- reported as such, not smoothed
over.

Run:
    python vqe/task29b_honest_1q_spectroscopy.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, native_target
from task2_fold_response_dataset import native_basis_change
from task27c_full_h4_folds import kept_slots_for_K
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from task28d_all_gate_zne import optimized_native_circuit
from qiskit_ionq.ionq_gates import GPIGate, GPI2Gate
from qiskit.quantum_info import Statevector, Pauli

K = 6
SHOTS = 100_000
N_SEEDS = 8
GATE_NAME = "zz"          # forte-1's family -- the one tied to the 89 kcal/mol native-vs-abstract gap
MODELS = ["ideal", "aria-1", "forte-1"]
REPRESENTATIVE_SLOTS = ["u_0", "u_1"]
N_PAD_LEVELS = [50, 100]   # how many EXTRA 1q gates to pad in, on top of each arm's own baseline
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task29b_honest_1q_spectroscopy.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task29b_honest_1q_spectroscopy_results.json")


def n1q_count(qc):
    return sum(1 for instr in qc.data if instr.operation.name in ("gpi", "gpi2"))


def pad_identical(qc, n_pairs):
    """Style (i): n_pairs of GPi(0), GPi(0) -- exact self-inverse, angle=0 always (28A's original method)."""
    padded = qc.copy()
    q0 = qc.qubits[0]
    for _ in range(n_pairs):
        padded.append(GPIGate(0), [q0])
        padded.append(GPIGate(0), [q0])
    return padded


def pad_varied(qc, n_pairs, rng):
    """Style (ii): n_pairs of GPi2(phi_k), GPi2(phi_k+0.5) -- exact identity PER PAIR (verified
    separately before use), but a DIFFERENT phi_k drawn per pair -- varied angles throughout,
    much closer to what a real synthesis compiler actually emits than repeating one gate."""
    padded = qc.copy()
    q0 = qc.qubits[0]
    for _ in range(n_pairs):
        phi = float(rng.uniform(0, 1))
        padded.append(GPI2Gate(phi), [q0])
        padded.append(GPI2Gate(phi + 0.5), [q0])
    return padded


def verify_exact(base_qc, padded_qc, label):
    sv_base = np.asarray(Statevector.from_instruction(base_qc))
    sv_pad = np.asarray(Statevector.from_instruction(padded_qc))
    idx = int(np.argmax(np.abs(sv_base)))
    phase = sv_pad[idx] / sv_base[idx] if abs(sv_base[idx]) > 1e-9 else 1.0
    err = float(np.max(np.abs(sv_pad / phase - sv_base)))
    assert err < 1e-10, f"{label}: padded circuit not unitary-equivalent to base, err={err:.3e}"
    return err


def main():
    print("\n" + "=" * 96)
    print("  task29b_honest_1q_spectroscopy.py -- two circuit arms x two padding styles, kcal/mol units")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)

    # -- energy sensitivity weight per label: |coefficient| summed over Hamiltonian terms
    #    touching that alpha label -- this is the explicit conversion from Pauli-expectation
    #    error to an energy-error contribution that 28A never applied --
    label_weight = {}
    for (a_label, b_label, coeff) in p["terms"]:
        w = abs(coeff) * HARTREE_TO_KCAL_MOL
        label_weight[a_label] = label_weight.get(a_label, 0.0) + w
    print(f"  built Hamiltonian-coefficient energy weight for {len(label_weight)} alpha labels "
          f"(sum={sum(label_weight.values()):.4f} kcal/mol total 'budget' if every label were maximally wrong)")

    # -- two circuit arms --
    arms = {}
    for name in REPRESENTATIVE_SLOTS:
        angles = fixed_solutions[name]["angles"]
        unopt = to_native(build_ansatz(angles), GATE_NAME)
        opt = optimized_native_circuit(angles, GATE_NAME)
        arms.setdefault("unoptimized", {})[name] = unopt
        arms.setdefault("optimized", {})[name] = opt
        print(f"  {name}: unoptimized N_1q={n1q_count(unopt)}  optimized N_1q={n1q_count(opt)}")

    # -- build every (arm, name, style, level) padded circuit, verify each exactly before use --
    variants = {}   # (arm, name, style, n1q_target) -> circuit
    rng_pad = np.random.default_rng(20290)
    for arm, circuits_by_name in arms.items():
        for name, base_qc in circuits_by_name.items():
            base_n1q = n1q_count(base_qc)
            variants[(arm, name, "baseline", base_n1q)] = base_qc
            for extra in N_PAD_LEVELS:
                target_n1q = base_n1q + extra
                n_pairs = extra // 2
                v_id = pad_identical(base_qc, n_pairs)
                verify_exact(base_qc, v_id, f"{arm}/{name}/identical/+{extra}")
                variants[(arm, name, "identical", target_n1q)] = v_id
                v_var = pad_varied(base_qc, n_pairs, rng_pad)
                verify_exact(base_qc, v_var, f"{arm}/{name}/varied/+{extra}")
                variants[(arm, name, "varied", target_n1q)] = v_var
    print(f"  built + verified {len(variants)} circuit variants (all unitary-equivalent to their own base, <1e-10)")

    exact_ideal = {}
    for name in REPRESENTATIVE_SLOTS:
        sv = Statevector.from_instruction(build_ansatz(fixed_solutions[name]["angles"]))
        exact_ideal[name] = {l: float(sv.expectation_value(Pauli(l)).real) for l in non_id_labels}

    if os.path.exists(CKPT_PATH):
        print("  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"  connected, backend={backend.name}")

        all_jobs = {}
        t0 = time.time()
        for model in MODELS:
            circuits, tags = [], []
            for (arm, name, style, target_n1q), qc in variants.items():
                for group in groups:
                    combined = effrag_mod.combined_basis_label(group)
                    basis_qc = native_basis_change(combined, GATE_NAME)
                    full_qc = qc.compose(basis_qc)
                    full_qc.measure_all()
                    circuits.append(full_qc)
                    tags.append((arm, name, style, target_n1q, list(group)))
            job = submit_job(circuits, backend, model, shots=SHOTS)
            all_jobs[model] = (job, tags)
            print(f"    submitted model={model}: {len(circuits)} circuits, job_id={job.job_id()}")
        t_submit = time.time() - t0

        t0 = time.time()
        all_counts = {}
        for model, (job, tags) in all_jobs.items():
            all_counts[model] = get_counts_list(job)
            print(f"    retrieved model={model}, {time.time()-t0:.1f}s elapsed")
        t_retrieve = time.time() - t0

        ck = {
            "tags": {m: [[t[0], t[1], t[2], t[3], t[4]] for t in v[1]] for m, v in all_jobs.items()},
            "counts": all_counts, "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        }
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- mean|delta| (Pauli-expectation units) AND energy-weighted mean|delta| (kcal/mol) vs N_1q --
    print(f"\n  -- dose-response: Pauli-expectation error AND energy-weighted (kcal/mol) error vs N_1q --")
    dose_response = {}
    for model in MODELS:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        per_key = {}
        for (arm, name, style, target_n1q, group), counts in zip(tags, counts_list):
            per_key.setdefault((arm, name, style, target_n1q), {}).setdefault(tuple(group), counts)

        by_arm_style_n1q_pauli = {}
        by_arm_style_n1q_energy = {}
        for (arm, name, style, target_n1q), group_counts in per_key.items():
            pauli_deltas = []
            energy_deltas = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task29b", model, arm, name, style, target_n1q, seed))
                for group_t, counts in group_counts.items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    for l in group_t:
                        m = expectation_from_counts(resampled, l)
                        d = abs(m - exact_ideal[name][l])
                        pauli_deltas.append(d)
                        energy_deltas.append(d * label_weight.get(l, 0.0))
            key = (arm, style, target_n1q)
            by_arm_style_n1q_pauli.setdefault(key, []).extend(pauli_deltas)
            by_arm_style_n1q_energy.setdefault(key, []).extend(energy_deltas)

        curve = {}
        for key in by_arm_style_n1q_pauli:
            curve[key] = {
                "pauli_mean": float(np.mean(by_arm_style_n1q_pauli[key])),
                "energy_kcal_mean": float(np.mean(by_arm_style_n1q_energy[key])),
            }
        dose_response[model] = curve
        print(f"\n    {model}:")
        for (arm, style, n1q), v in sorted(curve.items()):
            print(f"      arm={arm:<11} style={style:<9} N_1q={n1q:<4} "
                  f"pauli_err={v['pauli_mean']:.5f}  energy_err={v['energy_kcal_mean']:.5f} kcal/mol")

    # -- slopes, kcal/mol per 1q gate, per (model, arm, style) --
    print(f"\n  -- slopes: d(energy_kcal)/d(N_1q), per model x arm x style --")
    slopes = {}
    for model in MODELS:
        slopes[model] = {}
        for arm in ["unoptimized", "optimized"]:
            for style in ["identical", "varied"]:
                pts = [(n1q, v["energy_kcal_mean"]) for (a, s, n1q), v in dose_response[model].items()
                       if a == arm and s == style]
                # include the shared baseline point (style-independent) for this arm
                base_pts = [(n1q, v["energy_kcal_mean"]) for (a, s, n1q), v in dose_response[model].items()
                            if a == arm and s == "baseline"]
                pts = sorted(set(pts + base_pts))
                if len(pts) < 2:
                    continue
                x = np.array([q[0] for q in pts], dtype=float)
                y = np.array([q[1] for q in pts], dtype=float)
                slope, intercept = np.polyfit(x, y, 1)
                slopes[model][f"{arm}_{style}"] = {"slope_kcal_per_1q": float(slope),
                                                    "intercept_kcal": float(intercept), "n_points": len(pts)}
                print(f"    {model:<8} arm={arm:<11} style={style:<9}: "
                      f"slope={slope:.6f} kcal/mol per 1q gate  (n={len(pts)} points)")

    print(f"\n  -- DOES PADDING STYLE MATTER? identical vs varied, same arm/model --")
    style_diff_flags = {}
    for model in MODELS:
        for arm in ["unoptimized", "optimized"]:
            s_id = slopes[model].get(f"{arm}_identical")
            s_var = slopes[model].get(f"{arm}_varied")
            if s_id and s_var:
                ratio = (s_var["slope_kcal_per_1q"] / s_id["slope_kcal_per_1q"]
                         if abs(s_id["slope_kcal_per_1q"]) > 1e-9 else float("inf"))
                diff_kcal = abs(s_var["slope_kcal_per_1q"] - s_id["slope_kcal_per_1q"])
                print(f"    {model:<8} arm={arm:<11}: identical={s_id['slope_kcal_per_1q']:.6f}  "
                      f"varied={s_var['slope_kcal_per_1q']:.6f}  ratio={ratio:.2f}x  diff={diff_kcal:.6f} kcal/mol")
                style_diff_flags[f"{model}_{arm}"] = {"ratio": ratio, "diff_kcal": diff_kcal}

    results = {
        "gate_name": GATE_NAME, "n_pad_levels": N_PAD_LEVELS, "representative_slots": REPRESENTATIVE_SLOTS,
        "label_weight_kcal": label_weight, "dose_response": {m: {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in c.items()}
                                                               for m, c in dose_response.items()},
        "slopes": slopes, "style_diff": style_diff_flags,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
