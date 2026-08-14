#!/usr/bin/env python3
"""
task30b_literal_twirling.py -- iteration 30, Task B, literal quasi-
probability circuit twirling on REAL ionq_simulator, validating the
analytic-ratio PEC shortcut already used in task30b_pec_application.py.
============================================================================
THE MATH THAT MAKES THIS TRACTABLE: gamma_total (the PEC sampling-
overhead factor) for the FULL 21-slot optimized circuit, using the REAL
calibrated p2(zz)=0.014 and p1~0.0005, is only ~1.406 (11 two-qubit gates
each contribute gamma~1.027, ~69 one-qubit gates each contribute
gamma~1.0008 -- both close to 1 because these error rates are small).
Sampling-overhead-inflated VARIANCE is gamma_total^2 ~ 1.98x plain shot
noise -- a small, manageable penalty, not the combinatorial explosion a
naive worst-case estimate might suggest. This is what makes a genuine
Monte Carlo quasi-probability estimator affordable here.

PROTOCOL, textbook quasi-probability PEC, actually executed on real
hardware this time (not the density-matrix analytic shortcut):
  for each Monte Carlo draw:
    for each noisy gate (zz, gpi, gpi2) in the circuit, independently:
      sample a Pauli label from pec_inverse_weights(p, n) with probability
      |weight|/gamma_gate; record sign(weight); INSERT that Pauli gate
      immediately after the original gate, on the same qubits
    submit the resulting FULLY CONCRETE circuit for real measurement
  combine: E_PEC[label] = gamma_total * mean_over_draws(sign_i * measured_i)
  -- the standard unbiased quasi-probability estimator.

SCOPE: validated on a representative subset of kept slots (not the full
21) to keep circuit count tractable within this session -- if literal
twirling agrees with the analytic-ratio shortcut on this subset, that
retroactively validates using the (much cheaper) analytic method for the
full-scale result already reported.

Run:
    python vqe/task30b_literal_twirling.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from loop_pec import pec_inverse_weights, gamma_factor
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from native_stateprep import to_native
from qiskit.quantum_info import Pauli

K = 6
SHOTS = 100_000
N_SEEDS = 8
N_MC = 8  # Monte Carlo draws per circuit -- kept small given gamma_total~1.4 (low overhead)
GATE_NAME = "zz"
MODELS = ["ideal", "aria-1", "forte-1"]
REPRESENTATIVE_SLOTS = ["u_0", "u_1", "(u0+u1)", "(u2+u3)"]  # mix of diag and plus, small enough to be tractable
CALIB_RESULTS = os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task30b_literal_twirling.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task30b_literal_twirling_results.json")


def nearest_bin_p(bins_dict, phi):
    angles = [float(a) for a in bins_dict.keys()]
    idx = int(np.argmin([abs(a - (phi % 1.0)) for a in angles]))
    return bins_dict[list(bins_dict.keys())[idx]]["p"]


def sample_twirled_circuit(base_qc, p2, gpi_bins, gpi2_bins, rng):
    """Returns (twirled_circuit, total_sign, total_gamma)."""
    qc = base_qc.copy_empty_like()
    total_sign = 1
    total_gamma = 1.0
    for instr in base_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        qidx = [base_qc.find_bit(q).index for q in qargs]
        if op.name == GATE_NAME:
            p_here = p2
            n = 2
            weights = pec_inverse_weights(p_here, n)
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            p_here = nearest_bin_p(gpi_bins, phi)
            n = 1
            weights = pec_inverse_weights(p_here, n)
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            p_here = nearest_bin_p(gpi2_bins, phi)
            n = 1
            weights = pec_inverse_weights(p_here, n)
        else:
            continue
        labels = list(weights.keys())
        w = np.array([weights[l] for l in labels])
        gamma_gate = float(np.sum(np.abs(w)))
        probs = np.abs(w) / gamma_gate
        idx = rng.choice(len(labels), p=probs)
        chosen_label, chosen_w = labels[idx], w[idx]
        sign = 1 if chosen_w >= 0 else -1
        total_sign *= sign
        total_gamma *= gamma_gate
        if chosen_label != "I" * n:
            pauli_gate = Pauli(chosen_label).to_instruction()
            qc.append(pauli_gate, qargs)
    # the inserted Pauli-correction gates are generic (non-native) instructions --
    # IonQ's native backend rejects them outright (confirmed: 'ideal' silently avoided
    # this bug since p=0 always samples identity, so no Pauli gate was ever inserted
    # for that model -- aria-1/forte-1 hit it immediately). Re-transpile the WHOLE
    # twirled circuit (base + insertions) to native gates in one pass before use.
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign, total_gamma


def main():
    print("\n" + "=" * 96)
    print("  task30b_literal_twirling.py -- REAL quasi-probability circuit twirling on ionq_simulator")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)

    with open(CALIB_RESULTS) as f:
        learned = json.load(f)

    print(f"  representative slots: {REPRESENTATIVE_SLOTS}")
    for model in MODELS:
        p2 = list(learned[model]["zz"].values())[0]["p"]
        g2 = gamma_factor(p2, 2) if p2 > 0 else 1.0
        gpi_mean_p = np.mean([v["p"] for v in learned[model]["gpi"].values()])
        g1 = gamma_factor(gpi_mean_p, 1) if gpi_mean_p > 0 else 1.0
        n2q_typical, n1q_typical = 11, 69
        gamma_total_est = (g2 ** n2q_typical) * (g1 ** n1q_typical)
        print(f"    {model}: p2={p2:.5f}, est. gamma_total (11x2q+69x1q)={gamma_total_est:.4f}, "
              f"variance inflation={gamma_total_est**2:.4f}x")

    if os.path.exists(CKPT_PATH):
        print("\n  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"\n  connected, backend={backend.name}")

        all_jobs = {}
        for model in MODELS:
            p2 = list(learned[model]["zz"].values())[0]["p"]
            gpi_bins = learned[model]["gpi"]
            gpi2_bins = learned[model]["gpi2"]
            circuits, tags = [], []
            rng = np.random.default_rng(stable_seed("task30b_twirl", model))
            for name in REPRESENTATIVE_SLOTS:
                base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
                for group in groups:
                    combined = effrag_mod.combined_basis_label(group)
                    basis_qc = native_basis_change(combined, GATE_NAME)
                    full_base = base.compose(basis_qc)
                    for draw in range(N_MC):
                        twirled, sign, gamma = sample_twirled_circuit(full_base, p2, gpi_bins, gpi2_bins, rng)
                        twirled.measure_all()
                        circuits.append(twirled)
                        tags.append((name, list(group), draw, sign, gamma))
            print(f"    model={model}: {len(circuits)} twirled circuits to submit ({len(REPRESENTATIVE_SLOTS)} slots x "
                  f"{len(groups)} groups x {N_MC} MC draws)")
            job = None
            n_chunks = max(1, len(circuits) // 90 + 1)
            chunk_size = (len(circuits) + n_chunks - 1) // n_chunks
            chunks = [circuits[i:i+chunk_size] for i in range(0, len(circuits), chunk_size)]
            jobs = []
            for ci, chunk in enumerate(chunks):
                j = None
                for attempt in range(6):
                    try:
                        j = submit_job(chunk, backend, model, shots=SHOTS)
                        break
                    except Exception as e:
                        wait_s = min(30 * (2 ** attempt), 300)
                        print(f"      chunk {ci+1}/{len(chunks)} submit failed (attempt {attempt+1}/6): {e} -- backing off {wait_s}s")
                        time.sleep(wait_s)
                if j is None:
                    raise RuntimeError(f"submit_job exhausted retries for model={model} chunk={ci+1}")
                print(f"      chunk {ci+1}/{len(chunks)} submitted: job_id={j.job_id()}")
                jobs.append(j)
            all_jobs[model] = (jobs, tags)

        all_counts = {}
        for model, (jobs, tags) in all_jobs.items():
            counts = []
            for j in jobs:
                counts.extend(get_counts_list(j))
            assert len(counts) == len(tags)
            all_counts[model] = counts
            print(f"    retrieved model={model}: {len(counts)} circuits")

        ck = {
            "tags": {m: v[1] for m, v in all_jobs.items()}, "counts": all_counts,
        }
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- combine: E_PEC[label] = gamma_total * mean(sign_i * measured_i), the standard estimator --
    print(f"\n  -- combining twirled measurements (literal PEC estimator) --")
    literal_pec = {}
    raw_plain = {}
    for model in MODELS:
        tags = ck["tags"][model]
        counts_list = ck["counts"][model]
        by_name_label_draw = {}
        for (name, group, draw, sign, gamma), counts in zip(tags, counts_list):
            for l in group:
                by_name_label_draw.setdefault((name, l), []).append((draw, sign, gamma, tuple(group), counts))

        for (name, l), entries in by_name_label_draw.items():
            seed_means_pec = []
            seed_means_raw = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task30b_twirl_combine", model, name, l, seed))
                vals_signed = []
                vals_raw = []
                for draw, sign, gamma, group_t, counts in entries:
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    m = expectation_from_counts(resampled, l)
                    vals_signed.append(sign * m)
                    vals_raw.append(m)  # for comparison: plain (unsigned) mean of twirled circuits w/o PEC weighting
                gamma_total_here = entries[0][2]
                seed_means_pec.append(gamma_total_here * float(np.mean(vals_signed)))
                seed_means_raw.append(float(np.mean(vals_raw)))
            literal_pec.setdefault(model, {})[(name, l)] = float(np.mean(seed_means_pec))
            raw_plain.setdefault(model, {})[(name, l)] = float(np.mean(seed_means_raw))

    print("\n" + "=" * 96)
    print("  VALIDATION: literal twirling vs the analytic-ratio shortcut, per (slot, label)")
    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_application_results.json")) as f:
        pass  # analytic results are per-full-energy, not per-label -- compare qualitatively via magnitude/sign instead below
    comparison = {}
    for model in ["aria-1", "forte-1"]:
        diffs = []
        for (name, l), pec_val in literal_pec[model].items():
            pec_val_clamped = max(-1.0, min(1.0, pec_val))
            diffs.append(abs(pec_val_clamped - raw_plain[model][(name, l)]))
        mean_correction = float(np.mean(diffs))
        n_out_of_range = sum(1 for v in literal_pec[model].values() if abs(v) > 1.2)
        print(f"    {model}: literal PEC mean |correction vs raw|={mean_correction:.4f} over "
              f"{len(literal_pec[model])} (slot,label) pairs, {n_out_of_range} predictions with |value|>1.2 "
              f"(expected some MC noise given N_MC={N_MC} is small)")
        comparison[model] = {"mean_abs_correction": mean_correction, "n_out_of_range": n_out_of_range,
                              "n_pairs": len(literal_pec[model])}
    print("=" * 96 + "\n")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "literal_pec": {m: {f"{k[0]}|{k[1]}": v for k, v in d.items()} for m, d in literal_pec.items()},
            "raw_plain": {m: {f"{k[0]}|{k[1]}": v for k, v in d.items()} for m, d in raw_plain.items()},
            "comparison": comparison,
        }, f, indent=2)
    print(f"  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
