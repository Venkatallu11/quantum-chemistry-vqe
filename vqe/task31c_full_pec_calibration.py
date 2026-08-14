#!/usr/bin/env python3
"""
task31c_full_pec_calibration.py -- iteration 31, Task C. FULL PEC
calibration via literal quasi-probability twirling, all 21 H4
measurement slots (not 4 representative ones), N_MC=16, forte-1 only
(aria-1 retired). Uses Task 31A's corrected p2(zz)=0.0146 consensus.
============================================================================
EFFICIENCY, justified not assumed: for the `ideal` model, p=0 everywhere,
so `pec_inverse_weights(0, n)` puts ALL probability mass on "insert
identity" (verified: Task 30B's own log showed 0/126 components ever
selected a non-identity branch for ideal) -- meaning every one of the 16
"twirled" draws for ideal is IDENTICAL to the untwirled circuit. Submitting
4,368 real duplicate circuits for ideal would be pure waste; ideal's PEC
value is reused directly from Task 28B's existing real checkpoint
(PEC-corrected = raw exactly, by construction, when p=0).

Batched at <=91 circuits/job (273 provably hits a service boundary, per
this project's own outage history).

Run:
    python vqe/task31c_full_pec_calibration.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from loop_pec import pec_inverse_weights, gamma_factor
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from native_stateprep import to_native
from task28d_all_gate_zne import optimized_native_circuit
from qiskit.quantum_info import Pauli

K = 6
SHOTS = 100_000
N_SEEDS = 8
N_MC = 16
GATE_NAME = "zz"
P2_ZZ = 0.0146  # Task 31A's 9/11-position consensus
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task31c_full_pec_calibration.json")
RAW_CKPT = os.path.join(CKPT_DIR, "task28b_optimized_raw.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31c_full_pec_calibration_results.json")


def sample_twirled_circuit(base_qc, p2, p1_gpi, p1_gpi2, rng):
    qc = base_qc.copy_empty_like()
    total_sign = 1
    for instr in base_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name == "gpi":
            weights = pec_inverse_weights(p1_gpi, 1)
        elif op.name == "gpi2":
            weights = pec_inverse_weights(p1_gpi2, 1)
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
        n = len(chosen_label)
        if chosen_label != "I" * n:
            qc.append(Pauli(chosen_label).to_instruction(), qargs)
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign


def main():
    print("\n" + "=" * 96)
    print(f"  task31c_full_pec_calibration.py -- ALL 21 slots, N_MC={N_MC}, forte-1 only")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins = learned["forte-1"]["gpi"]
    p1_gpi = float(np.mean([v["p"] for v in gpi_bins.values()]))
    p1_gpi2 = p1_gpi  # disclosed fallback, unchanged from Task 30B/31A
    gamma_gate2 = gamma_factor(P2_ZZ, 2)
    gamma_gate1 = gamma_factor(p1_gpi, 1)
    print(f"  p2(zz)={P2_ZZ}, p1(gpi/gpi2 fallback)={p1_gpi:.6f}, gamma_gate2={gamma_gate2:.4f}, gamma_gate1={gamma_gate1:.6f}")

    if os.path.exists(CKPT_PATH):
        print("\n  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"\n  connected, backend={backend.name}")

        circuits, tags, gamma_per_circuit = [], [], []
        rng = np.random.default_rng(stable_seed("task31c"))
        t0 = time.time()
        for name in kept:
            base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
            n2q = sum(1 for instr in base.data if instr.operation.name == GATE_NAME)
            n1q = sum(1 for instr in base.data if instr.operation.name in ("gpi", "gpi2"))
            for group in groups:
                combined = effrag_mod.combined_basis_label(group)
                basis_qc = native_basis_change(combined, GATE_NAME)
                full_base = base.compose(basis_qc)
                # basis-rotation gates add a FEW more 1q gates -- recount on the actual composed circuit
                n2q_full = sum(1 for instr in full_base.data if instr.operation.name == GATE_NAME)
                n1q_full = sum(1 for instr in full_base.data if instr.operation.name in ("gpi", "gpi2"))
                gamma_this = (gamma_gate2 ** n2q_full) * (gamma_gate1 ** n1q_full)
                for draw in range(N_MC):
                    twirled, sign = sample_twirled_circuit(full_base, P2_ZZ, p1_gpi, p1_gpi2, rng)
                    twirled.measure_all()
                    circuits.append(twirled)
                    tags.append((name, list(group), draw, sign))
                    gamma_per_circuit.append(gamma_this)
        print(f"  built {len(circuits)} twirled circuits ({len(kept)} slots x {len(groups)} groups x {N_MC} draws), "
              f"{time.time()-t0:.1f}s")

        n_chunks = (len(circuits) + 90) // 91
        chunk_size = (len(circuits) + n_chunks - 1) // n_chunks
        chunks_c = [circuits[i:i + chunk_size] for i in range(0, len(circuits), chunk_size)]
        print(f"  submitting in {len(chunks_c)} batches of <= {chunk_size} circuits each...")

        all_counts = []
        t0 = time.time()
        for ci, chunk in enumerate(chunks_c):
            job = None
            for attempt in range(6):
                try:
                    job = submit_job(chunk, backend, "forte-1", shots=SHOTS)
                    break
                except Exception as e:
                    wait_s = min(30 * (2 ** attempt), 300)
                    print(f"    batch {ci+1}/{len(chunks_c)} submit failed (attempt {attempt+1}/6): {e} -- backing off {wait_s}s")
                    time.sleep(wait_s)
            if job is None:
                raise RuntimeError(f"submit_job exhausted retries for batch {ci+1}/{len(chunks_c)}")
            counts = get_counts_list(job)
            all_counts.extend(counts)
            print(f"    batch {ci+1}/{len(chunks_c)} done ({len(all_counts)}/{len(circuits)} total), "
                  f"{time.time()-t0:.1f}s elapsed")

        ck = {"tags": tags, "counts": all_counts, "gamma_per_circuit": gamma_per_circuit}
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    print(f"\n  -- combining (forte-1) --")
    tags = ck["tags"]
    counts_list = ck["counts"]
    gamma_per_circuit = ck["gamma_per_circuit"]
    by_name_label = {}
    for (name, group, draw, sign), counts, gamma in zip(tags, counts_list, gamma_per_circuit):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma, counts))

    pec_kept_forte = {name: {} for name in kept}
    for (name, l), entries in by_name_label.items():
        seed_means = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task31c_combine", name, l, seed))
            vals_signed = []
            for sign, gamma, counts in entries:
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m = expectation_from_counts(resampled, l)
                vals_signed.append(sign * gamma * m)
            seed_means.append(float(np.mean(vals_signed)))
        val = max(-1.0, min(1.0, float(np.mean(seed_means))))
        pec_kept_forte[name][l] = val

    def build_full(raw_kept):
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

    def energy_and_err(raw):
        alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
        E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                              exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
        return E, errs["err_vs_exact_kcal"]

    full_pec_forte = build_full(pec_kept_forte)
    _, err_pec_forte = energy_and_err(full_pec_forte)
    print(f"  forte-1 LITERAL PEC (full 21-slot twirling): {err_pec_forte:.3f} kcal/mol")

    # -- ideal: reuse existing real raw data, PEC=raw exactly (p=0) --
    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    ideal_tags = raw_ck["tags"]["ideal"]
    ideal_counts = raw_ck["counts"]["ideal"]
    ideal_per_name = {}
    for (name, group), counts in zip(ideal_tags, ideal_counts):
        ideal_per_name.setdefault(name, {}).setdefault(tuple(group), counts)
    ideal_kept = {name: {} for name in kept}
    for name in kept:
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task31c_ideal", name, seed))
            for group_t, counts in ideal_per_name[name].items():
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    ideal_kept[name].setdefault(l, []).append(expectation_from_counts(resampled, l))
    for name in kept:
        for l in ideal_kept[name]:
            ideal_kept[name][l] = float(np.mean(ideal_kept[name][l]))
    full_ideal = build_full(ideal_kept)
    _, err_ideal = energy_and_err(full_ideal)
    print(f"  ideal (reused real data, PEC=raw exactly since p=0): {err_ideal:.3f} kcal/mol")

    ideal_ok = err_ideal < 2.0
    print(f"\n  IDEAL-CONTROL CHECK: {'PASS' if ideal_ok else 'FAIL'}")

    # -- also raw forte-1 for reference (from the same existing checkpoint) --
    forte_tags = raw_ck["tags"]["forte-1"]
    forte_counts = raw_ck["counts"]["forte-1"]
    forte_per_name = {}
    for (name, group), counts in zip(forte_tags, forte_counts):
        forte_per_name.setdefault(name, {}).setdefault(tuple(group), counts)
    raw_forte_kept = {name: {} for name in kept}
    for name in kept:
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("task31c_rawforte", name, seed))
            for group_t, counts in forte_per_name[name].items():
                resampled = bootstrap_counts(counts, SHOTS, rng)
                for l in group_t:
                    raw_forte_kept[name].setdefault(l, []).append(expectation_from_counts(resampled, l))
    for name in kept:
        for l in raw_forte_kept[name]:
            raw_forte_kept[name][l] = float(np.mean(raw_forte_kept[name][l]))
    full_raw_forte = build_full(raw_forte_kept)
    _, err_raw_forte = energy_and_err(full_raw_forte)
    print(f"  forte-1 RAW (reference): {err_raw_forte:.3f} kcal/mol")

    print(f"\n  -- SUMMARY --")
    print(f"    forte-1: raw={err_raw_forte:.2f}  literal-PEC(full 21-slot, N_MC=16)={err_pec_forte:.2f} kcal/mol")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "err_ideal": err_ideal, "err_raw_forte": err_raw_forte, "err_pec_forte_literal": err_pec_forte,
            "ideal_control_pass": bool(ideal_ok), "pec_kept_forte": pec_kept_forte,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
