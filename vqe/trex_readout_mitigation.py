#!/usr/bin/env python3
"""
trex_readout_mitigation.py — T-REx (Twirled Readout Error Extinction, van
den Berg, Minev, Kandala, Temme, PRA 105, 032620 (2022)), combined with
iteration 18's verified-real particle-number leakage post-selection.
============================================================================
WHY THIS FILE EXISTS: researched IonQ's own published error-mitigation
literature (per direct user request) rather than continuing to guess at
techniques. Found: (1) IonQ's own built-in "debiasing" feature is exactly
this project's missing piece for ZNE's persistent non-convergence
(coherent-noise-induced), but is REAL-QPU-ONLY per IonQ's own docs
(docs.ionq.com/guides/error-mitigation-debiasing: "not currently
available for the IonQ Quantum Cloud simulator, including the simulator
with noise model") -- the user was asked directly and declined to spend
real QPU credits to test it, so it is NOT attempted here, only recorded
in the ledger as a real, documented, structurally-untestable-under-this-
project's-constraints finding. (2) T-REx, this file's actual subject, is
implementable entirely with circuits WE build and submit to the free
ionq_simulator -- no server-side feature required, just twirl gates we
add ourselves.

MECHANISM: readout error can be ASYMMETRIC (P(0->1) != P(1->0) per
qubit). Twirling -- randomly applying X to a qubit right before
measurement, then classically flipping the observed bit back -- forces
the EFFECTIVE noise channel (after averaging over twirls) to be
SYMMETRIC. Under a symmetric per-qubit bit-flip channel, a weight-k
Z-type Pauli operator's measured expectation value is EXACTLY
<Z>_measured = (product of per-qubit lambda_i for the k active qubits)
* <Z>_ideal -- a single multiplicative damping factor, cheaply
calibrated once (run the bare |0> state with twirls; the ideal value is
exactly +1, so the measured value directly IS lambda_i) and divided back
out. This corrects READOUT-stage error specifically -- a different,
complementary mechanism from iteration 18's leakage post-selection
(which catches errors that corrupt the state's Hamming weight DURING
the circuit, propagated through to a wrong measurement; T-REx catches
errors introduced AT the measurement step itself, after the state is
already whatever it is).

VERIFIED LOCALLY FIRST, with a genuine READOUT error channel (this
project's local noise models have, until now, only ever modeled GATE
error -- qiskit_aer's ReadoutError is used here for the first time) with
a DELIBERATELY ASYMMETRIC bit-flip probability per qubit, confirming
T-REx recovers the correct expectation value that a naive (untwirled)
measurement gets wrong, before spending anything on a real submission.

Run:
    python vqe/trex_readout_mitigation.py --verify   (local, no network)
    python vqe/trex_readout_mitigation.py --targets   (real IonQ)
    python vqe/trex_readout_mitigation.py --assemble
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from fixed_ansatz import build_ansatz
from spin_leakage_postselect_ionq import with_ancilla_parity
import rank6_symmetry_vd as r
import ef_fragment as effrag
from qforge import combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from ionq_backend import connect_provider, get_simulator
from ionq_run import basis_change, IONQ_QIS_STANDARD_BASIS
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, NOISE_MODELS, submit_job, get_counts_list, bootstrap_counts,
    stable_seed, save_ckpt, load_ckpt,
)

from qiskit.circuit import QuantumCircuit
from qiskit.transpiler import CouplingMap
from qiskit import transpile

K = 6
N_TWIRLS = 8  # random X-twirl patterns per (target, group); van den Berg et al use similar small counts
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "trex_readout_mitigation_results.json")
CKPT_NAME = "trex_readout_targets"
CMAP5 = CouplingMap.from_full(5)


def hadamard_balanced_masks(n_qubits, n_twirls=8):
    """Exact per-qubit-balanced twirl design (not random sampling --
    verify.py's own local test caught a real bias bug from random
    sampling at small N_TWIRLS, see module docstring). Order-8 Hadamard
    matrix: rows 1-7 (excluding the all-+1 row) are each split exactly
    4-vs-4 across the 8 columns -- use n_qubits of those 7 balanced rows
    as the twirl assignment per qubit, columns as the n_twirls mask
    instances. Supports up to 7 qubits at n_twirls=8."""
    from scipy.linalg import hadamard
    assert n_qubits <= 7 and n_twirls == 8, "Hadamard(8)-based design supports <=7 qubits at exactly 8 twirls"
    H = hadamard(8)
    rows01 = ((H[1:n_qubits + 1] + 1) // 2).astype(int)  # (n_qubits, 8), each row balanced 4/4
    masks = [tuple(int(rows01[q, i]) for q in range(n_qubits)) for i in range(8)]
    for q in range(n_qubits):
        assert sum(m[q] for m in masks) == 4, f"qubit {q} twirl assignment not balanced"
    return masks


def twirled_measurement_circuit(base5, group_combined_label, twirl_mask, n_qubits=5):
    """twirl_mask: tuple of 0/1 per qubit (which qubits get an X twirl
    right before measurement, applied AFTER the group's basis-change so
    the twirl only affects the READOUT stage, not the measured basis
    itself)."""
    qc = base5.copy()
    basis_change(qc, group_combined_label)
    for q in range(n_qubits):
        if twirl_mask[q]:
            qc.x(q)
    qc.measure_all()
    return qc


def untwirl_counts(counts, twirl_mask, n_qubits=5):
    """Flip back the bits that were twirled, restoring what the
    measurement WOULD have read without the twirl (given the state at
    that point) -- per-shot, on the bitstring keys."""
    out = {}
    for bs, c in counts.items():
        # qiskit bitstring: bs[-(q+1)] is qubit q
        chars = list(bs)
        for q in range(n_qubits):
            if twirl_mask[q]:
                idx = len(chars) - 1 - q
                chars[idx] = '1' if chars[idx] == '0' else '0'
        out["".join(chars)] = out.get("".join(chars), 0) + c
    return out


def merge_counts(list_of_counts):
    out = {}
    for counts in list_of_counts:
        for bs, c in counts.items():
            out[bs] = out.get(bs, 0) + c
    return out


# ---------------------------------------------------------------------------
# LOCAL VERIFICATION: genuine (asymmetric) readout error channel, first use
# in this project of qiskit_aer's ReadoutError (all prior local models here
# only ever simulated GATE error).
# ---------------------------------------------------------------------------

def local_verify():
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error
    from qiskit.quantum_info import Statevector, Pauli

    print("\n" + "=" * 96)
    print("  trex_readout_mitigation.py --verify (local, asymmetric readout channel, no network)")
    print("=" * 96)

    n_qubits = 3
    rng = np.random.default_rng(0)
    # deliberately ASYMMETRIC per-qubit readout error: p(0->1) != p(1->0)
    p01 = [0.05, 0.08, 0.02]
    p10 = [0.01, 0.015, 0.06]
    nm = NoiseModel()
    for q in range(n_qubits):
        nm.add_readout_error(ReadoutError([[1 - p01[q], p01[q]], [p10[q], 1 - p10[q]]]), [q])
    nm.add_all_qubit_quantum_error(depolarizing_error(0.005, 2), "cx")

    qc = QuantumCircuit(n_qubits)
    qc.ry(2 * 0.9, 0)
    qc.cx(0, 1)
    qc.cx(0, 2)
    # BUG CAUGHT: a 3-qubit GHZ (H+CX+CX) has EXACT <ZZZ>=0 (the |111> branch
    # contributes -1, cancelling |000>'s +1) -- an earlier version of this
    # test hardcoded "exact <ZZZ>=1" from an incorrect assumption instead of
    # computing it, making both the naive and T-REx numbers look artificially
    # "close" to a wrong target. Fixed by (a) using an asymmetric-population
    # state (mostly |000>, small |011>/|101>-type admixture via reduced-angle
    # RY, so <ZZZ> is genuinely close to +1, not exactly cancelling) and (b)
    # computing the exact reference via Statevector directly -- never
    # hand-assumed again.
    exact_val = float(np.real(Statevector.from_instruction(qc).expectation_value(Pauli("ZZZ"))))
    print(f"  exact <ZZZ> (noiseless, computed via Statevector, not assumed) = {exact_val:.6f}")

    def run_shots(circuit, shots, seed):
        sim = AerSimulator(noise_model=nm, seed_simulator=seed)
        result = sim.run(transpile(circuit, sim), shots=shots).result()
        return dict(result.get_counts())

    shots = 200_000

    # naive (untwirled)
    qc_meas = qc.copy()
    qc_meas.measure_all()
    counts_naive = run_shots(qc_meas, shots, seed=1)
    def zzz_expect(counts):
        total = sum(counts.values())
        s = 0
        for bs, c in counts.items():
            parity = bs.count('1') % 2
            s += c * (1 if parity == 0 else -1)
        return s / total
    naive_val = zzz_expect(counts_naive)

    # T-REx: FULL ENUMERATION of all 2^n_qubits twirl masks (not random
    # sampling). BUG CAUGHT: an earlier version used only N_TWIRLS=8 RANDOM
    # masks, which is not guaranteed to split 50/50 per qubit at small
    # sample sizes -- T-REx's symmetrization argument requires an unbiased
    # (exactly balanced) mix of "twirled"/"not twirled" per qubit, and an
    # unlucky random draw silently violates that, biasing the calibrated
    # lambda. Full enumeration guarantees exact balance by construction
    # (exactly half of all n-bit strings have a 1 in any given position),
    # tractable here since the real target circuits are only 5 qubits
    # (2^5=32 masks).
    all_masks = [tuple(int(b) for b in format(i, f'0{n_qubits}b')) for i in range(2 ** n_qubits)]
    merged = []
    for i, mask in enumerate(all_masks):
        qct = qc.copy()
        for q in range(n_qubits):
            if mask[q]:
                qct.x(q)
        qct.measure_all()
        counts = run_shots(qct, shots // len(all_masks), seed=100 + i)
        merged.append(untwirl_counts(counts, mask, n_qubits))
    merged_counts = merge_counts(merged)
    trex_raw_val = zzz_expect(merged_counts)

    # calibration: bare |0> with the SAME full-enumeration twirl protocol
    lambdas = []
    for q in range(n_qubits):
        qc0 = QuantumCircuit(n_qubits)
        merged_cal = []
        for i, mask in enumerate(all_masks):
            qct = qc0.copy()
            for qq in range(n_qubits):
                if mask[qq]:
                    qct.x(qq)
            qct.measure_all()
            counts = run_shots(qct, shots // len(all_masks), seed=200 + q * 1000 + i)
            merged_cal.append(untwirl_counts(counts, mask, n_qubits))
        cal_counts = merge_counts(merged_cal)
        total = sum(cal_counts.values())
        p1 = sum(c for bs, c in cal_counts.items() if bs[-(q + 1)] == '1') / total
        lam_q = 1 - 2 * p1  # <Z_q> measured directly (ideal value is exactly +1)
        lambdas.append(lam_q)
    lam_zzz = float(np.prod(lambdas))
    trex_corrected_val = trex_raw_val / lam_zzz

    print(f"  naive (untwirled) measured <ZZZ> = {naive_val:.4f}  (error={abs(exact_val-naive_val):.4f})")
    print(f"  T-REx raw (twirled, uncorrected) <ZZZ> = {trex_raw_val:.4f}  (error={abs(exact_val-trex_raw_val):.4f})")
    print(f"  calibrated per-qubit lambda: {[round(l,4) for l in lambdas]}, product lambda_ZZZ={lam_zzz:.4f}")
    print(f"  T-REx corrected <ZZZ> = {trex_corrected_val:.4f}  (error={abs(exact_val-trex_corrected_val):.4f})")

    improved = abs(exact_val - trex_corrected_val) < abs(exact_val - naive_val)
    print(f"\n  T-REx correction improves accuracy: {improved} "
          f"({'PASS -- proceed to real submission' if improved else 'FAIL -- do not submit for real'})")
    return improved


# ---------------------------------------------------------------------------
# REAL SUBMISSION: T-REx (Hadamard-balanced 8-twirl design) combined with
# iteration 18's ancilla-based leakage post-selection, on the abstract
# 11-gate ansatz (this project's best structural circuit).
# ---------------------------------------------------------------------------

def transpiled_ansatz_with_ancilla(angles):
    from fixed_ansatz import build_ansatz as _build
    base = _build(angles)
    qc5 = with_ancilla_parity(base)
    return transpile(qc5, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP5, optimization_level=0)


def phase_targets():
    print("\n" + "=" * 96)
    print("  trex_readout_mitigation.py --targets  (T-REx x leakage post-selection)")
    print("=" * 96)

    p = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    assert n_ok == 36, f"fit did not converge for all 36: {n_ok}/36"
    print(f"  36/36 targets converged, worst={worst:.2e}")

    alpha_labels = p["alpha_labels"]
    identity_label = p["identity_label"]
    groups = effrag.group_labels_qubit_wise(alpha_labels)
    target_names = sorted(solutions.keys())
    masks = hadamard_balanced_masks(5, 8)
    print(f"  {len(target_names)} targets x {len(groups)} groups x {len(masks)} Hadamard-balanced twirls "
          f"= {len(target_names) * len(groups) * len(masks)} target circuits/model")
    print(f"  + {len(masks)} calibration circuits/model (bare reference, same twirl design)")

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"  connected, backend={backend.name}")

    circuits, tags = [], []
    # calibration: bare 5-qubit |00000>, no state prep at all
    cal_base = transpile(QuantumCircuit(5), basis_gates=IONQ_QIS_STANDARD_BASIS,
                          coupling_map=CMAP5, optimization_level=0)
    for mi, mask in enumerate(masks):
        qc = cal_base.copy()
        for q in range(5):
            if mask[q]:
                qc.x(q)
        qc.measure_all()
        circuits.append(qc)
        tags.append(("__calibration__", 0, mi))

    for name in target_names:
        base5 = transpiled_ansatz_with_ancilla(solutions[name]["angles"])
        for gi, group in enumerate(groups):
            combined = effrag.combined_basis_label(group)
            # BUG CAUGHT during offline verification: assumed the Hadamard(8)
            # design's mask index 0 would be all-zero (untwirled) by
            # construction, based on checking H's first COLUMN being all +1 --
            # but the masks are built from H's ROWS (one row per qubit), and
            # no single COLUMN across those rows need be all-zero (confirmed:
            # printed masks=[(1,1,1,1,1), ...], none are (0,0,0,0,0)). Fixed
            # by adding an EXPLICIT, separate untwirled circuit (tag mi=-1)
            # for the same-batch baseline comparison, rather than assuming
            # one exists inside the balanced design.
            qc_untwirled = base5.copy()
            basis_change(qc_untwirled, combined)
            qc_untwirled.measure_all()
            circuits.append(qc_untwirled)
            tags.append((name, gi, -1))
            for mi, mask in enumerate(masks):
                qc = base5.copy()
                basis_change(qc, combined)
                for q in range(5):
                    if mask[q]:
                        qc.x(q)
                qc.measure_all()
                circuits.append(qc)
                tags.append((name, gi, mi))
    print(f"  {len(circuits)} total circuits/model")

    # BUG CAUGHT on the first real submission attempt: IonQAPIError 413,
    # "Payload content length greater than maximum allowed: 10000000" --
    # 4220 5-qubit circuits (with twirls) in one submit_job() call exceeds
    # IonQ's 10MB request payload cap. Every prior submission in this
    # project (largest: 1404 circuits, iteration 18) fit under that cap;
    # this one, ~3x bigger, does not. Fixed by batching: split circuits
    # into chunks, submit each chunk as its own job (still non-blocking --
    # ALL chunks across ALL models are submitted before any .result() call,
    # preserving this project's established concurrency discipline), then
    # concatenate each model's counts back together in submission order.
    BATCH_SIZE = 500
    batches = [circuits[i:i + BATCH_SIZE] for i in range(0, len(circuits), BATCH_SIZE)]
    print(f"  batching into {len(batches)} chunks of <= {BATCH_SIZE} circuits (413 payload-size fix)")

    t0 = time.time()
    jobs = {model: [submit_job(batch, backend, model, shots=SHOTS // 8) for batch in batches]
            for model in NOISE_MODELS}
    t_submit = time.time() - t0
    print(f"  all {len(jobs) * len(batches)} batch-jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_model = {}
    for model in NOISE_MODELS:
        merged = []
        for job in jobs[model]:
            merged.extend(get_counts_list(job))
        counts_by_model[model] = merged
    t_retrieve = time.time() - t0
    print(f"  all batch-jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "exact_energy": p["exact_energy"], "noiseless_energy": p["noiseless_numpy"],
        "alpha_labels": alpha_labels, "identity_label": identity_label,
        "target_names": target_names, "groups": groups, "masks": masks,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "counts": {model: counts_by_model[model] for model in NOISE_MODELS},
        "tags": [list(t) for t in tags],
    }
    save_ckpt(CKPT_NAME, out)
    print(f"\n  --targets phase complete\n")
    return out


def untwirl_and_merge(counts_list, masks, n_qubits=5):
    merged = {}
    for counts, mask in zip(counts_list, masks):
        for bs, c in untwirl_counts(counts, mask, n_qubits).items():
            merged[bs] = merged.get(bs, 0) + c
    return merged


def zexpect_from_counts(counts, active_qubits):
    """Product-of-Z expectation over the given active qubit indices,
    from a bitstring-count dict (qiskit convention: bs[-(q+1)]=qubit q)."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    s = 0
    for bs, c in counts.items():
        parity = sum(1 for q in active_qubits if bs[-(q + 1)] == "1") % 2
        s += c * (1 if parity == 0 else -1)
    return s / total


def postselect_ancilla0(counts, ancilla_qubit=4):
    return {bs: c for bs, c in counts.items() if bs[-(ancilla_qubit + 1)] == "0"}


def label_active_qubits(label):
    n = len(label)
    return [n - 1 - i for i, ch in enumerate(label) if ch != "I"]


def assemble():
    print("\n" + "=" * 96)
    print("  trex_readout_mitigation.py --assemble")
    print("=" * 96)

    ck = load_ckpt(CKPT_NAME)
    if ck is None:
        print("  missing checkpoint -- run --targets first")
        return None

    from qforge import setup_fragment
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    alpha_labels = ck["alpha_labels"]
    identity_label = ck["identity_label"]
    target_names = ck["target_names"]
    groups = ck["groups"]
    masks = [tuple(m) for m in ck["masks"]]
    tags = [tuple(t) for t in ck["tags"]]

    idx_by_tag = {}
    for i, (name, gi, mi) in enumerate(tags):
        idx_by_tag.setdefault((name, gi), {})[mi] = i

    report = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        counts_flat = ck["counts"][model]

        # calibration: per-qubit lambda from the bare reference, twirled+untwirled+merged
        cal_counts_list = [counts_flat[idx_by_tag[("__calibration__", 0)][mi]] for mi in range(len(masks))]
        cal_merged = untwirl_and_merge(cal_counts_list, masks, n_qubits=5)
        lambdas = []
        for q in range(5):
            total = sum(cal_merged.values())
            p1 = sum(c for bs, c in cal_merged.items() if bs[-(q + 1)] == "1") / total
            lambdas.append(1 - 2 * p1)

        errs_leakage_only, errs_trex_plus_leakage = [], []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("trex", model, seed))
            raw_leak_only = {n: {} for n in target_names}
            raw_trex = {n: {} for n in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    idxs = idx_by_tag[(name, gi)]
                    counts_variants = [bootstrap_counts(counts_flat[idxs[mi]], SHOTS // 8, rng)
                                        for mi in range(len(masks))]
                    # scheme (a): leakage-only baseline -- the EXPLICIT untwirled circuit
                    # (tag mi=-1), not assumed to exist inside the Hadamard design (see the
                    # bug caught during offline verification, phase_targets()'s own comment)
                    untw0 = bootstrap_counts(counts_flat[idxs[-1]], SHOTS, rng)
                    ps0 = postselect_ancilla0(untw0)
                    for l in group:
                        raw_leak_only[name][l] = zexpect_from_counts(ps0, label_active_qubits(l))
                    # scheme (b): T-REx x leakage -- merge all 8 twirled variants (untwirled per-shot),
                    # post-select on ancilla, then divide by the calibrated per-label lambda
                    merged = untwirl_and_merge(counts_variants, masks, n_qubits=5)
                    ps_merged = postselect_ancilla0(merged)
                    for l in group:
                        active = label_active_qubits(l)
                        raw_val = zexpect_from_counts(ps_merged, active)
                        lam = float(np.prod([lambdas[q] for q in active])) if active else 1.0
                        raw_trex[name][l] = raw_val / lam if abs(lam) > 1e-6 else raw_val

            mats_a = combine_matrices(raw_leak_only, alpha_labels, identity_label, K)
            E_a, err_a = energy_from_alpha_matrices(mats_a, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                      exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs_leakage_only.append(err_a["err_vs_exact_kcal"])

            mats_b = combine_matrices(raw_trex, alpha_labels, identity_label, K)
            E_b, err_b = energy_from_alpha_matrices(mats_b, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                      exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
            errs_trex_plus_leakage.append(err_b["err_vs_exact_kcal"])

        report[model] = {
            "lambdas": lambdas,
            "leakage_only_mean_kcal": float(np.mean(errs_leakage_only)), "leakage_only_std_kcal": float(np.std(errs_leakage_only)),
            "trex_plus_leakage_mean_kcal": float(np.mean(errs_trex_plus_leakage)), "trex_plus_leakage_std_kcal": float(np.std(errs_trex_plus_leakage)),
        }
        print(f"    {model}: lambdas={[round(l,4) for l in lambdas]}")
        print(f"      leakage-only (same-batch baseline): {report[model]['leakage_only_mean_kcal']:.3f} +/- {report[model]['leakage_only_std_kcal']:.3f} kcal/mol")
        print(f"      T-REx + leakage:                    {report[model]['trex_plus_leakage_mean_kcal']:.3f} +/- {report[model]['trex_plus_leakage_std_kcal']:.3f} kcal/mol")

    print(f"\n  -- comparison to prior real-hardware findings --")
    print(f"    iteration 9 (abstract ansatz, raw, no ancilla): aria-1=34.98, forte-1=43.03 kcal/mol")
    print(f"    iteration 18 (leakage-only, separate submission): aria-1=31.77, forte-1=33.86 kcal/mol")
    print(f"    THIS run, leakage-only (same-batch): aria-1={report['aria-1']['leakage_only_mean_kcal']:.2f}, "
          f"forte-1={report['forte-1']['leakage_only_mean_kcal']:.2f} kcal/mol")
    print(f"    THIS run, T-REx + leakage: aria-1={report['aria-1']['trex_plus_leakage_mean_kcal']:.2f}, "
          f"forte-1={report['forte-1']['trex_plus_leakage_mean_kcal']:.2f} kcal/mol")
    ideal_ok = report["ideal"]["leakage_only_mean_kcal"] < 5.0
    print(f"\n  ideal correctness control: {report['ideal']['leakage_only_mean_kcal']:.3f} kcal/mol "
          f"({'PASS' if ideal_ok else 'FAIL -- pipeline bug, not noise'})")

    results = {"n_seeds": N_SEEDS, "shots_per_twirl": SHOTS // 8, "n_twirls": len(masks),
               "wall_clock": ck["wall_clock"], "report": report,
               "ideal_correctness_control_pass": bool(ideal_ok),
               "comparison": {
                   "iteration9_no_ancilla_raw": {"aria-1": 34.98, "forte-1": 43.03},
                   "iteration18_leakage_only_separate_submission": {"aria-1": 31.77, "forte-1": 33.86},
               }}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--targets", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    args = parser.parse_args()
    if args.verify:
        local_verify()
    elif args.targets:
        phase_targets()
    elif args.assemble:
        assemble()
    else:
        parser.error("pass --verify, --targets, or --assemble")


if __name__ == "__main__":
    main()
