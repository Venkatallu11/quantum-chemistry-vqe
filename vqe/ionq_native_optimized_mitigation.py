#!/usr/bin/env python3
"""
ionq_native_optimized_mitigation.py — does the TrappedIonOptimizerPlugin
gate-count reduction (Task 3 of the prior iteration: mean 9.28 vs 11
two-qubit gates, K=6, per-target non-uniform 4-11) actually help when
combined with real ZNE/CDR/PEC on IonQ's free simulators, concurrently?
============================================================================
Task 3 verified the gate-count reduction LOCALLY only, and cited OLDER
real data using a DIFFERENT circuit (native_stateprep.py's hand-derived
K=5 tree, not this project's fixed 11-gate K=6 ansatz). This file is the
missing real test: the ACTUAL TrappedIonOptimizerPlugin-optimized fixed
ansatz, submitted for real to ideal/aria-1/forte-1 CONCURRENTLY, combined
with native ZNE (gate folding, correct order verified: optimize FIRST,
fold the already-optimized circuit SECOND, never re-optimize a folded
circuit -- re-optimizing was independently verified elsewhere in this
project to cancel about half of every fold), CDR, and PEC.

HONEST EXPECTATION, stated before running, not after: the gate-count
reduction is real but modest (~15% fewer two-qubit gates on average).
Every real-hardware number in this project so far (raw 35-135 kcal/mol,
CDR worse than raw, PEC's best real result 32 kcal/mol) sits 30-40x above
chemical accuracy. A 15% gate reduction is not expected to close that gap
on its own -- this run answers whether it helps at all and by how much,
not whether it reaches chemical accuracy.

A REAL COMPLICATION, handled explicitly, not glossed over: the optimized
circuit's gate count VARIES per target (4-11), which breaks CDR's usual
structural-identity premise (training circuits are supposed to be
structurally identical to targets). This file does NOT dodge that -- CDR
training draws are EACH individually optimized too (so each training
circuit has its own natural post-optimization structure, exactly like a
target would), and the per-basis scale fit is applied and reported AS
MEASURED, including if it fails to generalize. That failure, if it
happens, is itself the honest answer to "does CDR still work here."

SCOPE, reduced from iteration 9's full design for real-network-time
reasons, disclosed: single geometry (d=1.0), CDR training 8 seeds x 3
draws/seed (not 5), PEC real randomized-circuit protocol at 4 seeds
(not 8) -- both still genuinely independent real executions, just fewer
of them. Raw + ZNE (the primary "does gate count help" question) run at
full scope (all 36 K=6 targets x 13 qubit-wise groups x folds 1/3/5).

Run:
    python vqe/ionq_native_optimized_mitigation.py --calibrate
    python vqe/ionq_native_optimized_mitigation.py --targets
    python vqe/ionq_native_optimized_mitigation.py --pec
    python vqe/ionq_native_optimized_mitigation.py --assemble
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
import rank6_symmetry_vd as r
import loop_pec as pec
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, native_target
from ionq_fold_check import fold_native_2q
from ionq_native_forged_energy import native_basis_change
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, HARTREE_TO_KCAL_MOL, bootstrap_counts, stable_seed, save_ckpt, load_ckpt,
)

from qiskit.transpiler import PassManagerConfig
from qiskit_ionq import TrappedIonOptimizerPlugin
from qiskit_ionq.ionq_gates import MSGate, GPIGate, GPI2Gate
from qiskit.quantum_info import Statevector

K = 6
GATE = "ms"
NOISE_MODELS = ["ideal", "aria-1", "forte-1"]
CORRECTABLE_MODELS = ["aria-1", "forte-1"]
FOLDS = [1, 3, 5]
N_TRAIN_DRAWS = 3
PEC_SEEDS = 4
N_LIST_MS = [1, 3, 5, 7, 9]
N_LIST_1Q = [1, 3, 5, 7, 9]
CHEM_ACC_KCAL = 1.0
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "ionq_native_optimized_mitigation_results.json")


# ---------------------------------------------------------------------------
# Circuit construction: build -> native -> OPTIMIZE -> (fold) -> measure
# ---------------------------------------------------------------------------

def optimized_native_circuit(angles):
    native = to_native(build_ansatz(angles), GATE)
    tgt = native_target(native.num_qubits, GATE)
    pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt), optimization_level=3)
    optimized = pm.run(native)
    return optimized


def native_measurement_circuit_for_group(base, group):
    combined = effrag.combined_basis_label(group)
    qc = base.compose(native_basis_change(combined, GATE))
    qc.measure_all()
    return qc


def gate_counts(qc):
    ops = qc.count_ops()
    return ops.get(GATE, 0), sum(ops.get(g, 0) for g in ("gpi", "gpi2"))


def submit_job(circuits, backend, noise_model, shots=SHOTS):
    return backend.run(circuits, noise_model=noise_model, shots=shots)


def get_counts_list(job):
    result = job.result()
    counts = result.get_counts()
    if not isinstance(counts, list):
        counts = [counts]
    return [dict(c) for c in counts]


def counts_to_probs(counts):
    total = sum(counts.values())
    return {b: n / total for b, n in counts.items()} if total > 0 else {}


def expectation_from_counts(counts, label):
    return pauli_expectation(counts_to_probs(counts), label)


# ---------------------------------------------------------------------------
# PHASE: calibrate -- native MS/1q Clifford decay (real channel learning)
# + CDR training (each draw individually optimized)
# ---------------------------------------------------------------------------

def ms_decay_circuit(pair, N):
    """N copies of MS(0,0,0.25) (full strength -- verified in Task 3 to be
    the ONLY strength this optimizer emits) on `pair`, starting from
    |0000>. MS(0,0,theta) = exp(-i*theta*pi*XX) (verified by direct matrix
    comparison against the qiskit_ionq-published MS matrix in Task 3's
    XXPlusYY-vs-MS check), which preserves the |00>/|11> parity subspace
    -- so <ZZ> on `pair` is EXACTLY +1 for any N under ideal (noiseless)
    execution, checked below rather than assumed, giving the same
    "ideal=1 at every depth, real noise decays it" calibration signal
    iteration 9's CX-decay probe used, with a direct Z-basis measurement
    (no basis-change gates needed at all)."""
    from qiskit.circuit import QuantumCircuit
    a, b = pair
    qc = QuantumCircuit(4)
    for _ in range(N):
        qc.append(MSGate(0, 0, 0.25), [a, b])
    return qc


def one_q_decay_circuit(qubit, N):
    """N copies of GPI(0)-GPI(0) (=X.X=I, verified below) as the
    logically-identity repeated single-qubit probe -- same role as
    iteration 9's U3Gate(0,0,0) probe, native gates instead."""
    from qiskit.circuit import QuantumCircuit
    qc = QuantumCircuit(4)
    for _ in range(N):
        qc.append(GPIGate(0.0), [qubit])
        qc.append(GPIGate(0.0), [qubit])
    return qc


def verify_calibration_circuits():
    sv = np.asarray(Statevector.from_instruction(one_q_decay_circuit(0, 3)))
    assert abs(abs(sv[0]) - 1.0) < 1e-9, "one_q_decay_circuit is not logically identity"
    for N in (1, 2, 3, 4, 5):
        sv2 = np.asarray(Statevector.from_instruction(ms_decay_circuit((0, 1), N)))
        probs = np.abs(sv2) ** 2
        # <ZZ> = sum over bitstrings of probs * (+1 if qubit0==qubit1 else -1); qubits 0,1 are the two LSBs
        zz = 0.0
        for i, p in enumerate(probs):
            b0, b1 = (i >> 0) & 1, (i >> 1) & 1
            zz += p * (1 if b0 == b1 else -1)
        assert abs(zz - 1.0) < 1e-9, f"ideal <ZZ> at N={N} is {zz}, expected 1.0 -- MS calibration probe design is wrong"


def build_calibration_circuits():
    verify_calibration_circuits()
    circuits, tags = [], []
    for N in N_LIST_MS:
        qc = ms_decay_circuit((0, 1), N)
        qc.measure_all()
        circuits.append(qc)
        tags.append(("ms", N))
    for N in N_LIST_1Q:
        qc = one_q_decay_circuit(0, N)
        qc.measure_all()
        circuits.append(qc)
        tags.append(("1q", N))
    return circuits, tags


def build_cdr_training_circuits(non_id_labels, groups):
    circuits, tags = [], []
    exact_by_draw = {}
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed + 500)
        for draw in range(N_TRAIN_DRAWS):
            angles = rng.uniform(-np.pi, np.pi, 5).tolist()
            sv = Statevector.from_instruction(build_ansatz(angles))
            from qiskit.quantum_info import Pauli
            exact_vals = {l: float(sv.expectation_value(Pauli(l)).real) for l in non_id_labels}
            exact_by_draw[(seed, draw)] = exact_vals
            base = optimized_native_circuit(angles)
            for gi, group in enumerate(groups):
                circuits.append(native_measurement_circuit_for_group(base, group))
                tags.append((seed, draw, gi))
    return circuits, tags, exact_by_draw


def label_for_qubits(chars_at_qubit, n=4):
    chars = ["I"] * n
    for q, c in chars_at_qubit.items():
        chars[n - 1 - q] = c
    return "".join(chars)


def fit_decay_rate(N_list, measured_vals):
    Ns = np.array(N_list, dtype=float)
    lns = np.log(np.clip(np.abs(np.asarray(measured_vals)), 1e-6, None))
    slope, intercept = np.polyfit(Ns, lns, 1)
    fit_vals = slope * Ns + intercept
    residual_rms = float(np.sqrt(np.mean((lns - fit_vals) ** 2)))
    p_learned = float(1 - np.exp(slope))
    return p_learned, residual_rms


def phase_calibrate():
    print("\n" + "=" * 96)
    print("  ionq_native_optimized_mitigation.py --calibrate")
    print("=" * 96)

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"\n  connected, backend={backend.name} (native gateset)")

    p_ref = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p_ref["targets"])
    p_ref["solutions"] = solutions
    non_id_labels = [l for l in p_ref["alpha_labels"] if l != p_ref["identity_label"]]
    groups = effrag.group_labels_qubit_wise(p_ref["alpha_labels"])
    print(f"  {n_ok}/36 targets converged, {len(non_id_labels)} non-id labels, {len(groups)} groups")

    cal_circuits, cal_tags = build_calibration_circuits()
    train_circuits, train_tags, exact_by_draw = build_cdr_training_circuits(non_id_labels, groups)
    print(f"  calibration circuits: {len(cal_circuits)}, CDR training circuits: {len(train_circuits)} "
          f"({N_SEEDS} seeds x {N_TRAIN_DRAWS} draws x {len(groups)} groups)")

    t0 = time.time()
    jobs = {}
    for model in CORRECTABLE_MODELS:
        jobs[("cal", model)] = submit_job(cal_circuits, backend, model)
        jobs[("train", model)] = submit_job(train_circuits, backend, model)
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    results_counts = {key: get_counts_list(job) for key, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

    out = {"non_id_labels": non_id_labels, "groups": groups,
           "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve}, "per_model": {}}

    for model in CORRECTABLE_MODELS:
        cal_counts = results_counts[("cal", model)]
        ms_vals = [expectation_from_counts(cal_counts[i], label_for_qubits({0: "Z", 1: "Z"}))
                   for i, (kind, N) in enumerate(cal_tags) if kind == "ms"]
        q_vals = [expectation_from_counts(cal_counts[i], label_for_qubits({0: "Z"}))
                  for i, (kind, N) in enumerate(cal_tags) if kind == "1q"]
        p_ms_learned, p_ms_resid = fit_decay_rate(N_LIST_MS, ms_vals)
        p_1q_learned, p_1q_resid = fit_decay_rate(N_LIST_1Q, q_vals)
        gamma_ms = pec.gamma_factor(max(p_ms_learned, 1e-9), 2)
        gamma_1q = pec.gamma_factor(max(p_1q_learned, 1e-9), 1)
        print(f"\n  {model}: p_ms={p_ms_learned:.6f} (resid={p_ms_resid:.4f}), "
              f"p_1q={p_1q_learned:.6f} (resid={p_1q_resid:.4f})")

        train_counts = results_counts[("train", model)]
        training_by_seed = {seed: [] for seed in range(N_SEEDS)}
        for i, (seed, draw, gi) in enumerate(train_tags):
            group = groups[gi]
            counts = train_counts[i]
            for l in group:
                if l not in exact_by_draw[(seed, draw)]:
                    continue
                training_by_seed[seed].append({"slot": f"s{seed}d{draw}", "label": l,
                                                "exact": exact_by_draw[(seed, draw)][l],
                                                "noisy": expectation_from_counts(counts, l)})

        out["per_model"][model] = {
            "p_ms_learned": p_ms_learned, "p_ms_fit_residual": p_ms_resid,
            "p_1q_learned": p_1q_learned, "p_1q_fit_residual": p_1q_resid,
            "gamma_ms": gamma_ms, "gamma_1q": gamma_1q,
            "training_by_seed": training_by_seed,
        }

    save_ckpt("native_opt_calibrate", out)
    print(f"\n  --calibrate phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: targets -- raw + ZNE (native fold of the ALREADY-optimized
# circuit, never re-optimized after folding), all 3 models concurrent
# ---------------------------------------------------------------------------

def phase_targets():
    print("\n" + "=" * 96)
    print("  ionq_native_optimized_mitigation.py --targets")
    print("=" * 96)

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    p_d = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p_d["targets"])
    p_d["solutions"] = solutions
    assert n_ok == 36
    groups = effrag.group_labels_qubit_wise(p_d["alpha_labels"])
    target_names = sorted(solutions.keys())
    print(f"  36/36 converged, {len(groups)} groups, exact_energy={p_d['exact_energy']:.6f} Ha, "
          f"noiseless_numpy(K=6)={p_d['noiseless_numpy']:.6f} Ha")

    gate_count_by_target = {}
    circuits_by_fold = {f: [] for f in FOLDS}
    tags = []
    for name in target_names:
        optimized = optimized_native_circuit(solutions[name]["angles"])
        n2q, n1q = gate_counts(optimized)
        gate_count_by_target[name] = {"n2q": n2q, "n1q": n1q}
        for f in FOLDS:
            folded = fold_native_2q(optimized, f, GATE)
            for gi, group in enumerate(groups):
                circuits_by_fold[f].append(native_measurement_circuit_for_group(folded, group))
        for gi, group in enumerate(groups):
            tags.append((name, gi))
    print(f"  per-target native 2q gate count: min={min(v['n2q'] for v in gate_count_by_target.values())} "
          f"max={max(v['n2q'] for v in gate_count_by_target.values())} "
          f"mean={np.mean([v['n2q'] for v in gate_count_by_target.values()]):.2f} (abstract baseline: 11)")
    print(f"  {len(circuits_by_fold[1])} circuits/fold/model x {len(FOLDS)} folds")

    t0 = time.time()
    jobs = {}
    for model in NOISE_MODELS:
        for f in FOLDS:
            jobs[(model, f)] = submit_job(circuits_by_fold[f], backend, model)
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_key = {key: get_counts_list(job) for key, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "exact_energy": p_d["exact_energy"], "noiseless_numpy": p_d["noiseless_numpy"],
        "alpha_labels": p_d["alpha_labels"], "identity_label": p_d["identity_label"],
        "target_names": target_names, "groups": groups, "gate_count_by_target": gate_count_by_target,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "counts": {model: {f: counts_by_key[(model, f)] for f in FOLDS} for model in NOISE_MODELS},
        "tags": [list(t) for t in tags],
    }
    save_ckpt("native_opt_targets", out)
    print(f"\n  --targets phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: pec -- real randomized quasi-probability circuits on the learned
# native MS/1q channel, per-target gate structure (non-uniform, handled
# per-circuit not assumed uniform)
# ---------------------------------------------------------------------------

def sample_pauli_correction(weights, rng):
    labels = list(weights.keys())
    w = np.array([weights[l] for l in labels])
    absw = np.abs(w)
    probs = absw / absw.sum()
    idx = rng.choice(len(labels), p=probs)
    return labels[idx], float(np.sign(w[idx]) if w[idx] != 0 else 1.0)


def build_pec_variant_circuit(optimized_circuit, group, p_ms, p_1q, rng):
    """Walk the OPTIMIZED (per-target, non-uniform) native circuit
    gate-by-gate, insert a sampled Pauli correction (via native gates)
    after every ms/gpi/gpi2 instruction, track sign -- same real,
    physically-executable PEC protocol as iteration 9, adapted to
    per-target-varying native gate sequences instead of a fixed
    11-CX/51-1q abstract structure."""
    weights2 = pec.pec_inverse_weights(max(p_ms, 1e-9), 2)
    weights1 = pec.pec_inverse_weights(max(p_1q, 1e-9), 1)
    from qiskit.circuit import QuantumCircuit
    qc = QuantumCircuit(optimized_circuit.num_qubits)
    sign = 1.0
    for instr in optimized_circuit.data:
        op, qargs = instr.operation, instr.qubits
        if op.name in ("measure", "barrier"):
            continue
        idxs = [optimized_circuit.find_bit(q).index for q in qargs]
        qc.append(op, qargs)
        if op.num_qubits == 2:
            label, s = sample_pauli_correction(weights2, rng)
            sign *= s
            append_pauli_label_z_safe(qc, label, idxs)
        elif op.num_qubits == 1:
            label, s = sample_pauli_correction(weights1, rng)
            sign *= s
            append_pauli_label_z_safe(qc, label, idxs)
    combined = effrag.combined_basis_label(group)
    qc = qc.compose(native_basis_change(combined, GATE))
    qc.measure_all()
    return qc, sign


def append_pauli_label_z_safe(qc, label, qargs):
    """Verified by DIRECT matrix comparison (not derived from memory)
    before use here: GPI(0) == X exactly, GPI(0.25) == Y exactly. Z is
    built as GPI(0.25) then GPI(0) (matrix product X.Y = i*Z) -- matches
    Z up to a global phase of i, which is exactly irrelevant here since
    only RELATIVE signs from pec_inverse_weights matter for a quasi-
    probability correction's measurement statistics."""
    for ch, q in zip(label, qargs):
        if ch == "X":
            qc.append(GPIGate(0.0), [q])
        elif ch == "Y":
            qc.append(GPIGate(0.25), [q])
        elif ch == "Z":
            qc.append(GPIGate(0.25), [q])
            qc.append(GPIGate(0.0), [q])


def phase_pec():
    print("\n" + "=" * 96)
    print("  ionq_native_optimized_mitigation.py --pec")
    print("=" * 96)

    calib = load_ckpt("native_opt_calibrate")
    if calib is None:
        raise RuntimeError("run --calibrate first")

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    p_d = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p_d["targets"])
    p_d["solutions"] = solutions
    groups = effrag.group_labels_qubit_wise(p_d["alpha_labels"])
    target_names = sorted(solutions.keys())

    jobs = {}
    signs_by_key, gate_counts_by_key = {}, {}
    t0 = time.time()
    for model in CORRECTABLE_MODELS:
        p_ms = calib["per_model"][model]["p_ms_learned"]
        p_1q = calib["per_model"][model]["p_1q_learned"]
        for seed in range(PEC_SEEDS):
            rng = np.random.default_rng(stable_seed("native_pec", model, seed))
            circuits, tags, signs, n2q_list, n1q_list = [], [], [], [], []
            for name in target_names:
                optimized = optimized_native_circuit(solutions[name]["angles"])
                n2q, n1q = gate_counts(optimized)
                for gi, group in enumerate(groups):
                    qc, sign = build_pec_variant_circuit(optimized, group, p_ms, p_1q, rng)
                    circuits.append(qc)
                    tags.append((name, gi))
                    signs.append(sign)
                    n2q_list.append(n2q)
                    n1q_list.append(n1q)
            job = submit_job(circuits, backend, model)
            jobs[(model, seed)] = (job, tags)
            signs_by_key[(model, seed)] = signs
            gate_counts_by_key[(model, seed)] = list(zip(n2q_list, n1q_list))
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} PEC jobs submitted ({len(CORRECTABLE_MODELS)} models x {PEC_SEEDS} seeds), {t_submit:.1f}s")

    t0 = time.time()
    counts_by_key = {key: get_counts_list(job) for key, (job, tags) in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} PEC jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "exact_energy": p_d["exact_energy"], "noiseless_numpy": p_d["noiseless_numpy"],
        "alpha_labels": p_d["alpha_labels"], "identity_label": p_d["identity_label"],
        "target_names": target_names, "groups": groups,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "per_model_seed": {},
    }
    for model in CORRECTABLE_MODELS:
        gamma_ms = calib["per_model"][model]["gamma_ms"]
        gamma_1q = calib["per_model"][model]["gamma_1q"]
        for seed in range(PEC_SEEDS):
            job, tags = jobs[(model, seed)]
            counts_list = counts_by_key[(model, seed)]
            signs = signs_by_key[(model, seed)]
            gc_list = gate_counts_by_key[(model, seed)]
            raw = {name: {} for name in target_names}
            for i, (name, gi) in enumerate(tags):
                group = groups[gi]
                counts = counts_list[i]
                sign = signs[i]
                n2q, n1q = gc_list[i]
                gamma_total = gamma_ms ** n2q * gamma_1q ** n1q
                for l in group:
                    e = expectation_from_counts(counts, l)
                    raw[name][l] = gamma_total * sign * e
            out["per_model_seed"][f"{model}_seed{seed}"] = raw
    save_ckpt("native_opt_pec", out)
    print(f"\n  --pec phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: assemble
# ---------------------------------------------------------------------------

def assemble():
    print("\n" + "=" * 96)
    print("  ionq_native_optimized_mitigation.py --assemble")
    print("=" * 96)

    calib = load_ckpt("native_opt_calibrate")
    targets = load_ckpt("native_opt_targets")
    pec_ck = load_ckpt("native_opt_pec")
    if calib is None or targets is None:
        print("  missing checkpoints -- run --calibrate and --targets first")
        return None

    p_d = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p_d["targets"])
    p_d["solutions"] = solutions
    exact_energy = targets["exact_energy"]
    noiseless_numpy = targets["noiseless_numpy"]
    groups = targets["groups"]
    target_names = targets["target_names"]
    alpha_labels = targets["alpha_labels"]
    identity_label = targets["identity_label"]

    idx_map = []
    for name in target_names:
        for _ in groups:
            idx_map.append(name)

    def energies_for(model, fold):
        counts_flat = targets["counts"][model][str(fold)]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])
        results = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("native_opt_targets", model, fold, seed))
            raw = {name: {} for name in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    counts = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                    for l in group:
                        raw[name][l] = expectation_from_counts(counts, l)
            mats = r.combine_matrices(raw, alpha_labels, identity_label, K)
            E, _ = r.energy_from_alpha_matrices(mats, p_d, K)
            results.append(E)
        return results

    print(f"\n  -- RAW + ZNE, per model --")
    report = {}
    for model in NOISE_MODELS:
        energies_by_fold = {f: energies_for(model, f) for f in FOLDS}
        raw_errs = [abs(e - exact_energy) * HARTREE_TO_KCAL_MOL for e in energies_by_fold[1]]
        zne_lin_errs, zne_quad_errs = [], []
        for seed in range(N_SEEDS):
            es = [energies_by_fold[f][seed] for f in FOLDS]
            lin = np.polyval(np.polyfit(FOLDS, es, 1), 0)
            quad = np.polyval(np.polyfit(FOLDS, es, 2), 0)
            zne_lin_errs.append(abs(lin - exact_energy) * HARTREE_TO_KCAL_MOL)
            zne_quad_errs.append(abs(quad - exact_energy) * HARTREE_TO_KCAL_MOL)
        report[model] = {
            "raw_mean_kcal": float(np.mean(raw_errs)), "raw_std_kcal": float(np.std(raw_errs)),
            "zne_linear_mean_kcal": float(np.mean(zne_lin_errs)), "zne_linear_std_kcal": float(np.std(zne_lin_errs)),
            "zne_quadratic_mean_kcal": float(np.mean(zne_quad_errs)), "zne_quadratic_std_kcal": float(np.std(zne_quad_errs)),
        }
        print(f"    {model}: raw={report[model]['raw_mean_kcal']:.2f}+/-{report[model]['raw_std_kcal']:.2f}  "
              f"zne_lin={report[model]['zne_linear_mean_kcal']:.2f}+/-{report[model]['zne_linear_std_kcal']:.2f}  "
              f"zne_quad={report[model]['zne_quadratic_mean_kcal']:.2f}+/-{report[model]['zne_quadratic_std_kcal']:.2f}")

    print(f"\n  -- CDR, per correctable model --")
    for model in CORRECTABLE_MODELS:
        scales_by_seed = calib["per_model"][model]["training_by_seed"]
        counts_flat = targets["counts"][model]["1"]
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_flat[i])
        cdr_errs = []
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(stable_seed("native_opt_cdr", model, seed))
            raw = {name: {} for name in target_names}
            for name in target_names:
                for gi, group in enumerate(groups):
                    counts = bootstrap_counts(per_name[name][gi], SHOTS, rng)
                    for l in group:
                        raw[name][l] = expectation_from_counts(counts, l)
            training = scales_by_seed[str(seed)] if str(seed) in scales_by_seed else scales_by_seed[seed]
            non_id_labels = [l for l in alpha_labels if l != identity_label]
            scales = r.fit_all_scales(training, non_id_labels, K)
            mats = r.combine_matrices(raw, alpha_labels, identity_label, K, per_label_scale=scales["per_basis_scale"])
            E, _ = r.energy_from_alpha_matrices(mats, p_d, K)
            cdr_errs.append(abs(E - exact_energy) * HARTREE_TO_KCAL_MOL)
        report[model]["cdr_mean_kcal"] = float(np.mean(cdr_errs))
        report[model]["cdr_std_kcal"] = float(np.std(cdr_errs))
        print(f"    {model}: cdr={report[model]['cdr_mean_kcal']:.2f}+/-{report[model]['cdr_std_kcal']:.2f}")

    if pec_ck is not None:
        print(f"\n  -- PEC, per correctable model --")
        for model in CORRECTABLE_MODELS:
            pec_errs = []
            for seed in range(PEC_SEEDS):
                raw = pec_ck["per_model_seed"][f"{model}_seed{seed}"]
                mats = r.combine_matrices(raw, alpha_labels, identity_label, K)
                E, _ = r.energy_from_alpha_matrices(mats, p_d, K)
                pec_errs.append(abs(E - exact_energy) * HARTREE_TO_KCAL_MOL)
            report[model]["pec_mean_kcal"] = float(np.mean(pec_errs))
            report[model]["pec_std_kcal"] = float(np.std(pec_errs))
            print(f"    {model}: pec={report[model]['pec_mean_kcal']:.2f}+/-{report[model]['pec_std_kcal']:.2f} "
                  f"({PEC_SEEDS} seeds)")

    print(f"\n  -- comparison to prior real-hardware findings --")
    for model in ("aria-1", "forte-1"):
        print(f"    {model}: abstract-ansatz raw (iteration 9) vs native-optimized raw (this run) = "
              f"{report[model]['raw_mean_kcal']:.2f} kcal/mol")

    results = {"K": K, "n_seeds": N_SEEDS, "pec_seeds": PEC_SEEDS, "shots": SHOTS, "folds": FOLDS,
               "exact_energy_ha": exact_energy, "noiseless_numpy_ha": noiseless_numpy,
               "gate_count_by_target": targets["gate_count_by_target"],
               "report": report, "chemical_accuracy_kcal": CHEM_ACC_KCAL}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibrate", action="store_true")
    parser.add_argument("--targets", action="store_true")
    parser.add_argument("--pec", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    args = parser.parse_args()
    if args.calibrate:
        phase_calibrate()
    elif args.targets:
        phase_targets()
    elif args.pec:
        phase_pec()
    elif args.assemble:
        assemble()
    else:
        parser.error("pass --calibrate / --targets / --pec / --assemble")


if __name__ == "__main__":
    main()
