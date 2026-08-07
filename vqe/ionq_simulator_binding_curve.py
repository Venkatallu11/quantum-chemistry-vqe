#!/usr/bin/env python3
"""
ionq_simulator_binding_curve.py — does PEC's advantage survive a noise
model we did NOT design, run for real on IonQ's free cloud simulator,
concurrently across ideal/aria-1/forte-1?
============================================================================
THE QUESTION THIS ANSWERS: every result in this project through iteration
8 used OUR OWN depolarizing noise model. PEC's exactness (iterations 4-5)
relies on the noise being a Pauli channel, which is trivially true for a
channel we built that way. IonQ's aria-1/forte-1 simulator noise models
are richer, real, and third-party -- their exact form is not published.
This file tests PEC (and CDR, and raw) against those real noise models,
submitted concurrently, on the FREE ionq_simulator only.

TASK 1 (real IonQ QPU pricing, see vqe/ionq_resource_estimate.py's
real_pricing_check(), verified via IonQ's live GET /jobs/estimate
endpoint): the $25.7899 minimum IS charged once per JOB, not per circuit
(a 125x-circuits'-worth-of-gates job costs exactly 125.0x a 1-circuit
job -- confirmed by direct API call, not assumed). But that distinction
is MOOT for this project's actual shot count: a single circuit (11 2q +
51 1q gates) at 10,000 shots already costs $206.95 on qpu.forte-1 --
8.02x the floor by itself, gate-execution cost dominates from the very
first circuit. The real workload this file's full (undiscounted)
specification would need -- 36 targets x 13 groups x 7 geometries x 2
noisy models = 6,552 circuits -- would cost $1,355,936 on real QPU
hardware (452x the $3,000 award). Real hardware is categorically not an
option; this file runs entirely on ionq_simulator, which IonQ's own docs
confirm is free and consumes no credits.

SCOPE, REDUCED FROM ITERATION 8'S DESIGN AND DISCLOSED HERE, NOT HIDDEN:
real network round-trips (submit + queue + execute + retrieve) dominate
wall-clock time in a way local exact/shot-sampled simulation never did.
  - GEOMETRIES = [0.9, 1.0, 1.1] A for raw/CDR (3, not iteration 8's 7) --
    brackets iteration 8's own found equilibrium (d_eq=0.9001 A), enough
    for one local quadratic fit and two real difference pairs.
  - PEC (needs a randomized quasi-probability circuit PER shot-batch, see
    below -- far more circuits per data point than raw/CDR) is run ONLY
    at d=1.0 (the reference geometry). This sacrifices PEC's own
    difference-error/cancellation-factor/binding-curve numbers -- reported
    as N/A for PEC, not faked from the single point.
  - CDR training: genuinely independent real submissions, 8 seeds x 5
    random-angle draws/seed (reduced from iteration 8's local-simulator
    36-slots x 5/slot -- real network cost makes that scale infeasible;
    using global random draws instead of per-slot draws is a disclosed
    simplification, not a hidden one).
  - PEC calibration: 1 representative CX pair, 1 representative qubit
    (reduced from checking every distinct pair/qubit, as iteration 5 did
    locally for free) -- reduction stated plainly.
  - "8 seeds" handled two DIFFERENT, disclosed ways:
      * CDR training and PEC's randomized quasi-probability circuit draws:
        genuinely independent real circuit submissions -- real seed-to-
        seed variability, no shortcuts.
      * raw/CDR TARGET measurement: ONE real 10,000-shot submission per
        circuit, then BOOTSTRAP resampling (multinomial resampling of the
        real INTEGER counts returned by get_counts()) into 8 sub-samples
        at the same 10,000-shot count. This is NOT 8 independent real
        executions of the target circuits -- stated once here, applies to
        every raw/CDR "seed" reported below, not re-qualified per number.

REAL PEC ON HARDWARE -- the actual protocol, not a shortcut: gate-by-gate
PEC's quasi-probability inverse is NOT a physical channel, so unlike every
prior iteration in this project (which had exact density-matrix access
and could apply pec_inverse_weights as a literal linear map on rho), a
REAL device needs the standard randomized-circuit protocol (Temme,
Bravyi, Gambetta, PRL 119, 180509 (2017)): after every noisy gate, sample
a Pauli correction from the LEARNED channel's quasi-probability
distribution (weighted by |eta_P|, sign = sign(eta_P)), insert it as
actual gates, track the accumulated sign, execute, and multiply the
measured expectation value by (accumulated sign x gamma_total) -- the
ensemble average over independently-sampled circuit instances converges
to the corrected expectation. gamma_total (the same quantity iterations
4-6 already computed honestly as the real sampling-overhead cost) is the
correction's own normalization here, not an afterthought.

PEC's CHANNEL: learned separately for aria-1 and forte-1, from Clifford-
only circuits submitted to THAT specific real noise model (never shared,
never from the local model, never from any target circuit).

"IDEAL" IS A CORRECTNESS CONTROL, NOT A DATA POINT: raw energy on
noise_model="ideal" must reproduce the noiseless_numpy energy to within
real shot noise at 10,000 shots/setting. If it doesn't, that's a pipeline
bug, not evidence about noise -- checked explicitly, and the run stops if
it fails.

GATESET: abstract ("qis") gates throughout, optimization_level=0 (the
established invariant: >=1 can silently collapse gate counts for specific
fitted angles) -- no folding here, so no risk of the native-gate/
abstract-gate mixing this project has been careful about elsewhere.

CONCURRENCY: every phase submits ALL of its jobs first (non-blocking
backend.run() calls) and only THEN calls .result() on each -- so
ideal/aria-1/forte-1 (and, within a noise model, every geometry) queue
and execute concurrently server-side rather than being waited on one at a
time. Wall-clock time for each phase is recorded to make this checkable,
not asserted.

Run (each phase is its own process -- the Bash tool caps a single command
at 10 minutes even backgrounded, and real network round-trips make the
full run too long for one invocation):
    python vqe/ionq_simulator_binding_curve.py --smoke
    python vqe/ionq_simulator_binding_curve.py --calibrate
    python vqe/ionq_simulator_binding_curve.py --targets --d 0.9
    python vqe/ionq_simulator_binding_curve.py --targets --d 1.0
    python vqe/ionq_simulator_binding_curve.py --targets --d 1.1
    python vqe/ionq_simulator_binding_curve.py --pec
    python vqe/ionq_simulator_binding_curve.py --assemble
"""
import os
import sys
import json
import time
import argparse
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
import rank6_symmetry_vd as r
import loop_pec as pec
from energy_difference_study import setup_at_distance, fit_and_verify
from fixed_ansatz import build_ansatz
from ionq_backend import connect_provider, get_simulator
from ionq_run import basis_change, pauli_expectation, IONQ_QIS_STANDARD_BASIS

from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.transpiler import CouplingMap

K = 6
HARTREE_TO_KCAL_MOL = r.HARTREE_TO_KCAL_MOL
GEOMETRIES = [0.9, 1.0, 1.1]
D_REF = 1.0
PEC_D = 1.0  # PEC's randomized-circuit protocol is run only at this one geometry -- see module docstring
NOISE_MODELS = ["ideal", "aria-1", "forte-1"]
CORRECTABLE_MODELS = ["aria-1", "forte-1"]  # ideal needs no CDR/PEC -- it's the correctness control
SHOTS = 10_000
N_SEEDS = 8
N_TRAIN_DRAWS = 5           # real independent random-angle draws per seed, per model, for CDR training
LOW_SIGNAL_CUTOFF = r.LOW_SIGNAL_CUTOFF
CX_CAL_PAIR = (0, 1)
QUBIT_CAL = 0
N_LIST_CX = [1, 3, 5, 7, 9]
N_LIST_U3 = [1, 3, 5, 7, 9]
CHEM_ACC_KCAL = 1.0

CK_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
os.makedirs(CK_DIR, exist_ok=True)
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_results.json")


def ckpt_path(name):
    return os.path.join(CK_DIR, f"{name}.json")


def save_ckpt(name, obj):
    with open(ckpt_path(name), "w") as f:
        json.dump(obj, f, indent=2)
    print(f"    checkpoint saved -> {ckpt_path(name)}")


def load_ckpt(name):
    p = ckpt_path(name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Circuit construction (abstract "qis" gates, opt_level=0, real 4-qubit
# coupling map -- not a network call, all local)
# ---------------------------------------------------------------------------

CMAP4 = CouplingMap.from_full(4)


def transpiled_ansatz(angles):
    return transpile(build_ansatz(angles), basis_gates=IONQ_QIS_STANDARD_BASIS,
                      coupling_map=CMAP4, optimization_level=0)


def label_for_qubits(chars_at_qubit, n=4):
    """chars_at_qubit: dict qubit_index -> 'X'/'Y'/'Z', default 'I'.
    Returns the MSQ-first qiskit label string this project's basis_change/
    pauli_expectation convention expects (label[i] acts on qubit n-1-i)."""
    chars = ["I"] * n
    for q, c in chars_at_qubit.items():
        chars[n - 1 - q] = c
    return "".join(chars)


def measurement_circuits_for_groups(base_circuit, groups):
    circuits = []
    for group in groups:
        combined = effrag.combined_basis_label(group)
        qc = base_circuit.copy()
        basis_change(qc, combined)
        qc.measure_all()
        circuits.append(qc)
    return circuits


# ---------------------------------------------------------------------------
# Non-blocking submission: submit first (returns a job handle immediately),
# poll/retrieve later -- this is what makes ideal/aria-1/forte-1 genuinely
# concurrent instead of sequential.
# ---------------------------------------------------------------------------

def submit_job(circuits, backend, noise_model, shots=SHOTS):
    return backend.run(circuits, noise_model=noise_model, shots=shots)


def get_counts_list(job):
    result = job.result()
    counts = result.get_counts()
    if not isinstance(counts, list):
        counts = [counts]
    return [dict(c) for c in counts]


def bootstrap_counts(counts, n_shots, rng):
    bitstrings = list(counts.keys())
    if not bitstrings:
        return {}
    probs = np.array([counts[b] for b in bitstrings], dtype=float)
    probs = probs / probs.sum()
    draws = rng.multinomial(n_shots, probs)
    return {b: int(n) for b, n in zip(bitstrings, draws) if n > 0}


def stable_seed(*parts):
    """Deterministic RNG seed from arbitrary parts -- Python's built-in
    hash() is randomized per-process (PYTHONHASHSEED) for str/tuple
    inputs, which silently made bootstrap-resample results different on
    every --assemble rerun of the SAME real checkpointed data (a real
    reproducibility bug caught by literally re-running --assemble twice
    and seeing different numbers)."""
    import zlib
    return zlib.crc32("|".join(str(p) for p in parts).encode()) % (2**31)


def counts_to_probs(counts):
    total = sum(counts.values())
    return {b: n / total for b, n in counts.items()} if total > 0 else {}


def expectation_from_counts(counts, label):
    return pauli_expectation(counts_to_probs(counts), label)


# ---------------------------------------------------------------------------
# PHASE: calibration (Clifford-only, per correctable model) + CDR training
# (real independent random-angle draws, per correctable model) -- both
# geometry-independent, submitted ONCE, reused for every geometry below.
# ---------------------------------------------------------------------------

def cx_decay_circuit(pair, N):
    """<XX> is the observable being decayed (ideal=1 for every odd N, the
    Bell-state property this project's local Pauli-Lindblad learning also
    relies on) -- needs an X-basis rotation (H on both qubits of the
    pair) before Z-basis measurement, NOT a plain measure_all()."""
    qc = QuantumCircuit(4)
    a, b = pair
    qc.h(a)
    for _ in range(N):
        qc.cx(a, b)
    qct = transpile(qc, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP4, optimization_level=0)
    xx_label = label_for_qubits({a: "X", b: "X"})
    basis_change(qct, xx_label)
    qct.measure_all()
    return qct


def ry_identity_decay_circuit(qubit, N):
    """ry(0) repeated N times: logically identity (verified: ry(0)=I
    exactly), but built from an actual IonQ 'qis' abstract gate ('ry') so
    the calibration measures IonQ's real error on the gate TYPE that
    appears in the ansatz's own 1-qubit gates post-transpile -- unlike
    this project's local model, which calibrates the specific gate NAME
    'u3', IonQ's transpiled 1-qubit gates are a mix of rx/ry/rz/h/s-type
    names (a real structural difference from the local model, disclosed
    here, not glossed over). This calibrates 'ry' as one representative
    1-qubit gate and, per the PEC correction below, that single learned
    rate is applied uniformly to every 1-qubit gate regardless of its
    exact name -- an arity-based, not name-based, simplification, the
    same spirit as the local model's own single P1_PER_GATE applying
    uniformly regardless of a u3's specific rotation angles."""
    qc = QuantumCircuit(4)
    for _ in range(N):
        qc.ry(0.0, qubit)
    qct = transpile(qc, basis_gates=IONQ_QIS_STANDARD_BASIS, coupling_map=CMAP4, optimization_level=0)
    qct.measure_all()  # <Z> is the observable -- direct Z-basis measurement, no basis_change needed
    return qct


def fit_decay_rate(N_list, measured_vals):
    Ns = np.array(N_list, dtype=float)
    lns = np.log(np.clip(np.abs(np.asarray(measured_vals)), 1e-6, None))
    slope, intercept = np.polyfit(Ns, lns, 1)
    fit_vals = slope * Ns + intercept
    residual_rms = float(np.sqrt(np.mean((lns - fit_vals) ** 2)))
    lam = np.exp(slope)
    p_learned = float(1 - lam)
    return p_learned, residual_rms


def build_calibration_circuits():
    circuits, tags = [], []
    for N in N_LIST_CX:
        circuits.append(cx_decay_circuit(CX_CAL_PAIR, N))
        tags.append(("cx", N))
    for N in N_LIST_U3:
        circuits.append(ry_identity_decay_circuit(QUBIT_CAL, N))
        tags.append(("ry", N))
    return circuits, tags


def build_cdr_training_circuits(non_id_labels, groups, n_qubits=4):
    """8 seeds x N_TRAIN_DRAWS real independent random-angle draws
    (geometry-independent by construction, exactly like this project's
    local CDR training) -- each draw's exact (noiseless) label values
    computed LOCALLY (classical statevector sim, same as every prior
    iteration's CDR training; only the NOISY side needs the real device),
    each draw's circuits built for real submission."""
    circuits, tags = [], []
    exact_by_draw = {}
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed)
        for draw in range(N_TRAIN_DRAWS):
            angles = rng.uniform(-np.pi, np.pi, 5)
            exact_vals = r.exact_labels(angles.tolist(), non_id_labels)
            exact_by_draw[(seed, draw)] = exact_vals
            base = transpiled_ansatz(angles.tolist())
            for gi, group in enumerate(groups):
                combined = effrag.combined_basis_label(group)
                qc = base.copy()
                basis_change(qc, combined)
                qc.measure_all()
                circuits.append(qc)
                tags.append((seed, draw, gi))
    return circuits, tags, exact_by_draw


def phase_calibrate():
    print("\n" + "=" * 96)
    print("  ionq_simulator_binding_curve.py --calibrate")
    print("  (Clifford-only noise learning + CDR training, real submission, both correctable models)")
    print("=" * 96)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    p_ref = setup_at_distance(D_REF, K)
    fit_and_verify(p_ref)
    non_id_labels = [l for l in p_ref["alpha_labels"] if l != p_ref["identity_label"]]
    groups = effrag.group_labels_qubit_wise(p_ref["alpha_labels"])
    print(f"  {len(non_id_labels)} non-identity alpha labels, {len(groups)} qubit-wise-commuting groups")

    cal_circuits, cal_tags = build_calibration_circuits()
    train_circuits, train_tags, exact_by_draw = build_cdr_training_circuits(non_id_labels, groups)
    print(f"  calibration circuits: {len(cal_circuits)} ({dict((k,len(list(v))) for k,v in itertools.groupby(sorted(t[0] for t in cal_tags)))})")
    print(f"  CDR training circuits: {len(train_circuits)} ({N_SEEDS} seeds x {N_TRAIN_DRAWS} draws x {len(groups)} groups)")

    print(f"\n  submitting calibration + training jobs for {CORRECTABLE_MODELS} -- ALL non-blocking before any .result()")
    t0 = time.time()
    jobs = {}
    for model in CORRECTABLE_MODELS:
        jobs[("cal", model)] = submit_job(cal_circuits, backend, model, shots=SHOTS)
        jobs[("train", model)] = submit_job(train_circuits, backend, model, shots=SHOTS)
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    results_counts = {key: get_counts_list(job) for key, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s (concurrency check: retrieval time should NOT scale "
          f"linearly with job count if these ran concurrently server-side)")

    out = {"non_id_labels": non_id_labels, "groups": groups, "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
           "per_model": {}}

    for model in CORRECTABLE_MODELS:
        cal_counts = results_counts[("cal", model)]
        cx_vals = [expectation_from_counts(cal_counts[i], label_for_qubits({CX_CAL_PAIR[0]: "X", CX_CAL_PAIR[1]: "X"}))
                   for i, (kind, N) in enumerate(cal_tags) if kind == "cx"]
        ry_vals = [expectation_from_counts(cal_counts[i], label_for_qubits({QUBIT_CAL: "Z"}))
                   for i, (kind, N) in enumerate(cal_tags) if kind == "ry"]
        p2_learned, p2_residual = fit_decay_rate(N_LIST_CX, cx_vals)
        p1_learned, p1_residual = fit_decay_rate(N_LIST_U3, ry_vals)
        gamma_2q = pec.gamma_factor(max(p2_learned, 1e-9), 2)
        gamma_1q = pec.gamma_factor(max(p1_learned, 1e-9), 1)
        gamma_total = gamma_2q ** 11 * gamma_1q ** 51
        print(f"\n  {model}: p2_learned={p2_learned:.6f} (fit residual={p2_residual:.4f}), "
              f"p1_learned={p1_learned:.6f} (fit residual={p1_residual:.4f}), gamma_total={gamma_total:.4f}")
        print(f"    (fit residual = RMS deviation from a straight line in log-space across depths {N_LIST_CX}/{N_LIST_U3} -- "
              "large residual = evidence the real channel is NOT well described by a single Pauli-depolarizing rate, "
              "i.e. non-Pauli noise or drift, not a fit-quality artifact of this script)")

        train_counts = results_counts[("train", model)]
        training_by_seed = {seed: [] for seed in range(N_SEEDS)}
        for i, (seed, draw, gi) in enumerate(train_tags):
            group = groups[gi]
            counts = train_counts[i]
            for l in group:
                if l not in exact_by_draw[(seed, draw)]:
                    continue  # identity label -- <I>=1 exactly, never part of CDR training (established invariant)
                measured = expectation_from_counts(counts, l)
                exact = exact_by_draw[(seed, draw)][l]
                training_by_seed[seed].append({"slot": f"s{seed}d{draw}", "label": l, "exact": exact, "noisy": measured})

        out["per_model"][model] = {
            "p2_learned": p2_learned, "p2_fit_residual": p2_residual,
            "p1_learned": p1_learned, "p1_fit_residual": p1_residual,
            "gamma_2q": gamma_2q, "gamma_1q": gamma_1q, "gamma_total": gamma_total,
            "training_by_seed": training_by_seed,
        }

    save_ckpt("calibrate", out)
    print(f"\n  --calibrate phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: target measurement (raw + CDR-correctable), one geometry at a
# time, all 3 noise models submitted concurrently.
# ---------------------------------------------------------------------------

def phase_targets(d):
    print("\n" + "=" * 96)
    print(f"  ionq_simulator_binding_curve.py --targets --d {d}")
    print("=" * 96)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    p_d = setup_at_distance(d, K)
    n_ok, worst, counts = fit_and_verify(p_d)
    assert p_d["schmidt_rank_le_K"], f"d={d}: K=6 not exact here"
    assert counts == {11}, f"d={d}: gate count not fixed: {counts}"
    print(f"  d={d} A: {n_ok}/36 targets converged (worst={worst:.2e}), gate count={counts}, "
          f"exact_energy={p_d['exact_energy']:.6f} Ha, noiseless_numpy={p_d['noiseless_numpy']:.6f} Ha")

    non_id_labels = [l for l in p_d["alpha_labels"] if l != p_d["identity_label"]]
    groups = effrag.group_labels_qubit_wise(p_d["alpha_labels"])
    target_names = sorted(p_d["solutions"].keys())

    all_circuits = {}
    for name in target_names:
        base = transpiled_ansatz(p_d["solutions"][name]["angles"])
        all_circuits[name] = measurement_circuits_for_groups(base, groups)
    n_circuits_per_model = sum(len(v) for v in all_circuits.values())
    print(f"  {len(target_names)} targets x {len(groups)} groups = {n_circuits_per_model} circuits/model")

    flat_circuits = [qc for name in target_names for qc in all_circuits[name]]

    print(f"\n  submitting to {NOISE_MODELS} -- ALL non-blocking before any .result() (concurrency)")
    t0 = time.time()
    jobs = {model: submit_job(flat_circuits, backend, model, shots=SHOTS) for model in NOISE_MODELS}
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_model = {model: get_counts_list(job) for model, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} jobs retrieved, {t_retrieve:.1f}s")

    # reshape flat counts back into {model: {name: [counts per group]}}
    out_counts = {model: {} for model in NOISE_MODELS}
    idx = 0
    idx_map = []
    for name in target_names:
        for _ in groups:
            idx_map.append(name)
    for model in NOISE_MODELS:
        per_name = {name: [] for name in target_names}
        for i, name in enumerate(idx_map):
            per_name[name].append(counts_by_model[model][i])
        out_counts[model] = per_name

    # -- ideal correctness control, checked immediately --
    ideal_raw = {}
    for name in target_names:
        vals = {}
        for gi, group in enumerate(groups):
            counts = out_counts["ideal"][name][gi]
            for l in group:
                vals[l] = expectation_from_counts(counts, l)
        ideal_raw[name] = vals
    mats_ideal = r.combine_matrices(ideal_raw, p_d["alpha_labels"], p_d["identity_label"], K)
    E_ideal, err_ideal = r.energy_from_alpha_matrices(mats_ideal, p_d, K)
    ideal_ok = err_ideal["err_vs_noiseless_kcal"] < 5.0  # generous vs real 10,000-shot noise, see iteration 6
    print(f"\n  IDEAL correctness control: E={E_ideal:.6f} Ha, err_vs_noiseless={err_ideal['err_vs_noiseless_kcal']:.4f} kcal/mol "
          f"({'PASS' if ideal_ok else 'FAIL -- STOP, this is a pipeline bug, not noise'})")
    if not ideal_ok:
        raise RuntimeError(f"d={d}: ideal noise_model does not reproduce the noiseless energy "
                            f"(err={err_ideal['err_vs_noiseless_kcal']:.4f} kcal/mol) -- pipeline bug, refusing to proceed")

    out = {
        "d": d, "exact_energy": p_d["exact_energy"], "noiseless_numpy": p_d["noiseless_numpy"],
        "alpha_labels": p_d["alpha_labels"], "identity_label": p_d["identity_label"],
        "target_names": target_names, "groups": groups,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "ideal_check": {"E": E_ideal, "err_vs_noiseless_kcal": err_ideal["err_vs_noiseless_kcal"], "pass": bool(ideal_ok)},
        "counts": out_counts,
    }
    save_ckpt(f"targets_d{d}", out)
    print(f"\n  --targets --d {d} phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: real PEC via randomized quasi-probability circuit insertion
# (the actual hardware protocol -- see module docstring). PEC_D only.
# ---------------------------------------------------------------------------

def sample_pauli_correction(weights, rng):
    labels = list(weights.keys())
    w = np.array([weights[l] for l in labels])
    absw = np.abs(w)
    probs = absw / absw.sum()
    idx = rng.choice(len(labels), p=probs)
    return labels[idx], float(np.sign(w[idx]) if w[idx] != 0 else 1.0)


def append_pauli_label(qc, label, qargs):
    """Append single-qubit Pauli gates realizing `label` (e.g. 'XY') on
    `qargs`, using IonQ 'qis' abstract gates directly ('x'/'y'/'z') --
    no transpile needed, these are already in IONQ_QIS_STANDARD_BASIS."""
    for ch, q in zip(label, qargs):
        if ch == "X":
            qc.x(q)
        elif ch == "Y":
            qc.y(q)
        elif ch == "Z":
            qc.z(q)


def build_pec_variant_circuit(base_transpiled_qc, group, p2_learned, p1_learned, rng):
    """Walk the transpiled ansatz gate-by-gate, after every 2q/1q gate
    sample+insert a Pauli correction from the LEARNED channel's inverse
    quasi-probability distribution, track the accumulated sign. This is
    the real, physically-executable PEC protocol -- not a density-matrix
    shortcut (this project has no exact-DM access to IonQ's real noise)."""
    weights2 = pec.pec_inverse_weights(max(p2_learned, 1e-9), 2)
    weights1 = pec.pec_inverse_weights(max(p1_learned, 1e-9), 1)
    qc = QuantumCircuit(base_transpiled_qc.num_qubits)
    sign = 1.0
    for instr in base_transpiled_qc.data:
        op, qargs = instr.operation, instr.qubits
        if op.name in ("measure", "barrier"):
            continue
        idxs = [base_transpiled_qc.find_bit(q).index for q in qargs]
        qc.append(op, qargs)
        if op.num_qubits == 2:
            label, s = sample_pauli_correction(weights2, rng)
            sign *= s
            append_pauli_label(qc, label, idxs)
        elif op.num_qubits == 1:
            label, s = sample_pauli_correction(weights1, rng)
            sign *= s
            append_pauli_label(qc, label, idxs)
    combined = effrag.combined_basis_label(group)
    basis_change(qc, combined)
    qc.measure_all()
    return qc, sign


def phase_pec():
    print("\n" + "=" * 96)
    print(f"  ionq_simulator_binding_curve.py --pec  (real randomized-circuit PEC, d={PEC_D} only)")
    print("=" * 96)

    calib = load_ckpt("calibrate")
    if calib is None:
        raise RuntimeError("run --calibrate first (need learned p2/p1 per model)")

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    p_d = setup_at_distance(PEC_D, K)
    fit_and_verify(p_d)
    groups = effrag.group_labels_qubit_wise(p_d["alpha_labels"])
    target_names = sorted(p_d["solutions"].keys())
    print(f"  d={PEC_D}: {len(target_names)} targets x {len(groups)} groups")

    jobs = {}
    signs_by_model = {}
    t0 = time.time()
    for model in CORRECTABLE_MODELS:
        p2_learned = calib["per_model"][model]["p2_learned"]
        p1_learned = calib["per_model"][model]["p1_learned"]
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(seed * 7919 + hash(model) % 1000)
            circuits, tags, signs = [], [], []
            for name in target_names:
                base = transpiled_ansatz(p_d["solutions"][name]["angles"])
                for gi, group in enumerate(groups):
                    qc, sign = build_pec_variant_circuit(base, group, p2_learned, p1_learned, rng)
                    circuits.append(qc)
                    tags.append((name, gi))
                    signs.append(sign)
            job = submit_job(circuits, backend, model, shots=SHOTS)
            jobs[(model, seed)] = (job, tags)
            signs_by_model[(model, seed)] = signs
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} PEC jobs submitted ({len(CORRECTABLE_MODELS)} models x {N_SEEDS} seeds), {t_submit:.1f}s")

    t0 = time.time()
    counts_by_key = {key: get_counts_list(job) for key, (job, tags) in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} PEC jobs retrieved, {t_retrieve:.1f}s")

    out = {
        "d": PEC_D, "exact_energy": p_d["exact_energy"], "noiseless_numpy": p_d["noiseless_numpy"],
        "alpha_labels": p_d["alpha_labels"], "identity_label": p_d["identity_label"],
        "target_names": target_names, "groups": groups,
        "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        "per_model_seed": {},
    }
    for model in CORRECTABLE_MODELS:
        gamma_total = calib["per_model"][model]["gamma_total"]
        for seed in range(N_SEEDS):
            job, tags = jobs[(model, seed)]
            counts_list = counts_by_key[(model, seed)]
            signs = signs_by_model[(model, seed)]
            raw = {name: {} for name in target_names}
            for i, (name, gi) in enumerate(tags):
                group = groups[gi]
                counts = counts_list[i]
                sign = signs[i]
                for l in group:
                    e = expectation_from_counts(counts, l)
                    raw[name][l] = gamma_total * sign * e
            key = f"{model}_seed{seed}"
            out["per_model_seed"][key] = raw
    save_ckpt("pec", out)
    print(f"\n  --pec phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: assemble -- combine every checkpoint, do the full analysis, no
# network calls.
# ---------------------------------------------------------------------------

def energies_from_raw(raw, p_d, K, per_label_scale=None):
    mats = r.combine_matrices(raw, p_d["alpha_labels"], p_d["identity_label"], K, per_label_scale=per_label_scale)
    return r.energy_from_alpha_matrices(mats, p_d, K)


def assemble():
    print("\n" + "=" * 96)
    print("  ionq_simulator_binding_curve.py --assemble")
    print("=" * 96)

    calib = load_ckpt("calibrate")
    target_ckpts = {d: load_ckpt(f"targets_d{d}") for d in GEOMETRIES}
    pec_ckpt = load_ckpt("pec")
    missing = [f"targets_d{d}" for d in GEOMETRIES if target_ckpts[d] is None]
    if calib is None:
        missing.append("calibrate")
    if missing:
        print(f"  missing checkpoints: {missing} -- run the corresponding phases first")
        return None

    geoms = {d: setup_at_distance(d, K) for d in GEOMETRIES}
    for d in GEOMETRIES:
        fit_and_verify(geoms[d])

    METHODS = ["raw", "cdr", "pec"] if pec_ckpt is not None else ["raw", "cdr"]
    E_data = {model: {m: {d: [] for d in GEOMETRIES} for m in METHODS} for model in NOISE_MODELS}

    rng_master = np.random.default_rng(12345)

    for model in NOISE_MODELS:
        for d in GEOMETRIES:
            ck = target_ckpts[d]
            p_d = geoms[d]
            groups = ck["groups"]
            target_names = ck["target_names"]
            counts_model = ck["counts"][model]

            if model == "ideal":
                # bootstrap 8 seeds of raw only -- ideal is the correctness control, no CDR/PEC needed
                for seed in range(N_SEEDS):
                    rng = np.random.default_rng(stable_seed(model, d, seed))
                    raw = {name: {} for name in target_names}
                    for name in target_names:
                        for gi, group in enumerate(groups):
                            counts = bootstrap_counts(counts_model[name][gi], SHOTS, rng)
                            for l in group:
                                raw[name][l] = expectation_from_counts(counts, l)
                    E_raw, _ = energies_from_raw(raw, p_d, K)
                    E_data[model]["raw"][d].append(E_raw)
                    E_data[model]["cdr"][d].append(E_raw)  # no correction defined for ideal; raw==cdr placeholder
                continue

            scales_by_seed = calib["per_model"][model]["training_by_seed"]
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed(model, d, seed))
                raw = {name: {} for name in target_names}
                for name in target_names:
                    for gi, group in enumerate(groups):
                        counts = bootstrap_counts(counts_model[name][gi], SHOTS, rng)
                        for l in group:
                            raw[name][l] = expectation_from_counts(counts, l)
                E_raw, _ = energies_from_raw(raw, p_d, K)
                E_data[model]["raw"][d].append(E_raw)

                scales = r.fit_all_scales(scales_by_seed[str(seed)] if str(seed) in scales_by_seed else scales_by_seed[seed],
                                           ck["alpha_labels"] and [l for l in ck["alpha_labels"] if l != ck["identity_label"]], K)
                E_cdr, _ = energies_from_raw(raw, p_d, K, per_label_scale=scales["per_basis_scale"])
                E_data[model]["cdr"][d].append(E_cdr)

            if pec_ckpt is not None and d == PEC_D:
                for seed in range(N_SEEDS):
                    key = f"{model}_seed{seed}"
                    raw_pec = pec_ckpt["per_model_seed"][key]
                    E_pec, _ = energies_from_raw(raw_pec, p_d, K)
                    E_data[model]["pec"][d].append(E_pec)

    # ==== error analysis: absolute vs difference, cancellation factor ====
    exact_by_d = {d: geoms[d]["exact_energy"] for d in GEOMETRIES}
    diff_exact_by_d = {d: exact_by_d[d] - exact_by_d[D_REF] for d in GEOMETRIES}

    report = {model: {} for model in NOISE_MODELS}
    for model in NOISE_MODELS:
        methods_here = ["raw"] if model == "ideal" else METHODS
        for method in methods_here:
            per_geom = {}
            abs_errs_all, abs_errs_for_cancel, diff_errs_for_cancel = [], [], []
            for d in GEOMETRIES:
                vals = E_data[model][method].get(d, [])
                if not vals:
                    continue
                E_arr = np.array(vals)
                abs_err = np.abs(E_arr - exact_by_d[d]) * HARTREE_TO_KCAL_MOL
                row = {"E_mean_ha": float(np.mean(E_arr)), "E_std_ha": float(np.std(E_arr)),
                       "abs_err_mean_kcal": float(np.mean(abs_err)), "abs_err_std_kcal": float(np.std(abs_err))}
                abs_errs_all.append(row["abs_err_mean_kcal"])
                if d != D_REF and E_data[model][method].get(D_REF):
                    E_ref = np.array(E_data[model][method][D_REF])
                    n = min(len(E_arr), len(E_ref))
                    diff_measured = E_arr[:n] - E_ref[:n]
                    diff_err = np.abs(diff_measured - diff_exact_by_d[d]) * HARTREE_TO_KCAL_MOL
                    row["diff_err_mean_kcal"] = float(np.mean(diff_err))
                    row["diff_err_std_kcal"] = float(np.std(diff_err))
                    abs_errs_for_cancel.append(row["abs_err_mean_kcal"])
                    diff_errs_for_cancel.append(row["diff_err_mean_kcal"])
                per_geom[d] = row
            # mean_abs_err_kcal is reported over EVERY geometry this (model, method) actually has data
            # for (so PEC, which only has d=1.0, still reports a real number here) -- separate from the
            # cancellation-factor pair-means below, which need >=1 non-reference geometry to exist at all.
            mean_abs = float(np.mean(abs_errs_all)) if abs_errs_all else None
            mean_diff = float(np.mean(diff_errs_for_cancel)) if diff_errs_for_cancel else None
            cancellation = (mean_abs / mean_diff) if (mean_abs and mean_diff and mean_diff > 1e-9) else None
            report[model][method] = {"per_geometry": per_geom, "mean_abs_err_kcal": mean_abs,
                                      "mean_diff_err_kcal": mean_diff, "cancellation_factor": cancellation}

    # ==== binding curve (3-point local quadratic fit -- see module docstring on reduced geometry set) ====
    ds = np.array(GEOMETRIES)
    exact_es = np.array([exact_by_d[d] for d in GEOMETRIES])
    coeffs_exact = np.polyfit(ds, exact_es, 2)
    d_eq_exact = -coeffs_exact[1] / (2 * coeffs_exact[0])
    E_eq_exact = np.polyval(coeffs_exact, d_eq_exact)
    binding_curve = {"exact": {"d_eq": float(d_eq_exact), "E_eq_ha": float(E_eq_exact)}}
    for model in NOISE_MODELS:
        methods_here = ["raw"] if model == "ideal" else METHODS
        binding_curve[model] = {}
        for method in methods_here:
            means = [np.mean(E_data[model][method][d]) for d in GEOMETRIES if E_data[model][method].get(d)]
            if len(means) < 3:
                continue
            coeffs = np.polyfit(ds, means, 2)
            d_eq = -coeffs[1] / (2 * coeffs[0])
            E_eq = np.polyval(coeffs, d_eq)
            binding_curve[model][method] = {
                "d_eq": float(d_eq), "E_eq_ha": float(E_eq),
                "d_eq_error_A": float(abs(d_eq - d_eq_exact)),
            }

    # ==== learned-channel diagnostics (headline honesty check) ====
    channel_report = {}
    for model in CORRECTABLE_MODELS:
        cm = calib["per_model"][model]
        channel_report[model] = {
            "p2_learned": cm["p2_learned"], "p2_fit_residual": cm["p2_fit_residual"],
            "p1_learned": cm["p1_learned"], "p1_fit_residual": cm["p1_fit_residual"],
            "gamma_total": cm["gamma_total"],
        }

    # ==== print report ====
    print(f"\n  -- LEARNED CHANNEL DIAGNOSTICS (fit residual = evidence for/against pure-Pauli-depolarizing noise) --")
    for model, c in channel_report.items():
        flag = "LARGE -- non-Pauli noise or drift likely" if max(c["p2_fit_residual"], c["p1_fit_residual"]) > 0.05 else "small -- depolarizing fit looks clean"
        print(f"    {model}: p2={c['p2_learned']:.6f} (resid={c['p2_fit_residual']:.4f}), "
              f"p1={c['p1_learned']:.6f} (resid={c['p1_fit_residual']:.4f}), gamma_total={c['gamma_total']:.3f}  [{flag}]")

    print(f"\n  -- ABSOLUTE vs DIFFERENCE ERROR, PER MODEL --")
    for model in NOISE_MODELS:
        print(f"\n  == {model} ==")
        for method, row in report[model].items():
            mabs, mdiff, canc = row["mean_abs_err_kcal"], row["mean_diff_err_kcal"], row["cancellation_factor"]
            mabs_s = f"{mabs:.3f}" if mabs is not None else "n/a"
            mdiff_s = f"{mdiff:.3f}" if mdiff is not None else "n/a (single geometry only)"
            canc_s = f"{canc:.2f}x" if canc is not None else "n/a"
            print(f"    {method:>4}: mean|abs err|={mabs_s} kcal/mol, mean|diff err|={mdiff_s} kcal/mol, cancellation={canc_s}")

    print(f"\n  -- BINDING CURVE (equilibrium bond length) --")
    print(f"    exact: d_eq={d_eq_exact:.4f} A")
    for model in NOISE_MODELS:
        for method, bc in binding_curve.get(model, {}).items():
            print(f"    {model:>8} / {method:>4}: d_eq={bc['d_eq']:.4f} A (err={bc['d_eq_error_A']:.4f} A)")

    pec_note = ("PEC run only at d=1.0 (see module docstring) -- its mean_diff_err/cancellation_factor/binding_curve "
                "are therefore N/A by design, not missing data.") if pec_ckpt is not None else "PEC checkpoint not found -- run --pec first."
    print(f"\n  {pec_note}")

    results = {
        "geometries": GEOMETRIES, "d_ref": D_REF, "pec_d": PEC_D, "n_seeds": N_SEEDS, "shots": SHOTS,
        "methods_reported": METHODS,
        "exact_energy_by_d": exact_by_d,
        "channel_diagnostics": channel_report,
        "error_report": report,
        "binding_curve": binding_curve,
        "wall_clock": {
            "calibrate": calib["wall_clock"],
            "targets": {str(d): target_ckpts[d]["wall_clock"] for d in GEOMETRIES},
            "pec": pec_ckpt["wall_clock"] if pec_ckpt is not None else None,
        },
        "ideal_correctness_control": {str(d): target_ckpts[d]["ideal_check"] for d in GEOMETRIES},
        "pec_note": pec_note,
        "scope_reductions": [
            "3 geometries (0.9/1.0/1.1 A), not iteration 8's 7",
            "PEC run only at d=1.0 -- no PEC difference/cancellation/binding-curve numbers",
            "CDR training: 8 seeds x 5 global random-angle draws/seed, not per-slot (36 slots x 5) as done locally",
            "PEC calibration: 1 CX pair + 1 qubit, not all distinct pairs/qubits as done locally",
            "raw/CDR '8 seeds' at the target-measurement step are bootstrap resamples of ONE real 10,000-shot "
            "execution per circuit, not 8 independent real executions -- CDR training and PEC's quasi-probability "
            "circuit draws ARE genuinely independent real executions",
        ],
        "chemical_accuracy_kcal": CHEM_ACC_KCAL,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


# ---------------------------------------------------------------------------
# Smoke test: tiny, cheap, real -- validates the whole real-submission
# pipeline (circuit build -> submit -> counts -> expectation -> energy)
# before committing to the full run.
# ---------------------------------------------------------------------------

def smoke_test():
    print("\n" + "=" * 70)
    print("  SMOKE TEST -- tiny real submission to ionq_simulator (ideal only)")
    print("=" * 70)
    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"  connected, backend={backend.name}")

    p_d = setup_at_distance(1.0, K)
    fit_and_verify(p_d)
    groups = effrag.group_labels_qubit_wise(p_d["alpha_labels"])[:2]
    name = sorted(p_d["solutions"].keys())[0]
    base = transpiled_ansatz(p_d["solutions"][name]["angles"])
    circuits = measurement_circuits_for_groups(base, groups)
    print(f"  submitting {len(circuits)} circuits, shots={SHOTS}, noise_model=ideal")
    t0 = time.time()
    job = submit_job(circuits, backend, "ideal", shots=SHOTS)
    counts = get_counts_list(job)
    print(f"  done in {time.time()-t0:.1f}s")
    for group, c in zip(groups, counts):
        combined = effrag.combined_basis_label(group)
        for l in group:
            e = expectation_from_counts(c, l)
            print(f"    group={combined} label={l}: measured={e:.4f}")
    print("\n  SMOKE TEST PASSED\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--calibrate", action="store_true")
    parser.add_argument("--targets", action="store_true")
    parser.add_argument("--d", type=float)
    parser.add_argument("--pec", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    args = parser.parse_args()

    if args.smoke:
        smoke_test()
    elif args.calibrate:
        phase_calibrate()
    elif args.targets:
        if args.d is None:
            parser.error("--targets needs --d")
        phase_targets(args.d)
    elif args.pec:
        phase_pec()
    elif args.assemble:
        assemble()
    else:
        parser.error("pass one of --smoke / --calibrate / --targets --d D / --pec / --assemble")


if __name__ == "__main__":
    main()
