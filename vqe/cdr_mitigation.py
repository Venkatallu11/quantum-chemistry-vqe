#!/usr/bin/env python3
"""
cdr_mitigation.py — Clifford-Data-Regression-style noise-scale learning for
the H4 forged energy, using the FIXED-STRUCTURE ansatz from fixed_ansatz.py
(v2: 11 two-qubit gates, cut from the original 33 -- see fixed_ansatz.py's
docstring) so training and target circuits are structurally IDENTICAL.
SIMULATOR ONLY. Nothing here touches IonQ.
============================================================================
WHY THE OLD ATTEMPT FAILED (context, not re-derived here): the earliest CDR
attempt used qiskit's `StatePreparation`, which compiles a different circuit
per input vector, so training and target noise profiles never matched.
fixed_ansatz.py fixes that: every training circuit and every target circuit
is `build_ansatz(angles)`, the identical gate sequence, angles the only
thing that varies.

TWO CHANGES FROM THE PRECEDING VERSION OF THIS FILE, BOTH REQUIRED:

  1. PER-GATE CALIBRATION, NOT PER-CIRCUIT. The preceding version root-found
     a single 2-qubit depolarizing rate so the FULL 33-gate ansatz's forged
     energy reproduced f=0.842805 (the real aria-1 fold=1 measurement).
     That's circular: it implies a per-gate error of 0.842805^(1/33) =
     0.517%, when the number f=0.842805 was actually measured on the
     14-GATE native circuit, implying 0.842805^(1/14) = 1.214% per gate --
     2.35x higher. Forcing a fixed TOTAL f onto however many gates a new
     ansatz happens to have quietly dilutes the per-gate rate as gate count
     grows, flattering any circuit that got smaller. Fixed here: apply the
     universal per-gate rate P2_PER_GATE=0.01214 directly (via
     qiskit-aer's depolarizing_error, whose parameter IS the fractional
     Pauli-expectation shrink per application -- verified empirically
     below, not assumed from a hand-derived formula: an earlier derivation
     assuming a Pauli-twirl-over-nonidentity-Paulis convention predicted
     shrink=1-(16/15)*param and was WRONG; qiskit's actual convention
     gives shrink=1-param exactly, confirmed via a direct Aer probe before
     use). No more root-finding -- whatever f the ansatz's own gate count
     produces is what gets reported and used.
  2. SINGLE-QUBIT GATES COUNTED TOO, at P1_PER_GATE = P2_PER_GATE/40 (per
     fixed_ansatz.py's honest native-gate accounting: 163-397 single-qubit
     GPi/GPi2 gates after native transpile + IonQ's own
     TrappedIonOptimizerPlugin is not negligible at that rate). Applied
     here as a depolarizing_error on every 'u3' gate in the LOCAL
     simulator's own u3/cx-transpiled circuit (verified 'u3' really is
     the gate name qiskit's transpiler emits for basis_gates=["u3","cx"],
     not assumed).
  3. SEED SWEEP, not one seed. A single training-seed run overstates what
     CDR reliably achieves. Every scheme below is run across N_SEEDS
     independent training-data draws; the headline numbers are mean +/-
     std over seeds, not a cherry-picked best case.

THIRD FIX, found while building zne_vs_cdr.py: every transpile call here
uses optimization_level=0, not 1. Verified (not assumed): at
optimization_level>=1, qiskit's two-qubit gate synthesis is numerically
ADAPTIVE -- for specific fitted-angle solutions that happen to land near
periodic special values (multiples of pi, which several of the 25 targets'
least_squares solutions did), it silently collapses the circuit to FEWER
CX gates (7-9 instead of 11 for 4 of the 25 targets, measured directly).
That means 4 target circuits were getting LESS noise than every training
circuit (which almost never lands on a special angle), a real violation of
the "structurally identical" premise this whole file exists to satisfy.
optimization_level=0 (pure rule-based substitution, no adaptive synthesis)
was confirmed to give a constant 11 CX across all 25 targets AND 10 random
training-style angle draws before adopting it here -- see fixed_ansatz.py's
own now-added self-check (abstract_cx_count_fixed_across_25_targets).

TWO FIXES CARRIED FORWARD FROM EARLIER (already-debugged) ATTEMPTS -- both
silent-corruption bugs if skipped:
  1. NEVER rescale the identity Pauli label. <I>=1 exactly for any
     trace-preserving channel -- dividing it by a learned scale corrupts
     the forged energy (multiplies straight into alpha*beta for every term
     carrying that label). First attempt at this overshot to -2.67 Ha.
  2. DROP training pairs whose EXACT expectation value has |value| < 0.05
     as uninformative about the noise scale.

Run:
    python vqe/cdr_mitigation.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
from fixed_ansatz import build_ansatz, P2_PER_GATE, P1_PER_GATE, fidelity_from_counts

from qiskit import transpile
from qiskit.quantum_info import Pauli, Statevector
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error

HARTREE_TO_KCAL_MOL = 627.5094740631
ATOMS = [0, 1, 2, 3]
NELEC = 4
K = 5
BASIS_GATES = ["u3", "cx"]
LOW_SIGNAL_CUTOFF = 0.05
N_TRAIN_PER_SLOT = 5
N_SEEDS = 8
SEEDS = list(range(N_SEEDS))

# Context from the earliest (StatePreparation-based) CDR attempt -- reported
# purely for comparison, NOT reproduced/recomputed here (that code no
# longer exists).
PRIOR_STATEPREP_RAW_KCAL = 199
PRIOR_STATEPREP_GLOBAL_KCAL = 116
PRIOR_STATEPREP_GAIN = PRIOR_STATEPREP_RAW_KCAL / PRIOR_STATEPREP_GLOBAL_KCAL

ANSATZ_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "fixed_ansatz_results.json")
NATIVE_ZNE_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "native_forged_zne_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "fixed_ansatz_v2_results.json")


# ---------------------------------------------------------------------------
# Setup: same fragment physics as ionq_native_forged_energy.py, but the
# angle solutions come from fixed_ansatz.py's fixed circuit.
# ---------------------------------------------------------------------------

def setup():
    qop_bare, qop_pen, enuc = effrag.build_fragment_qop(ATOMS, NELEC, d=1.0)
    n_qubits = qop_bare.num_qubits
    e_elec, psi = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    psi_real, residual = effrag.real_gauge(psi)
    lambdas, u_vecs, v_vecs = effrag.schmidt_decompose_real(psi_real, n_qubits)
    terms = effrag.decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    assert set(alpha_labels) == set(beta_labels)

    u_top, v_top = u_vecs[:K], v_vecs[:K]
    signs = np.array([1.0 if np.dot(v_top[n], u_top[n]) >= 0 else -1.0 for n in range(K)])
    sign_residual = max(float(np.max(np.abs(v_top[n] - signs[n] * u_top[n]))) for n in range(K))
    assert sign_residual < 1e-8, f"beta=sign*alpha residual {sign_residual:.3e} -- shortcut doesn't hold here"

    identity_label = "I" * (n_qubits // 2)
    identity_coeff = 0.0
    for label, coeff in qop_bare.to_list():
        if label == "I" * n_qubits:
            identity_coeff = float(coeff.real)
            break
    e_mixed = identity_coeff + enuc

    alpha_cache, beta_cache = effrag.precompute_exact_matrices(terms, u_top, v_top)
    noiseless_numpy = effrag.ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)

    with open(ANSATZ_RESULTS_PATH) as f:
        ansatz_results = json.load(f)
    solutions = ansatz_results["solutions"]

    expected_names = [f"u_{n}" for n in range(K)]
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        expected_names += [f"(u{n}+u{m})", f"(u{n}-u{m})"]
    assert set(expected_names) == set(solutions.keys()), "fixed_ansatz_results.json solutions don't match expected 25 slot names"
    assert all(solutions[name]["verified"] for name in expected_names), (
        "fixed_ansatz_results.json contains an unverified solution -- refusing to build CDR on top of it"
    )
    assert ansatz_results["two_qubit_gate_count"] == 11, (
        f"expected the v2 (11-gate) ansatz, found two_qubit_gate_count="
        f"{ansatz_results['two_qubit_gate_count']} -- run fixed_ansatz.py again"
    )

    with open(NATIVE_ZNE_RESULTS_PATH) as f:
        zne_results = json.load(f)
    quadratic_zne_kcal = zne_results["zne_report"]["aria-1"]["zne_quadratic_kcal"]
    classical_floor_kcal = zne_results["classical_floor_kcal"]
    chemical_accuracy_kcal = zne_results["chemical_accuracy_kcal"]

    return {
        "terms": terms, "alpha_labels": alpha_labels, "identity_label": identity_label,
        "signs": signs, "lambdas": lambdas, "u_vecs": u_top, "enuc": enuc,
        "e_mixed": e_mixed, "noiseless_numpy": noiseless_numpy, "exact_energy": exact_energy,
        "alpha_cache": alpha_cache, "solutions": solutions,
        "quadratic_zne_kcal": quadratic_zne_kcal,
        "classical_floor_kcal": classical_floor_kcal, "chemical_accuracy_kcal": chemical_accuracy_kcal,
        "gate_count": ansatz_results["two_qubit_gate_count"],
        "gate_count_v1": ansatz_results["two_qubit_gate_count_v1"],
        "baseline_native_14": ansatz_results["baseline_native_14"],
        "baseline_cx_11": ansatz_results["baseline_cx_11"],
        "native_gate_accounting": ansatz_results["native_gate_accounting"],
    }


def derive_beta_matrices(alpha_matrices, signs, K):
    S = np.diag(signs)
    return {label: S @ mat @ S for label, mat in alpha_matrices.items()}


# ---------------------------------------------------------------------------
# Pure-numpy injection sanity test (unchanged from the preceding version):
# confirms injecting a KNOWN uniform scale and dividing it back out is an
# exact algebraic no-op on the forged energy -- isolates "does the
# correction math work" from "can CDR LEARN the scale".
# ---------------------------------------------------------------------------

def injection_sanity_test(p):
    terms, lambdas = p["terms"], p["lambdas"]
    alpha_cache, identity_label, signs, enuc = p["alpha_cache"], p["identity_label"], p["signs"], p["enuc"]
    out = {}
    for f in (0.90, 0.95):
        alpha_scaled = {}
        for label, mat in alpha_cache.items():
            m = mat[:K, :K].copy()
            if label != identity_label:
                m = (m * f) / f
            alpha_scaled[label] = m
        beta_scaled = derive_beta_matrices(alpha_scaled, signs, K)
        E = effrag.ef_energy_from_noisy_matrices(terms, lambdas, alpha_scaled, beta_scaled, enuc, K)
        out[f] = E
    return out


# ---------------------------------------------------------------------------
# Local noise model: depolarizing_error(param, n) applied per gate. The
# qiskit-aer convention verified empirically (not assumed): a single
# application shrinks any non-identity n-qubit Pauli expectation by EXACTLY
# a factor of (1-param) -- confirmed via a direct Aer probe (1 cx gate,
# param=0.02 -> measured shrink 0.980..0.980 across XX/ZZ/YY; 1 u3 gate,
# param=0.02 -> measured shrink 0.980 on X/Z) before trusting it here. An
# earlier hand-derived Pauli-twirl formula (shrink=1-(16/15)*param) was
# WRONG and is not used.
# ---------------------------------------------------------------------------

def build_noise_model():
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(P2_PER_GATE, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(P1_PER_GATE, 1), "u3")
    return nm


def make_noisy_estimator():
    nm = build_noise_model()
    return AerEstimatorV2(options={"backend_options": {"noise_model": nm, "method": "density_matrix"}})


def exact_labels(angles, labels):
    sv = Statevector.from_instruction(build_ansatz(angles))
    return {l: float(sv.expectation_value(Pauli(l)).real) for l in labels}


def noisy_labels(angles, labels, estimator):
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    obs = [Pauli(l) for l in labels]
    result = estimator.run([(qc, obs)]).result()
    evs = np.atleast_1d(result[0].data.evs)
    return {l: float(np.real(v)) for l, v in zip(labels, evs)}


def abstract_circuit_gate_counts():
    """The gate counts the LOCAL simulator actually sees (u3/cx-transpiled
    ansatz) -- reported for honesty, since these differ from both the
    abstract-CX count (11) and the native-optimized counts in
    fixed_ansatz.py's own report (which vary per target)."""
    qc = transpile(build_ansatz([0.1, 0.2, 0.3, 0.4, 0.5]), basis_gates=BASIS_GATES, optimization_level=0)
    counts = qc.count_ops()
    return {"n2q": counts.get("cx", 0), "n1q": counts.get("u3", 0)}


# ---------------------------------------------------------------------------
# The 25 named target slots (5 diagonal Schmidt states + 20 phase states,
# 2 per pair via the real gauge) exactly match fixed_ansatz_results.json's
# "solutions" keys.
# ---------------------------------------------------------------------------

def slot_names(K):
    names = [f"u_{n}" for n in range(K)]
    for n in range(K):
        for m in range(K):
            if n < m:
                names += [f"(u{n}+u{m})", f"(u{n}-u{m})"]
    return names


def measure_raw_per_slot(solutions, non_id_labels, measure_fn):
    return {name: measure_fn(solutions[name]["angles"], non_id_labels) for name in solutions}


def combine_matrices(raw, alpha_labels, identity_label, K, per_slot_scale=None, per_label_scale=None):
    """Build the K x K alpha matrices from 25 raw per-slot label values,
    optionally dividing each raw value by a per-slot and/or per-label scale
    BEFORE the (E0-E2)/2 cross-term reconstruction. Identity label is NEVER
    scaled -- always exactly 1 on the diagonal, 0 off-diagonal."""
    def corrected(name, label):
        v = raw[name][label]
        if per_slot_scale is not None:
            v = v / per_slot_scale.get(name, 1.0)
        if per_label_scale is not None:
            v = v / per_label_scale.get(label, 1.0)
        return v

    mats = {l: np.zeros((K, K)) for l in alpha_labels}
    for n in range(K):
        name = f"u_{n}"
        for l in alpha_labels:
            mats[l][n, n] = 1.0 if l == identity_label else corrected(name, l)
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            name0, name2 = f"(u{n}+u{m})", f"(u{n}-u{m})"
            for l in alpha_labels:
                if l == identity_label:
                    re = 0.0
                else:
                    re = (corrected(name0, l) - corrected(name2, l)) / 2
                mats[l][n, m] = re
                mats[l][m, n] = re
    return mats


def energy_from_alpha_matrices(alpha_mats, p):
    beta_mats = derive_beta_matrices(alpha_mats, p["signs"], K)
    E = effrag.ef_energy_from_noisy_matrices(p["terms"], p["lambdas"], alpha_mats, beta_mats, p["enuc"], K)
    err_kcal = abs(E - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    f_signal = (E - p["e_mixed"]) / (p["noiseless_numpy"] - p["e_mixed"])
    return E, err_kcal, f_signal


# ---------------------------------------------------------------------------
# Training data + scale fitting
# ---------------------------------------------------------------------------

def generate_training_data(non_id_labels, noisy_estimator, n_per_slot, seed):
    rng = np.random.default_rng(seed)
    names = slot_names(K)
    rows = []
    for name in names:
        for _ in range(n_per_slot):
            angles = rng.uniform(-np.pi, np.pi, 5)
            exact_vals = exact_labels(angles, non_id_labels)
            noisy_vals = noisy_labels(angles, non_id_labels, noisy_estimator)
            for l in non_id_labels:
                rows.append({"slot": name, "label": l, "exact": exact_vals[l], "noisy": noisy_vals[l]})
    return rows


def fit_scale(pairs):
    """Weighted least squares through the origin: f = sum(exact*noisy)/sum(exact^2)."""
    if len(pairs) < 3:
        return None
    exact = np.array([e for e, _ in pairs])
    noisy = np.array([n for _, n in pairs])
    denom = float(np.sum(exact * exact))
    if denom < 1e-12:
        return None
    return float(np.sum(exact * noisy) / denom)


def filtered_pairs(training, label=None, slot=None):
    out = []
    for row in training:
        if label is not None and row["label"] != label:
            continue
        if slot is not None and row["slot"] != slot:
            continue
        if abs(row["exact"]) < LOW_SIGNAL_CUTOFF:
            continue
        out.append((row["exact"], row["noisy"]))
    return out


def fit_all_scales(training, non_id_labels):
    global_scale = fit_scale(filtered_pairs(training))

    per_basis, per_basis_fallback = {}, []
    for l in non_id_labels:
        f = fit_scale(filtered_pairs(training, label=l))
        if f is None:
            f = global_scale
            per_basis_fallback.append(l)
        per_basis[l] = f

    per_circuit, per_circuit_fallback = {}, []
    for name in slot_names(K):
        f = fit_scale(filtered_pairs(training, slot=name))
        if f is None:
            f = global_scale
            per_circuit_fallback.append(name)
        per_circuit[name] = f

    return {
        "global_scale": global_scale,
        "per_basis_scale": per_basis, "per_basis_fallback_to_global": per_basis_fallback,
        "per_circuit_scale": per_circuit, "per_circuit_fallback_to_global": per_circuit_fallback,
    }


def run_one_seed(p, non_id_labels, raw_noisy, noisy_estimator, seed):
    training = generate_training_data(non_id_labels, noisy_estimator, N_TRAIN_PER_SLOT, seed)
    n_dropped = sum(1 for row in training if abs(row["exact"]) < LOW_SIGNAL_CUTOFF)
    scales = fit_all_scales(training, non_id_labels)

    global_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                    per_slot_scale=None if scales["global_scale"] is None else
                                    {n: scales["global_scale"] for n in slot_names(K)})
    E_global, err_global_kcal, f_global = energy_from_alpha_matrices(global_mats, p)

    basis_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                   per_label_scale=scales["per_basis_scale"])
    E_basis, err_basis_kcal, f_basis = energy_from_alpha_matrices(basis_mats, p)

    circuit_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                     per_slot_scale=scales["per_circuit_scale"])
    E_circuit, err_circuit_kcal, f_circuit = energy_from_alpha_matrices(circuit_mats, p)

    return {
        "seed": seed, "n_training_pairs_dropped_low_signal": n_dropped,
        "global_scale": scales["global_scale"],
        "per_basis_scale_mean": float(np.mean(list(scales["per_basis_scale"].values()))),
        "per_circuit_scale_mean": float(np.mean(list(scales["per_circuit_scale"].values()))),
        "err_global_kcal": err_global_kcal, "err_basis_kcal": err_basis_kcal, "err_circuit_kcal": err_circuit_kcal,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 70)
    print("  CDR noise-scale learning, v2 ansatz (11 2q gates), seed-swept (simulator only)")
    print("=" * 70)

    p = setup()
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    print(f"\n  {len(p['alpha_labels'])} alpha labels ({len(non_id_labels)} non-identity)")
    print(f"  exact fragment energy       = {p['exact_energy']:.6f} Ha")
    print(f"  exact numpy forged (K={K})    = {p['noiseless_numpy']:.6f} Ha")
    print(f"  ansatz 2-qubit gate count    = {p['gate_count']} "
          f"(v1 was {p['gate_count_v1']}; baselines: {p['baseline_native_14']} native, {p['baseline_cx_11']} CX)")

    nga = p["native_gate_accounting"]
    f_before_opt = fidelity_from_counts(nga["n2q_before_optimizer"], nga["n1q_before_optimizer"])
    print(f"\n  native gate accounting (aria-1/ms, from fixed_ansatz.py):")
    print(f"    before TrappedIonOptimizerPlugin: {nga['n2q_before_optimizer']} 2q, "
          f"{nga['n1q_before_optimizer']} 1q (FIXED across all 25 targets) -> f={f_before_opt:.4f}")
    print(f"    after  TrappedIonOptimizerPlugin: {nga['n2q_after_optimizer']} 2q, "
          f"{nga['n1q_after_optimizer']} 1q (this ONE representative target; varies by target -- "
          f"see fixed_ansatz_results.json/native_gate_accounting, "
          f"counts_consistent_across_angle_sets={nga['counts_consistent_across_angle_sets']})")
    print(f"    f = (1-{P2_PER_GATE})^n2q * (1-{P1_PER_GATE:.6f})^n1q, p2 given, p1=p2/40")

    abstract_counts = abstract_circuit_gate_counts()
    f_abstract_single_gate_estimate = fidelity_from_counts(abstract_counts["n2q"], abstract_counts["n1q"])
    print(f"\n  LOCAL SIMULATOR circuit (u3/cx-transpiled, what CDR actually simulates): "
          f"{abstract_counts['n2q']} cx, {abstract_counts['n1q']} u3 "
          f"-> naive single-gate-type f-estimate {f_abstract_single_gate_estimate:.4f} "
          "(actual simulated f measured below, not assumed from this count alone since different "
          "labels sit behind different numbers of upstream gates)")

    inj = injection_sanity_test(p)
    print("\n  Injection sanity test (inject known uniform f, divide it back out --")
    print("  should be an exact algebraic no-op):")
    for f, E in inj.items():
        print(f"    f={f}: recovered E = {E:.6f} Ha (diff {abs(E - p['noiseless_numpy']):.2e})")
    injection_ok = all(abs(E - p["noiseless_numpy"]) < 1e-9 for E in inj.values())
    print(f"  {'PASS' if injection_ok else 'FAIL'}")

    print(f"\n  Building local noise model: depolarizing_error(cx)={P2_PER_GATE}, "
          f"depolarizing_error(u3)={P1_PER_GATE:.6f} (P2_PER_GATE/40) -- applied directly, no calibration/root-find.")
    noisy_estimator = make_noisy_estimator()

    print("  Measuring the 25 real target circuits under this noise model...")
    raw_noisy = measure_raw_per_slot(p["solutions"], non_id_labels, lambda a, l: noisy_labels(a, l, noisy_estimator))
    raw_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K)
    E_raw, err_raw_kcal, f_raw = energy_from_alpha_matrices(raw_mats, p)
    print(f"  raw (uncorrected): E={E_raw:.6f} Ha, error={err_raw_kcal:.3f} kcal/mol, "
          f"measured full-forged-energy f={f_raw:.6f}")

    print(f"\n  Seed sweep: {N_SEEDS} seeds, {len(slot_names(K))} slots x {N_TRAIN_PER_SLOT} "
          f"random-angle training circuits each...")
    seed_results = []
    for seed in SEEDS:
        r = run_one_seed(p, non_id_labels, raw_noisy, noisy_estimator, seed)
        seed_results.append(r)
        print(f"    seed={seed}: global={r['err_global_kcal']:.3f}  per_basis={r['err_basis_kcal']:.3f}  "
              f"per_circuit={r['err_circuit_kcal']:.3f} kcal/mol")

    def stats(key):
        vals = [r[key] for r in seed_results]
        return float(np.mean(vals)), float(np.std(vals)), float(np.min(vals)), float(np.max(vals))

    g_mean, g_std, g_min, g_max = stats("err_global_kcal")
    b_mean, b_std, b_min, b_max = stats("err_basis_kcal")
    c_mean, c_std, c_min, c_max = stats("err_circuit_kcal")

    def n_reached(key, threshold):
        return sum(1 for r in seed_results if r[key] < threshold)

    chem_acc = p["chemical_accuracy_kcal"]
    g_reached = n_reached("err_global_kcal", chem_acc)
    b_reached = n_reached("err_basis_kcal", chem_acc)
    c_reached = n_reached("err_circuit_kcal", chem_acc)

    print("\n" + "-" * 78)
    print(f"  {'scheme':<14}{'mean':>10}{'std':>10}{'min':>10}{'max':>10}{'reach chem.acc.':>20}")
    print(f"  {'raw':<14}{err_raw_kcal:>10.3f}{'--':>10}{'--':>10}{'--':>10}{'n/a':>20}")
    print(f"  {'global':<14}{g_mean:>10.3f}{g_std:>10.3f}{g_min:>10.3f}{g_max:>10.3f}{f'{g_reached}/{N_SEEDS} seeds':>20}")
    print(f"  {'per-basis':<14}{b_mean:>10.3f}{b_std:>10.3f}{b_min:>10.3f}{b_max:>10.3f}{f'{b_reached}/{N_SEEDS} seeds':>20}")
    print(f"  {'per-circuit':<14}{c_mean:>10.3f}{c_std:>10.3f}{c_min:>10.3f}{c_max:>10.3f}{f'{c_reached}/{N_SEEDS} seeds':>20}")
    print("-" * 78)
    print(f"  (kcal/mol; gain over raw = mean: global {err_raw_kcal / g_mean:.2f}x, "
          f"per-basis {err_raw_kcal / b_mean:.2f}x, per-circuit {err_raw_kcal / c_mean:.2f}x)")

    print(f"\n  reference: quadratic ZNE (real aria-1 hardware) = {p['quadratic_zne_kcal']:.2f} kcal/mol")
    print(f"  reference: classical floor                      = {p['classical_floor_kcal']:.4f} kcal/mol")
    print(f"  reference: chemical accuracy target              = {chem_acc:.2f} kcal/mol")
    print(f"  context: prior StatePreparation-based CDR attempt: "
          f"{PRIOR_STATEPREP_RAW_KCAL} -> {PRIOR_STATEPREP_GLOBAL_KCAL} kcal/mol "
          f"({PRIOR_STATEPREP_GAIN:.2f}x, NOT reproduced here -- that code no longer exists)")

    best_scheme, best_mean, best_reached = max(
        [("global", g_mean, g_reached), ("per_basis", b_mean, b_reached), ("per_circuit", c_mean, c_reached)],
        key=lambda t: -t[1],
    )
    print("\n  HONEST DIAGNOSIS:")
    print(f"  Best scheme by mean error: {best_scheme} ({best_mean:.3f} kcal/mol mean over {N_SEEDS} seeds).")
    if best_reached >= N_SEEDS - 1:
        print(f"  Chemical accuracy ({chem_acc} kcal/mol) is reached RELIABLY: {best_reached}/{N_SEEDS} seeds.")
    elif best_reached >= N_SEEDS // 2:
        print(f"  Chemical accuracy is reached on a MAJORITY but not all seeds: {best_reached}/{N_SEEDS}. "
              "Treat as 'usually gets there', not guaranteed.")
    else:
        print(f"  Chemical accuracy is NOT reached reliably: only {best_reached}/{N_SEEDS} seeds. "
              "A single lucky seed is not evidence CDR reliably hits chemical accuracy here -- "
              "the mean and std above are the honest summary, not the best individual seed.")

    results = {
        "gate_count": p["gate_count"], "gate_count_v1": p["gate_count_v1"],
        "baseline_native_14": p["baseline_native_14"], "baseline_cx_11": p["baseline_cx_11"],
        "native_gate_accounting": nga,
        "f_before_optimizer_fixed": f_before_opt,
        "p2_per_gate": P2_PER_GATE, "p1_per_gate": P1_PER_GATE,
        "abstract_local_sim_gate_counts": abstract_counts,
        "injection_sanity_test": {str(f): E for f, E in inj.items()},
        "injection_sanity_test_pass": injection_ok,
        "raw_measured_f": f_raw,
        "n_seeds": N_SEEDS, "n_train_per_slot": N_TRAIN_PER_SLOT, "low_signal_cutoff": LOW_SIGNAL_CUTOFF,
        "seed_results": seed_results,
        "results_kcal_mol": {
            "raw": round(err_raw_kcal, 4),
            "global": {"mean": round(g_mean, 4), "std": round(g_std, 4), "min": round(g_min, 4), "max": round(g_max, 4),
                       "reached_chem_acc": f"{g_reached}/{N_SEEDS}"},
            "per_basis": {"mean": round(b_mean, 4), "std": round(b_std, 4), "min": round(b_min, 4), "max": round(b_max, 4),
                          "reached_chem_acc": f"{b_reached}/{N_SEEDS}"},
            "per_circuit": {"mean": round(c_mean, 4), "std": round(c_std, 4), "min": round(c_min, 4), "max": round(c_max, 4),
                            "reached_chem_acc": f"{c_reached}/{N_SEEDS}"},
        },
        "gain_over_raw_mean": {
            "global": round(err_raw_kcal / g_mean, 4),
            "per_basis": round(err_raw_kcal / b_mean, 4),
            "per_circuit": round(err_raw_kcal / c_mean, 4),
        },
        "best_scheme_by_mean": best_scheme,
        "reaches_chemical_accuracy_reliably": bool(best_reached >= N_SEEDS - 1),
        "references": {
            "quadratic_zne_kcal_real_aria1": p["quadratic_zne_kcal"],
            "classical_floor_kcal": p["classical_floor_kcal"],
            "chemical_accuracy_kcal": chem_acc,
            "prior_stateprep_cdr_attempt_raw_kcal": PRIOR_STATEPREP_RAW_KCAL,
            "prior_stateprep_cdr_attempt_global_kcal": PRIOR_STATEPREP_GLOBAL_KCAL,
            "prior_stateprep_cdr_attempt_gain": round(PRIOR_STATEPREP_GAIN, 4),
            "prior_stateprep_cdr_attempt_note": "context from an earlier attempt in this project; that code no longer exists and is not reproduced here",
        },
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
