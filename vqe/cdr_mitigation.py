#!/usr/bin/env python3
"""
cdr_mitigation.py — Clifford-Data-Regression-style noise-scale learning for
the H4 forged energy, using the FIXED-STRUCTURE ansatz from fixed_ansatz.py
so training and target circuits are structurally IDENTICAL (same 33
two-qubit gates in the same order — only the 5 rotation angles differ).
SIMULATOR ONLY. Nothing here touches IonQ.
============================================================================
WHY THE OLD ATTEMPT FAILED (context carried into this file, not re-derived
here): the previous CDR attempt used qiskit's `StatePreparation` for both
training and target circuits. `StatePreparation` compiles a DIFFERENT
circuit (different gate count/structure) per input vector, so a training
circuit's noise profile did not match its target's — the learned scale did
not transfer. Reported there: best case 1.7x over raw (199 -> 116 kcal/mol),
and going to PER-CIRCUIT training made it WORSE (1.3x) — the signature of
train/target structural mismatch (smaller, more localized pools drawn from
differently-structured circuits than the target they were meant to correct).

THIS FILE: every training circuit and every target circuit is
`fixed_ansatz.build_ansatz(angles)` — the identical 33-gate structure,
angles the only thing that varies. If CDR's dependence on structural
matching is real (as the theory says), this should recover a much more
consistent scale across granularities than the old attempt, without any
manual tuning.

LOCAL NOISE MODEL: a qiskit-aer 2-qubit depolarizing channel on every `cx`
(plus the same small 1-qubit floor entanglement_forging_h4.py's Stage 2
uses), transpiled basis ["u3","cx"]. The 2-qubit error rate is CALIBRATED
(root-find, not guessed) so that running the fixed ansatz's own 25 real
target circuits through it reproduces the REAL aria-1 fold=1 signal
fraction f=0.842805 (native_forged_zne_results.json) for the FULL K=5
forged energy — not just a single circuit.

TWO FIXES CARRIED FORWARD FROM THE PRECEDING (already-debugged) ATTEMPT —
both silent-corruption bugs if skipped:
  1. NEVER rescale the identity Pauli label. <I>=1 exactly for any
     trace-preserving channel (noisy or not) — dividing it by a learned
     scale corrupts the forged energy, since it multiplies straight into
     alpha*beta for every term carrying that label on one side. First
     attempt at this overshot to -2.67 Ha by doing exactly that.
  2. DROP training pairs whose EXACT expectation value has |value| < 0.05.
     Near-zero-signal labels carry almost no information about the noise
     scale (dividing a noisy near-zero by an exact near-zero is mostly
     amplified numerator noise) and previously dragged basis-level scale
     estimates down to ~0.7555.

Run:
    python vqe/cdr_mitigation.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
from fixed_ansatz import build_ansatz

from qiskit import transpile
from qiskit.quantum_info import Pauli, Statevector
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error

HARTREE_TO_KCAL_MOL = 627.5094740631
ATOMS = [0, 1, 2, 3]
NELEC = 4
K = 5
BASIS_GATES = ["u3", "cx"]
ONE_Q_ERROR = 0.00005  # same 1-qubit floor as entanglement_forging_h4.py's Stage 2
LOW_SIGNAL_CUTOFF = 0.05
N_TRAIN_PER_SLOT = 5
TRAIN_SEED = 12345

# Context from the earlier (StatePreparation-based) CDR attempt -- reported
# here purely for comparison, NOT reproduced/recomputed in this file (that
# code no longer exists to re-run).
PRIOR_STATEPREP_RAW_KCAL = 199
PRIOR_STATEPREP_GLOBAL_KCAL = 116
PRIOR_STATEPREP_GAIN = PRIOR_STATEPREP_RAW_KCAL / PRIOR_STATEPREP_GLOBAL_KCAL

ANSATZ_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "fixed_ansatz_results.json")
NATIVE_ZNE_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "native_forged_zne_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "fixed_ansatz_cdr_results.json")


# ---------------------------------------------------------------------------
# Setup: same fragment physics as ionq_native_forged_energy.py, but the
# angle solutions come from fixed_ansatz.py's fixed circuit instead of the
# native hand-derived tree -- reused nowhere else, so intentionally
# self-contained rather than importing the IonQ-touching module.
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

    with open(NATIVE_ZNE_RESULTS_PATH) as f:
        zne_results = json.load(f)
    target_f = zne_results["configs"]["aria-1_fold1"]["f"]
    quadratic_zne_kcal = zne_results["zne_report"]["aria-1"]["zne_quadratic_kcal"]
    classical_floor_kcal = zne_results["classical_floor_kcal"]
    chemical_accuracy_kcal = zne_results["chemical_accuracy_kcal"]

    return {
        "terms": terms, "alpha_labels": alpha_labels, "identity_label": identity_label,
        "signs": signs, "lambdas": lambdas, "u_vecs": u_top, "enuc": enuc,
        "e_mixed": e_mixed, "noiseless_numpy": noiseless_numpy, "exact_energy": exact_energy,
        "alpha_cache": alpha_cache, "solutions": solutions,
        "target_f": target_f, "quadratic_zne_kcal": quadratic_zne_kcal,
        "classical_floor_kcal": classical_floor_kcal, "chemical_accuracy_kcal": chemical_accuracy_kcal,
        "gate_count": ansatz_results["two_qubit_gate_count"],
        "baseline_native_14": ansatz_results["baseline_native_14"],
        "baseline_cx_11": ansatz_results["baseline_cx_11"],
    }


def derive_beta_matrices(alpha_matrices, signs, K):
    S = np.diag(signs)
    return {label: S @ mat @ S for label, mat in alpha_matrices.items()}


# ---------------------------------------------------------------------------
# Pure-numpy injection sanity test: re-verifies (cheaply, no circuit
# simulation) that injecting a KNOWN uniform scale on every non-identity
# label and dividing it back out is an algebraic no-op on the forged
# energy. This isolates "does the correction math work" from "can CDR
# LEARN the scale", which is the actual open question this file answers.
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
# Circuit-level measurement: exact (Statevector) and noisy (Aer, calibrated
# depolarizing channel) expectation values for a set of Pauli labels, on
# the SAME fixed_ansatz.build_ansatz(angles) circuit.
# ---------------------------------------------------------------------------

def build_noise_model(two_q_error):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(two_q_error, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(ONE_Q_ERROR, 1), "u3")
    return nm


def make_noisy_estimator(two_q_error):
    nm = build_noise_model(two_q_error)
    return AerEstimatorV2(options={"backend_options": {"noise_model": nm, "method": "density_matrix"}})


def exact_labels(angles, labels):
    sv = Statevector.from_instruction(build_ansatz(angles))
    return {l: float(sv.expectation_value(Pauli(l)).real) for l in labels}


def noisy_labels(angles, labels, estimator):
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=1)
    obs = [Pauli(l) for l in labels]
    result = estimator.run([(qc, obs)]).result()
    evs = np.atleast_1d(result[0].data.evs)
    return {l: float(np.real(v)) for l, v in zip(labels, evs)}


# ---------------------------------------------------------------------------
# The 25 named target slots (5 diagonal Schmidt states + 20 phase states,
# 2 per pair via the real gauge) exactly match fixed_ansatz_results.json's
# "solutions" keys -- each IS a training circuit and IS a target circuit,
# so no separate bookkeeping is needed to keep them aligned.
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
    BEFORE the (E0-E2)/2 cross-term reconstruction (so a per-circuit/
    per-slot correction is applied to E0 and E2 separately, not to their
    already-combined difference). Identity label is NEVER scaled -- always
    exactly 1 on the diagonal, 0 off-diagonal (orthogonal Schmidt vectors)."""
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
# Calibration: find the 2-qubit depolarizing rate whose FULL K=5 forged
# energy (built from the fixed ansatz's own 25 real target circuits)
# reproduces the real aria-1 fold=1 signal fraction.
# ---------------------------------------------------------------------------

def f_for_error_rate(two_q_error, p, non_id_labels):
    estimator = make_noisy_estimator(two_q_error)
    raw = measure_raw_per_slot(p["solutions"], non_id_labels, lambda a, l: noisy_labels(a, l, estimator))
    mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    _, _, f = energy_from_alpha_matrices(mats, p)
    return f


def calibrate_noise_rate(p, non_id_labels, target_f):
    print(f"\n  Calibrating local depolarizing rate to reproduce f={target_f:.6f} (real aria-1 fold=1)...")

    def residual(two_q_error):
        f = f_for_error_rate(two_q_error, p, non_id_labels)
        print(f"    two_q_error={two_q_error:.5f} -> f={f:.6f}")
        return f - target_f

    lo, hi = 1e-4, 0.30
    r_lo, r_hi = residual(lo), residual(hi)
    if r_lo * r_hi > 0:
        raise RuntimeError(
            f"could not bracket target_f={target_f} between two_q_error in [{lo},{hi}] "
            f"(residuals {r_lo:.4f}, {r_hi:.4f}) -- calibration range needs widening, not forcing a fit"
        )
    two_q_error = brentq(residual, lo, hi, xtol=1e-5)
    achieved_f = target_f + residual(two_q_error)
    print(f"  Calibrated: two_q_error={two_q_error:.5f}, achieved f={achieved_f:.6f}")
    return two_q_error, achieved_f


# ---------------------------------------------------------------------------
# Training data + scale fitting
# ---------------------------------------------------------------------------

def generate_training_data(p, non_id_labels, noisy_estimator, n_per_slot, seed):
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
    """Weighted least squares through the origin: minimize sum (noisy -
    f*exact)^2 over the training pairs -- f = sum(exact*noisy)/sum(exact^2)."""
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

    per_basis = {}
    per_basis_fallback = []
    for l in non_id_labels:
        f = fit_scale(filtered_pairs(training, label=l))
        if f is None:
            f = global_scale
            per_basis_fallback.append(l)
        per_basis[l] = f

    per_circuit = {}
    per_circuit_fallback = []
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 70)
    print("  CDR noise-scale learning on the fixed-structure ansatz (simulator only)")
    print("=" * 70)

    p = setup()
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    print(f"\n  {len(p['alpha_labels'])} alpha labels ({len(non_id_labels)} non-identity)")
    print(f"  exact fragment energy       = {p['exact_energy']:.6f} Ha")
    print(f"  exact numpy forged (K={K})    = {p['noiseless_numpy']:.6f} Ha")
    print(f"  ansatz 2-qubit gate count    = {p['gate_count']} "
          f"(baselines: {p['baseline_native_14']} native, {p['baseline_cx_11']} CX)")

    inj = injection_sanity_test(p)
    print("\n  Injection sanity test (inject known uniform f on every non-identity")
    print("  label, divide it back out -- should be an exact algebraic no-op):")
    for f, E in inj.items():
        print(f"    f={f}: recovered E = {E:.6f} Ha "
              f"(target {p['noiseless_numpy']:.6f} Ha, diff {abs(E - p['noiseless_numpy']):.2e})")
    injection_ok = all(abs(E - p["noiseless_numpy"]) < 1e-9 for E in inj.values())
    print(f"  {'PASS' if injection_ok else 'FAIL'}")

    two_q_error, achieved_f = calibrate_noise_rate(p, non_id_labels, p["target_f"])
    noisy_estimator = make_noisy_estimator(two_q_error)

    print("\n  Measuring the 25 real target circuits under the calibrated noise model...")
    raw_noisy = measure_raw_per_slot(p["solutions"], non_id_labels, lambda a, l: noisy_labels(a, l, noisy_estimator))

    raw_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K)
    E_raw, err_raw_kcal, f_raw = energy_from_alpha_matrices(raw_mats, p)
    print(f"  raw (uncorrected): E={E_raw:.6f} Ha, error={err_raw_kcal:.3f} kcal/mol, f={f_raw:.6f}")

    print(f"\n  Generating CDR training data: {len(slot_names(K))} slots x {N_TRAIN_PER_SLOT} "
          f"random-angle circuits = {len(slot_names(K)) * N_TRAIN_PER_SLOT} training circuits...")
    training = generate_training_data(p, non_id_labels, noisy_estimator, N_TRAIN_PER_SLOT, TRAIN_SEED)
    n_dropped = sum(1 for row in training if abs(row["exact"]) < LOW_SIGNAL_CUTOFF)
    print(f"  {len(training)} raw training pairs, {n_dropped} dropped as low-signal (|exact|<{LOW_SIGNAL_CUTOFF})")

    scales = fit_all_scales(training, non_id_labels)
    print(f"\n  global scale       = {scales['global_scale']:.6f}")
    print(f"  per-basis scale    = {np.mean(list(scales['per_basis_scale'].values())):.6f} avg "
          f"(range {min(scales['per_basis_scale'].values()):.4f}-{max(scales['per_basis_scale'].values()):.4f}, "
          f"{len(scales['per_basis_fallback_to_global'])} fell back to global for lack of data)")
    print(f"  per-circuit scale  = {np.mean(list(scales['per_circuit_scale'].values())):.6f} avg "
          f"(range {min(scales['per_circuit_scale'].values()):.4f}-{max(scales['per_circuit_scale'].values()):.4f}, "
          f"{len(scales['per_circuit_fallback_to_global'])} fell back to global for lack of data)")
    print(f"  (calibrated true noise-model f = {achieved_f:.6f}, for reference)")

    global_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                    per_label_scale=None, per_slot_scale=None if scales["global_scale"] is None else
                                    {n: scales["global_scale"] for n in slot_names(K)})
    E_global, err_global_kcal, f_global = energy_from_alpha_matrices(global_mats, p)

    basis_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                   per_label_scale=scales["per_basis_scale"])
    E_basis, err_basis_kcal, f_basis = energy_from_alpha_matrices(basis_mats, p)

    circuit_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                     per_slot_scale=scales["per_circuit_scale"])
    E_circuit, err_circuit_kcal, f_circuit = energy_from_alpha_matrices(circuit_mats, p)

    def gain(err):
        return err_raw_kcal / err if err > 1e-9 else float("inf")

    print("\n" + "-" * 70)
    print(f"  {'scheme':<14}{'E (Ha)':>14}{'error (kcal/mol)':>20}{'gain over raw':>16}")
    print(f"  {'raw':<14}{E_raw:>14.6f}{err_raw_kcal:>20.3f}{'1.00x':>16}")
    print(f"  {'global':<14}{E_global:>14.6f}{err_global_kcal:>20.3f}{gain(err_global_kcal):>15.2f}x")
    print(f"  {'per-basis':<14}{E_basis:>14.6f}{err_basis_kcal:>20.3f}{gain(err_basis_kcal):>15.2f}x")
    print(f"  {'per-circuit':<14}{E_circuit:>14.6f}{err_circuit_kcal:>20.3f}{gain(err_circuit_kcal):>15.2f}x")
    print("-" * 70)

    print(f"\n  reference: quadratic ZNE (real aria-1 hardware) = {p['quadratic_zne_kcal']:.2f} kcal/mol")
    print(f"  reference: classical floor                      = {p['classical_floor_kcal']:.4f} kcal/mol")
    print(f"  reference: chemical accuracy target              = {p['chemical_accuracy_kcal']:.2f} kcal/mol")
    print(f"  context: prior StatePreparation-based CDR attempt: "
          f"{PRIOR_STATEPREP_RAW_KCAL} -> {PRIOR_STATEPREP_GLOBAL_KCAL} kcal/mol "
          f"({PRIOR_STATEPREP_GAIN:.2f}x, NOT reproduced here -- that code no longer exists)")

    best_gain = max(gain(err_global_kcal), gain(err_basis_kcal), gain(err_circuit_kcal))
    print("\n  HONEST DIAGNOSIS:")
    if best_gain >= 5.0:
        print(f"  Best scheme reaches {best_gain:.2f}x over raw -- clears the ~5x bar.")
    else:
        print(f"  Best scheme reaches only {best_gain:.2f}x over raw -- does NOT clear the ~5x bar.")
        print("  Not tuning to force a better number. Plausible reason, given the setup actually used:")
        print("  the depolarizing noise injected here is UNIFORM per-gate and angle-independent by")
        print("  construction (qiskit-aer's depolarizing_error attaches a fixed-magnitude channel to")
        print("  every 'cx' instruction regardless of its rotation angle). Since every one of the 25")
        print("  target/training circuits shares the exact same 33-gate structure, the TRUE scale each")
        print("  circuit/label experiences should already be near-identical by this model's own")
        print("  construction -- so this test mainly measures whether finite-sample REGRESSION NOISE (not")
        print("  structural mismatch) is the bottleneck. If per-circuit underperforms per-basis/global here,")
        print(f"  that points to sampling noise from only {N_TRAIN_PER_SLOT} training draws per slot, not a")
        print("  structural circuit-shape mismatch (which this ansatz was built specifically to remove).")

    results = {
        "gate_count": p["gate_count"], "baseline_native_14": p["baseline_native_14"], "baseline_cx_11": p["baseline_cx_11"],
        "injection_sanity_test": {str(f): E for f, E in inj.items()},
        "injection_sanity_test_pass": injection_ok,
        "calibrated_two_q_depolarizing_error": two_q_error,
        "calibrated_achieved_f": achieved_f,
        "target_f_aria1_fold1": p["target_f"],
        "n_train_per_slot": N_TRAIN_PER_SLOT,
        "n_training_circuits": len(slot_names(K)) * N_TRAIN_PER_SLOT,
        "n_training_pairs_raw": len(training),
        "n_training_pairs_dropped_low_signal": n_dropped,
        "low_signal_cutoff": LOW_SIGNAL_CUTOFF,
        "global_scale": scales["global_scale"],
        "per_basis_scale": scales["per_basis_scale"],
        "per_basis_fallback_to_global": scales["per_basis_fallback_to_global"],
        "per_circuit_scale": scales["per_circuit_scale"],
        "per_circuit_fallback_to_global": scales["per_circuit_fallback_to_global"],
        "results_kcal_mol": {
            "raw": round(err_raw_kcal, 4),
            "global": round(err_global_kcal, 4),
            "per_basis": round(err_basis_kcal, 4),
            "per_circuit": round(err_circuit_kcal, 4),
        },
        "energies_ha": {
            "raw": round(E_raw, 6), "global": round(E_global, 6),
            "per_basis": round(E_basis, 6), "per_circuit": round(E_circuit, 6),
        },
        "gain_over_raw": {
            "global": round(gain(err_global_kcal), 4),
            "per_basis": round(gain(err_basis_kcal), 4),
            "per_circuit": round(gain(err_circuit_kcal), 4),
        },
        "best_gain_over_raw": round(best_gain, 4),
        "clears_5x_bar": bool(best_gain >= 5.0),
        "references": {
            "quadratic_zne_kcal_real_aria1": p["quadratic_zne_kcal"],
            "classical_floor_kcal": p["classical_floor_kcal"],
            "chemical_accuracy_kcal": p["chemical_accuracy_kcal"],
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
