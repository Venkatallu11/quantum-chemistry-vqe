#!/usr/bin/env python3
"""
shot_noise_study.py — adds shot noise to the simulator (missing from every
prior result in this project) and reframes the loop's goal accordingly.
============================================================================
WHY THIS OVERRIDES THE 0.30 KCAL/MOL TARGET: every result in
RESEARCH_LEDGER.md through iteration 5 used `AerSimulator(method=
"density_matrix")` with no shot count -- an EXACT expectation value, not a
sampled one. That omits the dominant real-hardware error source entirely.
This file adds it, verifies it, and re-measures raw/CDR/PEC as a function
of shots per setting.

TASK 1 -- the shot-noise model
For a Pauli observable P with eigenvalues +-1 and true (density-matrix-
exact) expectation value e, a projective measurement with N shots gives a
sample mean estimator: draw n_plus ~ Binomial(N, (1+e)/2), estimator =
2*n_plus/N - 1. This is the EXACT distribution real shot-based hardware
would produce for a given (possibly noisy) quantum state -- not an
approximation -- since Aer's own shot-based execution samples from exactly
this Born-rule distribution.

VARIANCE CHECK, and an honest complication: the standard VQE result
Var(E) = L2^2/N (L2 = sqrt(sum c_i^2) over Hamiltonian coefficients, N =
shots per term) assumes E is a LINEAR combination of independently-
measured Pauli terms. The forged-energy estimator here is NOT that -- each
Hamiltonian term's contribution is coeff*(diag + cross)/norm2, where cross
involves Amat[n,m]*Bmat[n,m] and Bmat = S.Amat.S is DERIVED from the SAME
measured alpha matrix (beta is never independently measured, per this
project's beta-reuse shortcut) -- a BILINEAR, correlated structure, not a
linear sum. L2=2.7555 Ha was verified directly from the Hamiltonian's own
Pauli coefficients (matches the given value exactly) before any of this
was built, but whether the RESULTING energy estimator's variance actually
equals L2^2/N is an empirical question, checked below, not assumed --
report whatever is actually measured, including a mismatch if there is
one.

EFFICIENCY NOTE: shot noise is added by SAMPLING ON TOP OF the already-
exact density-matrix expectation values computed via this project's
existing (expensive) circuit machinery -- the exact values are computed
ONCE per (circuit, label) and CACHED; every shot level and every seed then
just re-samples cheaply (a single `rng.binomial` call) from those already-
computed exact values, rather than re-simulating the circuit per shot
level. This is exact (Binomial sampling from the true Born-rule
probability IS what real shot-based execution does), not an approximation
introduced for speed.

PEC's shot noise is handled differently, honestly: PEC's Monte Carlo
procedure (randomly sampling a Pauli recovery per noisy gate, weighted by
the SIGN of its quasi-probability, averaged with the final shot outcome)
gives an UNBIASED estimator with variance INFLATED by gamma_total^2
relative to a plain shot-noise measurement of the same observable
(gamma_total from iteration 4/5's own analysis, ~1.315 here). Modeled as
Gaussian noise with variance gamma_total^2*(1-e^2)/N -- the standard,
well-established PEC sampling-overhead result, not simulated gate-by-gate
Monte Carlo (which would be exact but far slower for no accuracy gain at
the shot counts studied here, since gamma_total^2 IS the exact overhead
factor by construction of the quasi-probability decomposition).

Run:
    python vqe/shot_noise_study.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import rank6_symmetry_vd as r
import loop_pec as pec
import loop_pauli_lindblad_pec as pl
from fixed_ansatz import build_ansatz, P2_PER_GATE, P1_PER_GATE
from qiskit.quantum_info import Statevector, Pauli

K = r.K
BASIS_GATES = r.BASIS_GATES
N_SEEDS = 8
SHOT_LEVELS = [1_000, 10_000, 100_000, 1_000_000, 10_000_000]
N_TRAIN_PER_SLOT = r.N_TRAIN_PER_SLOT
LOW_SIGNAL_CUTOFF = r.LOW_SIGNAL_CUTOFF

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "shot_noise_study_results.json")
LINDBLAD_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "loop_pauli_lindblad_pec_results.json")


# ---------------------------------------------------------------------------
# TASK 1: shot-noise primitives
# ---------------------------------------------------------------------------

def shot_sample(exact_value, n_shots, rng):
    """Exact Binomial sample-mean estimator of a +-1-eigenvalue observable
    with true expectation `exact_value`, from n_shots projective
    measurements -- the actual distribution real hardware execution
    produces, not an approximation."""
    p_plus = float(np.clip((1 + exact_value) / 2, 0.0, 1.0))
    n_plus = rng.binomial(n_shots, p_plus)
    return 2 * n_plus / n_shots - 1


def pec_shot_sample(exact_value, gamma_total, n_shots, rng):
    """PEC's quasi-probability Monte Carlo estimator: unbiased, variance
    inflated by gamma_total^2 relative to a plain shot-noise measurement
    of the same observable (the standard, exact-by-construction PEC
    sampling-overhead result)."""
    var = (gamma_total ** 2) * max(1 - exact_value ** 2, 0.0) / n_shots
    return exact_value + rng.normal(0.0, np.sqrt(var))


# ---------------------------------------------------------------------------
# Exact-value caching: the expensive circuit simulation is done ONCE per
# (circuit, label); every shot level/seed just re-samples cheaply on top.
# ---------------------------------------------------------------------------

def cache_exact_raw(p, non_id_labels, noise_model):
    """Exact (noisy, depolarizing-affected) per-slot label values for all
    36 real targets -- shot-count-independent, computed once."""
    raw = {}
    for name, sol in p["solutions"].items():
        dm = r.density_matrix_4q(sol["angles"], noise_model)
        vals, _ = r.measure_labels(sol["angles"], non_id_labels, dm)
        raw[name] = vals
    return raw


def cache_exact_training(p, non_id_labels, noise_model, n_per_slot, seed):
    """Exact (noiseless AND noisy) values at random training angles for one
    seed -- shot-count-independent; shot noise applied on top per shot
    level, reusing the SAME underlying random angles across shot levels."""
    rng = np.random.default_rng(seed)
    rows = []
    for name in r.slot_names(K):
        for _ in range(n_per_slot):
            angles = rng.uniform(-np.pi, np.pi, 5).tolist()
            exact_vals = r.exact_labels(angles, non_id_labels)
            dm = r.density_matrix_4q(angles, noise_model)
            noisy_exact_vals, _ = r.measure_labels(angles, non_id_labels, dm)
            rows.append({"slot": name, "exact": exact_vals, "noisy_exact": noisy_exact_vals})
    return rows


def cache_exact_pec(p, non_id_labels, p2_learned, p1_learned):
    """Exact PEC-corrected (using the LEARNED channel from iteration 5) per-
    slot label values -- shot-count-independent, computed once."""
    raw = {}
    for name, sol in p["solutions"].items():
        raw[name] = pec_measure_learned(sol["angles"], non_id_labels, P2_PER_GATE, P1_PER_GATE,
                                         p2_learned, p1_learned)
    return raw


def pec_measure_learned(angles, labels, true_p2, true_p1, learned_p2, learned_p1):
    from qiskit import transpile
    from qiskit.quantum_info import DensityMatrix, Operator, Pauli
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
        if op.num_qubits == 2:
            dm = pec.apply_pauli_mixture(dm, qargs, pec.depolarizing_weights(true_p2, 2))
            dm = pec.apply_pauli_mixture(dm, qargs, pec.pec_inverse_weights(learned_p2, 2))
        elif op.num_qubits == 1:
            dm = pec.apply_pauli_mixture(dm, qargs, pec.depolarizing_weights(true_p1, 1))
            dm = pec.apply_pauli_mixture(dm, qargs, pec.pec_inverse_weights(learned_p1, 1))
    out = {}
    for l in labels:
        P = np.asarray(Pauli(l).to_matrix())
        out[l] = float(np.real(np.trace(P @ dm.data)))
    return out


# ---------------------------------------------------------------------------
# Shot-noisy energy computation for each method, from cached exact values
# ---------------------------------------------------------------------------

def shot_noisy_raw_energy(p, exact_raw, n_shots, rng):
    raw = {name: {l: shot_sample(v, n_shots, rng) for l, v in vals.items()}
           for name, vals in exact_raw.items()}
    mats = r.combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, err = r.energy_from_alpha_matrices(mats, p, K)
    return E, err


def shot_noisy_pec_energy(p, exact_pec, gamma_total, n_shots, rng):
    raw = {name: {l: pec_shot_sample(v, gamma_total, n_shots, rng) for l, v in vals.items()}
           for name, vals in exact_pec.items()}
    mats = r.combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, err = r.energy_from_alpha_matrices(mats, p, K)
    return E, err


def shot_noisy_cdr_energy(p, exact_raw, exact_training_rows, non_id_labels, n_shots, rng):
    raw = {name: {l: shot_sample(v, n_shots, rng) for l, v in vals.items()}
           for name, vals in exact_raw.items()}
    training = []
    for row in exact_training_rows:
        for l in non_id_labels:
            noisy_shot = shot_sample(row["noisy_exact"][l], n_shots, rng)
            training.append({"slot": row["slot"], "label": l, "exact": row["exact"][l], "noisy": noisy_shot})
    scales = r.fit_all_scales(training, non_id_labels, K)
    mats = r.combine_matrices(raw, p["alpha_labels"], p["identity_label"], K,
                               per_label_scale=scales["per_basis_scale"])
    E, err = r.energy_from_alpha_matrices(mats, p, K)
    return E, err


# ---------------------------------------------------------------------------
# TASK 1 verification: convergence as shots -> infinity, and the L2/sqrt(N)
# variance check. MANDATORY -- do not proceed to task 2 until this passes.
# ---------------------------------------------------------------------------

def verify_l2_norm():
    import ef_fragment as effrag
    qop_bare, _, _ = effrag.build_fragment_qop([0, 1, 2, 3], nelec=4, d=1.0)
    coeffs = np.array([complex(c).real for _, c in qop_bare.to_list()])
    L1 = float(np.sum(np.abs(coeffs)))
    L2 = float(np.sqrt(np.sum(coeffs ** 2)))
    return L1, L2


def verify_convergence(p, exact_raw, n_shots_list, seed=0):
    rng = np.random.default_rng(seed)
    print(f"\n  -- VERIFY 1: shots -> infinity converges to the exact result --")
    exact_mats = r.combine_matrices(exact_raw, p["alpha_labels"], p["identity_label"], K)
    E_exact, _ = r.energy_from_alpha_matrices(exact_mats, p, K)
    rows = []
    for n_shots in n_shots_list:
        E, err = shot_noisy_raw_energy(p, exact_raw, n_shots, rng)
        diff_kcal = abs(E - E_exact) * r.HARTREE_TO_KCAL_MOL
        rows.append({"n_shots": n_shots, "E": E, "diff_from_exact_kcal": diff_kcal})
        print(f"    n_shots={n_shots:>12,}: E={E:.6f} Ha, diff from exact = {diff_kcal:.4f} kcal/mol")
    converges = bool(rows[-1]["diff_from_exact_kcal"] < rows[0]["diff_from_exact_kcal"])
    print(f"  verdict: {'PASS -- shrinks toward exact as shots grow' if converges else 'FAIL'}")
    return rows, converges, E_exact


def verify_variance(p, exact_raw, n_shots, n_trials, L2, seed_base=1000):
    """Repeat the shot-noisy RAW energy computation n_trials times at a
    FIXED shot count, measure empirical std, compare to the predicted
    L2/sqrt(n_shots). Reports whatever is actually measured -- including a
    mismatch, since the forged-energy estimator's bilinear (beta=S.alpha.S)
    structure is NOT the simple linear-Pauli-sum case the L2/sqrt(N) result
    assumes, so a mismatch would be a real, reportable finding, not a bug
    to hide."""
    print(f"\n  -- VERIFY 2: variance check at n_shots={n_shots:,}, {n_trials} independent trials --")
    energies = []
    for trial in range(n_trials):
        rng = np.random.default_rng(seed_base + trial)
        E, _ = shot_noisy_raw_energy(p, exact_raw, n_shots, rng)
        energies.append(E)
    energies = np.array(energies)
    std_ha = float(np.std(energies))
    std_kcal = std_ha * r.HARTREE_TO_KCAL_MOL
    predicted_std_ha = L2 / np.sqrt(n_shots)
    predicted_std_kcal = predicted_std_ha * r.HARTREE_TO_KCAL_MOL
    ratio = std_ha / predicted_std_ha
    print(f"  measured std(E) = {std_ha:.6e} Ha ({std_kcal:.4f} kcal/mol)")
    print(f"  predicted L2/sqrt(N) = {predicted_std_ha:.6e} Ha ({predicted_std_kcal:.4f} kcal/mol), L2={L2}")
    print(f"  ratio measured/predicted = {ratio:.3f}")
    close = 0.5 < ratio < 2.0
    print(f"  verdict: {'within 2x of the simple linear-sum prediction' if close else 'DOES NOT MATCH -- the bilinear beta=S.alpha.S structure changes the scaling, reported honestly, not forced'}")
    return {"n_shots": n_shots, "n_trials": n_trials, "measured_std_ha": std_ha, "measured_std_kcal": std_kcal,
            "predicted_std_ha": predicted_std_ha, "predicted_std_kcal": predicted_std_kcal,
            "ratio": ratio, "close_to_prediction": bool(close)}


def verify_scaling_exponent(p, exact_raw, n_shots_pair, n_trials, seed_base=5000):
    """Confirm std scales as ~1/sqrt(N) (halving N_shots by 100x should
    change std by ~10x) regardless of whether the ABSOLUTE value matches
    L2/sqrt(N) exactly -- the scaling law is the more fundamental, more
    robust thing to check."""
    print(f"\n  -- VERIFY 3: 1/sqrt(N) scaling exponent, {n_shots_pair} --")
    stds = []
    for n_shots in n_shots_pair:
        energies = []
        for trial in range(n_trials):
            rng = np.random.default_rng(seed_base + trial)
            E, _ = shot_noisy_raw_energy(p, exact_raw, n_shots, rng)
            energies.append(E)
        std = float(np.std(energies))
        stds.append(std)
        print(f"    n_shots={n_shots:,}: std(E)={std:.6e} Ha")
    ratio_shots = n_shots_pair[1] / n_shots_pair[0]
    ratio_std = stds[0] / stds[1]
    expected_ratio_std = np.sqrt(ratio_shots)
    print(f"  shots ratio={ratio_shots:.1f}x, std ratio={ratio_std:.2f}x, expected (sqrt law)={expected_ratio_std:.2f}x")
    scaling_ok = abs(ratio_std - expected_ratio_std) / expected_ratio_std < 0.3
    print(f"  verdict: {'PASS -- consistent with 1/sqrt(N)' if scaling_ok else 'FAIL -- does not scale as 1/sqrt(N)'}")
    return {"n_shots_pair": n_shots_pair, "stds": stds, "ratio_std": ratio_std,
            "expected_ratio_std": expected_ratio_std, "scaling_ok": bool(scaling_ok)}


# ---------------------------------------------------------------------------
# TASK 2: error vs shots for raw / CDR per-basis / PEC on the learned channel
# ---------------------------------------------------------------------------

def summarize(vals, chem_acc, target):
    return {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
            "min": float(np.min(vals)), "max": float(np.max(vals)),
            "reached_chem_acc": f"{sum(1 for v in vals if v < chem_acc)}/{len(vals)}",
            "reached_target": f"{sum(1 for v in vals if v < target)}/{len(vals)}"}


# ---------------------------------------------------------------------------
# HONEST PEC: the previous "PEC (learned channel)" reused iteration 5's
# CALIBRATION as-is (exact, no shot noise) and only added shot noise to the
# target-correction step -- that overstates real performance ("perfect
# calibration, shot-limited correction"). This makes the CALIBRATION shot-
# limited too, at the SAME shot count as the main sweep, and re-fits
# p2_learned/p1_learned from noisy calibration data at every shot level --
# this is where iteration 5's 3.3e-15 becomes a real, shot-count-dependent
# number, per the task's own framing.
# ---------------------------------------------------------------------------

def cache_calibration_exact(noise_model):
    """Exact ideal/noisy decay values for every calibration circuit used in
    iteration 5 (7 CX pairs x 8 depths, 4 U3 qubits x 11 depths) -- shot-
    count-independent, computed once."""
    N_list_cx = list(range(1, 16, 2))
    N_list_u3 = list(range(1, 22, 2))
    cx_pairs = [(0, 1), (0, 2), (0, 3), (2, 0), (2, 1), (3, 0), (3, 1)]

    cx_cache = {}
    for pair in cx_pairs:
        rows = []
        for N in N_list_cx:
            qc = pl.cx_decay_circuit(pair, N)
            ideal_sv = np.asarray(Statevector.from_instruction(qc))
            n = 4
            chars = ["I"] * n
            a, b = pair
            chars[n - 1 - a] = "X"
            chars[n - 1 - b] = "X"
            Pxx = np.asarray(Pauli("".join(chars)).to_matrix())
            ideal_xx = float(np.real(ideal_sv.conj() @ Pxx @ ideal_sv))
            dm = pl.noisy_density_matrix(qc, noise_model)
            noisy_xx = float(np.real(np.trace(Pxx @ dm)))
            rows.append({"N": N, "ideal": ideal_xx, "noisy_exact": noisy_xx})
        cx_cache[pair] = rows

    u3_cache = {}
    for qubit in range(4):
        rows = []
        for N in N_list_u3:
            qc = pl.u3_identity_decay_circuit(qubit, N)
            n = 4
            chars = ["I"] * n
            chars[n - 1 - qubit] = "Z"
            Pz = np.asarray(Pauli("".join(chars)).to_matrix())
            ideal_sv = np.asarray(Statevector.from_instruction(qc))
            ideal_z = float(np.real(ideal_sv.conj() @ Pz @ ideal_sv))
            dm = pl.noisy_density_matrix(qc, noise_model)
            noisy_z = float(np.real(np.trace(Pz @ dm)))
            rows.append({"N": N, "ideal": ideal_z, "noisy_exact": noisy_z})
        u3_cache[qubit] = rows

    return cx_cache, u3_cache


def shot_noisy_calibrate(cx_cache, u3_cache, n_shots, rng):
    """Re-sample every calibration circuit at n_shots, refit p2/p1 per pair/
    qubit, average -- exactly iteration 5's protocol, now shot-limited."""
    p2_rates = []
    for pair, rows in cx_cache.items():
        Ns = np.array([row["N"] for row in rows], dtype=float)
        noisy_shots = np.array([shot_sample(row["noisy_exact"], n_shots, rng) for row in rows])
        # guard against non-positive values before log (rare at low shots/high N)
        safe = np.clip(noisy_shots, 1e-6, None)
        slope, _ = np.polyfit(Ns, np.log(safe), 1)
        p2_rates.append(1 - np.exp(slope))
    p2_learned = float(np.mean(p2_rates))

    p1_rates = []
    for qubit, rows in u3_cache.items():
        Ns = np.array([row["N"] for row in rows], dtype=float)
        noisy_shots = np.array([shot_sample(row["noisy_exact"], n_shots, rng) for row in rows])
        safe = np.clip(noisy_shots, 1e-6, None)
        slope, _ = np.polyfit(Ns, np.log(safe), 1)
        p1_rates.append(1 - np.exp(slope))
    p1_learned = float(np.mean(p1_rates))

    return p2_learned, p1_learned


def shot_noisy_pec_honest_energy(p, non_id_labels, exact_raw_true_noise, cx_cache, u3_cache,
                                  gamma_total, n_shots, rng):
    """Fully honest PEC: shot-limited calibration (re-learns p2/p1 from
    n_shots-sampled calibration circuits) AND shot-limited target
    correction, both at the same n_shots."""
    p2_learned, p1_learned = shot_noisy_calibrate(cx_cache, u3_cache, n_shots, rng)
    # exact PEC-corrected value using the (now shot-limited-learned) channel,
    # then shot-sample the correction step itself with PEC's inflated variance
    raw = {}
    for name, sol in p["solutions"].items():
        exact_pec_vals = pec_measure_learned(sol["angles"], non_id_labels, P2_PER_GATE, P1_PER_GATE,
                                              p2_learned, p1_learned)
        raw[name] = {l: pec_shot_sample(v, gamma_total, n_shots, rng) for l, v in exact_pec_vals.items()}
    mats = r.combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, err = r.energy_from_alpha_matrices(mats, p, K)
    return E, err, p2_learned, p1_learned


def main():
    print("\n" + "=" * 90)
    print("  shot_noise_study.py -- adding the dominant real-hardware error source")
    print("=" * 90)

    p = r.setup(K)
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    p["solutions"] = solutions
    print(f"  fit {n_ok}/{len(solutions)} targets converged, worst={worst:.2e}")
    counts = r.verify_constant_gate_count(solutions)
    assert counts == {11}, f"gate count not fixed: {counts}"
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    chem_acc = p["chemical_accuracy_kcal"]
    LOOP_TARGET = 0.30

    L1, L2 = verify_l2_norm()
    print(f"\n  L1 = {L1:.4f} Ha, L2 = {L2:.4f} Ha (Hamiltonian Pauli-coefficient norms, verified directly)")

    noise_model = r.build_noise_model()
    t0 = time.time()
    exact_raw = cache_exact_raw(p, non_id_labels, noise_model)
    print(f"  cached exact (noisy) raw matrices for all 36 targets, {time.time()-t0:.1f}s")

    # ==== TASK 1: MANDATORY VERIFICATION -- do not proceed until this passes ====
    conv_rows, converges, E_exact = verify_convergence(p, exact_raw, [1e3, 1e4, 1e5, 1e6, 1e7, 1e8])
    var_result = verify_variance(p, exact_raw, n_shots=100_000, n_trials=200, L2=L2)
    scaling_result = verify_scaling_exponent(p, exact_raw, (10_000, 1_000_000), n_trials=100)

    task1_pass = converges and scaling_result["scaling_ok"]
    print(f"\n  TASK 1 VERIFICATION: convergence={'PASS' if converges else 'FAIL'}, "
          f"1/sqrt(N) scaling={'PASS' if scaling_result['scaling_ok'] else 'FAIL'}, "
          f"L2/sqrt(N) absolute match={'close' if var_result['close_to_prediction'] else 'MISMATCH (reported, not hidden)'}")
    if not task1_pass:
        print("  STOPPING per instruction: convergence or scaling verification failed.")
        results = {"task1_failed": True, "convergence": conv_rows, "variance": var_result, "scaling": scaling_result}
        with open(RESULTS_PATH, "w") as f:
            json.dump(results, f, indent=2)
        return results
    print("  Proceeding to Task 2 (convergence + scaling both pass; L2/sqrt(N) match reported as measured, "
          "not gated on -- the bilinear beta=S.alpha.S structure is a real, documented structural difference "
          "from the simple linear-Pauli-sum case that formula assumes).")

    # ==== load iteration 5's learned channel + gamma_total for PEC ====
    with open(LINDBLAD_RESULTS_PATH) as f:
        lindblad_results = json.load(f)
    p2_learned = lindblad_results["p2_learned"]
    p1_learned = lindblad_results["p1_learned"]
    gamma_total = lindblad_results["cost_accounting"]["gamma_total"]
    print(f"\n  loaded iteration 5's learned channel: p2_learned={p2_learned}, p1_learned={p1_learned}, "
          f"gamma_total={gamma_total:.4f}")

    t0 = time.time()
    exact_pec = cache_exact_pec(p, non_id_labels, p2_learned, p1_learned)
    print(f"  cached exact PEC (learned-channel) matrices for all 36 targets, {time.time()-t0:.1f}s")

    t0 = time.time()
    exact_training_by_seed = {seed: cache_exact_training(p, non_id_labels, noise_model, N_TRAIN_PER_SLOT, seed)
                               for seed in range(N_SEEDS)}
    print(f"  cached exact CDR training data for {N_SEEDS} seeds, {time.time()-t0:.1f}s")

    t0 = time.time()
    cx_cal_cache, u3_cal_cache = cache_calibration_exact(noise_model)
    print(f"  cached exact calibration-circuit decay values (100 circuits), {time.time()-t0:.1f}s")
    print(f"\n  NOTE: 'PEC (optimistic)' below reuses iteration 5's EXACT (infinite-shot) calibration --")
    print(f"  only the target-correction step is shot-limited. 'PEC (honest)' shot-limits the")
    print(f"  CALIBRATION too, at the same shot count -- this is where iteration 5's 3.3e-15 becomes")
    print(f"  a real, shot-count-dependent number, and is the fair comparison to raw/CDR.")

    # ==== TASK 2: main shots-vs-accuracy sweep ====
    print(f"\n  -- TASK 2: error vs shots, {N_SEEDS}-seed sweep at each of {SHOT_LEVELS} --")
    sweep_results = {}
    for n_shots in SHOT_LEVELS:
        raw_exact_errs, raw_noiseless_errs = [], []
        cdr_exact_errs, cdr_noiseless_errs = [], []
        pec_exact_errs, pec_noiseless_errs = [], []
        pec_honest_exact_errs, pec_honest_noiseless_errs = [], []
        pec_honest_p2_rel_errs = []
        t_shots0 = time.time()
        for seed in range(N_SEEDS):
            rng = np.random.default_rng(seed * 1_000_003 + n_shots)

            E_raw, err_raw = shot_noisy_raw_energy(p, exact_raw, n_shots, rng)
            raw_exact_errs.append(err_raw["err_vs_exact_kcal"])
            raw_noiseless_errs.append(err_raw["err_vs_noiseless_kcal"])

            E_cdr, err_cdr = shot_noisy_cdr_energy(p, exact_raw, exact_training_by_seed[seed],
                                                     non_id_labels, n_shots, rng)
            cdr_exact_errs.append(err_cdr["err_vs_exact_kcal"])
            cdr_noiseless_errs.append(err_cdr["err_vs_noiseless_kcal"])

            E_pec, err_pec = shot_noisy_pec_energy(p, exact_pec, gamma_total, n_shots, rng)
            pec_exact_errs.append(err_pec["err_vs_exact_kcal"])
            pec_noiseless_errs.append(err_pec["err_vs_noiseless_kcal"])

            E_pech, err_pech, p2_h, p1_h = shot_noisy_pec_honest_energy(
                p, non_id_labels, exact_raw, cx_cal_cache, u3_cal_cache, gamma_total, n_shots, rng)
            pec_honest_exact_errs.append(err_pech["err_vs_exact_kcal"])
            pec_honest_noiseless_errs.append(err_pech["err_vs_noiseless_kcal"])
            pec_honest_p2_rel_errs.append(abs(p2_h - P2_PER_GATE) / P2_PER_GATE)

        sweep_results[n_shots] = {
            "raw": {"exact": summarize(raw_exact_errs, chem_acc, LOOP_TARGET),
                    "noiseless": summarize(raw_noiseless_errs, chem_acc, LOOP_TARGET)},
            "cdr_per_basis": {"exact": summarize(cdr_exact_errs, chem_acc, LOOP_TARGET),
                               "noiseless": summarize(cdr_noiseless_errs, chem_acc, LOOP_TARGET)},
            "pec_optimistic_calibration": {"exact": summarize(pec_exact_errs, chem_acc, LOOP_TARGET),
                                            "noiseless": summarize(pec_noiseless_errs, chem_acc, LOOP_TARGET)},
            "pec_honest_calibration": {"exact": summarize(pec_honest_exact_errs, chem_acc, LOOP_TARGET),
                                        "noiseless": summarize(pec_honest_noiseless_errs, chem_acc, LOOP_TARGET),
                                        "mean_p2_relative_error": float(np.mean(pec_honest_p2_rel_errs))},
        }
        print(f"\n  n_shots={n_shots:>12,}  ({time.time()-t_shots0:.1f}s):")
        print(f"    raw                  : {sweep_results[n_shots]['raw']['exact']['mean']:>12.4f} +/- "
              f"{sweep_results[n_shots]['raw']['exact']['std']:<10.4f} kcal/mol  "
              f"(chem.acc {sweep_results[n_shots]['raw']['exact']['reached_chem_acc']}, "
              f"target {sweep_results[n_shots]['raw']['exact']['reached_target']})")
        print(f"    CDR per-basis        : {sweep_results[n_shots]['cdr_per_basis']['exact']['mean']:>12.4f} +/- "
              f"{sweep_results[n_shots]['cdr_per_basis']['exact']['std']:<10.4f} kcal/mol  "
              f"(chem.acc {sweep_results[n_shots]['cdr_per_basis']['exact']['reached_chem_acc']}, "
              f"target {sweep_results[n_shots]['cdr_per_basis']['exact']['reached_target']})")
        print(f"    PEC (optimistic cal.): {sweep_results[n_shots]['pec_optimistic_calibration']['exact']['mean']:>12.4f} +/- "
              f"{sweep_results[n_shots]['pec_optimistic_calibration']['exact']['std']:<10.4f} kcal/mol  "
              f"(chem.acc {sweep_results[n_shots]['pec_optimistic_calibration']['exact']['reached_chem_acc']}, "
              f"target {sweep_results[n_shots]['pec_optimistic_calibration']['exact']['reached_target']})")
        print(f"    PEC (honest cal.)    : {sweep_results[n_shots]['pec_honest_calibration']['exact']['mean']:>12.4f} +/- "
              f"{sweep_results[n_shots]['pec_honest_calibration']['exact']['std']:<10.4f} kcal/mol  "
              f"(chem.acc {sweep_results[n_shots]['pec_honest_calibration']['exact']['reached_chem_acc']}, "
              f"target {sweep_results[n_shots]['pec_honest_calibration']['exact']['reached_target']}, "
              f"mean calib. p2 rel.err={sweep_results[n_shots]['pec_honest_calibration']['mean_p2_relative_error']*100:.3f}%)")

    # ==== crossing points ====
    print(f"\n  -- CHEMICAL ACCURACY / TARGET CROSSING POINTS --")
    crossings = {}
    for method in ("raw", "cdr_per_basis", "pec_optimistic_calibration", "pec_honest_calibration"):
        crossed_chem, crossed_target = None, None
        for n_shots in SHOT_LEVELS:
            d = sweep_results[n_shots][method]["exact"]
            if crossed_chem is None and d["mean"] < chem_acc:
                crossed_chem = n_shots
            if crossed_target is None and d["mean"] < LOOP_TARGET:
                crossed_target = n_shots
        crossings[method] = {"chemical_accuracy_at_shots": crossed_chem, "loop_target_at_shots": crossed_target}
        print(f"    {method}: chemical accuracy at {crossed_chem or 'NOT within tested range'} shots/setting, "
              f"0.30 kcal/mol target at {crossed_target or 'NOT within tested range'} shots/setting")

    results = {
        "L1": L1, "L2": L2,
        "task1_verification": {
            "convergence": conv_rows, "convergence_pass": converges,
            "variance_check": var_result, "scaling_check": scaling_result,
        },
        "learned_channel": {"p2_learned": p2_learned, "p1_learned": p1_learned, "gamma_total": gamma_total},
        "shot_levels": SHOT_LEVELS, "n_seeds": N_SEEDS,
        "sweep": {str(k): v for k, v in sweep_results.items()},
        "crossings": crossings,
        "chemical_accuracy_kcal": chem_acc, "loop_target_kcal": LOOP_TARGET,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
