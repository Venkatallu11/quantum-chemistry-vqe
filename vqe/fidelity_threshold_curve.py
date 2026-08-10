#!/usr/bin/env python3
"""
fidelity_threshold_curve.py — Task 5: THE DELIVERABLE. Sweeps two-qubit
gate fidelity and asks, precisely, for each mitigation method: at what
fidelity does its error cross chemical accuracy (1.0 kcal/mol)?
============================================================================
This converts "it didn't reproduce on IonQ" (iteration 9's real-hardware
finding, raw ~35-43 kcal/mol; Task 2's finding that the ORIGINAL
StatePreparation circuit is worse still, ~123-135 kcal/mol) into a
quantitative statement: H4 forged VQE (this fixed 11-gate ansatz, K=6)
needs two-qubit fidelity >= X%, evidenced by sweeping the SAME
depolarizing-channel parameter every method in this ledger already uses,
not a new noise model invented for this file.

FIDELITY CONVENTION, stated explicitly (this is the task's own framing,
not re-derived from the more precise average-gate-fidelity formula
F = 1 - p*(d^2-1)/d^2): fidelity = 1 - p2. p1 = p2/40 throughout,
matching this project's established ratio (fixed_ansatz.P2_PER_GATE /
P1_PER_GATE).

FIVE METHODS, ALL LOCAL (no network, no real hardware -- this is a sweep
over a CLASSICAL parameter of a noise MODEL, not a hardware measurement):
  - raw: no mitigation
  - ZNE-linear / ZNE-quadratic: noise scaled x1/x2/x3, fit + extrapolate
    to scale=0 (matching entanglement_forging_zne.py's own SCALES)
  - CDR: per-basis scale, fresh training data generated AT EACH p2 (a
    method trained on the WRONG noise level is not evidence about the
    RIGHT one)
  - PEC: gate-by-gate EXACT quasi-probability inverse (the channel is
    exactly known at every swept p2 by construction -- this is loop_pec.py's
    original, best-case framing, not iteration 9's real-hardware-learned-
    channel framing; stated explicitly, not conflated)

Run:
    python vqe/fidelity_threshold_curve.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, HARTREE_TO_KCAL_MOL,
    combine_matrices, energy_from_alpha_matrices, CDRStrategy, build_ansatz, slot_names,
)
import loop_pec as pec
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.quantum_info import Pauli, DensityMatrix, Operator

K = 6
BASIS_GATES = ["u3", "cx"]
SCALES = [1.0, 2.0, 3.0]
N_TRAIN_PER_SLOT = 5
LOW_SIGNAL_CUTOFF = 0.05
CHEM_ACC_KCAL = 1.0

# CORRECTED iteration 24, Task 0 (was WRONG before this fix -- see
# RESEARCH_LEDGER.md Task 0 write-up). This project's own local-model
# constant, used throughout as "real aria-1/forte-1" -- NOT a live device
# reading, just this project's fixed depolarizing-model parameter.
LOCAL_MODEL_CONSTANT_FIDELITY = 0.98786

# Queried LIVE from IonQ's real /backends/<name>/characterizations API,
# 2026-08-10 (see task0_fidelity_correction.py). "2q" field is the
# device-wide MEDIAN two-qubit gate fidelity IonQ's own API reports (not
# a best-pair number, not a mean -- the API exposes only median+stderr).
IONQ_FORTE1_REAL_FIDELITY = 0.9952   # qpu.forte-1, AVAILABLE, characterization dated 2026-08-09
# qpu.aria-1 is RETIRED (confirmed live via /backends). Its own MOST
# RECENT characterization record (2026-02-19) has 2q=null and a
# nonsensical 1q=0.4745 -- clearly a stale/non-representative end-of-life
# record, not used. Its LAST VALID reading with a real 2q number:
IONQ_ARIA1_LAST_VALID_FIDELITY = 0.9820  # qpu.aria-1, characterization dated 2025-09-02, RETIRED since

# FIXED (was 0.9782/0.9891 -- those were (0.998)^11-type 11-GATE CIRCUIT
# fidelities mislabeled as single-gate reference fidelities, matching this
# project's own 11-CX ansatz gate count by coincidence, not derived from a
# real per-gate number). Real per-gate values, from this project's own
# QUANTINUUM_TWO_Q_ERROR=0.002 / QUANTINUUM_ONE_Q_ERROR=0.00005 constants
# (entanglement_forging_h4.py), i.e. genuinely 1 - per_gate_error:
QUANTINUUM_H1_FIDELITY = 0.998    # 1 - QUANTINUUM_TWO_Q_ERROR, per 2-qubit gate
QUANTINUUM_H2_FIDELITY = 0.999    # H2 characterized as modestly better than H1 in this project's own prior notes

# P2 sweep: extended to 0.02 (was 0.015) to cover aria-1's real last-valid
# p2=0.018, which fell OUTSIDE the original swept range
P2_VALUES = sorted(set(np.round(np.geomspace(0.0001, 0.02, 24), 6).tolist()), reverse=True)
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "fidelity_threshold_curve_results.json")


def build_noise_model(p2, p1):
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), "u3")
    return nm


def noisy_density_matrix(angles, noise_model):
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    qc2 = qc.copy()
    qc2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc2).result()
    return np.asarray(result.data(0)["density_matrix"])


def measure_exact_noisy_raw(p, non_id_labels, noise_model):
    raw = {}
    for name, sol in p["solutions"].items():
        dm = noisy_density_matrix(sol["angles"], noise_model)
        raw[name] = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm))) for l in non_id_labels}
    return raw


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def generate_cdr_training(p, non_id_labels, noise_model, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for name in slot_names(K):
        for _ in range(N_TRAIN_PER_SLOT):
            angles = rng.uniform(-np.pi, np.pi, 5)
            sv_exact = {}
            from qiskit.quantum_info import Statevector
            sv = Statevector.from_instruction(build_ansatz(angles.tolist()))
            for l in non_id_labels:
                sv_exact[l] = float(sv.expectation_value(Pauli(l)).real)
            dm = noisy_density_matrix(angles.tolist(), noise_model)
            for l in non_id_labels:
                noisy_v = float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm)))
                rows.append({"slot": name, "label": l, "exact": sv_exact[l], "noisy": noisy_v})
    return rows


def measure_pec_exact(p, non_id_labels, true_p2, true_p1):
    """Gate-by-gate EXACT quasi-probability inverse, channel exactly
    known (loop_pec.py's original framing) -- deterministic, no shot
    noise, no learning."""
    raw = {}
    for name, sol in p["solutions"].items():
        qc = transpile(build_ansatz(sol["angles"]), basis_gates=BASIS_GATES, optimization_level=0)
        dm = DensityMatrix.from_label("0" * qc.num_qubits)
        for instr in qc.data:
            op = instr.operation
            if op.name in ("measure", "barrier"):
                continue
            qargs = [qc.find_bit(q).index for q in instr.qubits]
            dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
            if op.num_qubits == 2:
                dm = pec.apply_pauli_mixture(dm, qargs, pec.depolarizing_weights(true_p2, 2))
                dm = pec.apply_pauli_mixture(dm, qargs, pec.pec_inverse_weights(true_p2, 2))
            elif op.num_qubits == 1:
                dm = pec.apply_pauli_mixture(dm, qargs, pec.depolarizing_weights(true_p1, 1))
                dm = pec.apply_pauli_mixture(dm, qargs, pec.pec_inverse_weights(true_p1, 1))
        raw[name] = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm.data))) for l in non_id_labels}
    return raw


def crossing_fidelity(p2_values, errs, target_kcal=CHEM_ACC_KCAL):
    """Linear-interpolate (in log(p2)) the p2 at which err crosses
    target_kcal, scanning from HIGH fidelity (low p2) to LOW fidelity
    (high p2) -- returns None if never crosses in the swept range."""
    pairs = sorted(zip(p2_values, errs))  # ascending p2 (= descending fidelity)
    for i in range(len(pairs) - 1):
        p2_a, e_a = pairs[i]
        p2_b, e_b = pairs[i + 1]
        if (e_a - target_kcal) * (e_b - target_kcal) <= 0 and e_a != e_b:
            frac = (target_kcal - e_a) / (e_b - e_a)
            log_p2 = np.log(p2_a) + frac * (np.log(p2_b) - np.log(p2_a))
            p2_cross = float(np.exp(log_p2))
            return p2_cross, 1 - p2_cross
    return None, None


def main():
    print("\n" + "=" * 96)
    print("  fidelity_threshold_curve.py -- Task 5, THE DELIVERABLE")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K)
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    assert n_ok == 36
    counts = verify_constant_gate_count(solutions)
    assert counts == {11}
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    print(f"  setup OK: 36/36 converged, gate count={counts}")
    print(f"  sweeping p2 over {len(P2_VALUES)} values from {P2_VALUES[0]} (fidelity="
          f"{1-P2_VALUES[0]:.4%}) to {P2_VALUES[-1]} (fidelity={1-P2_VALUES[-1]:.4%})")

    rows = []
    for i, p2 in enumerate(P2_VALUES):
        p1 = p2 / 40
        fidelity = 1 - p2

        energies_by_scale = []
        for s in SCALES:
            nm = build_noise_model(p2 * s, p1 * s)
            raw = measure_exact_noisy_raw(p, non_id_labels, nm)
            E, _ = energy_and_err(p, raw, K)
            energies_by_scale.append(E)
        raw_err = abs(energies_by_scale[0] - p["exact_energy"]) * HARTREE_TO_KCAL_MOL

        lin_coeffs = np.polyfit(SCALES, energies_by_scale, 1)
        E_lin0 = float(np.polyval(lin_coeffs, 0))
        zne_lin_err = abs(E_lin0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL

        quad_coeffs = np.polyfit(SCALES, energies_by_scale, 2)
        E_quad0 = float(np.polyval(quad_coeffs, 0))
        zne_quad_err = abs(E_quad0 - p["exact_energy"]) * HARTREE_TO_KCAL_MOL

        nm1 = build_noise_model(p2, p1)
        training = generate_cdr_training(p, non_id_labels, nm1, seed=12345)
        from qforge.mitigation import fit_all_scales
        scales_fit = fit_all_scales(training, non_id_labels, K, LOW_SIGNAL_CUTOFF)
        raw1 = measure_exact_noisy_raw(p, non_id_labels, nm1)
        alpha_mats_cdr = combine_matrices(raw1, p["alpha_labels"], p["identity_label"], K,
                                           per_label_scale=scales_fit["per_basis_scale"])
        E_cdr, errs_cdr = energy_from_alpha_matrices(alpha_mats_cdr, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                                      exact_energy=p["exact_energy"])
        cdr_err = errs_cdr["err_vs_exact_kcal"]

        raw_pec = measure_pec_exact(p, non_id_labels, p2, p1)
        E_pec, pec_err = energy_and_err(p, raw_pec, K)

        rows.append({"p2": p2, "p1": p1, "fidelity": fidelity, "raw_kcal": raw_err,
                     "zne_linear_kcal": zne_lin_err, "zne_quadratic_kcal": zne_quad_err,
                     "cdr_kcal": cdr_err, "pec_kcal": pec_err})
        print(f"    [{i+1:>2}/{len(P2_VALUES)}] p2={p2:.5f} (F={fidelity:.4%}): raw={raw_err:7.2f}  "
              f"zne_lin={zne_lin_err:7.2f}  zne_quad={zne_quad_err:7.2f}  cdr={cdr_err:7.2f}  pec={pec_err:7.3f}")

    print(f"\n  -- CROSSING FIDELITY (error crosses {CHEM_ACC_KCAL} kcal/mol chemical accuracy) --")
    p2_vals = [r["p2"] for r in rows]
    crossings = {}
    for method in ("raw_kcal", "zne_linear_kcal", "zne_quadratic_kcal", "cdr_kcal", "pec_kcal"):
        errs = [r[method] for r in rows]
        p2_cross, fid_cross = crossing_fidelity(p2_vals, errs)
        crossings[method] = {"p2_crossing": p2_cross, "fidelity_crossing": fid_cross}
        if fid_cross is not None:
            print(f"    {method:>18}: crosses at F={fid_cross:.4%} (p2={p2_cross:.5f})")
        else:
            print(f"    {method:>18}: does NOT cross {CHEM_ACC_KCAL} kcal/mol in the swept range "
                  f"[{1-P2_VALUES[0]:.4%}, {1-P2_VALUES[-1]:.4%}]")

    print(f"\n  -- reference fidelities (CORRECTED, iteration 24 Task 0) --")
    print(f"    this project's local-model constant : {LOCAL_MODEL_CONSTANT_FIDELITY:.4%}  "
          f"(NOT a live reading -- fixed_ansatz.P2_PER_GATE, applied identically to both aria-1 and forte-1 local sims)")
    print(f"    IonQ forte-1 REAL (live, 2026-08-09) : {IONQ_FORTE1_REAL_FIDELITY:.4%}  (AVAILABLE, device-wide median)")
    print(f"    IonQ aria-1 REAL (last valid, 2025-09-02): {IONQ_ARIA1_LAST_VALID_FIDELITY:.4%}  (RETIRED since; device-wide median)")
    print(f"    Quantinuum H1 (per-gate, FIXED)      : {QUANTINUUM_H1_FIDELITY:.4%}")
    print(f"    Quantinuum H2 (per-gate, FIXED)      : {QUANTINUUM_H2_FIDELITY:.4%}")

    # np.interp requires xp ASCENDING -- p2_vals comes out of `rows` in
    # DESCENDING order (P2_VALUES was built sorted reverse=True), so
    # interpolate against an explicitly ascending-sorted copy (same fix
    # crossing_fidelity() already applies via sorted(zip(...)) above;
    # this loop originally skipped it, silently returning nonsense).
    p2_vals_asc, idx_asc = np.unique(p2_vals, return_index=True)
    reference_points = {}
    for name, fid in (("local_model_constant", LOCAL_MODEL_CONSTANT_FIDELITY),
                       ("ionq_forte1_real", IONQ_FORTE1_REAL_FIDELITY),
                       ("ionq_aria1_real_last_valid", IONQ_ARIA1_LAST_VALID_FIDELITY),
                       ("quantinuum_h1", QUANTINUUM_H1_FIDELITY), ("quantinuum_h2", QUANTINUUM_H2_FIDELITY)):
        p2_ref = 1 - fid
        row = {"fidelity": fid, "p2": p2_ref}
        row_str = f"    at {name:<28} (F={fid:.4%}, p2={p2_ref:.5f}): "
        for method in ("raw_kcal", "zne_linear_kcal", "zne_quadratic_kcal", "cdr_kcal", "pec_kcal"):
            errs_asc = np.array([r[method] for r in rows])[idx_asc]
            interp = float(np.interp(p2_ref, p2_vals_asc, errs_asc))
            row[method] = interp
            row_str += f"{method}={interp:.2f}  "
        reference_points[name] = row
        print(row_str)

    # gap decomposition: how much of the real-hardware raw gap (iteration
    # 9's actual submitted numbers) is explained by the LOCAL depolarizing
    # model's own miscalibration, vs everything this parameter sweep
    # cannot capture at all (SPAM, leakage, coherent noise, crosstalk,
    # shot noise -- already independently investigated elsewhere in this
    # ledger, e.g. iteration 22's coherent-noise RC study)
    REAL_HARDWARE_RAW_KCAL = {"aria-1": 34.98, "forte-1": 43.03}  # iteration 9's actual submitted result
    raw_at_constant = reference_points["local_model_constant"]["raw_kcal"]
    raw_at_forte1_real = reference_points["ionq_forte1_real"]["raw_kcal"]
    raw_at_aria1_real = reference_points["ionq_aria1_real_last_valid"]["raw_kcal"]
    print(f"\n  -- gap decomposition: how much of the real-hardware raw gap was miscalibration? --")
    print(f"    forte-1: real submitted={REAL_HARDWARE_RAW_KCAL['forte-1']:.2f} kcal/mol; "
          f"local depolarizing-only model predicts {raw_at_constant:.2f} kcal/mol at the assumed constant, "
          f"{raw_at_forte1_real:.2f} kcal/mol at forte-1's REAL fidelity")
    print(f"    aria-1:  real submitted={REAL_HARDWARE_RAW_KCAL['aria-1']:.2f} kcal/mol; "
          f"local depolarizing-only model predicts {raw_at_constant:.2f} kcal/mol at the assumed constant, "
          f"{raw_at_aria1_real:.2f} kcal/mol at aria-1's REAL last-valid fidelity")
    print(f"    HONEST CAVEAT: this depolarizing-only sweep was never a full hardware model (no SPAM, no "
          f"leakage, no coherent noise, no shot noise -- iteration 22 already found real IonQ noise has a "
          f"genuine coherent component this exact sweep cannot see). It answers a narrower question: given "
          f"ONLY the depolarizing-rate assumption, how much of the raw gap traces to picking the wrong rate, "
          f"vs everything else. It is not a claim that the full 33-43 kcal/mol gap is explained.")

    results = {
        "K": K, "p2_p1_ratio": 40, "fidelity_convention": "fidelity = 1 - p2 (not the average-gate-fidelity formula)",
        "sweep": rows,
        "crossings": crossings,
        "reference_fidelities": {
            "local_model_constant": LOCAL_MODEL_CONSTANT_FIDELITY,
            "ionq_forte1_real_live_2026_08_09": IONQ_FORTE1_REAL_FIDELITY,
            "ionq_aria1_real_last_valid_2025_09_02_RETIRED": IONQ_ARIA1_LAST_VALID_FIDELITY,
            "quantinuum_h1_per_gate_FIXED": QUANTINUUM_H1_FIDELITY,
            "quantinuum_h2_per_gate_FIXED": QUANTINUUM_H2_FIDELITY,
        },
        "reference_points": reference_points,
        "real_hardware_raw_kcal_iteration9": REAL_HARDWARE_RAW_KCAL,
        "correction_note": (
            "iteration 24 Task 0: fixed WRONG quantinuum_h1/h2 reference fidelities "
            "(were 0.9782/0.9891, an 11-gate CIRCUIT fidelity mislabeled as a per-gate "
            "number; real per-gate values are ~0.998/0.999). Added REAL live-queried "
            "IonQ calibration (forte-1=99.52% available; aria-1's last valid reading "
            "before retirement=98.20%, its current record has a null/broken 2q entry). "
            "The local_model_constant (98.786%) was never a live reading -- it is this "
            "project's own fixed depolarizing-model parameter, applied identically to "
            "both backends despite forte-1 being meaningfully better."
        ),
        "chemical_accuracy_kcal": CHEM_ACC_KCAL,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
