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
IONQ_ARIA_FORTE_FIDELITY = 0.98786
QUANTINUUM_H1_FIDELITY = 0.9782  # 1 - QUANTINUUM_TWO_Q_ERROR-derived, matches entanglement_forging_h4.py's own numbers
QUANTINUUM_H2_FIDELITY = 0.9891
# P2 sweep: 0.015 (98.5%) down to 0.0001 (99.99%), log-spaced for even coverage of the crossing region
P2_VALUES = sorted(set(np.round(np.geomspace(0.0001, 0.015, 22), 6).tolist()), reverse=True)
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

    print(f"\n  -- reference fidelities --")
    print(f"    IonQ Aria/Forte : {IONQ_ARIA_FORTE_FIDELITY:.4%}")
    print(f"    Quantinuum H1   : {QUANTINUUM_H1_FIDELITY:.4%}")
    print(f"    Quantinuum H2   : {QUANTINUUM_H2_FIDELITY:.4%}")
    for name, fid in (("IonQ Aria/Forte", IONQ_ARIA_FORTE_FIDELITY),
                       ("Quantinuum H1", QUANTINUUM_H1_FIDELITY), ("Quantinuum H2", QUANTINUUM_H2_FIDELITY)):
        p2_ref = 1 - fid
        # interpolate each method's error at this reference fidelity
        row_str = f"    at {name} (F={fid:.4%}, p2={p2_ref:.5f}): "
        for method in ("raw_kcal", "zne_quadratic_kcal", "pec_kcal"):
            errs = [r[method] for r in rows]
            interp = float(np.interp(p2_ref, p2_vals, errs))
            row_str += f"{method}={interp:.2f}  "
        print(row_str)

    results = {
        "K": K, "p2_p1_ratio": 40, "fidelity_convention": "fidelity = 1 - p2 (not the average-gate-fidelity formula)",
        "sweep": rows,
        "crossings": crossings,
        "reference_fidelities": {
            "ionq_aria_forte": IONQ_ARIA_FORTE_FIDELITY,
            "quantinuum_h1": QUANTINUUM_H1_FIDELITY,
            "quantinuum_h2": QUANTINUUM_H2_FIDELITY,
        },
        "chemical_accuracy_kcal": CHEM_ACC_KCAL,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
