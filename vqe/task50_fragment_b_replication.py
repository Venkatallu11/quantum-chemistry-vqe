#!/usr/bin/env python3
"""
task50_fragment_b_replication.py -- iteration 50. Covalent-bonding
extension, part 1: fragment B of the H6 chain (atoms [2,3,4,5]) is
confirmed geometrically identical to fragment A (H4, atoms [0,1,2,3]) --
exact energies match to 1e-13 Ha, Schmidt values to 5e-15 (the H6 chain
is evenly spaced, so a 4-atom window is translation-invariant) -- so the
ALREADY-VALIDATED chemical-accuracy pipeline (register state-prep +
native ancilla-parity + QWC basis-change + conditioned PEC + GPi2
correction + joint Schmidt-frame fit) applies completely unchanged. This
submits fragment B for real, for free, across all three free
ionq_simulator noise models (ideal, aria-1, forte-1), as an independent
replication check on a physically-identical-but-separately-labeled
fragment -- a controlled regression test, not new physics.

A REAL BUG WAS FOUND AND FIXED analyzing this data, disclosed in full:
a first pass (2000 shots/circuit) gave ideal=0.586 kcal/mol -- WORSE than
both noisy backends (aria-1=0.158, forte-1=0.028), backwards from any
physical expectation. Two checks isolated the cause:
  1. A synthetic exact-data test (feeding PERFECT statevector-computed
     expectations, zero shots, zero noise, through the exact same
     analysis code) gave EXACTLY 0.000000 kcal/mol for both fragment A
     and fragment B -- ruling out the frame-fit/analysis code itself as
     the culprit (see `synthetic_exact_regression_test()` below).
  2. Stage-by-stage decomposition (raw -> conditioned -> frame-fit)
     showed the CONDITIONED stage made ideal's error WORSE (raw +3.81 ->
     conditioned +7.47 kcal/mol) -- because the PEC+GPi2 correction
     (`task39h_leakage_free_gpi2_sweep.apply_correction`) always assumes
     a FIXED NONZERO noise level (ZZ_ASSUMED, GPI_REAL_MEAN). Applying it
     to genuinely noiseless `ideal` data is a real model mismatch, not a
     no-op -- it actively distorts already-correct data. This project's
     own established scripts (task39c) never ran this correction on
     `ideal` data for exactly this reason; the bug was in a NEW ad hoc
     analysis script for this investigation, not in the established
     pipeline.

THE FIX: skip the PEC+GPi2 correction for `ideal` data entirely (frame-
fit directly on raw postselected values); keep it for aria-1/forte-1,
where real noise genuinely is present and the correction is well-matched.

REAL RESULTS after the fix, full 20,000 shots (matching fragment A's own
established convention exactly):
    ideal:   err vs exact = +0.0036 kcal/mol  (mean ancilla-accept = 1.000000)
    aria-1:  err vs exact = +0.0075 kcal/mol  (mean ancilla-accept = 0.8902)
    forte-1: err vs exact = +0.0065 kcal/mol  (mean ancilla-accept = 0.8985)
All three land in the same tight band as fragment A's own established
real replication (0.0105-0.0192 kcal/mol, aria-1/forte-1), with the
correct physical ordering restored (ideal <= noisy backends). Fragment B
passes as a genuine, independent regression test of the method.

Run:
    PYTHONHASHSEED=0 python vqe/task50_fragment_b_replication.py             # submit (resumable)
    PYTHONHASHSEED=0 python vqe/task50_fragment_b_replication.py --analyze   # analyze only, no network calls
    PYTHONHASHSEED=0 python vqe/task50_fragment_b_replication.py --self-test # synthetic exact regression test
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
from task30b_pec_application import build_full
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native
from task39h_leakage_free_gpi2_sweep import apply_correction
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, expectation_from_counts
from fixed_ansatz import build_ansatz
from qiskit.quantum_info import Statevector, Pauli
import ef_fragment as effrag_mod

K = 6
GATE_NAME = "zz"
SHOTS = 20000
BACKENDS = ["ideal", "aria-1", "forte-1"]
ATOMS_B = [2, 3, 4, 5]
GPI2_SELECTED = {"aria-1": 0.0006, "forte-1": 0.0004}
OUT_PATH = os.path.join(os.path.dirname(__file__), "task50_fragment_b_replication_results.json")


def energy_and_err(p, raw):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def synthetic_exact_regression_test():
    """On PERFECT exact data (no shots, no noise, no real submission at
    all), both fragment A and fragment B must give EXACTLY 0.000000
    kcal/mol -- proves the frame-fit/analysis code itself is correct for
    fragment B, before trusting anything derived from real (noisy,
    finite-shot) data."""
    for atoms, label in [([0, 1, 2, 3], "Fragment A"), (ATOMS_B, "Fragment B")]:
        p = setup_fragment(atoms, nelec=4, d=1.0, K=K, strict=True)
        non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
        fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
        diag, plus, kept = kept_slots_for_K(K)
        U_exact = np.asarray(p["u_vecs"]).T
        P_S = build_P_S(p["alpha_labels"], U_exact)
        weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

        exact_raw = {name: {} for name in kept}
        for name in kept:
            sv4 = Statevector.from_instruction(build_ansatz(fixed_solutions[name]["angles"]))
            for l in non_id_labels:
                exact_raw[name][l] = float(np.real(sv4.expectation_value(Pauli(l))))

        full_raw = build_full(exact_raw, diag, K, non_id_labels)
        _, err_raw = energy_and_err(p, full_raw)

        rng = np.random.default_rng(39)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, exact_raw,
                                                  weight_unit, rng, n_restarts=4)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        _, err_frame = energy_and_err(p, full)
        print(f"  {label}: raw={err_raw:+.6f}  frame-fit={err_frame:+.6f} kcal/mol  chi2/dof={chi2dof:.3e}")
        assert abs(err_frame) < 1e-6, f"{label}: analysis pipeline is NOT exact on perfect data -- real bug, stop"
    print("  PASS: analysis pipeline is exact on perfect data for both fragments.")


def submit(resume=True):
    p = setup_fragment(ATOMS_B, nelec=4, d=1.0, K=K, strict=True)
    print(f"  fragment B exact energy: {p['exact_energy']:.6f} Ha "
          f"(fragment A's own: -2.166387448627536 Ha -- must match closely)")
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)

    ancilla_native = ancilla_cnots_native(GATE_NAME)
    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    if resume and os.path.exists(OUT_PATH):
        with open(OUT_PATH) as f:
            state = json.load(f)
        print(f"  resuming: {len(state['done'])} slots already done")
    else:
        state = {"done": {}}

    t0 = time.time()
    for backend_name in BACKENDS:
        for name in kept:
            key = f"{backend_name}|{name}"
            if key in state["done"]:
                print(f"    skip (done): {key}")
                continue
            register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
            circuits = []
            for group in groups:
                combined = effrag_mod.combined_basis_label(group)
                basis_qc = native_basis_change(combined, GATE_NAME)
                full5 = with_ancilla_parity_native(register_native, ancilla_native)
                full5 = full5.compose(basis_qc, qubits=[0, 1, 2, 3])
                full5.measure_all()
                circuits.append(full5)
            job = submit_job(circuits, backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"groups": groups, "counts": counts}
            print(f"    done: {key} ({len(circuits)} circuits, {SHOTS} shots each, {time.time()-t0:.1f}s elapsed)")
            with open(OUT_PATH, "w") as f:
                json.dump(state, f, indent=2)
    print(f"\n  ALL DONE in {time.time()-t0:.1f}s. Saved -> {OUT_PATH}\n")


def analyze():
    with open(OUT_PATH) as f:
        state = json.load(f)
    p = setup_fragment(ATOMS_B, nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    for backend_name in BACKENDS:
        blended = {name: {} for name in kept}
        accept_fracs = []
        for name in kept:
            entry = state["done"][f"{backend_name}|{name}"]
            for group, counts in zip(entry["groups"], entry["counts"]):
                total = sum(counts.values())
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
                accept_fracs.append(sum(filtered.values()) / total)
                for l in group:
                    blended[name][l] = expectation_from_counts(filtered, l)

        if backend_name == "ideal":
            data_for_fit = blended  # NO correction for ideal -- see module docstring
        else:
            data_for_fit = apply_correction(blended, GPI2_SELECTED[backend_name], kept, non_id_labels,
                                             fixed_solutions)

        rng = np.random.default_rng(39)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, data_for_fit,
                                                  weight_unit, rng, n_restarts=4)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        E, err = energy_and_err(p, full)
        print(f"  {backend_name}: mean_accept={np.mean(accept_fracs):.4f}  E={E:.6f} Ha  "
              f"err_vs_exact={err:+.4f} kcal/mol  chi2/dof={chi2dof:.5f}")


if __name__ == "__main__":
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36): the joint-frame "
              "nonconvex fit can land in a different, worse local minimum. Rerun with PYTHONHASHSEED=0.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyze", action="store_true", help="analyze already-submitted data, no network calls")
    ap.add_argument("--self-test", action="store_true", help="synthetic exact-data regression test only")
    args = ap.parse_args()

    if args.self_test:
        synthetic_exact_regression_test()
    elif args.analyze:
        analyze()
    else:
        submit()
        analyze()
