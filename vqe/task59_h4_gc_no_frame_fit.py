#!/usr/bin/env python3
"""
task59_h4_gc_no_frame_fit.py -- iteration 59. Real H4 no-frame-fit
measurement using the LOWEST-CIRCUIT-COUNT design available (general-
commuting measurement grouping, 4 groups instead of the original 13 QWC
groups -- 84 circuits/backend instead of 273) plus the already-reduced
native gate count (optimized_native_circuit, IonQ's own
TrappedIonOptimizerPlugin applied). This is a full, complete 21-slot
sweep -- the earlier GC real-hardware check (task49) only ever covered 3
slots as a physical-correctness spot-check, never a full energy
measurement. Every circuit is still built by fitting angles to the
EXACT classical Schmidt vectors -- this is the same oracle-informed
circuit construction task40's ablation and Vadim's Sep 22 email both
already account for; only the MEASUREMENT GROUPING changes here, purely
a real, already-validated efficiency improvement, unrelated to the
oracle-basis question.

Reconstruction: same as task40's own no-frame-fit ablation --
QED (ancilla-parity postselection) + conditioned PEC + GPi2 correction,
raw corrected Pauli values placed DIRECTLY into the K x K matrices (no
P_S, no frame fit, no reference to U_exact anywhere in reconstruction).

Free simulators only (ideal/aria-1/forte-1), no real hardware spend.

Run:
    PYTHONHASHSEED=0 python vqe/task59_h4_gc_no_frame_fit.py --canary
    PYTHONHASHSEED=0 python vqe/task59_h4_gc_no_frame_fit.py
    PYTHONHASHSEED=0 python vqe/task59_h4_gc_no_frame_fit.py --analyze
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native, verify_native_ancilla
from general_commuting_measurements import build_general_commuting_measurement_plan, expectations_from_counts
from native_stateprep import to_native
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list
from task37c_extended_forward_model import analytic_A_and_B_5param
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from task30b_pec_application import build_full
from qiskit.circuit import QuantumCircuit

HARTREE_TO_KCAL_MOL = 627.5094740631
GATE_NAME = "zz"
K = 6
ZZ_ASSUMED = 0.014593
GPI2_SELECTED = {"aria-1": 0.0006, "forte-1": 0.0004}
SHOTS = 20_000
BACKENDS = ["ideal", "aria-1", "forte-1"]
CKPT_PATH = os.path.join(os.path.dirname(__file__), "task59_h4_gc_no_frame_fit.partial.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task59_h4_gc_no_frame_fit_results.json")


def build_gc_measurement_circuits(register_native_4q, ancilla_native_5q, diagonalizers, diag_natives):
    circuits = []
    for diag, diag_native in zip(diagonalizers, diag_natives):
        full5 = with_ancilla_parity_native(register_native_4q, ancilla_native_5q)
        # compose the (already-native) diagonalizer onto qubits [0,1,2,3] of
        # the 5-qubit (4 register + 1 ancilla) circuit, measurement LAST --
        # same native-then-measure ordering fix verified in task58.
        qc = full5.compose(diag_native, qubits=[0, 1, 2, 3])
        qc.measure_all()
        circuits.append(qc)
    return circuits


def load_partial():
    if os.path.exists(CKPT_PATH):
        with open(CKPT_PATH) as f:
            return json.load(f)
    return {"done": {}}


def save_partial(state):
    with open(CKPT_PATH, "w") as f:
        json.dump(state, f, indent=2)


def setup():
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    diag_natives = [to_native(dg.to_circuit(), GATE_NAME) for dg in diagonalizers]
    print(f"  {len(non_id_labels)} labels -> {len(groups)} GC groups -> "
          f"{len(kept)*len(groups)} circuits/backend (vs 273 for the original QWC design)")
    return p, non_id_labels, fixed_solutions, kept, groups, diagonalizers, diag_natives


def submit(args):
    p, non_id_labels, fixed_solutions, kept, groups, diagonalizers, diag_natives = setup()

    ancilla_native = ancilla_cnots_native(GATE_NAME)
    diff = verify_native_ancilla(ancilla_native)
    print(f"  ancilla circuit regression check: diff={diff:.3e}  {'PASS' if diff < 1e-8 else 'FAIL -- STOP'}")
    assert diff < 1e-8, "ancilla circuit failed its own regression check -- stop before submitting"

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    backends_to_run = ["ideal"] if args.canary else BACKENDS
    slots_to_run = [kept[0]] if args.canary else kept

    state = load_partial()
    for backend_name in backends_to_run:
        for name in slots_to_run:
            key = f"{backend_name}|{name}"
            if key in state["done"]:
                print(f"    skip (already done): {key}")
                continue
            register_native = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
            circuits = build_gc_measurement_circuits(register_native, ancilla_native, diagonalizers, diag_natives)
            job = submit_job(circuits, backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"counts": counts}
            save_partial(state)
            print(f"    done: {key} ({len(circuits)} circuits, {SHOTS} shots each)")

    print(f"\n  {'CANARY PASSED' if args.canary else 'FULL SWEEP COMPLETE'} -- saved -> {CKPT_PATH}")
    return p, non_id_labels, fixed_solutions, kept, diagonalizers


_CACHE = {}
def get_ab(name, p_gpi2, fixed_solutions, non_id_labels):
    key = (name, round(p_gpi2, 6))
    if key not in _CACHE:
        A, B, _retA, _retB = analytic_A_and_B_conditioned(fixed_solutions[name]["angles"], GATE_NAME, ZZ_ASSUMED,
                                                             GPI_REAL_MEAN, p_gpi2, 0.0, 0.0, non_id_labels)
        _CACHE[key] = (A, B)
    return _CACHE[key]


def analyze(p, non_id_labels, fixed_solutions, kept, diagonalizers):
    with open(CKPT_PATH) as f:
        state = json.load(f)
    diag, plus, _ = kept_slots_for_K(K)

    print(f"\n  H4 exact_energy={p['exact_energy']:.6f} Ha")
    results = {}
    for backend_name in BACKENDS:
        # postselected (ancilla=0) raw expectations, via GC decode
        blended = {name: {} for name in kept}
        accept_fracs = []
        for name in kept:
            entry = state["done"][f"{backend_name}|{name}"]
            for dg, counts in zip(diagonalizers, entry["counts"]):
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
                total = sum(counts.values())
                accept_fracs.append(sum(filtered.values()) / total if total else 0.0)
                exp = expectations_from_counts(filtered, dg)
                blended[name].update(exp)

        if backend_name == "ideal":
            raw_kept = {name: {l: max(-1.0, min(1.0, v)) for l, v in blended[name].items()} for name in kept}
        else:
            p_gpi2_sel = GPI2_SELECTED[backend_name]
            raw_kept = {}
            for name in kept:
                raw_kept[name] = {}
                A, B = get_ab(name, p_gpi2_sel, fixed_solutions, non_id_labels)
                for l, m_raw in blended[name].items():
                    ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
                    raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))

        full = build_full(raw_kept, diag, K, non_id_labels)
        from qforge import energy_from_alpha_matrices
        # build_full already gives alpha-matrix-shaped dict; reuse combine path via qforge
        from qforge import combine_matrices
        alpha_mats = combine_matrices(full, p["alpha_labels"], p["identity_label"], K)
        E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                              exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
        print(f"    {backend_name}: mean_accept={np.mean(accept_fracs):.4f}  E={E:.6f} Ha  "
              f"err_vs_exact={errs['err_vs_exact_kcal']:+.4f} kcal/mol")
        results[backend_name] = {"E": E, "err_kcal": errs["err_vs_exact_kcal"], "mean_accept": float(np.mean(accept_fracs))}

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved -> {RESULTS_PATH}")
    return results


if __name__ == "__main__":
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    args = ap.parse_args()

    if args.analyze:
        p, non_id_labels, fixed_solutions, kept, groups, diagonalizers, diag_natives = setup()
        analyze(p, non_id_labels, fixed_solutions, kept, diagonalizers)
    else:
        p, non_id_labels, fixed_solutions, kept, diagonalizers = submit(args)
        if not args.canary:
            analyze(p, non_id_labels, fixed_solutions, kept, diagonalizers)
