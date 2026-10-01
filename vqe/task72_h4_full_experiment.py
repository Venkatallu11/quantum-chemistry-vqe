#!/usr/bin/env python3
"""
task72_h4_full_experiment.py -- iteration 72. The actual dress-run of the
task71-approved 54-circuit real-hardware plan, executed for real against
IonQ's FREE ionq_simulator (ideal/aria-1/forte-1 noise models) -- zero
real-hardware dollars spent -- before committing the approved hardware
budget. Unlike task71 (which reused task59's OLD checkpoint, bootstrap-
resampled down to the target shot count), this submits FRESH real jobs,
built with exactly the same 54-circuit design (6-slot full-coverage
panel measured through all 4 GC groups, the other 15 slots measured
through GC0+GC1 only), at the real proposed shot count (1,100
shots/circuit), and analyzes the REAL returned counts directly -- no
bootstrap step anywhere in this script.

Reuses, unchanged: task59's native ancilla-parity + ancilla regression
check, task70's GC-aware conditioned correction, task60's per-slot
independent-state fit + build_full_from_independent_states, task36's
joint Schmidt-frame fit. Same three locked analyses as every hardware-
plan rehearsal in this project: RAW (uncorrected), NO-FRAME (corrected,
no shared frame), SHARED-FRAME (corrected, joint 15-parameter frame fit).

Run (resumable -- checkpoints after every slot):
    PYTHONHASHSEED=0 python vqe/task72_h4_full_experiment.py --submit
    PYTHONHASHSEED=0 python vqe/task72_h4_full_experiment.py --analyze
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task28d_all_gate_zne import optimized_native_circuit
from task39b_native_ancilla_parity import ancilla_cnots_native, with_ancilla_parity_native, verify_native_ancilla
from general_commuting_measurements import build_general_commuting_measurement_plan, expectations_from_counts
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from task59_h4_gc_no_frame_fit import GATE_NAME, GPI2_SELECTED, BACKENDS
from task60_ionq_no_frame_h4 import fit_all_slots, build_full_from_independent_states, energy_report, variance_weights
from task70_gc_aware_correction import analytic_A_and_B_conditioned_gc
from task71_hw_plan_rehearsal_v2 import PANEL_SLOTS, measured_groups_for_slot
from task37b_h4_noise_model import GPI_REAL_MEAN
from native_stateprep import to_native
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list

K = 6
ZZ_ASSUMED = 0.014593
PROPOSED_SHOTS = 1100


def ckpt_path(trial):
    suffix = "" if trial == 1 else f"_trial{trial}"
    return os.path.join(os.path.dirname(__file__), f"task72_h4_full_experiment{suffix}.partial.json")


def results_path(trial):
    suffix = "" if trial == 1 else f"_trial{trial}"
    return os.path.join(os.path.dirname(__file__), f"task72_h4_full_experiment{suffix}_results.json")


# kept for backward compatibility with anything importing these directly (trial 1 paths)
CKPT_PATH = ckpt_path(1)
RESULTS_PATH = results_path(1)


def setup():
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    diag_natives = [to_native(dg.to_circuit(), GATE_NAME) for dg in diagonalizers]
    n_groups = len(groups)
    measured_groups_map = {name: measured_groups_for_slot(name, n_groups) for name in kept}
    n_circuits = sum(len(v) for v in measured_groups_map.values())
    print(f"  design: {n_circuits} circuits/backend (panel: {PANEL_SLOTS})")
    return dict(p=p, non_id_labels=non_id_labels, fixed_solutions=fixed_solutions, diag=diag, kept=kept,
                P_S=P_S, groups=groups, diagonalizers=diagonalizers, diag_natives=diag_natives,
                n_groups=n_groups, measured_groups_map=measured_groups_map)


def load_partial(trial):
    path = ckpt_path(trial)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"done": {}}


def save_partial(state, trial):
    with open(ckpt_path(trial), "w") as f:
        json.dump(state, f, indent=2)


def build_circuits_for_slot(register_native, ancilla_native, diagonalizers, diag_natives, group_idxs):
    circuits = []
    for gi in group_idxs:
        full5 = with_ancilla_parity_native(register_native, ancilla_native)
        qc = full5.compose(diag_natives[gi], qubits=[0, 1, 2, 3])
        qc.measure_all()
        circuits.append(qc)
    return circuits


def submit(ctx, args, trial):
    ancilla_native = ancilla_cnots_native(GATE_NAME)
    diff = verify_native_ancilla(ancilla_native)
    print(f"  ancilla circuit regression check: diff={diff:.3e}  {'PASS' if diff < 1e-8 else 'FAIL -- STOP'}")
    assert diff < 1e-8, "ancilla circuit failed its own regression check -- stop before submitting"

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}, shots={PROPOSED_SHOTS}/circuit, trial={trial}")

    backends_to_run = ["ideal"] if args.canary else BACKENDS
    slots_to_run = [ctx["kept"][0]] if args.canary else ctx["kept"]

    state = load_partial(trial)
    for backend_name in backends_to_run:
        for name in slots_to_run:
            key = f"{backend_name}|{name}"
            if key in state["done"]:
                print(f"    skip (already done): {key}")
                continue
            group_idxs = ctx["measured_groups_map"][name]
            register_native = optimized_native_circuit(ctx["fixed_solutions"][name]["angles"], GATE_NAME)
            circuits = build_circuits_for_slot(register_native, ancilla_native, ctx["diagonalizers"],
                                                 ctx["diag_natives"], group_idxs)
            job = submit_job(circuits, backend, backend_name, shots=PROPOSED_SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"counts": counts, "group_idxs": group_idxs}
            save_partial(state, trial)
            print(f"    done: {key} ({len(circuits)} circuits, groups={group_idxs}, {PROPOSED_SHOTS} shots each)")

    print(f"\n  {'CANARY PASSED' if args.canary else 'FULL SWEEP COMPLETE'} -- saved -> {ckpt_path(trial)}")


_AB_CACHE = {}
def corrected_observables_gc_aware(postselected, p_gpi2, kept, fixed_solutions, diag_natives, diagonalizers,
                                     groups, measured_groups_map):
    corrected = {name: {} for name in kept}
    for name in kept:
        for gi in measured_groups_map[name]:
            key = (name, gi, round(float(p_gpi2), 7))
            if key not in _AB_CACHE:
                A, B, _, _ = analytic_A_and_B_conditioned_gc(
                    fixed_solutions[name]["angles"], GATE_NAME, ZZ_ASSUMED, GPI_REAL_MEAN, float(p_gpi2), 0.0, 0.0,
                    diag_natives[gi], diagonalizers[gi].transformed)
                _AB_CACHE[key] = (A, B)
            A, B = _AB_CACHE[key]
            for label in groups[gi]:
                raw = float(postselected[name][label])
                denom = float(A[label])
                ratio = float(B[label] / denom) if abs(denom) > 1e-6 else 1.0
                corrected[name][label] = float(np.clip(raw * ratio, -1.0, 1.0))
    return corrected


def analyze_backend(ctx, backend_name, state):
    kept = ctx["kept"]
    diagonalizers = ctx["diagonalizers"]
    measured_groups_map = ctx["measured_groups_map"]

    postselected = {name: {} for name in kept}
    kept_shots = {name: {} for name in kept}
    accept_fracs = []
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for local_i, gi in enumerate(entry["group_idxs"]):
            dg = diagonalizers[gi]
            full_counts = entry["counts"][local_i]
            filtered = {bs[1:]: c for bs, c in full_counts.items() if bs[0] == "0"}
            total = sum(full_counts.values())
            n_kept = sum(filtered.values())
            accept_fracs.append(n_kept / total if total else 0.0)
            exp = expectations_from_counts(filtered, dg) if filtered else {l: 0.0 for l in ctx["groups"][gi]}
            for l in ctx["groups"][gi]:
                postselected[name][l] = max(-1.0, min(1.0, exp.get(l, 0.0)))
                kept_shots[name][l] = n_kept

    weights_raw = variance_weights(postselected, kept_shots)
    fits_raw = fit_all_slots(postselected, weights_raw, ctx["P_S"], kept)
    full_raw = build_full_from_independent_states(fits_raw, ctx["diag"], K, ctx["P_S"], ctx["non_id_labels"])
    _, errs_raw = energy_report(ctx["p"], full_raw, K)

    if backend_name == "ideal":
        corrected = postselected
    else:
        corrected = corrected_observables_gc_aware(postselected, GPI2_SELECTED[backend_name], kept,
                                                      ctx["fixed_solutions"], ctx["diag_natives"], diagonalizers,
                                                      ctx["groups"], measured_groups_map)
    weights_c = variance_weights(corrected, kept_shots)
    fits_nf = fit_all_slots(corrected, weights_c, ctx["P_S"], kept)
    full_nf = build_full_from_independent_states(fits_nf, ctx["diag"], K, ctx["P_S"], ctx["non_id_labels"])
    _, errs_nf = energy_report(ctx["p"], full_nf, K)

    weight_unit = {name: {l: 1.0 for l in ctx["non_id_labels"]} for name in kept}
    rng_frame = np.random.default_rng(7200)
    U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), ctx["P_S"], K, kept, ctx["non_id_labels"], corrected,
                                             weight_unit, rng_frame, n_restarts=4)
    full_sf = build_full_from_frame(U_hat, ctx["P_S"], K, ctx["non_id_labels"], kept)
    _, errs_sf = energy_report(ctx["p"], full_sf, K)

    return {
        "mean_accept": float(np.mean(accept_fracs)),
        "raw_err_kcal": errs_raw["err_vs_exact_kcal"],
        "no_frame_err_kcal": errs_nf["err_vs_exact_kcal"],
        "shared_frame_err_kcal": errs_sf["err_vs_exact_kcal"],
        "chi2_dof": float(chi2dof),
    }


def _run_analysis_once(ctx, state):
    results = {}
    for backend_name in BACKENDS:
        results[backend_name] = analyze_backend(ctx, backend_name, state)
    return results


def analyze(ctx, trial, verify=True):
    """Computes the three locked analyses from the trial's real checkpoint.

    SAFEGUARD (added after task72's first run printed a number that never
    reproduced again -- root cause never pinned down, see RESEARCH_LEDGER):
    this reloads the checkpoint from disk and recomputes EVERYTHING a
    SECOND time, from a completely independent in-memory state, and
    refuses to report anything unless both passes agree exactly. A live
    run's own immediate printout is never trusted on its own again."""
    path = ckpt_path(trial)
    with open(path) as f:
        state = json.load(f)
    missing = [bn for bn in BACKENDS if any(f"{bn}|{name}" not in state["done"] for name in ctx["kept"])]
    if missing:
        print(f"  ERROR: real checkpoint incomplete for backends {missing} -- run --submit first.")
        return None

    print(f"\n  H4 exact_energy={ctx['p']['exact_energy']:.6f} Ha")
    print(f"  data: REAL fresh submission, trial {trial}, 54-circuit design, 3 backends, "
          f"{PROPOSED_SHOTS} shots/circuit")
    results = _run_analysis_once(ctx, state)

    if verify:
        with open(path) as f:
            state2 = json.load(f)
        results2 = _run_analysis_once(ctx, state2)
        for bn in BACKENDS:
            for k in ("raw_err_kcal", "no_frame_err_kcal", "shared_frame_err_kcal"):
                if abs(results[bn][k] - results2[bn][k]) > 1e-9:
                    raise RuntimeError(
                        f"REPRODUCIBILITY CHECK FAILED for trial {trial}, {bn}.{k}: "
                        f"{results[bn][k]!r} vs {results2[bn][k]!r} on two independent passes over the SAME "
                        f"checkpoint -- refusing to report an unverified number. Do not trust this trial.")
        print("  reproducibility check: two independent passes agree exactly -- PASS")

    for backend_name in BACKENDS:
        entry = results[backend_name]
        print(f"    {backend_name}: accept={entry['mean_accept']:.4f}  "
              f"raw={entry['raw_err_kcal']:+.4f}  no_frame={entry['no_frame_err_kcal']:+.4f}  "
              f"shared_frame={entry['shared_frame_err_kcal']:+.4f} kcal/mol  chi2/dof={entry['chi2_dof']:.4f}")

    with open(results_path(trial), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved -> {results_path(trial)}\n")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--canary", action="store_true")
    ap.add_argument("--trial", type=int, default=1, help="independent trial number (separate checkpoint/results)")
    args = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- rerun with PYTHONHASHSEED=0")

    ctx = setup()
    if args.submit or args.canary:
        submit(ctx, args, args.trial)
        if not args.canary:
            analyze(ctx, args.trial)
    elif args.analyze:
        analyze(ctx, args.trial)
    else:
        print("  pass --submit (full run), --canary (1 slot, ideal only), or --analyze (add --trial N for trial N)")


if __name__ == "__main__":
    main()
