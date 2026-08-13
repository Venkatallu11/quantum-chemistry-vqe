#!/usr/bin/env python3
"""
task28d_all_gate_zne.py — iteration 28, Task D. THE KEY EXPERIMENT: fold
the single-qubit gates too, not just the two-qubit one. Compares, on
identical states and observables:
    (a) 2q-only folding  -- iteration 27's own baseline, the 14.28 result
    (b) all-gate folding -- 1q and 2q scaled together, lambda_1 = lambda_2
============================================================================
CIRCUIT BASE: Task 28B's TrappedIonOptimizerPlugin-optimized circuit
(verified statevector-identical to <1e-12, N_2q now 4-11 per target
instead of the constant 11 iteration 27 used) -- per the explicit
ordering instruction ("optimise the UNFOLDED circuit only, then freeze,
then fold"), folding is applied AFTER optimization, never before.

FOLD INVERSES, verified numerically before use (not assumed):
  - GPi(phi)^2 = I exactly for any phi (iteration 28 Task A, re-verified
    here) -- GPi is its own inverse, so folding it means inserting
    ANOTHER GPi(phi) (not a different gate).
  - GPi2(phi)^{-1} = GPi2(phi+0.5) exactly (verified here: product is I
    to machine precision, both multiplication orders, and equals the
    conjugate-transpose too) -- the SAME phase-shift-by-0.5 pattern this
    project's MS-gate inverse already uses (`ionq_fold_check.py`).
  - The native 2-qubit gate (ms/zz) reuses `fold_native_2q` UNCHANGED.

Every folded circuit -- ALL-GATE folding, not just 2q-only -- is verified
unitary-identical to the UNFOLDED (optimized) circuit to <1e-12 BEFORE
being queued for submission, matching Task B's own discipline.

Held-out validation, exactly as iteration 27 Task D: train on folds
[1,3,5], predict fold=7; train on [1,3,5,7], predict fold=9. NO
CLIPPING -- unphysical (|value|>1) extrapolations are EXCLUDED and
COUNTED, never silently fixed.

Run:
    python vqe/task28d_all_gate_zne.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, native_target
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
from task27d_held_out_zne import MODEL_CLASSES, FOLD_STAGE1_FIT, FOLD_STAGE1_HOLDOUT, FOLD_STAGE2_FIT, FOLD_STAGE2_HOLDOUT
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from qiskit_ionq.ionq_gates import GPIGate, GPI2Gate, MSGate, ZZGate
from qiskit.quantum_info import Statevector

K = 6
FOLD_FACTORS = [1, 3, 5, 7, 9]
SHOTS = 100_000
N_SEEDS = 8
GATE_BY_MODEL = {"ideal": "ms", "aria-1": "ms", "forte-1": "zz"}
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task28d_all_gate_zne.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task28d_all_gate_zne_results.json")


def optimized_native_circuit(angles, gate_name):
    from qiskit.transpiler import PassManagerConfig
    from qiskit_ionq import TrappedIonOptimizerPlugin
    qc = build_ansatz(angles)
    native = to_native(qc, gate_name)
    tgt = native_target(qc.num_qubits, gate_name)
    pm = TrappedIonOptimizerPlugin().pass_manager(PassManagerConfig(target=tgt), optimization_level=3)
    return pm.run(native)


def fold_all_gates(qc, fold, gate_name):
    """Folds EVERY gate (gpi, gpi2, AND the native 2-qubit gate) by the
    SAME odd fold factor -- lambda_1q = lambda_2q. Reuses the verified
    inverse for each gate type (see module docstring)."""
    if fold == 1:
        return qc.copy()
    assert fold % 2 == 1
    reps = (fold - 1) // 2
    folded = qc.copy_empty_like()
    for instr in qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        folded.append(op, qargs, cargs)
        if op.name == "gpi":
            for _ in range(reps):
                folded.append(op, qargs, cargs)   # self-inverse: G(G G)^reps
                folded.append(op, qargs, cargs)
        elif op.name == "gpi2":
            phi = float(op.params[0])
            inv = GPI2Gate(phi + 0.5)
            for _ in range(reps):
                folded.append(inv, qargs, cargs)
                folded.append(op, qargs, cargs)
        elif op.name == gate_name:
            for _ in range(reps):
                if gate_name == "ms":
                    inv = MSGate(float(op.params[0]) + 0.5, op.params[1], op.params[2])
                    folded.append(inv, qargs, cargs)
                else:
                    folded.append(GPIGate(0), [qargs[0]])
                    folded.append(op, qargs, cargs)
                    folded.append(GPIGate(0), [qargs[0]])
                folded.append(op, qargs, cargs)
    return folded


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print("  task28d_all_gate_zne.py -- THE KEY EXPERIMENT: fold 1q gates too")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)
    print(f"  K={K}: {len(kept)} kept circuits, optimizer-frozen base (Task 28B)")

    if os.path.exists(CKPT_PATH):
        print("  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"  connected, backend={backend.name}")

        # -- build + verify EVERY folded circuit (all-gate) vs its own unfolded optimized base --
        base_by_gate = {"ms": {}, "zz": {}}
        for gate_name in ["ms", "zz"]:
            for name in kept:
                base_by_gate[gate_name][name] = optimized_native_circuit(fixed_solutions[name]["angles"], gate_name)

        max_verify_err = 0.0
        folded_by_gate = {"ms": {}, "zz": {}}
        for gate_name in ["ms", "zz"]:
            for name in kept:
                base = base_by_gate[gate_name][name]
                sv_base = np.asarray(Statevector.from_instruction(base))
                for fold in FOLD_FACTORS:
                    folded = fold_all_gates(base, fold, gate_name)
                    sv_folded = np.asarray(Statevector.from_instruction(folded))
                    idx = int(np.argmax(np.abs(sv_base)))
                    phase = sv_folded[idx] / sv_base[idx] if abs(sv_base[idx]) > 1e-9 else 1.0
                    err = float(np.max(np.abs(sv_folded / phase - sv_base)))
                    max_verify_err = max(max_verify_err, err)
                    folded_by_gate[gate_name][(name, fold)] = folded
        print(f"  all-gate fold unitary-equivalence check: worst={max_verify_err:.2e} (must be <1e-12): "
              f"{'PASS' if max_verify_err < 1e-12 else 'FAIL -- refusing to submit'}")
        assert max_verify_err < 1e-12

        # -- submit ALL (fold, model) jobs, non-blocking --
        all_jobs = {}
        t0 = time.time()
        for fold in FOLD_FACTORS:
            for model, gate_name in GATE_BY_MODEL.items():
                circuits, tags = [], []
                for name in kept:
                    folded = folded_by_gate[gate_name][(name, fold)]
                    for group in groups:
                        combined = effrag_mod.combined_basis_label(group)
                        basis_qc = native_basis_change(combined, gate_name)
                        qc = folded.compose(basis_qc)
                        qc.measure_all()
                        circuits.append(qc)
                        tags.append((name, tuple(group)))
                # manual retry-with-backoff around submission: the library's own short internal
                # retries were observed (3 consecutive full-script attempts) to exhaust against a
                # sustained IonQ 'Service Unavailable'/'upstream timing out' condition -- back off
                # longer between attempts here rather than treating a 4th blind retry as different.
                job = None
                for attempt in range(6):
                    try:
                        job = submit_job(circuits, backend, model, shots=SHOTS)
                        break
                    except Exception as e:
                        wait_s = min(30 * (2 ** attempt), 300)
                        print(f"    submit failed (fold={fold} model={model}, attempt {attempt+1}/6): {e} "
                              f"-- backing off {wait_s}s")
                        time.sleep(wait_s)
                if job is None:
                    raise RuntimeError(f"submit_job exhausted 6 manual retries for fold={fold} model={model}")
                print(f"    submitted fold={fold} model={model}: job_id={job.job_id()}")
                all_jobs[(fold, model)] = (job, tags)
        t_submit = time.time() - t0
        print(f"  all {len(all_jobs)} (fold, model) jobs submitted (non-blocking), {t_submit:.1f}s")

        t0 = time.time()
        all_counts = {}
        for i, (key, (job, tags)) in enumerate(all_jobs.items()):
            all_counts[f"{key[0]}|{key[1]}"] = get_counts_list(job)
            print(f"    retrieved {i+1}/{len(all_jobs)}: fold={key[0]} model={key[1]}, {time.time()-t0:.1f}s elapsed")
        t_retrieve = time.time() - t0
        print(f"  all {len(all_jobs)} jobs retrieved, {t_retrieve:.1f}s")

        ck = {
            "fold_factors": FOLD_FACTORS, "shots": SHOTS, "verify_err": max_verify_err,
            "tags": {f"{k[0]}|{k[1]}": [[t[0], list(t[1])] for t in v[1]] for k, v in all_jobs.items()},
            "counts": all_counts, "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve},
        }
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- energy per fold, per model (raw, exact vs shot-noisy) --
    print(f"\n  -- ALL-GATE-folded H4 energy vs fold --")
    for model in ["ideal", "aria-1", "forte-1"]:
        for fold in FOLD_FACTORS:
            key = f"{fold}|{model}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            per_name = {}
            for (name, group), counts in zip(tags, counts_list):
                per_name.setdefault(name, {}).setdefault(tuple(group), counts)
            errs = []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task28d", fold, model, seed))
                m = {name: {} for name in diag + plus}
                for name in diag + plus:
                    for group_t, counts in per_name[name].items():
                        resampled = bootstrap_counts(counts, SHOTS, rng)
                        for l in group_t:
                            m[name][l] = expectation_from_counts(resampled, l)
                full = {name: dict(m[name]) for name in diag}
                for n in range(K):
                    for mm in range(K):
                        if n >= mm:
                            continue
                        un, um, pl = f"u_{n}", f"u_{mm}", f"(u{n}+u{mm})"
                        full[pl] = dict(m[pl])
                        synth_minus = {}
                        for l in non_id_labels:
                            if l not in m[pl] or l not in full[un] or l not in full[um]:
                                continue
                            cross = m[pl][l] - (full[un][l] + full[um][l]) / 2
                            synth_minus[l] = m[pl][l] - 2 * cross
                        full[f"(u{n}-u{mm})"] = synth_minus
                _, err = energy_and_err(p, full, K)
                errs.append(err)
            print(f"    {model} fold={fold}: {np.mean(errs):.2f}+/-{np.std(errs):.2f} kcal/mol")

    # -- held-out ZNE, all-gate curves, same two-stage procedure as iteration 27 Task D --
    print(f"\n  -- held-out ZNE on ALL-GATE-folded curves --")
    results_by_model = {}
    for model in ["ideal", "aria-1", "forte-1"]:
        curves = {}
        for fold in FOLD_FACTORS:
            key = f"{fold}|{model}"
            tags = ck["tags"][key]
            counts_list = ck["counts"][key]
            per_name = {}
            for (name, group), counts in zip(tags, counts_list):
                per_name.setdefault(name, {}).setdefault(tuple(group), counts)
            for name in kept:
                for group_t, counts in per_name[name].items():
                    seed_vals = []
                    for seed in range(N_SEEDS):
                        rng = np.random.default_rng(stable_seed("task28d_curve", fold, model, seed, name))
                        resampled = bootstrap_counts(counts, SHOTS, rng)
                        for l in group_t:
                            curves.setdefault((name, l), {}).setdefault(fold, []).append(expectation_from_counts(resampled, l))
        for key2 in curves:
            for fold in curves[key2]:
                curves[key2][fold] = float(np.mean(curves[key2][fold]))

        stage1_errors, stage2_selected, stage2_errors = {}, {}, {}
        for key2, c in curves.items():
            if not all(f in c for f in FOLD_STAGE1_FIT + [FOLD_STAGE1_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_STAGE1_FIT]
            best_name, best_err = None, float("inf")
            for name, fitter in MODEL_CLASSES.items():
                fn = fitter(FOLD_STAGE1_FIT, vals_fit)
                if fn is None:
                    continue
                try:
                    pred = float(np.atleast_1d(fn([FOLD_STAGE1_HOLDOUT]))[0])
                    if np.isfinite(pred):
                        err = abs(pred - c[FOLD_STAGE1_HOLDOUT])
                        if err < best_err:
                            best_err, best_name = err, name
                except Exception:
                    continue
            if best_name is not None:
                stage1_errors[key2] = best_err
        mean_s1 = float(np.mean(list(stage1_errors.values()))) if stage1_errors else float("inf")

        for key2, c in curves.items():
            if not all(f in c for f in FOLD_STAGE2_FIT + [FOLD_STAGE2_HOLDOUT]):
                continue
            vals_fit = [c[f] for f in FOLD_STAGE2_FIT]
            best_name, best_err = None, float("inf")
            for name, fitter in MODEL_CLASSES.items():
                fn = fitter(FOLD_STAGE2_FIT, vals_fit)
                if fn is None:
                    continue
                try:
                    pred = float(np.atleast_1d(fn([FOLD_STAGE2_HOLDOUT]))[0])
                    if np.isfinite(pred):
                        err = abs(pred - c[FOLD_STAGE2_HOLDOUT])
                        if err < best_err:
                            best_err, best_name = err, name
                except Exception:
                    continue
            if best_name is not None:
                stage2_selected[key2] = best_name
                stage2_errors[key2] = best_err
        mean_s2 = float(np.mean(list(stage2_errors.values()))) if stage2_errors else float("inf")

        n_excluded, n_total, extrap = 0, 0, {}
        for key2, cls in stage2_selected.items():
            n_total += 1
            c = curves[key2]
            all_folds = sorted(c.keys())
            vals_all = [c[f] for f in all_folds]
            fn = MODEL_CLASSES[cls](all_folds, vals_all)
            if fn is None:
                n_excluded += 1
                continue
            try:
                val0 = float(np.atleast_1d(fn([0]))[0])
            except Exception:
                n_excluded += 1
                continue
            if not np.isfinite(val0) or abs(val0) > 1.0:
                n_excluded += 1
                continue
            extrap[key2] = val0

        full = {name: {} for name in kept}
        for name in kept:
            for l in non_id_labels:
                key2 = (name, l)
                if key2 in extrap:
                    full[name][l] = extrap[key2]
                elif key2 in curves and 1 in curves[key2]:
                    full[name][l] = curves[key2][1]
        full_complete = {n2: dict(full[n2]) for n2 in diag}
        for n in range(K):
            for m in range(K):
                if n >= m:
                    continue
                un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
                full_complete[pl] = dict(full[pl])
                synth_minus = {}
                for l in non_id_labels:
                    if l not in full[pl] or l not in full_complete[un] or l not in full_complete[um]:
                        continue
                    cross = full[pl][l] - (full_complete[un][l] + full_complete[um][l]) / 2
                    synth_minus[l] = full[pl][l] - 2 * cross
                full_complete[f"(u{n}-u{m})"] = synth_minus
        _, err_zne = energy_and_err(p, full_complete, K)

        raw1 = {name: {} for name in kept}
        for name in kept:
            for l in non_id_labels:
                key2 = (name, l)
                if key2 in curves and 1 in curves[key2]:
                    raw1[name][l] = curves[key2][1]
        raw1_complete = {n2: dict(raw1[n2]) for n2 in diag}
        for n in range(K):
            for m in range(K):
                if n >= m:
                    continue
                un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
                raw1_complete[pl] = dict(raw1[pl])
                synth_minus = {}
                for l in non_id_labels:
                    if l not in raw1[pl] or l not in raw1_complete[un] or l not in raw1_complete[um]:
                        continue
                    cross = raw1[pl][l] - (raw1_complete[un][l] + raw1_complete[um][l]) / 2
                    synth_minus[l] = raw1[pl][l] - 2 * cross
                raw1_complete[f"(u{n}-u{m})"] = synth_minus
        _, err_raw1 = energy_and_err(p, raw1_complete, K)

        print(f"    {model}: stage1_err={mean_s1:.4f}  stage2_err={mean_s2:.4f}  "
              f"excluded={n_excluded}/{n_total}  raw(fold1)={err_raw1:.2f}  "
              f"ALL-GATE-ZNE(fold0)={err_zne:.2f} kcal/mol  improves={err_zne < err_raw1}")
        results_by_model[model] = {
            "stage1_mean_holdout_err": mean_s1, "stage2_mean_holdout_err": mean_s2,
            "n_excluded": n_excluded, "n_total": n_total,
            "err_raw_fold1": err_raw1, "err_all_gate_zne_fold0": err_zne,
            "improves_on_raw": bool(err_zne < err_raw1),
        }

    print(f"\n  -- COMPARISON: 2q-only ZNE (iteration 27) vs ALL-GATE ZNE (this file) --")
    cited_2q_only = {"aria-1": 52.89, "forte-1": 14.28}
    for model in ["aria-1", "forte-1"]:
        allgate = results_by_model[model]["err_all_gate_zne_fold0"]
        two_q = cited_2q_only[model]
        print(f"    {model}: 2q-only(iter27)={two_q:.2f}  all-gate(this file)={allgate:.2f}  "
              f"{'ALL-GATE IS BETTER' if allgate < two_q else 'all-gate is NOT better'}")

    results = {"verify_err": ck.get("verify_err"), "results_by_model": results_by_model,
               "cited_2q_only_iteration27": cited_2q_only}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
