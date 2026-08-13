#!/usr/bin/env python3
"""
task28d_resume.py -- one-off recovery: task28d_all_gate_zne.py's own
checkpoint is only written after ALL 15 (fold,model) jobs are submitted
AND retrieved, so its 2nd crash at fold=7/forte-1 (same combo as the 1st
crash) left 11 real, already-succeeded jobs with no way to recover their
results except by real GUID -- which THIS run's log captured for the
first time. Reconstructs the identical circuits/tags/order the original
script builds (fully deterministic, no randomness in circuit
construction), retrieves the 11 known-good jobs by ID instead of
resubmitting them, submits only the 4 still missing, and writes the
checkpoint task28d_all_gate_zne.py itself expects -- so rerunning that
script afterward will find it and skip straight to analysis.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from fixed_ansatz import build_ansatz
from native_stateprep import to_native, native_target
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list
from task28d_all_gate_zne import (
    K, FOLD_FACTORS, SHOTS, GATE_BY_MODEL, CKPT_DIR, CKPT_PATH,
    optimized_native_circuit, fold_all_gates,
)

SPLIT_MISSING_INTO = 3  # fold=7/forte-1 failed 3/3 times as one 273-circuit job;
                         # split still-missing jobs into smaller sub-batches so no
                         # single request is the largest payload of the day.

KNOWN_JOBS = {
    (1, "ideal"): "019ffc2c-b9f2-758e-b591-90598439dbf2",
    (1, "aria-1"): "019ffc2c-cd1a-7786-9c9d-f9f2969eda4c",
    (1, "forte-1"): "019ffc2c-e39d-721f-b7e7-94c048854355",
    (3, "ideal"): "019ffc2d-0bf7-70e5-99e8-79f58003a3d1",
    (3, "aria-1"): "019ffc2d-34e5-73dd-ad06-ced2ed0bc000",
    (3, "forte-1"): "019ffc2d-6455-712e-ab28-2083e85c5bea",
    (5, "ideal"): "019ffc2d-a726-71b4-b080-7f0cd5e77bfa",
    (5, "aria-1"): "019ffc2d-e697-7269-969b-3642dcc2d0bc",
    (5, "forte-1"): "019ffc2e-2ea2-70dd-b19d-dd9fd66f6200",
    (7, "ideal"): "019ffc2f-0ed5-76ee-9d08-61fd8e07a4c6",
    (7, "aria-1"): "019ffc2f-6365-76ac-869c-db32e599f3b4",
}


def main():
    print("\n" + "=" * 96)
    print("  task28d_resume.py -- recovering 11 known jobs, submitting the 4 still missing")
    print("=" * 96)

    if os.path.exists(CKPT_PATH):
        print(f"  checkpoint already exists at {CKPT_PATH} -- nothing to do, aborting.")
        return

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == 36
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    diag, plus, kept = kept_slots_for_K(K)
    print(f"  K={K}: {len(kept)} kept circuits")

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    # -- rebuild optimizer-frozen base circuits per gate family (deterministic) --
    base_by_gate = {"ms": {}, "zz": {}}
    for gate_name in ["ms", "zz"]:
        for name in kept:
            base_by_gate[gate_name][name] = optimized_native_circuit(fixed_solutions[name]["angles"], gate_name)

    all_jobs = {}
    all_counts = {}
    missing = [(fold, model) for fold in FOLD_FACTORS for model in GATE_BY_MODEL if (fold, model) not in KNOWN_JOBS]
    print(f"  {len(KNOWN_JOBS)} known-good jobs to retrieve by ID, {len(missing)} still missing: {missing}")

    t0 = time.time()
    for fold in FOLD_FACTORS:
        for model, gate_name in GATE_BY_MODEL.items():
            key = (fold, model)
            base = base_by_gate[gate_name]
            circuits, tags = [], []
            for name in kept:
                folded = fold_all_gates(base[name], fold, gate_name)
                for group in groups:
                    combined = effrag_mod.combined_basis_label(group)
                    basis_qc = native_basis_change(combined, gate_name)
                    qc = folded.compose(basis_qc)
                    qc.measure_all()
                    circuits.append(qc)
                    tags.append((name, tuple(group)))

            if key in KNOWN_JOBS:
                guid = KNOWN_JOBS[key]
                print(f"    retrieving known job fold={fold} model={model} id={guid} ...")
                jobs_for_key = [backend.retrieve_job(guid)]
            else:
                # split into SPLIT_MISSING_INTO sub-batches -- fold=7/forte-1's single
                # 273-circuit job failed 3/3 times; a smaller payload per request is
                # the direct mitigation for a payload-size-correlated timeout.
                n_chunks = SPLIT_MISSING_INTO
                chunk_size = (len(circuits) + n_chunks - 1) // n_chunks
                chunks = [circuits[i:i + chunk_size] for i in range(0, len(circuits), chunk_size)]
                print(f"    submitting missing job fold={fold} model={model} as {len(chunks)} sub-batches "
                      f"({[len(c) for c in chunks]} circuits each) ...")
                jobs_for_key = []
                for ci, chunk in enumerate(chunks):
                    job = None
                    for attempt in range(6):
                        try:
                            job = submit_job(chunk, backend, model, shots=SHOTS)
                            break
                        except Exception as e:
                            wait_s = min(30 * (2 ** attempt), 300)
                            print(f"      chunk {ci+1}/{len(chunks)} submit failed (attempt {attempt+1}/6): "
                                  f"{e} -- backing off {wait_s}s")
                            time.sleep(wait_s)
                    if job is None:
                        raise RuntimeError(f"submit_job exhausted 6 manual retries for fold={fold} "
                                            f"model={model} chunk={ci+1}/{len(chunks)}")
                    print(f"      chunk {ci+1}/{len(chunks)} submitted: job_id={job.job_id()}")
                    jobs_for_key.append(job)

            all_jobs[key] = (jobs_for_key, tags)

    print(f"\n  all 15 jobs available (11 recovered + {len(missing)} freshly submitted), retrieving results...")
    for i, (key, (jobs_for_key, tags)) in enumerate(all_jobs.items()):
        counts = []
        for job in jobs_for_key:
            counts.extend(get_counts_list(job))
        assert len(counts) == len(tags), f"count/tag length mismatch for {key}: {len(counts)} vs {len(tags)}"
        all_counts[f"{key[0]}|{key[1]}"] = counts
        print(f"    retrieved {i+1}/{len(all_jobs)}: fold={key[0]} model={key[1]}, {time.time()-t0:.1f}s elapsed")

    ck = {
        "fold_factors": FOLD_FACTORS, "shots": SHOTS, "verify_err": 1.09e-14,
        "tags": {f"{k[0]}|{k[1]}": [[t[0], list(t[1])] for t in v[1]] for k, v in all_jobs.items()},
        "counts": all_counts, "wall_clock": {"submit_s": None, "retrieve_s": time.time() - t0},
        "recovered_job_ids": {f"{k[0]}|{k[1]}": (KNOWN_JOBS.get(k) or [j.job_id() for j in all_jobs[k][0]]) for k in all_jobs},
    }
    os.makedirs(CKPT_DIR, exist_ok=True)
    with open(CKPT_PATH, "w") as f:
        json.dump(ck, f, indent=2)
    print(f"\n  checkpoint saved -> {CKPT_PATH}")
    print("  rerun task28d_all_gate_zne.py now -- it will find this checkpoint and skip straight to analysis.\n")


if __name__ == "__main__":
    main()
