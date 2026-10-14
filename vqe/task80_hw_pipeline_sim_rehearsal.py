#!/usr/bin/env python3
"""
task80_hw_pipeline_sim_rehearsal.py -- iteration 80. End-to-end rehearsal of
the exact hardware pipeline (task74 build + checks -> submit -> fetch ->
task73 locked analysis) on IonQ's FREE ionq_simulator, once per noise model
(ideal, aria-1, forte-1). NOT hardware, no cost.

What is identical to the planned Forte Enterprise 1 run:
  - the 81 new circuits, built and checked by task74.prepare("full")
    (all 6 pre-submission checks, including the debiasing payload gate)
  - one bundled job, with task74.RUN_OPTIONS (1,100 shots, debiasing off)
  - fetch-side assert that every circuit's counts total 1,100
  - counts written in task73's hardware format, analysed by task73 --counts
    (so the 3 reused task49 GC3 hardware circuits are filled in exactly as
    they will be for the real run, and the locked forte-1 GPi2 correction
    is applied to every noise model, as it will be to the hardware data)
What differs: backend = ionq_simulator with noise_model set.

Each noise model's job id is checkpointed before waiting, so a rerun
re-fetches instead of resubmitting.

Run:
    PYTHONHASHSEED=0 python vqe/task80_hw_pipeline_sim_rehearsal.py
    PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py --counts vqe/task80_sim_<model>_counts.json \
        --out vqe/task80_sim_<model>_analysis.json
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
import task74_hw_submission as t74
from ionq_backend import connect_provider, get_native_simulator

HERE = os.path.dirname(__file__)
NOISE_MODELS = ["ideal", "aria-1", "forte-1"]
CKPT_PATH = os.path.join(HERE, "task80_sim_jobs.json")


def counts_path(model):
    return os.path.join(HERE, f"task80_sim_{model}_counts.json")


def main():
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0")
    ctx, circuits, submit_keys, counts, est = t74.prepare("full", t74.BUDGET_TOTAL)
    print(f"  (hardware estimate for reference only: ${est:,.2f}; this rehearsal is free)")

    backend = get_native_simulator(connect_provider())
    assert backend.name == "ionq_simulator", f"refusing to run on {backend.name} -- simulator only"
    print(f"  connected to {backend.name}")

    jobs = {}
    if os.path.exists(CKPT_PATH):
        with open(CKPT_PATH) as f:
            jobs = json.load(f)

    circuit_list = [circuits[k] for k in submit_keys]
    for model in NOISE_MODELS:
        if os.path.exists(counts_path(model)):
            print(f"  skip (counts already saved): {model}")
            continue
        if model in jobs:
            job = backend.retrieve_job(jobs[model])
            print(f"  {model}: re-fetching job {jobs[model]}")
        else:
            job = backend.run(circuit_list, noise_model=model, **t74.RUN_OPTIONS)
            jobs[model] = job.job_id()
            with open(CKPT_PATH, "w") as f:
                json.dump(jobs, f, indent=2)
            print(f"  {model}: submitted {len(circuit_list)} circuits as one job, id={jobs[model]}")
        job.wait_for_final_state(timeout=4 * 3600)
        assert job.done(), f"{model}: job not done: {job.status()}"
        counts_list = job.result().get_counts()
        if not isinstance(counts_list, list):
            counts_list = [counts_list]
        assert len(counts_list) == len(submit_keys), f"{model}: {len(counts_list)} results for {len(submit_keys)} circuits"
        totals = {k: sum(c.values()) for k, c in zip(submit_keys, counts_list)}
        bad = {k: t for k, t in totals.items() if t != t74.SHOTS}
        assert not bad, f"{model}: per-circuit shot totals != {t74.SHOTS}: {bad}"
        raw = backend.client.retrieve_job(jobs[model])
        print(f"  {model}: done, all {len(counts_list)} circuits have {t74.SHOTS} shots; "
              f"API settings={raw.get('settings')}")
        out = {"job_id": jobs[model], "backend": f"ionq_simulator/{model}", "shots": t74.SHOTS,
               "design": "full", "is_hardware": False, "submitted_at": raw.get("submitted_at"),
               "completed_at": raw.get("completed_at"), "billed_cost_usd": raw.get("cost"),
               "counts": {k: dict(c) for k, c in zip(submit_keys, counts_list)},
               "raw_job_metadata": raw}
        with open(counts_path(model), "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"  saved -> {counts_path(model)}")


if __name__ == "__main__":
    main()
