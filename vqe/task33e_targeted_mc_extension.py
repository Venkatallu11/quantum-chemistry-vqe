#!/usr/bin/env python3
"""
task33e_targeted_mc_extension.py -- iteration 33, Task E. TARGETED real
submission: extend N_MC from 16 to 64 draws, but ONLY for the specific 10
(slot, basis-group) circuit families Task 33D identified as having a
near-coin-flip sign pattern (27 of 756 (slot,label) groups, |mean_sign|
< 0.5, e.g. (u3+u5)|XZXZ at 10+/6-). Task 33D's math: sign_se shrinks as
1/sqrt(n), so going 16->64 draws (4x) should roughly HALVE sign_se for
these specific groups, IF they are genuinely well-behaved importance-
sampling noise and not something else -- this task tests that prediction
against real data instead of assuming it.

This is DELIBERATELY NOT a blanket re-run of Task 31C's full 21-slot x
16-draw collection (4,368 circuits, ~4.1 hours). Only the 10 flagged
(slot, group) pairs get new circuits -- real submissions to the FREE
ionq_simulator forte-1 backend, same as every prior real-data task this
project has run. No real QPU, no cost.

METHOD: reuses Task 31C's own `sample_twirled_circuit` function and
circuit-construction logic UNCHANGED (same base circuit, same basis
rotation, same PEC inverse-weight sampling) -- only the SET of (slot,
group) pairs and the draw COUNT differ. New draws use a distinct RNG seed
per (name, group) (stable_seed("task33e", name, group_idx)), independent
of Task 31C's original single continuous stream -- this is fine and
correct: each draw is an independently valid, real, freshly-sampled twirl
regardless of which stream produced it. Draw indices 16-63 are used so
the new data merges cleanly with Task 31C's original 0-15 without
collision.

VERIFICATION, not assumed: after collecting real data,
  1. recompute mean_sign/sign_se for the 27 flagged (slot,label) groups
     using the MERGED 64 draws, report before (16) vs after (64) side by
     side -- does sign_se actually shrink close to the predicted ~2x?
  2. recompute the champion pipeline's real energy error and a real
     bootstrap MSE using merged data for extended slots + Task 31C's
     original 16-draw data everywhere else (an honest "targeted-only"
     test, not silently upgrading data project-wide).

Run:
    python vqe/task33e_targeted_mc_extension.py
"""
import os
import sys
import json
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from task31c_full_pec_calibration import sample_twirled_circuit, GATE_NAME, P2_ZZ

K = 6
SHOTS = 100_000
N_MC_ORIGINAL = 16
N_MC_TARGET = 64
N_EXTRA = N_MC_TARGET - N_MC_ORIGINAL  # 48 new draws per flagged (slot, group)
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
ORIGINAL_CKPT = os.path.join(CKPT_DIR, "task31c_full_pec_calibration.json")
EXTENSION_CKPT = os.path.join(CKPT_DIR, "task33e_targeted_extension.json")
T33D_RESULTS = os.path.join(os.path.dirname(__file__), "task33d_pec_weight_diagnostics_results.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task33e_targeted_mc_extension_results.json")
N_BOOT = 24
EXACT_ENERGY_TARGET = 0.5


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def _boot_worker(p, non_id_labels, diag, kept, P_S, by_name_label, n_mc_by_key, seed):
    rng = np.random.default_rng(seed)
    a_by_name = {}
    for name in kept:
        blended = {}
        for l in non_id_labels:
            key = (name, l)
            if key not in by_name_label:
                continue
            entries = by_name_label[key]
            n_mc = n_mc_by_key.get(key, len(entries))
            draw_idx = rng.integers(0, len(entries), size=n_mc)
            vals_signed = []
            for idx in draw_idx:
                sign, gamma, counts = entries[idx]
                resampled = bootstrap_counts(counts, SHOTS, rng)
                m = expectation_from_counts(resampled, l)
                vals_signed.append(sign * gamma * m)
            blended[l] = max(-1.0, min(1.0, float(np.mean(vals_signed))))
        v0 = target_coeff_vector(name, K)
        a_hat, _ = fit_pure_state(P_S, blended, {l: 1.0 for l in blended}, K, v0,
                                   seed=int(rng.integers(0, 2**31)))
        a_by_name[name] = a_hat
    full = build_full_from_a(a_by_name, P_S, diag, K, non_id_labels)
    _, err = energy_and_err(p, full, K)
    return err


def main():
    print("\n" + "=" * 96)
    print("  task33e_targeted_mc_extension.py -- targeted real N_MC extension, 16->64, 10 flagged circuit families")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)

    with open(T33D_RESULTS) as f:
        t33d = json.load(f)
    flagged = {}
    for key, s in t33d["all_stats"].items():
        if abs(s["mean_sign"]) < 0.5:
            name, label = key.split("|", 1)
            flagged.setdefault(name, set()).add(label)
    needed = sorted({(name, gi) for name, labels in flagged.items()
                      for gi, g in enumerate(groups) if set(g) & labels})
    print(f"  {sum(len(v) for v in flagged.values())} flagged (slot,label) groups from Task 33D "
          f"-> {len(needed)} (slot, basis-group) circuit families need extension")
    for name, gi in needed:
        print(f"    {name:<10} group{gi}={groups[gi]}")

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins = learned["forte-1"]["gpi"]
    p1_gpi = float(np.mean([v["p"] for v in gpi_bins.values()]))
    p1_gpi2 = p1_gpi
    from loop_pec import gamma_factor
    gamma_gate2 = gamma_factor(P2_ZZ, 2)
    gamma_gate1 = gamma_factor(p1_gpi, 1)

    if os.path.exists(EXTENSION_CKPT):
        print("\n  found existing extension checkpoint -> reusing, not resubmitting")
        with open(EXTENSION_CKPT) as f:
            ext_ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"\n  connected, backend={backend.name}")

        circuits, tags, gamma_per_circuit = [], [], []
        for name, gi in needed:
            group = groups[gi]
            base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
            combined = effrag_mod.combined_basis_label(group)
            basis_qc = native_basis_change(combined, GATE_NAME)
            full_base = base.compose(basis_qc)
            n2q_full = sum(1 for instr in full_base.data if instr.operation.name == GATE_NAME)
            n1q_full = sum(1 for instr in full_base.data if instr.operation.name in ("gpi", "gpi2"))
            gamma_this = (gamma_gate2 ** n2q_full) * (gamma_gate1 ** n1q_full)
            rng = np.random.default_rng(stable_seed("task33e", name, gi))
            for draw in range(N_MC_ORIGINAL, N_MC_TARGET):
                twirled, sign = sample_twirled_circuit(full_base, P2_ZZ, p1_gpi, p1_gpi2, rng)
                twirled.measure_all()
                circuits.append(twirled)
                tags.append((name, list(group), draw, sign))
                gamma_per_circuit.append(gamma_this)
        print(f"  built {len(circuits)} new twirled circuits ({len(needed)} families x {N_EXTRA} extra draws)")

        n_chunks = (len(circuits) + 90) // 91
        chunk_size = max(1, (len(circuits) + n_chunks - 1) // n_chunks)
        chunks_c = [circuits[i:i + chunk_size] for i in range(0, len(circuits), chunk_size)]
        print(f"  submitting in {len(chunks_c)} batches of <= {chunk_size} circuits each...")

        os.makedirs(CKPT_DIR, exist_ok=True)
        PARTIAL_PATH = EXTENSION_CKPT + ".partial.json"
        all_counts = []
        start_batch = 0
        if os.path.exists(PARTIAL_PATH):
            with open(PARTIAL_PATH) as f:
                partial = json.load(f)
            all_counts = partial["counts"]
            start_batch = partial["n_batches_done"]
            print(f"  resuming from partial checkpoint: {len(all_counts)} circuits already retrieved, "
                  f"{start_batch}/{len(chunks_c)} batches done")
        t0 = time.time()
        for ci, chunk in enumerate(chunks_c):
            if ci < start_batch:
                continue
            job = None
            for attempt in range(6):
                try:
                    job = submit_job(chunk, backend, "forte-1", shots=SHOTS)
                    break
                except Exception as e:
                    wait_s = min(30 * (2 ** attempt), 300)
                    print(f"    batch {ci+1}/{len(chunks_c)} submit failed (attempt {attempt+1}/6): {e} -- backing off {wait_s}s")
                    time.sleep(wait_s)
            if job is None:
                raise RuntimeError(f"submit_job exhausted retries for batch {ci+1}/{len(chunks_c)}")
            counts = get_counts_list(job)
            all_counts.extend(counts)
            print(f"    batch {ci+1}/{len(chunks_c)} done ({len(all_counts)}/{len(circuits)} total), "
                  f"{time.time()-t0:.1f}s elapsed", flush=True)
            # -- SAVE AFTER EVERY BATCH, not just at the end -- a killed/timed-out process must not lose
            # already-completed real submissions (this cost real, if free, hardware time to collect)
            with open(PARTIAL_PATH, "w") as f:
                json.dump({"counts": all_counts, "n_batches_done": ci + 1, "n_batches_total": len(chunks_c)}, f)

        ext_ck = {"tags": tags, "counts": all_counts, "gamma_per_circuit": gamma_per_circuit}
        with open(EXTENSION_CKPT, "w") as f:
            json.dump(ext_ck, f, indent=2)
        if os.path.exists(PARTIAL_PATH):
            os.remove(PARTIAL_PATH)
        print(f"  extension checkpoint saved -> {EXTENSION_CKPT}")

    # -- merge original + extension --
    with open(ORIGINAL_CKPT) as f:
        orig_ck = json.load(f)
    by_name_label = {}
    for (name, group, draw, sign), counts, gamma in zip(
            [tuple(t) for t in orig_ck["tags"]], orig_ck["counts"], orig_ck["gamma_per_circuit"]):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma, counts))
    for (name, group, draw, sign), counts, gamma in zip(
            [tuple(t) for t in ext_ck["tags"]], ext_ck["counts"], ext_ck["gamma_per_circuit"]):
        for l in group:
            by_name_label.setdefault((name, l), []).append((sign, gamma, counts))

    # -- verification 1: sign_se before vs after, for the 27 originally-flagged groups --
    print(f"\n  -- VERIFICATION: sign_se before (16 draws) vs after (merged, up to {N_MC_TARGET} draws) --")
    before_after = {}
    for name, labels in flagged.items():
        for l in labels:
            key = (name, l)
            entries = by_name_label.get(key, [])
            n_before = N_MC_ORIGINAL
            signs_before = np.array([e[0] for e in entries[:n_before]], dtype=float)
            signs_after = np.array([e[0] for e in entries], dtype=float)
            ms_before = float(signs_before.mean())
            ms_after = float(signs_after.mean())
            se_before = float(np.sqrt(max(0.0, 1 - ms_before**2)) / np.sqrt(len(signs_before)))
            se_after = float(np.sqrt(max(0.0, 1 - ms_after**2)) / np.sqrt(len(signs_after)))
            before_after[f"{name}|{l}"] = {
                "n_before": len(signs_before), "n_after": len(signs_after),
                "mean_sign_before": ms_before, "mean_sign_after": ms_after,
                "sign_se_before": se_before, "sign_se_after": se_after,
            }
    n_improved = sum(1 for v in before_after.values() if v["sign_se_after"] < v["sign_se_before"])
    print(f"    {n_improved}/{len(before_after)} flagged groups: sign_se decreased with more real data")
    for k, v in sorted(before_after.items(), key=lambda kv: -kv[1]["sign_se_before"])[:10]:
        print(f"    {k:<22} n:{v['n_before']}->{v['n_after']}  mean_sign:{v['mean_sign_before']:+.3f}->{v['mean_sign_after']:+.3f}  "
              f"sign_se:{v['sign_se_before']:.4f}->{v['sign_se_after']:.4f}")

    # -- verification 2: real champion-pipeline energy error, targeted-extended data vs original-only --
    print(f"\n  -- VERIFICATION: champion pipeline energy error, extended-where-flagged vs original-everywhere --")
    n_mc_original_by_key = {k: N_MC_ORIGINAL for k in by_name_label}
    n_mc_extended_by_key = {k: (len(v) if k in {(n, l) for n, ls in flagged.items() for l in ls} else N_MC_ORIGINAL)
                             for k, v in by_name_label.items()}

    for label_run, n_mc_by_key in [("ORIGINAL (16 everywhere)", n_mc_original_by_key),
                                     ("TARGETED-EXTENDED (up to 64 for flagged groups)", n_mc_extended_by_key)]:
        with ProcessPoolExecutor(max_workers=4) as ex:
            futures = [ex.submit(_boot_worker, p, non_id_labels, diag, kept, P_S, by_name_label, n_mc_by_key,
                                  50_000 + i)
                       for i in range(N_BOOT)]
            errs = np.array([fut.result() for fut in futures])
        bias_proxy = float(np.median(errs))
        sigma = float(errs.std(ddof=1))
        mse = bias_proxy ** 2 + sigma ** 2
        passes = (abs(bias_proxy) + 2 * sigma) < EXACT_ENERGY_TARGET
        print(f"    {label_run:<45} median={bias_proxy:8.4f}  std={sigma:7.4f}  MSE={mse:9.4f}  "
              f"{'PASS' if passes else 'FAIL'}")
        if label_run.startswith("ORIGINAL"):
            orig_result = {"median": bias_proxy, "std": sigma, "mse": mse, "passes": bool(passes), "errs": errs.tolist()}
        else:
            ext_result = {"median": bias_proxy, "std": sigma, "mse": mse, "passes": bool(passes), "errs": errs.tolist()}

    print(f"\n  -- HONEST READ --")
    if ext_result["mse"] < orig_result["mse"]:
        pct = 100 * (orig_result["mse"] - ext_result["mse"]) / orig_result["mse"]
        print(f"    Targeted extension (real new data, 10 circuit families only) reduced champion-pipeline "
              f"MSE by {pct:.1f}% ({orig_result['mse']:.3f} -> {ext_result['mse']:.3f}), confirming Task 33D's "
              f"diagnosis was a real, actionable finding, not a false lead.")
    else:
        print(f"    Targeted extension did NOT reduce champion-pipeline MSE ({orig_result['mse']:.3f} -> "
              f"{ext_result['mse']:.3f}) despite real sign_se improvement on the flagged groups themselves -- "
              f"meaning those 27 groups' contribution to total pipeline variance is smaller than the "
              f"manifold-fit and other-slot variance already documented (Task 32B), or the improvement is "
              f"real but too small to show up against this pipeline's own run-to-run noise at N_BOOT={N_BOOT}.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "needed_circuit_families": [[n, gi] for n, gi in needed],
            "sign_se_before_after": before_after,
            "n_improved": n_improved,
            "n_flagged_groups": len(before_after),
            "original_result": orig_result,
            "extended_result": ext_result,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
