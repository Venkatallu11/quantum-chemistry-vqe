#!/usr/bin/env python3
"""
task31f_convergence_study.py -- iteration 31, Task F, part 2 (scoped).
128 independent full-pipeline real submissions is not achievable in this
session's timeframe (Task C's single N_MC=16 run took ~4.1 hours; 128
would take ~22 days). Scoped to N_SUBMISSIONS=8 independent REAL
submissions of the cheaper analytic-PEC+manifold pipeline (raw 273-
circuit measurement per submission; calibration parameters already
established in Task A/B, reused unchanged, not re-measured each time).
============================================================================
Records E_1..E_8 (forte-1, analytic-PEC+manifold energy error vs exact),
fits log(SE_N) = c + alpha*log(N) across the running-mean standard error
at each cumulative N (want alpha ~ -0.5 for genuine sqrt(N) convergence),
and reports the lag-1 autocorrelation of the raw E_i sequence -- if the
"85 submissions" reasoning behind chemical-accuracy-by-averaging depends
on i.i.d. draws, real autocorrelation would invalidate that reasoning and
requires N_eff = N/(1+2*sum(rho_k)) instead of raw N.

Run:
    python vqe/task31f_convergence_study.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices, HARTREE_TO_KCAL_MOL
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state, build_full_from_a
from phys_constrained_reconstruction import build_P_S
from task30b_pec_application import analytic_A_and_B
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod

K = 6
SHOTS = 100_000
N_SEEDS = 8
GATE_NAME = "zz"
N_SUBMISSIONS = 8  # scoped down from 128 -- see module docstring
P2_ZZ = 0.0146
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH_TMPL = os.path.join(CKPT_DIR, "task31f_convergence_submission_{}.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31f_convergence_study_results.json")


def build_full(raw_kept, diag, K, non_id_labels):
    full = {name: dict(raw_kept[name]) for name in diag}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            full[pl] = dict(raw_kept[pl])
            synth_minus = {}
            for l in non_id_labels:
                if l not in raw_kept[pl] or l not in full[un] or l not in full[um]:
                    continue
                synth_minus[l] = full[un][l] + full[um][l] - raw_kept[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print(f"  task31f_convergence_study.py -- {N_SUBMISSIONS} independent real submissions (scoped from 128)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    U = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U)

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins = learned["forte-1"]["gpi"]

    provider = None
    E_list = []
    for sub in range(N_SUBMISSIONS):
        ckpt_path = CKPT_PATH_TMPL.format(sub)
        if os.path.exists(ckpt_path):
            print(f"  submission {sub+1}/{N_SUBMISSIONS}: found existing checkpoint, reusing")
            with open(ckpt_path) as f:
                ck = json.load(f)
        else:
            if provider is None:
                provider = connect_provider()
                backend = get_native_simulator(provider)
                print(f"  connected, backend={backend.name}")
            circuits, tags = [], []
            for name in kept:
                base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
                for group in groups:
                    combined = effrag_mod.combined_basis_label(group)
                    basis_qc = native_basis_change(combined, GATE_NAME)
                    qc = base.compose(basis_qc)
                    qc.measure_all()
                    circuits.append(qc)
                    tags.append((name, list(group)))
            n_chunks = (len(circuits) + 90) // 91
            chunk_size = (len(circuits) + n_chunks - 1) // n_chunks
            chunks = [circuits[i:i + chunk_size] for i in range(0, len(circuits), chunk_size)]
            all_counts = []
            t0 = time.time()
            for ci, chunk in enumerate(chunks):
                job = None
                for attempt in range(6):
                    try:
                        job = submit_job(chunk, backend, "forte-1", shots=SHOTS)
                        break
                    except Exception as e:
                        wait_s = min(30 * (2 ** attempt), 300)
                        print(f"    submission {sub+1} chunk {ci+1}/{len(chunks)} failed (attempt {attempt+1}/6): "
                              f"{e} -- backing off {wait_s}s")
                        time.sleep(wait_s)
                if job is None:
                    raise RuntimeError(f"exhausted retries: submission {sub+1} chunk {ci+1}")
                all_counts.extend(get_counts_list(job))
            print(f"  submission {sub+1}/{N_SUBMISSIONS}: {len(circuits)} circuits, {len(chunks)} batches, "
                  f"{time.time()-t0:.1f}s")
            ck = {"tags": tags, "counts": all_counts}
            os.makedirs(CKPT_DIR, exist_ok=True)
            with open(ckpt_path, "w") as f:
                json.dump(ck, f, indent=2)

        # -- analytic PEC + manifold on this submission's own real data --
        tags = ck["tags"]
        counts_list = ck["counts"]
        per_name = {}
        for (name, group), counts in zip(tags, counts_list):
            per_name.setdefault(name, {}).setdefault(tuple(group), counts)

        manifold_kept = {}
        for name in kept:
            m_dict = {}
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task31f_conv", sub, name, seed))
                for group_t, counts in per_name[name].items():
                    resampled = bootstrap_counts(counts, SHOTS, rng)
                    for l in group_t:
                        m_dict.setdefault(l, []).append(expectation_from_counts(resampled, l))
            m_dict = {l: float(np.mean(v)) for l, v in m_dict.items()}
            labels_here = list(m_dict.keys())
            A, B = analytic_A_and_B(fixed_solutions[name]["angles"], GATE_NAME, P2_ZZ, gpi_bins, gpi_bins, labels_here)
            pec_dict = {}
            for l in labels_here:
                ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
                pec_dict[l] = max(-1.0, min(1.0, m_dict[l] * ratio))
            v0 = target_coeff_vector(name, K)
            a_hat, _ = fit_pure_state(P_S, pec_dict, {l: 1.0 for l in pec_dict}, K, v0,
                                       seed=stable_seed("task31f_conv_fit", sub, name))
            manifold_kept[name] = a_hat

        full = build_full_from_a(manifold_kept, P_S, diag, K, non_id_labels)
        E, err = energy_and_err(p, full, K)
        E_list.append(err)
        print(f"  submission {sub+1}/{N_SUBMISSIONS}: E_err = {err:.4f} kcal/mol")

    E_arr = np.array(E_list)
    print(f"\n  E_1..E_{N_SUBMISSIONS}: {[round(e, 4) for e in E_list]}")
    print(f"  mean={E_arr.mean():.4f}  std={E_arr.std():.4f}")

    # -- running mean/SE, fit log(SE_N) = c + alpha*log(N) --
    running_se = []
    for n in range(2, N_SUBMISSIONS + 1):
        running_se.append(float(np.std(E_arr[:n], ddof=1) / np.sqrt(n)))
    Ns = np.arange(2, N_SUBMISSIONS + 1)
    if len(running_se) >= 2:
        log_se = np.log(np.array(running_se))
        log_n = np.log(Ns)
        alpha, c = np.polyfit(log_n, log_se, 1)
        print(f"\n  fit: log(SE_N) = {c:.4f} + {alpha:.4f}*log(N)  (want alpha ~ -0.5)")
    else:
        alpha = None

    # -- lag-1 autocorrelation --
    if N_SUBMISSIONS > 2:
        x = E_arr - E_arr.mean()
        rho1 = float(np.sum(x[:-1] * x[1:]) / np.sum(x ** 2)) if np.sum(x ** 2) > 0 else 0.0
    else:
        rho1 = None
    print(f"  lag-1 autocorrelation rho_1 = {rho1}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"E_list": E_list, "mean": float(E_arr.mean()), "std": float(E_arr.std()),
                   "alpha_fit": float(alpha) if alpha is not None else None, "rho1": rho1,
                   "n_submissions": N_SUBMISSIONS, "scoped_down_from": 128}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
