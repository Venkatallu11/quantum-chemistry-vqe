#!/usr/bin/env python3
"""
task32f_alpha_determination.py -- iteration 32, Task F. Task 31F's
alpha=-0.154 (want -0.5) came from only 8 real submissions -- far too few
to call it a law. Extends to N=32 real independent submissions (reusing
Task 31F's own 8 checkpointed submissions unchanged, adding 24 new ones),
records real wall-clock TIMESTAMPS at submission and retrieval (not
previously saved anywhere in this project), and fits

    SE(N) = a * N^alpha + b

by nonlinear least squares (3 free parameters, not the old 2-parameter
log-log linear fit) -- the INTERCEPT b matters more than alpha: b~0 means
averaging eventually reaches chemical accuracy; b>0.25 means it structurally
cannot, no matter how many submissions are averaged.

Reports autocorrelation rho_k for lags 1-5 (not just lag-1) and effective
sample size N_eff = N / (1 + 2*sum(rho_k)) -- if submissions are not i.i.d.
(plausible: IonQ's free simulator may share underlying random state or
load-dependent behavior across temporally close submissions), raw N
overstates the real independence this averaging argument needs.

Circuits, calibration parameters, and per-submission analytic-PEC+manifold
pipeline are UNCHANGED from Task 31F (same production circuit, same
P2_ZZ=0.0146, same gpi_bins fallback) -- this task only extends N and adds
timing/analysis, it does not re-derive the estimator.

Run:
    python vqe/task32f_alpha_determination.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.optimize import curve_fit

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
N_SUBMISSIONS = 32
P2_ZZ = 0.0146
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH_TMPL = os.path.join(CKPT_DIR, "task31f_convergence_submission_{}.json")  # REUSE Task 31F's 0-7
TIMING_PATH = os.path.join(os.path.dirname(__file__), "task32f_submission_timestamps.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32f_alpha_determination_results.json")


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def main():
    print("\n" + "=" * 96)
    print(f"  task32f_alpha_determination.py -- {N_SUBMISSIONS} independent real submissions "
          f"(reusing Task 31F's 8, adding 24), with timestamps")
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

    timestamps = {}
    if os.path.exists(TIMING_PATH):
        with open(TIMING_PATH) as f:
            timestamps = json.load(f)

    provider = None
    E_list = []
    for sub in range(N_SUBMISSIONS):
        ckpt_path = CKPT_PATH_TMPL.format(sub)
        if os.path.exists(ckpt_path):
            print(f"  submission {sub+1}/{N_SUBMISSIONS}: found existing checkpoint, reusing")
            with open(ckpt_path) as f:
                ck = json.load(f)
            if str(sub) not in timestamps:
                # Task 31F's original 8 predate timestamp recording -- disclosed as missing, not fabricated
                timestamps[str(sub)] = {"submit_time": None, "retrieve_time": None,
                                         "note": "predates timestamp recording (Task 31F original)"}
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
            partial_path = ckpt_path + ".partial.json"
            start_chunk = 0
            if os.path.exists(partial_path):
                with open(partial_path) as f:
                    partial = json.load(f)
                all_counts = partial["counts"]
                start_chunk = len(all_counts) // chunk_size  # whole chunks already retrieved
                print(f"  submission {sub+1}/{N_SUBMISSIONS}: resuming from partial checkpoint, "
                      f"{len(all_counts)}/{len(circuits)} circuits already retrieved")
            submit_t0 = time.time()
            t0 = time.time()
            chunks = chunks[start_chunk:]
            for ci, chunk in enumerate(chunks):
                # retry as FRESH submit+retrieve pairs, not re-polling the same job's result():
                # Task 32A found (and this exact submission run just reproduced, at submission 21)
                # that re-polling the SAME job after a "invalid literal for int()" malformed-response
                # error reproduces identically -- a fresh job can package its response differently.
                counts = None
                for attempt in range(6):
                    job = None
                    for sub_attempt in range(4):
                        try:
                            job = submit_job(chunk, backend, "forte-1", shots=SHOTS)
                            break
                        except Exception as e:
                            wait_s = min(30 * (2 ** sub_attempt), 300)
                            print(f"    submission {sub+1} chunk {ci+1}/{len(chunks)} submit failed "
                                  f"(attempt {sub_attempt+1}/4): {e} -- backing off {wait_s}s")
                            time.sleep(wait_s)
                    if job is None:
                        continue
                    try:
                        counts = get_counts_list(job)
                        break
                    except Exception as e:
                        wait_s = min(20 * (2 ** attempt), 180)
                        print(f"    submission {sub+1} chunk {ci+1}/{len(chunks)} RESULT RETRIEVAL failed "
                              f"(fresh-job attempt {attempt+1}/6): {e} -- backing off {wait_s}s and resubmitting fresh")
                        time.sleep(wait_s)
                if counts is None:
                    raise RuntimeError(f"exhausted fresh-submit retries: submission {sub+1} chunk {ci+1}")
                all_counts.extend(counts)
                # per-chunk checkpoint: a submission-level checkpoint alone (written only after ALL
                # chunks succeed) would lose earlier chunks' real work if a LATER chunk fails --
                # exactly what happened when submission 21 failed on some chunk after 20 clean
                # submissions in a row.
                with open(partial_path, "w") as f:
                    json.dump({"tags": tags[:len(all_counts)], "counts": all_counts}, f)
            retrieve_t0 = time.time()
            print(f"  submission {sub+1}/{N_SUBMISSIONS}: {len(circuits)} circuits total "
                  f"({len(all_counts)} retrieved), {time.time()-t0:.1f}s this run")
            ck = {"tags": tags, "counts": all_counts}
            os.makedirs(CKPT_DIR, exist_ok=True)
            with open(ckpt_path, "w") as f:
                json.dump(ck, f, indent=2)
            if os.path.exists(partial_path):
                os.remove(partial_path)
            timestamps[str(sub)] = {"submit_time": submit_t0, "retrieve_time": retrieve_t0,
                                     "wall_clock_s": retrieve_t0 - submit_t0}
            with open(TIMING_PATH, "w") as f:
                json.dump(timestamps, f, indent=2)

        tags = ck["tags"]
        counts_list = ck["counts"]
        per_name = {}
        for (name, group), counts in zip(tags, counts_list):
            per_name.setdefault(name, {}).setdefault(tuple(group), counts)

        manifold_kept = {}
        for name in kept:
            m_dict = {}
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("task32f_conv", sub, name, seed))
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
                                       seed=stable_seed("task32f_conv_fit", sub, name))
            manifold_kept[name] = a_hat

        full = build_full_from_a(manifold_kept, P_S, diag, K, non_id_labels)
        E, err = energy_and_err(p, full, K)
        E_list.append(err)
        print(f"  submission {sub+1}/{N_SUBMISSIONS}: E_err = {err:.4f} kcal/mol")

    E_arr = np.array(E_list)
    print(f"\n  E_1..E_{N_SUBMISSIONS}: {[round(e, 4) for e in E_list]}")
    print(f"  mean={E_arr.mean():.4f}  std={E_arr.std(ddof=1):.4f}")

    # -- running SE, fit SE(N) = a*N^alpha + b (3-param nonlinear, not 2-param log-log) --
    running_se = []
    for n in range(2, N_SUBMISSIONS + 1):
        running_se.append(float(np.std(E_arr[:n], ddof=1) / np.sqrt(n)))
    Ns = np.arange(2, N_SUBMISSIONS + 1, dtype=float)
    running_se = np.array(running_se)

    def se_model(N, a, alpha, b):
        return a * N ** alpha + b

    try:
        popt, pcov = curve_fit(se_model, Ns, running_se, p0=[running_se[0], -0.5, 0.05],
                                bounds=([0, -3, 0], [np.inf, 0, 1.0]), maxfev=10000)
        a_fit, alpha_fit, b_fit = popt
        perr = np.sqrt(np.diag(pcov))
        print(f"\n  3-PARAMETER FIT: SE(N) = {a_fit:.4f}*N^{alpha_fit:.4f} + {b_fit:.4f}")
        print(f"    parameter SEs: a={perr[0]:.4f}, alpha={perr[1]:.4f}, b={perr[2]:.4f}")
        print(f"    INTERCEPT b={b_fit:.4f}: {'averaging CAN eventually reach chemical accuracy (b<0.25)' if b_fit < 0.25 else 'averaging CANNOT reach chemical accuracy no matter how many submissions (b>=0.25)'}")
    except Exception as e:
        print(f"\n  3-parameter fit failed to converge: {e}")
        a_fit = alpha_fit = b_fit = None
        perr = [None, None, None]

    # -- autocorrelation, lags 1-5, and N_eff --
    x = E_arr - E_arr.mean()
    denom = np.sum(x ** 2)
    rhos = {}
    for k in range(1, 6):
        if k >= N_SUBMISSIONS:
            break
        rhos[k] = float(np.sum(x[:-k] * x[k:]) / denom) if denom > 0 else 0.0
    sum_rho = sum(rhos.values())
    n_eff = N_SUBMISSIONS / (1 + 2 * sum_rho) if (1 + 2 * sum_rho) > 0 else N_SUBMISSIONS
    print(f"\n  autocorrelation rho_1..rho_5: {rhos}")
    print(f"  N_eff = N/(1+2*sum(rho_k)) = {N_SUBMISSIONS}/(1+2*{sum_rho:.4f}) = {n_eff:.2f} "
          f"(vs raw N={N_SUBMISSIONS})")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "E_list": E_list, "mean": float(E_arr.mean()), "std": float(E_arr.std(ddof=1)),
            "n_submissions": N_SUBMISSIONS,
            "fit_a": float(a_fit) if a_fit is not None else None,
            "fit_alpha": float(alpha_fit) if alpha_fit is not None else None,
            "fit_b": float(b_fit) if b_fit is not None else None,
            "fit_param_se": [float(x) if x is not None else None for x in perr],
            "autocorrelation": rhos, "n_eff": float(n_eff),
            "timestamps": timestamps,
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
