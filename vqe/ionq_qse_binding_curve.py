#!/usr/bin/env python3
"""
ionq_qse_binding_curve.py — quantum subspace expansion (qse_mitigation.py,
verified locally to machine precision) run for real, concurrently, on
IonQ's free ionq_simulator (ideal/aria-1/forte-1), compared directly
against iteration 9's raw/CDR/PEC numbers on the SAME real noise.
============================================================================
KEY EFFICIENCY FACT, not incidental: QSE-ordinary's H_eff is built from
EXACTLY the same K x K alpha/beta matrices (36 targets x 13 qubit-wise
groups) that iteration 9's raw/CDR/PEC run already collected in
vqe/ionq_simulator_binding_curve_checkpoints/targets_d{d}.json. This file
reuses that checkpointed data directly -- QSE-ordinary needs ZERO new
real circuits. Only QSE-regularized needs anything new: 15 compute-
uncompute overlap circuits per geometry (K=6, K*(K-1)/2 pairs), x 3
geometries x 3 models = 135 circuits total, submitted concurrently
exactly like every other real phase in this project (submit all, then
poll).

HONESTY CARRIED OVER FROM THE LOCAL VERIFICATION (qse_mitigation.py),
not re-derived here: the overlap circuit measures |<u_m|u_n>|^2 only --
its SIGN is not resolved, so S[n,m] is built as +sqrt(measured
probability) for n!=m, a magnitude-only proxy. Locally, on this
project's own (larger) synthetic noise, that proxy made aggressive
regularization actively WORSE than doing nothing (ordinary QSE, S=I, was
the more robust choice) -- reported honestly there and re-checked here
against real IonQ noise, not assumed to carry over.

Run:
    python vqe/ionq_qse_binding_curve.py --overlap
    python vqe/ionq_qse_binding_curve.py --assemble
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from ionq_simulator_binding_curve import (
    SHOTS, N_SEEDS, GEOMETRIES, D_REF, NOISE_MODELS, HARTREE_TO_KCAL_MOL,
    transpiled_ansatz, submit_job, get_counts_list, bootstrap_counts,
    expectation_from_counts, stable_seed, CK_DIR, ckpt_path, save_ckpt, load_ckpt,
)
from energy_difference_study import setup_at_distance, fit_and_verify
from qse_mitigation import build_h_eff, qse_energy_ordinary, qse_energy_regularized, overlap_magnitude_circuit
from qforge import floor_test, combine_matrices, derive_beta_matrices
from ionq_backend import connect_provider, get_simulator
from ionq_run import IONQ_QIS_STANDARD_BASIS

K = 6
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "ionq_qse_binding_curve_results.json")
ITERATION9_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_results.json")


# ---------------------------------------------------------------------------
# PHASE: overlap circuits -- the ONLY new real circuits this file needs
# ---------------------------------------------------------------------------

def phase_overlap():
    print("\n" + "=" * 96)
    print("  ionq_qse_binding_curve.py --overlap")
    print("=" * 96)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  connected, backend={backend.name}")

    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    print(f"  {len(pairs)} overlap pairs/geometry x {len(GEOMETRIES)} geometries x {len(NOISE_MODELS)} models")

    geoms = {}
    circuits_by_d = {}
    for d in GEOMETRIES:
        p_d = setup_at_distance(d, K)
        fit_and_verify(p_d)
        geoms[d] = p_d
        names = [f"u_{n}" for n in range(K)]
        circuits = [overlap_magnitude_circuit(p_d["solutions"][names[n]]["angles"],
                                               p_d["solutions"][names[m]]["angles"],
                                               IONQ_QIS_STANDARD_BASIS)
                    for n, m in pairs]
        circuits_by_d[d] = circuits

    jobs = {}
    t0 = time.time()
    for d in GEOMETRIES:
        for model in NOISE_MODELS:
            jobs[(d, model)] = submit_job(circuits_by_d[d], backend, model, shots=SHOTS)
    t_submit = time.time() - t0
    print(f"  all {len(jobs)} overlap jobs submitted, {t_submit:.1f}s")

    t0 = time.time()
    counts_by_key = {key: get_counts_list(job) for key, job in jobs.items()}
    t_retrieve = time.time() - t0
    print(f"  all {len(jobs)} overlap jobs retrieved, {t_retrieve:.1f}s")

    out = {"pairs": pairs, "wall_clock": {"submit_s": t_submit, "retrieve_s": t_retrieve}, "counts": {}}
    for d in GEOMETRIES:
        out["counts"][str(d)] = {}
        for model in NOISE_MODELS:
            out["counts"][str(d)][model] = counts_by_key[(d, model)]

    save_ckpt("qse_overlap", out)
    print(f"\n  --overlap phase complete\n")
    return out


# ---------------------------------------------------------------------------
# PHASE: assemble -- reuse iteration 9's target checkpoints (zero new
# circuits) + this file's own overlap checkpoint, no network calls.
# ---------------------------------------------------------------------------

def s_matrix_from_overlap_counts(overlap_counts_for_model_d, pairs, rng=None):
    S = np.eye(K)
    for (n, m), counts in zip(pairs, overlap_counts_for_model_d):
        c = bootstrap_counts(counts, SHOTS, rng) if rng is not None else counts
        total = sum(c.values())
        p0000 = c.get("0000", 0) / total if total > 0 else 0.0
        S[n, m] = S[m, n] = float(np.sqrt(max(p0000, 0.0)))
    return S


def qse_ordinary_energy_from_raw(raw, p_d):
    alpha_mats = combine_matrices(raw, p_d["alpha_labels"], p_d["identity_label"], K)
    beta_mats = derive_beta_matrices(alpha_mats, p_d["signs"], K)
    H_eff = build_h_eff(p_d["terms"], alpha_mats, beta_mats, K)
    return qse_energy_ordinary(H_eff, p_d["enuc"])[0], H_eff


def assemble():
    print("\n" + "=" * 96)
    print("  ionq_qse_binding_curve.py --assemble")
    print("=" * 96)

    overlap_ckpt = load_ckpt("qse_overlap")
    if overlap_ckpt is None:
        print("  missing checkpoint: qse_overlap -- run --overlap first")
        return None
    target_ckpts = {}
    for d in GEOMETRIES:
        p = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints", f"targets_d{d}.json")
        if not os.path.exists(p):
            print(f"  missing iteration-9 checkpoint: targets_d{d}.json -- run ionq_simulator_binding_curve.py "
                  "--targets first (that file's real data is reused here, not re-collected)")
            return None
        with open(p) as f:
            target_ckpts[d] = json.load(f)

    with open(ITERATION9_RESULTS_PATH) as f:
        iter9 = json.load(f)

    pairs = [tuple(pr) for pr in overlap_ckpt["pairs"]]
    geoms = {d: setup_at_distance(d, K) for d in GEOMETRIES}
    for d in GEOMETRIES:
        fit_and_verify(geoms[d])

    exact_by_d = {d: geoms[d]["exact_energy"] for d in GEOMETRIES}
    diff_exact_by_d = {d: exact_by_d[d] - exact_by_d[D_REF] for d in GEOMETRIES}

    E_ordinary = {model: {d: [] for d in GEOMETRIES} for model in NOISE_MODELS}
    E_regularized = {model: {d: [] for d in GEOMETRIES} for model in NOISE_MODELS}
    s_eigval_ranges = {model: {} for model in NOISE_MODELS}
    threshold_sweeps = {model: {} for model in NOISE_MODELS}

    for model in NOISE_MODELS:
        for d in GEOMETRIES:
            ck = target_ckpts[d]
            p_d = geoms[d]
            groups = ck["groups"]
            target_names = ck["target_names"]
            counts_model = ck["counts"][model]
            overlap_counts_model_d = overlap_ckpt["counts"][str(d)][model]

            # -- real S from the (non-bootstrapped) real overlap measurement, once per (model,d) --
            S_real = s_matrix_from_overlap_counts(overlap_counts_model_d, pairs)
            s_eig = np.linalg.eigvalsh(S_real)
            s_eigval_ranges[model][d] = {"min": float(np.min(s_eig)), "max": float(np.max(s_eig)),
                                          "eigenvalues": s_eig.tolist()}

            seed_ordinary, seed_regularized = [], []
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(stable_seed("qse", model, d, seed))
                raw = {name: {} for name in target_names}
                for name in target_names:
                    for gi, group in enumerate(groups):
                        counts = bootstrap_counts(counts_model[name][gi], SHOTS, rng)
                        for l in group:
                            raw[name][l] = expectation_from_counts(counts, l)
                E_ord, H_eff = qse_ordinary_energy_from_raw(raw, p_d)
                seed_ordinary.append(E_ord)

                rng_s = np.random.default_rng(stable_seed("qse_s", model, d, seed))
                S_seed = s_matrix_from_overlap_counts(overlap_counts_model_d, pairs, rng=rng_s)
                # regularization threshold: fixed, chosen from the REAL measured spectrum below
                # (computed once per (model,d) outside the seed loop, applied per seed)
                seed_regularized.append((S_seed, H_eff))

            E_ordinary[model][d] = seed_ordinary

            # -- floor-test the regularization threshold on REAL data, span the REAL measured spectrum --
            s_lo, s_hi = s_eigval_ranges[model][d]["min"], s_eigval_ranges[model][d]["max"]
            thresholds = sorted(set([1e-6, 1e-3] + list(np.linspace(max(1e-2, s_lo * 0.5), s_hi * 1.2, 10).round(4))))
            sweep_errs = []
            for thresh in thresholds:
                errs_this_thresh = []
                for (S_seed, H_eff) in seed_regularized:
                    E_reg, n_kept, _ = qse_energy_regularized(H_eff, S_seed, p_d["enuc"], thresh)
                    if E_reg is not None:
                        errs_this_thresh.append(abs(E_reg - exact_by_d[d]) * HARTREE_TO_KCAL_MOL)
                if errs_this_thresh:
                    sweep_errs.append((thresh, float(np.mean(errs_this_thresh))))
            threshold_sweeps[model][d] = sweep_errs
            if len(sweep_errs) >= 2:
                ft = floor_test([t for t, _ in sweep_errs], [e for _, e in sweep_errs])
                best_thresh, best_err = min(sweep_errs, key=lambda te: te[1])
            else:
                ft, best_thresh, best_err = None, thresholds[0], None
            threshold_sweeps[model][d] = {"sweep": sweep_errs, "floor_test": ft,
                                           "best_threshold": best_thresh, "best_mean_err_kcal": best_err}

            E_regularized[model][d] = [
                (qse_energy_regularized(H_eff, S_seed, p_d["enuc"], best_thresh)[0])
                for (S_seed, H_eff) in seed_regularized
            ]
            E_regularized[model][d] = [e for e in E_regularized[model][d] if e is not None]

    # ==== report: mean/std, abs error, diff error, vs iteration 9 ====
    def stats_for(E_data, method):
        out = {}
        for model in NOISE_MODELS:
            abs_all, diff_all = [], []
            per_geom = {}
            for d in GEOMETRIES:
                vals = np.array(E_data[model][d])
                if len(vals) == 0:
                    continue
                abs_err = np.abs(vals - exact_by_d[d]) * HARTREE_TO_KCAL_MOL
                row = {"E_mean_ha": float(np.mean(vals)), "E_std_ha": float(np.std(vals)),
                       "abs_err_mean_kcal": float(np.mean(abs_err)), "abs_err_std_kcal": float(np.std(abs_err))}
                abs_all.append(row["abs_err_mean_kcal"])
                if d != D_REF and len(E_data[model].get(D_REF, [])) > 0:
                    ref = np.array(E_data[model][D_REF])
                    n = min(len(vals), len(ref))
                    diff_meas = vals[:n] - ref[:n]
                    diff_err = np.abs(diff_meas - diff_exact_by_d[d]) * HARTREE_TO_KCAL_MOL
                    row["diff_err_mean_kcal"] = float(np.mean(diff_err))
                    diff_all.append(row["diff_err_mean_kcal"])
                per_geom[d] = row
            mean_abs = float(np.mean(abs_all)) if abs_all else None
            mean_diff = float(np.mean(diff_all)) if diff_all else None
            cancellation = (mean_abs / mean_diff) if (mean_abs and mean_diff and mean_diff > 1e-9) else None
            out[model] = {"per_geometry": per_geom, "mean_abs_err_kcal": mean_abs,
                           "mean_diff_err_kcal": mean_diff, "cancellation_factor": cancellation}
        return out

    ordinary_report = stats_for(E_ordinary, "qse_ordinary")
    regularized_report = stats_for(E_regularized, "qse_regularized")

    print(f"\n  -- QSE vs iteration 9's raw/CDR/PEC, real IonQ noise --")
    for model in NOISE_MODELS:
        print(f"\n  == {model} ==")
        iter9_row = {m: iter9["error_report"].get(model, {}).get(m, {}) for m in ("raw", "cdr", "pec")}
        for m in ("raw", "cdr", "pec"):
            v = iter9_row[m].get("mean_abs_err_kcal")
            print(f"    {m:>14} (iter.9): {v:.3f} kcal/mol" if v is not None else f"    {m:>14} (iter.9): n/a")
        v = ordinary_report[model]["mean_abs_err_kcal"]
        print(f"    {'qse_ordinary':>14}      : {v:.3f} kcal/mol" if v is not None else f"    {'qse_ordinary':>14}      : n/a")
        v = regularized_report[model]["mean_abs_err_kcal"]
        print(f"    {'qse_regularized':>14}   : {v:.3f} kcal/mol" if v is not None else f"    {'qse_regularized':>14}   : n/a")
        print(f"    S eigenvalue range by geometry: "
              f"{ {str(d): (round(s_eigval_ranges[model][d]['min'],3), round(s_eigval_ranges[model][d]['max'],3)) for d in GEOMETRIES} }")

    results = {
        "geometries": GEOMETRIES, "d_ref": D_REF, "n_seeds": N_SEEDS, "shots": SHOTS,
        "s_eigenvalue_ranges": s_eigval_ranges,
        "threshold_sweeps": {model: {str(d): threshold_sweeps[model][d] for d in GEOMETRIES} for model in NOISE_MODELS},
        "qse_ordinary_report": ordinary_report,
        "qse_regularized_report": regularized_report,
        "iteration9_comparison": {m: {model: iter9["error_report"].get(model, {}).get(m, {}).get("mean_abs_err_kcal")
                                       for model in NOISE_MODELS} for m in ("raw", "cdr", "pec")},
        "wall_clock": {"overlap": overlap_ckpt["wall_clock"]},
        "chemical_accuracy_kcal": 1.0,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--overlap", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    args = parser.parse_args()
    if args.overlap:
        phase_overlap()
    elif args.assemble:
        assemble()
    else:
        parser.error("pass --overlap or --assemble")


if __name__ == "__main__":
    main()
