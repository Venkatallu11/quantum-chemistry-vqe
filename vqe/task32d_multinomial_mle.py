#!/usr/bin/env python3
"""
task32d_multinomial_mle.py -- iteration 32, Task D. MULTINOMIAL LIKELIHOOD
ON RAW COUNTS. Every prior estimator in this project converts each real
circuit's raw bitstring counts to a per-label EXPECTATION VALUE first
(losing information: two circuits with very different count DISTRIBUTIONS
can give the same expectation value), then fits a physical state to those
756 derived numbers. This fits the physical state DIRECTLY against the raw
counts instead:

    p_{s,x}(a) = |<x| R_s |psi(a)>|^2,   maximize sum_{s,x} n_{s,x} log p_{s,x}(a)

subject to |a|=1 (parametrized a=v/|v|, same as the existing estimator).
R_s is the REAL basis-rotation unitary this project's own
`qforge.forging.basis_change_h` applies before measurement (reused
unchanged, exact match to how the real counts were generated -- not a
re-derived convention). |<P_l>| <= 1 holds by construction here for a much
more basic reason than the existing estimator's unit-sphere parametrization:
p_model is a genuine probability distribution from a normalized quantum
state, so no clipping is EVER needed, at any stage.
============================================================================
PEC AS A SEPARATE CORRECTION LAYER, not forced into the multinomial: PEC's
literal-twirling correction (Task 31C) produces SIGNED, individually
unphysical intermediate quantities (quasi-probability reweighting) that
cannot be multinomial counts themselves -- forcing them in would be
incoherent. Instead: fit `a` by multinomial MLE on RAW (untwirled) counts,
then apply Task 31C's own already-computed literal-PEC value as a
MULTIPLICATIVE per-label correction ratio to the MLE-predicted expectation
-- exactly the same "ratio correction" pattern task30b_pec_application.py
and task31f already use elsewhere in this project (corrected = predicted *
(PEC_value / raw_simple_value)), just applied to a multinomial-fit
prediction instead of a raw bootstrap mean.

Run:
    python vqe/task32d_multinomial_mle.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from qforge import combine_matrices, energy_from_alpha_matrices
from qforge.forging import basis_change_h, pauli_expectation
import ef_fragment as effrag_mod
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Pauli
from ionq_simulator_binding_curve import stable_seed

K = 6
RAW_CKPT = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints",
                         "task28b_optimized_raw.json")
PEC_RESULTS = os.path.join(os.path.dirname(__file__), "task31c_full_pec_calibration_results.json")
BASELINE_FORTE = 0.31659421376349306
BASELINE_IDEAL = 0.046613622763887454
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32d_multinomial_mle_results.json")


def probs_for_state(psi16, combined_label):
    qc = QuantumCircuit(4)
    basis_change_h(qc, combined_label)
    sv = Statevector(psi16).evolve(qc)
    return sv.probabilities_dict()


def multinomial_nll(v, U, groups_counts, K):
    a = v / np.linalg.norm(v)
    psi16 = U @ a
    nll = 0.0
    for combined_label, counts in groups_counts:
        probs = probs_for_state(psi16, combined_label)
        for bitstring, n in counts.items():
            p = probs.get(bitstring, 1e-15)
            nll -= n * np.log(max(p, 1e-15))
    return nll


def fit_mle(U, groups_counts, K, v0, seed, n_restarts=12):
    """L-BFGS-B, not Nelder-Mead: a direct reproduction test (main() reported
    36.43 kcal/mol on IDEAL data; re-running the identical function call in a
    fresh process gave 0.22) found Nelder-Mead lands on genuinely different
    local optima across runs for this 6D nonconvex problem -- the SAME
    cross-process degenerate-optimum failure mode Task 32B found and
    quantified for the existing L-BFGS-B-based estimator, now confirmed to
    also afflict Nelder-Mead, and evidently worse for it (no gradient
    information, weaker convergence guarantees). The multinomial NLL is
    smooth away from p=0 (regularized via max(p,1e-15)), so L-BFGS-B is a
    valid, and typically more reliable, choice here -- but per Task 32B's
    own lesson, robustness must be VERIFIED across real separate processes,
    not assumed from switching optimizers alone (see run_mle_for_model's
    median-of-independent-process-runs wrapper)."""
    rng = np.random.default_rng(seed)
    inits = [v0] + [v0 + rng.normal(0, s, K) for s in [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5]]
    best_val, best_v = float("inf"), None
    for v_init in inits[:n_restarts]:
        res = minimize(multinomial_nll, v_init, args=(U, groups_counts, K), method="L-BFGS-B")
        if res.fun < best_val:
            best_val, best_v = res.fun, res.x
    a_hat = best_v / np.linalg.norm(best_v)
    return a_hat, best_val


def target_coeff_vector(name, K):
    a = np.zeros(K)
    if name.startswith("u_"):
        a[int(name[2:])] = 1.0
        return a
    inner = name.strip("()")
    sign = "+" if "+" in inner else "-"
    n_str, m_str = inner.split(sign)
    n, m = int(n_str[1:]), int(m_str[1:])
    a[n] = 1.0 / np.sqrt(2)
    a[m] = (1.0 if sign == "+" else -1.0) / np.sqrt(2)
    return a


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return E, errs["err_vs_exact_kcal"]


def build_full_from_a(a_by_name, U, diag, K, non_id_labels):
    full = {}
    for name in diag:
        psi16 = U @ a_by_name[name]
        sv = Statevector(psi16)
        full[name] = {l: float(np.real(sv.expectation_value(Pauli(l)))) for l in non_id_labels}
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
            psi16 = U @ a_by_name[pl]
            sv = Statevector(psi16)
            full[pl] = {l: float(np.real(sv.expectation_value(Pauli(l)))) for l in non_id_labels}
            synth_minus = {}
            for l in non_id_labels:
                synth_minus[l] = full[un][l] + full[um][l] - full[pl][l]
            full[f"(u{n}-u{m})"] = synth_minus
    return full


def mle_worker(U, groups_counts, K, v0, seed):
    a_hat, nll = fit_mle(U, groups_counts, K, v0, seed)
    return a_hat.tolist(), nll


def run_mle_for_model(p, non_id_labels, diag, kept, U, model_tags, model_counts, seed_tag, n_proc_check=0):
    per_name_groups = {}
    for (name, group), counts in zip(model_tags, model_counts):
        combined = effrag_mod.combined_basis_label(group)
        per_name_groups.setdefault(name, []).append((combined, counts))

    a_by_name = {}
    nlls = {}
    cross_process_check = None
    for name in kept:
        groups_counts = per_name_groups[name]
        v0 = target_coeff_vector(name, K)
        a_hat, nll = fit_mle(U, groups_counts, K, v0, seed=stable_seed(seed_tag, name))
        a_by_name[name] = a_hat
        nlls[name] = nll
        if name == kept[0] and n_proc_check > 0:
            with ProcessPoolExecutor(max_workers=n_proc_check) as ex:
                futures = [ex.submit(mle_worker, U, groups_counts, K, v0, stable_seed(seed_tag, name))
                           for _ in range(n_proc_check)]
                reruns = [fut.result() for fut in futures]
            worst = max(abs(nll - r[1]) for r in reruns)
            cross_process_check = {"worst_nll_diff": float(worst), "reproducible": bool(worst < 1e-6)}
    return a_by_name, nlls, cross_process_check


def _full_pipeline_worker(p, non_id_labels, diag, kept, U, model_tags, model_counts, seed_tag,
                           pec_kept_forte=None, raw_simple=None):
    """Picklable top-level worker: the FULL per-model pipeline (fit every
    slot, reconstruct, compute energy, and -- for forte-1 -- apply the
    PEC-ratio correction layer) in one call, for use across genuinely
    separate OS processes -- a direct reproduction test found the ideal
    model's reported 36.43 kcal/mol was a bad local optimum on ONE process
    launch (a fresh rerun of the identical call gave 0.22), so ANY single
    run of this pipeline must be treated as one draw from a real
    run-to-run distribution, not a definitive answer -- report median and
    range across independent processes instead."""
    a_by_name, nlls, _ = run_mle_for_model(p, non_id_labels, diag, kept, U, model_tags, model_counts,
                                            seed_tag, n_proc_check=0)
    full_mle = build_full_from_a(a_by_name, U, diag, K, non_id_labels)
    _, err_mle = energy_and_err(p, full_mle, K)

    err_mle_pec = None
    if pec_kept_forte is not None and raw_simple is not None:
        full_mle_pec = {}
        for name in kept:
            full_mle_pec[name] = {}
            for l in non_id_labels:
                if l not in full_mle[name]:
                    continue
                raw_s = raw_simple[name].get(l, None)
                pec_v = pec_kept_forte.get(name, {}).get(l, None)
                if raw_s is None or pec_v is None or abs(raw_s) < 1e-6:
                    full_mle_pec[name][l] = full_mle[name][l]
                    continue
                ratio = pec_v / raw_s
                full_mle_pec[name][l] = max(-1.0, min(1.0, full_mle[name][l] * ratio))
        for n in range(K):
            for m in range(K):
                if n >= m:
                    continue
                un, um, pl = f"u_{n}", f"u_{m}", f"(u{n}+u{m})"
                for l in non_id_labels:
                    if l not in full_mle_pec[un] or l not in full_mle_pec[um] or l not in full_mle_pec[pl]:
                        continue
                    full_mle_pec.setdefault(f"(u{n}-u{m})", {})[l] = (
                        full_mle_pec[un][l] + full_mle_pec[um][l] - full_mle_pec[pl][l])
        _, err_mle_pec = energy_and_err(p, full_mle_pec, K)
    return err_mle, err_mle_pec


def main():
    print("\n" + "=" * 96)
    print("  task32d_multinomial_mle.py -- fit the physical state to RAW COUNTS directly (no expectation-value step)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    diag, plus, kept = kept_slots_for_K(K)
    U = np.asarray(p["u_vecs"])
    if U.shape[0] != 16:
        U = U.T
    print(f"  U shape (should be 16 x {K}): {U.shape}")

    with open(RAW_CKPT) as f:
        raw_ck = json.load(f)
    with open(PEC_RESULTS) as f:
        pec_results = json.load(f)
    pec_kept_forte = pec_results["pec_kept_forte"]

    # -- PEC as a SEPARATE correction layer: precompute the (deterministic) raw-simple expectation
    # values once, outside the per-process fit, so every worker uses the identical correction ratio --
    tags_forte = [tuple(t) for t in raw_ck["tags"]["forte-1"]]
    counts_forte = raw_ck["counts"]["forte-1"]
    per_name_counts = {}
    for (name, group), counts in zip(tags_forte, counts_forte):
        per_name_counts.setdefault(name, {}).setdefault(tuple(group), counts)
    raw_simple = {}
    for name in kept:
        raw_simple[name] = {}
        for group_t, counts in per_name_counts[name].items():
            total = sum(counts.values())
            probs = {b: n / total for b, n in counts.items()}
            for l in group_t:
                raw_simple[name][l] = pauli_expectation(probs, l)

    N_PROC_TRIALS = 8
    results_by_model = {}
    for model in ["ideal", "forte-1"]:
        print(f"\n  -- {model}: multinomial MLE fit, {N_PROC_TRIALS} independent process launches "
              f"(per Task 32B's precedent -- a single run is not a trustworthy point estimate for this "
              f"class of nonconvex fit) --")
        tags = [tuple(t) for t in raw_ck["tags"][model]]
        counts_list = raw_ck["counts"][model]
        pec_arg = pec_kept_forte if model == "forte-1" else None
        raw_arg = raw_simple if model == "forte-1" else None

        with ProcessPoolExecutor(max_workers=4) as ex:
            futures = [ex.submit(_full_pipeline_worker, p, non_id_labels, diag, kept, U, tags, counts_list,
                                  f"t32d_{model}_trial{i}", pec_arg, raw_arg)
                       for i in range(N_PROC_TRIALS)]
            trials = [fut.result() for fut in futures]
        err_raw_trials = np.array([t[0] for t in trials])
        print(f"    raw MLE energy across {N_PROC_TRIALS} independent processes: "
              f"{[round(e, 3) for e in err_raw_trials]}")
        print(f"    median={np.median(err_raw_trials):.4f}  min={err_raw_trials.min():.4f}  "
              f"max={err_raw_trials.max():.4f}  ({len(np.unique(np.round(err_raw_trials,4)))}/{N_PROC_TRIALS} distinct)")
        results_by_model[model] = {
            "err_raw_trials": err_raw_trials.tolist(), "err_raw_median": float(np.median(err_raw_trials)),
            "err_raw_min": float(err_raw_trials.min()), "err_raw_max": float(err_raw_trials.max()),
        }
        if model == "forte-1":
            err_pec_trials = np.array([t[1] for t in trials])
            print(f"    +PEC-ratio-layer energy across trials: {[round(e, 3) for e in err_pec_trials]}")
            print(f"    median={np.median(err_pec_trials):.4f}  min={err_pec_trials.min():.4f}  "
                  f"max={err_pec_trials.max():.4f}")
            results_by_model[model]["err_pec_trials"] = err_pec_trials.tolist()
            results_by_model[model]["err_pec_median"] = float(np.median(err_pec_trials))
            results_by_model[model]["err_pec_min"] = float(err_pec_trials.min())
            results_by_model[model]["err_pec_max"] = float(err_pec_trials.max())

    print(f"\n" + "=" * 96)
    print(f"  SUMMARY (median across {N_PROC_TRIALS} independent process launches, NOT a single run)")
    ideal_med = results_by_model["ideal"]["err_raw_median"]
    print(f"    ideal (sanity control): multinomial-MLE median = {ideal_med:.4f} "
          f"[range {results_by_model['ideal']['err_raw_min']:.4f}-{results_by_model['ideal']['err_raw_max']:.4f}] "
          f"(Task 31D uniform-fit baseline: {BASELINE_IDEAL:.4f})")
    print(f"    forte-1 RAW: median = {results_by_model['forte-1']['err_raw_median']:.4f} kcal/mol "
          f"[range {results_by_model['forte-1']['err_raw_min']:.4f}-{results_by_model['forte-1']['err_raw_max']:.4f}]")
    print(f"    forte-1 + PEC-ratio layer: median = {results_by_model['forte-1']['err_pec_median']:.4f} kcal/mol "
          f"[range {results_by_model['forte-1']['err_pec_min']:.4f}-{results_by_model['forte-1']['err_pec_max']:.4f}]  "
          f"(Task 31D nonlinear-fit-on-PEC-data baseline: {BASELINE_FORTE:.4f})")
    ideal_ok = ideal_med <= BASELINE_IDEAL * 3 + 0.2
    print(f"    IDEAL CONTROL (on median): {'PASS' if ideal_ok else 'FAIL'}")
    print(f"    NOTE: even the MEDIAN across {N_PROC_TRIALS} runs is a materially uncertain summary of a "
          f"heavy-tailed, run-dependent quantity -- Task 32B's lesson applies here just as much as to the "
          f"existing estimator; see the range reported alongside every median above.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"results_by_model": results_by_model, "ideal_control_pass": bool(ideal_ok),
                   "n_proc_trials": N_PROC_TRIALS,
                   "baseline_forte_task31d": BASELINE_FORTE, "baseline_ideal_task31d": BASELINE_IDEAL}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
