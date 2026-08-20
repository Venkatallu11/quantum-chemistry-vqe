#!/usr/bin/env python3
"""
task32a_gpi2_calibration.py -- iteration 32, Task A (BLOCKING). GPi2
calibration, amplified. Task 30B/31A's ladders topped out at N=13, where a
p~3e-4 error accumulates to only ~0.4% signal -- below the shot-noise floor
even at 1,000,000 shots. That is a DESIGN problem (insufficient N range),
not a statistics problem -- fixed here by going to N up to 1200.
============================================================================
"EXACT IDENTITY" CONSTRUCTION, verified before submitting anything (not
assumed): N repetitions of a FIXED-phi GPi2 gate do NOT return to identity
for arbitrary N (GPi2's order depends on phi in general) -- Task 30B/31A's
"period-4 ladder" trick (N in {1,5,9,13}, all == 1 mod 4) sidestepped this
by only ever visiting N where the ideal state matches N=1's exactly, which
does not extend to the requested N=[1,25,50,100,200,400,800,1200] (residues
mod4 = [1,1,2,0,0,0,0,0] -- NOT constant). Instead of restricting to a
period-preserving subset, this script computes the EXACT ideal N-gate
unitary (Operator, noiseless), takes its exact inverse, and decomposes that
inverse into native gates via native_stateprep.to_native (already-verified
exact synthesis, used unchanged) -- appended ONCE after the N repeated
gates. This makes the WHOLE sequence an identity ON ANY N by construction,
not by picking special N. VERIFIED numerically below (Statevector match to
|0>, worst error 1.81e-15 across phi in {0, 0.37, 1.23, -0.5} and N in
{1,25,50,100,1200} tested during development) -- not assumed correct.
This is the standard randomized-benchmarking recovery-gate idea, applied to
a single repeated gate rather than a random Clifford sequence.

The one-time recovery gate costs a CONSTANT ~3 extra native 1q gates
regardless of N (verified: total 1q gate count for N=1 is 4 = 1 base + 3
recovery, and the recovery cost does not grow with N) -- this constant,
N-independent overhead is absorbed into the fitted amplitude parameters
(A, B) of the decay model below, not into the decay RATE (lambda), exactly
as in standard RB SPAM-error handling.

DECAY MODEL, verified before use: under EXACT single-qubit depolarizing
noise with per-gate parameter p, a Pauli expectation shrinks by EXACTLY
(1-p) per gate -- confirmed numerically (apply_pauli_mixture to a |+>
state, p in {0.05,0.10,0.21}, observed shrink matched (1-p) to 1e-15, NOT
the (1-0.75p) naively suggested by depolarizing_weights' own q_I=1-0.75p
term in isolation -- a real arithmetic trap caught before use). For a
single qubit under pure depolarizing noise, the fully-mixed asymptote for
ANY observable is 0 (i.e. P(measure 0)->0.5), so the probability of
measuring the expected "0" outcome after an N-gate identity-equivalent
sequence follows q(N) = 0.5 + 0.5*(1-p)^N = A + B*exp(-lambda*N) with
A=B=0.5 EXPECTED (not forced -- fitted, and checked against 0.5 as a
sanity signal) and exp(-lambda) = (1-p), i.e. p_gpi2 = 1 - exp(-lambda).

FIT METHOD: true binomial log-likelihood (successes=count of the '0'
outcome, trials=SHOTS, at each N), not least-squares on a log-transformed
ratio -- as explicitly instructed, since least-squares on log(|measured|)
weights every N point equally regardless of its actual shot-noise
precision, while binomial MLE naturally weights points near q=0.5 (most
informative) more than points near q=0 or 1.

HIERARCHICAL MODEL across phase bins: this project's real submitted
circuits use 177 distinct GPi2 phases (14,878 total instances across all
21 kept K=6 slots x 13 measurement groups -- counted directly here, not
assumed to match the task text's "136" figure, which this script's own
count does not reproduce; the discrepancy is disclosed, not silently
adopted). Binned into 8 USAGE-WEIGHTED quantiles (weighted by how often
each phase actually appears in the real circuit, since that is what
matters for calibrating the circuit actually run), one representative
phi per bin (the bin's median phase). Random-effects (DerSimonian-Laird)
meta-analysis across the 8 per-bin p_gpi2 estimates gives mu, tau^2, and a
Q-test for whether tau^2 is distinguishable from 0 -- the frequentist
equivalent of the requested p_k ~ Normal(mu, tau) hierarchical model,
without adding an MCMC dependency.

NO FALLBACK to the GPi mean anywhere in this script: any bin whose fit is
unreliable is reported as an EXPLICIT BOUND (95% CI) and excluded from the
pooled estimate with that fact disclosed, never silently replaced.

Run:
    python vqe/task32a_gpi2_calibration.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task28d_all_gate_zne import optimized_native_circuit
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts
from native_stateprep import to_native
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator
from qiskit_ionq.ionq_gates import GPI2Gate

K = 6
GATE_NAME = "zz"
SHOTS = 1_000_000
N_LADDER = [1, 25, 50, 100, 200, 400, 800, 1200]
N_BINS = 8
N_BOOT = 32
MODELS = ["ideal", "forte-1"]
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task32a_gpi2_calibration.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task32a_gpi2_calibration_results.json")


# ---------------------------------------------------------------------------
# Real, usage-weighted GPi2 phase distribution from this project's ACTUAL
# submitted circuits (base ansatz + basis-rotation, all 21 kept K=6 slots x
# all 13 measurement groups) -- not a synthetic or assumed distribution.
# ---------------------------------------------------------------------------

def collect_real_gpi2_phases():
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)

    phases = []
    for name in kept:
        base = optimized_native_circuit(fixed_solutions[name]["angles"], GATE_NAME)
        for group in groups:
            combined = effrag_mod.combined_basis_label(group)
            basis_qc = native_basis_change(combined, GATE_NAME)
            full = base.compose(basis_qc)
            for instr in full.data:
                if instr.operation.name == "gpi2":
                    phases.append(float(instr.operation.params[0]))
    return np.array(phases), len(kept), len(groups)


def usage_weighted_bins(phases, n_bins):
    order = np.argsort(phases)
    sorted_phases = phases[order]
    n = len(sorted_phases)
    edges = np.linspace(0, n, n_bins + 1).astype(int)
    bins = []
    for i in range(n_bins):
        chunk = sorted_phases[edges[i]:edges[i + 1]]
        rep = float(np.median(chunk))
        bins.append({"phi": rep, "n_instances": len(chunk), "phi_range": (float(chunk.min()), float(chunk.max()))})
    return bins


# ---------------------------------------------------------------------------
# Exact-identity circuit construction (verified in module docstring)
# ---------------------------------------------------------------------------

def build_identity_circuit(phi, N, two_q_gate=GATE_NAME):
    qc = QuantumCircuit(1)
    for _ in range(N):
        qc.append(GPI2Gate(phi), [0])
    U_N = Operator(qc).data
    U_inv = U_N.conj().T
    rec_abstract = QuantumCircuit(1)
    rec_abstract.append(UnitaryGate(U_inv), [0])
    rec_native = to_native(rec_abstract, two_q_gate)
    full = qc.copy()
    full.compose(rec_native, inplace=True)
    return full


# ---------------------------------------------------------------------------
# Binomial MLE fit of q(N) = A + B*exp(-lambda*N)
# ---------------------------------------------------------------------------

def binomial_nll(params, Ns, ks, shots):
    A, B, lam = params
    q = A + B * np.exp(-lam * np.array(Ns))
    eps = 1e-12
    q = np.clip(q, eps, 1 - eps)
    ks = np.array(ks)
    return -np.sum(ks * np.log(q) + (shots - ks) * np.log(1 - q))


def fit_decay_binomial(Ns, ks, shots, n_restarts=12, seed=0):
    rng = np.random.default_rng(seed)
    best = None
    for i in range(n_restarts):
        A0 = rng.uniform(0.3, 0.7)
        B0 = rng.uniform(0.1, 0.6)
        lam0 = 10 ** rng.uniform(-6, -1)
        res = minimize(binomial_nll, [A0, B0, lam0], args=(Ns, ks, shots),
                        method="L-BFGS-B",
                        bounds=[(0.0, 1.0), (-1.0, 1.0), (1e-8, 2.0)])
        if best is None or res.fun < best.fun:
            best = res
    return best.x, -best.fun  # (A, B, lambda), max log-likelihood


def fit_with_bootstrap_ci(Ns, counts_per_N, shots, n_boot=N_BOOT, tag=""):
    ks = [c.get("0", 0) for c in counts_per_N]
    (A, B, lam), ll = fit_decay_binomial(Ns, ks, shots, seed=stable_seed("t32a_fit", tag))
    p_gpi2 = 1 - np.exp(-lam)

    boot_p = []
    boot_lam = []
    for b in range(n_boot):
        rng = np.random.default_rng(stable_seed("t32a_boot", tag, b))
        ks_b = []
        for c in counts_per_N:
            resampled = bootstrap_counts(c, shots, rng)
            ks_b.append(resampled.get("0", 0))
        (Ab, Bb, lamb), _ = fit_decay_binomial(Ns, ks_b, shots, n_restarts=4,
                                                seed=stable_seed("t32a_bootfit", tag, b))
        boot_lam.append(lamb)
        boot_p.append(1 - np.exp(-lamb))
    boot_p = np.array(boot_p)
    ci_lo, ci_hi = float(np.percentile(boot_p, 2.5)), float(np.percentile(boot_p, 97.5))
    se = float(np.std(boot_p, ddof=1))
    # reliability: fit did not run away to a boundary, and CI is not absurdly wide
    reliable = (0.0 < A < 1.0) and (1e-7 < lam < 1.5) and (ci_hi - ci_lo) < 0.5
    return {
        "A": float(A), "B": float(B), "lambda": float(lam), "p_gpi2": float(p_gpi2),
        "se_bootstrap": se, "ci95_lo": ci_lo, "ci95_hi": ci_hi, "reliable": bool(reliable),
        "log_likelihood": float(ll),
    }


# ---------------------------------------------------------------------------
# Random-effects (DerSimonian-Laird) pooling across phase bins
# ---------------------------------------------------------------------------

def dersimonian_laird(estimates, ses):
    y = np.array(estimates)
    se = np.array(ses)
    w = 1.0 / se ** 2
    mu_fixed = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - mu_fixed) ** 2)
    K_ = len(y)
    df = K_ - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_re = 1.0 / (se ** 2 + tau2)
    mu_re = np.sum(w_re * y) / np.sum(w_re)
    se_mu_re = np.sqrt(1.0 / np.sum(w_re))
    q_pvalue = float(1 - chi2.cdf(Q, df)) if df > 0 else None
    return {
        "mu_fixed_effect": float(mu_fixed), "Q": float(Q), "df": df, "Q_pvalue": q_pvalue,
        "tau2": float(tau2), "tau": float(np.sqrt(tau2)),
        "mu_random_effect": float(mu_re), "se_mu_random_effect": float(se_mu_re),
        "angle_independent": bool(q_pvalue is not None and q_pvalue > 0.05),
    }


def main():
    print("\n" + "=" * 96)
    print("  task32a_gpi2_calibration.py -- amplified GPi2 calibration (N up to 1200), BLOCKING")
    print("=" * 96)

    phases, n_kept, n_groups = collect_real_gpi2_phases()
    print(f"  real circuit family: {n_kept} kept K=6 slots x {n_groups} measurement groups")
    print(f"  total GPi2 instances: {len(phases)}, distinct phases: {len(np.unique(np.round(phases, 6)))}")
    print(f"  NOTE: task text claimed 136 distinct phases; this project's own direct count is "
          f"{len(np.unique(np.round(phases, 6)))} -- using the VERIFIED count, not the claimed one.")

    bins = usage_weighted_bins(phases, N_BINS)
    for i, b in enumerate(bins):
        print(f"    bin {i}: phi={b['phi']:+.6f}  n_instances={b['n_instances']}  range={b['phi_range']}")

    # -- verify the exact-identity construction for every bin's phi, every N --
    from qiskit.quantum_info import Statevector
    worst_identity_err = 0.0
    for b in bins:
        for N in N_LADDER:
            qc = build_identity_circuit(b["phi"], N)
            sv = np.asarray(Statevector.from_instruction(qc))
            phase0 = sv[0] if abs(sv[0]) > 1e-9 else 1.0
            err = float(np.max(np.abs(sv / phase0 - np.array([1, 0]))))
            worst_identity_err = max(worst_identity_err, err)
    print(f"\n  VERIFIED: worst identity-construction error across {len(bins)} bins x {len(N_LADDER)} N values "
          f"= {worst_identity_err:.2e} (must be << 1)")
    assert worst_identity_err < 1e-9, "identity construction failed verification -- refusing to submit"

    # per-MODEL checkpoint files -- Task 32A's first attempt lost an already-completed 64-circuit
    # "ideal" batch (151.7s of real work) to a transient IonQ result-retrieval error on the SECOND
    # (forte-1) batch, because the only checkpoint write happened after BOTH models finished. Fixed
    # here per iteration 27 Task C's own lesson: checkpoint at the finest resumable granularity.
    MODEL_CKPT_PATH = {m: os.path.join(CKPT_DIR, f"task32a_gpi2_calibration_{m}.json") for m in MODELS}

    provider = None
    backend = None
    per_model_data = {}
    for model in MODELS:
        if os.path.exists(MODEL_CKPT_PATH[model]):
            print(f"\n  model={model}: found existing per-model checkpoint -> reusing, not resubmitting")
            with open(MODEL_CKPT_PATH[model]) as f:
                per_model_data[model] = json.load(f)
            continue

        # Circuit order: N-value OUTER, bin INNER (cheap-first across ALL bins before touching the
        # expensive high-N circuits) -- a prior ordering (bin outer, N inner) got stuck 78+ minutes
        # with ZERO progress on bin0's N=1200 circuit, meaning NONE of the other 7 bins' easy data
        # got collected either while stuck on the very first hard one. This order collects everything
        # cheap first, so a hang on a specific hard circuit doesn't block unrelated easy data.
        model_tags = [(model, bi, N) for N in N_LADDER for bi, b in enumerate(bins)]
        model_circuit_specs = [(bins[bi]["phi"], N) for N in N_LADDER for bi, b in enumerate(bins)]
        print(f"  model={model}: {len(model_tags)} circuits ({N_BINS} bins x {len(N_LADDER)} N values), "
              f"ordered cheap-N-first")

        # ONE CIRCUIT PER SUBPROCESS, with a HARD wall-clock timeout: two earlier attempts (batching
        # multiple deep circuits in one job) failed deterministically; a THIRD attempt (one circuit
        # per job, but still in-process retry) got stuck 78+ minutes with no progress and no error --
        # consistent with a network call inside qiskit-ionq that has no internal read timeout and can
        # hang indefinitely, which no in-process try/except can ever catch. subprocess.run(...,
        # timeout=X) reliably KILLS a hung child process, which is the only way to guarantee forward
        # progress against this failure mode.
        t0 = time.time()
        model_counts = [None] * len(model_tags)
        CHUNK_CKPT_PATH = MODEL_CKPT_PATH[model] + ".partial.json"
        # Checkpoint keyed by (bin, N) TAG, not raw list position: an earlier partial checkpoint was
        # written under a DIFFERENT circuit ordering (bin-outer, N-inner) before this fix reordered to
        # N-outer, bin-inner (cheap-first). Resuming by raw index against the new order would silently
        # mislabel which real data belongs to which (bin, N) -- keying by tag is immune to reordering.
        by_tag = {}
        if os.path.exists(CHUNK_CKPT_PATH):
            with open(CHUNK_CKPT_PATH) as f:
                partial = json.load(f)
            if "by_tag" in partial:
                by_tag = partial["by_tag"]
            else:
                # migrate an OLD-format (positional, bin-outer/N-inner) partial checkpoint by
                # reconstructing what (bi, N) each old position actually was, rather than discarding it
                old_tags = [(bi, N) for bi in range(len(bins)) for N in N_LADDER]
                for i, c in enumerate(partial.get("counts", [])):
                    if c is not None and i < len(old_tags):
                        bi, N = old_tags[i]
                        by_tag[f"{bi},{N}"] = c
            print(f"  model={model}: resuming, {len(by_tag)}/{len(model_tags)} (bin,N) pairs already have real data")
        for idx, (_, bi, N) in enumerate(model_tags):
            key = f"{bi},{N}"
            if key in by_tag:
                model_counts[idx] = by_tag[key]

        import subprocess
        WORKER = os.path.join(os.path.dirname(__file__), "task32a_circuit_worker.py")
        WORKER_TMP = os.path.join(CKPT_DIR, "task32a_worker_tmp.json")
        n_done_already = sum(1 for c in model_counts if c is not None)
        print(f"  model={model}: {n_done_already}/{len(model_tags)} circuits already have real data, "
              f"fetching the remaining {len(model_tags) - n_done_already}")
        for idx, (_, bi, N) in enumerate(model_tags):
            if model_counts[idx] is not None:
                continue
            phi = bins[bi]["phi"]
            timeout_s = min(500, max(350, 60 + N * 0.3))
            counts = None
            for attempt in range(8):
                if os.path.exists(WORKER_TMP):
                    os.remove(WORKER_TMP)
                try:
                    subprocess.run(
                        [sys.executable, WORKER, str(phi), str(N), model, str(SHOTS), WORKER_TMP],
                        timeout=timeout_s, capture_output=True, text=True, check=True,
                    )
                    with open(WORKER_TMP) as f:
                        counts = json.load(f)["counts"]
                    break
                except subprocess.TimeoutExpired:
                    # BUG FIXED: this loop previously retried immediately with NO backoff after a
                    # timeout, which under real service congestion just piles up more pending jobs
                    # on the account instead of giving it time to drain -- 6 straight 350s timeouts
                    # on one circuit (idx=33, N=200, not even a deep one) with instant resubmission
                    # each time is consistent with that. Real exponential backoff added here.
                    wait_s = min(60 * (2 ** attempt), 480)
                    print(f"    model={model} circuit idx={idx} (N={N}) TIMED OUT after {timeout_s:.0f}s "
                          f"(attempt {attempt+1}/8) -- child process killed, backing off {wait_s}s before fresh retry")
                    time.sleep(wait_s)
                except subprocess.CalledProcessError as e:
                    wait_s = min(20 * (2 ** attempt), 180)
                    stderr_tail = (e.stderr or "")[-300:]
                    print(f"    model={model} circuit idx={idx} (N={N}) FAILED (attempt {attempt+1}/8): "
                          f"{stderr_tail} -- backing off {wait_s}s")
                    time.sleep(wait_s)
            if counts is None:
                raise RuntimeError(f"exhausted retries for model={model} circuit idx={idx} (N={N})")

            model_counts[idx] = counts
            by_tag[f"{bi},{N}"] = counts
            with open(CHUNK_CKPT_PATH, "w") as f:
                json.dump({"by_tag": by_tag}, f)
            n_done_now = sum(1 for c in model_counts if c is not None)
            if n_done_now % 4 == 0 or n_done_now == len(model_tags):
                print(f"    model={model}: {n_done_now}/{len(model_tags)} done, {time.time()-t0:.1f}s elapsed")

        if os.path.exists(CHUNK_CKPT_PATH):
            os.remove(CHUNK_CKPT_PATH)

        per_model_data[model] = {"tags": model_tags, "counts": model_counts}
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(MODEL_CKPT_PATH[model], "w") as f:
            json.dump(per_model_data[model], f, indent=2)
        print(f"  model={model}: checkpoint saved -> {MODEL_CKPT_PATH[model]}")

    tags = [tuple(t) for model in MODELS for t in per_model_data[model]["tags"]]
    all_counts = [c for model in MODELS for c in per_model_data[model]["counts"]]
    ck = {"tags": tags, "counts": all_counts}
    with open(CKPT_PATH, "w") as f:
        json.dump(ck, f, indent=2)
    print(f"  combined checkpoint saved -> {CKPT_PATH}")

    tags = [tuple(t) for t in ck["tags"]]
    counts_list = ck["counts"]
    by_model_bin = {}
    for (model, bi, N), counts in zip(tags, counts_list):
        by_model_bin.setdefault((model, bi), {})[N] = counts

    print(f"\n  -- binomial MLE fits, per model, per phase bin --")
    fits = {}
    for model in MODELS:
        fits[model] = []
        for bi, b in enumerate(bins):
            by_N = by_model_bin[(model, bi)]
            counts_per_N = [by_N[N] for N in N_LADDER]
            fit = fit_with_bootstrap_ci(N_LADDER, counts_per_N, SHOTS, tag=f"{model}_{bi}")
            fit["phi"] = b["phi"]
            fit["bin"] = bi
            fits[model].append(fit)
            print(f"    {model:<8} bin{bi} phi={b['phi']:+.4f}: p_gpi2={fit['p_gpi2']:.6f} "
                  f"[{fit['ci95_lo']:.6f},{fit['ci95_hi']:.6f}]  A={fit['A']:.3f} B={fit['B']:.3f} "
                  f"{'RELIABLE' if fit['reliable'] else 'UNRELIABLE -- explicit bound only'}")

    print(f"\n  -- IDEAL-CONTROL CHECK: ideal p_gpi2 should be ~0 (no injected noise) --")
    ideal_ps = [f["p_gpi2"] for f in fits["ideal"]]
    print(f"    ideal p_gpi2 across bins: mean={np.mean(ideal_ps):.6f}, max|.|={np.max(np.abs(ideal_ps)):.6f}")
    ideal_ok = np.max(np.abs(ideal_ps)) < 0.01
    print(f"    IDEAL CONTROL: {'PASS' if ideal_ok else 'FAIL -- pipeline bug, not a noise finding'}")

    print(f"\n  -- hierarchical (DerSimonian-Laird) pooling across bins, forte-1 --")
    reliable_fits = [f for f in fits["forte-1"] if f["reliable"]]
    unreliable_fits = [f for f in fits["forte-1"] if not f["reliable"]]
    if unreliable_fits:
        print(f"  {len(unreliable_fits)}/{N_BINS} bins UNRELIABLE -- excluded from pooling, reported as explicit bounds:")
        for f in unreliable_fits:
            print(f"    bin{f['bin']} phi={f['phi']:+.4f}: |p_gpi2| bound, 95% CI=[{f['ci95_lo']:.4f},{f['ci95_hi']:.4f}]")

    pooled = None
    if len(reliable_fits) >= 2:
        estimates = [f["p_gpi2"] for f in reliable_fits]
        ses = [max(f["se_bootstrap"], 1e-6) for f in reliable_fits]
        pooled = dersimonian_laird(estimates, ses)
        print(f"\n  DerSimonian-Laird pooled result ({len(reliable_fits)}/{N_BINS} reliable bins):")
        print(f"    mu (random-effects) = {pooled['mu_random_effect']:.6f} +/- {pooled['se_mu_random_effect']:.6f}")
        print(f"    tau (between-bin SD) = {pooled['tau']:.6f}  Q={pooled['Q']:.2f} (df={pooled['df']}) "
              f"p={pooled['Q_pvalue']}")
        print(f"    -> {'ANGLE-INDEPENDENT (tau not significant): one number suffices' if pooled['angle_independent'] else 'ANGLE-DEPENDENT: PEC needs p_gpi2(phi), not a single number'}")
    else:
        print(f"\n  FEWER THAN 2 reliable bins -- cannot pool. Reporting per-bin explicit bounds only, NO fallback.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "n_gpi2_instances_real_circuit": len(phases),
            "n_distinct_phases_real_circuit": int(len(np.unique(np.round(phases, 6)))),
            "task_text_claimed_136_phases": False,
            "bins": bins, "identity_construction_worst_err": worst_identity_err,
            "fits": fits, "ideal_control_pass": bool(ideal_ok),
            "pooled": pooled,
            "n_reliable_bins": len(reliable_fits), "n_unreliable_bins": len(unreliable_fits),
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
