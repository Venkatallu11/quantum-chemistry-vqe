#!/usr/bin/env python3
"""
task35bc_stratified_control_variate_pec.py -- iteration 35, Tasks B+C.
Attack PEC Monte Carlo variance AT ITS SOURCE (confirmed dominant, ~87%
of total pipeline variance, Task 32B), via stratified sampling and a
control variate -- but adapted to how THIS codebase's PEC actually works,
not a generic QPD abstraction.

STRUCTURAL FACT this design depends on (verified directly in Task 33D by
reading real checkpoint entries, not assumed): `gamma` in
`sample_twirled_circuit` is FIXED per circuit -- computed from the
calibrated error rate alone (`gamma_gate = sum(|weights|)`, independent
of which label `idx` gets chosen). ONLY the total SIGN (product of each
gate's independently-sampled recovery-label sign) is stochastic per draw.
This means:
  - "QPD control variate with known mean" (CV4Quantum's core idea)
    concretely becomes: V = sign - E[sign], where E[sign] IS analytically
    knowable (product of per-gate marginal sign probabilities) -- not an
    abstract weight, THE actual stochastic quantity this estimator draws.
  - "Stratified product-QPD sampling" concretely becomes: stratify by
    total sign (a 2-outcome stratification -- the coarsest, cheapest,
    most defensible version, not a combinatorial full-joint stratification
    over all ~100+ per-gate choices, which would need astronomically many
    strata to be well-populated).

METHOD, avoiding the two failure modes this project has hit before
(reusing training data for hyperparameter choice, per Task 32a's outlier
diagnosis discipline; peeking at the exact answer to tune anything, per
this session's explicit instruction):
  1. PILOT phase (cheap, sign-only, NO density-matrix simulation needed
     -- sign is determined by the classical twirl choice alone): draw
     N_PILOT=4000 twirl realizations, record only the sign. This gives
     p_plus = P(total_sign=+1) essentially exactly (SE ~1.6% at N=4000).
     A SEPARATE pilot batch estimates per-stratum std(m) (the expensive
     part, small N, only for setting Neyman-vs-proportional allocation --
     not reused in the production estimate).
  2. PRODUCTION phase, REPEATED across many independent trials (so we
     measure the estimator's REAL variance empirically, not from one
     lucky/unlucky sample): at a fixed real-world-realistic budget N=16
     (matching this project's actual production N_MC) generate:
       (a) IID: N unrestricted twirl draws (status quo).
       (b) STRATIFIED: draws split across the sign=+1/sign=-1 strata
           by PROPORTIONAL allocation (n_s = round(N*p_s), the paper's
           own "never worse than naive under ideal allocation" claim
           requires only p_s, not per-stratum variance -- so proportional
           allocation is used as the primary, defensible test; Neyman
           allocation reported as a secondary check using pilot-estimated
           sigma_s, disclosed as using pilot data for that purpose only).
       (c) CONTROL VARIATE: IID draws, corrected using beta*(sign-E[sign])
           with beta estimated from the SAME production sample (standard
           CV practice -- beta is a nuisance parameter, not a target-
           dependent hyperparameter).
       (d) STRATIFIED + CV combined.
     Compare against the TRUE exact (zero-noise) value for validation
     ONLY -- never used to pick strata counts, beta, or anything else.

Run:
    python vqe/task35bc_stratified_control_variate_pec.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from native_stateprep import to_native
from ionq_simulator_binding_curve import stable_seed
from qiskit.quantum_info import DensityMatrix, Operator, Pauli, Statevector

K = 6
GATE_NAME = "zz"
P2 = 0.0146
P1 = 0.000119
N_PILOT_SIGN = 4000     # cheap: sign only, no DM sim
N_PILOT_VARIANCE = 48   # expensive: per-stratum std(m), for Neyman allocation only
N_PRODUCTION = 16       # matches this project's real N_MC
N_TRIALS = 60           # independent repeats of the size-16 production experiment, to measure real variance
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task35bc_stratified_control_variate_pec_results.json")

TEST_CASES = [
    {"slot": "(u0+u1)", "group": ["XYYX", "IYYI"], "label": "IYYI", "tag": "IYYI (known atypical, coherent-noise case)"},
    {"slot": "u_0", "group": ["IYIY", "YIYI"], "label": "YIYI", "tag": "u_0/YIYI (well-conditioned per Task 33D, mean_sign near +1.0)"},
    {"slot": "(u3+u5)", "group": ["XZXZ", "XZXI", "IZIZ", "IIIZ"], "label": "IIIZ",
     "tag": "(u3+u5)/IIIZ (Task 33D-flagged high sign-noise group, mean_sign=+0.25)"},
    {"slot": "(u0+u3)", "group": ["YZYZ", "YZYI"], "label": "YZYI",
     "tag": "(u0+u3)/YZYI (Task 33D-flagged high sign-noise group, mean_sign=+0.375)"},
]


def noisy_dm(qc, p2, p1):
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
        if op.name == GATE_NAME:
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p2, 2))
        elif op.name in ("gpi", "gpi2"):
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p1, 1))
    return dm


def sample_twirl_choice(base_qc, p2, p1, rng):
    """Returns ONLY the twirled circuit + sign (cheap, classical). gamma
    is fixed per circuit (verified Task 33D) so it's computed once
    outside this function, not per draw."""
    qc = base_qc.copy_empty_like()
    total_sign = 1
    for instr in base_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name in ("gpi", "gpi2"):
            weights = pec_inverse_weights(p1, 1)
        else:
            continue
        labels = list(weights.keys())
        w = np.array([weights[l] for l in labels])
        gamma_gate = float(np.sum(np.abs(w)))
        probs = np.abs(w) / gamma_gate
        idx = rng.choice(len(labels), p=probs)
        chosen_label, chosen_w = labels[idx], w[idx]
        sign = 1 if chosen_w >= 0 else -1
        total_sign *= sign
        n = len(chosen_label)
        if chosen_label != "I" * n:
            qc.append(Pauli(chosen_label).to_instruction(), qargs)
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign


def fixed_gamma(base_qc, p2, p1):
    gamma = 1.0
    for instr in base_qc.data:
        op = instr.operation
        if op.name == GATE_NAME:
            gamma *= float(np.sum(np.abs(list(pec_inverse_weights(p2, 2).values()))))
        elif op.name in ("gpi", "gpi2"):
            gamma *= float(np.sum(np.abs(list(pec_inverse_weights(p1, 1).values()))))
    return gamma


def measure(qc, p2, p1, Pmat):
    dm = noisy_dm(qc, p2, p1)
    return float(np.real(np.trace(Pmat @ dm.data)))


def run_case(case, p, fixed_solutions):
    slot, group, label = case["slot"], case["group"], case["label"]
    base = optimized_native_circuit(fixed_solutions[slot]["angles"], GATE_NAME)
    basis_qc = native_basis_change(effrag_mod.combined_basis_label(group), GATE_NAME)
    full_base = base.compose(basis_qc)
    Pmat = np.asarray(Pauli(label).to_matrix())
    gamma = fixed_gamma(full_base, P2, P1)
    exact_val = float(np.real(Statevector.from_instruction(full_base).expectation_value(Pauli(label))))

    print(f"\n  === {case['tag']} ===")
    print(f"    exact={exact_val:+.4f}  gamma={gamma:.4f}")

    # -- PILOT 1: sign distribution, cheap, no DM sim --
    rng_pilot = np.random.default_rng(stable_seed("task35_pilot_sign", slot, label))
    signs_pilot = np.zeros(N_PILOT_SIGN, dtype=int)
    for i in range(N_PILOT_SIGN):
        _, s = sample_twirl_choice(full_base, P2, P1, rng_pilot)
        signs_pilot[i] = s
    p_plus = float((signs_pilot > 0).mean())
    p_minus = 1.0 - p_plus
    se_p_plus = float(np.sqrt(p_plus * p_minus / N_PILOT_SIGN))
    print(f"    PILOT (N={N_PILOT_SIGN}, sign-only, no DM sim): p_plus={p_plus:.4f} +/- {se_p_plus:.4f}")

    # -- PILOT 2: per-stratum std(m), expensive (DM sim), used ONLY for Neyman allocation --
    m_plus_pilot, m_minus_pilot = [], []
    rng_pv = np.random.default_rng(stable_seed("task35_pilot_var", slot, label))
    tries = 0
    while (len(m_plus_pilot) < N_PILOT_VARIANCE or len(m_minus_pilot) < N_PILOT_VARIANCE) and tries < 20 * N_PILOT_VARIANCE:
        tries += 1
        qc, s = sample_twirl_choice(full_base, P2, P1, rng_pv)
        if s > 0 and len(m_plus_pilot) < N_PILOT_VARIANCE:
            m_plus_pilot.append(measure(qc, P2, P1, Pmat))
        elif s < 0 and len(m_minus_pilot) < N_PILOT_VARIANCE:
            m_minus_pilot.append(measure(qc, P2, P1, Pmat))
    sigma_plus = float(np.std(m_plus_pilot, ddof=1)) if len(m_plus_pilot) > 1 else 0.0
    sigma_minus = float(np.std(m_minus_pilot, ddof=1)) if len(m_minus_pilot) > 1 else 0.0
    print(f"    PILOT variance: sigma(m|+)={sigma_plus:.4f} (n={len(m_plus_pilot)})  "
          f"sigma(m|-)={sigma_minus:.4f} (n={len(m_minus_pilot)})")

    # -- production allocations (frozen BEFORE any production trial, per this project's own
    # anti-adaptive-bias discipline) --
    n_plus_prop = max(1, round(N_PRODUCTION * p_plus)) if p_minus > 0 else N_PRODUCTION
    n_minus_prop = N_PRODUCTION - n_plus_prop
    total_w = p_plus * sigma_plus + p_minus * sigma_minus
    if total_w > 0 and p_minus > 1e-6:
        n_plus_ney = max(1, round(N_PRODUCTION * (p_plus * sigma_plus) / total_w))
        n_minus_ney = N_PRODUCTION - n_plus_ney
    else:
        n_plus_ney, n_minus_ney = n_plus_prop, n_minus_prop
    print(f"    allocation (N={N_PRODUCTION}): proportional n+={n_plus_prop} n-={n_minus_prop}  "
          f"Neyman n+={n_plus_ney} n-={n_minus_ney}")

    def draw_n_with_sign(target_sign, n, rng, max_tries=2000):
        got = []
        tries = 0
        while len(got) < n and tries < max_tries * n:
            tries += 1
            qc, s = sample_twirl_choice(full_base, P2, P1, rng)
            if (s > 0) == (target_sign > 0):
                got.append(measure(qc, P2, P1, Pmat))
        return got

    def one_trial(seed):
        rng = np.random.default_rng(seed)
        # (a) IID
        iid_signed = []
        for _ in range(N_PRODUCTION):
            qc, s = sample_twirl_choice(full_base, P2, P1, rng)
            m = measure(qc, P2, P1, Pmat)
            iid_signed.append(s * gamma * m)
        E_iid = float(np.mean(iid_signed))

        # (c) control variate: V = sign - E[sign] = sign - (p_plus-p_minus), on a FRESH iid sample
        # (paired sign/m tracked directly, avoiding any sign-from-value ambiguity)
        e_sign = p_plus - p_minus
        rng_cv = np.random.default_rng(seed + 500_000)
        pairs = []
        for _ in range(N_PRODUCTION):
            qc, s = sample_twirl_choice(full_base, P2, P1, rng_cv)
            m = measure(qc, P2, P1, Pmat)
            pairs.append((s, m))
        vals = np.array([s * gamma * m for s, m in pairs])
        V = np.array([s - e_sign for s, m in pairs])
        if V.var() > 1e-12:
            beta = float(np.cov(vals, V, ddof=1)[0, 1] / V.var(ddof=1))
        else:
            beta = 0.0
        E_cv = float(np.mean(vals - beta * V))

        # (b) stratified, proportional allocation
        rng_s1 = np.random.default_rng(seed + 1_000_000)
        m_plus = draw_n_with_sign(+1, n_plus_prop, rng_s1) if n_plus_prop > 0 else []
        m_minus = draw_n_with_sign(-1, n_minus_prop, rng_s1) if n_minus_prop > 0 else []
        mean_plus = np.mean(m_plus) if m_plus else 0.0
        mean_minus = np.mean(m_minus) if m_minus else 0.0
        E_strat = gamma * (p_plus * mean_plus - p_minus * mean_minus)

        # (d) stratified + control variate: apply CV correction WITHIN each stratum using the
        # analytic within-stratum mean (trivial here since stratified draws are already conditioned
        # on their sign -- the CV correction on sign is a no-op within a pure stratum by construction;
        # report as identical to (b) for this design, disclosed rather than fabricated as new)
        E_strat_cv = E_strat

        return E_iid, E_cv, E_strat, E_strat_cv

    print(f"    running {N_TRIALS} independent production trials (N={N_PRODUCTION} each)...")
    results = {"iid": [], "cv": [], "strat": [], "strat_cv": []}
    for t in range(N_TRIALS):
        e_iid, e_cv, e_strat, e_strat_cv = one_trial(70_000 + t)
        results["iid"].append(e_iid)
        results["cv"].append(e_cv)
        results["strat"].append(e_strat)
        results["strat_cv"].append(e_strat_cv)

    summary = {}
    for method, vals in results.items():
        arr = np.array(vals)
        bias = float(arr.mean() - exact_val)
        std = float(arr.std(ddof=1))
        mse = bias ** 2 + std ** 2
        summary[method] = {"bias": bias, "std": std, "mse": mse}
        print(f"    {method:<10} bias={bias:+.4f}  std={std:.4f}  MSE={mse:.4f}")

    std_iid = summary["iid"]["std"]
    for method in ["cv", "strat"]:
        pct = 100 * (std_iid - summary[method]["std"]) / std_iid if std_iid > 0 else 0.0
        print(f"    {method} std reduction vs IID: {pct:+.1f}%")
        summary[method]["std_reduction_vs_iid_pct"] = pct

    return {
        "exact_val": exact_val, "gamma": gamma, "p_plus": p_plus, "p_minus": p_minus,
        "sigma_plus": sigma_plus, "sigma_minus": sigma_minus,
        "allocation_proportional": [n_plus_prop, n_minus_prop], "allocation_neyman": [n_plus_ney, n_minus_ney],
        "summary": summary,
    }


def main():
    print("\n" + "=" * 96)
    print("  task35bc_stratified_control_variate_pec.py -- sign-stratified PEC + sign control variate")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)

    all_results = {}
    for case in TEST_CASES:
        all_results[case["tag"]] = run_case(case, p, fixed_solutions)

    print(f"\n  -- CROSS-CASE SUMMARY --")
    for tag, r in all_results.items():
        print(f"    {tag}: p_plus={r['p_plus']:.4f}  "
              f"CV std reduction={r['summary']['cv']['std_reduction_vs_iid_pct']:+.1f}%  "
              f"Strat std reduction={r['summary']['strat']['std_reduction_vs_iid_pct']:+.1f}%")

    with open(RESULTS_PATH, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
