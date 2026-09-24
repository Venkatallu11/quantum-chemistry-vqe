#!/usr/bin/env python3
"""
task58_blind_h4_vqe.py -- iteration 58. A genuinely BLIND H4 forged-energy
estimator, built in direct response to an IonQ reviewer's Sep 22
2026 methodology review: the project's existing 0.0105-0.0192 kcal/mol
headline result relies on state-preparation circuits whose angles are
fit directly to the classically pre-computed EXACT FCI Schmidt vectors --
"the preparation angles are also fitted to the exact Schmidt vectors, so
there is no variational step." Verified directly against this project's
own code (task40_certification_ablation_adversarial.py,
phys_constrained_reconstruction.build_P_S) -- the critique is accurate.

THE RULE, enforced throughout: no exact FCI/Schmidt information enters
circuit construction, the reconstruction basis, or the optimization
objective. FCI is used ONLY offline, after a result is produced, to
report the error -- never fed back into circuit-building or parameter
selection.

ARCHITECTURE (v2 -- see module history below for what changed and why):
  1. Two independent 15-parameter orthogonal frames, U_alpha(theta_a) and
     U_beta(theta_b) = U0 * expm(A(theta)), applied to the PLAIN
     COMPUTATIONAL weight-2 basis (U0 = identity), NOT the exact Schmidt
     vectors. 15 parameters fully spans the space of 6-dim orthonormal
     frames (a mathematical fact, not something needing FCI verification).
     Both registers are optimized independently -- this project's existing
     beta=sign*alpha shortcut is itself derived from the exact answer
     (verified true for it, never assumed for an unknown one), so a
     genuinely blind estimator cannot use it; the real, disclosed cost is
     roughly double the circuits of the oracle-informed pipeline.
  2. Each frame vector (and pairwise sums, for cross-term phase circuits)
     becomes a REAL quantum circuit via native_stateprep's already-
     validated (<1e-15) generic real-state-prep machinery -- nothing here
     is fit to a known target.
  3. Measurement grouping: GENERAL-COMMUTING (GC) groups, not qubit-wise
     commuting -- reuses this project's own already-verified GC machinery
     (iteration 46-49, general_commuting_measurements.py) to cut circuits
     per register from 273 (13 QWC groups x 21 slots) to 84 (4 GC groups
     x 21 slots), a real ~3.25x reduction, unrelated to the oracle
     question (purely a measurement-efficiency choice already proven safe).
  4. Energy reconstruction: raw measured Pauli expectations go directly
     into the K x K alpha/beta matrices (entanglement_forging_h4
     .ef_energy_from_noisy_matrices) -- NO P_S, NO reference to U_exact.
  5. Optimizer: SPSA (Simultaneous Perturbation Stochastic Approximation),
     not COBYLA. This is a genuine, disclosed engineering fix, not a
     shortcut around the science: COBYLA needs >= n_params+2 = 32
     evaluations just to build its initial model before any real
     optimization begins, and scipy's COBYLA implementation exposes no
     resumable internal state -- on this machine (repeated background-
     process kills from memory pressure, observed constantly this
     session), a killed COBYLA run loses 100% of its progress and must
     restart from theta=0. SPSA needs exactly 2 real evaluations per
     step, and its ENTIRE state is the current theta vector + iteration
     count -- trivially checkpointed to disk after every step and resumed
     exactly where it left off. A first real COBYLA attempt (this
     project's own dev history, not repeated here) confirmed the problem
     empirically: 10 real evaluations / 108 minutes / 5,460 circuits
     completed before a memory-kill, still inside COBYLA's un-converged
     warm-up phase, and none of that progress was resumable.

VALIDATION, mandatory before trusting anything: the architecture's math
(15-param frame expressivity, real-circuit fidelity, blind-reconstruction
correctness) was verified OFFLINE ONLY in scratchpad -- fitting theta* to
reproduce the exact answer and confirming a 0.000000 kcal/mol round-trip
energy -- BEFORE any real optimizer run. That offline theta*/exact-answer
information is not visible to anything in this file.

Run:
    PYTHONHASHSEED=0 python vqe/task58_blind_h4_vqe.py --backend ideal --iters 40
    PYTHONHASHSEED=0 python vqe/task58_blind_h4_vqe.py --backend ideal --iters 40   # resumes automatically
"""
import os
import sys
import json
import time
import argparse
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment
from entanglement_forging_h4 import ef_energy_from_noisy_matrices
from task27c_full_h4_folds import kept_slots_for_K
from native_stateprep import build_real_state_prep_circuit, to_native
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list
from general_commuting_measurements import (
    build_general_commuting_measurement_plan, measurement_circuit_general, expectations_from_counts,
)

HARTREE_TO_KCAL_MOL = 627.5094740631
K = 6
N_THETA = K * (K - 1) // 2  # 15
N_REGISTER = 4
WEIGHT2_INDICES = [3, 5, 6, 9, 10, 12]
GATE_NAME = "zz"
SHOTS = 2000
CKPT_PATH = os.path.join(os.path.dirname(__file__), "task58_blind_h4_vqe_ckpt.json")
LOG_PATH = os.path.join(os.path.dirname(__file__), "task58_blind_h4_vqe_log.json")

# SPSA standard hyperparameters (Spall's conventions)
SPSA_a = 0.15
SPSA_c = 0.15
SPSA_A = 5.0
SPSA_alpha = 0.602
SPSA_gamma = 0.101


def skew_from_theta(theta, K):
    A = np.zeros((K, K))
    idx = 0
    for i in range(K):
        for j in range(i + 1, K):
            A[i, j] = theta[idx]
            A[j, i] = -theta[idx]
            idx += 1
    return A


def U_from_theta(theta, K):
    return expm(skew_from_theta(theta, K))  # U0 = identity, the plain computational basis


def embed_weight2(vec6):
    full = np.zeros(16)
    full[WEIGHT2_INDICES] = vec6
    return full


def frame_targets(theta):
    U6 = U_from_theta(theta, K)
    targets = {f"u_{n}": embed_weight2(U6[:, n]) for n in range(K)}
    for n in range(K):
        for m in range(K):
            if n < m:
                targets[f"(u{n}+u{m})"] = embed_weight2((U6[:, n] + U6[:, m]) / np.sqrt(2))
    return targets


def setup_problem():
    """Only reads the Hamiltonian's own terms/lambdas/enuc/alpha_labels/
    identity_label -- NEVER touches p['u_vecs'] (exact Schmidt vectors) or
    p['signs'] anywhere in this file."""
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    return {
        "terms": p["terms"], "lambdas": p["lambdas"], "enuc": p["enuc"],
        "alpha_labels": p["alpha_labels"], "identity_label": p["identity_label"],
        "exact_energy_OFFLINE_ONLY": p["exact_energy"],
    }


def measure_register(vecs_by_slot, kept, groups, diagonalizers, diag_natives, backend, backend_name, shots):
    out = {name: {} for name in kept}
    n_circuits = 0
    for name in kept:
        base = build_real_state_prep_circuit(vecs_by_slot[name])
        base_native = to_native(base, GATE_NAME)
        circuits = []
        for diag, diag_native in zip(diagonalizers, diag_natives):
            # Compose two ALREADY-native circuits, then add measurement LAST --
            # never call to_native() on a circuit that already carries a
            # measure instruction (the native gate Target has no synthesis
            # rule for "measure", which raised TranspilerError when this was
            # tried the naive way -- verified and fixed here, not assumed).
            qc = base_native.compose(diag_native, qubits=list(range(N_REGISTER)))
            qc.measure_all()
            circuits.append(qc)
        job = submit_job(circuits, backend, backend_name, shots=shots)
        counts = get_counts_list(job)
        n_circuits += len(circuits)
        for diag, c in zip(diagonalizers, counts):
            exp = expectations_from_counts(c, diag)
            out[name].update(exp)
    return out, n_circuits


def energy_from_raw(problem, raw_alpha, raw_beta, kept):
    K_ = K
    alpha_mats = {l: np.zeros((K_, K_)) for l in problem["alpha_labels"]}
    beta_mats = {l: np.zeros((K_, K_)) for l in problem["alpha_labels"]}
    for l in problem["alpha_labels"]:
        if l == problem["identity_label"]:
            for n in range(K_):
                alpha_mats[l][n, n] = 1.0
                beta_mats[l][n, n] = 1.0
            continue
        for n in range(K_):
            name = f"u_{n}"
            alpha_mats[l][n, n] = raw_alpha[name].get(l, 0.0)
            beta_mats[l][n, n] = raw_beta[name].get(l, 0.0)
        for n in range(K_):
            for m in range(K_):
                if n < m:
                    name = f"(u{n}+u{m})"
                    va = raw_alpha[name].get(l, 0.0)
                    vb = raw_beta[name].get(l, 0.0)
                    alpha_mats[l][n, m] = alpha_mats[l][m, n] = va - 0.5 * (alpha_mats[l][n, n] + alpha_mats[l][m, m])
                    beta_mats[l][n, m] = beta_mats[l][m, n] = vb - 0.5 * (beta_mats[l][n, n] + beta_mats[l][m, m])
    E = ef_energy_from_noisy_matrices(problem["terms"], problem["lambdas"], alpha_mats, beta_mats, problem["enuc"], K_)
    return E


def load_checkpoint():
    if os.path.exists(CKPT_PATH):
        with open(CKPT_PATH) as f:
            d = json.load(f)
        return np.array(d["theta"]), d["k"], d.get("history", [])
    return np.zeros(2 * N_THETA), 0, []


def save_checkpoint(theta, k, history):
    with open(CKPT_PATH, "w") as f:
        json.dump({"theta": theta.tolist(), "k": k, "history": history}, f, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="ideal")
    ap.add_argument("--iters", type=int, default=40, help="additional SPSA iterations to run this invocation")
    args = ap.parse_args()

    problem = setup_problem()
    non_id_labels = sorted(l for l in problem["alpha_labels"] if l != problem["identity_label"])
    diag, plus, kept = kept_slots_for_K(K)
    groups, diagonalizers = build_general_commuting_measurement_plan(non_id_labels)
    diag_natives = [to_native(dg.to_circuit(), GATE_NAME) for dg in diagonalizers]
    print(f"  {len(non_id_labels)} labels -> {len(groups)} GC groups, {len(kept)} slots -> "
          f"{len(kept)*len(groups)} circuits PER REGISTER, x2 per SPSA eval, x2 evals/step "
          f"= {4*len(kept)*len(groups)} circuits/SPSA iteration")

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}, target noise_model={args.backend}")

    theta, k0, history = load_checkpoint()
    if k0 > 0:
        print(f"  RESUMING from checkpoint: k={k0}, last E={history[-1]['E_ha']:.6f} Ha")
    else:
        print(f"  Starting fresh: theta=0 (no oracle info)")

    t0 = time.time()

    def objective(th):
        theta_a = th[:N_THETA]
        theta_b = th[N_THETA:]
        targets_a = frame_targets(theta_a)
        targets_b = frame_targets(theta_b)
        raw_a, na = measure_register(targets_a, kept, groups, diagonalizers, diag_natives, backend, args.backend, SHOTS)
        raw_b, nb = measure_register(targets_b, kept, groups, diagonalizers, diag_natives, backend, args.backend, SHOTS)
        E = energy_from_raw(problem, raw_a, raw_b, kept)
        return E, na + nb

    rng = np.random.default_rng(1234 + k0)  # advance the stream on resume, not identical replay
    n_circ_total = 0
    for step in range(args.iters):
        k = k0 + step
        c_k = SPSA_c / (k + 1) ** SPSA_gamma
        a_k = SPSA_a / (k + 1 + SPSA_A) ** SPSA_alpha
        delta = rng.choice([-1.0, 1.0], size=2 * N_THETA)

        E_plus, n1 = objective(theta + c_k * delta)
        E_minus, n2 = objective(theta - c_k * delta)
        n_circ_total += n1 + n2
        ghat = (E_plus - E_minus) / (2 * c_k) * (1.0 / delta)
        theta = theta - a_k * ghat

        E_center, n3 = objective(theta)  # honest checkpoint of current best, real data
        n_circ_total += n3
        elapsed = time.time() - t0
        print(f"    k={k}: E+={E_plus:.6f} E-={E_minus:.6f}  E(theta_new)={E_center:.6f} Ha  "
              f"[{n1+n2+n3} circuits, {elapsed:.1f}s elapsed]")
        history.append({"k": k, "E_plus": E_plus, "E_minus": E_minus, "E_ha": E_center,
                         "n_circuits": n1 + n2 + n3, "elapsed_s": elapsed})
        save_checkpoint(theta, k + 1, history)
        with open(LOG_PATH, "w") as f:
            json.dump({"history": history}, f, indent=2)

    print(f"\n  {args.iters} SPSA iterations done this invocation, {n_circ_total} circuits, "
          f"{time.time()-t0:.1f}s. Total k so far: {k0 + args.iters}")

    exact = problem["exact_energy_OFFLINE_ONLY"]
    last_E = history[-1]["E_ha"] if history else None
    if last_E is not None:
        err_kcal = abs(last_E - exact) * HARTREE_TO_KCAL_MOL
        print(f"  OFFLINE comparison (never used during optimization): exact={exact:.6f} Ha, "
              f"current E={last_E:.6f} Ha, err={err_kcal:.4f} kcal/mol")


if __name__ == "__main__":
    main()
