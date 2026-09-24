#!/usr/bin/env python3
"""
task61_h4_no_frame_gls.py -- covariance-aware no-frame H4 reconstruction.

This is intentionally narrower than a new mitigation method:
  * no joint Schmidt-frame fit
  * no exact-energy optimization/selection
  * same ancilla-parity + conditioned PEC + GPi2 correction
  * same exact-FCI Schmidt-coordinate basis disclosure as Task 60
  * same 5-parameter independent physical state per slot

NEW STATISTICAL PIECE:
  Pauli expectation values measured in one commuting group are correlated
  because they come from the same accepted bitstrings. Task 60 treated them
  as independent. Here we use the full multinomial covariance block for each
  group and generalized least squares (GLS).

The GPi2 value passed on the command line must be selected by the previously
frozen measurement-only validation, not by final exact energy.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(__file__))

from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from phys_constrained_reconstruction import build_P_S
from ionq_simulator_binding_curve import expectation_from_counts
from task60_vadim_no_frame_h4 import (
    K,
    GATE_NAME,
    ZZ_ASSUMED,
    angles_to_vector,
    vector_to_angles,
)

FIT_RESTARTS = 8
FIT_FLOOR = 1e-12
JITTER_REL = 1e-8


def pauli_sign(bs, label):
    """Sign of a Pauli label on a measured Z-basis bitstring."""
    bits = bs[::-1]
    sign = 1
    for q in range(len(label)):
        ch = label[len(label) - 1 - q]
        if ch != "I" and q < len(bits) and bits[q] == "1":
            sign = -sign
    return sign


def pooled_kept_counts(checkpoints, backend, name, gi):
    pooled = {}
    for path in checkpoints:
        with open(path, "r", encoding="utf-8") as f:
            state = json.load(f)
        entry = state["done"][f"{backend}|{name}"]
        for bs, count in entry["counts"][gi].items():
            pooled[bs] = pooled.get(bs, 0) + int(count)

    # Ancilla=0 is the same physical postselection used by Task 60.
    return {bs[1:]: c for bs, c in pooled.items() if bs[0] == "0"}


def group_stats(kept_counts, labels):
    total = sum(kept_counts.values())
    if total <= 0:
        raise RuntimeError("empty postselected group")

    probs = {bs: c / total for bs, c in kept_counts.items()}
    means = np.array([
        sum(p * pauli_sign(bs, label) for bs, p in probs.items())
        for label in labels
    ], dtype=float)

    X = np.array([
        [pauli_sign(bs, label) for label in labels]
        for bs in probs
    ], dtype=float)
    pvec = np.array([probs[bs] for bs in probs], dtype=float)
    second = (X.T * pvec) @ X
    cov = (second - np.outer(means, means)) / total
    cov = 0.5 * (cov + cov.T)
    return means, cov, int(total)


def whiten_block(cov):
    cov = 0.5 * (cov + cov.T)
    scale = max(float(np.max(np.diag(cov))), 1e-12)
    reg = JITTER_REL * scale
    cov_reg = cov + reg * np.eye(cov.shape[0])
    vals, vecs = np.linalg.eigh(cov_reg)
    vals = np.maximum(vals, reg)
    return vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.T


def load_corrected_gls_data(checkpoints, backend, kept, fixed_solutions, non_id, p_gpi2, P_S):
    corrected = {name: {} for name in kept}
    whitener = {name: {} for name in kept}
    group_labels = {name: {} for name in kept}

    for name in kept:
        # Each submission contains the same group ordering.
        with open(checkpoints[0], "r", encoding="utf-8") as f:
            state0 = json.load(f)
        entry0 = state0["done"][f"{backend}|{name}"]

        A, B, _, _ = analytic_A_and_B_conditioned(
            fixed_solutions[name]["angles"],
            GATE_NAME,
            ZZ_ASSUMED,
            GPI_REAL_MEAN,
            float(p_gpi2),
            0.0,
            0.0,
            non_id,
        )

        for gi, group in enumerate(entry0["groups"]):
            labels = [l for l in group if l != "IIII"]
            if not labels:
                continue

            kept_counts = pooled_kept_counts(checkpoints, backend, name, gi)
            raw_mean, raw_cov, _ = group_stats(kept_counts, labels)

            ratios = np.array([
                float(B[l] / A[l]) if abs(A[l]) > 1e-6 else 1.0
                for l in labels
            ])
            corrected_mean = raw_mean * ratios
            corrected_cov = np.diag(ratios) @ raw_cov @ np.diag(ratios)

            corrected[name].update({
                l: float(m) for l, m in zip(labels, corrected_mean)
            })
            whitener[name][gi] = whiten_block(corrected_cov)
            group_labels[name][gi] = labels

    return corrected, whitener, group_labels


def spectral_initializer(P_S, measured):
    M = np.zeros((K, K), dtype=float)
    for l, m in measured.items():
        P = np.real_if_close(np.asarray(P_S[l])).astype(float)
        M += float(m) * P
    M = 0.5 * (M + M.T)
    vals, vecs = np.linalg.eigh(M)
    v = vecs[:, int(np.argmax(vals))]
    v = np.asarray(v, dtype=float)
    v /= max(np.linalg.norm(v), FIT_FLOOR)
    return v


def fit_slot_gls(P_S, corrected, whiteners, group_labels, seed):
    rng = np.random.default_rng(seed)
    theta0 = vector_to_angles(spectral_initializer(P_S, corrected))

    def residual(theta):
        a = angles_to_vector(theta)
        chunks = []
        for gi, labels in group_labels.items():
            r = np.array([
                float(a @ np.real_if_close(np.asarray(P_S[l])).astype(float) @ a)
                - float(corrected[l])
                for l in labels
            ])
            chunks.append(whiteners[gi] @ r)
        return np.concatenate(chunks)

    inits = [theta0, np.zeros(K - 1, dtype=float)]
    inits += [rng.uniform(-0.5 * np.pi, 0.5 * np.pi, K - 1)
              for _ in range(max(0, FIT_RESTARTS - 2))]

    best = None
    for x0 in inits:
        try:
            sol = least_squares(
                residual,
                x0,
                method="trf",
                bounds=(-0.5 * np.pi, 0.5 * np.pi),
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
                max_nfev=1800,
            )
        except Exception:
            continue
        sse = float(np.sum(sol.fun ** 2))
        if np.isfinite(sse) and (best is None or sse < best[0]):
            best = (sse, np.asarray(sol.x, dtype=float))

    if best is None:
        raise RuntimeError("GLS fit failed")
    return angles_to_vector(best[1]), best[0]


def build_full_from_states(states, diag, K, non_id):
    full = {name: {} for name in diag}
    for name in diag:
        a = states[name]
        for label in non_id:
            P = np.real_if_close(np.asarray(P_S_GLOBAL[label])).astype(float)
            full[name][label] = float(a @ P @ a)

    for n in range(K):
        for m in range(n + 1, K):
            plus = f"(u{n}+u{m})"
            minus = f"(u{n}-u{m})"
            ap = states[plus]
            full[plus] = {}
            full[minus] = {}
            for label in non_id:
                P = np.real_if_close(np.asarray(P_S_GLOBAL[label])).astype(float)
                full[plus][label] = float(ap @ P @ ap)
                full[minus][label] = (
                    full[f"u_{n}"][label] + full[f"u_{m}"][label] - full[plus][label]
                )
    return full


P_S_GLOBAL = {}


def run(backend, checkpoints, p_gpi2):
    global P_S_GLOBAL

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, _ = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == len(p["targets"])
    diag, _, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S_GLOBAL = build_P_S(p["alpha_labels"], U_exact)

    corrected, whiteners, group_labels = load_corrected_gls_data(
        checkpoints, backend, kept, fixed_solutions, non_id, p_gpi2, P_S_GLOBAL
    )

    states = {}
    chi2 = []
    for i, name in enumerate(kept):
        a, sse = fit_slot_gls(
            P_S_GLOBAL,
            corrected[name],
            whiteners[name],
            group_labels[name],
            seed=6100 + i,
        )
        states[name] = a
        n_obs = sum(len(v) for v in group_labels[name].values())
        chi2.append(sse / max(1, n_obs - (K - 1)))

    full = build_full_from_states(states, diag, K, non_id)
    mats = combine_matrices(
        full,
        p["alpha_labels"],
        p["identity_label"],
        K,
    )
    E, errs = energy_from_alpha_matrices(
        mats,
        p["terms"],
        p["lambdas"],
        p["enuc"],
        p["signs"],
        K,
        exact_energy=p["exact_energy"],
        noiseless_energy=p["noiseless_energy"],
    )

    result = {
        "backend": backend,
        "p_gpi2": float(p_gpi2),
        "pooled_draws": len(checkpoints),
        "shots_per_circuit": 20000 * len(checkpoints),
        "energy_ha": float(E),
        "err_vs_exact_kcal": float(errs["err_vs_exact_kcal"]),
        "mean_gls_chi2_dof": float(np.mean(chi2)),
        "max_gls_chi2_dof": float(np.max(chi2)),
        "uses_joint_schmidt_frame": False,
        "uses_exact_energy_for_selection": False,
        "uses_exact_schmidt_basis": True,
        "uses_exact_target_preparation": True,
    }
    out = os.path.join(
        os.path.dirname(__file__),
        f"task61_h4_no_frame_gls_{backend}_result.json",
    )
    with open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["aria-1", "forte-1"], required=True)
    ap.add_argument("--p-gpi2", type=float, required=True)
    ap.add_argument("--checkpoint", required=True)
    args = ap.parse_args()
    checkpoints = [x.strip() for x in args.checkpoint.split(",") if x.strip()]
    for path in checkpoints:
        if not os.path.exists(path):
            raise SystemExit(f"checkpoint not found: {path}")
    run(args.backend, checkpoints, args.p_gpi2)


if __name__ == "__main__":
    main()
