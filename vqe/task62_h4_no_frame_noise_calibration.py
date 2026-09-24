#!/usr/bin/env python3
"""
Task 62 -- H4 no-frame, measurement-only recalibration of (p_ZZ, p_GPi2).

Scope is deliberately narrow:
  - NO shared Schmidt-frame fit.
  - NO exact-energy objective.
  - NO new quantum circuits.
  - SAME ancilla-parity + conditioned PEC model.
  - SAME exact-FCI Schmidt-coordinate operators P_S (disclosed oracle-assisted
    basis, but no fitted frame).
  - Each of the 21 H4 slots is fitted independently to a 5-angle real pure
    state.
  - Only p_ZZ and p_GPi2 are scanned; selection uses held-out Pauli residuals.

p_ZZ prior comes from the real calibration: 0.014593 +/- 0.000124.
p_GPi2 remains locally unresolved, so the search stays in the empirically
relevant 0.0003--0.0007 region identified by prior measurement-only tests.

Final exact energy is printed ONLY as an informational validation number.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from phys_constrained_reconstruction import build_P_S
from task60_ionq_no_frame_h4 import (
    K,
    GATE_NAME,
    fit_slot_data_only,
    fit_and_score_candidate,
    split_labels,
    load_pooled_postselected,
    variance_weights,
)

ZZ_MEAN = 0.014593
ZZ_STD = 0.000124

# Fine but finite, pre-registered grids. No exact energy enters selection.
PZZ_GRID = np.array([ZZ_MEAN + ZZ_STD * s for s in (-3,-2,-1,0,1,2,3)], dtype=float)
GPI2_GRID = np.arange(0.00030, 0.000701, 0.000050, dtype=float)

VAL_FRACTION = 0.30
FIT_RESTARTS = 4

_AB_CACHE = {}


def corrected_for_candidate(postselected, pzz, pgpi2, kept, non_id, fixed_solutions):
    out = {name: {} for name in kept}
    key_base = (round(float(pzz), 8), round(float(pgpi2), 8))
    for name in kept:
        key = (name,) + key_base
        if key not in _AB_CACHE:
            A, B, _, _ = analytic_A_and_B_conditioned(
                fixed_solutions[name]["angles"],
                GATE_NAME,
                float(pzz),
                GPI_REAL_MEAN,
                float(pgpi2),
                0.0,
                0.0,
                non_id,
            )
            _AB_CACHE[key] = (A, B)
        A, B = _AB_CACHE[key]
        for label in non_id:
            denom = float(A[label])
            ratio = float(B[label] / denom) if abs(denom) > 1e-6 else 1.0
            out[name][label] = float(
                np.clip(postselected[name][label] * ratio, -1.0, 1.0)
            )
    return out


def score_candidate(corrected, weights, P_S, kept, train_labels, val_labels):
    total = 0.0
    n = 0
    # Fit each slot independently on train labels only, then predict held-out labels.
    for i, name in enumerate(kept):
        train_m = {l: corrected[name][l] for l in train_labels[name]}
        train_w = {l: weights[name][l] for l in train_labels[name]}
        fit = fit_slot_data_only(
            P_S,
            train_m,
            train_w,
            seed=62000 + i,
            n_restarts=FIT_RESTARTS,
        )
        a = fit.vector
        for label in val_labels[name]:
            P = np.real_if_close(np.asarray(P_S[label])).astype(float)
            r = float(a @ P @ a) - float(corrected[name][label])
            total += float(weights[name][label]) * r * r
            n += 1
    dof = max(1, n - len(kept) * (K - 1))
    return float(total / dof)


def build_full(fits, diag, non_id, P_S):
    full = {name: {} for name in diag}
    for name in diag:
        a = fits[name].vector
        for l in non_id:
            P = np.real_if_close(np.asarray(P_S[l])).astype(float)
            full[name][l] = float(a @ P @ a)

    for n in range(K):
        for m in range(n + 1, K):
            plus = f"(u{n}+u{m})"
            minus = f"(u{n}-u{m})"
            ap = fits[plus].vector
            full[plus] = {}
            full[minus] = {}
            for l in non_id:
                P = np.real_if_close(np.asarray(P_S[l])).astype(float)
                full[plus][l] = float(ap @ P @ ap)
                full[minus][l] = (
                    full[f"u_{n}"][l] + full[f"u_{m}"][l] - full[plus][l]
                )
    return full


def run(backend, checkpoints):
    p = setup_fragment([0,1,2,3], nelec=4, d=1.0, K=K, strict=True)
    non_id = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, _ = fit_all_targets(p["targets"], tol=1e-10)
    assert n_ok == len(p["targets"])

    diag, _, kept = kept_slots_for_K(K)
    P_S = build_P_S(
        p["alpha_labels"],
        np.asarray(p["u_vecs"]).T,
    )

    train_labels, val_labels = split_labels(
        kept,
        non_id,
        seed=62,
        val_fraction=VAL_FRACTION,
    )

    postselected, kept_shots = load_pooled_postselected(
        checkpoints, kept, backend
    )
    weights = variance_weights(postselected, kept_shots)

    rows = []
    best = None
    for pzz in PZZ_GRID:
        for pgpi2 in GPI2_GRID:
            corrected = corrected_for_candidate(
                postselected, pzz, pgpi2, kept, non_id, fixed_solutions
            )
            val_score = score_candidate(
                corrected, weights, P_S, kept, train_labels, val_labels
            )
            row = {
                "p_zz": float(pzz),
                "p_gpi2": float(pgpi2),
                "heldout_chi2_dof": float(val_score),
            }
            rows.append(row)
            if best is None or val_score < best["heldout_chi2_dof"]:
                best = row
                print(
                    f"  backend={backend} NEW BEST pZZ={pzz:.8f} "
                    f"pGPi2={pgpi2:.7f} heldout_chi2/dof={val_score:.6g}",
                    flush=True,
                )

    # Final fit after parameter selection, using ALL available labels.
    corrected_best = corrected_for_candidate(
        postselected, best["p_zz"], best["p_gpi2"],
        kept, non_id, fixed_solutions
    )
    fits = {}
    for i, name in enumerate(kept):
        fits[name] = fit_slot_data_only(
            P_S, corrected_best[name], weights[name],
            seed=63000 + i, n_restarts=FIT_RESTARTS
        )

    full = build_full(fits, diag, non_id, P_S)
    mats = combine_matrices(
        full, p["alpha_labels"], p["identity_label"], K
    )
    E, errs = energy_from_alpha_matrices(
        mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
        exact_energy=p["exact_energy"],
        noiseless_energy=p["noiseless_energy"],
    )

    result = {
        "backend": backend,
        "pooled_draws": len(checkpoints),
        "p_zz": best["p_zz"],
        "p_gpi2": best["p_gpi2"],
        "heldout_chi2_dof": best["heldout_chi2_dof"],
        "final_energy_ha": float(E),
        "err_vs_exact_kcal": float(errs["err_vs_exact_kcal"]),
        "strict_chemical_accuracy_0p25": bool(abs(errs["err_vs_exact_kcal"]) <= 0.25),
        "stretch_0p02": bool(abs(errs["err_vs_exact_kcal"]) <= 0.02),
        "n_candidates": len(rows),
        "pzz_grid": PZZ_GRID.tolist(),
        "gpi2_grid": GPI2_GRID.tolist(),
        "uses_joint_schmidt_frame": False,
        "uses_exact_energy_for_selection": False,
        "uses_exact_schmidt_basis": True,
        "uses_exact_target_preparation": True,
        "all_candidates": rows,
    }

    out = os.path.join(
        os.path.dirname(__file__),
        f"task62_h4_no_frame_noise_calibration_{backend}_result.json",
    )
    with open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print("\nFINAL INFORMATIONAL ENERGY CHECK")
    print(json.dumps({
        k: result[k] for k in (
            "backend","p_zz","p_gpi2","heldout_chi2_dof",
            "final_energy_ha","err_vs_exact_kcal",
            "strict_chemical_accuracy_0p25","stretch_0p02"
        )
    }, indent=2))

    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["aria-1","forte-1"], required=True)
    ap.add_argument("--checkpoint", required=True)
    args = ap.parse_args()
    checkpoints = [p.strip() for p in args.checkpoint.split(",") if p.strip()]
    for path in checkpoints:
        if not os.path.exists(path):
            raise SystemExit(f"checkpoint not found: {path}")
    run(args.backend, checkpoints)


if __name__ == "__main__":
    main()
