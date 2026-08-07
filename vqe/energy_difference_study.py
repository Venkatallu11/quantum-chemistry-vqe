#!/usr/bin/env python3
"""
energy_difference_study.py — does CDR's/raw's systematic BIAS cancel in
energy DIFFERENCES between neighboring geometries, the way chemistry
actually needs (reaction energies, binding curves, barrier heights), even
though it doesn't shrink in the absolute energy (iteration 6)?
============================================================================
THE INSIGHT (established in iteration 6, not re-derived): raw and CDR's
residual error is BIAS, not variance — flat across four orders of
magnitude of shots (CDR: 2.70-2.92 kcal/mol from 1e4-1e7 shots/setting;
raw: pinned ~104). A bias that is SIMILAR at two nearby geometries cancels
in their difference; chemical accuracy is defined on energy differences
(reaction energies, binding curves, barrier heights), not absolute
energies. This file tests whether that cancellation is real here.

WHY THIS MIGHT NOT WORK, STATED UP FRONT: vqe/difference_cancellation_results.json
(pre-dates this ledger's iterations) found fragment errors ADD across
DIFFERENT molecules with DIFFERENT circuits — no cancellation there. This
study is a much more favorable case for cancellation in principle (same
molecule, same FIXED 11-gate circuit STRUCTURE, only the target angles
differ across geometries) — but it is not guaranteed, and the honest
result is whatever the cancellation factor (absolute error / difference
error) actually comes out to, including if it's ~1.0 (no cancellation).

EFFICIENCY, JUSTIFIED NOT ASSUMED: verified directly that the Hamiltonian's
alpha-label SET (37 labels) is IDENTICAL across all 7 geometries tested
(d=0.8 through d=2.0) before relying on it. That means:
  - CDR training (random angles, geometry-independent by construction) and
    its resulting per-basis scales are reused across all 7 geometries at a
    given (shot_level, seed) -- computed once per (shot_level, seed), not
    once per (shot_level, seed, geometry).
  - PEC's Clifford-only calibration (gate noise, not target-Hamiltonian-
    dependent) is likewise reused across all 7 geometries at a given
    (shot_level, seed).
Only the per-geometry target-angle refit and the per-geometry exact/noisy
matrix measurement are actually geometry-specific -- everything else is
shared, matching how a real experimentalist would scan a binding curve
(calibrate once, reuse for every point), not re-deriving the noise model
from scratch at each geometry.

K=6 EXACTNESS: re-verified (not assumed) at every geometry that the
Schmidt rank stays <=6 (the weight-2 x weight-2 physical-subspace bound is
a particle-counting argument, geometry-independent in principle, but
checked directly here rather than taken on faith).

Run:
    python vqe/energy_difference_study.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
import rank6_symmetry_vd as r
import loop_pec as pec
import loop_pauli_lindblad_pec as pl
import shot_noise_study as sns
from fixed_ansatz import build_ansatz, P2_PER_GATE, P1_PER_GATE

K = 6
ATOMS, NELEC = [0, 1, 2, 3], 4
HARTREE_TO_KCAL_MOL = r.HARTREE_TO_KCAL_MOL
BASIS_GATES = r.BASIS_GATES
GEOMETRIES = [0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]
D_REF = 1.0
SHOT_LEVELS = [1_000, 10_000, 100_000, 1_000_000, 10_000_000]
N_SEEDS = 8
N_TRAIN_PER_SLOT = r.N_TRAIN_PER_SLOT
GAMMA_TOTAL = None  # loaded from iteration 5's results in main()

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "energy_difference_study_results.json")


# ---------------------------------------------------------------------------
# Setup at an arbitrary bond distance
# ---------------------------------------------------------------------------

def setup_at_distance(d, K=6):
    qop_bare, qop_pen, enuc = effrag.build_fragment_qop(ATOMS, NELEC, d=d)
    n_qubits = qop_bare.num_qubits
    e_elec, psi = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    psi_real, residual = effrag.real_gauge(psi)
    lambdas, u_vecs, v_vecs = effrag.schmidt_decompose_real(psi_real, n_qubits)
    terms = effrag.decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    assert set(alpha_labels) == set(beta_labels)

    max_tail = float(np.max(np.abs(lambdas[K:])))
    schmidt_rank_le_K = bool(max_tail < 1e-9)

    u_top, v_top = u_vecs[:K], v_vecs[:K]
    signs = np.array([1.0 if np.dot(v_top[n], u_top[n]) >= 0 else -1.0 for n in range(K)])
    sign_residual = max(float(np.max(np.abs(v_top[n] - signs[n] * u_top[n]))) for n in range(K))
    assert sign_residual < 1e-8, f"d={d}: beta=sign*alpha residual {sign_residual:.3e}"

    identity_label = "I" * (n_qubits // 2)
    identity_coeff = 0.0
    for label, coeff in qop_bare.to_list():
        if label == "I" * n_qubits:
            identity_coeff = float(coeff.real)
            break
    e_mixed = identity_coeff + enuc

    alpha_cache, beta_cache = effrag.precompute_exact_matrices(terms, u_top, v_top)
    noiseless_numpy = effrag.ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)
    truncation_floor_kcal = abs(exact_energy - noiseless_numpy) * HARTREE_TO_KCAL_MOL

    targets = {}
    for n in range(K):
        targets[f"u_{n}"] = u_top[n]
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        targets[f"(u{n}+u{m})"] = (u_top[n] + u_top[m]) / np.sqrt(2)
        targets[f"(u{n}-u{m})"] = (u_top[n] - u_top[m]) / np.sqrt(2)
    assert set(targets.keys()) == set(r.slot_names(K))

    return {
        "d": d, "n_qubits": n_qubits, "terms": terms, "alpha_labels": alpha_labels,
        "identity_label": identity_label, "signs": signs, "lambdas": lambdas, "u_vecs": u_top,
        "enuc": enuc, "e_mixed": e_mixed, "noiseless_numpy": noiseless_numpy, "exact_energy": exact_energy,
        "truncation_floor_kcal": truncation_floor_kcal, "targets": targets,
        "schmidt_rank_le_K": schmidt_rank_le_K, "max_schmidt_tail": max_tail, "K": K,
    }


def fit_and_verify(p):
    solutions, n_ok, worst = r.fit_all_targets(p["targets"])
    p["solutions"] = solutions
    counts = r.verify_constant_gate_count(solutions)
    return n_ok, worst, counts


def energy_errors_at(E, p):
    err_vs_exact = abs(E - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    err_vs_noiseless = abs(E - p["noiseless_numpy"]) * HARTREE_TO_KCAL_MOL
    return err_vs_exact, err_vs_noiseless


# ---------------------------------------------------------------------------
# Per-geometry exact-value caching (raw noisy matrices; PEC needs the
# per-shot-level learned channel applied per geometry, done in the sweep)
# ---------------------------------------------------------------------------

def cache_exact_raw_at(p, non_id_labels, noise_model):
    raw = {}
    for name, sol in p["solutions"].items():
        dm = r.density_matrix_4q(sol["angles"], noise_model)
        vals, _ = r.measure_labels(sol["angles"], non_id_labels, dm)
        raw[name] = vals
    return raw


PARTIAL_PATH_TMPL = os.path.join(os.path.dirname(__file__), "energy_difference_study_partial_{}.json")


def setup_common():
    """Everything that does NOT depend on shot level: per-geometry setup +
    angle fits, CDR training cache, PEC calibration cache, exact raw
    matrices. Computed once, reused by whichever shot level is being run
    in this process (the Bash tool caps a single command at 10 minutes,
    even backgrounded, so the full 5-shot-level x 8-seed x 7-geometry
    sweep is run as one process PER SHOT LEVEL, checkpointed to disk, and
    combined in --assemble -- this setup work is small enough (~1-2 min)
    to just redo per chunk rather than build cross-process caching for it)."""
    print("\n" + "=" * 96)
    print("  energy_difference_study.py -- does systematic bias cancel between geometries?")
    print("=" * 96)

    noise_model = r.build_noise_model()

    global GAMMA_TOTAL
    with open(os.path.join(os.path.dirname(__file__), "loop_pauli_lindblad_pec_results.json")) as f:
        lindblad_results = json.load(f)
    GAMMA_TOTAL = lindblad_results["cost_accounting"]["gamma_total"]
    print(f"\n  loaded PEC sampling-overhead factor gamma_total={GAMMA_TOTAL:.4f} from iteration 5")

    # ==== setup + verification at every geometry ====
    print(f"\n  -- setting up and fitting the fixed 11-gate ansatz at {len(GEOMETRIES)} geometries --")
    geoms = {}
    for d in GEOMETRIES:
        t0 = time.time()
        p_d = setup_at_distance(d, K)
        n_ok, worst, counts = fit_and_verify(p_d)
        geoms[d] = p_d
        print(f"    d={d} A: Schmidt rank<=6: {p_d['schmidt_rank_le_K']} (tail={p_d['max_schmidt_tail']:.2e}), "
              f"fit {n_ok}/36 converged (worst={worst:.2e}), gate count={counts}, "
              f"exact_energy={p_d['exact_energy']:.6f} Ha, {time.time()-t0:.1f}s")
        assert p_d["schmidt_rank_le_K"], f"d={d}: K=6 is NOT exact here -- do not proceed on that assumption"
        assert counts == {11}, f"d={d}: gate count not fixed: {counts}"

    p_ref = geoms[D_REF]
    non_id_labels = [l for l in p_ref["alpha_labels"] if l != p_ref["identity_label"]]
    for d in GEOMETRIES:
        assert set(geoms[d]["alpha_labels"]) == set(p_ref["alpha_labels"]), \
            f"alpha-label set differs at d={d} -- CDR-scale/calibration reuse across geometries is invalid"
    print(f"\n  verified: alpha-label set identical across all {len(GEOMETRIES)} geometries "
          "(justifies reusing CDR training and PEC calibration across geometries, not assumed)")

    # ==== cache exact raw matrices per geometry (shot-count independent) ====
    exact_raw_by_d = {}
    for d in GEOMETRIES:
        exact_raw_by_d[d] = cache_exact_raw_at(geoms[d], non_id_labels, noise_model)

    # ==== cache CDR training (geometry-independent) and PEC calibration (geometry-independent) ====
    print(f"\n  caching CDR training data ({N_SEEDS} seeds, geometry-independent) and PEC calibration "
          "(geometry-independent) -- reused across all 7 geometries below")
    t0 = time.time()
    exact_training_by_seed = {seed: sns.cache_exact_training(p_ref, non_id_labels, noise_model,
                                                              N_TRAIN_PER_SLOT, seed)
                               for seed in range(N_SEEDS)}
    print(f"    CDR training cached, {time.time()-t0:.1f}s")
    t0 = time.time()
    cx_cal_cache, u3_cal_cache = sns.cache_calibration_exact(noise_model)
    print(f"    PEC calibration circuits cached, {time.time()-t0:.1f}s")

    return {
        "geoms": geoms, "non_id_labels": non_id_labels, "noise_model": noise_model,
        "exact_raw_by_d": exact_raw_by_d, "exact_training_by_seed": exact_training_by_seed,
        "cx_cal_cache": cx_cal_cache, "u3_cal_cache": u3_cal_cache,
    }


METHODS = ("raw", "cdr", "pec_honest")


def run_shot_level(n_shots, common):
    """Full 8-seed x 7-geometry sweep at ONE shot level -- kept as its own
    process invocation (see PARTIAL_PATH_TMPL) because the full 5-level
    sweep does not fit in a single Bash-tool command's 10-minute cap, even
    backgrounded."""
    geoms, non_id_labels, noise_model = common["geoms"], common["non_id_labels"], common["noise_model"]
    exact_raw_by_d, exact_training_by_seed = common["exact_raw_by_d"], common["exact_training_by_seed"]
    cx_cal_cache, u3_cal_cache = common["cx_cal_cache"], common["u3_cal_cache"]

    E_data = {m: {d: [] for d in GEOMETRIES} for m in METHODS}
    p2_rel_errs = []

    t_shots0 = time.time()
    if True:
        for seed in range(N_SEEDS):
            # --- CDR scales for this (shot_level, seed), shared across geometries ---
            rng_cdr = np.random.default_rng(seed * 1_000_003 + n_shots)
            training = []
            for row in exact_training_by_seed[seed]:
                for l in non_id_labels:
                    noisy_shot = sns.shot_sample(row["noisy_exact"][l], n_shots, rng_cdr)
                    training.append({"slot": row["slot"], "label": l, "exact": row["exact"][l], "noisy": noisy_shot})
            scales = r.fit_all_scales(training, non_id_labels, K)

            # --- PEC calibration for this (shot_level, seed), shared across geometries ---
            rng_calib = np.random.default_rng(seed * 2_000_003 + n_shots)
            p2_learned, p1_learned = sns.shot_noisy_calibrate(cx_cal_cache, u3_cal_cache, n_shots, rng_calib)
            p2_rel_errs.append(abs(p2_learned - P2_PER_GATE) / P2_PER_GATE)

            # --- per-geometry measurement, independent shot draws per geometry ---
            for d_idx, d in enumerate(GEOMETRIES):
                p_d = geoms[d]
                rng_raw = np.random.default_rng(seed * 3_000_003 + n_shots * 100 + d_idx)
                raw_d = {name: {l: sns.shot_sample(v, n_shots, rng_raw) for l, v in vals.items()}
                         for name, vals in exact_raw_by_d[d].items()}
                mats_raw = r.combine_matrices(raw_d, p_d["alpha_labels"], p_d["identity_label"], K)
                E_raw, _ = r.energy_from_alpha_matrices(mats_raw, p_d, K)
                E_data["raw"][d].append(E_raw)

                mats_cdr = r.combine_matrices(raw_d, p_d["alpha_labels"], p_d["identity_label"], K,
                                               per_label_scale=scales["per_basis_scale"])
                E_cdr, _ = r.energy_from_alpha_matrices(mats_cdr, p_d, K)
                E_data["cdr"][d].append(E_cdr)

                rng_pec = np.random.default_rng(seed * 4_000_003 + n_shots * 100 + d_idx)
                raw_pec_d = {}
                for name, sol in p_d["solutions"].items():
                    exact_pec_vals = sns.pec_measure_learned(sol["angles"], non_id_labels, P2_PER_GATE, P1_PER_GATE,
                                                              p2_learned, p1_learned)
                    raw_pec_d[name] = {l: sns.pec_shot_sample(v, GAMMA_TOTAL, n_shots, rng_pec)
                                        for l, v in exact_pec_vals.items()}
                mats_pec = r.combine_matrices(raw_pec_d, p_d["alpha_labels"], p_d["identity_label"], K)
                E_pec, _ = r.energy_from_alpha_matrices(mats_pec, p_d, K)
                E_data["pec_honest"][d].append(E_pec)

    print(f"    n_shots={n_shots:>12,}: done, {time.time()-t_shots0:.1f}s, "
          f"mean calib p2 rel.err={np.mean(p2_rel_errs)*100:.3f}%")

    partial = {
        "n_shots": n_shots, "n_seeds": N_SEEDS, "geometries": GEOMETRIES,
        "E_data": {m: {str(d): E_data[m][d] for d in GEOMETRIES} for m in METHODS},
        "mean_calib_p2_rel_err": float(np.mean(p2_rel_errs)),
    }
    with open(PARTIAL_PATH_TMPL.format(n_shots), "w") as f:
        json.dump(partial, f, indent=2)
    print(f"    checkpoint saved -> {PARTIAL_PATH_TMPL.format(n_shots)}")
    return partial


def assemble(geoms):
    """Load every available per-shot-level checkpoint and do the full
    analysis: absolute vs difference error, cancellation factor, binding
    curve."""
    print("\n" + "=" * 96)
    print("  energy_difference_study.py --assemble")
    print("=" * 96)

    available = []
    for n_shots in SHOT_LEVELS:
        path = PARTIAL_PATH_TMPL.format(n_shots)
        if os.path.exists(path):
            with open(path) as f:
                available.append(json.load(f))
    found_shots = [p["n_shots"] for p in available]
    print(f"\n  found checkpoints for shot levels: {found_shots} (expected {SHOT_LEVELS})")
    E_data = {m: {d: {} for d in GEOMETRIES} for m in METHODS}
    mean_calib_err_by_shots = {}
    for part in available:
        n_shots = part["n_shots"]
        mean_calib_err_by_shots[n_shots] = part["mean_calib_p2_rel_err"]
        for m in METHODS:
            for d in GEOMETRIES:
                E_data[m][d][n_shots] = part["E_data"][m][str(d)]

    # ==== analysis: absolute error, difference error, cancellation factor ====
    print(f"\n  -- ABSOLUTE vs DIFFERENCE error, d_ref={D_REF} A --")
    exact_by_d = {d: geoms[d]["exact_energy"] for d in GEOMETRIES}
    noiseless_by_d = {d: geoms[d]["noiseless_numpy"] for d in GEOMETRIES}
    diff_exact_by_d = {d: exact_by_d[d] - exact_by_d[D_REF] for d in GEOMETRIES}

    per_geometry = {m: {} for m in METHODS}
    cancellation_summary = {m: {} for m in METHODS}

    for method in METHODS:
        print(f"\n  == {method} ==")
        for n_shots in found_shots:
            abs_err_all_d = []
            diff_err_all_d = []
            per_geom_row = {}
            for d in GEOMETRIES:
                E_vals = np.array(E_data[method][d][n_shots])
                abs_errs = np.abs(E_vals - exact_by_d[d]) * HARTREE_TO_KCAL_MOL
                noiseless_errs = np.abs(E_vals - noiseless_by_d[d]) * HARTREE_TO_KCAL_MOL
                per_geom_row[d] = {
                    "abs_err_vs_exact": {"mean": float(np.mean(abs_errs)), "std": float(np.std(abs_errs))},
                    "abs_err_vs_noiseless": {"mean": float(np.mean(noiseless_errs)), "std": float(np.std(noiseless_errs))},
                }
                if d != D_REF:
                    E_ref_vals = np.array(E_data[method][D_REF][n_shots])
                    diff_measured = E_vals - E_ref_vals
                    diff_errs = np.abs(diff_measured - diff_exact_by_d[d]) * HARTREE_TO_KCAL_MOL
                    per_geom_row[d]["diff_err_vs_exact"] = {"mean": float(np.mean(diff_errs)), "std": float(np.std(diff_errs))}
                    abs_err_all_d.append(np.mean(abs_errs))
                    diff_err_all_d.append(np.mean(diff_errs))
            per_geometry[method][n_shots] = per_geom_row

            mean_abs = float(np.mean(abs_err_all_d))
            mean_diff = float(np.mean(diff_err_all_d))
            cancellation_factor = mean_abs / mean_diff if mean_diff > 1e-9 else float("inf")
            cancellation_summary[method][n_shots] = {
                "mean_abs_error_kcal": mean_abs, "mean_diff_error_kcal": mean_diff,
                "cancellation_factor": cancellation_factor,
            }
            print(f"    n_shots={n_shots:>10,}: mean|abs err|={mean_abs:>8.3f} kcal/mol, "
                  f"mean|diff err|={mean_diff:>8.3f} kcal/mol, cancellation factor={cancellation_factor:.2f}x "
                  f"({'REAL CANCELLATION' if cancellation_factor > 1.5 else 'little/no cancellation'})")

    # ==== binding curve: equilibrium bond length + well depth, exact vs noisy ====
    # NOTE, a real methodological fix: a single quadratic across the FULL
    # 0.8-2.0 A range is badly distorted by the anharmonic dissociation tail
    # (H4's energy flattens at large d, it does not keep curving like a
    # parabola) -- an all-7-point fit gave a nonsensical d_eq around -4.4 A
    # on the first run here. Fixed by fitting the parabola only to a LOCAL
    # WINDOW of points actually near the minimum (the standard way to
    # extract equilibrium geometry/well depth from a sampled curve), chosen
    # from the EXACT curve's own lowest point, not assumed in advance.
    d_dissoc = max(GEOMETRIES)
    ds_all = np.array(GEOMETRIES)
    exact_es_all = np.array([exact_by_d[d] for d in GEOMETRIES])
    min_idx = int(np.argmin(exact_es_all))
    lo = max(0, min_idx - 2)
    hi = min(len(GEOMETRIES), min_idx + 3)
    local_ds = [GEOMETRIES[i] for i in range(lo, hi)]
    print(f"  local fit window around the exact minimum (d={GEOMETRIES[min_idx]} A): {local_ds}")

    ds = np.array(local_ds)
    exact_es = np.array([exact_by_d[d] for d in local_ds])
    coeffs_exact = np.polyfit(ds, exact_es, 2)
    d_eq_exact = -coeffs_exact[1] / (2 * coeffs_exact[0])
    E_eq_exact = np.polyval(coeffs_exact, d_eq_exact)
    well_depth_exact = (exact_by_d[d_dissoc] - E_eq_exact) * HARTREE_TO_KCAL_MOL
    print(f"    exact: d_eq={d_eq_exact:.4f} A, E(d_eq)={E_eq_exact:.6f} Ha, "
          f"well depth (vs d={d_dissoc})={well_depth_exact:.3f} kcal/mol")

    binding_curve_results = {"exact": {"d_eq": float(d_eq_exact), "E_eq_ha": float(E_eq_exact),
                                        "well_depth_kcal": float(well_depth_exact), "fit_window": local_ds}}
    representative_shots = 100_000 if 100_000 in found_shots else found_shots[len(found_shots) // 2]
    for method in METHODS:
        noisy_means_local = np.array([np.mean(E_data[method][d][representative_shots]) for d in local_ds])
        coeffs_noisy = np.polyfit(ds, noisy_means_local, 2)
        d_eq_noisy = -coeffs_noisy[1] / (2 * coeffs_noisy[0])
        E_eq_noisy = np.polyval(coeffs_noisy, d_eq_noisy)
        noisy_dissoc = float(np.mean(E_data[method][d_dissoc][representative_shots]))
        well_depth_noisy = (noisy_dissoc - E_eq_noisy) * HARTREE_TO_KCAL_MOL
        d_eq_err = abs(d_eq_noisy - d_eq_exact)
        well_depth_err = abs(well_depth_noisy - well_depth_exact)
        print(f"    {method} (at {representative_shots:,} shots/setting): d_eq={d_eq_noisy:.4f} A "
              f"(err={d_eq_err:.4f} A), well depth={well_depth_noisy:.3f} kcal/mol (err={well_depth_err:.3f} kcal/mol)")
        binding_curve_results[method] = {"d_eq": float(d_eq_noisy), "E_eq_ha": float(E_eq_noisy),
                                          "well_depth_kcal": float(well_depth_noisy),
                                          "d_eq_error_A": float(d_eq_err), "well_depth_error_kcal": float(well_depth_err),
                                          "representative_shots": representative_shots, "fit_window": local_ds}

    results = {
        "geometries": GEOMETRIES, "d_ref": D_REF, "shot_levels_intended": SHOT_LEVELS,
        "shot_levels_found": found_shots, "n_seeds": N_SEEDS,
        "mean_calib_p2_rel_err_by_shots": mean_calib_err_by_shots,
        "exact_energy_by_d": exact_by_d, "noiseless_energy_by_d": noiseless_by_d,
        "schmidt_rank_check": {str(d): {"le_K": geoms[d]["schmidt_rank_le_K"], "max_tail": geoms[d]["max_schmidt_tail"]}
                                for d in GEOMETRIES},
        "per_geometry": {m: {str(k): v for k, v in per_geometry[m].items()} for m in METHODS},
        "cancellation_summary": {m: {str(k): v for k, v in cancellation_summary[m].items()} for m in METHODS},
        "binding_curve": binding_curve_results,
        "chemical_accuracy_kcal": 1.0,
    }
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--shots", type=int, help="run the full 8-seed x 7-geometry sweep for this one shot level")
    parser.add_argument("--assemble", action="store_true", help="combine all available checkpoints and analyze")
    args = parser.parse_args()

    if args.assemble:
        # geometry setup only (fast, no CDR training/calibration needed for analysis)
        common = setup_common()
        return assemble(common["geoms"])
    if args.shots:
        common = setup_common()
        return run_shot_level(args.shots, common)
    parser.error("pass --shots N (one of 1000/10000/100000/1000000/10000000) or --assemble")


if __name__ == "__main__":
    main()
