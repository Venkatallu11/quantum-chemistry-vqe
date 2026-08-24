#!/usr/bin/env python3
"""
task40_certification_ablation_adversarial.py -- iteration 40, Tasks
40G+40H+40I (partial). Certification of the 0.011-0.018 kcal/mol result
(Task 39H/I/J, 3-for-3 real replication) via cheap, LOCAL checks on the
ALREADY-COLLECTED real data (zero new real submissions) -- deliberately
run BEFORE spending any more real shot budget, since these could
immediately falsify the result if something is wrong.

ABLATION MATRIX (40G): decompose the pipeline into its components on the
SAME real draw-0 data, to see whether the 0.01-scale result comes from
the genuine COMBINATION or from one component doing nearly everything:
  raw            -- no ancilla marginalized out, no correction at all
  QED only       -- postselected, no PEC correction (ratio=1)
  PEC only       -- not postselected, WITH the standard (unconditional)
                    analytic correction, GPi2=0
  GPi2 only      -- not postselected, correction WITH GPi2 (no QED)
  QED+PEC        -- postselected + conditioned correction, GPi2=0
  QED+PEC+GPi2   -- the full candidate (Task 39H's selected point)

ADVERSARIAL / NEGATIVE CONTROLS (40H): the estimator MUST fail loudly on
these, or the good result is suspect --
  shuffled ancilla bit   -- randomly reassign which shots are "kept",
                            same retained fraction, destroys any real
                            correlation between the ancilla and the
                            actual leaked state
  wrong-parity condition -- postselect on ODD weight instead of even
                            (the WRONG physical condition)
  wrong-sign GPi2         -- use -p_gpi2_selected instead of the
                            selected value
  shuffled labels         -- real, established adversarial-rejection
                            check (Task 36/37's own convention): shuffle
                            which (slot,label) pair each real value is
                            assigned to before fitting the frame

Run:
    PYTHONHASHSEED=0 python vqe/task40_certification_ablation_adversarial.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task30b_pec_application import build_full
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame, slot_vector
from task37c_extended_forward_model import analytic_A_and_B_5param
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from ionq_simulator_binding_curve import expectation_from_counts

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
GPI2_SELECTED = {"aria-1": 0.0006, "forte-1": 0.0004}  # Task 39H/I/J's consistently-selected candidates
NEW_CKPT = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission.partial.json")

_CACHE_UNC, _CACHE_COND = {}, {}


def energy_and_err(p, raw, K):
    alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
    E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                          exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
    return errs["err_vs_exact_kcal"]


def get_ab(name, p_gpi2, delta_zz, delta_gpi2, fixed_solutions, non_id_labels, conditioned):
    cache = _CACHE_COND if conditioned else _CACHE_UNC
    key = (name, round(p_gpi2, 6), round(delta_zz, 6), round(delta_gpi2, 6))
    if key not in cache:
        if conditioned:
            A, B, retA, retB = analytic_A_and_B_conditioned(fixed_solutions[name]["angles"], GATE_NAME,
                                                               ZZ_ASSUMED, GPI_REAL_MEAN, p_gpi2, delta_zz,
                                                               delta_gpi2, non_id_labels)
        else:
            A, B = analytic_A_and_B_5param(fixed_solutions[name]["angles"], GATE_NAME, ZZ_ASSUMED,
                                             GPI_REAL_MEAN, p_gpi2, delta_zz, delta_gpi2, non_id_labels)
        cache[key] = (A, B)
    return cache[key]


def corrected_energy(raw_by_slot_label, p_gpi2, kept, non_id_labels, fixed_solutions, diag, p, conditioned,
                      no_correction=False):
    """no_correction=True: TRUE raw, ratio identically 1 for every label --
    the genuine 'no correction at all' ablation baseline. Without this
    explicit path, p_gpi2=0 with the unconditional correction still
    applies the ZZ_ASSUMED ratio (Task 38C's own real, if tiny, ZZ
    correction) -- fine for a 'PEC-only-at-p_gpi2=0' ablation row, but
    NOT the same thing as genuinely no correction, so 'raw' needs this
    separate, explicit path to mean what it says."""
    raw_kept = {}
    for name in kept:
        raw_kept[name] = {}
        if no_correction:
            for l, m_raw in raw_by_slot_label[name].items():
                raw_kept[name][l] = max(-1.0, min(1.0, m_raw))
            continue
        A, B = get_ab(name, p_gpi2, 0.0, 0.0, fixed_solutions, non_id_labels, conditioned)
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    full = build_full(raw_kept, diag, K, non_id_labels)
    return energy_and_err(p, full, K)


def load_raw(kept, backend_name, mode, rng=None):
    """mode: 'raw_no_ancilla' (marginalize, no postselection), 'postselected'
    (real ancilla=0), 'shuffled_ancilla' (randomly reassign keep/discard,
    same retained count), 'wrong_parity' (postselect on ODD weight)."""
    with open(NEW_CKPT) as f:
        state = json.load(f)
    blended = {name: {} for name in kept}
    for name in kept:
        entry = state["done"][f"{backend_name}|{name}"]
        for group, counts in zip(entry["groups"], entry["counts"]):
            if mode == "raw_no_ancilla":
                filtered = {}
                for bs, c in counts.items():
                    filtered[bs[1:]] = filtered.get(bs[1:], 0) + c
            elif mode == "postselected":
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
            elif mode == "wrong_parity":
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "1"}
            elif mode == "shuffled_ancilla":
                keys = list(counts.keys())
                anc_bits = [k[0] for k in keys for _ in range(counts[k])]
                rng.shuffle(anc_bits)
                filtered = {}
                idx = 0
                for k in keys:
                    reg = k[1:]
                    for _ in range(counts[k]):
                        if anc_bits[idx] == "0":
                            filtered[reg] = filtered.get(reg, 0) + 1
                        idx += 1
            else:
                raise ValueError(mode)
            for l in group:
                blended[name][l] = expectation_from_counts(filtered, l)
    return blended


def main():
    print("\n" + "=" * 96)
    print("  task40_certification_ablation_adversarial.py -- ablation matrix + adversarial controls, real data")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}
    rng = np.random.default_rng(40)

    print("\n" + "=" * 60 + "  ABLATION MATRIX (40G)  " + "=" * 60)
    for backend_name in ["aria-1", "forte-1"]:
        p_gpi2_sel = GPI2_SELECTED[backend_name]
        raw_no_anc = load_raw(kept, backend_name, "raw_no_ancilla")
        ps = load_raw(kept, backend_name, "postselected")

        err_raw = corrected_energy(raw_no_anc, 0.0, kept, non_id_labels, fixed_solutions, diag, p,
                                    conditioned=False, no_correction=True)
        err_qed_only = corrected_energy(ps, 0.0, kept, non_id_labels, fixed_solutions, diag, p,
                                          conditioned=False, no_correction=True)
        err_pec_only = corrected_energy(raw_no_anc, 0.0, kept, non_id_labels, fixed_solutions, diag, p, conditioned=False)
        err_gpi2_only = corrected_energy(raw_no_anc, p_gpi2_sel, kept, non_id_labels, fixed_solutions, diag, p, conditioned=False)
        err_qed_pec = corrected_energy(ps, 0.0, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)
        err_full = corrected_energy(ps, p_gpi2_sel, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)

        print(f"\n  === {backend_name} (p_gpi2_selected={p_gpi2_sel}) ===")
        print(f"    raw (no ancilla, NO correction at all):  {err_raw:+.4f} kcal/mol")
        print(f"    QED only (postselected, NO correction):  {err_qed_only:+.4f} kcal/mol")
        print(f"    PEC only (no QED, ZZ-only correction):   {err_pec_only:+.4f} kcal/mol")
        print(f"    GPi2 only (no QED, +GPi2 correction):    {err_gpi2_only:+.4f} kcal/mol")
        print(f"    QED+PEC (postselect+conditioned, no GPi2): {err_qed_pec:+.4f} kcal/mol")
        print(f"    QED+PEC+GPi2 (FULL CANDIDATE):            {err_full:+.4f} kcal/mol")

    print("\n" + "=" * 55 + "  ADVERSARIAL / NEGATIVE CONTROLS (40H)  " + "=" * 55)
    for backend_name in ["aria-1", "forte-1"]:
        p_gpi2_sel = GPI2_SELECTED[backend_name]
        ps = load_raw(kept, backend_name, "postselected")
        wrong_parity = load_raw(kept, backend_name, "wrong_parity")
        shuffled_anc = load_raw(kept, backend_name, "shuffled_ancilla", rng=rng)

        err_full = corrected_energy(ps, p_gpi2_sel, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)
        err_wrong_parity = corrected_energy(wrong_parity, p_gpi2_sel, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)
        err_shuffled_anc = corrected_energy(shuffled_anc, p_gpi2_sel, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)
        err_wrong_sign_gpi2 = corrected_energy(ps, -p_gpi2_sel, kept, non_id_labels, fixed_solutions, diag, p, conditioned=True)

        print(f"\n  === {backend_name} ===")
        print(f"    REAL candidate (correct everything):        {err_full:+.4f} kcal/mol")
        print(f"    WRONG parity (postselect on ODD weight):    {err_wrong_parity:+.4f} kcal/mol   "
              f"{'FAILS LOUDLY as expected' if abs(err_wrong_parity) > 1.0 else 'DID NOT FAIL -- suspicious, investigate'}")
        print(f"    SHUFFLED ancilla (destroy real correlation): {err_shuffled_anc:+.4f} kcal/mol   "
              f"{'FAILS LOUDLY as expected' if abs(err_shuffled_anc) > 1.0 else 'DID NOT FAIL -- suspicious, investigate'}")
        print(f"    WRONG-SIGN GPi2 correction:                  {err_wrong_sign_gpi2:+.4f} kcal/mol   "
              f"{'FAILS LOUDLY as expected' if abs(err_wrong_sign_gpi2) > 1.0 else 'DID NOT FAIL -- suspicious, investigate'}")

    print("\n" + "=" * 50 + "  SHUFFLED-LABEL ADVERSARIAL FRAME FIT (40H, established convention)  " + "=" * 30)
    for backend_name in ["aria-1", "forte-1"]:
        p_gpi2_sel = GPI2_SELECTED[backend_name]
        ps = load_raw(kept, backend_name, "postselected")
        raw_kept = {}
        for name in kept:
            A, B = get_ab(name, p_gpi2_sel, 0.0, 0.0, fixed_solutions, non_id_labels, conditioned=True)
            raw_kept[name] = {}
            for l, m_raw in ps[name].items():
                ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
                raw_kept[name][l] = max(-1.0, min(1.0, m_raw * ratio))
        shuffled = {}
        for name in kept:
            labels_here = list(raw_kept[name].keys())
            vals_here = list(raw_kept[name].values())
            perm = rng.permutation(len(vals_here))
            shuffled[name] = {l: vals_here[perm[i]] for i, l in enumerate(labels_here)}
        _, _, chi2dof_real = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, raw_kept, weight_unit, rng)
        _, _, chi2dof_adv = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, shuffled, weight_unit, rng)
        ratio = chi2dof_adv / max(chi2dof_real, 1e-12)
        print(f"    {backend_name}: chi2/dof REAL={chi2dof_real:.5f}  SHUFFLED={chi2dof_adv:.5f}  "
              f"ratio={ratio:.1f}x  {'PASS (real data far more consistent than scrambled)' if ratio > 3 else 'FAIL -- suspicious'}")


if __name__ == "__main__":
    main()
