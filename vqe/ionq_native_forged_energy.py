#!/usr/bin/env python3
"""
ionq_native_forged_energy.py — the FULL H4 entanglement-forging energy
(K=5, 185 terms, 37 alpha labels), measured on real IonQ circuits built
entirely from NATIVE gates (GPi/GPi2/MS or GPi/GPi2/ZZ), at folds 1/3/5.
============================================================================
ionq_fold_check.py --native proved gate-folding survives on a simple
2-qubit probe when submitted natively (f decays cleanly 0.986->0.339 on
aria-1, 0.985->0.312 on forte-1, folds 1..81). But that probe measured one
Schmidt vector's worth of signal through a single Z-type observable. The
real forged energy sums 185 Pauli terms across 5 Schmidt indices and 10
pairs, each through its own basis rotation -- more circuits, more
opportunities for the per-circuit signal to fall through differently.
This measures the real thing: does the clean per-gate decay seen on the
probe actually show up in the quantity that matters (the reassembled
energy), and does ZNE recover anything useful from it.

ONE GATESET PER CIRCUIT, NEVER MIXED: every circuit here is built as
state-prep (native, from native_stateprep.py) -> FOLDED (native, before
anything else touches it) -> basis-change (built abstractly, then
translated to native and appended -- never re-transpiled together with
the already-folded state-prep, so the fold is never at risk of being
optimized away by a second pass). Reuses, unchanged:
  - native_stateprep.build_real_state_prep_circuit / to_native
    (hand-derived real-only tree, translated via IonQ's own tested CX
    equivalence library -- 14 two-qubit gates, verified on real IonQ)
  - ionq_fold_check.fold_native_2q (explicit, numerically-verified
    inverse: phase-shift for ms, X-conjugation sandwich for zz -- NOT
    Gate.inverse(), which returns a mis-parametrized, unusable gate)
  - ef_fragment.group_labels_qubit_wise / combined_basis_label
    (qubit-wise commuting grouping, 37 labels -> 13 circuits/Schmidt
    index -- the SAFE grouping this repo chose over a hand-derived
    Clifford-diagonalization synthesizer)
  - the real-gauge 2-phase-circuit shortcut and beta-register sign reuse
    (verified to hold for fragment A / H4 in ionq_tailoring.py)

Fragment: H4 (atoms [0,1,2,3], d=1.0 Ang) -- exactly fragment A from
ionq_tailoring.py, same physical system, same K=5.

Run (ideal MUST be run and confirmed before any noisy config -- refuses
otherwise, same gate as ionq_tailoring.py):
    python vqe/ionq_native_forged_energy.py --noise-model ideal --fold 1
    python vqe/ionq_native_forged_energy.py --noise-model aria-1 --fold 1
    python vqe/ionq_native_forged_energy.py --noise-model aria-1 --fold 3
    python vqe/ionq_native_forged_energy.py --noise-model aria-1 --fold 5
    python vqe/ionq_native_forged_energy.py --noise-model forte-1 --fold 1
    python vqe/ionq_native_forged_energy.py --noise-model forte-1 --fold 3
    python vqe/ionq_native_forged_energy.py --noise-model forte-1 --fold 5
    python vqe/ionq_native_forged_energy.py --assemble
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
from native_stateprep import build_real_state_prep_circuit, to_native
from ionq_fold_check import fold_native_2q
from ionq_backend import connect_provider, get_native_simulator
from ionq_run import pauli_expectation

HARTREE_TO_KCAL_MOL = 627.5094740631
ATOMS = [0, 1, 2, 3]
NELEC = 4
K = 5
GATE_BY_MODEL = {"ideal": "ms", "aria-1": "ms", "forte-1": "zz"}
RAW_BASELINE_KCAL = {"aria-1": 125.07, "forte-1": 134.62}
CLASSICAL_FLOOR_KCAL = 0.5655
CHEMICAL_ACCURACY_KCAL = 1.0

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "native_forged_zne_results.json")


# ---------------------------------------------------------------------------
# Setup (local, no network) -- identical physics to ionq_tailoring.py fragment A
# ---------------------------------------------------------------------------

def setup():
    qop_bare, qop_pen, enuc = effrag.build_fragment_qop(ATOMS, NELEC, d=1.0)
    n_qubits = qop_bare.num_qubits
    e_elec, psi = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    psi_real, residual = effrag.real_gauge(psi)
    lambdas, u_vecs, v_vecs = effrag.schmidt_decompose_real(psi_real, n_qubits)
    terms = effrag.decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    assert set(alpha_labels) == set(beta_labels)

    u_top, v_top = u_vecs[:K], v_vecs[:K]
    signs = np.array([1.0 if np.dot(v_top[n], u_top[n]) >= 0 else -1.0 for n in range(K)])
    sign_residual = max(float(np.max(np.abs(v_top[n] - signs[n] * u_top[n]))) for n in range(K))
    assert sign_residual < 1e-8, (
        f"beta=sign*alpha residual {sign_residual:.3e} -- shortcut doesn't hold, "
        "would need independent beta measurement (not implemented here)"
    )

    identity_label = "I" * n_qubits
    identity_coeff = 0.0
    for label, coeff in qop_bare.to_list():
        if label == identity_label:
            identity_coeff = float(coeff.real)
            break
    e_mixed = identity_coeff + enuc

    alpha_cache, beta_cache = effrag.precompute_exact_matrices(terms, u_top, v_top)
    noiseless_numpy = effrag.ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)

    groups = effrag.group_labels_qubit_wise(alpha_labels)

    return {
        "enuc": enuc, "exact_energy": exact_energy, "real_gauge_residual": residual,
        "lambdas": lambdas, "u_vecs": u_top, "terms": terms, "alpha_labels": alpha_labels,
        "signs": signs, "e_mixed": e_mixed, "noiseless_numpy": noiseless_numpy,
        "groups": groups,
    }


# ---------------------------------------------------------------------------
# Native circuit construction: state prep -> fold -> basis change (native, appended)
# ---------------------------------------------------------------------------

def native_basis_change(label, gate_name):
    """Abstract h/sdg circuit, translated to native via IonQ's own tested
    equivalence library (same safe pattern as native_stateprep.py's
    CX->native translation) -- verified (scratchpad, development) to use
    ZERO 2-qubit gates and match the abstract version to <1e-15."""
    from qiskit.circuit import QuantumCircuit
    n = len(label)
    qc = QuantumCircuit(n)
    for i, ch in enumerate(label):
        qubit = n - 1 - i
        if ch == "X":
            qc.h(qubit)
        elif ch == "Y":
            qc.sdg(qubit)
            qc.h(qubit)
    return to_native(qc, gate_name)


def native_measurement_circuit(folded_stateprep, label, gate_name):
    """Compose the (already folded) state-prep circuit with a SEPARATELY
    native-translated basis-change circuit -- never re-transpiled
    together, so the fold survives untouched."""
    qc = folded_stateprep.compose(native_basis_change(label, gate_name))
    qc.measure_all()
    return qc


def submit_batch(circuits, backend, noise_model, shots=100_000):
    job = backend.run(circuits, noise_model=noise_model, shots=shots)
    result = job.result()
    probs = result.get_probabilities()
    if not isinstance(probs, list):
        probs = [probs]
    return [dict(p) for p in probs]


def measure_alpha_matrices_native(u_vecs, groups, alpha_labels, backend, noise_model, gate_name, fold, K):
    matrices = {label: np.zeros((K, K), dtype=complex) for label in alpha_labels}
    n_circuits = 0

    # --- diagonal: one job per Schmidt index, one circuit per group ---
    for n in range(K):
        base = to_native(build_real_state_prep_circuit(u_vecs[n]), gate_name)
        folded = fold_native_2q(base, fold, gate_name)
        circuits, metas = [], []
        for group in groups:
            combined = effrag.combined_basis_label(group)
            circuits.append(native_measurement_circuit(folded, combined, gate_name))
            metas.append(group)
        probs_list = submit_batch(circuits, backend, noise_model)
        n_circuits += len(circuits)
        for group, probs in zip(metas, probs_list):
            for label in group:
                matrices[label][n, n] = pauli_expectation(probs, label)

    # --- cross: one job per pair, k=0 and k=2 phase states, one circuit per group each ---
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        circuits, metas = [], []
        for k, sign in ((0, 1.0), (2, -1.0)):
            vec_k = (u_vecs[n] + sign * u_vecs[m]) / np.sqrt(2)
            base = to_native(build_real_state_prep_circuit(vec_k), gate_name)
            folded = fold_native_2q(base, fold, gate_name)
            for group in groups:
                combined = effrag.combined_basis_label(group)
                circuits.append(native_measurement_circuit(folded, combined, gate_name))
                metas.append((k, group))
        probs_list = submit_batch(circuits, backend, noise_model)
        n_circuits += len(circuits)
        e_by_label_k = {label: {} for label in alpha_labels}
        for (k, group), probs in zip(metas, probs_list):
            for label in group:
                e_by_label_k[label][k] = pauli_expectation(probs, label)
        for label in alpha_labels:
            e0, e2 = e_by_label_k[label][0], e_by_label_k[label][2]
            re = (e0 - e2) / 2
            matrices[label][n, m] = re
            matrices[label][m, n] = re

    return matrices, n_circuits


def derive_beta_matrices(alpha_matrices, signs, K):
    S = np.diag(signs)
    return {label: S @ mat @ S for label, mat in alpha_matrices.items()}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def load_results():
    if os.path.exists(RESULTS_PATH):
        with open(RESULTS_PATH) as f:
            return json.load(f)
    return {}


def save_results(results):
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


def run_one(noise_model, fold):
    print("\n" + "=" * 70)
    print(f"  Native forged H4 energy -- noise_model={noise_model}, fold={fold}")
    print("=" * 70)

    if noise_model != "ideal":
        results = load_results()
        ideal = results.get("configs", {}).get("ideal_fold1")
        if not ideal:
            print("\n  REFUSING: run --noise-model ideal --fold 1 first and confirm it "
                  "matches the classical/numpy reference before any noisy run.")
            sys.exit(1)

    gate_name = GATE_BY_MODEL[noise_model]
    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"\n  Connected. Backend: {backend.name} (gateset=native, gate={gate_name})")

    p = setup()
    print(f"  real-gauge residual: {p['real_gauge_residual']:.3e}")
    print(f"  {len(p['alpha_labels'])} alpha labels -> {len(p['groups'])} qubit-wise groups")
    print(f"  exact fragment energy = {p['exact_energy']:.6f} Ha")
    print(f"  exact numpy EF energy (K={K}) = {p['noiseless_numpy']:.6f} Ha")

    alpha_mats, n_circuits = measure_alpha_matrices_native(
        p["u_vecs"], p["groups"], p["alpha_labels"], backend, noise_model, gate_name, fold, K,
    )
    beta_mats = derive_beta_matrices(alpha_mats, p["signs"], K)
    E = effrag.ef_energy_from_noisy_matrices(p["terms"], p["lambdas"], alpha_mats, beta_mats, p["enuc"], K)

    err_ha = abs(E - p["exact_energy"])
    err_kcal = err_ha * HARTREE_TO_KCAL_MOL
    f_signal = (E - p["e_mixed"]) / (p["noiseless_numpy"] - p["e_mixed"])

    print(f"\n  E = {E:.6f} Ha, {n_circuits} circuits, error = {err_kcal:.3f} kcal/mol, f = {f_signal:.6f}")

    results = load_results()
    results.setdefault("configs", {})
    key = f"{noise_model}_fold{fold}"
    results["configs"][key] = {
        "noise_model": noise_model, "gate": gate_name, "fold": fold, "K": K,
        "n_circuits": n_circuits,
        "exact_energy_ha": round(p["exact_energy"], 6),
        "exact_numpy_ef_k_ha": round(p["noiseless_numpy"], 6),
        "e_mixed_ha": round(p["e_mixed"], 6),
        "measured_energy_ha": round(E, 6),
        "error_kcal": round(err_kcal, 4),
        "f": round(f_signal, 6),
    }
    save_results(results)
    return results


def assemble():
    results = load_results()
    configs = results.get("configs", {})
    print("\n  Configs completed:", sorted(configs.keys()))

    ideal_key = "ideal_fold1"
    if ideal_key not in configs:
        print("  Missing ideal_fold1 -- run that first.")
        return results
    ideal = configs[ideal_key]
    print(f"\n  ideal (fold=1): E={ideal['measured_energy_ha']} Ha, "
          f"error vs exact numpy EF-K={K} = "
          f"{abs(ideal['measured_energy_ha'] - ideal['exact_numpy_ef_k_ha']) * HARTREE_TO_KCAL_MOL:.3f} kcal/mol")

    zne_report = {}
    for model in ("aria-1", "forte-1"):
        rows = []
        for fold in (1, 3, 5):
            key = f"{model}_fold{fold}"
            if key not in configs:
                print(f"  Missing {key}")
                continue
            rows.append(configs[key])
        if len(rows) < 3:
            print(f"  {model}: incomplete ({len(rows)}/3 folds)")
            continue

        folds = [r["fold"] for r in rows]
        energies = [r["measured_energy_ha"] for r in rows]
        f_vals = [r["f"] for r in rows]
        exact = rows[0]["exact_energy_ha"]

        eff_rate = [1 - f_vals[i] ** (1 / folds[i]) for i in range(len(folds))]
        rate_spread = max(eff_rate) - min(eff_rate)
        rate_consistent = rate_spread < 0.01  # generous but explicit threshold

        print(f"\n  -- {model} --")
        print(f"  {'fold':>6} {'E (Ha)':>12} {'error (kcal/mol)':>18} {'f':>10} {'1-f^(1/L)':>12}")
        for i, r in enumerate(rows):
            print(f"  {r['fold']:>6} {r['measured_energy_ha']:>12.6f} "
                  f"{r['error_kcal']:>18.3f} {r['f']:>10.6f} {eff_rate[i]:>12.6f}")
        print(f"  per-fold effective error rate spread: {rate_spread:.6f} "
              f"({'roughly constant -- ' if rate_consistent else 'NOT constant -- '}"
              f"{'extrapolation trustworthy by this check' if rate_consistent else 'STOP: something differs from the probe, extrapolation is NOT trustworthy'})")

        no_mit_kcal = abs(energies[0] - exact) * HARTREE_TO_KCAL_MOL
        lin_fit = np.polyfit(folds, energies, 1)
        E_lin0 = float(np.polyval(lin_fit, 0))
        lin_kcal = abs(E_lin0 - exact) * HARTREE_TO_KCAL_MOL

        quad_fit = np.polyfit(folds, energies, 2)
        E_quad0 = float(np.polyval(quad_fit, 0))
        quad_kcal = abs(E_quad0 - exact) * HARTREE_TO_KCAL_MOL

        # exponential fit: E(L) = A*exp(-k*L) + E_inf, via linearizing on
        # (E - E_last) if monotonic, else report as unavailable rather than fake it
        try:
            diffs = np.array(energies) - energies[-1]
            if np.all(diffs > 0) or np.all(diffs < 0):
                slope, intercept = np.polyfit(folds, np.log(np.abs(diffs)), 1)
                A = np.sign(diffs[0]) * np.exp(intercept)
                E_exp0 = energies[-1] + A
                exp_kcal = abs(E_exp0 - exact) * HARTREE_TO_KCAL_MOL
                exp_note = None
            else:
                exp_kcal = None
                exp_note = "energies not monotonic across folds -- exponential fit not well-posed, not faked"
        except Exception as e:
            exp_kcal = None
            exp_note = f"exponential fit failed: {e}"

        print(f"\n  no mitigation (fold=1) : {no_mit_kcal:.2f} kcal/mol")
        print(f"  ZNE linear              : {lin_kcal:.2f} kcal/mol")
        print(f"  ZNE quadratic           : {quad_kcal:.2f} kcal/mol")
        print(f"  ZNE exponential         : {exp_kcal:.2f} kcal/mol" if exp_kcal is not None
              else f"  ZNE exponential         : n/a ({exp_note})")

        zne_report[model] = {
            "folds": folds, "energies_ha": energies, "f": f_vals,
            "effective_error_rate": eff_rate, "rate_spread": round(rate_spread, 6),
            "rate_consistent": rate_consistent,
            "no_mitigation_kcal": round(no_mit_kcal, 4),
            "zne_linear_kcal": round(lin_kcal, 4),
            "zne_quadratic_kcal": round(quad_kcal, 4),
            "zne_exponential_kcal": round(exp_kcal, 4) if exp_kcal is not None else None,
            "zne_exponential_note": exp_note,
            "raw_baseline_kcal_from_abstract_gate_run": RAW_BASELINE_KCAL[model],
        }

    results["zne_report"] = zne_report
    results["classical_floor_kcal"] = CLASSICAL_FLOOR_KCAL
    results["chemical_accuracy_kcal"] = CHEMICAL_ACCURACY_KCAL
    save_results(results)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--noise-model", choices=["ideal", "aria-1", "forte-1"])
    parser.add_argument("--fold", type=int, choices=[1, 3, 5])
    parser.add_argument("--assemble", action="store_true")
    args = parser.parse_args()

    if args.assemble:
        return assemble()
    if not args.noise_model or not args.fold:
        parser.error("pass --noise-model and --fold, or --assemble")
    if args.noise_model == "ideal" and args.fold != 1:
        parser.error("ideal is only run at fold=1 (correctness control)")
    return run_one(args.noise_model, args.fold)


if __name__ == "__main__":
    main()
