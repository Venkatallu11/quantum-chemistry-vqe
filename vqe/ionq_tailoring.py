#!/usr/bin/env python3
"""
ionq_tailoring.py — Multi-fragment molecular tailoring on real IonQ
circuits: H6 chain ground state reconstructed from overlapping fragments,
each solved by entanglement forging.
============================================================================
Generalizes ionq_run.py's single-H4-fragment approach to THREE fragments
of the covalent_fragment.py H6 layout (4-atom blocks, d=1.0 Ang):

  fragment A = atoms [0,1,2,3]   (4 electrons, 8 qubits -> two 4-qubit EF registers)
  fragment B = atoms [2,3,4,5]   (4 electrons, 8 qubits -> two 4-qubit EF registers)
  overlap    = atoms [2,3]       (2 electrons, 4 qubits -> two 2-qubit EF registers)
  E_tailored = E(A) + E(B) - E(overlap)     [inclusion-exclusion, covalent_fragment.py]

Classical references (already in the repo, covalent_fragment_results.json):
  H6 full exact              = -3.236066 Ha
  H6 tailored, 4-atom blocks = -3.231625 Ha   (method floor: 2.79 kcal/mol)

Three techniques cut real circuit count by ~10x vs the naive approach,
EACH independently verified against exact local math before use here (see
ef_fragment.py's docstring and the scratchpad verify_*.py scripts from
development):

  1. Real gauge (ef_fragment.real_gauge): eigsh's ground state is real up
     to a global phase; rotating it real means every Pauli label's cross
     term is provably either purely real (even Y-count) or purely
     imaginary (odd Y-count) -- never both -- so only 2 of the 4 phase
     circuits (k=0,2 for this Hamiltonian, where every label turned out
     even-Y) are needed per Schmidt pair, not 4.

  2. Beta-register sign reuse (verified: v_n = s_n * u_n exactly, to
     ~1e-13, for H4): the beta register's Schmidt vectors are just
     signed copies of the alpha ones for this real-orbital, closed-shell
     Hamiltonian. B_nm = s_n * s_m * A_nm exactly -- so ONLY the alpha
     register needs real measurement; beta is derived classically for
     free. This is checked per-fragment (not assumed) before being relied
     on -- see verify_fragment_symmetries() below.

  3. Qubit-wise commuting Pauli grouping (verified: max error 3e-14 vs
     exact Statevector math across 222 sampled combinations): labels that
     agree on every qubit's measurement axis share ONE circuit instead of
     one each (37 alpha labels -> 13 groups on H4-shaped fragments).

Only fold=1 is run -- gate folding was proven to produce zero noise
response on IonQ's simulator (ionq_fold_check.py: a Bell state folded 1x
to 81x gave a bit-identical result every time), so ZNE would waste
circuits for no benefit here.

Reporting includes the "surviving signal fraction" f per fragment:
  E(f) = E_mixed + f*(E_ideal - E_mixed)
  => f = (E_measured - E_mixed) / (E_ideal - E_mixed)
E_mixed = (coefficient of the identity Pauli term) + enuc -- the energy a
fully-decohered state would give (every non-identity expectation -> 0).
f=1 means the real measurement fully preserved the ideal signal; f=0 means
it collapsed all the way to the mixed-state value.

Run (each --fragment is independent and resumable, matching ionq_run.py's
pattern -- merges into vqe/ionq_tailoring_results*.json, namespaced by
noise_model so different profiles never overwrite each other):
    python vqe/ionq_tailoring.py --fragment A --noise-model ideal --K 3
    python vqe/ionq_tailoring.py --fragment B --noise-model ideal --K 3
    python vqe/ionq_tailoring.py --fragment overlap --noise-model ideal --K 3
    python vqe/ionq_tailoring.py --assemble --noise-model ideal --K 3
    # ^ MUST match the classical tailored energy before running any noisy config.
    python vqe/ionq_tailoring.py --fragment A --noise-model aria-1 --K 3
    ...
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
from ionq_run import (
    transpiled_state_prep, fold_circuit, basis_change, pauli_expectation,
    submit_and_get_probabilities,
)
from ionq_backend import connect_provider, get_simulator

HARTREE_TO_KCAL_MOL = 627.5094740631

FRAGMENTS = {
    "A": {"atoms": [0, 1, 2, 3], "nelec": 4},
    "B": {"atoms": [2, 3, 4, 5], "nelec": 4},
    "overlap": {"atoms": [2, 3], "nelec": 2},
}
ASSEMBLY_SIGN = {"A": +1, "B": +1, "overlap": -1}  # inclusion-exclusion

CLASSICAL_H6_FULL_EXACT_HA = -3.236066
CLASSICAL_H6_TAILORED_4ATOM_HA = -3.231625

DEFAULT_K = 3


# ---------------------------------------------------------------------------
# Local (no-network) setup + verification for one fragment
# ---------------------------------------------------------------------------

def setup_fragment(frag_key, K):
    spec = FRAGMENTS[frag_key]
    atoms, nelec = spec["atoms"], spec["nelec"]

    qop_bare, qop_pen, enuc = effrag.build_fragment_qop(atoms, nelec, d=1.0)
    n_qubits = qop_bare.num_qubits
    e_elec, psi = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    psi_real, residual_imag = effrag.real_gauge(psi)

    K_max = 2 ** (n_qubits // 2)
    K_run = min(K, K_max)

    lambdas, u_vecs, v_vecs = effrag.schmidt_decompose_real(psi_real, n_qubits)
    terms = effrag.decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))

    # --- verify the two symmetry shortcuts hold for THIS fragment, not just H4 ---
    assert set(alpha_labels) == set(beta_labels), (
        f"fragment {frag_key}: alpha/beta label sets differ -- beta-reuse shortcut "
        "does not apply, would need independent beta-register measurement"
    )
    non_even_labels = [l for l in alpha_labels if not effrag.label_is_real(l)]
    assert not non_even_labels, (
        f"fragment {frag_key}: labels with odd Y-count found ({non_even_labels}) -- "
        "the 'always k=0,2' shortcut doesn't hold, some labels need k=1,3 instead"
    )

    u_top, v_top = u_vecs[:K_run], v_vecs[:K_run]
    signs = np.array([1.0 if np.dot(v_top[n], u_top[n]) >= 0 else -1.0 for n in range(K_run)])
    sign_residual = max(
        float(np.max(np.abs(v_top[n] - signs[n] * u_top[n]))) for n in range(K_run)
    )
    use_beta_reuse = sign_residual < 1e-8
    if not use_beta_reuse:
        print(f"  NOTE: fragment {frag_key} does NOT satisfy beta=sign*alpha "
              f"(residual {sign_residual:.3e}) -- measuring beta register "
              "independently instead of reusing alpha (verified per-fragment, "
              "not assumed from H4/fragment A).")

    # identity-term coefficient, for the E_mixed signal-fraction baseline
    identity_label = "I" * n_qubits
    identity_coeff = 0.0
    for label, coeff in qop_bare.to_list():
        if label == identity_label:
            identity_coeff = float(coeff.real)
            break
    e_mixed = identity_coeff + enuc

    alpha_cache, beta_cache = effrag.precompute_exact_matrices(terms, u_top, v_top)
    noiseless_numpy = effrag.ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K_run)

    groups = effrag.group_labels_qubit_wise(alpha_labels)

    return {
        "atoms": atoms, "nelec": nelec, "n_qubits": n_qubits,
        "enuc": enuc, "exact_energy": exact_energy, "real_gauge_residual": residual_imag,
        "lambdas": lambdas, "u_vecs": u_top, "v_vecs": v_top,
        "terms": terms, "alpha_labels": alpha_labels, "beta_labels": beta_labels,
        "signs": signs, "use_beta_reuse": use_beta_reuse,
        "K_run": K_run, "e_mixed": e_mixed, "noiseless_numpy": noiseless_numpy,
        "groups": groups,
    }


# ---------------------------------------------------------------------------
# Real IonQ measurement (one register; the caller decides whether the OTHER
# register needs its own call or can be derived via sign reuse)
# ---------------------------------------------------------------------------

def measure_register_matrices_real(vecs, groups, labels, backend, noise_model, K):
    """Mirrors ionq_run.py's build_matrices_real, but: (a) only THIS
    register (caller decides whether the other one is reused via signs or
    measured independently), (b) only k=0,2 phase circuits (real gauge,
    all-even-Y labels, verified per-fragment), (c) grouped measurement
    circuits (qubit-wise commuting) instead of one circuit per label.
    fold=1 always -- folding was proven not to scale noise on this
    simulator."""
    matrices = {label: np.zeros((K, K), dtype=complex) for label in labels}
    n_circuits = 0

    # --- diagonal: one job per Schmidt index, one circuit per group ---
    for n in range(K):
        base = fold_circuit(transpiled_state_prep(vecs[n], backend), fold=1)
        circuits, metas = [], []
        for group in groups:
            combined = effrag.combined_basis_label(group)
            qc = base.copy()
            basis_change(qc, combined)
            qc.measure_all()
            circuits.append(qc)
            metas.append(group)
        probs_list = submit_and_get_probabilities(circuits, backend, noise_model)
        n_circuits += len(circuits)
        for group, probs in zip(metas, probs_list):
            for label in group:
                matrices[label][n, n] = pauli_expectation(probs, label)

    # --- cross: one job per pair, k=0 and k=2 state preps, one circuit per group each ---
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        circuits, metas = [], []
        for k, sign in ((0, 1.0), (2, -1.0)):
            vec_k = (vecs[n] + sign * vecs[m]) / np.sqrt(2)
            base = fold_circuit(transpiled_state_prep(vec_k, backend), fold=1)
            for group in groups:
                combined = effrag.combined_basis_label(group)
                qc = base.copy()
                basis_change(qc, combined)
                qc.measure_all()
                circuits.append(qc)
                metas.append((k, group))
        probs_list = submit_and_get_probabilities(circuits, backend, noise_model)
        n_circuits += len(circuits)
        e_by_label_k = {label: {} for label in labels}
        for (k, group), probs in zip(metas, probs_list):
            for label in group:
                e_by_label_k[label][k] = pauli_expectation(probs, label)
        for label in labels:
            e0, e2 = e_by_label_k[label][0], e_by_label_k[label][2]
            re = (e0 - e2) / 2  # always real for these labels (verified even-Y)
            matrices[label][n, m] = re
            matrices[label][m, n] = re

    return matrices, n_circuits


def derive_beta_matrices(alpha_matrices, signs, K):
    """B_nm = s_n * s_m * A_nm (verified exactly for these fragments)."""
    S = np.diag(signs)
    return {label: S @ mat @ S for label, mat in alpha_matrices.items()}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def results_path(noise_model):
    suffix = "" if noise_model == "ideal" else f"_{noise_model}"
    return os.path.join(os.path.dirname(__file__), f"ionq_tailoring_results{suffix}.json")


def load_results(noise_model):
    path = results_path(noise_model)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_results(results, noise_model):
    path = results_path(noise_model)
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {path}\n")


def run_one_fragment(frag_key, noise_model, K):
    print("\n" + "=" * 70)
    print(f"  IonQ tailoring -- fragment={frag_key}, noise_model={noise_model}, K={K}")
    print("=" * 70)

    if noise_model != "ideal":
        ideal_results = load_results("ideal")
        ideal_frag = ideal_results.get("fragments", {}).get(frag_key)
        if not ideal_frag or ideal_frag.get("requested_K") != K:
            print(f"\n  REFUSING: no matching --noise-model ideal --K {K} result for "
                  f"fragment {frag_key} yet. Run and confirm ideal first -- if ideal "
                  "doesn't match the classical reference, the pipeline is broken and "
                  "a noisy run would just be measuring a bug.")
            sys.exit(1)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  Connected. Target backend: {backend.name}")

    p = setup_fragment(frag_key, K)
    print(f"  atoms={p['atoms']}, nelec={p['nelec']}, {p['n_qubits']} qubits, "
          f"K={p['K_run']} (requested {K})")
    print(f"  real-gauge residual imag: {p['real_gauge_residual']:.3e}")
    print(f"  {len(p['alpha_labels'])} alpha labels -> {len(p['groups'])} qubit-wise groups")
    print(f"  exact fragment energy = {p['exact_energy']:.6f} Ha")
    print(f"  exact numpy EF energy (K={p['K_run']}) = {p['noiseless_numpy']:.6f} Ha")

    alpha_mats, n_circuits = measure_register_matrices_real(
        p["u_vecs"], p["groups"], p["alpha_labels"], backend, noise_model, p["K_run"],
    )
    if p["use_beta_reuse"]:
        beta_mats = derive_beta_matrices(alpha_mats, p["signs"], p["K_run"])
    else:
        beta_mats, n_circuits_beta = measure_register_matrices_real(
            p["v_vecs"], p["groups"], p["beta_labels"], backend, noise_model, p["K_run"],
        )
        n_circuits += n_circuits_beta
    E = effrag.ef_energy_from_noisy_matrices(
        p["terms"], p["lambdas"], alpha_mats, beta_mats, p["enuc"], p["K_run"],
    )

    f_signal = (E - p["e_mixed"]) / (p["noiseless_numpy"] - p["e_mixed"]) \
        if abs(p["noiseless_numpy"] - p["e_mixed"]) > 1e-9 else None

    err_vs_own_exact_kcal = abs(E - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    print(f"\n  E({frag_key}) = {E:.6f} Ha, {n_circuits} circuits, "
          f"error vs own exact = {err_vs_own_exact_kcal:.3f} kcal/mol, "
          f"f = {f_signal}")

    results = load_results(noise_model)
    results.setdefault("fragments", {})
    results["fragments"][frag_key] = {
        "atoms": p["atoms"], "nelec": p["nelec"], "n_qubits": p["n_qubits"],
        "requested_K": K, "K": p["K_run"], "n_circuits": n_circuits,
        "exact_energy_ha": round(p["exact_energy"], 6),
        "exact_numpy_ef_k_ha": round(p["noiseless_numpy"], 6),
        "e_mixed_ha": round(p["e_mixed"], 6),
        "measured_energy_ha": round(E, 6),
        "error_vs_own_exact_kcal": round(err_vs_own_exact_kcal, 4),
        "surviving_signal_fraction": round(f_signal, 6) if f_signal is not None else None,
        "noise_model": noise_model, "backend": backend.name,
    }
    save_results(results, noise_model)
    return results


def assemble(noise_model, K):
    results = load_results(noise_model)
    fragments = results.get("fragments", {})
    print(f"\n  Fragments completed for noise_model={noise_model}: {sorted(fragments.keys())}")

    missing = [k for k in FRAGMENTS if k not in fragments or fragments[k]["requested_K"] != K]
    if missing:
        print(f"  Missing (or wrong K) fragments: {missing} -- run those first.")
        return results

    E_tailored = sum(ASSEMBLY_SIGN[k] * fragments[k]["measured_energy_ha"] for k in FRAGMENTS)
    E_tailored_numpy_ef_k = sum(ASSEMBLY_SIGN[k] * fragments[k]["exact_numpy_ef_k_ha"] for k in FRAGMENTS)
    total_circuits = sum(fragments[k]["n_circuits"] for k in FRAGMENTS)

    err_vs_tailored_kcal = abs(E_tailored - CLASSICAL_H6_TAILORED_4ATOM_HA) * HARTREE_TO_KCAL_MOL
    err_vs_full_kcal = abs(E_tailored - CLASSICAL_H6_FULL_EXACT_HA) * HARTREE_TO_KCAL_MOL
    # pipeline-correctness check: does the REAL measurement match what this
    # exact SAME K-truncated method predicts with zero hardware noise? This is
    # the honest "is the pipeline broken" signal -- separate from the EF-K
    # method's own truncation error vs the fully-exact classical reference,
    # which is a real, expected property of a low Schmidt rank, not a bug.
    err_vs_numpy_ef_k_kcal = abs(E_tailored - E_tailored_numpy_ef_k) * HARTREE_TO_KCAL_MOL

    print(f"\n  {'fragment':<10} {'E (Ha)':>12} {'err vs own exact (kcal/mol)':>28} {'f':>10} {'circuits':>10}")
    for k in FRAGMENTS:
        fr = fragments[k]
        f_str = f"{fr['surviving_signal_fraction']:.4f}" if fr["surviving_signal_fraction"] is not None else "n/a"
        print(f"  {k:<10} {fr['measured_energy_ha']:>12.6f} {fr['error_vs_own_exact_kcal']:>28.3f} "
              f"{f_str:>10} {fr['n_circuits']:>10}")

    print(f"\n  E_tailored (assembled)        = {E_tailored:.6f} Ha")
    print(f"  exact-numpy EF-K={K} tailored   = {E_tailored_numpy_ef_k:.6f} Ha "
          f"(pipeline-correctness gap: {err_vs_numpy_ef_k_kcal:.3f} kcal/mol)")
    print(f"  classical tailored (exact frags) = {CLASSICAL_H6_TAILORED_4ATOM_HA:.6f} Ha "
          f"(EF-K={K} method-truncation gap: {err_vs_tailored_kcal:.3f} kcal/mol)")
    print(f"  classical full exact             = {CLASSICAL_H6_FULL_EXACT_HA:.6f} Ha "
          f"(error {err_vs_full_kcal:.3f} kcal/mol)")
    print(f"  total circuits submitted across all fragments: {total_circuits}")

    if noise_model == "ideal":
        # the real gate: does the real hardware round-trip reproduce what the
        # SAME K-truncated method predicts noiselessly? (not: does the method
        # itself match full exact diagonalization -- that's a separate, honest
        # truncation-error report, not a pipeline-correctness signal)
        ok = err_vs_numpy_ef_k_kcal < 2.0
        print(f"\n  IDEAL SANITY CHECK (real vs exact-numpy EF-K={K}, same method): "
              f"{'PASS' if ok else 'FAIL -- pipeline is broken, stop here'}"
              f" ({err_vs_numpy_ef_k_kcal:.3f} kcal/mol)")
        print(f"  (separately, EF-K={K} truncation vs fully-exact classical tailoring: "
              f"{err_vs_tailored_kcal:.3f} kcal/mol -- a real method limitation, not a bug)")

    results["assembly"] = {
        "K": K,
        "E_tailored_ha": round(E_tailored, 6),
        "exact_numpy_ef_k_tailored_ha": round(E_tailored_numpy_ef_k, 6),
        "pipeline_correctness_gap_kcal": round(err_vs_numpy_ef_k_kcal, 4),
        "classical_tailored_ha": CLASSICAL_H6_TAILORED_4ATOM_HA,
        "classical_full_exact_ha": CLASSICAL_H6_FULL_EXACT_HA,
        "error_vs_classical_tailored_kcal": round(err_vs_tailored_kcal, 4),
        "error_vs_classical_full_exact_kcal": round(err_vs_full_kcal, 4),
        "total_circuits": total_circuits,
    }
    save_results(results, noise_model)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragment", choices=list(FRAGMENTS.keys()),
                         help="run exactly one fragment and merge into the results JSON")
    parser.add_argument("--assemble", action="store_true",
                         help="no network calls -- assemble E_tailored from saved fragments")
    parser.add_argument("--noise-model", default="ideal",
                         help="ideal | aria-1 | forte-1 (fold=1 always -- folding doesn't "
                              "scale noise on this simulator, proven separately)")
    parser.add_argument("--K", type=int, default=DEFAULT_K)
    args = parser.parse_args()

    if args.assemble:
        assemble(args.noise_model, args.K)
    elif args.fragment:
        run_one_fragment(args.fragment, args.noise_model, args.K)
    else:
        parser.error("pass --fragment <A|B|overlap> or --assemble")


if __name__ == "__main__":
    main()
