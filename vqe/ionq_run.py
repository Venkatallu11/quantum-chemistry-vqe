#!/usr/bin/env python3
"""
ionq_run.py — Entanglement Forging (K=5) for the H4 fragment, measured via
real circuit submission to IonQ's cloud simulator (free tier only:
ionq_simulator -- never ionq_qpu, which is real paid trapped-ion hardware
and deliberately out of scope for this script).
============================================================================
Reuses entanglement_forging_h4.py's Hamiltonian / Schmidt decomposition /
Pauli-term machinery unchanged -- H4 fragment (atoms 0-3, d=1.0 Ang), 185
Pauli terms, 37 unique alpha-register labels, 37 unique beta-register
labels. The only thing that's new here is HOW the alpha/beta matrix
elements get measured:

  entanglement_forging_h4.py (Stage 2): local qiskit-aer density-matrix
  simulator + a hand-set depolarizing noise model standing in for
  Quantinuum-grade hardware.

  ionq_run.py (this file): real circuits, actually submitted to IonQ's
  cloud simulator, either with noise_model="ideal" (exact probabilities --
  a sanity check that the real round-trip matches the known numpy answer)
  or noise_model="aria-1" (IonQ's own noise emulation of their real Aria
  trapped-ion device). Because it's a SIMULATOR, IonQ's API returns exact
  probabilities either way (see get_probabilities()), not shot-sampled
  counts -- so these numbers carry no shot noise, only the noise_model's
  physics, same determinism as the local Aer density-matrix approach.

ZNE here can't reuse entanglement_forging_zne.py's trick of scaling a
LOCAL noise model's error rate by a continuous factor -- IonQ's simulator
only exposes a fixed, named noise_model. Instead this uses the standard
real-hardware ZNE technique: local unitary gate folding. Every 2-qubit
gate G in the (already-transpiled) circuit is replaced with
G (G^-1 G)^((fold-1)/2), which physically stretches the circuit to
1x/3x/5x depth -- and therefore 1x/3x/5x noise-channel applications --
without changing the ideal unitary. Folding is applied AFTER
transpilation, right before submission, so no further transpiler pass
gets a chance to cancel the inserted gate pairs back out (IonQ's own
circuit converter does not re-optimize submitted circuits).

Circuit count is real: the full 37-label-per-register K=5 reconstruction
is 5 diagonal + 10x4=40 cross state preparations per register, x 37
labels each = ~1665 circuits per register x 2 registers = ~3330 circuits
per (noise_model, fold) configuration. Submitted as 15 batched jobs per
register (one per Schmidt index for diagonal, one per pair for cross)
rather than one giant job, so a single bad job doesn't blow up the run.

Run:
    python vqe/ionq_run.py --stage ideal --smoke   # tiny pipeline test
    python vqe/ionq_run.py --stage ideal           # full sanity check vs exact numpy
    python vqe/ionq_run.py --stage zne             # aria-1, fold=1,3,5 + extrapolation
    python vqe/ionq_run.py --stage all             # both
"""
import os
import sys
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit import transpile
from qiskit.transpiler import CouplingMap

from ionq_backend import connect_provider, get_simulator

HARTREE_TO_KCAL_MOL = 627.5094740631
K = 5
FOLD_FACTORS = [1, 3, 5]
DEFAULT_SHOTS = 100_000  # matters! IonQ's simulator "probabilities" carry real
# shot noise for anything beyond trivial 2-outcome distributions -- confirmed
# empirically: shots=100 gave ~0.01 error on a near-zero expectation value,
# shots>=1000 matched the exact numpy answer to machine precision. 100k keeps
# per-term noise well under the ~0.001 Ha (~1 kcal/mol) precision this project
# cares about.
NOISY_MODEL = "aria-1"

# Standard qiskit gate names (a subset of IonQ's 'qis' gateset -- excludes
# IonQ-only aliases like "cnot"/"mcx" that qiskit's transpiler basis_gates
# parser doesn't recognize as real gate classes). rx/ry/rz + cx is already
# universal, so this is sufficient for arbitrary 4-qubit state preparation.
IONQ_QIS_STANDARD_BASIS = [
    "rx", "ry", "rz", "h", "x", "y", "z",
    "s", "sdg", "sx", "sxdg", "t", "tdg",
    "cx", "cy", "cz", "swap",
]


# ---------------------------------------------------------------------------
# Circuit construction: state prep -> transpile -> fold -> basis change -> measure
# ---------------------------------------------------------------------------

def state_prep_circuit(vec):
    qc = QuantumCircuit(4)
    qc.append(StatePreparation(np.asarray(vec)), range(4))
    return qc


def transpiled_state_prep(vec, backend):
    """Decompose StatePreparation into the backend's own basis gates
    (IonQ's 'qis' gateset), on a clean 4-qubit fully-connected coupling
    map -- NOT `transpile(..., backend=backend)`, which lays the circuit
    out across the backend's full qubit count (29 for ionq_simulator) and
    can place the 4 logical qubits on arbitrary physical indices, silently
    breaking any code (like basis_change()) that assumes qubit i is at
    position i."""
    cmap = CouplingMap.from_full(4)
    return transpile(state_prep_circuit(vec), basis_gates=IONQ_QIS_STANDARD_BASIS,
                      coupling_map=cmap, optimization_level=1)


def fold_circuit(qc, fold):
    """Local unitary folding on every 2-qubit gate: G -> G (G^-1 G)^reps,
    physically stretching depth (and real device noise) by `fold`x while
    leaving the ideal unitary unchanged. fold must be odd; fold=1 is a
    no-op copy. Folding after transpilation means IonQ's circuit
    converter (a straight instruction-by-instruction dict conversion,
    no optimization passes) submits it exactly as folded."""
    if fold == 1:
        return qc.copy()
    assert fold % 2 == 1, "fold factor must be odd"
    reps = (fold - 1) // 2
    folded = qc.copy_empty_like()
    for instr in qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        folded.append(op, qargs, cargs)
        if op.num_qubits == 2 and op.name not in ("measure", "barrier"):
            for _ in range(reps):
                folded.append(op.inverse(), qargs, cargs)
                folded.append(op, qargs, cargs)
    return folded


def basis_change(qc, label):
    """Append per-qubit basis-change gates so a Z-basis measurement reads
    out the eigenbasis of Pauli string `label` (qiskit label convention:
    label[i] acts on qubit n-1-i, i.e. leftmost char = highest qubit)."""
    n = len(label)
    for i, ch in enumerate(label):
        qubit = n - 1 - i
        if ch == "X":
            qc.h(qubit)
        elif ch == "Y":
            qc.sdg(qubit)
            qc.h(qubit)


def measurement_circuit(base_circuit, label):
    qc = base_circuit.copy()
    basis_change(qc, label)
    qc.measure_all()
    return qc


# ---------------------------------------------------------------------------
# Real-circuit expectation values from IonQ probabilities
# ---------------------------------------------------------------------------

def pauli_expectation(probs, label):
    """probs: dict[bitstring -> probability] (already normalized to 1.0).
    Standard qiskit bit ordering: bitstring[-1] is qubit 0."""
    n = len(label)
    exp = 0.0
    for bitstring, p in probs.items():
        bits = bitstring[::-1]
        sign = 1
        for q in range(n):
            ch = label[n - 1 - q]
            if ch != "I" and q < len(bits) and bits[q] == "1":
                sign = -sign
        exp += sign * p
    return exp


def submit_and_get_probabilities(circuits, backend, noise_model, shots=DEFAULT_SHOTS):
    """Submit one batched real job to IonQ and return a list of
    probability dicts, one per circuit. Real API call -- no local
    fallback, no faking a result if this errors."""
    job = backend.run(circuits, noise_model=noise_model, shots=shots)
    result = job.result()
    probs = result.get_probabilities()
    if not isinstance(probs, list):
        probs = [probs]
    return [dict(p) for p in probs]


# ---------------------------------------------------------------------------
# Matrix-element reconstruction (mirrors entanglement_forging_h4.build_noisy_matrices)
# ---------------------------------------------------------------------------

def build_matrices_real(vecs, unique_labels, backend, noise_model, fold, K):
    unique_labels = sorted(set(unique_labels))
    matrices = {label: np.zeros((K, K), dtype=complex) for label in unique_labels}

    # --- diagonal: one job per Schmidt index n ---
    for n in range(K):
        base = fold_circuit(transpiled_state_prep(vecs[n], backend), fold)
        circuits = [measurement_circuit(base, label) for label in unique_labels]
        probs_list = submit_and_get_probabilities(circuits, backend, noise_model)
        for label, probs in zip(unique_labels, probs_list):
            matrices[label][n, n] = pauli_expectation(probs, label)

    # --- cross: one job per (n, m) pair, 4 phase circuits x n_labels ---
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        circuits, metas = [], []
        for k in range(4):
            vec_k = (vecs[n] + (1j ** k) * vecs[m]) / np.sqrt(2)
            base = fold_circuit(transpiled_state_prep(vec_k, backend), fold)
            for label in unique_labels:
                circuits.append(measurement_circuit(base, label))
                metas.append((k, label))
        probs_list = submit_and_get_probabilities(circuits, backend, noise_model)
        vals = {label: [None] * 4 for label in unique_labels}
        for (k, label), probs in zip(metas, probs_list):
            vals[label][k] = pauli_expectation(probs, label)
        for label in unique_labels:
            e0, e1, e2, e3 = vals[label]
            re = (e0 - e2) / 2
            im = (e3 - e1) / 2
            x = re + 1j * im
            matrices[label][n, m] = x
            matrices[label][m, n] = np.conj(x)

    return matrices


def ef_energy_from_matrices_real(terms, lambdas, alpha_mats, beta_mats, enuc, K):
    return ef.ef_energy_from_noisy_matrices(terms, lambdas, alpha_mats, beta_mats, enuc, K)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_stage(name, backend, noise_model, fold, terms, lambdas, u_vecs, v_vecs,
              enuc, exact_energy, alpha_labels, beta_labels, K):
    print(f"\n  -- {name}: noise_model={noise_model}, fold={fold} --")
    alpha_mats = build_matrices_real(u_vecs, alpha_labels, backend, noise_model, fold, K)
    print(f"    alpha register done ({len(alpha_labels)} labels)")
    beta_mats = build_matrices_real(v_vecs, beta_labels, backend, noise_model, fold, K)
    print(f"    beta register done ({len(beta_labels)} labels)")
    E = ef_energy_from_matrices_real(terms, lambdas, alpha_mats, beta_mats, enuc, K)
    err_ha = abs(E - exact_energy)
    err_kcal = err_ha * HARTREE_TO_KCAL_MOL
    print(f"    E = {E:.6f} Ha, error = {err_ha:.6f} Ha ({err_kcal:.3f} kcal/mol)")
    return E, err_ha, err_kcal


def results_path(smoke):
    suffix = "_smoke" if smoke else ""
    return os.path.join(os.path.dirname(__file__), f"ionq_run_results{suffix}.json")


def load_results(smoke):
    path = results_path(smoke)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_results(results, smoke):
    path = results_path(smoke)
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {path}\n")


def setup_problem(smoke):
    """Shared, cheap (local, no network) setup: Hamiltonian, Schmidt
    decomposition, term/label lists, and the exact-numpy reference. Every
    --config invocation redoes this (fast) so each one is a fully
    independent, resumable process -- no shared state between runs."""
    qop_bare, qop_penalized, enuc = ef.build_h4_qop(d=1.0)
    e_elec, psi = ef.exact_ground_state(qop_penalized)
    exact_energy = e_elec + enuc
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))

    K_run = K
    if smoke:
        K_run = 2
        alpha_labels = alpha_labels[:3]
        beta_labels = beta_labels[:3]
        alpha_set, beta_set = set(alpha_labels), set(beta_labels)
        terms = [(a, b, c) for (a, b, c) in terms if a in alpha_set and b in beta_set]

    alpha_cache, beta_cache = ef.precompute_exact_matrices(terms, u_vecs[:K_run], v_vecs[:K_run])
    noiseless_numpy = ef.ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K_run)

    return {
        "enuc": enuc, "exact_energy": exact_energy, "lambdas": lambdas,
        "u_vecs": u_vecs, "v_vecs": v_vecs, "terms": terms,
        "alpha_labels": alpha_labels, "beta_labels": beta_labels,
        "K_run": K_run, "noiseless_numpy": noiseless_numpy,
    }


CONFIGS = {
    "ideal": ("ideal", 1),
    "fold1": (NOISY_MODEL, 1),
    "fold3": (NOISY_MODEL, 3),
    "fold5": (NOISY_MODEL, 5),
}


def run_one_config(config_key, smoke):
    """Run exactly ONE (noise_model, fold) configuration end to end and
    merge its result into the persistent results JSON. Small, bounded
    runtime by design -- meant to be invoked repeatedly (once per config)
    rather than as one long-running process, so a kill/interruption only
    loses one config's worth of work, not the whole run."""
    noise_model, fold = CONFIGS[config_key]

    print("\n" + "=" * 70)
    print(f"  IonQ EF run -- config={config_key} (noise_model={noise_model}, fold={fold})")
    print("=" * 70)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  Connected. Target backend: {backend.name}")

    p = setup_problem(smoke)
    numpy_err_kcal = abs(p["noiseless_numpy"] - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    print(f"  exact ground state = {p['exact_energy']:.6f} Ha")
    print(f"  K = {p['K_run']}, {len(p['alpha_labels'])} alpha labels, "
          f"{len(p['beta_labels'])} beta labels, {len(p['terms'])} terms")
    print(f"  exact numpy EF energy (K={p['K_run']}) = {p['noiseless_numpy']:.6f} Ha "
          f"(error {numpy_err_kcal:.3f} kcal/mol)")

    E, err_ha, err_kcal = run_stage(
        config_key, backend, noise_model, fold,
        p["terms"], p["lambdas"], p["u_vecs"], p["v_vecs"], p["enuc"], p["exact_energy"],
        p["alpha_labels"], p["beta_labels"], p["K_run"],
    )

    results = load_results(smoke)
    results.update({
        "molecule": "H4 fragment (atoms 0-3, d=1.0 Ang)",
        "K": p["K_run"],
        "n_alpha_labels": len(p["alpha_labels"]),
        "n_beta_labels": len(p["beta_labels"]),
        "n_terms": len(p["terms"]),
        "exact_energy_ha": round(p["exact_energy"], 6),
        "exact_numpy_ef_k_ha": round(p["noiseless_numpy"], 6),
        "exact_numpy_ef_k_err_kcal": round(numpy_err_kcal, 4),
        "smoke_test": smoke,
        "backend": backend.name,
    })
    results.setdefault("configs", {})
    results["configs"][config_key] = {
        "noise_model": noise_model, "fold": fold,
        "energy_ha": round(E, 6), "error_ha": round(err_ha, 6),
        "error_kcal_mol": round(err_kcal, 4),
    }
    save_results(results, smoke)
    return results


def assemble(smoke):
    """No network calls -- just load whatever configs have been saved so
    far and, if all three aria-1 folds are present, fit + extrapolate."""
    results = load_results(smoke)
    configs = results.get("configs", {})
    print("\n  Configs completed so far:", sorted(configs.keys()))

    fold_energies = {}
    for key, fold in (("fold1", 1), ("fold3", 3), ("fold5", 5)):
        if key in configs:
            fold_energies[fold] = configs[key]["energy_ha"]

    if len(fold_energies) < 2:
        print("  Not enough aria-1 fold results yet for ZNE extrapolation.")
        return results

    folds = sorted(fold_energies)
    energies = [fold_energies[f] for f in folds]
    exact_energy = results["exact_energy_ha"]

    no_mit_err = abs(energies[0] - exact_energy) * HARTREE_TO_KCAL_MOL
    zne_result = {
        "fold_factors": folds,
        "energies_ha": energies,
        "no_mitigation_err_kcal": round(no_mit_err, 4),
    }
    print(f"  no mitigation (fold={folds[0]}) : {no_mit_err:.2f} kcal/mol")

    lin_fit = np.polyfit(folds, energies, 1)
    E_lin0 = float(np.polyval(lin_fit, 0))
    zne_lin_err = abs(E_lin0 - exact_energy) * HARTREE_TO_KCAL_MOL
    zne_result["zne_linear_energy_ha"] = round(E_lin0, 6)
    zne_result["zne_linear_err_kcal"] = round(zne_lin_err, 4)
    print(f"  ZNE linear extrapolation      : {zne_lin_err:.2f} kcal/mol")

    if len(folds) >= 3:
        quad_fit = np.polyfit(folds, energies, 2)
        E_quad0 = float(np.polyval(quad_fit, 0))
        zne_quad_err = abs(E_quad0 - exact_energy) * HARTREE_TO_KCAL_MOL
        zne_result["zne_quadratic_energy_ha"] = round(E_quad0, 6)
        zne_result["zne_quadratic_err_kcal"] = round(zne_quad_err, 4)
        print(f"  ZNE quadratic extrapolation   : {zne_quad_err:.2f} kcal/mol")

    results["zne"] = zne_result
    save_results(results, smoke)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", choices=list(CONFIGS.keys()),
                         help="run exactly one (noise_model, fold) configuration "
                              "and merge into the results JSON")
    parser.add_argument("--assemble", action="store_true",
                         help="no network calls -- fit/extrapolate ZNE from "
                              "whatever configs are already saved")
    parser.add_argument("--smoke", action="store_true",
                         help="tiny K=2, 3 labels/register -- pipeline validation only")
    args = parser.parse_args()

    if args.assemble:
        assemble(args.smoke)
    elif args.config:
        run_one_config(args.config, args.smoke)
    else:
        parser.error("pass --config <ideal|fold1|fold3|fold5> or --assemble")


if __name__ == "__main__":
    main()
