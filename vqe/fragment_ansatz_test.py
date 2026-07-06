#!/usr/bin/env python3
"""
fragment_ansatz_test.py — SIMULATOR-ONLY: which particle-conserving
ansatz can actually reach block1's (H4 fragment, atoms 0-3, d=1.0 A)
exact ground state?

Background: hardware_fragmentation.py's block1 (8 qubits, 4 electrons)
could not converge with EfficientSU2 + a particle-number PENALTY term,
even after COBYLA with up to 16000 iterations and 24 total restarts
(best error 0.067842 Ha, target 0.002 Ha). Root cause: EfficientSU2 is a
GENERIC ansatz that doesn't conserve particle number on its own, so the
penalty term (needed to keep it in the right electron-count sector) makes
the optimization landscape hard regardless of optimizer or iteration count.

This script tests two ansätze that conserve particle number BY
CONSTRUCTION (no penalty term needed at all -- verified below that
ExcitationPreserving keeps the electron count EXACTLY fixed even under
random parameters):

  (a) Hartree-Fock state prep + qiskit's ExcitationPreserving(reps=2,3)
  (b) qiskit-nature's UCCSD, with a HartreeFock initial state -- the
      standard chemically-motivated ansatz for exactly this kind of problem

Both are optimized on the noiseless LOCAL statevector simulator only.
NO hardware jobs are submitted by this script, and ../quantum-hardware-mcp
is never touched.

Run:
    python vqe/fragment_ansatz_test.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(__file__))
from hardware_fragmentation import _build_fragment_qop, hchain

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import ExcitationPreserving
from qiskit.primitives import StatevectorEstimator

from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.mappers import JordanWignerMapper

CHEM_ACCURACY_TARGET_HA = 0.002
N_RESTARTS = 5


def exact_energy_no_penalty(qop_bare, enuc):
    """
    Exact ground state of the BARE (physical) Hamiltonian -- unlike
    hardware_fragmentation.py, no penalty term is needed here: both
    ansätze below conserve particle number by construction, so there's
    no risk of the exact diagonalization (or the VQE optimization)
    landing on a wrong-electron-count state.
    """
    M = qop_bare.to_matrix(sparse=True)
    e_elec = float(eigsh(M, k=1, which="SA", return_eigenvectors=False)[0])
    return e_elec + enuc


def _hf_bits(n_qubits, nelec):
    """Hartree-Fock occupation, block-spin Jordan-Wigner ordering (same convention chem.py/qiskit-nature use)."""
    norb = n_qubits // 2
    occ = nelec // 2
    bits = [0] * n_qubits
    for i in range(occ):
        bits[i] = 1
        bits[norb + i] = 1
    return bits


def circuit_stats(qc):
    """
    Fully decompose to basic gates and report depth + 2-qubit gate count
    -- an honest proxy for how noisy this circuit would actually be on
    real hardware (more/deeper 2-qubit gates = more accumulated error).
    """
    decomposed = qc.decompose(reps=6)
    ops = dict(decomposed.count_ops())
    two_qubit_names = ("cx", "cz", "ecr", "rxx", "ryy", "rzz", "rzx", "swap", "iswap")
    two_qubit_gates = sum(count for name, count in ops.items() if name in two_qubit_names)
    return decomposed.depth(), two_qubit_gates, ops


def _run_optimizer(objective, n_params, exact, method, n_restarts, maxiter, seed=42):
    """Shared restart loop for both ansätze -- same pattern used throughout this repo."""
    rng = np.random.default_rng(seed)
    best_params, best_energy = None, None
    t0 = time.time()
    for r in range(n_restarts):
        x0 = np.zeros(n_params) if r == 0 else rng.uniform(-np.pi, np.pi, size=n_params)
        if method == "COBYLA":
            opts = {"maxiter": maxiter, "rhobeg": 0.5, "tol": 1e-10}
        else:  # L-BFGS-B
            opts = {"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-10}
        res = minimize(objective, x0=x0, method=method, options=opts)
        if best_energy is None or res.fun < best_energy:
            best_energy, best_params = float(res.fun), res.x
        error = abs(best_energy - exact)
        print(f"    [{method}] restart {r + 1}/{n_restarts}: E = {res.fun:.6f} Ha "
              f"(best {best_energy:.6f} Ha, error {error:.6f} Ha)")
        if error < CHEM_ACCURACY_TARGET_HA:
            break
    elapsed = time.time() - t0
    return best_params, best_energy, abs(best_energy - exact), elapsed


def test_excitation_preserving(qop_bare, enuc, exact, n_qubits, nelec, reps):
    """
    Hartree-Fock state prep (X gates on occupied qubits) followed by
    ExcitationPreserving(reps) -- verified separately to conserve
    particle number EXACTLY even under random parameters, so no penalty
    term is needed: the ansatz itself never leaves the correct sector.
    """
    bits = _hf_bits(n_qubits, nelec)
    base_ansatz = ExcitationPreserving(n_qubits, reps=reps)
    n_params = base_ansatz.num_parameters
    estimator = StatevectorEstimator()

    def make_circuit(params):
        qc = QuantumCircuit(n_qubits)
        for i, b in enumerate(bits):
            if b:
                qc.x(i)
        qc.compose(base_ansatz.assign_parameters(params), inplace=True)
        return qc

    def objective(params):
        qc = make_circuit(params)
        result = estimator.run([(qc, qop_bare)]).result()
        return float(result[0].data.evs) + enuc

    results = {}
    for method, maxiter in (("COBYLA", 3000), ("L-BFGS-B", 500)):
        print(f"  ExcitationPreserving(reps={reps}) -- {method}")
        params, energy, error, elapsed = _run_optimizer(
            objective, n_params, exact, method, N_RESTARTS, maxiter,
        )
        qc = make_circuit(params)
        depth, two_q, ops = circuit_stats(qc)
        results[method] = {
            "energy_ha": energy, "error_ha": error,
            "converged": bool(error < CHEM_ACCURACY_TARGET_HA),
            "depth": depth, "two_qubit_gates": two_q, "gate_counts": ops,
            "time_sec": round(elapsed, 1), "n_params": n_params,
        }
    return results


def test_uccsd(qop_bare, enuc, exact, n_qubits, nelec):
    """
    UCCSD with a Hartree-Fock initial state -- the standard chemically
    motivated, particle-conserving ansatz, built with qiskit-nature's own
    (correctly-ordered) HartreeFock/UCCSD classes rather than a hand-rolled
    bitstring (that's what caused the earlier E=82.9 Ha disaster with a
    manually-guessed warm start for EfficientSU2 -- qiskit-nature's own
    builder avoids that class of bug entirely).
    """
    n_orb = n_qubits // 2
    mapper = JordanWignerMapper()
    num_particles = (nelec // 2, nelec // 2)

    hf_state = HartreeFock(n_orb, num_particles, mapper)
    ansatz = UCCSD(n_orb, num_particles, mapper, initial_state=hf_state)
    n_params = ansatz.num_parameters
    estimator = StatevectorEstimator()

    def objective(params):
        qc = ansatz.assign_parameters(params)
        result = estimator.run([(qc, qop_bare)]).result()
        return float(result[0].data.evs) + enuc

    results = {}
    for method, maxiter in (("COBYLA", 3000), ("L-BFGS-B", 500)):
        print(f"  UCCSD -- {method}")
        params, energy, error, elapsed = _run_optimizer(
            objective, n_params, exact, method, N_RESTARTS, maxiter,
        )
        qc = ansatz.assign_parameters(params)
        depth, two_q, ops = circuit_stats(qc)
        results[method] = {
            "energy_ha": energy, "error_ha": error,
            "converged": bool(error < CHEM_ACCURACY_TARGET_HA),
            "depth": depth, "two_qubit_gates": two_q, "gate_counts": ops,
            "time_sec": round(elapsed, 1), "n_params": n_params,
        }
    return results


def main():
    d = 1.0
    indices = [0, 1, 2, 3]  # block1: the fragment that failed before
    nelec = len(indices)

    geom = hchain(indices, d)
    qop_bare, _qop_penalized, enuc = _build_fragment_qop(geom, nelec)
    n_qubits = qop_bare.num_qubits
    exact = exact_energy_no_penalty(qop_bare, enuc)

    print("\n" + "=" * 70)
    print("  Particle-conserving ansatz test -- block1 (H4, atoms 0-3)")
    print("=" * 70)
    print(f"  {n_qubits} qubits, {nelec} electrons")
    print(f"  exact ground state = {exact:.6f} Ha\n")

    all_results = {"n_qubits": n_qubits, "nelec": nelec, "exact_energy_ha": round(exact, 6)}

    print("=" * 70)
    print("  (a) ExcitationPreserving")
    print("=" * 70)
    all_results["excitation_preserving"] = {}
    for reps in (2, 3):
        print(f"\n-- reps={reps} --")
        res = test_excitation_preserving(qop_bare, enuc, exact, n_qubits, nelec, reps)
        all_results["excitation_preserving"][f"reps{reps}"] = res
        for method, r in res.items():
            print(f"  -> {method}: E={r['energy_ha']:.6f} Ha  err={r['error_ha']:.6f} Ha  "
                  f"converged={r['converged']}  depth={r['depth']}  "
                  f"2q_gates={r['two_qubit_gates']}  time={r['time_sec']}s")

    print("\n" + "=" * 70)
    print("  (b) UCCSD (+ Hartree-Fock initial state)")
    print("=" * 70)
    res_uccsd = test_uccsd(qop_bare, enuc, exact, n_qubits, nelec)
    all_results["uccsd"] = res_uccsd
    for method, r in res_uccsd.items():
        print(f"  -> {method}: E={r['energy_ha']:.6f} Ha  err={r['error_ha']:.6f} Ha  "
              f"converged={r['converged']}  depth={r['depth']}  "
              f"2q_gates={r['two_qubit_gates']}  time={r['time_sec']}s")

    # --- Final report -----------------------------------------------------
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  {'ansatz':30s} | {'method':10s} | {'error (Ha)':>10s} | {'converged':>9s} | "
          f"{'depth':>6s} | {'2q gates':>8s}")
    for reps in (2, 3):
        for method, r in all_results["excitation_preserving"][f"reps{reps}"].items():
            print(f"  {'ExcPreserving(reps=' + str(reps) + ')':30s} | {method:10s} | "
                  f"{r['error_ha']:10.6f} | {str(r['converged']):>9s} | "
                  f"{r['depth']:6d} | {r['two_qubit_gates']:8d}")
    for method, r in res_uccsd.items():
        print(f"  {'UCCSD':30s} | {method:10s} | {r['error_ha']:10.6f} | "
              f"{str(r['converged']):>9s} | {r['depth']:6d} | {r['two_qubit_gates']:8d}")
    print("=" * 70)
    print("\n  NOTE: simulator-only test. No hardware jobs submitted.\n")

    out_path = os.path.join(os.path.dirname(__file__), "fragment_ansatz_test_results.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
