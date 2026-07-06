#!/usr/bin/env python3
"""
hardware_fragmentation.py — covalent-bond fragmentation (molecular
tailoring), measured on REAL IBM hardware, for a bonded H6 chain.

Same idea as covalent_fragment.py (overlapping fragments + inclusion-
exclusion), but this time each fragment's optimized VQE circuit is
measured on ACTUAL hardware instead of a simulator:

    E(H6) = E_hw(atoms 0-3) + E_hw(atoms 2-5) - E_hw(atoms 2-3)

Pipeline per fragment:
  1. chem.py integrals -> RHF -> qiskit-nature qubit Hamiltonian
     (exact same pipeline as molecules_real.py / covalent_fragment.py).
  2. Optimize an EfficientSU2 ansatz on the LOCAL statevector simulator
     until it matches that fragment's EXACT (eigsh) ground state within
     0.002 Ha -- this is a real correctness gate, not a formality: if a
     fragment can't get there even after escalating reps/restarts, this
     script stops and says which fragment failed, rather than silently
     measuring a wrong circuit on hardware.
  3. Measure that ONE optimized circuit ONCE on real hardware (no loop,
     no re-optimization on hardware -- classical parameters, one real
     measurement, same as h2_hardware_direct.py).

Three fragments -> three hardware jobs total, no more.

Run:
    python vqe/hardware_fragmentation.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import eigsh

from qiskit.circuit.library import EfficientSU2
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import StatevectorEstimator
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

sys.path.insert(0, os.path.dirname(__file__))
import chem
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper

from covalent_fragment import hchain               # reuse the same geometry builder
from molecules_real import _constrain_particle_number  # reuse the same penalty trick

from dotenv import load_dotenv
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "..", ".env"))

HARTREE_TO_KCAL_MOL = 627.5094740631
CHEM_ACCURACY_TARGET_HA = 0.002  # how close local VQE must get to exact before trusting it

# Reference numbers from covalent_fragment.py, for comparison at the end.
FULL_EXACT_H6_HA = -3.236066
SIM_TAILORING_H6_HA = -3.231625


# --------------------------------------------------------------------------
# Building each fragment's Hamiltonian
# --------------------------------------------------------------------------

def _build_fragment_qop(geom, nelec):
    """
    Build the qubit Hamiltonian for one fragment, using the identical
    integral -> RHF -> qiskit-nature pipeline as molecules_real.py.

    Returns:
        qop_bare:      plain electronic SparsePauliOp (no penalty term).
                       This is what gets MEASURED on hardware and is the
                       physically meaningful electronic energy.
        qop_penalized: qop_bare + a particle-number penalty. Used ONLY to
                       steer classical VQE and exact diagonalization into
                       the correct-electron-count subspace -- never
                       measured on hardware (its penalty term would just
                       amplify noise sensitivity for no physical reason).
        enuc: nuclear repulsion (classical constant). Add this to the
              electronic expectation value to get the fragment's total
              energy.
    """
    S, T, V, eri, enuc = chem.integrals(geom)
    _ehf, C, Hc = chem.rhf(S, T, V, eri, enuc, nelec=nelec)
    h1 = C.T @ Hc @ C
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)
    ee = ElectronicEnergy.from_raw_integrals(h1, h2)
    qop_bare = JordanWignerMapper().map(ee.second_q_op())
    qop_penalized = _constrain_particle_number(qop_bare, nelec)
    return qop_bare, qop_penalized, enuc


def _exact_energy(qop_penalized, enuc):
    """Exact ground state via sparse diagonalization -- the reference answer."""
    M = qop_penalized.to_matrix(sparse=True)
    e_elec = float(eigsh(M, k=1, which="SA", return_eigenvectors=False)[0])
    return e_elec + enuc


# --------------------------------------------------------------------------
# Local (simulator) VQE optimization, with escalation if it doesn't converge
# --------------------------------------------------------------------------

def _optimize_ansatz(qop_penalized, enuc, exact_energy, n_qubits, reps, n_restarts, maxiter):
    """
    Optimize an EfficientSU2(n_qubits, reps) ansatz on the noiseless local
    statevector simulator, trying `n_restarts` different starting points.

    Back to COBYLA (same optimizer as h2_vqe.py / mcp_energy.py) -- an
    L-BFGS-B + hand-built Hartree-Fock warm-start attempt was tried and
    abandoned: EfficientSU2's entangling CX gates are unconditional (they
    fire regardless of rotation angles), so a bitstring carefully prepared
    in the first rotation layer gets scrambled by the fixed entanglers
    before the optimizer even takes a step -- that produced a nonsense
    result (E = 82.9 Ha), not a real one. COBYLA doesn't need that kind of
    careful state construction. "Enhanced" here means: more restarts and
    much higher maxiter than the original attempt (COBYLA's own earlier
    best was error 0.288 Ha at maxiter=2000/3 restarts -- this pushes both
    much further before concluding it can't converge).

    Returns (best_params, best_energy, best_error).
    """
    estimator = StatevectorEstimator()
    ansatz = EfficientSU2(n_qubits, reps=reps)
    n_params = ansatz.num_parameters

    def objective(params):
        qc = ansatz.assign_parameters(params)
        result = estimator.run([(qc, qop_penalized)]).result()
        return float(result[0].data.evs) + enuc

    rng = np.random.default_rng(42)
    best_params, best_energy = None, None
    for r in range(n_restarts):
        x0 = np.full(n_params, 0.1) if r == 0 else rng.uniform(-np.pi, np.pi, size=n_params)
        result = minimize(
            objective, x0=x0, method="COBYLA",
            options={"maxiter": maxiter, "rhobeg": 0.5, "tol": 1e-10},
        )
        if best_energy is None or result.fun < best_energy:
            best_energy, best_params = float(result.fun), result.x

        error = abs(best_energy - exact_energy)
        print(f"    restart {r + 1}/{n_restarts}: E = {result.fun:.6f} Ha "
              f"(best so far: {best_energy:.6f} Ha, error {error:.6f} Ha)")
        if error < CHEM_ACCURACY_TARGET_HA:
            break

    return best_params, best_energy, abs(best_energy - exact_energy)


def process_fragment(name, indices, d):
    """
    Full local pipeline for one fragment: build its Hamiltonian, then
    optimize + verify VQE reaches the exact ground state, escalating
    reps/restarts if needed. Raises RuntimeError (stop-and-report) if it
    still can't get there.
    """
    geom = hchain(indices, d)
    nelec = len(indices)
    print(f"\n--- Fragment '{name}' (atoms {indices}, {nelec} electrons) ---")

    qop_bare, qop_penalized, enuc = _build_fragment_qop(geom, nelec)
    n_qubits = qop_bare.num_qubits
    exact = _exact_energy(qop_penalized, enuc)
    print(f"  qubits = {n_qubits}, exact ground state = {exact:.6f} Ha")

    # COBYLA (same optimizer as h2_vqe.py / mcp_energy.py), enhanced with
    # many more restarts and a much higher maxiter than the original
    # attempt (which used maxiter=2000/3 restarts and got to 0.288 Ha).
    reps, restarts, maxiter = 2, 5, 8000
    max_attempts = 3  # escalate restarts/maxiter this many times before giving up
    for attempt in range(1, max_attempts + 1):
        print(f"  Optimizing EfficientSU2(reps={reps}) with COBYLA, {restarts} restarts, "
              f"maxiter={maxiter} (attempt {attempt}/{max_attempts})...")
        params, sim_energy, sim_error = _optimize_ansatz(
            qop_penalized, enuc, exact, n_qubits, reps=reps, n_restarts=restarts, maxiter=maxiter,
        )
        print(f"  -> local VQE best: E = {sim_energy:.6f} Ha, error = {sim_error:.6f} Ha")

        if sim_error < CHEM_ACCURACY_TARGET_HA:
            break
        if attempt == max_attempts:
            raise RuntimeError(
                f"Fragment '{name}' could NOT reach the exact ground state within "
                f"{CHEM_ACCURACY_TARGET_HA} Ha even after {max_attempts} escalation attempts "
                f"with COBYLA (best error achieved: {sim_error:.6f} Ha). STOPPING."
            )
        print(f"  Not converged -- escalating: restarts {restarts}->{restarts + 3}, "
              f"maxiter {maxiter}->{maxiter + 4000}")
        restarts += 3
        maxiter += 4000

    return {
        "name": name, "indices": indices, "nelec": nelec, "n_qubits": n_qubits,
        "reps": reps, "exact_energy_ha": exact, "sim_energy_ha": sim_energy,
        "sim_error_ha": sim_error, "params": params, "qop_bare": qop_bare, "enuc": enuc,
    }


# --------------------------------------------------------------------------
# Real hardware measurement
# --------------------------------------------------------------------------

def _pick_least_busy_backend(service):
    """Check every accessible backend's queue depth directly and return the least-busy operational one."""
    candidates = []
    for b in service.backends():
        try:
            s = b.status()
            if s.operational:
                candidates.append((b, s.pending_jobs))
        except Exception:
            pass
    backend, queue_depth = min(candidates, key=lambda x: x[1])
    print(f"  Least-busy device right now: {backend.name} ({queue_depth} pending jobs)")
    return backend


def measure_on_hardware(fragment, backend):
    """
    Submit ONE real hardware job for this fragment's already-optimized
    ansatz, measuring the BARE (non-penalized) electronic Hamiltonian,
    then add back the nuclear repulsion to report the total fragment energy.
    """
    ansatz = EfficientSU2(fragment["n_qubits"], reps=fragment["reps"])
    qc = ansatz.assign_parameters(fragment["params"])

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)
    isa_observable = fragment["qop_bare"].apply_layout(isa_circuit.layout)

    estimator = EstimatorV2(mode=backend)
    estimator.options.resilience_level = 1
    estimator.options.default_shots = 4096

    job = estimator.run([(isa_circuit, isa_observable)])
    job_id = job.job_id()
    print(f"  Submitted job {job_id} on {backend.name} for fragment '{fragment['name']}'")

    result = job.result()
    electronic_energy = float(result[0].data.evs)
    hardware_energy = electronic_energy + fragment["enuc"]
    return hardware_energy, job_id


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    d = 1.0  # Angstrom, same bonded H6 chain as covalent_fragment.py

    print("\n" + "=" * 70)
    print("  Covalent-bond fragmentation on REAL hardware -- H6 chain")
    print("=" * 70)

    # --- Step 1: build + locally verify each fragment -------------------
    fragments = {
        "block1":  process_fragment("block1 (atoms 0-3)",  [0, 1, 2, 3], d),
        "block2":  process_fragment("block2 (atoms 2-5)",  [2, 3, 4, 5], d),
        "overlap": process_fragment("overlap (atoms 2-3)", [2, 3],       d),
    }

    # --- Step 2: pick ONE backend for all 3 hardware measurements -------
    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        sys.exit(1)
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = _pick_least_busy_backend(service)

    # --- Step 3: measure each fragment ONCE on real hardware ------------
    print("\n" + "=" * 70)
    print("  Hardware measurements (3 jobs total)")
    print("=" * 70)
    for frag in fragments.values():
        hw_energy, job_id = measure_on_hardware(frag, backend)
        frag["hardware_energy_ha"] = hw_energy
        frag["job_id"] = job_id
        err = abs(hw_energy - frag["exact_energy_ha"])
        print(f"  {frag['name']}: hardware = {hw_energy:.6f} Ha, "
              f"exact = {frag['exact_energy_ha']:.6f} Ha, error = {err:.6f} Ha")

    # --- Step 4: reassemble via molecular tailoring ----------------------
    e_hw_h6 = (fragments["block1"]["hardware_energy_ha"]
               + fragments["block2"]["hardware_energy_ha"]
               - fragments["overlap"]["hardware_energy_ha"])

    err_vs_exact = abs(e_hw_h6 - FULL_EXACT_H6_HA)
    err_vs_sim_tailor = abs(e_hw_h6 - SIM_TAILORING_H6_HA)

    # --- Summary table ----------------------------------------------------
    print("\n" + "=" * 70)
    print("  FRAGMENT SUMMARY")
    print("=" * 70)
    print(f"  {'fragment':22s} | {'qubits':>6s} | {'exact (Ha)':>12s} | {'hardware (Ha)':>14s} | {'error (Ha)':>10s}")
    for frag in fragments.values():
        err = abs(frag["hardware_energy_ha"] - frag["exact_energy_ha"])
        print(f"  {frag['name']:22s} | {frag['n_qubits']:6d} | {frag['exact_energy_ha']:12.6f} | "
              f"{frag['hardware_energy_ha']:14.6f} | {err:10.6f}")

    print("\n" + "=" * 70)
    print("  REASSEMBLED H6 (molecular tailoring on REAL hardware)")
    print("=" * 70)
    print("  E_hardware(H6) = E_hw(block1) + E_hw(block2) - E_hw(overlap)")
    print(f"                 = {fragments['block1']['hardware_energy_ha']:.6f} "
          f"+ {fragments['block2']['hardware_energy_ha']:.6f} "
          f"- {fragments['overlap']['hardware_energy_ha']:.6f}")
    print(f"                 = {e_hw_h6:.6f} Ha\n")
    print(f"  Full exact H6           = {FULL_EXACT_H6_HA:.6f} Ha")
    print(f"  Simulator tailoring H6  = {SIM_TAILORING_H6_HA:.6f} Ha")
    print(f"  Hardware tailoring H6   = {e_hw_h6:.6f} Ha\n")
    print(f"  Error vs full exact     = {err_vs_exact:.6f} Ha  ({err_vs_exact * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print(f"  Error vs sim tailoring  = {err_vs_sim_tailor:.6f} Ha  ({err_vs_sim_tailor * HARTREE_TO_KCAL_MOL:.3f} kcal/mol)")
    print("=" * 70)

    # --- Save results -------------------------------------------------------
    results = {
        "backend": backend.name,
        "fragments": {
            key: {
                "name": frag["name"],
                "indices": frag["indices"],
                "n_qubits": frag["n_qubits"],
                "reps": frag["reps"],
                "exact_energy_ha": round(frag["exact_energy_ha"], 6),
                "sim_vqe_energy_ha": round(frag["sim_energy_ha"], 6),
                "sim_vqe_error_ha": round(frag["sim_error_ha"], 6),
                "hardware_energy_ha": round(frag["hardware_energy_ha"], 6),
                "hardware_error_ha": round(abs(frag["hardware_energy_ha"] - frag["exact_energy_ha"]), 6),
                "job_id": frag["job_id"],
            }
            for key, frag in fragments.items()
        },
        "reassembled_hardware_h6_ha": round(e_hw_h6, 6),
        "full_exact_h6_ha": FULL_EXACT_H6_HA,
        "sim_tailoring_h6_ha": SIM_TAILORING_H6_HA,
        "error_vs_full_exact_ha": round(err_vs_exact, 6),
        "error_vs_full_exact_kcal_mol": round(err_vs_exact * HARTREE_TO_KCAL_MOL, 4),
        "error_vs_sim_tailoring_ha": round(err_vs_sim_tailor, 6),
        "error_vs_sim_tailoring_kcal_mol": round(err_vs_sim_tailor * HARTREE_TO_KCAL_MOL, 4),
    }
    out_path = os.path.join(os.path.dirname(__file__), "hardware_fragmentation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
