#!/usr/bin/env python3
"""
h2_vqe.py  —  Ground State Energy of H₂ via VQE
=================================================
The first real step toward simulating the FeMoco enzyme
(the molecule that could solve nitrogen fixation).

What this does
--------------
VQE (Variational Quantum Eigensolver) finds the lowest energy a molecule
can have — its "ground state." This matters because:
  - Ground state energy tells you if a chemical reaction will happen
  - It tells you how molecules bond and break apart
  - Classical computers can only approximate it for big molecules
  - Quantum computers can calculate it EXACTLY

H₂ only needs 2 qubits. We start here, then scale up:
  H₂ (2 qubits) → LiH (4 qubits) → H₂O (8 qubits) → FeMoco (54 qubits)

How to run
----------
  # Local simulation (fast, no IBM token needed, always try this first)
  python vqe/h2_vqe.py

  # Real IBM quantum hardware (needs IBM_QUANTUM_TOKEN in .env)
  python vqe/h2_vqe.py --real

  # Real hardware on a specific machine
  python vqe/h2_vqe.py --real --device ibm_brisbane
"""

import os
import sys
import json
import argparse
import numpy as np
from scipy.optimize import minimize

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import StatevectorEstimator

# Load .env from project root (same token the MCP server uses)
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from dotenv import load_dotenv
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "..", ".env"))


# --------------------------------------------------------------------------
# The H₂ Hamiltonian
# --------------------------------------------------------------------------
# To simulate a molecule on a quantum computer, we first convert it into
# a "Hamiltonian" — a mathematical description of its energy.
#
# For H₂ at its natural bond length (0.735 Å), in the STO-3G basis set,
# using Jordan-Wigner mapping, we get this 2-qubit Hamiltonian.
# Each key is a Pauli operator; each value is its energy coefficient (Hartree).
#
# Think of it like: E = a*(IZ measurement) + b*(ZI measurement) + ...
# The quantum computer measures each term; we multiply by coefficients and sum.
#
# The XX and YY terms are the interesting ones — they represent quantum
# correlations between electrons that classical computers can't handle at scale.
H2_HAMILTONIAN = {
    "II":  -0.8105479,   # constant offset (nuclear repulsion + mean field)
    "IZ":   0.1712076,   # energy of electron in orbital 0
    "ZI":  -0.2257534,   # energy of electron in orbital 1
    "ZZ":   0.1209057,   # electron-electron repulsion
    "XX":   0.1686013,   # quantum correlation ← classically hard at scale
    "YY":   0.1686013,   # quantum correlation ← classically hard at scale
}

# How close to the exact answer is "good enough" for chemistry?
# 1 kcal/mol = 0.0016 Hartree. If we're within this, we can predict
# whether reactions happen, how drugs bind, how enzymes work.
CHEMICAL_ACCURACY_HA = 0.0016  # Hartree


def exact_ground_state_energy() -> float:
    """
    Compute the exact ground state energy by diagonalizing the Hamiltonian.
    This is what VQE is trying to find. We use this as the reference answer.

    Classical computers can do this for small molecules like H₂ —
    but the cost grows exponentially with molecule size. For FeMoco
    (54 qubits), this classical diagonalization is impossible.
    """
    # Build the full Hamiltonian matrix from our Pauli terms
    hamiltonian_op = SparsePauliOp.from_list(list(H2_HAMILTONIAN.items()))
    matrix = hamiltonian_op.to_matrix()

    # eigvalsh finds all eigenvalues of a Hermitian matrix.
    # The smallest eigenvalue IS the ground state energy — this is a
    # fundamental result from quantum mechanics (variational principle).
    eigenvalues = np.linalg.eigvalsh(matrix)
    return float(eigenvalues[0])


# --------------------------------------------------------------------------
# The Ansatz Circuit
# --------------------------------------------------------------------------

def make_ansatz(theta: float) -> QuantumCircuit:
    """
    Build the particle-conserving ansatz for H₂.

    For H₂ with 2 electrons, the ground state is always a superposition of
    exactly two states: |01⟩ and |10⟩ (one electron in each spin orbital).
    The ansatz must stay in that subspace — this is the "particle conservation"
    constraint that mirrors real physics (electrons don't disappear).

    The circuit: X(q0) → RY(θ, q1) → CX(q1→q0)
    produces:    cos(θ/2)|01⟩ + sin(θ/2)|10⟩

    At the optimal θ, this IS the exact ground state — no approximation.

    Args:
        theta: Single rotation angle (radians). The optimizer tunes this.

    Returns:
        A 2-qubit circuit. No measurements — the Estimator handles that.
    """
    qc = QuantumCircuit(2)

    # Start in the Hartree-Fock state |01⟩ (qubit 0 has an electron, qubit 1 doesn't)
    # This is the classical "best guess" ground state — VQE improves on it.
    qc.x(0)

    # Rotate qubit 1 — this creates a superposition of the two electron configs:
    #   θ=0  → stays in |01⟩ (pure Hartree-Fock, no quantum correlation)
    #   optimal θ → the true quantum ground state (electrons correlated)
    qc.ry(theta, 1)

    # Entangle with qubit 0: when qubit 1 rotates to |1⟩, qubit 0 flips to |0⟩.
    # This keeps the total electron count exactly 1 — physically correct.
    # Result: cos(θ/2)|01⟩ + sin(θ/2)|10⟩
    qc.cx(1, 0)

    return qc


# --------------------------------------------------------------------------
# Energy measurement
# --------------------------------------------------------------------------

def measure_energy(theta: float, estimator) -> float:
    """
    Compute ⟨ψ(θ)|H|ψ(θ)⟩ — the energy of our trial state.

    This is the core of VQE:
    1. Build the circuit for this angle
    2. For each Pauli term in H, run the circuit and measure that observable
    3. Weight each measurement by its Hamiltonian coefficient
    4. Sum everything up → total energy of the molecule at this angle

    Args:
        theta: Current rotation angle
        estimator: Quantum estimator (simulator or real IBM hardware)

    Returns:
        Energy in Hartree. The optimizer will try to make this as small as possible.
    """
    qc = make_ansatz(theta)
    energy = 0.0

    for pauli_str, coefficient in H2_HAMILTONIAN.items():
        # "II" = identity on both qubits = just a constant, no measurement needed
        if pauli_str == "II":
            energy += coefficient
            continue

        # For all other Pauli terms: measure the expectation value on our circuit.
        # The estimator runs the circuit and returns a number in [-1, +1].
        observable = SparsePauliOp(pauli_str)
        result = estimator.run([(qc, observable)]).result()
        expectation = float(result[0].data.evs)

        # Add this term's contribution to the total energy
        energy += coefficient * expectation

    return energy


# --------------------------------------------------------------------------
# VQE on local simulator
# --------------------------------------------------------------------------

def run_vqe_local() -> dict:
    """
    Run VQE using a local statevector simulator.

    The statevector simulator computes quantum mechanics EXACTLY (no noise).
    It's limited to ~30 qubits on a normal laptop, but it's instant and
    perfect for verifying our approach before using real hardware.
    """
    print("\n" + "="*60)
    print("  VQE — Local Statevector Simulator")
    print("  (Exact quantum simulation on your CPU)")
    print("="*60 + "\n")

    exact = exact_ground_state_energy()
    print(f"  Target (exact classical answer): {exact:.6f} Hartree\n")

    estimator = StatevectorEstimator()
    iteration_log = []

    def objective(params):
        theta = params[0]
        energy = measure_energy(theta, estimator)
        n = len(iteration_log) + 1
        iteration_log.append({
            "iteration": n,
            "theta": round(float(theta), 5),
            "energy": round(float(energy), 8),
        })
        # Print progress every 10 iterations
        if n % 10 == 0 or n == 1:
            gap = abs(energy - exact)
            marker = " ✓" if gap < CHEMICAL_ACCURACY_HA else ""
            print(f"  Iter {n:3d} | θ = {theta:+.4f} | E = {energy:.6f} Ha | gap = {gap:.6f} Ha{marker}")
        return energy

    # COBYLA: gradient-free optimizer — works well with quantum noise too.
    # We start near the Hartree-Fock state (theta=0) and let it explore.
    result = minimize(
        objective,
        x0=[0.1],
        method="COBYLA",
        options={"maxiter": 300, "rhobeg": 0.5},
    )

    final_energy = result.fun
    error = abs(final_energy - exact)

    return {
        "method": "local_statevector_simulator",
        "optimal_theta": float(result.x[0]),
        "vqe_energy_ha": round(final_energy, 8),
        "exact_energy_ha": round(exact, 8),
        "error_ha": round(error, 8),
        "chemical_accuracy_threshold_ha": CHEMICAL_ACCURACY_HA,
        "reached_chemical_accuracy": bool(error < CHEMICAL_ACCURACY_HA),
        "total_iterations": len(iteration_log),
        "optimizer_converged": bool(result.success),
        "history": iteration_log,
    }


# --------------------------------------------------------------------------
# VQE on real IBM quantum hardware
# --------------------------------------------------------------------------

def run_vqe_real_hardware(device_name: str | None = None) -> dict:
    """
    Run VQE on real IBM quantum hardware via the MCP server's connection.

    On real hardware:
    - There is quantum noise (gates aren't perfect)
    - Each energy evaluation = a real job submitted to IBM's queue
    - More shots per measurement = less statistical noise, but slower
    - We use fewer optimizer iterations than the simulator (to save queue time)

    This is the real deal — actual electrons, actual quantum mechanics.
    """
    from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2 as RealEstimator
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    print("\n" + "="*60)
    print("  VQE — Real IBM Quantum Hardware")
    print("  (Actual quantum computation — may take time in queue)")
    print("="*60 + "\n")

    token = os.getenv("IBM_QUANTUM_TOKEN")
    if not token:
        print("ERROR: IBM_QUANTUM_TOKEN not found in .env")
        print("Get your token at https://quantum.ibm.com/account")
        sys.exit(1)

    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)

    # Auto-select the least-busy operational backend if none specified
    if not device_name:
        print("  Finding least-busy IBM quantum computer...")
        backends = service.backends()
        operational = []
        for b in backends:
            try:
                s = b.status()
                if s.operational:
                    operational.append((b, s.pending_jobs))
            except Exception:
                pass
        if not operational:
            print("ERROR: No operational backends found.")
            sys.exit(1)
        best, queue_depth = min(operational, key=lambda x: x[1])
        device_name = best.name
        print(f"  Selected: {device_name} ({queue_depth} jobs in queue)\n")

    backend = service.backend(device_name)
    exact = exact_ground_state_energy()
    print(f"  Target (exact classical answer): {exact:.6f} Hartree\n")

    # Transpile once — compile our circuit to this backend's native gates.
    # We reuse the transpiled circuit for all optimizer iterations to save time.
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)

    estimator = RealEstimator(mode=backend)
    # 2048 shots: more measurements = less noise in each energy estimate.
    # 1024 is the default; we use more for better accuracy.
    estimator.options.default_shots = 2048

    n_backend_qubits = backend.num_qubits
    iteration_log = []

    def objective(params):
        theta = params[0]
        qc = make_ansatz(theta)
        isa_circuit = pm.run(qc)  # transpile to this backend's native gate set

        energy = 0.0
        pubs = []
        coefficients = []

        for pauli_str, coefficient in H2_HAMILTONIAN.items():
            if pauli_str == "II":
                energy += coefficient
                continue
            # IBM backends have many qubits; our observable must match the
            # full transpiled circuit width. Pad with identity qubits on the left.
            padded = "I" * (n_backend_qubits - len(pauli_str)) + pauli_str
            pubs.append((isa_circuit, SparsePauliOp(padded)))
            coefficients.append(coefficient)

        # Submit all 5 observables as one job — efficient use of queue time
        result = estimator.run(pubs).result()
        for i, coeff in enumerate(coefficients):
            energy += coeff * float(result[i].data.evs)

        n = len(iteration_log) + 1
        iteration_log.append({
            "iteration": n,
            "theta": round(float(theta), 5),
            "energy": round(float(energy), 8),
        })
        gap = abs(energy - exact)
        marker = " ✓" if gap < CHEMICAL_ACCURACY_HA else ""
        print(f"  Iter {n:3d} | θ = {theta:+.4f} | E = {energy:.6f} Ha | gap = {gap:.6f} Ha{marker}")
        return energy

    # Fewer iterations on real hardware — each one is a real quantum job
    result = minimize(
        objective,
        x0=[0.1],
        method="COBYLA",
        options={"maxiter": 50, "rhobeg": 0.5},
    )

    final_energy = result.fun
    error = abs(final_energy - exact)

    return {
        "method": "real_ibm_hardware",
        "device": device_name,
        "optimal_theta": float(result.x[0]),
        "vqe_energy_ha": round(final_energy, 8),
        "exact_energy_ha": round(exact, 8),
        "error_ha": round(error, 8),
        "chemical_accuracy_threshold_ha": CHEMICAL_ACCURACY_HA,
        "reached_chemical_accuracy": bool(error < CHEMICAL_ACCURACY_HA),
        "total_iterations": len(iteration_log),
        "history": iteration_log,
    }


# --------------------------------------------------------------------------
# Pretty results summary
# --------------------------------------------------------------------------

def print_results(r: dict) -> None:
    print("\n" + "="*60)
    print("  RESULTS")
    print("="*60)
    print(f"  Method:          {r['method']}")
    if r.get("device"):
        print(f"  Device:          {r['device']}")
    print()
    print(f"  VQE energy:      {r['vqe_energy_ha']:.6f} Hartree")
    print(f"  Exact energy:    {r['exact_energy_ha']:.6f} Hartree")
    print(f"  Error:           {r['error_ha']:.6f} Hartree")
    print(f"  Threshold:       {r['chemical_accuracy_threshold_ha']} Hartree  (1 kcal/mol)")
    print()

    if r["reached_chemical_accuracy"]:
        print("  ✓ CHEMICAL ACCURACY REACHED")
        print("  ✓ This result is accurate enough to predict real chemistry!")
    else:
        times_off = r["error_ha"] / r["chemical_accuracy_threshold_ha"]
        print(f"  ✗ {times_off:.1f}x above chemical accuracy")
        if r.get("device"):
            print("    Hardware noise is likely the cause.")
            print("    Try --device with a lower-error backend (use compare_devices).")
        else:
            print("    Try increasing maxiter in the optimizer.")

    print()
    print("  Why H₂ matters:")
    print("  ─────────────────────────────────────────────────────")
    print("  H₂ = 2 qubits. We solved it.")
    print("  LiH = 4 qubits.  Next step.")
    print("  H₂O = 8 qubits.  The step after.")
    print("  FeMoco (nitrogen fixation enzyme) = ~54 qubits.")
    print("  IBM hardware today: 127–433 qubits.")
    print("  The door is open.")
    print("="*60)


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="VQE for H₂ — ground state energy on simulator or real IBM hardware",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python vqe/h2_vqe.py                        # local simulation (recommended first)
  python vqe/h2_vqe.py --real                 # real IBM hardware (auto-selects best machine)
  python vqe/h2_vqe.py --real --device ibm_fez  # real hardware, specific machine
        """,
    )
    parser.add_argument(
        "--real",
        action="store_true",
        help="Run on real IBM quantum hardware instead of local simulator",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="IBM backend name (e.g. ibm_brisbane). Auto-selects least-busy if omitted.",
    )
    args = parser.parse_args()

    if args.real:
        results = run_vqe_real_hardware(device_name=args.device)
    else:
        results = run_vqe_local()

    print_results(results)

    # Save full results (including iteration history) to JSON
    out_path = os.path.join(os.path.dirname(__file__), "h2_vqe_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Full results saved → {out_path}\n")
