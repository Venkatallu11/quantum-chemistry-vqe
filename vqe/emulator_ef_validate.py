#!/usr/bin/env python3
"""
emulator_ef_validate.py — validate the local Quantinuum-like noise model
used in entanglement_forging_h4.py / entanglement_forging_zne.py against
a REAL Quantinuum emulator job.
============================================================================
Takes ONE real circuit from the EF pipeline: state-preparation of the
alpha register's leading Schmidt vector (u_vecs[0], 4 qubits), then
measures a diagonal Pauli observable -- the same kind of measurement
build_noisy_matrices() does locally via AerEstimatorV2.

Label choice: "ZZZZ" is NOT actually one of this Hamiltonian's Pauli
terms (checked directly against decompose_pauli_terms() output) -- using
it would validate a term that isn't part of the real energy expression.
"ZIII" is the real diagonal alpha-register term with the largest
coefficient, so that's what's measured here.

Runs the SAME observable/state three ways:
  1. Exact (numpy)                                    -- ground truth
  2. Local Aer, Quantinuum-like noise model (scale=1)  -- what this repo's
     EF/ZNE results have been assuming Quantinuum-grade noise looks like
  3. Real quantinuum.sim.h2-1e emulator (Azure Quantum) -- what Quantinuum's
     own emulator actually predicts

If (2) and (3) disagree substantially, the local noise model used
throughout entanglement_forging_h4.py / entanglement_forging_zne.py is
not a trustworthy stand-in for real Quantinuum-grade noise.

Run:
    python vqe/emulator_ef_validate.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
from azure_backend import connect_workspace
from azure.quantum.qiskit import AzureQuantumProvider
from qiskit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Pauli
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2

LABEL = "ZIII"
TARGET = "quantinuum.sim.h2-1e"
SHOTS = 500


def zlabel_expectation_from_counts(counts, label, shots):
    """<P> for a diagonal I/Z-only label from measured Z-basis counts:
    sum_b (-1)^(parity of bits where label is Z) * count(b) / shots."""
    z_positions = [i for i, c in enumerate(reversed(label)) if c == "Z"]
    total = 0.0
    for bitstring, n in counts.items():
        bits = bitstring.replace(" ", "")
        parity = sum(int(bits[len(bits) - 1 - i]) for i in z_positions) % 2
        total += ((-1) ** parity) * n
    return total / shots


def main():
    print("\n" + "=" * 70)
    print("  Validating the local Quantinuum-like noise model against a")
    print(f"  REAL Quantinuum emulator job -- single EF circuit, {LABEL}")
    print("=" * 70)

    qop_bare, qop_penalized, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_penalized)
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)

    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = set(a for a, _, _ in terms)
    if LABEL not in alpha_labels:
        print(f"\n  {LABEL} is not one of this Hamiltonian's alpha-register "
              f"Pauli terms -- pick a real one from decompose_pauli_terms().")
        sys.exit(1)

    u0 = u_vecs[0]  # leading Schmidt vector, alpha register (4 qubits)

    # 1. Exact (numpy)
    Pz = Pauli(LABEL).to_matrix()
    exact_val = float(np.real(u0.conj() @ Pz @ u0))
    print(f"\n  1. Exact <psi_0|{LABEL}|psi_0>              = {exact_val:.6f}")

    # 2. Local Aer, Quantinuum-like noise
    circ_t = ef.transpiled_state_prep_circuit(u0)
    nm = ef.build_noise_model()
    aer_est = AerEstimatorV2(
        options={"backend_options": {"noise_model": nm, "method": "density_matrix"}}
    )
    pub = (circ_t, [Pauli(LABEL)])
    aer_val = float(aer_est.run([pub]).result()[0].data.evs[0])
    print(f"  2. Local Aer (Quantinuum-like noise model)   = {aer_val:.6f}  "
          f"(diff from exact: {abs(aer_val - exact_val):.6f})")

    # 3. Real Quantinuum emulator
    print(f"\n  Connecting to Azure Quantum, submitting to {TARGET} ...")
    ws = connect_workspace()
    backend = AzureQuantumProvider(ws).get_backend(TARGET)

    meas_circ = QuantumCircuit(4, 4)
    meas_circ.append(StatePreparation(np.asarray(u0)), range(4))
    meas_circ.measure(range(4), range(4))

    job = backend.run(meas_circ, shots=SHOTS)
    print(f"  job id: {job.id()} -- waiting ...")
    counts = job.result().get_counts()
    real_val = zlabel_expectation_from_counts(counts, LABEL, SHOTS)
    print(f"  3. Real {TARGET} emulator                    = {real_val:.6f}  "
          f"(diff from exact: {abs(real_val - exact_val):.6f})")

    print("\n" + "-" * 70)
    print(f"  Local-model vs real-emulator gap: {abs(aer_val - real_val):.6f}")
    print("-" * 70)

    results = {
        "label": LABEL,
        "exact": round(exact_val, 6),
        "local_aer_quantinuum_like": round(aer_val, 6),
        "local_aer_error_vs_exact": round(abs(aer_val - exact_val), 6),
        "real_quantinuum_emulator": round(real_val, 6),
        "real_emulator_error_vs_exact": round(abs(real_val - exact_val), 6),
        "local_vs_real_gap": round(abs(aer_val - real_val), 6),
        "shots": SHOTS,
        "job_id": job.id(),
    }
    out = os.path.join(os.path.dirname(__file__), "emulator_ef_validate_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
