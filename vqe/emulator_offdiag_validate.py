#!/usr/bin/env python3
"""
emulator_offdiag_validate.py — validate the local Quantinuum-like noise
model on an OFF-DIAGONAL (X/Y-containing) Hamiltonian term, via the EF
cross-term superposition trick, against a REAL Quantinuum emulator.
============================================================================
Every prior validation script (emulator_ef_validate.py) only checked
DIAGONAL (I/Z-only) terms, because those can be read off directly from a
single computational-basis measurement. Off-diagonal terms like "XXII"
appear in the EF energy reconstruction as CROSS terms <u_n|P|u_m> (n!=m)
between different Schmidt vectors, and need the superposition trick from
entanglement_forging_h4.py's measure_cross_noisy(): prepare
    phi_k = (u_n + i^k * u_m) / sqrt(2),  k = 0,1,2,3
measure <phi_k|P|phi_k> for each, then
    Re(<u_n|P|u_m>) = (E0 - E2) / 2
    Im(<u_n|P|u_m>) = (E3 - E1) / 2

Label + pair: the first version of this script used "XXII" (the
largest-*Hamiltonian-coefficient* off-diagonal term), but its actual
cross-term magnitude |<u_0|XXII|u_1>| turned out to be tiny (~0.0029) --
too small to resolve against ordinary shot noise at 1000 shots/phase,
making that run statistically inconclusive. Hamiltonian-coefficient size
and cross-term magnitude are NOT the same thing. Checked directly (not
assumed): among all off-diagonal alpha labels, "IXXI" has by far the
largest actual cross-term magnitude between the two leading Schmidt
vectors, |<u_0|IXXI|u_1>| = 0.961330 -- a strong, easily resolvable
signal at 1000 shots. (n, m) = (0, 1): the two leading Schmidt vectors.

REAL-HARDWARE CONSTRAINT: Azure's Quantinuum backend only supports ONE
circuit per job (confirmed directly against the SDK source --
AzureQirBackend.run() raises NotImplementedError for a circuit list of
length != 1). Unlike the diagonal validation (one job, many labels
derived from one shot record), the 4 cross-term phase circuits each need
their own real job -- there is no way to batch them into one submission
with this backend.

Run:
    python vqe/emulator_offdiag_validate.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
from azure_backend import connect_workspace
from azure.quantum.qiskit import AzureQuantumProvider
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Pauli
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2

LABEL = "IXXI"
PAIR = (0, 1)
TARGET = "quantinuum.sim.h2-1e"
SHOTS = 1000
K_PHASES = [0, 1, 2, 3]


def build_measurement_circuit(vec, label):
    """State-prep + basis rotation for `label` + Z-basis measurement.
    Same left-to-right = qubit(n-1)...qubit0 convention as the rest of
    this repo (decompose_pauli_terms, measure_diagonal_noisy)."""
    n = len(label)
    qc = QuantumCircuit(n, n)
    qc.append(StatePreparation(np.asarray(vec)), range(n))
    for i, ch in enumerate(reversed(label)):
        if ch == "X":
            qc.h(i)
        elif ch == "Y":
            qc.sdg(i)
            qc.h(i)
    qc.measure(range(n), range(n))
    return qc


def pauli_expectation_from_counts(counts, label, shots):
    """After the basis-rotation in build_measurement_circuit, every
    non-identity factor of `label` behaves like a Z in the measured
    frame -- same parity formula as the diagonal-only validation."""
    positions = [i for i, c in enumerate(reversed(label)) if c != "I"]
    total = 0.0
    for bitstring, n in counts.items():
        bits = bitstring.replace(" ", "")
        parity = sum(int(bits[len(bits) - 1 - i]) for i in positions) % 2
        total += ((-1) ** parity) * n
    return total / shots


def main():
    print("\n" + "=" * 70)
    print("  Validating the local Quantinuum-like noise model on an")
    print(f"  OFF-DIAGONAL term ({LABEL}, Schmidt pair {PAIR}) against a")
    print("  REAL Quantinuum emulator -- EF cross-term superposition trick")
    print("=" * 70)

    qop_bare, qop_penalized, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_penalized)
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)
    alpha_labels = set(a for a, _, _ in terms)
    if LABEL not in alpha_labels:
        print(f"\n  {LABEL} is not one of this Hamiltonian's alpha-register "
              f"Pauli terms -- pick a real one.")
        sys.exit(1)

    n, m = PAIR
    u_n, u_m = u_vecs[n], u_vecs[m]
    print(f"\n  Schmidt coefficients: lambda_{n}={lambdas[n]:.6f}, lambda_{m}={lambdas[m]:.6f}")

    # 1. Exact (numpy)
    Pmat = Pauli(LABEL).to_matrix()
    exact_cross = complex(u_n.conj() @ Pmat @ u_m)
    print(f"\n  1. Exact <u_{n}|{LABEL}|u_{m}>  = {exact_cross.real:+.6f} {exact_cross.imag:+.6f}i")

    phis = [(u_n + (1j ** k) * u_m) / np.sqrt(2) for k in K_PHASES]

    # 2. Local Aer, Quantinuum-like noise -- all 4 phases in one batched call
    nm = ef.build_noise_model()
    aer_est = AerEstimatorV2(
        options={"backend_options": {"noise_model": nm, "method": "density_matrix"}}
    )
    circs_t = [ef.transpiled_state_prep_circuit(phi) for phi in phis]
    pubs = [(c, [Pauli(LABEL)]) for c in circs_t]
    aer_result = aer_est.run(pubs).result()
    E_local = [float(aer_result[k].data.evs[0]) for k in K_PHASES]
    re_local = (E_local[0] - E_local[2]) / 2
    im_local = (E_local[3] - E_local[1]) / 2
    local_cross = complex(re_local, im_local)
    print(f"  2. Local Aer (Quantinuum-like)   = {local_cross.real:+.6f} {local_cross.imag:+.6f}i  "
          f"(gap from exact: {abs(local_cross - exact_cross):.6f})")

    # 3. Real Quantinuum emulator -- 4 separate jobs, one per phase (backend
    #    constraint: single circuit per job, cannot batch)
    print(f"\n  Connecting to Azure Quantum, submitting 4 jobs (one per phase) to {TARGET} ...")
    ws = connect_workspace()
    backend = AzureQuantumProvider(ws).get_backend(TARGET)

    E_real = []
    job_ids = []
    for k, phi in zip(K_PHASES, phis):
        circ = build_measurement_circuit(phi, LABEL)
        circ_t = transpile(circ, backend=backend)
        job = backend.run(circ_t, shots=SHOTS)
        print(f"    phase k={k}: job id {job.id()} -- waiting ...")
        counts = job.result().get_counts()
        E_real.append(pauli_expectation_from_counts(counts, LABEL, SHOTS))
        job_ids.append(job.id())

    re_real = (E_real[0] - E_real[2]) / 2
    im_real = (E_real[3] - E_real[1]) / 2
    real_cross = complex(re_real, im_real)
    print(f"\n  3. Real {TARGET} emulator        = {real_cross.real:+.6f} {real_cross.imag:+.6f}i  "
          f"(gap from exact: {abs(real_cross - exact_cross):.6f})")

    local_err = abs(local_cross - exact_cross)
    real_err = abs(real_cross - exact_cross)
    local_vs_real_gap = abs(local_cross - real_cross)

    print("\n" + "-" * 70)
    print(f"  local err (|local - exact|) = {local_err:.6f}")
    print(f"  real err  (|real - exact|)  = {real_err:.6f}")
    print(f"  local-vs-real gap           = {local_vs_real_gap:.6f}")
    print("-" * 70)

    results = {
        "label": LABEL,
        "pair": list(PAIR),
        "shots_per_job": SHOTS,
        "job_ids": job_ids,
        "exact": {"real": round(exact_cross.real, 6), "imag": round(exact_cross.imag, 6)},
        "local_aer": {"real": round(local_cross.real, 6), "imag": round(local_cross.imag, 6)},
        "real_emulator": {"real": round(real_cross.real, 6), "imag": round(real_cross.imag, 6)},
        "local_err": round(local_err, 6),
        "real_err": round(real_err, 6),
        "local_vs_real_gap": round(local_vs_real_gap, 6),
    }
    out = os.path.join(os.path.dirname(__file__), "emulator_offdiag_validate_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
