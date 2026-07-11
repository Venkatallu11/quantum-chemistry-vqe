#!/usr/bin/env python3
"""
emulator_ef_validate.py — validate the local Quantinuum-like noise model
used in entanglement_forging_h4.py / entanglement_forging_zne.py against
a REAL Quantinuum emulator job, across the 10 largest diagonal terms of
the H4 Hamiltonian's alpha register (not just one).
============================================================================
The single-term version of this script (ZIII only) found the local noise
model underestimates real Quantinuum-grade noise by ~3x on that one term.
One term isn't enough to know if that's a general pattern or a fluke, but
submitting 10 separate real jobs would mean paying for 10x the shots.

Instead: state-prepare the alpha register's leading Schmidt vector
(u_vecs[0], 4 qubits) and measure ONCE in the computational (Z) basis
with 4000 shots. That single shot record contains the full 4-qubit
bitstring distribution, from which the expectation value of ANY diagonal
(I/Z-only) Pauli observable can be computed classically -- so one real
job validates all 10 terms, not just one.

Runs the SAME 10 diagonal observables three ways:
  1. Exact (numpy)                                    -- ground truth
  2. Local Aer, Quantinuum-like noise model (scale=1)  -- what this repo's
     EF/ZNE results have been assuming Quantinuum-grade noise looks like
  3. Real quantinuum.sim.h2-1e emulator, ONE 4000-shot job -- what
     Quantinuum's own emulator actually predicts

Run:
    python vqe/emulator_ef_validate.py
"""
import os
import sys
import json
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import entanglement_forging_h4 as ef
from azure_backend import connect_workspace
from azure.quantum.qiskit import AzureQuantumProvider
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Pauli
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2

N_TERMS = 10
TARGET = "quantinuum.sim.h2-1e"
SHOTS = 4000


def top_diagonal_labels(terms, n):
    """The n alpha-register I/Z-only (diagonal) labels with the largest
    total |coefficient| weight in the Hamiltonian -- ranked by physical
    importance, not chosen arbitrarily."""
    weight = defaultdict(float)
    for a, _, coeff in terms:
        if set(a) <= {"I", "Z"} and a != "IIII":
            weight[a] += abs(coeff.real)
    ranked = sorted(weight.items(), key=lambda x: -x[1])
    return [label for label, _ in ranked[:n]]


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
    print(f"  REAL Quantinuum emulator job -- {N_TERMS} diagonal EF observables,")
    print(f"  ONE {SHOTS}-shot job")
    print("=" * 70)

    qop_bare, qop_penalized, enuc = ef.build_h4_qop(1.0)
    e_elec, psi = ef.exact_ground_state(qop_penalized)
    lambdas, u_vecs, v_vecs = ef.schmidt_decompose(psi)
    terms = ef.decompose_pauli_terms(qop_bare)

    labels = top_diagonal_labels(terms, N_TERMS)
    print(f"\n  Top {len(labels)} diagonal alpha-register terms by |coefficient| weight:")
    print(f"    {labels}")

    u0 = u_vecs[0]  # leading Schmidt vector, alpha register (4 qubits)

    # 1. Exact (numpy), all labels
    exact_vals = {}
    for label in labels:
        Pz = Pauli(label).to_matrix()
        exact_vals[label] = float(np.real(u0.conj() @ Pz @ u0))

    # 2. Local Aer, Quantinuum-like noise, all labels in one batched call
    circ_t = ef.transpiled_state_prep_circuit(u0)
    nm = ef.build_noise_model()
    aer_est = AerEstimatorV2(
        options={"backend_options": {"noise_model": nm, "method": "density_matrix"}}
    )
    pub = (circ_t, [Pauli(l) for l in labels])
    aer_evs = aer_est.run([pub]).result()[0].data.evs
    aer_vals = {label: float(v) for label, v in zip(labels, aer_evs)}

    # 3. Real Quantinuum emulator -- ONE job, all labels derived from it
    print(f"\n  Connecting to Azure Quantum, submitting ONE {SHOTS}-shot job to {TARGET} ...")
    ws = connect_workspace()
    backend = AzureQuantumProvider(ws).get_backend(TARGET)

    meas_circ = QuantumCircuit(4, 4)
    meas_circ.append(StatePreparation(np.asarray(u0)), range(4))
    meas_circ.measure(range(4), range(4))
    meas_circ_t = transpile(meas_circ, backend=backend)

    job = backend.run(meas_circ_t, shots=SHOTS)
    print(f"  job id: {job.id()} -- waiting ...")
    counts = job.result().get_counts()
    real_vals = {label: zlabel_expectation_from_counts(counts, label, SHOTS) for label in labels}

    # --- report ---
    print(f"\n  {'label':<8} {'exact':>10} {'local Aer':>12} {'real emu':>12} "
          f"{'local err':>10} {'real err':>10} {'local-vs-real gap':>18}")
    print("  " + "-" * 82)
    per_term = {}
    for label in labels:
        e, a, r = exact_vals[label], aer_vals[label], real_vals[label]
        local_err = abs(a - e)
        real_err = abs(r - e)
        gap = abs(a - r)
        print(f"  {label:<8} {e:>10.6f} {a:>12.6f} {r:>12.6f} "
              f"{local_err:>10.6f} {real_err:>10.6f} {gap:>18.6f}")
        per_term[label] = {
            "exact": round(e, 6), "local_aer": round(a, 6), "real_emulator": round(r, 6),
            "local_err": round(local_err, 6), "real_err": round(real_err, 6),
            "local_vs_real_gap": round(gap, 6),
        }

    mean_local_err = float(np.mean([per_term[l]["local_err"] for l in labels]))
    mean_real_err = float(np.mean([per_term[l]["real_err"] for l in labels]))
    mean_gap = float(np.mean([per_term[l]["local_vs_real_gap"] for l in labels]))
    ratio = mean_real_err / mean_local_err if mean_local_err > 0 else float("inf")

    print("\n" + "-" * 70)
    print(f"  mean |local err|  = {mean_local_err:.6f}")
    print(f"  mean |real err|   = {mean_real_err:.6f}")
    print(f"  mean local-vs-real gap = {mean_gap:.6f}")
    print(f"  real/local error ratio = {ratio:.2f}x")
    print("-" * 70)

    results = {
        "n_terms": len(labels),
        "labels": labels,
        "shots": SHOTS,
        "job_id": job.id(),
        "per_term": per_term,
        "mean_local_err": round(mean_local_err, 6),
        "mean_real_err": round(mean_real_err, 6),
        "mean_local_vs_real_gap": round(mean_gap, 6),
        "real_over_local_error_ratio": round(ratio, 4),
    }
    out = os.path.join(os.path.dirname(__file__), "emulator_ef_validate_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
