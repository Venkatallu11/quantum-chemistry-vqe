#!/usr/bin/env python3
"""ONE circuit on the real Quantinuum emulator — cost probe before the full run."""
import os
import sys
from dotenv import load_dotenv
load_dotenv()

sys.path.insert(0, os.path.dirname(__file__))
from azure_backend import connect_workspace
from azure.quantum.qiskit import AzureQuantumProvider
from qiskit import QuantumCircuit

TARGET = "quantinuum.sim.h2-1e"
SHOTS = 100


def build_probe_circuit():
    qc = QuantumCircuit(4, 4)
    qc.h(0); qc.cx(0, 1); qc.cx(1, 2); qc.cx(2, 3)
    qc.measure(range(4), range(4))
    return qc


def main():
    ws = connect_workspace()
    backend = AzureQuantumProvider(ws).get_backend(TARGET)
    qc = build_probe_circuit()
    try:
        c = backend.estimate_cost(qc, shots=SHOTS)
        print(f"\n  Cost for THIS 1 job ({SHOTS} shots): {c.estimated_total} {c.currency_code}")
        print(f"  => full ~270-circuit run rough cost: {c.estimated_total*270:.2f} {c.currency_code}")
    except Exception as e:
        print(f"\n  (cost preview not available: {e})")
    print(f"\n  Submitting 1 circuit to {TARGET} ...")
    job = backend.run(qc, shots=SHOTS)
    print(f"  job id: {job.id()} — waiting ...")
    for b, n in sorted(job.result().get_counts().items(), key=lambda x: -x[1]):
        print(f"    {b}: {n}")


if __name__ == "__main__":
    main()
