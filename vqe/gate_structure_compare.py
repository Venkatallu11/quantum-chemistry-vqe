#!/usr/bin/env python3
"""
gate_structure_compare.py — diagnostic: compare the abstract 11-gate
ansatz's and the Z2-tapered circuit's actual gate STRUCTURE (not just
2-qubit gate count), to understand the real-hardware gap the user asked
about (abstract 34.98/43.03 kcal/mol vs tapered-StatePreparation
47.78/51.25, despite fewer 2-qubit gates).
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

print("importing...", flush=True)
t0 = time.time()
from fixed_ansatz import build_ansatz
from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import StatePreparation
import rank6_symmetry_vd as r
print(f"imports done, {time.time()-t0:.1f}s", flush=True)

BASIS = ["u3", "cx"]

t0 = time.time()
p = r.setup(6)
print(f"setup done, {time.time()-t0:.1f}s", flush=True)

t0 = time.time()
solutions, n_ok, worst = r.fit_all_targets(p["targets"])
print(f"fit_all_targets done ({n_ok}/36), {time.time()-t0:.1f}s", flush=True)

t0 = time.time()
n2q_list, n1q_list = [], []
for name, sol in solutions.items():
    qc = build_ansatz(sol["angles"])
    t = transpile(qc, basis_gates=BASIS, optimization_level=0)
    ops = t.count_ops()
    n2q_list.append(ops.get("cx", 0))
    n1q_list.append(ops.get("u3", 0))
print(f"abstract ansatz gate counts done, {time.time()-t0:.1f}s", flush=True)
print(f"ABSTRACT 11-gate ansatz: 2q={set(n2q_list)}  1q min={min(n1q_list)} max={max(n1q_list)} mean={np.mean(n1q_list):.2f}")

t0 = time.time()
from z2_tapered_zne import build_reduced_problem
problem = build_reduced_problem()
print(f"build_reduced_problem done, {time.time()-t0:.1f}s", flush=True)

t0 = time.time()
n2q_list_r, n1q_list_r = [], []
for name, vec in problem["reduced_targets"].items():
    qc = QuantumCircuit(3)
    qc.append(StatePreparation(vec), range(3))
    t = transpile(qc, basis_gates=BASIS, optimization_level=0)
    ops = t.count_ops()
    n2q_list_r.append(ops.get("cx", 0))
    n1q_list_r.append(ops.get("u3", 0))
print(f"tapered circuit gate counts done, {time.time()-t0:.1f}s", flush=True)
print(f"TAPERED (StatePreparation) circuit: 2q min={min(n2q_list_r)} max={max(n2q_list_r)} mean={np.mean(n2q_list_r):.2f}  "
      f"1q min={min(n1q_list_r)} max={max(n1q_list_r)} mean={np.mean(n1q_list_r):.2f}")

def print_gate_list(t, title):
    print(f"\n-- {title} ({len(t.data)} instructions) --")
    for instr in t.data:
        op = instr.operation
        qargs = [t.find_bit(q).index for q in instr.qubits]
        params = [round(float(p), 4) for p in op.params] if op.params else []
        print(f"  {op.name:<6} qubits={qargs} params={params}")


print_gate_list(transpile(build_ansatz(solutions["u_0"]["angles"]), basis_gates=BASIS, optimization_level=0),
                 "abstract ansatz circuit (u_0)")

qc2 = QuantumCircuit(3)
qc2.append(StatePreparation(problem["reduced_targets"]["u_0"]), range(3))
print_gate_list(transpile(qc2, basis_gates=BASIS, optimization_level=0), "tapered circuit (u_0)")

# also show a target with the MAX 2q gate count in the tapered set, for the worst case
gc_by_name = {}
for name, vec in problem["reduced_targets"].items():
    qc3 = QuantumCircuit(3)
    qc3.append(StatePreparation(vec), range(3))
    t3 = transpile(qc3, basis_gates=BASIS, optimization_level=0)
    gc_by_name[name] = (t3.count_ops().get("cx", 0), t3)
worst_name = max(gc_by_name, key=lambda k: gc_by_name[k][0])
print_gate_list(gc_by_name[worst_name][1], f"tapered circuit, MAX-2q-gate target ({worst_name}, {gc_by_name[worst_name][0]} CX)")
