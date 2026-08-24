#!/usr/bin/env python3
"""
task39b_native_ancilla_parity.py -- iteration 39, Task B. Adapts
iteration 18's real-hardware-validated ancilla-parity leakage detector
(`spin_leakage_postselect_ionq.py`) from its original abstract-gate/CX
construction to the CURRENT native-gate (GPi/GPi2/ZZ) circuit family
Tasks 36-38 actually use -- the first step toward combining it with PEC
+ the joint Schmidt frame, which iteration 18 never had access to.

DESIGN, matching the project's own established fusion-builder pattern
(`gate_fusion.build_fused_measurement_circuit`, Task 38A): build the
register's native state-prep circuit UNCHANGED (whatever this project's
current pipeline already produces), and the ancilla-CNOT sub-circuit
SEPARATELY, transpiled to native gates IN ISOLATION at optimization_
level=1 (matching `to_native`'s own established convention, NOT the
higher levels this project has repeatedly found unsafe for IonQ's own
compiler) -- then compose the two, never letting the transpiler see
both pieces at once (the same reasoning iteration 38A's fusion builder
already established: composing separately-verified native pieces is
safer than re-transpiling a combined circuit and trusting the compiler
not to do something unexpected to it).

VERIFICATION, matching iteration 18's own `verify_ancilla_scheme`
exactly (same two checks, same tolerance), now for the native version:
(1) tracing out the ancilla must leave the original 4-qubit register's
marginal state EXACTLY unchanged; (2) for the noiseless ideal state, the
ancilla must read 0 with probability exactly 1 (certain, no leakage).
Both checked to machine precision before this circuit is trusted for
anything downstream.

Run:
    PYTHONHASHSEED=0 python vqe/task39b_native_ancilla_parity.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from native_stateprep import to_native, native_target
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator, partial_trace

GATE_NAME = "zz"


def ancilla_cnots_abstract():
    qc = QuantumCircuit(5)
    qc.cx(0, 4)
    qc.cx(1, 4)
    qc.cx(2, 4)
    qc.cx(3, 4)
    return qc


def ancilla_cnots_native(gate_name=GATE_NAME):
    return to_native(ancilla_cnots_abstract(), gate_name)


def with_ancilla_parity_native(register_native_4q, ancilla_native_5q):
    """register_native_4q: this project's own existing native state-prep
    circuit (4 qubits, unchanged, built exactly as everywhere else in
    Tasks 27-38). ancilla_native_5q: the verified, natively-transpiled
    CNOT-parity sub-circuit (built once, reused for every slot -- it
    does not depend on the state-prep angles at all)."""
    qc5 = QuantumCircuit(5)
    qc5.compose(register_native_4q, qubits=[0, 1, 2, 3], inplace=True)
    qc5.compose(ancilla_native_5q, inplace=True)
    return qc5


def verify_native_ancilla(ancilla_native_5q):
    """Regression check: does the NATIVELY-TRANSPILED ancilla circuit
    reproduce the ABSTRACT CX version's unitary exactly?"""
    U_native = Operator(ancilla_native_5q).data
    U_abstract = Operator(ancilla_cnots_abstract()).data
    mask = np.abs(U_abstract) > 1e-9
    idx = np.argwhere(mask)[0]
    phase = U_native[tuple(idx)] / U_abstract[tuple(idx)]
    diff = float(np.max(np.abs(U_native - phase * U_abstract)))
    return diff


def verify_ancilla_scheme_native(angles, gate_name=GATE_NAME):
    """Matches iteration 18's verify_ancilla_scheme exactly, for the
    native version: (1) marginal-state preservation, (2) ideal ancilla=0
    certainty."""
    from fixed_ansatz import build_ansatz
    base_native = to_native(build_ansatz(angles), gate_name)
    ancilla_native = ancilla_cnots_native(gate_name)
    qc5 = with_ancilla_parity_native(base_native, ancilla_native)

    sv5 = Statevector.from_instruction(qc5)
    sv4 = Statevector.from_instruction(base_native)
    rho4_direct = np.outer(np.asarray(sv4), np.asarray(sv4).conj())
    rho4_marginal = np.asarray(partial_trace(sv5, [4]))
    marginal_err = float(np.max(np.abs(rho4_direct - rho4_marginal)))

    probs = sv5.probabilities_dict()
    p_ancilla1 = sum(v for k, v in probs.items() if k[0] == "1")
    return marginal_err, p_ancilla1


def _self_test():
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(__file__))
    from qforge import setup_fragment, fit_all_targets

    print("\n" + "=" * 96)
    print("  task39b_native_ancilla_parity.py -- native ancilla-parity circuit, verified before use")
    print("=" * 96)

    ancilla_native = ancilla_cnots_native()
    print(f"\n  ancilla sub-circuit native gate counts: {ancilla_native.count_ops()}")
    diff = verify_native_ancilla(ancilla_native)
    print(f"  regression vs abstract CX circuit: unitary diff = {diff:.3e}  {'PASS' if diff < 1e-8 else 'FAIL -- STOP'}")
    if diff >= 1e-8:
        raise RuntimeError("native ancilla circuit does not match the abstract CX circuit -- do not trust it")

    K = 6
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)

    print(f"\n  -- verify_ancilla_scheme_native on {len(fixed_solutions)} state-prep angle sets --")
    worst_marginal, worst_p1 = 0.0, 0.0
    for name, sol in fixed_solutions.items():
        marginal_err, p_ancilla1 = verify_ancilla_scheme_native(sol["angles"])
        worst_marginal = max(worst_marginal, marginal_err)
        worst_p1 = max(worst_p1, p_ancilla1)
    print(f"    worst marginal-state error across all slots: {worst_marginal:.3e}  "
          f"{'PASS' if worst_marginal < 1e-9 else 'FAIL'}")
    print(f"    worst ideal ancilla=1 probability across all slots: {worst_p1:.3e}  "
          f"{'PASS (certain ancilla=0)' if worst_p1 < 1e-9 else 'FAIL'}")
    print(f"\n  OVERALL: {'VERIFIED, safe to use downstream' if diff < 1e-8 and worst_marginal < 1e-9 and worst_p1 < 1e-9 else 'DO NOT TRUST -- investigate before proceeding'}")


if __name__ == "__main__":
    _self_test()
