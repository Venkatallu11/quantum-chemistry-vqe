#!/usr/bin/env python3
"""
gate_fusion.py -- iteration 37 (circuit-level branch). Fuses consecutive
single-qubit native gates (gpi/gpi2) into the minimal native sequence,
using Qiskit's OWN unitary-synthesis machinery (via a `UnitaryGate` +
`transpile(target=...)`) rather than hand-derived Euler-angle math --
reuses well-tested infrastructure, applied SURGICALLY only where runs of
consecutive 1-qubit gates actually exist, not as a blanket re-transpile
of the whole circuit (already checked: `transpile(..., optimization_level
=0..3)` on the FULL composed circuit does NOT fuse this project's
ansatz/basis-rotation seam -- levels 0/1 leave it untouched, levels 2/3
make it WORSE (76 vs 48 gpi2 gates for one representative circuit),
directly confirming the IonQTranspileLevelWarning seen in every run log
this session: IonQ's own recommendation against optimization_level 2+ is
correct for this gate set).

WHY THIS EXISTS: found by direct inspection that `base.compose(basis_qc)`
concatenates two INDEPENDENTLY-compiled native circuits without ever
re-optimizing across the seam -- every qubit shows a clean `gpi2-gpi-gpi2`
Euler-style block from the ansatz immediately followed by another
`gpi2-gpi-gpi2` block from the basis rotation, when mathematically the
two rotations should be multiplied into ONE net unitary and decomposed
into at most 3 native gates, not 6. This is a real, free, zero-physics-
risk gate-count reduction -- purely a compilation fix, verified via exact
unitary comparison (not just "looks shorter"), not a new physics
assumption or a new calibration requirement.

Run (self-test):
    python vqe/gate_fusion.py
"""
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator


def fuse_single_qubit_runs(qc, native_target):
    """Returns a new circuit with every maximal run of consecutive
    single-qubit gates on the SAME qubit (with no intervening multi-qubit
    gate) replaced by its combined unitary, re-synthesized into the
    minimal native sequence via Qiskit's own UnitarySynthesis (through
    `transpile(target=native_target)`). Gates on classical bits
    (measure/barrier) end any pending run. Multi-qubit gates end pending
    runs on ALL qubits they touch."""
    n = qc.num_qubits
    fused = qc.copy_empty_like()
    pending = [None] * n  # per-qubit accumulated Operator, or None if no run active

    def flush(q):
        if pending[q] is None:
            return
        U = pending[q]
        pending[q] = None
        # skip re-synthesis entirely if it's (numerically) identity -- nothing to emit
        if np.allclose(U.data, np.eye(2), atol=1e-12):
            return
        sub = QuantumCircuit(1)
        sub.append(UnitaryGate(U.data), [0])
        sub_native = transpile(sub, target=native_target, optimization_level=1)
        for instr in sub_native.data:
            if instr.operation.name in ("id",):
                continue
            fused.append(instr.operation, [fused.qubits[q]], [])

    for instr in qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qidxs = [qc.find_bit(q).index for q in qargs]
        if len(qidxs) == 1 and not cargs and op.name not in ("barrier", "measure"):
            q = qidxs[0]
            U_gate = Operator(op).data
            if pending[q] is None:
                pending[q] = Operator(U_gate)
            else:
                pending[q] = Operator(U_gate) @ pending[q]  # apply order: existing THEN new
        else:
            for q in qidxs:
                flush(q)
            fused.append(op, qargs, cargs)
    for q in range(n):
        flush(q)
    return fused


def _self_test():
    import sys
    import os
    from collections import Counter
    sys.path.insert(0, os.path.dirname(__file__))
    from qforge import setup_fragment, fit_all_targets
    from task28d_all_gate_zne import optimized_native_circuit
    from task2_fold_response_dataset import native_basis_change
    import ef_fragment as effrag_mod
    from native_stateprep import native_target

    K = 6
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    GATE_NAME = "zz"

    print("=" * 96)
    print("  gate_fusion.py self-test: fuse the ansatz/basis-rotation seam, verify exact equivalence")
    print("=" * 96)

    base = optimized_native_circuit(fixed_solutions["u_0"]["angles"], GATE_NAME)
    combined = effrag_mod.combined_basis_label(["XYYX", "IYYI"])
    basis_qc = native_basis_change(combined, GATE_NAME)
    full = base.compose(basis_qc)
    tgt = native_target(full.num_qubits, GATE_NAME)

    counts_before = Counter(instr.operation.name for instr in full.data)
    fused = fuse_single_qubit_runs(full, tgt)
    counts_after = Counter(instr.operation.name for instr in fused.data)

    U_before = Operator(full).data
    U_after = Operator(fused).data
    mask = np.abs(U_before) > 1e-9
    idx = np.argwhere(mask)[0]
    phase = U_after[tuple(idx)] / U_before[tuple(idx)]
    diff = np.max(np.abs(U_after - phase * U_before))

    print(f"  BEFORE: {dict(counts_before)}  depth={full.depth()}")
    print(f"  AFTER:  {dict(counts_after)}  depth={fused.depth()}")
    print(f"  unitary diff (up to global phase): {diff:.3e}  (must be ~1e-9 or smaller)")
    print(f"  EXACT EQUIVALENCE: {'PASS' if diff < 1e-8 else 'FAIL -- DO NOT TRUST THIS FUSION'}")

    gpi2_before, gpi2_after = counts_before.get("gpi2", 0), counts_after.get("gpi2", 0)
    gpi_before, gpi_after = counts_before.get("gpi", 0), counts_after.get("gpi", 0)
    pct_gpi2 = 100 * (gpi2_before - gpi2_after) / gpi2_before if gpi2_before else 0
    pct_gpi = 100 * (gpi_before - gpi_after) / gpi_before if gpi_before else 0
    print(f"\n  GPi2: {gpi2_before} -> {gpi2_after}  ({pct_gpi2:+.1f}%)")
    print(f"  GPi:  {gpi_before} -> {gpi_after}  ({pct_gpi:+.1f}%)")


if __name__ == "__main__":
    _self_test()
