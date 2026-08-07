#!/usr/bin/env python3
"""
qforge.ansatz — the fixed-structure, particle-number-conserving 5-angle
ansatz this whole project's forged-energy pipeline is built on, extracted
from fixed_ansatz.py (all physics/derivation unchanged, see that file's
own docstring for the full derivation history and the opt_level=0 bug
that motivated GATE_TRANSPILE_OPT_LEVEL being hardcoded, not a caller
option, below).

INVARIANT ENFORCED IN CODE, not left as a comment: transpile_fixed(...)
does not accept an optimization_level argument at all -- it is hardcoded
to 0. optimization_level>=1 was found (see fixed_ansatz.py) to silently
collapse the 2-qubit gate count for specific fitted angles landing near
periodic special values, which breaks CDR's core assumption (training
and target circuits must be structurally identical). Removing the
parameter, rather than defaulting it, makes that mistake unrepresentable
through this API.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from fixed_ansatz import (  # noqa: F401  -- re-exported, verified physics unchanged
    build_ansatz,
    fit_angles,
    statevector_of,
    max_abs_error,
    N_PARAMS,
    P2_PER_GATE,
    P1_PER_GATE,
)

from qiskit import transpile
from qiskit.transpiler import CouplingMap

GATE_TRANSPILE_OPT_LEVEL = 0  # invariant -- see module docstring


def transpile_fixed(qc, basis_gates, n_qubits=4, coupling_map="full"):
    """The ONLY transpile entry point this package exposes for the fixed
    ansatz -- optimization_level is not a parameter, it is always
    GATE_TRANSPILE_OPT_LEVEL (0)."""
    cmap = CouplingMap.from_full(n_qubits) if coupling_map == "full" else coupling_map
    return transpile(qc, basis_gates=basis_gates, coupling_map=cmap,
                      optimization_level=GATE_TRANSPILE_OPT_LEVEL)


def two_qubit_gate_count(angles, basis_gates=("u3", "cx")):
    qc = build_ansatz(angles)
    t = transpile_fixed(qc, list(basis_gates))
    return t.count_ops().get("cx", 0)


def verify_constant_gate_count(solutions, basis_gates=("u3", "cx")):
    """Re-verify (not assume) that every solved target's angles produce
    the SAME 2-qubit gate count under GATE_TRANSPILE_OPT_LEVEL -- the
    condition CDR structurally depends on. Returns the observed set of
    counts; callers should assert len(...) == 1."""
    counts = set()
    for name in solutions:
        angles = solutions[name]["angles"] if isinstance(solutions[name], dict) else solutions[name]
        counts.add(two_qubit_gate_count(angles, basis_gates))
    return counts


def fit_all_targets(targets, tol=1e-10):
    """Solve fit_angles() for every named target vector, verified <1e-10
    (fixed_ansatz.fit_angles' own default tolerance) -- the same
    all-36-targets verification iteration ran at K=6 (rank6_symmetry_vd.py
    setup()), extracted here as a standalone, reusable function rather
    than inline in a script's main()."""
    solutions = {}
    worst = 0.0
    for name, vec in targets.items():
        angles, err = fit_angles(vec, tol=tol)
        solutions[name] = {"angles": angles.tolist() if hasattr(angles, "tolist") else list(angles),
                            "max_abs_error": err, "verified": bool(err < tol)}
        worst = max(worst, err)
    n_ok = sum(s["verified"] for s in solutions.values())
    return solutions, n_ok, worst
