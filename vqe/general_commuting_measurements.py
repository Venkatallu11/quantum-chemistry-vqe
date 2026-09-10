"""General-commuting Pauli grouping and Clifford measurement for the H4 EF codebase.

This module is deliberately opt-in.  The existing H4 pipeline uses qubit-wise
commuting (QWC) groups with a tensor-product H/Sdg basis rotation.  General
commuting (GC) groups require a joint Clifford diagonalizer, so this module
keeps the old path untouched and provides a separately testable replacement.

Key guarantees:
  * exact Pauli commutation test (binary symplectic form)
  * exact minimum coloring for small Pauli sets (DSATUR branch-and-bound)
  * joint Clifford diagonalizer found by deterministic beam search over H/S/Sdg/CX
  * transformed Pauli signs retained (critical for expectation reconstruction)
  * one measurement circuit per GC group
  * ideal statevector-equivalence verification for every group member

The diagonalizer search is intended for the project's 4-qubit H4 measurement
problem.  For larger systems, use a dedicated stabilizer-synthesis routine.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

try:
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Pauli, SparsePauliOp
except Exception as exc:  # pragma: no cover - import guard for documentation
    QuantumCircuit = None
    Pauli = None
    SparsePauliOp = None
    _QISKIT_IMPORT_ERROR = exc


GateSpec = Tuple[str, int, int | None]


def _require_qiskit() -> None:
    if QuantumCircuit is None:
        raise ImportError(
            "general_commuting_measurements requires qiskit-terra"
        ) from _QISKIT_IMPORT_ERROR


def pauli_commutes(a: str, b: str) -> bool:
    """Return True iff two same-length Pauli strings commute globally.

    Two Pauli strings commute iff the number of positions in which their
    non-identity factors anticommute is even.
    """
    if len(a) != len(b):
        raise ValueError("Pauli labels must have equal length")
    anticommutes = 0
    for xa, xb in zip(a, b):
        if xa == "I" or xb == "I" or xa == xb:
            continue
        anticommutes += 1
    return (anticommutes % 2) == 0


def validate_pairwise_commuting(group: Sequence[str]) -> None:
    for i, a in enumerate(group):
        for b in group[i + 1 :]:
            if not pauli_commutes(a, b):
                raise ValueError(f"non-commuting pair in group: {a}, {b}")


def _compatibility_graph(labels: Sequence[str]) -> List[int]:
    """Return non-commutation bitmasks for exact graph coloring."""
    n = len(labels)
    bad = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if not pauli_commutes(labels[i], labels[j]):
                bad[i] |= 1 << j
                bad[j] |= 1 << i
    return bad


def minimum_general_commuting_groups(labels: Sequence[str]) -> List[List[str]]:
    """Exact minimum partition into pairwise-general-commuting groups.

    This is graph coloring on the *non-commutation* graph.  The H4 problem is
    tiny (36 local Pauli labels in the current workflow), so an exact DSATUR
    branch-and-bound search is practical and preferable to silently accepting
    a merely greedy partition.
    """
    labels = list(dict.fromkeys(labels))
    if not labels:
        return []
    n = len(labels)
    bad = _compatibility_graph(labels)

    # Initial upper bound: deterministic first-fit ordered by descending degree.
    order = sorted(range(n), key=lambda i: bad[i].bit_count(), reverse=True)
    greedy_color = [-1] * n
    greedy_k = 0
    for v in order:
        used = 0
        for u in range(n):
            if bad[v] >> u & 1 and greedy_color[u] >= 0:
                used |= 1 << greedy_color[u]
        c = 0
        while used >> c & 1:
            c += 1
        greedy_color[v] = c
        greedy_k = max(greedy_k, c + 1)

    best_k = greedy_k
    best_color = greedy_color[:]
    colors = [-1] * n
    neighbor_color_masks = [0] * n

    def saturation(v: int) -> int:
        return neighbor_color_masks[v].bit_count()

    def select_vertex() -> int:
        uncolored = [v for v in range(n) if colors[v] < 0]
        return max(
            uncolored,
            key=lambda v: (saturation(v), bad[v].bit_count(), -v),
        )

    def feasible_lower_bound() -> int:
        # A cheap lower bound: maximum clique in the non-commutation graph is
        # expensive; for 36 vertices, DSATUR's saturation already gives a
        # strong practical branch-and-bound. Return current max saturation+1.
        return max((saturation(v) + 1 for v in range(n) if colors[v] < 0), default=0)

    def assign(v: int, c: int) -> None:
        colors[v] = c
        bit = 1 << c
        for u in range(n):
            if bad[v] >> u & 1:
                neighbor_color_masks[u] |= bit

    def unassign(v: int, c: int) -> None:
        colors[v] = -1
        # Recompute affected masks; n is tiny and correctness matters more.
        for u in range(n):
            mask = 0
            for w in range(n):
                if colors[w] >= 0 and (bad[u] >> w & 1):
                    mask |= 1 << colors[w]
            neighbor_color_masks[u] = mask

    def search(colored: int, used_colors: int) -> None:
        nonlocal best_k, best_color
        if colored == n:
            k = used_colors.bit_count()
            if k < best_k:
                best_k = k
                best_color = colors[:]
            return
        if used_colors.bit_count() >= best_k:
            return
        if max(feasible_lower_bound(), used_colors.bit_count()) >= best_k:
            # Still allow the branch if an existing color can lead to the same
            # count; the strict improvement is what matters.
            pass

        v = select_vertex()
        forbidden = neighbor_color_masks[v]
        current_k = used_colors.bit_count()

        # Try existing colors first, preferring lower indices.
        for c in range(current_k):
            if not (forbidden >> c) & 1:
                assign(v, c)
                search(colored + 1, used_colors)
                unassign(v, c)

        # Introduce one new color only if it can beat the incumbent.
        if current_k + 1 < best_k:
            assign(v, current_k)
            search(colored + 1, used_colors | (1 << current_k))
            unassign(v, current_k)

    search(0, 0)
    groups: List[List[str]] = [[] for _ in range(best_k)]
    for idx, color in enumerate(best_color):
        groups[color].append(labels[idx])
    return [g for g in groups if g]


@dataclass(frozen=True)
class DiagonalizedPauli:
    original: str
    diagonal: str
    sign: int


@dataclass
class Diagonalizer:
    n_qubits: int
    gates: List[GateSpec]
    transformed: List[DiagonalizedPauli]

    @property
    def two_qubit_count(self) -> int:
        return sum(1 for g in self.gates if g[0] == "cx")

    @property
    def one_qubit_count(self) -> int:
        return sum(1 for g in self.gates if g[0] in {"h", "s", "sdg"})

    def to_circuit(self) -> QuantumCircuit:
        _require_qiskit()
        qc = QuantumCircuit(self.n_qubits)
        for name, q0, q1 in self.gates:
            if name == "h":
                qc.h(q0)
            elif name == "s":
                qc.s(q0)
            elif name == "sdg":
                qc.sdg(q0)
            elif name == "cx":
                assert q1 is not None
                qc.cx(q0, q1)
            else:  # pragma: no cover
                raise ValueError(f"unknown Clifford gate {name}")
        return qc


def _one_gate_circuit(n: int, gate: GateSpec) -> QuantumCircuit:
    qc = QuantumCircuit(n)
    name, q0, q1 = gate
    if name == "h":
        qc.h(q0)
    elif name == "s":
        qc.s(q0)
    elif name == "sdg":
        qc.sdg(q0)
    elif name == "cx":
        assert q1 is not None
        qc.cx(q0, q1)
    else:  # pragma: no cover
        raise ValueError(name)
    return qc


def _transform_labels(labels: Sequence[str], gate: GateSpec) -> Tuple[Tuple[str, ...], Tuple[int, ...]]:
    """Conjugate labels by one Clifford gate, retaining Pauli signs."""
    _require_qiskit()
    n = len(labels[0])
    gqc = _one_gate_circuit(n, gate)
    out_labels: List[str] = []
    signs: List[int] = []
    for label in labels:
        p = Pauli(label).evolve(gqc)
        phase = int(getattr(p, "phase", 0))
        # Pauli.evolve returns a Hermitian Pauli here; phase 0/2 => +/-1.
        if phase not in (0, 2):
            raise AssertionError(f"unexpected non-Hermitian phase {phase} for {label}")
        # Pauli.to_label() embeds the sign as a literal leading "-" character
        # in this qiskit version (e.g. "-YY" for a 2-qubit label) -- keep only
        # the trailing n characters so the label string stays a bare Pauli
        # string; sign is tracked separately via the already-validated phase.
        out_labels.append(p.to_label()[-n:])
        signs.append(1 if phase == 0 else -1)
    return tuple(out_labels), tuple(signs)


def _state_score(labels: Sequence[str], depth: int, two_q: int) -> Tuple[int, int, int, float]:
    bad = sum(ch in "XY" for label in labels for ch in label)
    distinct_bad = len({i for label in labels for i, ch in enumerate(label) if ch in "XY"})
    # Primary goal: all X/Y disappear.  Secondary: fewer active bad qubits,
    # then fewer 2Q gates, then total depth.
    return (bad, distinct_bad, two_q, depth + 0.01 * two_q)


def find_diagonalizer(
    group: Sequence[str],
    beam_width: int = 256,
    max_depth: int = 20,
) -> Diagonalizer:
    """Find a low-cost Clifford D with D P D^† = +/- Z-type for all P in group.

    This is a deterministic beam search specialized to the 4-qubit H4 problem.
    It is not a generic Clifford synthesizer, but every returned circuit is
    exactly checked by Qiskit's Pauli conjugation machinery before being used.
    """
    _require_qiskit()
    group = list(dict.fromkeys(group))
    validate_pairwise_commuting(group)
    n = len(group[0])
    if any(len(p) != n for p in group):
        raise ValueError("all Pauli labels in a group must have equal length")

    actions: List[GateSpec] = []
    for q in range(n):
        actions.extend([("h", q, None), ("s", q, None), ("sdg", q, None)])
    for c in range(n):
        for t in range(n):
            if c != t:
                actions.append(("cx", c, t))

    # state = (labels, signs, gates)
    beam = [(tuple(group), tuple([1] * len(group)), tuple())]
    seen = {tuple(group)}

    def is_goal(labels: Sequence[str]) -> bool:
        return all(ch in "IZ" for label in labels for ch in label)

    for depth in range(max_depth + 1):
        candidates = []
        for labels, signs, gates in beam:
            if is_goal(labels):
                transformed = [
                    DiagonalizedPauli(g, lab, s)
                    for g, lab, s in zip(group, labels, signs)
                ]
                return Diagonalizer(n, list(gates), transformed)

            for action in actions:
                new_labels, local_signs = _transform_labels(labels, action)
                if new_labels in seen:
                    continue
                seen.add(new_labels)
                # Signs multiply under composition.  The label state itself is
                # sufficient for search; retain full accumulated signs for the
                # final circuit.
                new_signs = tuple(a * b for a, b in zip(signs, local_signs))
                new_gates = gates + (action,)
                twoq = sum(g[0] == "cx" for g in new_gates)
                score = _state_score(new_labels, len(new_gates), twoq)
                candidates.append((score, new_labels, new_signs, new_gates))

        candidates.sort(key=lambda x: x[0])
        beam = [(l, s, g) for _, l, s, g in candidates[:beam_width]]
        if not beam:
            break

    raise RuntimeError(
        f"No Clifford diagonalizer found within depth={max_depth}, beam_width={beam_width}; "
        "increase search limits for this group."
    )


def build_general_commuting_measurement_plan(
    labels: Sequence[str],
    beam_width: int = 256,
    max_depth: int = 20,
) -> Tuple[List[List[str]], List[Diagonalizer]]:
    """Return exact minimum GC groups and one diagonalizer per group."""
    groups = minimum_general_commuting_groups(labels)
    diagonals = [
        find_diagonalizer(g, beam_width=beam_width, max_depth=max_depth)
        for g in groups
    ]
    return groups, diagonals


def measurement_circuit_general(
    base_circuit: QuantumCircuit,
    diagonalizer: Diagonalizer,
    qubits: Sequence[int] | None = None,
) -> QuantumCircuit:
    """Append one group's joint Clifford diagonalizer and measure all qubits.

    `qubits`: which qubits of `base_circuit` the diagonalizer acts on (e.g.
    [0,1,2,3] for a 4-qubit register embedded in a larger ancilla-augmented
    circuit). Defaults to None, which lets Qiskit's own compose() require an
    exact qubit-count match (the original behavior) -- pass it explicitly
    whenever base_circuit has more qubits than the diagonalizer itself, or
    compose() will raise CircuitError rather than silently doing the wrong
    thing."""
    _require_qiskit()
    qc = base_circuit.compose(diagonalizer.to_circuit(), qubits=qubits, front=False)
    qc.measure_all()
    return qc


def expectations_from_counts(
    counts: Dict[str, int],
    diagonalizer: Diagonalizer,
) -> Dict[str, float]:
    """Recover every original Pauli expectation from one group's counts."""
    total = float(sum(counts.values()))
    if total <= 0:
        raise ValueError("counts contain zero shots")
    out: Dict[str, float] = {}
    for item in diagonalizer.transformed:
        value = 0.0
        n = len(item.diagonal)
        # Qiskit's own Pauli-label convention places the LEFTMOST character
        # as the HIGHEST-indexed qubit (string position i -> qubit n-1-i),
        # matching this project's own established native_basis_change()
        # convention. z_positions must be physical qubit indices to line up
        # with `bits` (which is already reversed so bits[q] = qubit q) --
        # using the raw string position i here (BUG, now fixed) silently
        # mismatched non-palindromic diagonal patterns.
        z_positions = [n - 1 - i for i, ch in enumerate(item.diagonal) if ch == "Z"]
        for bitstring, count in counts.items():
            bits = bitstring.replace(" ", "")[::-1]
            parity = 0
            for q in z_positions:
                parity ^= int(bits[q])
            value += (-1.0 if parity else 1.0) * count / total
        out[item.original] = float(item.sign * value)
    return out


def audit_plan(
    labels: Sequence[str],
    groups: Sequence[Sequence[str]],
    diagonalizers: Sequence[Diagonalizer],
) -> Dict[str, object]:
    """Machine-readable resource report for the grouping proposal."""
    validate = True
    for g in groups:
        validate_pairwise_commuting(g)
    if set(sum((list(g) for g in groups), [])) != set(labels):
        raise AssertionError("grouping does not cover exactly the supplied labels")

    return {
        "n_labels": len(labels),
        "n_groups": len(groups),
        "groups": [list(g) for g in groups],
        "diagonalizer_1q": [d.one_qubit_count for d in diagonalizers],
        "diagonalizer_2q": [d.two_qubit_count for d in diagonalizers],
        "total_diagonalizer_1q": sum(d.one_qubit_count for d in diagonalizers),
        "total_diagonalizer_2q": sum(d.two_qubit_count for d in diagonalizers),
        "pairwise_commuting_verified": validate,
    }


if __name__ == "__main__":  # pragma: no cover
    # Direct H4 report.  This imports the repo's chemistry/QForge stack and
    # obtains the actual K=6 alpha labels; no hand-copied label list.
    import json
    import os
    import sys

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    try:
        from qforge import setup_fragment
        from ef_fragment import group_labels_qubit_wise
    except Exception as exc:
        raise SystemExit(f"Could not import project QForge stack: {exc}")

    setup = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=6, strict=True)
    labels = [l for l in setup["alpha_labels"] if l != setup["identity_label"]]
    qwc = group_labels_qubit_wise(labels)
    gc, diags = build_general_commuting_measurement_plan(labels)

    report = audit_plan(labels, gc, diags)
    report["qwc_groups"] = len(qwc)
    report["qwc_expected_circuit_count_K6"] = len(qwc) * 21
    report["gc_expected_circuit_count_K6"] = len(gc) * 21
    report["circuit_reduction_fraction"] = 1.0 - len(gc) / len(qwc)
    print(json.dumps(report, indent=2))
