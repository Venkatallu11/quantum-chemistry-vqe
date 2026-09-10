import pytest

qiskit = pytest.importorskip("qiskit")

from general_commuting_measurements import (
    Diagonalizer,
    build_general_commuting_measurement_plan,
    expectations_from_counts,
    minimum_general_commuting_groups,
    pauli_commutes,
)


def test_general_commutation_rules():
    assert pauli_commutes("XX", "YY")
    assert pauli_commutes("XX", "ZZ")
    assert pauli_commutes("YY", "ZZ")
    assert not pauli_commutes("XX", "ZI")


def test_exact_minimum_for_bell_stabilizer_group():
    groups = minimum_general_commuting_groups(["XX", "YY", "ZZ"])
    assert len(groups) == 1
    assert set(groups[0]) == {"XX", "YY", "ZZ"}


def test_diagonalizer_maps_group_to_z_type():
    groups, diagonals = build_general_commuting_measurement_plan(["XX", "YY", "ZZ"], beam_width=128, max_depth=12)
    assert len(groups) == 1
    d = diagonals[0]
    assert all(set(x.diagonal) <= {"I", "Z"} for x in d.transformed)


def test_measurement_reconstruction_on_bell_state():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector

    groups, diagonals = build_general_commuting_measurement_plan(["XX", "YY", "ZZ"], beam_width=128, max_depth=12)
    d = diagonals[0]

    # Bell state |Phi+> has XX=+1, YY=-1, ZZ=+1.
    base = QuantumCircuit(2)
    base.h(0)
    base.cx(0, 1)
    meas = base.compose(d.to_circuit())
    sv = Statevector.from_instruction(meas)
    probs = sv.probabilities_dict()
    counts = {k: int(round(v * 1_000_000)) for k, v in probs.items()}
    vals = expectations_from_counts(counts, d)

    assert vals["XX"] == pytest.approx(1.0, abs=1e-6)
    assert vals["YY"] == pytest.approx(-1.0, abs=1e-6)
    assert vals["ZZ"] == pytest.approx(1.0, abs=1e-6)
