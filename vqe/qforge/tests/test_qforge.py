#!/usr/bin/env python3
"""
qforge/tests/test_qforge.py — exercises the extracted library end-to-end
against known values already established in RESEARCH_LEDGER.md, so the
extraction itself is verified, not just re-stated. Run:
    python vqe/qforge/tests/test_qforge.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from qforge import (
    setup_fragment, fit_all_targets, verify_constant_gate_count, combine_matrices,
    energy_from_alpha_matrices, RawStrategy, CDRStrategy, floor_test, LOW_SIGNAL_CUTOFF,
    build_ansatz, transpile_fixed,
)
from qforge.floor_test import _self_test as floor_test_self_test


def test_floor_test_catches_iteration2():
    assert floor_test_self_test()
    print("  [PASS] floor_test disqualifies iteration 2's known-bad PERTURB_RADIUS sweep")


def test_setup_and_gate_count():
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=6)
    assert abs(p["max_schmidt_tail"]) < 1e-9, "K=6 must be exact for H4 at d=1.0"
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    assert n_ok == 36, f"expected all 36 K=6 targets to converge, got {n_ok}"
    assert worst < 1e-10
    counts = verify_constant_gate_count(solutions)
    assert counts == {11}, f"expected constant 11 CX gates, got {counts}"
    print(f"  [PASS] setup_fragment + fit_all_targets: 36/36 converged, gate count={counts}, "
          f"exact_energy={p['exact_energy']:.6f} Ha")
    return p, solutions


def test_identity_never_rescaled(p, solutions):
    """The identity Pauli's diagonal must be exactly 1.0 in combine_matrices()
    regardless of any per_label_scale passed in -- even a deliberately
    absurd scale must not touch it."""
    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    fake_raw = {name: {l: 0.5 for l in non_id_labels} for name in solutions}
    absurd_scale = {l: 1e-9 for l in non_id_labels}  # would blow up any rescaled value
    mats = combine_matrices(fake_raw, p["alpha_labels"], p["identity_label"], p["K"],
                             per_label_scale=absurd_scale)
    identity_diag = np.diag(mats[p["identity_label"]])
    assert np.allclose(identity_diag, 1.0), f"identity Pauli was rescaled: {identity_diag}"
    print("  [PASS] identity Pauli diagonal stays exactly 1.0 even under an absurd per_label_scale")


def test_low_signal_cutoff_enforced():
    from qforge.mitigation import filtered_pairs
    training = [
        {"slot": "a", "label": "XX", "exact": 0.01, "noisy": 0.005},   # below cutoff, must be dropped
        {"slot": "b", "label": "XX", "exact": 0.5, "noisy": 0.4},
        {"slot": "c", "label": "XX", "exact": -0.6, "noisy": -0.5},
        {"slot": "d", "label": "XX", "exact": 0.7, "noisy": 0.6},
    ]
    pairs = filtered_pairs(training, label="XX")
    assert len(pairs) == 3, f"expected the |exact|<{LOW_SIGNAL_CUTOFF} row dropped, got {len(pairs)} pairs"
    assert all(abs(e) >= LOW_SIGNAL_CUTOFF for e, _ in pairs)
    print(f"  [PASS] |exact|<{LOW_SIGNAL_CUTOFF} training pairs are dropped (3/4 kept, as expected)")


def test_opt_level_hardcoded():
    """transpile_fixed has no optimization_level parameter at all."""
    import inspect
    sig = inspect.signature(transpile_fixed)
    assert "optimization_level" not in sig.parameters, \
        "transpile_fixed must not expose optimization_level as a caller-settable parameter"
    print("  [PASS] transpile_fixed has no optimization_level parameter (hardcoded to 0)")


def test_raw_vs_cdr_local_direction():
    """Sanity check the extracted CDRStrategy at least has the right
    SHAPE of interface (not a full noisy re-run, which belongs in the
    heavier iteration scripts, not this fast test) -- construct it from
    a tiny synthetic training set and confirm correct() runs and returns
    K x K matrices with the identity invariant intact."""
    from qforge.mitigation import CDRStrategy
    K = 2
    alpha_labels = ["II", "XX", "ZZ"]
    identity_label = "II"
    non_id_labels = ["XX", "ZZ"]
    training = []
    for lbl, scale in (("XX", 0.9), ("ZZ", 0.8)):
        for e in (0.2, 0.4, -0.3, 0.6, -0.5):
            training.append({"slot": "s", "label": lbl, "exact": e, "noisy": e * scale})
    strat = CDRStrategy(training, non_id_labels, K)
    raw = {
        "u_0": {"XX": 0.9 * 0.5, "ZZ": 0.8 * 0.3}, "u_1": {"XX": 0.9 * 0.4, "ZZ": 0.8 * 0.2},
        "(u0+u1)": {"XX": 0.9 * 0.1, "ZZ": 0.8 * 0.1}, "(u0-u1)": {"XX": 0.9 * 0.05, "ZZ": 0.8 * 0.05},
    }
    mats = strat.correct(raw, alpha_labels, identity_label, K)
    assert np.allclose(np.diag(mats["II"]), 1.0)
    assert abs(mats["XX"][0, 0] - 0.5) < 1e-6, "CDR should recover ~exact after dividing out the known scale"
    print("  [PASS] CDRStrategy recovers the known injected scale on a synthetic example")


if __name__ == "__main__":
    print("qforge test suite")
    test_floor_test_catches_iteration2()
    p, solutions = test_setup_and_gate_count()
    test_identity_never_rescaled(p, solutions)
    test_low_signal_cutoff_enforced()
    test_opt_level_hardcoded()
    test_raw_vs_cdr_local_direction()
    print("\nALL TESTS PASSED")
