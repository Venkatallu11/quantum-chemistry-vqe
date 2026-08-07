#!/usr/bin/env python3
"""
qforge — a clean, importable library for this project's validated
entanglement-forging + noise-mitigation pipeline (H4 fragment, K=6 fixed
ansatz, raw/CDR/PEC strategies, floor testing).

No qforge package existed before this extraction; everything here is
consolidated from vqe/*.py scripts (fixed_ansatz.py, ef_fragment.py,
rank6_symmetry_vd.py, loop_pec.py, shot_noise_study.py,
loop_local_perturbation_cdr.py) with physics and verification logic
unchanged -- see each qforge submodule's docstring for exactly which
source file it consolidates and which invariant it enforces in code
rather than in a comment.

    from qforge import setup_fragment, fit_all_targets, RawStrategy, CDRStrategy, PECStrategy, floor_test
"""
from qforge.ansatz import (  # noqa: F401
    build_ansatz, fit_angles, fit_all_targets, transpile_fixed, two_qubit_gate_count,
    verify_constant_gate_count, GATE_TRANSPILE_OPT_LEVEL, P2_PER_GATE, P1_PER_GATE, N_PARAMS,
)
from qforge.forging import (  # noqa: F401
    setup_fragment, beta_signs, derive_beta_matrices, combine_matrices,
    energy_from_alpha_matrices, measurement_circuit, basis_change_h, pauli_expectation,
    group_labels_qubit_wise, combined_basis_label, target_states, slot_names,
    build_fragment_qop, exact_ground_state, real_gauge, schmidt_decompose_real,
    HARTREE_TO_KCAL_MOL,
)
from qforge.mitigation import (  # noqa: F401
    MitigationStrategy, RawStrategy, CDRStrategy, PECStrategy,
    fit_scale, fit_all_scales, filtered_pairs, LOW_SIGNAL_CUTOFF,
    depolarizing_weights, pec_inverse_weights, gamma_factor, apply_pauli_mixture,
)
from qforge.shot_noise import (  # noqa: F401
    shot_sample, pec_shot_sample, verify_convergence, verify_scaling_exponent, summarize,
)
from qforge.floor_test import floor_test  # noqa: F401

__all__ = [
    "build_ansatz", "fit_angles", "fit_all_targets", "transpile_fixed", "two_qubit_gate_count",
    "verify_constant_gate_count", "GATE_TRANSPILE_OPT_LEVEL", "P2_PER_GATE", "P1_PER_GATE", "N_PARAMS",
    "setup_fragment", "beta_signs", "derive_beta_matrices", "combine_matrices",
    "energy_from_alpha_matrices", "measurement_circuit", "basis_change_h", "pauli_expectation",
    "group_labels_qubit_wise", "combined_basis_label", "target_states", "slot_names",
    "build_fragment_qop", "exact_ground_state", "real_gauge", "schmidt_decompose_real",
    "HARTREE_TO_KCAL_MOL",
    "MitigationStrategy", "RawStrategy", "CDRStrategy", "PECStrategy",
    "fit_scale", "fit_all_scales", "filtered_pairs", "LOW_SIGNAL_CUTOFF",
    "depolarizing_weights", "pec_inverse_weights", "gamma_factor", "apply_pauli_mixture",
    "shot_sample", "pec_shot_sample", "verify_convergence", "verify_scaling_exponent", "summarize",
    "floor_test",
]
