#!/usr/bin/env python3
"""
qforge.mitigation — raw / CDR (per-basis) / PEC (learned-channel) as
interchangeable strategies sharing one interface: `correct(raw, ...) ->
alpha_matrices`, consumed by qforge.forging.energy_from_alpha_matrices.
Extracted from rank6_symmetry_vd.py (CDR fitting) and loop_pec.py /
loop_pauli_lindblad_pec.py (PEC) -- physics unchanged.

INVARIANT ENFORCED IN CODE: LOW_SIGNAL_CUTOFF is not a per-call default a
caller can silently omit -- CDRStrategy's constructor requires it
explicitly have a value (defaults to the project-established 0.05, but
is a real constructor argument, not buried inside fit logic) and
filtered_pairs() enforces the drop unconditionally, every call, not just
when convenient. Every training pair with |exact| below the cutoff is
excluded from every fit this module performs -- there is no code path
that fits a CDR scale to a near-zero-signal pair.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from qforge.forging import combine_matrices, slot_names
import loop_pec as _pec

LOW_SIGNAL_CUTOFF = 0.05  # invariant: established across every iteration in RESEARCH_LEDGER.md

depolarizing_weights = _pec.depolarizing_weights
pec_inverse_weights = _pec.pec_inverse_weights
gamma_factor = _pec.gamma_factor
apply_pauli_mixture = _pec.apply_pauli_mixture


class MitigationStrategy:
    """Common interface every method in this module implements: turn
    per-slot raw label measurements into corrected alpha matrices, ready
    for qforge.forging.energy_from_alpha_matrices. Nothing about this
    interface assumes local exact simulation vs real shot-based
    measurement -- `raw` is just {slot_name: {label: float}}, produced
    however the caller measured it."""
    name = "base"

    def correct(self, raw, alpha_labels, identity_label, K):
        raise NotImplementedError


class RawStrategy(MitigationStrategy):
    """No correction -- the floor every mitigation method is measured
    against."""
    name = "raw"

    def correct(self, raw, alpha_labels, identity_label, K):
        return combine_matrices(raw, alpha_labels, identity_label, K)


def filtered_pairs(training, label=None, slot=None, low_signal_cutoff=LOW_SIGNAL_CUTOFF):
    out = []
    for row in training:
        if label is not None and row["label"] != label:
            continue
        if slot is not None and row["slot"] != slot:
            continue
        if abs(row["exact"]) < low_signal_cutoff:
            continue  # invariant: low-signal training pairs are never fit
        out.append((row["exact"], row["noisy"]))
    return out


def fit_scale(pairs):
    """exact = scale * noisy, least-squares through the origin (verified
    across every CDR variant in this project: no real intercept exists
    to fit -- see RESEARCH_LEDGER.md iteration 1)."""
    if len(pairs) < 3:
        return None
    exact = np.array([e for e, _ in pairs])
    noisy = np.array([n for _, n in pairs])
    denom = float(np.sum(exact * exact))
    if denom < 1e-12:
        return None
    return float(np.sum(exact * noisy) / denom)


def fit_all_scales(training, non_id_labels, K, low_signal_cutoff=LOW_SIGNAL_CUTOFF):
    global_scale = fit_scale(filtered_pairs(training, low_signal_cutoff=low_signal_cutoff))
    per_basis, per_basis_fallback = {}, []
    for l in non_id_labels:
        f = fit_scale(filtered_pairs(training, label=l, low_signal_cutoff=low_signal_cutoff))
        if f is None:
            f, per_basis_fallback = global_scale, per_basis_fallback + [l]
        per_basis[l] = f
    per_circuit, per_circuit_fallback = {}, []
    for name in slot_names(K):
        f = fit_scale(filtered_pairs(training, slot=name, low_signal_cutoff=low_signal_cutoff))
        if f is None:
            f, per_circuit_fallback = global_scale, per_circuit_fallback + [name]
        per_circuit[name] = f
    return {"global_scale": global_scale, "per_basis_scale": per_basis, "per_circuit_scale": per_circuit,
            "per_basis_fallback": per_basis_fallback, "per_circuit_fallback": per_circuit_fallback}


class CDRStrategy(MitigationStrategy):
    """Clifford Data Regression, per-basis scale (the best surviving-
    scrutiny noiseless-estimator result in this project, 2.850+/-0.490
    kcal/mol locally -- see RESEARCH_LEDGER.md; also the method iteration
    9 found makes things 2.1-2.6x WORSE on real IonQ noise, reported
    honestly there, not hidden here)."""
    name = "cdr"

    def __init__(self, training, non_id_labels, K, low_signal_cutoff=LOW_SIGNAL_CUTOFF):
        self.scales = fit_all_scales(training, non_id_labels, K, low_signal_cutoff)

    def correct(self, raw, alpha_labels, identity_label, K):
        return combine_matrices(raw, alpha_labels, identity_label, K,
                                 per_label_scale=self.scales["per_basis_scale"])


class PECStrategy(MitigationStrategy):
    """Probabilistic error cancellation on a LEARNED channel. `raw` here
    is expected to already be gamma-weighted, sign-corrected expectation
    values (from either the exact density-matrix shortcut this project's
    local simulator uses, or the real randomized-quasi-probability-
    circuit protocol iteration 9 used against real IonQ hardware) --
    PECStrategy.correct() is then just combine_matrices() with no further
    scale, since the correction already happened at measurement time.
    This mirrors the fact that PEC corrects noise WHERE IT HAPPENS
    (per-gate), not as an end-of-circuit rescale like CDR."""
    name = "pec"

    def __init__(self, learned_p2, learned_p1, gamma_total=None, n2q_gates=11, n1q_gates=51):
        self.learned_p2 = learned_p2
        self.learned_p1 = learned_p1
        if gamma_total is None:
            gamma_total = gamma_factor(max(learned_p2, 1e-12), 2) ** n2q_gates * \
                          gamma_factor(max(learned_p1, 1e-12), 1) ** n1q_gates
        self.gamma_total = gamma_total

    def correct(self, raw, alpha_labels, identity_label, K):
        return combine_matrices(raw, alpha_labels, identity_label, K)
