#!/usr/bin/env python3
"""
task34a2_gate_aligned_crn.py -- iteration 34, Task A (second attempt).
Task 34A's naive "shared RNG seed across folds" gave essentially ZERO
correlation (rho=0.02-0.07 for fold 1v3/1v5, RED per the proposal's own
thresholds). Before concluding the whole CRN idea is dead, check whether
that was a REAL negative result or an ARTIFACT of how the seed was
shared: `fold_all_gates` (Task 28d) inserts fold-repeat gates INLINE,
immediately after each original gate, not appended at the end. This means
a single sequentially-advancing RNG stream visits gates in a completely
different ORDER at fold=3 (gate1, gate1_repeat, gate1_repeat, gate2, ...)
than at fold=1 (gate1, gate2, gate3, ...) -- so "same seed" does NOT mean
"same twirl choice for the same physical gate" once folding changes the
circuit structure. That's a plausible, checkable, FIXABLE bug in the
naive construction, not evidence against the underlying physical idea.

THE FIX, a fundamentally different (and arguably more principled)
construction: instead of sharing an RNG STREAM, share the TWIRL DECISION
itself, keyed by GATE IDENTITY. Every fold-repeat copy of "the same
original gate" (fold_all_gates always keeps the ORIGINAL op as the first
copy, then appends repeat-echoes) gets the SAME PEC recovery-Pauli choice
as that gate's other copies, at every fold level and every repeat
instance -- not independently re-rolled. Gates introduced PURELY by the
folding construction itself (the intermediate GPI(0) gates used to build
the 2Q-gate inverse) have no fold=1 analog and get their own reproducible
but fold-local draws.

VALIDATION BEFORE TRUSTING: `build_folded_with_parent_tags` reimplements
`fold_all_gates`'s exact gate-insertion logic (same gate types/qargs/
params in the same order) so it can additionally tag each instruction
with its original-gate parent -- this is checked instruction-for-
instruction against the real `fold_all_gates` output before any energy
number is trusted, exactly this project's own house rule (never trust a
new local reimplementation without checking it against the established
one first).

Run:
    python vqe/task34a2_gate_aligned_crn.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit, fold_all_gates
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from native_stateprep import to_native
from ionq_simulator_binding_curve import stable_seed
from qiskit.quantum_info import DensityMatrix, Operator, Pauli, Statevector
from qiskit_ionq.ionq_gates import GPIGate, GPI2Gate

K = 6
SLOT = "(u0+u1)"
GROUP = ["XYYX", "IYYI"]
LABEL = "IYYI"
GATE_NAME = "zz"
P2 = 0.0146
P1 = 0.000119
FOLDS = [1, 3, 5]
M = 256
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task34a2_gate_aligned_crn_results.json")


def build_folded_with_parent_tags(qc, fold, gate_name):
    """Mirrors task28d's fold_all_gates EXACTLY (same gate objects, same
    order), additionally returning a parallel list of parent tags:
    ('primary', orig_idx) for direct copies of an original gate,
    ('extra', orig_idx, k) for a fold-construction-only gate (unique
    per repeat k, not present at fold=1)."""
    folded = qc.copy_empty_like()
    tags = []
    if fold == 1:
        for i, instr in enumerate(qc.data):
            folded.append(instr.operation, instr.qubits, instr.clbits)
            tags.append(("primary", i))
        return folded, tags
    assert fold % 2 == 1
    reps = (fold - 1) // 2
    for i, instr in enumerate(qc.data):
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        folded.append(op, qargs, cargs)
        tags.append(("primary", i))
        if op.name == "gpi":
            for k in range(reps):
                folded.append(op, qargs, cargs)
                tags.append(("extra", i, k, 0))
                folded.append(op, qargs, cargs)
                tags.append(("extra", i, k, 1))
        elif op.name == "gpi2":
            phi = float(op.params[0])
            inv = GPI2Gate(phi + 0.5)
            for k in range(reps):
                folded.append(inv, qargs, cargs)
                tags.append(("extra", i, k, 0))
                folded.append(op, qargs, cargs)
                tags.append(("extra", i, k, 1))
        elif op.name == gate_name:
            for k in range(reps):
                folded.append(GPIGate(0), [qargs[0]])
                tags.append(("extra", i, k, 0))
                folded.append(op, qargs, cargs)
                tags.append(("extra", i, k, 1))
                folded.append(GPIGate(0), [qargs[0]])
                tags.append(("extra", i, k, 2))
                folded.append(op, qargs, cargs)
                tags.append(("extra", i, k, 3))
    return folded, tags


def circuits_identical(qc_a, qc_b):
    if len(qc_a.data) != len(qc_b.data):
        return False, f"length mismatch {len(qc_a.data)} vs {len(qc_b.data)}"
    for i, (ia, ib) in enumerate(zip(qc_a.data, qc_b.data)):
        if ia.operation.name != ib.operation.name:
            return False, f"instr {i}: name {ia.operation.name} vs {ib.operation.name}"
        if list(ia.operation.params) != list(ib.operation.params):
            return False, f"instr {i}: params {ia.operation.params} vs {ib.operation.params}"
        if [qc_a.find_bit(q).index for q in ia.qubits] != [qc_b.find_bit(q).index for q in ib.qubits]:
            return False, f"instr {i}: qargs mismatch"
    return True, "identical"


def noisy_dm(qc, p2, p1):
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
        if op.name == GATE_NAME:
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p2, 2))
        elif op.name in ("gpi", "gpi2"):
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p1, 1))
    return dm


def gate_aligned_twirl(folded_qc, tags, p2, p1, draw_k):
    """Every ('primary', i) instruction's twirl choice is keyed ONLY by
    (SLOT, LABEL, draw_k, i) -- identical across every fold level, since
    fold_all_gates always keeps the original gate as the first copy.
    Every ('extra', i, k, sub) instruction is keyed by its own full tag
    -- unique, reproducible, but fold-local (fold=5's second repeat has
    no fold=3 analog and is not forced to correlate with anything)."""
    qc = folded_qc.copy_empty_like()
    total_sign = 1
    total_gamma = 1.0
    for instr, tag in zip(folded_qc.data, tags):
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name in ("gpi", "gpi2"):
            weights = pec_inverse_weights(p1, 1)
        else:
            continue
        seed = stable_seed("task34a2", SLOT, LABEL, draw_k, *tag)
        rng_local = np.random.default_rng(seed)
        labels = list(weights.keys())
        w = np.array([weights[l] for l in labels])
        gamma_gate = float(np.sum(np.abs(w)))
        probs = np.abs(w) / gamma_gate
        idx = rng_local.choice(len(labels), p=probs)
        chosen_label, chosen_w = labels[idx], w[idx]
        sign = 1 if chosen_w >= 0 else -1
        total_sign *= sign
        total_gamma *= gamma_gate
        n = len(chosen_label)
        if chosen_label != "I" * n:
            qc.append(Pauli(chosen_label).to_instruction(), qargs)
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign, total_gamma


def independent_twirl(folded_qc, p2, p1, seed_base):
    """Status-quo control: every instruction gets its own independent draw."""
    rng = np.random.default_rng(seed_base)
    qc = folded_qc.copy_empty_like()
    total_sign = 1
    total_gamma = 1.0
    for instr in folded_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name in ("gpi", "gpi2"):
            weights = pec_inverse_weights(p1, 1)
        else:
            continue
        labels = list(weights.keys())
        w = np.array([weights[l] for l in labels])
        gamma_gate = float(np.sum(np.abs(w)))
        probs = np.abs(w) / gamma_gate
        idx = rng.choice(len(labels), p=probs)
        chosen_label, chosen_w = labels[idx], w[idx]
        sign = 1 if chosen_w >= 0 else -1
        total_sign *= sign
        total_gamma *= gamma_gate
        n = len(chosen_label)
        if chosen_label != "I" * n:
            qc.append(Pauli(chosen_label).to_instruction(), qargs)
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign, total_gamma


def main():
    print("\n" + "=" * 96)
    print("  task34a2_gate_aligned_crn.py -- gate-identity-keyed CRN (second attempt, LOCAL, zero shots)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    base = optimized_native_circuit(fixed_solutions[SLOT]["angles"], GATE_NAME)
    combined = effrag_mod.combined_basis_label(GROUP)
    basis_qc = native_basis_change(combined, GATE_NAME)
    Pmat = np.asarray(Pauli(LABEL).to_matrix())

    # -- VALIDATION: my tagged reimplementation must produce IDENTICAL circuits to the established fold_all_gates.
    # IMPORTANT: fold_all_gates is applied to the ANSATZ ONLY, before composing the basis-rotation gates --
    # confirmed against task28d_all_gate_zne.py's and task28d_resume.py's own usage (`fold_all_gates(base, ...)`,
    # never the basis-composed circuit). Folding the basis-rotation gates too (an earlier draft of this script)
    # would test a DIFFERENT, non-standard construction -- fixed before running anything expensive. --
    print("\n  -- VALIDATION: build_folded_with_parent_tags vs task28d's fold_all_gates, must be gate-identical --")
    full_folded, tags_by_fold = {}, {}
    for fold in FOLDS:
        folded_ref = fold_all_gates(base, fold, GATE_NAME).compose(basis_qc)
        folded_ansatz_tagged, tags = build_folded_with_parent_tags(base, fold, GATE_NAME)
        ok, msg = circuits_identical(fold_all_gates(base, fold, GATE_NAME), folded_ansatz_tagged)
        print(f"    fold={fold}: {msg} -- {'PASS' if ok else 'FAIL, STOP'}")
        if not ok:
            raise RuntimeError(f"Tagged folding diverges from the established fold_all_gates at fold={fold}: {msg}")
        # basis-rotation gates appended AFTER folding, untagged (not part of the ansatz being fold-scaled) --
        # each gets its own independent twirl draw, same as any other untagged instruction downstream
        n_basis_gates = len(basis_qc.data)
        full_folded[fold] = folded_ansatz_tagged.compose(basis_qc)
        tags_by_fold[fold] = tags + [("basis", i) for i in range(n_basis_gates)]
        assert circuits_identical(full_folded[fold], folded_ref)[0]

    exact_val = float(np.real(Statevector.from_instruction(full_folded[1]).expectation_value(Pauli(LABEL))))
    print(f"\n  TRUE exact <{LABEL}> (zero noise): {exact_val:.4f}")

    print(f"\n  -- generating {M} draws per fold, GATE-ALIGNED (shared twirl decision per original gate) "
          f"and INDEPENDENT (control) --")
    E_aligned = {fold: np.zeros(M) for fold in FOLDS}
    E_indep = {fold: np.zeros(M) for fold in FOLDS}
    for k in range(M):
        for fold in FOLDS:
            twirled, sign, gamma = gate_aligned_twirl(full_folded[fold], tags_by_fold[fold], P2, P1, k)
            dm_t = noisy_dm(twirled, P2, P1)
            m = float(np.real(np.trace(Pmat @ dm_t.data)))
            E_aligned[fold][k] = sign * gamma * m
        for fold in FOLDS:
            seed = stable_seed("task34a2_indep", SLOT, LABEL, fold, k)
            twirled, sign, gamma = independent_twirl(full_folded[fold], P2, P1, seed)
            dm_t = noisy_dm(twirled, P2, P1)
            m = float(np.real(np.trace(Pmat @ dm_t.data)))
            E_indep[fold][k] = sign * gamma * m
        if (k + 1) % 64 == 0:
            print(f"    {k+1}/{M} draws done")

    def corr(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    print(f"\n  -- CORRELATION: GATE-ALIGNED CRN vs INDEPENDENT (control) --")
    results = {}
    for fold in [3, 5]:
        rho_a = corr(E_aligned[1], E_aligned[fold])
        rho_i = corr(E_indep[1], E_indep[fold])
        var1_a, varf_a = E_aligned[1].var(ddof=1), E_aligned[fold].var(ddof=1)
        var_diff_a = (E_aligned[1] - E_aligned[fold]).var(ddof=1)
        ratio_a = var_diff_a / (var1_a + varf_a)
        print(f"    fold 1 vs {fold}:")
        print(f"      GATE-ALIGNED: rho={rho_a:+.4f}   Var(E1-E{fold})/[Var(E1)+Var(E{fold})]={ratio_a:.4f}")
        print(f"      INDEPENDENT:  rho={rho_i:+.4f}  (control baseline)")
        print(f"      CRN GAP (aligned - independent) in rho: {rho_a - rho_i:+.4f}")
        results[f"1_vs_{fold}"] = {
            "rho_aligned": rho_a, "rho_independent": rho_i, "crn_gap": rho_a - rho_i,
            "var_ratio_aligned": ratio_a,
        }

    print(f"\n  -- VERDICT (thresholds: green rho>0.9, yellow 0.6-0.9, red rho<0.3) --")
    for fold in [3, 5]:
        rho = results[f"1_vs_{fold}"]["rho_aligned"]
        verdict = "GREEN" if rho > 0.9 else ("YELLOW" if rho >= 0.6 else ("RED" if rho < 0.3 else "borderline"))
        print(f"    fold 1 vs {fold}: rho_aligned={rho:+.4f} -> {verdict}")

    # -- THE BOTTOM-LINE QUESTION: does the confirmed correlation actually tighten a real zero-noise
    # extrapolation, or does it get swallowed the way Task 33E's confirmed-but-real fix did? Simplest
    # possible test: 2-point linear Richardson on folds (1,3), E0_est = 2*E1 - E3 (exact for a noise
    # response linear in fold number), applied per-draw k, compared MATCHED (CRN) vs INDEPENDENT --
    # against the TRUE exact value we already know (no noise), not an assumption. --
    print(f"\n  -- BOTTOM-LINE TEST: does CRN actually tighten the zero-noise (fold->0) extrapolation? --")
    print(f"     2-point linear Richardson, folds x=(1,3): solving E1=E0+b, E3=E0+3b gives "
          f"E0_est = 1.5*E1 - 0.5*E3 (CORRECTED -- an earlier draft of this script used the wrong "
          f"formula 2*E1-E3, caught by hand-deriving the 2-point linear system before trusting the "
          f"result: that bug alone produced sign-flipped, ~2x-magnitude-wrong extrapolated values, "
          f"a math error, not a finding about CRN), per draw, {M} draws")
    E0_aligned = 1.5 * E_aligned[1] - 0.5 * E_aligned[3]
    E0_indep = 1.5 * E_indep[1] - 0.5 * E_indep[3]
    bias_aligned = float(E0_aligned.mean() - exact_val)
    bias_indep = float(E0_indep.mean() - exact_val)
    std_aligned = float(E0_aligned.std(ddof=1))
    std_indep = float(E0_indep.std(ddof=1))
    std_raw_e1 = float(E_aligned[1].std(ddof=1))  # same draws feed both conditions' E1 population-wise, for reference
    print(f"     RAW fold-1 alone (no extrapolation, reference): std={std_raw_e1:.4f}")
    print(f"     ZNE, INDEPENDENT draws: E0_mean={E0_indep.mean():.4f}  bias={bias_indep:+.4f}  std={std_indep:.4f}")
    print(f"     ZNE, MATCHED (CRN) draws: E0_mean={E0_aligned.mean():.4f}  bias={bias_aligned:+.4f}  std={std_aligned:.4f}")
    if std_indep > 0:
        std_reduction_pct = 100 * (std_indep - std_aligned) / std_indep
        print(f"     CRN std reduction vs independent ZNE: {std_reduction_pct:+.1f}%")
    else:
        std_reduction_pct = None
    zne_vs_raw_pct = 100 * (std_raw_e1 - std_aligned) / std_raw_e1 if std_raw_e1 > 0 else None
    print(f"     CRN-ZNE std vs raw-fold-1-alone: {zne_vs_raw_pct:+.1f}% "
          f"({'ZNE is TIGHTER than just using fold 1' if (zne_vs_raw_pct or 0) > 0 else 'ZNE is WORSE than just using fold 1 -- extrapolation amplified noise, exactly the old per-Pauli failure mode, now at the scalar level'})")
    print(f"\n     HONEST READ: bias must also stay reasonable -- |bias_aligned|={abs(bias_aligned):.4f} vs "
          f"|bias_indep|={abs(bias_indep):.4f} (both should be small if 2-point linear Richardson is a "
          f"reasonable local model; a big bias here means the fold-1-to-3 response isn't well-approximated "
          f"as linear for this label, a separate problem from variance).")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "slot": SLOT, "label": LABEL, "folds": FOLDS, "M": M, "exact_val": exact_val,
            "fold_validation": "build_folded_with_parent_tags verified gate-identical to fold_all_gates for all folds",
            "correlations": results,
            "raw_draws": {
                "E_aligned": {str(f): E_aligned[f].tolist() for f in FOLDS},
                "E_indep": {str(f): E_indep[f].tolist() for f in FOLDS},
            },
            "bottom_line_zne_test": {
                "method": "2-point linear Richardson on folds (1,3): E0_est = 2*E1 - E3",
                "std_raw_fold1_alone": std_raw_e1,
                "bias_aligned": bias_aligned, "std_aligned": std_aligned,
                "bias_independent": bias_indep, "std_independent": std_indep,
                "crn_std_reduction_pct": std_reduction_pct,
                "zne_vs_raw_fold1_pct": zne_vs_raw_pct,
            },
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
