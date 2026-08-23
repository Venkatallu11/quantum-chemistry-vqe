#!/usr/bin/env python3
"""
task37a_iyyi_convention_lock.py -- iteration 37, Task A. Resolve the
IYYI sign puzzle Task 37c's real submission surfaced (sign-flipped vs
exact expectation, ON THE IDEAL/ZERO-NOISE BACKEND TOO) BEFORE it can
contaminate any calibration or PEC conclusion. Pure exact computation --
no noise, no PEC, no hardware, no shots. Compute <psi|P|psi> FOUR
independent ways for IYYI (and control labels for contrast) and check
they all agree:
  1. Statevector.expectation_value(Pauli(l)) on the unrotated ansatz.
  2. Explicit matrix contraction: psi.conj() @ P.to_matrix() @ psi,
     same state, different code path (rules out an expectation_value
     API quirk).
  3. Infinite-shot limit of the ACTUAL measurement pipeline this project
     uses for real data: basis-rotate (native_basis_change), take the
     exact probability distribution over computational-basis outcomes,
     apply expectation_from_counts' OWN sign/parity logic analytically
     (not sampled) -- this is the path most likely to diverge, since
     it's the one Task 37c's real submission actually used.
  4. The EF alpha/beta reduced-matrix path: a @ P_S[l] @ a using the
     ideal target_coeff_vector and the SAME P_S = build_P_S(...) matrix
     this project's whole PEC/manifold pipeline is built on.
If all four agree, IYYI's "wrong sign" in prior iterations (31B/32I) was
purely a hardware/coherent-noise effect, as originally diagnosed -- no
correction needed anywhere upstream. If they disagree, this pins down
EXACTLY which code path breaks the convention, before that break can
contaminate a self-consistent calibration model that would otherwise
"explain" it away as noise.

Run:
    PYTHONHASHSEED=0 python vqe/task37a_iyyi_convention_lock.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector
from phys_constrained_reconstruction import build_P_S
from qiskit.quantum_info import Statevector, Pauli

K = 6
GATE_NAME = "zz"
TEST_CASES = [
    {"slot": "(u0+u1)", "group": ["XYYX", "IYYI"], "focus": "IYYI"},
    {"slot": "(u0+u1)", "group": ["XYYX", "IYYI"], "focus": "XYYX"},  # control: same group, contrast case
    {"slot": "(u3+u5)", "group": ["XZXZ", "XZXI", "IZIZ", "IIIZ"], "focus": "IIIZ"},  # control: unrelated slot
]


def method1_statevector(base_qc, label):
    return float(np.real(Statevector.from_instruction(base_qc).expectation_value(Pauli(label))))


def method2_explicit_matrix(base_qc, label):
    psi = Statevector.from_instruction(base_qc).data
    Pmat = Pauli(label).to_matrix()
    return float(np.real(np.conj(psi) @ Pmat @ psi))


def method3_measurement_pipeline_exact(base_qc, basis_qc, group, label):
    """Exact (infinite-shot) version of the REAL measurement pipeline:
    rotate into measurement basis, take the exact outcome probability
    distribution, apply expectation_from_counts' OWN sign/parity
    convention analytically (sum over ALL 2^n outcomes weighted by
    their exact probability and the label's +-1 parity on that bitstring
    restricted to the label's non-identity qubits) -- exactly what
    `expectation_from_counts` computes in the N_shots->infinity limit."""
    full = base_qc.compose(basis_qc)
    sv = Statevector.from_instruction(full)
    probs = sv.probabilities_dict()  # {bitstring: probability}, Qiskit's own little-endian convention
    n = full.num_qubits
    # which qubit positions does this label act on non-trivially (matches expectation_from_counts'
    # own convention: Qiskit bitstrings are little-endian, qubit 0 = rightmost character)
    non_id_positions = [i for i, c in enumerate(reversed(label)) if c != "I"]
    total = 0.0
    for bitstring, prob in probs.items():
        bits = bitstring.zfill(n)
        parity = 1
        for pos in non_id_positions:
            bit = bits[n - 1 - pos]
            if bit == "1":
                parity *= -1
        total += parity * prob
    return float(total)


def method4_ef_reduced_matrix(slot, label, P_S):
    a = target_coeff_vector(slot, K)  # ideal K=6 coefficient vector for this slot
    P = P_S[label]
    return float(np.real(a @ P @ a))


def main():
    print("\n" + "=" * 96)
    print("  task37a_iyyi_convention_lock.py -- four independent <psi|P|psi> computations, must agree")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- results may not be reproducible (Task 36's known issue).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    all_pass = True
    for case in TEST_CASES:
        slot, group, label = case["slot"], case["group"], case["focus"]
        base = optimized_native_circuit(fixed_solutions[slot]["angles"], GATE_NAME)
        basis_qc = native_basis_change(effrag_mod.combined_basis_label(group), GATE_NAME)

        m1 = method1_statevector(base, label)
        m2 = method2_explicit_matrix(base, label)
        m3 = method3_measurement_pipeline_exact(base, basis_qc, group, label)
        m4 = method4_ef_reduced_matrix(slot, label, P_S)

        vals = {"1_statevector": m1, "2_explicit_matrix": m2, "3_measurement_pipeline": m3, "4_ef_reduced": m4}
        spread = max(vals.values()) - min(vals.values())
        agree = spread < 1e-6
        all_pass &= agree

        print(f"\n  === slot={slot}  label={label} ===")
        for name, v in vals.items():
            print(f"    {name:<24} {v:+.6f}")
        print(f"    spread = {spread:.6f}   {'AGREE (PASS)' if agree else 'DISAGREE -- REAL CONVENTION BUG'}")
        if not agree:
            # pinpoint which pair(s) actually diverge
            names = list(vals.keys())
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    d = abs(vals[names[i]] - vals[names[j]])
                    if d > 1e-6:
                        print(f"      DIVERGES: {names[i]} vs {names[j]}  diff={d:.6f}")

    print(f"\n  -- OVERALL: {'ALL CASES AGREE across all 4 methods -- no convention bug found here' if all_pass else 'AT LEAST ONE CASE DISAGREES -- see above for exactly which method(s) diverge'} --")


if __name__ == "__main__":
    main()
