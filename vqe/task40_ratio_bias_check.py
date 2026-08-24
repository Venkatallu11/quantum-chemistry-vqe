#!/usr/bin/env python3
"""
task40_ratio_bias_check.py -- iteration 40, ratio-estimator bias check
(user's own point 11: E[A/B] != E[A]/E[B] for ratio estimators -- does
the finite-shot postselected estimator have a bias that vanishes as
N->infinity, or a hidden plateau?).

METHOD: build the EXACT 5-qubit noisy density matrix (register +
ancilla, real calibration: ZZ_ASSUMED, GPI_REAL_MEAN, p_gpi2=0.0006) for
one representative (slot, label) case, get its EXACT probability
distribution over 5-qubit bitstrings. Sample MULTINOMIAL counts at
increasing shot counts N (1e3, 1e4, 2e4 [this task's real shot budget],
1e5, 1e6), postselect on ancilla=0 exactly as the real pipeline does,
compute the resulting label expectation, and compare its FINITE-SHOT
estimate against the TRUE (infinite-shot / exact) conditioned
expectation value. Repeat each N with 50 independent multinomial draws
to separate genuine bias (systematic offset of the mean) from ordinary
shot noise (spread around the mean).

Run:
    PYTHONHASHSEED=0 python vqe/task40_ratio_bias_check.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task39b_native_ancilla_parity import ancilla_cnots_abstract
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, apply_pauli_mixture
from ionq_simulator_binding_curve import expectation_from_counts
from qiskit.quantum_info import DensityMatrix, Operator, partial_trace

GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
P_GPI2 = 0.0006
LABEL = "IYYI"  # this project's own long-established validation case
N_LIST = [1000, 10000, 20000, 100000, 1000000]
N_TRIALS = 50


def build_4q_dm(angles, p_zz, p_gpi2, apply_noise=True):
    """4-qubit noisy (or, with apply_noise=False, exactly noiseless)
    state-prep density matrix. GPi's OWN real (tiny) depolarizing noise
    (GPI_REAL_MEAN) is applied unconditionally whenever apply_noise=True
    -- disclosed explicitly, since it means a 'p_zz=0, p_gpi2=0' call is
    NOT the same as a genuinely zero-noise circuit unless apply_noise is
    also set False."""
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native
    qc = to_native(build_ansatz(angles), GATE_NAME)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == GATE_NAME:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, 0.0))
            p_here, n_here = p_zz, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))
            p_here, n_here = GPI_REAL_MEAN, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, 0.0))
            p_here, n_here = p_gpi2, 1
        else:
            dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm = dm.evolve(U, qargs=qargs)
        if apply_noise:
            dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p_here, n_here))
    return dm


def append_ancilla(dm4):
    """Exact 5-qubit joint state: dm4 (register) tensor |0><0| (ancilla,
    placed at the HIGHEST qubit index by qiskit's own .expand()) then the
    (noiseless) ancilla-parity CNOTs."""
    dm5 = DensityMatrix(dm4.data).expand(DensityMatrix.from_label("0"))
    ancilla_qc = ancilla_cnots_abstract()
    for instr in ancilla_qc.data:
        op = instr.operation
        qargs = [ancilla_qc.find_bit(q).index for q in instr.qubits]
        dm5 = dm5.evolve(Operator(op.to_matrix()), qargs=qargs)
    return dm5


def build_true_5q_dm(angles, p_zz, p_gpi2):
    return append_ancilla(build_4q_dm(angles, p_zz, p_gpi2, apply_noise=True))


def main():
    print("\n" + "=" * 96)
    print("  task40_ratio_bias_check.py -- finite-shot postselection bias, does it vanish as N->infinity?")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    K = 6
    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, _, _ = fit_all_targets(p["targets"], tol=1e-10)
    angles = fixed_solutions["u_0"]["angles"]

    print("\n  -- sanity check 1: expand()+CNOT marginal-state preservation, using the SAME (real-noise) dm --")
    dm4_real = build_4q_dm(angles, p_zz=ZZ_ASSUMED, p_gpi2=P_GPI2, apply_noise=True)
    dm5_real = append_ancilla(dm4_real)
    dm4_marginal = partial_trace(dm5_real, [4])
    marginal_err = float(np.max(np.abs(np.asarray(dm4_marginal) - dm4_real.data)))
    print(f"    marginal-state error: {marginal_err:.3e}   {'PASS' if marginal_err < 1e-9 else 'FAIL -- STOP, convention bug'}")
    if marginal_err >= 1e-9:
        raise RuntimeError("expand()/CNOT convention check failed -- do not trust the bias numbers below")

    print("\n  -- sanity check 2: TRUE zero-noise case must give ancilla=0 with certainty --")
    dm4_true_ideal = build_4q_dm(angles, p_zz=0.0, p_gpi2=0.0, apply_noise=False)
    dm5_true_ideal = append_ancilla(dm4_true_ideal)
    probs_ideal = dm5_true_ideal.probabilities_dict()
    p_anc1_ideal = sum(v for k, v in probs_ideal.items() if k[0] == "1")
    print(f"    p(ancilla=1) at TRUE zero noise: {p_anc1_ideal:.3e}   {'PASS' if p_anc1_ideal < 1e-9 else 'FAIL -- STOP, convention bug'}")
    if p_anc1_ideal >= 1e-9:
        raise RuntimeError("expand()/CNOT convention check failed -- do not trust the bias numbers below")

    dm5 = build_true_5q_dm(angles, p_zz=ZZ_ASSUMED, p_gpi2=P_GPI2)
    probs = dm5.probabilities_dict()
    bitstrings = list(probs.keys())
    prob_arr = np.array([probs[bs] for bs in bitstrings])
    prob_arr = prob_arr / prob_arr.sum()  # renormalize away any tiny numerical residual

    # TRUE (infinite-shot) conditioned expectation -- exact reference
    true_filtered = {bs[1:]: probs[bs] for bs in bitstrings if bs[0] == "0"}
    true_val = expectation_from_counts(true_filtered, LABEL)
    print(f"\n  slot=u_0 label={LABEL}  TRUE (exact, infinite-shot) conditioned expectation = {true_val:+.6f}")

    rng = np.random.default_rng(40)
    print(f"\n  {'N shots':>10}{'mean estimate':>16}{'bias (mean-true)':>20}{'std across trials':>20}")
    for N in N_LIST:
        estimates = []
        for _ in range(N_TRIALS):
            counts_idx = rng.multinomial(N, prob_arr)
            counts = {bitstrings[i]: int(c) for i, c in enumerate(counts_idx) if c > 0}
            filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == "0"}
            if sum(filtered.values()) == 0:
                continue
            estimates.append(expectation_from_counts(filtered, LABEL))
        mean_est = float(np.mean(estimates))
        bias = mean_est - true_val
        std_est = float(np.std(estimates, ddof=1))
        print(f"  {N:>10}{mean_est:>16.6f}{bias:>20.6f}{std_est:>20.6f}")

    print(f"\n  -- HONEST READ --")
    print(f"  If bias shrinks toward 0 as N grows (roughly monotonically, well below the std at each N), "
          f"this is ordinary, expected shot noise with no hidden ratio-estimator bias floor. If bias stays "
          f"roughly CONSTANT (does not shrink) even at N=1,000,000, that would indicate a real, hidden "
          f"systematic problem -- worth flagging explicitly either way.")


if __name__ == "__main__":
    main()
