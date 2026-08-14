#!/usr/bin/env python3
"""
task31b_noisy_sim_validation.py -- iteration 31, Task B, the VALID
version of the unbiasedness test. My earlier attempt (pure exact-
statevector twirling, no underlying noise) was methodologically invalid:
PEC's correction mixture is only unbiased when it's undoing an ACTUAL
depolarizing channel of the same p -- testing the correction with no
noise present to correct isn't a meaningful test of anything.
============================================================================
THE VALID TEST: simulate REAL depolarizing noise via density-matrix
propagation (`loop_pec.py`'s own `apply_pauli_mixture` +
`depolarizing_weights`, matching what real IonQ execution would produce
if the calibrated p2/p1 accurately describe the channel), for BOTH (a)
the raw circuit and (b) EVERY twirled-circuit realization (the inserted
Pauli-correction gates are real physical gates too -- they experience the
SAME noise, exactly as real hardware execution would). This directly
tests: does literal twirling, applied on top of the SAME noise model the
analytic method assumes, actually recover the exact (noiseless) value
better than raw? -- entirely local, zero new real submissions, but this
time a mathematically valid test.

Run:
    python vqe/task31b_noisy_sim_validation.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture, gamma_factor
from qiskit.quantum_info import DensityMatrix, Operator, Pauli, Statevector

K = 6
SLOT = "(u0+u1)"
GROUP = ["XYYX", "IYYI"]
LABEL = "IYYI"
GATE_NAME = "zz"
P2 = 0.0146
P1 = 0.000119
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31b_noisy_sim_validation_results.json")


def noisy_dm(qc, p2, p1):
    """Propagate qc through REAL simulated depolarizing noise (no correction)."""
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


def sample_twirled_circuit_local(base_qc, p2, p1, rng):
    """Same logic as task30b_literal_twirling's version, inlined so this
    file has no import-time dependency on that module's own GATE_NAME
    global (kept local and explicit for this test)."""
    from native_stateprep import to_native
    qc = base_qc.copy_empty_like()
    total_sign = 1
    total_gamma = 1.0
    for instr in base_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name == "gpi":
            weights = pec_inverse_weights(p1, 1)
        elif op.name == "gpi2":
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
    print("  task31b_noisy_sim_validation.py -- the VALID unbiasedness test (local, no new submission)")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    base = optimized_native_circuit(fixed_solutions[SLOT]["angles"], GATE_NAME)
    combined = effrag_mod.combined_basis_label(GROUP)
    basis_qc = native_basis_change(combined, GATE_NAME)
    full_base = base.compose(basis_qc)

    exact_val = float(np.real(Statevector.from_instruction(full_base).expectation_value(Pauli(LABEL))))
    print(f"  TRUE exact <{LABEL}> (zero noise): {exact_val:.4f}")

    # -- raw, under the SAME simulated noise model --
    dm_raw = noisy_dm(full_base, P2, P1)
    Pmat = np.asarray(Pauli(LABEL).to_matrix())
    raw_val = float(np.real(np.trace(Pmat @ dm_raw.data)))
    print(f"  RAW under simulated noise (p2={P2}, p1={P1}): {raw_val:.4f}  "
          f"(|err|={abs(raw_val-exact_val):.4f})")

    # -- analytic PEC, for reference (same noise model, gate-by-gate correction) --
    dm_A = DensityMatrix.from_label("0" * full_base.num_qubits)
    dm_B = DensityMatrix.from_label("0" * full_base.num_qubits)
    for instr in full_base.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [full_base.find_bit(q).index for q in instr.qubits]
        U = Operator(op.to_matrix())
        dm_A = dm_A.evolve(U, qargs=qargs)
        dm_B = dm_B.evolve(U, qargs=qargs)
        if op.name == GATE_NAME:
            dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(P2, 2))
            dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(P2, 2))
            dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(P2, 2))
        elif op.name in ("gpi", "gpi2"):
            dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(P1, 1))
            dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(P1, 1))
            dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(P1, 1))
    analytic_val = float(np.real(np.trace(Pmat @ dm_B.data)))
    print(f"  ANALYTIC PEC under this SAME noise model: {analytic_val:.4f}  "
          f"(|err|={abs(analytic_val-exact_val):.4f}) -- should be very close to exact by PEC's own design")

    # -- literal twirling: EVERY sampled twirled circuit experiences the SAME real noise --
    print(f"\n  -- literal twirling, correction gates ALSO experience the same simulated noise --")
    rng = np.random.default_rng(2024)
    N_MC = 400
    vals = []
    for i in range(N_MC):
        twirled, sign, gamma = sample_twirled_circuit_local(full_base, P2, P1, rng)
        dm_t = noisy_dm(twirled, P2, P1)
        m_noisy = float(np.real(np.trace(Pmat @ dm_t.data)))
        vals.append(sign * gamma * m_noisy)
        if (i + 1) in [8, 16, 32, 64, 128, 200, 400]:
            running = float(np.mean(vals[: i + 1]))
            se = float(np.std(vals[: i + 1]) / np.sqrt(i + 1))
            print(f"    N_MC={i+1:>4}: literal_twirl_estimate={running:+.4f} +/- {se:.4f}  "
                  f"|err vs exact|={abs(running-exact_val):.4f}")

    final = float(np.mean(vals))
    print(f"\n  FINAL (N_MC={N_MC}): literal={final:.4f}  analytic={analytic_val:.4f}  "
          f"raw={raw_val:.4f}  exact={exact_val:.4f}")
    print(f"  literal error vs exact: {abs(final-exact_val):.4f}")
    print(f"  analytic error vs exact: {abs(analytic_val-exact_val):.4f}")
    verdict = "literal twirling matches analytic (both correct)" if abs(final - analytic_val) < 0.05 else \
              "literal twirling and analytic DISAGREE even under the SAME local noise model -- a real bug in one of the two implementations"
    print(f"  VERDICT: {verdict}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"exact": exact_val, "raw": raw_val, "analytic": analytic_val, "literal_final": final,
                    "n_mc": N_MC, "verdict": verdict}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
