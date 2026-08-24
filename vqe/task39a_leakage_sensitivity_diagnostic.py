#!/usr/bin/env python3
"""
task39a_leakage_sensitivity_diagnostic.py -- iteration 39, Task A. Before
building any new circuit (this project already has a real, hardware-
validated particle-number-leakage detector -- `spin_leakage_postselect_
ionq.py`'s ancilla-parity trick, iteration 18: real IonQ postselection
cut error from 34.98/43.03 to 31.77/33.86 kcal/mol on the RAW ansatz,
BEFORE PEC or the joint Schmidt frame existed, and has never been
combined with either), test the actual premise directly: does the
calibration-sensitive noise Task 38C found dominant (p_gpi2 75.1%,
p_readout 24.8% of V_cal) actually show up as LEAKAGE out of the
register's exact weight-2 physical sector (H4's own real symmetry,
verified exactly in `rank6_symmetry_vd.py`: Schmidt rank exactly 6,
living entirely in the 6-dim weight-2 x 6-dim weight-2 subspace) -- or
does it stay hidden WITHIN that sector, where a parity check cannot see
it at all? This is answerable locally, cheaply, with existing verified
machinery (Task 37C's own `biased_zz_matrix`/`biased_gpi2_matrix` +
`loop_pec`'s Pauli-mixture primitives), before committing to any new
circuit design or real submission.

METHOD: for each kept slot's native state-prep circuit (BEFORE basis
rotation, matching iteration 18's own ancilla-must-precede-rotation
requirement), build the noisy density matrix under a perturbation in
EACH of the 5 theta directions (one at a time, same abs_step convention
as Task 38C/38D, learned from the rel_step bug caught there), and
compute the WEIGHT-2 POPULATION FRACTION -- the probability a parity
check would read "even" (weight in {0,2,4}) and accept the shot. The
IDEAL circuit gives exactly 1.0 (verified, not assumed -- iteration 18's
own `verify_ancilla_scheme` established this). The SENSITIVITY of this
fraction to each theta direction is a direct, real answer to "would a
parity check actually catch this error" -- compared side by side with
Task 38C's own energy-sensitivity gradient g, to see whether the
calibration-DANGEROUS directions and the LEAKAGE-DETECTABLE directions
are the same directions or different ones.

READOUT is treated separately (analytic, not new machinery): a
symmetric per-qubit bit-flip readout channel with probability p_ro
produces a DETECTABLE (odd total weight) outcome with probability
1-(1-2p_ro)^n/2-ish only for an ODD number of flips among the n=4
register qubits -- computed exactly below, not assumed, since this
determines whether the SAME parity check that helps with GPi2 also
partially helps with readout error.

Run:
    PYTHONHASHSEED=0 python vqe/task39a_leakage_sensitivity_diagnostic.py
"""
import os
import sys
import numpy as np
from itertools import product

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, apply_pauli_mixture
from qiskit.quantum_info import DensityMatrix, Operator

K = 6
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
# same real, disclosed step sizes used throughout Task 38 (abs_step, never rel_step -- Task 37D's bug)
STEPS = {"p_zz": 2e-4, "p_gpi2": 2e-4, "delta_zz": 2e-3, "delta_gpi2": 2e-3}


def weight2_fraction(dm):
    diag = np.real(np.diag(dm.data))
    n = int(np.log2(len(diag)))
    total_w2 = sum(diag[i] for i in range(len(diag)) if bin(i).count("1") == 2)
    return float(total_w2)


def noisy_dm(angles, gate_name, p_zz, p_gpi2, delta_zz, delta_gpi2):
    from fixed_ansatz import build_ansatz
    from native_stateprep import to_native
    qc = to_native(build_ansatz(angles), gate_name)
    n = qc.num_qubits
    dm = DensityMatrix.from_label("0" * n)
    for instr in qc.data:
        op = instr.operation
        if op.name in ("measure", "barrier"):
            continue
        qargs = [qc.find_bit(q).index for q in instr.qubits]
        if op.name == gate_name:
            theta = float(op.params[0])
            U = Operator(biased_zz_matrix(theta, delta_zz))
            p_here, n_here = p_zz, 2
        elif op.name == "gpi":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi_matrix(phi, 0.0))
            p_here, n_here = GPI_REAL_MEAN, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2))
            p_here, n_here = p_gpi2, 1
        else:
            dm = dm.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm = dm.evolve(U, qargs=qargs)
        dm = apply_pauli_mixture(dm, qargs, depolarizing_weights(p_here, n_here))
    return dm


def readout_detection_rate(p_ro, n=4):
    """Exact: probability a symmetric per-qubit bit-flip channel with rate
    p_ro produces an ODD number of flips among n qubits (detected by the
    same even/odd parity check) -- direct combinatorial sum, not assumed."""
    p_odd = 0.0
    for k in range(1, n + 1, 2):
        from math import comb
        p_odd += comb(n, k) * (p_ro ** k) * ((1 - p_ro) ** (n - k))
    return p_odd


def main():
    print("\n" + "=" * 96)
    print("  task39a_leakage_sensitivity_diagnostic.py -- does the parity check see the dangerous noise?")
    print("=" * 96)
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36).")

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    diag, plus, kept = kept_slots_for_K(K)

    theta0 = {"p_zz": ZZ_ASSUMED, "p_gpi2": 0.0, "delta_zz": 0.0, "delta_gpi2": 0.0}

    print(f"\n  -- baseline weight-2 population fraction at theta_0 (should be ~1.0 exactly for GPi/ZZ, "
          f"non-trivial for GPi2 since it's genuinely noisy at theta_0 too) --")
    frac0_by_slot = {}
    for name in kept:
        dm0 = noisy_dm(fixed_solutions[name]["angles"], GATE_NAME, **theta0)
        frac0_by_slot[name] = weight2_fraction(dm0)
    mean_frac0 = float(np.mean(list(frac0_by_slot.values())))
    print(f"    mean weight-2 fraction across {len(kept)} slots at theta_0: {mean_frac0:.6f}")

    print(f"\n  -- leakage sensitivity: d(weight-2 fraction)/d(theta_i), central difference, all {len(kept)} slots, mean --")
    sensitivities = {}
    for name_perturb, step in STEPS.items():
        d_fracs = []
        for name in kept:
            theta_plus = dict(theta0); theta_plus[name_perturb] += step
            theta_minus = dict(theta0); theta_minus[name_perturb] -= step
            f_plus = weight2_fraction(noisy_dm(fixed_solutions[name]["angles"], GATE_NAME, **theta_plus))
            f_minus = weight2_fraction(noisy_dm(fixed_solutions[name]["angles"], GATE_NAME, **theta_minus))
            d_fracs.append((f_plus - f_minus) / (2 * step))
        sensitivities[name_perturb] = float(np.mean(d_fracs))

    print(f"\n  {'param':<12}{'d(weight-2 frac)/d(theta)':>28}{'|sensitivity|':>16}")
    for name, s in sensitivities.items():
        print(f"  {name:<12}{s:>28.4f}{abs(s):>16.4f}")

    print(f"\n  -- readout: exact ODD-flip detection rate (n=4 register qubits) at a few real p_readout values --")
    for p_ro in [0.001, 0.002, 0.005, 0.01, 0.02]:
        rate = readout_detection_rate(p_ro)
        print(f"    p_readout={p_ro:.3f}  P(detected as leakage) = {rate:.4f}  "
              f"({'meaningful fraction caught' if rate > 0.01 else 'negligible catch rate'})")

    print(f"\n  -- HONEST READ -- compare against Task 38C's own energy-sensitivity finding "
          f"(p_gpi2 75.1% of V_cal, p_readout 24.8%, p_zz/delta_zz/delta_gpi2 each ~0%) --")
    dominant_energy_dirs = {"p_gpi2"}
    max_leak_dir = max(sensitivities, key=lambda k: abs(sensitivities[k]))
    print(f"    Direction with the LARGEST leakage sensitivity: {max_leak_dir} (|d(frac)/d(theta)|={abs(sensitivities[max_leak_dir]):.4f})")
    if max_leak_dir in dominant_energy_dirs and abs(sensitivities["p_gpi2"]) > 0.01:
        print(f"    p_gpi2 IS both the dominant energy-sensitivity direction AND shows real leakage sensitivity -- "
              f"a parity check has a real chance of catching exactly the error that matters most. Worth building "
              f"the real combined circuit test.")
    else:
        print(f"    p_gpi2 (the dominant ENERGY-sensitivity direction, 75.1% of V_cal) does NOT show the largest "
              f"(or a meaningfully large) LEAKAGE sensitivity here -- if this holds, it means the calibration-"
              f"dangerous GPi2 error mostly stays WITHIN the weight-2 sector, where a parity check structurally "
              f"cannot see it, regardless of how well the circuit/detection machinery is engineered. This would be "
              f"a real, decisive reason NOT to invest further in the symmetry-detection direction for THIS specific "
              f"noise source, before any new circuit is built or any real submission is spent.")


if __name__ == "__main__":
    main()
