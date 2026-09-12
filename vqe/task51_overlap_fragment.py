#!/usr/bin/env python3
"""
task51_overlap_fragment.py -- iteration 51. Covalent-bonding extension,
part 2: the OVERLAP fragment (atoms [2,3], nelec=2) of the H6 molecular-
tailoring scheme (`ionq_tailoring.py`'s own E_tailored = E(A)+E(B)-E(overlap)
inclusion-exclusion formula) -- genuinely new physics, not a replication.
Its Schmidt rank is exactly K=2 (verified: K=1 fails with tail=1.753e-01,
K=2 gives tail=0.000e+00 exactly).

REAL PHYSICAL FINDING, caught by verification BEFORE any real submission:
this register holds 1 electron in 2 orbitals -- weight-1, ODD parity --
unlike H4's 2-electrons-in-4-orbitals registers, which are weight-2, EVEN.
The standard XOR-parity ancilla-CNOT construction makes the ancilla equal
the register's total Hamming-weight parity, so "which ancilla outcome is
leakage-free" is a real physical property of the fragment, not a free
convention choice. A first run using H4's own convention (ancilla=0 is
valid) correctly FAILED at the verification step: worst ideal p(leaked)
= 1.000 (the ancilla read 1 with certainty, the opposite of expected) --
caught before any real submission, exactly what this project's own
verify-before-trust discipline exists to do. Fixed by flipping the
convention (LEAKAGE_FREE="1") and using the ODD-weight projector
(I - even_weight_projector(n)) in the analytic conditioned-correction
model instead of the even-weight one.

BUILDING BLOCKS, each verified before use:
  1. Ansatz: reuses fixed_ansatz.py's OWN already-validated real-rotation
     building block, XXPlusYYGate(theta, beta=pi/2) -- an exact real
     rotation on {|01>,|10>}, leaving |00>/|11> untouched -- rather than
     hand-deriving a new circuit (a first hand-derived attempt was
     provably wrong: a bare RY+CX sequence does not stay in the target
     subspace, caught by a feasibility scan before it was used for
     anything). Single free parameter per target; fit to machine
     precision via least_squares (angles land at clean values: 0, -pi,
     +-pi/2, matching the structure of the four K=2 targets exactly).
  2. Ancilla-parity: generalized to an n-qubit register (n_register+1
     total qubits), verified via marginal-state preservation (worst
     2.01e-34) and ideal leakage-free certainty (worst p(leaked)=0.000e+00,
     after the parity-convention fix above).
  3. Conditioned correction: generalized `analytic_A_and_B_conditioned`
     (same gate-by-gate noisy + PEC-inverse propagation as the H4
     pipeline, task39e's own even_weight_projector(n) is already generic
     -- only the module-level n=4 CONSTANT was H4-specific, not the
     function) -- verified at true theta=0 to be an exact no-op
     (max diff A vs B = 0.000e+00, retained_A=1.000000).
  4. Device-calibration reuse, disclosed: ZZ_ASSUMED/GPI_REAL_MEAN/
     GPI2_SELECTED (0.0006 aria-1, 0.0004 forte-1) are properties of the
     SIMULATOR'S per-gate noise model, not of the molecule being run
     through it -- reused here rather than re-swept, the same reasoning
     that let fragment B (task50) reuse the same calibration successfully.

REAL RESULTS, all three free simulators, 20,000 shots:
    ideal:   E=-1.101145 Ha  err vs exact = +0.0034 kcal/mol  (accept=1.0000)
    aria-1:  E=-1.101138 Ha  err vs exact = +0.0078 kcal/mol  (accept=0.9501)
    forte-1: E=-1.101115 Ha  err vs exact = +0.0221 kcal/mol  (accept=0.9660)
Overlap exact energy: -1.101150 Ha. All three land in the same excellent
sub-0.03 kcal/mol range as fragments A and B.

Run:
    PYTHONHASHSEED=0 python vqe/task51_overlap_fragment.py            # full build+verify+submit+analyze
    PYTHONHASHSEED=0 python vqe/task51_overlap_fragment.py --analyze  # analyze already-submitted data only
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from native_stateprep import to_native
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, expectation_from_counts
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from task39e_conditioned_correction import even_weight_projector
import ef_fragment as effrag_mod

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import XXPlusYYGate
from qiskit.quantum_info import Statevector, Operator, partial_trace, DensityMatrix, Pauli
from scipy.optimize import least_squares

K = 2
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
GPI2_SELECTED = {"aria-1": 0.0006, "forte-1": 0.0004}
SHOTS = 20000
BACKENDS = ["ideal", "aria-1", "forte-1"]
ATOMS_OVERLAP = [2, 3]
LEAKAGE_FREE = "1"  # weight-1 (odd) register -- see module docstring
OUT_PATH = os.path.join(os.path.dirname(__file__), "task51_overlap_fragment_results.json")


def build_overlap_ansatz(theta):
    qc = QuantumCircuit(2)
    qc.x(0)
    qc.append(XXPlusYYGate(theta[0], np.pi / 2), [0, 1])
    return qc


def _residual(theta, target):
    sv = np.asarray(Statevector.from_instruction(build_overlap_ansatz(theta)))
    idx = int(np.argmax(np.abs(target)))
    phase = sv[idx] / target[idx] if abs(target[idx]) > 1e-9 else 1.0
    diff = sv / phase - target
    return np.concatenate([diff.real, diff.imag])


def fit_overlap_angle(target, n_attempts=10, tol=1e-12):
    best_err, best_x = None, None
    for attempt in range(n_attempts):
        rng = np.random.default_rng(attempt)
        x0 = rng.uniform(-np.pi, np.pi, 1)
        res = least_squares(_residual, x0, args=(target,), method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
        sv = np.asarray(Statevector.from_instruction(build_overlap_ansatz(res.x)))
        idx = int(np.argmax(np.abs(target)))
        phase = sv[idx] / target[idx] if abs(target[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv / phase - target)))
        if best_err is None or err < best_err:
            best_err, best_x = err, res.x
        if best_err < tol:
            break
    return best_x, best_err


def ancilla_cnots_abstract_n(n_register):
    qc = QuantumCircuit(n_register + 1)
    for q in range(n_register):
        qc.cx(q, n_register)
    return qc


def verify_ancilla_scheme_n(build_ansatz_fn, angles, n_register, leakage_free_value):
    base = build_ansatz_fn(angles)
    qc5 = QuantumCircuit(n_register + 1)
    qc5.compose(base, qubits=list(range(n_register)), inplace=True)
    qc5.compose(ancilla_cnots_abstract_n(n_register), inplace=True)

    sv_full = Statevector.from_instruction(qc5)
    sv_reg = Statevector.from_instruction(base)
    rho_direct = np.outer(np.asarray(sv_reg), np.asarray(sv_reg).conj())
    rho_marginal = np.asarray(partial_trace(sv_full, [n_register]))
    marginal_err = float(np.max(np.abs(rho_direct - rho_marginal)))

    probs = sv_full.probabilities_dict()
    leaked_value = "1" if leakage_free_value == "0" else "0"
    p_leaked = sum(v for k, v in probs.items() if k[0] == leaked_value)
    return marginal_err, p_leaked


def analytic_A_and_B_conditioned_n(build_ansatz_fn, angles, gate_name, p_zz, p_gpi, p_gpi2,
                                    delta_zz, delta_gpi2, labels, n_register, leakage_free_value="0"):
    qc = to_native(build_ansatz_fn(angles), gate_name)
    n = qc.num_qubits
    assert n == n_register
    dm_A = DensityMatrix.from_label("0" * n)
    dm_B = DensityMatrix.from_label("0" * n)
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
            p_here, n_here = p_gpi, 1
        elif op.name == "gpi2":
            phi = float(op.params[0]) % 1.0
            U = Operator(biased_gpi2_matrix(phi, delta_gpi2))
            p_here, n_here = p_gpi2, 1
        else:
            dm_A = dm_A.evolve(Operator(op.to_matrix()), qargs=qargs)
            dm_B = dm_B.evolve(Operator(op.to_matrix()), qargs=qargs)
            continue
        dm_A = dm_A.evolve(U, qargs=qargs)
        dm_B = dm_B.evolve(U, qargs=qargs)
        dm_A = apply_pauli_mixture(dm_A, qargs, depolarizing_weights(p_here, n_here))
        dm_B = apply_pauli_mixture(dm_B, qargs, depolarizing_weights(p_here, n_here))
        dm_B = apply_pauli_mixture(dm_B, qargs, pec_inverse_weights(p_here, n_here))

    P_even = even_weight_projector(n_register)
    P_valid = P_even if leakage_free_value == "0" else (np.eye(P_even.shape[0]) - P_even)

    def condition(dm):
        rho = dm.data
        rho_cond = P_valid @ rho @ P_valid
        norm = float(np.real(np.trace(rho_cond)))
        if norm < 1e-12:
            return dm, 0.0
        return DensityMatrix(rho_cond / norm), norm

    dm_A_cond, retained_A = condition(dm_A)
    dm_B_cond, retained_B = condition(dm_B)
    A = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_A_cond.data))) for l in labels}
    B = {l: float(np.real(np.trace(np.asarray(Pauli(l).to_matrix()) @ dm_B_cond.data))) for l in labels}
    return A, B, retained_A, retained_B


def apply_correction_n(build_ansatz_fn, raw_by_slot_label, p_gpi2, kept, non_id_labels, fixed_solutions,
                        n_register, leakage_free_value="0"):
    corrected = {name: {} for name in kept}
    for name in kept:
        A, B, retA, retB = analytic_A_and_B_conditioned_n(
            build_ansatz_fn, fixed_solutions[name], GATE_NAME, ZZ_ASSUMED, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0,
            non_id_labels, n_register, leakage_free_value)
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            corrected[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    return corrected


def setup_and_verify():
    p = setup_fragment(ATOMS_OVERLAP, nelec=2, d=1.0, K=K, strict=True)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    diag, plus, kept = kept_slots_for_K(K)
    print(f"  kept slots: {kept}  non_id_labels: {non_id_labels}")
    print(f"  register physical subspace: weight-1 (odd) -> leakage-free ancilla value = {LEAKAGE_FREE!r}")

    fixed_solutions = {}
    worst = 0.0
    for name in kept:
        theta, err = fit_overlap_angle(p["targets"][name])
        fixed_solutions[name] = theta
        worst = max(worst, err)
        print(f"    {name}: theta={theta[0]:+.6f}  max_abs_error={err:.3e}  {'OK' if err < 1e-10 else 'FAIL'}")
    assert worst < 1e-10, "angle fit did not converge to machine precision -- STOP"

    worst_marginal, worst_p_leaked = 0.0, 0.0
    for name in kept:
        marginal_err, p_leaked = verify_ancilla_scheme_n(build_overlap_ansatz, fixed_solutions[name], 2, LEAKAGE_FREE)
        worst_marginal = max(worst_marginal, marginal_err)
        worst_p_leaked = max(worst_p_leaked, p_leaked)
    print(f"    worst marginal-state error: {worst_marginal:.3e}  {'PASS' if worst_marginal < 1e-9 else 'FAIL -- STOP'}")
    print(f"    worst ideal p(leaked): {worst_p_leaked:.3e}  {'PASS' if worst_p_leaked < 1e-9 else 'FAIL -- STOP'}")
    assert worst_marginal < 1e-9 and worst_p_leaked < 1e-9, "ancilla-parity verification failed -- STOP"

    A0, B0, retA0, retB0 = analytic_A_and_B_conditioned_n(
        build_overlap_ansatz, fixed_solutions["u_0"], GATE_NAME, 0.0, 0.0, 0.0, 0.0, 0.0, non_id_labels, 2, LEAKAGE_FREE)
    max_diff0 = max(abs(A0[l] - B0[l]) for l in non_id_labels)
    print(f"    max diff A vs B at true theta=0: {max_diff0:.3e}  retained_A={retA0:.6f}  "
          f"{'PASS' if max_diff0 < 1e-9 and abs(retA0 - 1.0) < 1e-9 else 'FAIL -- STOP'}")
    assert max_diff0 < 1e-9 and abs(retA0 - 1.0) < 1e-9, "conditioning sanity check failed -- STOP"
    return p, non_id_labels, kept, fixed_solutions


def submit(resume=True):
    p, non_id_labels, kept, fixed_solutions = setup_and_verify()
    ancilla_native = to_native(ancilla_cnots_abstract_n(2), GATE_NAME)
    groups = effrag_mod.group_labels_qubit_wise(non_id_labels)
    print(f"  QWC groups for overlap: {groups}")

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    if resume and os.path.exists(OUT_PATH):
        with open(OUT_PATH) as f:
            state = json.load(f)
        print(f"  resuming: {len(state['done'])} already done")
    else:
        state = {"done": {}}

    t0 = time.time()
    for backend_name in BACKENDS:
        for name in kept:
            key = f"{backend_name}|{name}"
            if key in state["done"]:
                print(f"    skip (done): {key}")
                continue
            register_native = to_native(build_overlap_ansatz(fixed_solutions[name]), GATE_NAME)
            circuits = []
            for group in groups:
                combined = effrag_mod.combined_basis_label(group)
                basis_abstract = QuantumCircuit(2)
                for i, ch in enumerate(combined):
                    qubit = 2 - 1 - i
                    if ch == "X":
                        basis_abstract.h(qubit)
                    elif ch == "Y":
                        basis_abstract.sdg(qubit)
                        basis_abstract.h(qubit)
                basis_qc = to_native(basis_abstract, GATE_NAME)
                full3 = QuantumCircuit(3)
                full3.compose(register_native, qubits=[0, 1], inplace=True)
                full3.compose(ancilla_native, inplace=True)
                full3 = full3.compose(basis_qc, qubits=[0, 1])
                full3.measure_all()
                circuits.append(full3)
            job = submit_job(circuits, backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"groups": groups, "counts": counts}
            print(f"    done: {key} ({len(circuits)} circuits, {SHOTS} shots each, {time.time()-t0:.1f}s elapsed)")
            with open(OUT_PATH, "w") as f:
                json.dump(state, f, indent=2)
    print(f"\n  ALL SUBMISSIONS DONE in {time.time()-t0:.1f}s. Saved -> {OUT_PATH}")
    return p, non_id_labels, kept, fixed_solutions


def analyze(p=None, non_id_labels=None, kept=None, fixed_solutions=None):
    if p is None:
        p, non_id_labels, kept, fixed_solutions = setup_and_verify()
    with open(OUT_PATH) as f:
        state = json.load(f)

    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)
    weight_unit = {name: {l: 1.0 for l in non_id_labels} for name in kept}

    def energy_and_err(raw):
        alpha_mats = combine_matrices(raw, p["alpha_labels"], p["identity_label"], K)
        E, errs = energy_from_alpha_matrices(alpha_mats, p["terms"], p["lambdas"], p["enuc"], p["signs"], K,
                                              exact_energy=p["exact_energy"], noiseless_energy=p["noiseless_energy"])
        return E, errs["err_vs_exact_kcal"]

    print(f"\n  overlap exact energy: {p['exact_energy']:.6f} Ha")
    results = {}
    for backend_name in BACKENDS:
        blended = {name: {} for name in kept}
        accept_fracs = []
        for name in kept:
            entry = state["done"][f"{backend_name}|{name}"]
            for group, counts in zip(entry["groups"], entry["counts"]):
                total = sum(counts.values())
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == LEAKAGE_FREE}
                accept_fracs.append(sum(filtered.values()) / total)
                for l in group:
                    blended[name][l] = expectation_from_counts(filtered, l)

        if backend_name == "ideal":
            data_for_fit = blended
        else:
            data_for_fit = apply_correction_n(build_overlap_ansatz, blended, GPI2_SELECTED[backend_name],
                                               kept, non_id_labels, fixed_solutions, 2, LEAKAGE_FREE)

        rng = np.random.default_rng(39)
        U_hat, cost, chi2dof = fit_joint_frame(np.eye(K), P_S, K, kept, non_id_labels, data_for_fit,
                                                  weight_unit, rng, n_restarts=4)
        full = build_full_from_frame(U_hat, P_S, K, non_id_labels, kept)
        E, err = energy_and_err(full)
        print(f"    {backend_name}: mean_accept={np.mean(accept_fracs):.4f}  E={E:.6f} Ha  "
              f"err_vs_exact={err:+.4f} kcal/mol  chi2/dof={chi2dof:.5f}")
        results[backend_name] = {"E": E, "err_kcal": err, "mean_accept": float(np.mean(accept_fracs))}
    return results


if __name__ == "__main__":
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("  WARNING: PYTHONHASHSEED != 0 -- known nondeterminism risk (Task 36): the joint-frame "
              "nonconvex fit can land in a different, worse local minimum. Rerun with PYTHONHASHSEED=0.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyze", action="store_true", help="analyze already-submitted data only, no network calls")
    args = ap.parse_args()

    if args.analyze:
        analyze()
    else:
        p, non_id_labels, kept, fixed_solutions = submit()
        analyze(p, non_id_labels, kept, fixed_solutions)
