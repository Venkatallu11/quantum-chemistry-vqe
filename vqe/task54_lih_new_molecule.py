#!/usr/bin/env python3
"""
task54_lih_new_molecule.py -- iteration 54. Testing the validated forging +
ancilla-parity + conditioned-PEC/GPi2 + joint-Schmidt-frame pipeline (the
one that reached 0.0105-0.0192 kcal/mol on H4's own fragment, and
0.0025-0.0094 kcal/mol on the real covalent-bonding tailored result) on a
GENUINELY DIFFERENT molecule -- not another H-chain fragment. LiH: a real
heteronuclear molecule (Li has actual p-orbitals, unlike a bare H chain),
4 electrons, 6 spatial orbitals -> 12 qubits, same qubit-count class as
the full H6 problem (already shown tractable on this machine).

CROSS-CHECK AGAINST AN INDEPENDENT PRIOR RESULT: this repo's own
molecules_real.py / molecules_real_results.json already has a trusted LiH
reference computed independently (a totally different code path -- no
Schmidt decomposition, no forging, just direct dense FCI): exact_energy_ha
= -7.882324, hf_energy_ha = -7.861865. Rebuilding the SAME Hamiltonian here
via ef_fragment's own real-gauge + Schmidt-decompose pipeline reproduces
both numbers exactly (diff 0.0002 kcal/mol, just float rounding) --
confirms this new molecule's setup is correct before anything else.

REAL, HONEST FINDING (the interesting part): LiH's Schmidt rank is EXACTLY
K=15 -- the full theoretical maximum for its 6-qubit, 2-alpha-electron
weight-2 register (C(6,2)=15). Unlike H4 (also full-rank at K=6, but out of
a much SMALLER 4-qubit register), LiH does not compress the way an H-chain
fragment does; its natural register is just bigger. A full-rank K=15 real
campaign would need ~5,850 alpha-only circuits (26 qubit-wise groups x 225
slots, beta-reuse holds, sign_residual=6.3e-12) -- large but not
unreasonable, NOT attempted here.

INSTEAD: a K=6 TRUNCATION was chosen deliberately -- classically verified
(no quantum measurement involved) to already give 0.0176 kcal/mol vs the
full exact energy, comparable to H4's OWN real hardware noise floor
(0.0105-0.0192 kcal/mol). This makes K=6 the natural, honest choice for a
real test of whether the METHOD generalizes to a new molecule: same K,
same n_slots=36 (6 diagonal + 30 phase-pair) as H4's own real campaign, so
the real submission is the SAME SCALE (36 job submissions per backend, one
per named slot, each containing this Hamiltonian's own 26 qubit-wise-
commuting circuits) as an already-proven-tractable run -- not the ~10x
scale-up the full-H6/bigger-overlap paths would need.

ANSATZ: unlike H4/overlap's hand-derived weight-preserving circuits, this
uses native_stateprep.build_real_state_prep_circuit(vec) -- ALREADY
GENERIC and ALREADY VALIDATED (native_stateprep.py's own docstring:
"Verified correct via Statevector to <1e-15 across 60 random real
vectors"), so no new circuit derivation is needed; each of the 6 required
Schmidt-vector/phase-state targets gets an EXACT (to machine precision by
construction, not by numerical fitting) real-amplitude state prep. This is
NOT gate-count-optimized (a hand-derived weight-2-preserving ansatz, like
H4's, would likely need far fewer 2-qubit gates) -- correctness on a new
molecule was the priority for this first test, not circuit efficiency.

ANCILLA-PARITY CONVENTION: weight-2 (EVEN), same convention as H4/fragment
B -- LEAKAGE_FREE="0" (unlike the overlap fragment's odd-weight case).

Run:
    PYTHONHASHSEED=0 python vqe/task54_lih_new_molecule.py            # full build+verify+submit+analyze
    PYTHONHASHSEED=0 python vqe/task54_lih_new_molecule.py --analyze  # analyze already-submitted data only
    python vqe/task54_lih_new_molecule.py --classical-only            # setup+verify only, no network calls
"""
import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import chem
import ef_fragment as effrag
from entanglement_forging_h4 import precompute_exact_matrices, ef_energy_from_matrices
from qforge.forging import (
    group_labels_qubit_wise, combined_basis_label, combine_matrices, energy_from_alpha_matrices,
)
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task36_joint_schmidt_frame import fit_joint_frame, build_full_from_frame
from native_stateprep import build_real_state_prep_circuit, to_native
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, expectation_from_counts
from task37c_extended_forward_model import biased_zz_matrix, biased_gpi_matrix, biased_gpi2_matrix
from task37b_h4_noise_model import GPI_REAL_MEAN
from loop_pec import depolarizing_weights, pec_inverse_weights, apply_pauli_mixture
from task39e_conditioned_correction import even_weight_projector

from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper
from molecules_real import _constrain_particle_number
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator, partial_trace, DensityMatrix, Pauli

HARTREE_TO_KCAL_MOL = 627.5094740631
GEOM = [("Li", (0, 0, 0)), ("H", (0, 0, 1.6))]
NELEC = 4
K = 6
N_REGISTER = 6
LEAKAGE_FREE = "0"  # weight-2, even -- same convention as H4/fragment B
GATE_NAME = "zz"
ZZ_ASSUMED = 0.014593
GPI2_SELECTED = {"aria-1": 0.0006, "forte-1": 0.0004}
SHOTS = 20000
BACKENDS = ["ideal", "aria-1", "forte-1"]
TRUSTED_EXACT_HA = -7.882324  # molecules_real_results.json, independent code path
TRUSTED_HF_HA = -7.861865
OUT_PATH = os.path.join(os.path.dirname(__file__), "task54_lih_new_molecule_results.json")


def build_lih_qop():
    S, T, V, eri, enuc = chem.integrals(GEOM)
    ehf, C, Hc = chem.rhf(S, T, V, eri, enuc, nelec=NELEC)
    h1 = C.T @ Hc @ C
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)
    ee = ElectronicEnergy.from_raw_integrals(h1, h2)
    qop_bare = JordanWignerMapper().map(ee.second_q_op())
    qop_pen = _constrain_particle_number(qop_bare, NELEC)
    return qop_bare, qop_pen, enuc, ehf


def setup_lih_fragment(K):
    qop_bare, qop_pen, enuc, ehf = build_lih_qop()
    n_qubits = qop_bare.num_qubits
    assert n_qubits == 2 * N_REGISTER

    hf_err_kcal = abs(ehf - TRUSTED_HF_HA) * HARTREE_TO_KCAL_MOL
    assert hf_err_kcal < 0.01, f"HF energy does not match trusted repo reference: diff={hf_err_kcal:.4f} kcal/mol"

    e_elec, psi = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    exact_err_kcal = abs(exact_energy - TRUSTED_EXACT_HA) * HARTREE_TO_KCAL_MOL
    assert exact_err_kcal < 0.01, f"exact energy does not match trusted repo reference: diff={exact_err_kcal:.4f} kcal/mol"

    psi_real, real_gauge_residual = effrag.real_gauge(psi)
    lambdas, u_vecs, v_vecs = effrag.schmidt_decompose_real(psi_real, n_qubits)
    max_tail_at_K = float(np.max(np.abs(lambdas[K:])))

    terms = effrag.decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    assert set(alpha_labels) == set(beta_labels)
    non_even_labels = [l for l in alpha_labels if not effrag.label_is_real(l)]
    assert not non_even_labels, f"non-even-Y labels found: {non_even_labels}"
    identity_label = "I" * N_REGISTER

    u_top, v_top = u_vecs[:K], v_vecs[:K]
    signs = np.array([1.0 if np.dot(v_top[n], u_top[n]) >= 0 else -1.0 for n in range(K)])
    sign_residual = max(float(np.max(np.abs(v_top[n] - signs[n] * u_top[n]))) for n in range(K))
    assert sign_residual < 1e-8, f"beta-reuse shortcut failed: residual={sign_residual:.3e}"

    alpha_cache, beta_cache = precompute_exact_matrices(terms, u_top, v_top)
    noiseless_energy = ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)
    truncation_err_kcal = abs(noiseless_energy - exact_energy) * HARTREE_TO_KCAL_MOL

    diag = [f"u_{n}" for n in range(K)]
    plus_minus = []
    for n in range(K):
        for m in range(K):
            if n < m:
                plus_minus += [f"(u{n}+u{m})", f"(u{n}-u{m})"]
    targets = {f"u_{n}": u_top[n] for n in range(K)}
    for n in range(K):
        for m in range(K):
            if n < m:
                targets[f"(u{n}+u{m})"] = (u_top[n] + u_top[m]) / np.sqrt(2)
                targets[f"(u{n}-u{m})"] = (u_top[n] - u_top[m]) / np.sqrt(2)

    return {
        "n_qubits": n_qubits, "terms": terms, "alpha_labels": alpha_labels,
        "identity_label": identity_label, "signs": signs, "lambdas": lambdas, "u_vecs": u_top,
        "enuc": enuc, "exact_energy": exact_energy, "noiseless_energy": noiseless_energy,
        "hf_energy": ehf, "max_schmidt_tail_at_K": max_tail_at_K,
        "truncation_err_vs_exact_kcal": truncation_err_kcal,
        "real_gauge_residual": real_gauge_residual, "targets": targets, "K": K,
        "sign_residual": sign_residual,
    }


def build_lih_ansatz(vec):
    return build_real_state_prep_circuit(np.asarray(vec, dtype=float))


def ancilla_cnots_abstract_n(n_register):
    qc = QuantumCircuit(n_register + 1)
    for q in range(n_register):
        qc.cx(q, n_register)
    return qc


def verify_ancilla_scheme_n(build_ansatz_fn, vec, n_register, leakage_free_value):
    base = build_ansatz_fn(vec)
    qc = QuantumCircuit(n_register + 1)
    qc.compose(base, qubits=list(range(n_register)), inplace=True)
    qc.compose(ancilla_cnots_abstract_n(n_register), inplace=True)

    sv_full = Statevector.from_instruction(qc)
    sv_reg = Statevector.from_instruction(base)
    rho_direct = np.outer(np.asarray(sv_reg), np.asarray(sv_reg).conj())
    rho_marginal = np.asarray(partial_trace(sv_full, [n_register]))
    marginal_err = float(np.max(np.abs(rho_direct - rho_marginal)))

    probs = sv_full.probabilities_dict()
    leaked_value = "1" if leakage_free_value == "0" else "0"
    p_leaked = sum(v for k, v in probs.items() if k[0] == leaked_value)
    return marginal_err, p_leaked


def analytic_A_and_B_conditioned_n(build_ansatz_fn, vec, gate_name, p_zz, p_gpi, p_gpi2,
                                    delta_zz, delta_gpi2, labels, n_register, leakage_free_value="0"):
    qc = to_native(build_ansatz_fn(vec), gate_name)
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


def apply_correction_n(build_ansatz_fn, raw_by_slot_label, p_gpi2, kept, non_id_labels, targets,
                        n_register, leakage_free_value="0"):
    corrected = {name: {} for name in kept}
    for name in kept:
        A, B, retA, retB = analytic_A_and_B_conditioned_n(
            build_ansatz_fn, targets[name], GATE_NAME, ZZ_ASSUMED, GPI_REAL_MEAN, p_gpi2, 0.0, 0.0,
            non_id_labels, n_register, leakage_free_value)
        for l, m_raw in raw_by_slot_label[name].items():
            ratio = B[l] / A[l] if abs(A[l]) > 1e-6 else 1.0
            corrected[name][l] = max(-1.0, min(1.0, m_raw * ratio))
    return corrected


def setup_and_verify():
    p = setup_lih_fragment(K)
    non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
    diag, plus, kept = kept_slots_for_K(K)
    print(f"  LiH: n_qubits={p['n_qubits']}, K={K} (true rank=15, deliberately truncated)")
    print(f"  HF energy matches trusted reference (diff<0.01 kcal/mol): PASS")
    print(f"  exact energy matches trusted reference (diff<0.01 kcal/mol): PASS")
    print(f"  exact_energy={p['exact_energy']:.6f} Ha, noiseless K={K} EF energy={p['noiseless_energy']:.6f} Ha")
    print(f"  classical K={K} truncation error vs full exact: {p['truncation_err_vs_exact_kcal']:.4f} kcal/mol")
    print(f"  beta-reuse sign_residual={p['sign_residual']:.3e} (OK)")
    print(f"  {len(non_id_labels)} non-identity alpha labels, {len(kept)} kept slots")

    groups = group_labels_qubit_wise(non_id_labels)
    print(f"  {len(groups)} qubit-wise-commuting groups")

    worst_prep_err = 0.0
    for name in kept:
        sv = np.asarray(Statevector.from_instruction(build_lih_ansatz(p["targets"][name])))
        target = p["targets"][name]
        idx = int(np.argmax(np.abs(target)))
        phase = sv[idx] / target[idx] if abs(target[idx]) > 1e-9 else 1.0
        err = float(np.max(np.abs(sv / phase - target)))
        worst_prep_err = max(worst_prep_err, err)
    print(f"  worst ansatz statevector error (exact by construction): {worst_prep_err:.3e}  "
          f"{'PASS' if worst_prep_err < 1e-9 else 'FAIL -- STOP'}")
    assert worst_prep_err < 1e-9, "ansatz state prep is not exact -- STOP"

    worst_marginal, worst_p_leaked = 0.0, 0.0
    for name in kept:
        marginal_err, p_leaked = verify_ancilla_scheme_n(build_lih_ansatz, p["targets"][name], N_REGISTER, LEAKAGE_FREE)
        worst_marginal = max(worst_marginal, marginal_err)
        worst_p_leaked = max(worst_p_leaked, p_leaked)
    print(f"  worst marginal-state error: {worst_marginal:.3e}  {'PASS' if worst_marginal < 1e-9 else 'FAIL -- STOP'}")
    print(f"  worst ideal p(leaked): {worst_p_leaked:.3e}  {'PASS' if worst_p_leaked < 1e-9 else 'FAIL -- STOP'}")
    assert worst_marginal < 1e-9 and worst_p_leaked < 1e-9, "ancilla-parity verification failed -- STOP"

    A0, B0, retA0, retB0 = analytic_A_and_B_conditioned_n(
        build_lih_ansatz, p["targets"]["u_0"], GATE_NAME, 0.0, 0.0, 0.0, 0.0, 0.0, non_id_labels, N_REGISTER, LEAKAGE_FREE)
    max_diff0 = max(abs(A0[l] - B0[l]) for l in non_id_labels)
    print(f"  max diff A vs B at zero injected noise: {max_diff0:.3e}  retained_A={retA0:.6f}  "
          f"{'PASS' if max_diff0 < 1e-9 and abs(retA0 - 1.0) < 1e-9 else 'FAIL -- STOP'}")
    assert max_diff0 < 1e-9 and abs(retA0 - 1.0) < 1e-9, "conditioning sanity check failed -- STOP"

    n2q_counts = [to_native(build_lih_ansatz(p["targets"][name]), GATE_NAME).count_ops().get(GATE_NAME, 0) for name in kept]
    print(f"  native 2-qubit gate count per slot: min={min(n2q_counts)} max={max(n2q_counts)} "
          f"(unoptimized generic real-state-prep tree, not hand-tuned)")

    return p, non_id_labels, kept, groups


def submit(resume=True):
    p, non_id_labels, kept, groups = setup_and_verify()

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")
    ancilla_native = to_native(ancilla_cnots_abstract_n(N_REGISTER), GATE_NAME)

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
            register_native = to_native(build_lih_ansatz(p["targets"][name]), GATE_NAME)
            circuits = []
            for group in groups:
                combined = combined_basis_label(group)
                basis_abstract = QuantumCircuit(N_REGISTER)
                for i, ch in enumerate(combined):
                    qubit = N_REGISTER - 1 - i
                    if ch == "X":
                        basis_abstract.h(qubit)
                    elif ch == "Y":
                        basis_abstract.sdg(qubit)
                        basis_abstract.h(qubit)
                basis_qc = to_native(basis_abstract, GATE_NAME)
                full = QuantumCircuit(N_REGISTER + 1)
                full.compose(register_native, qubits=list(range(N_REGISTER)), inplace=True)
                full.compose(ancilla_native, inplace=True)
                full = full.compose(basis_qc, qubits=list(range(N_REGISTER)))
                full.measure_all()
                circuits.append(full)
            job = submit_job(circuits, backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            state["done"][key] = {"groups": groups, "counts": counts}
            print(f"    done: {key} ({len(circuits)} circuits, {SHOTS} shots each, {time.time()-t0:.1f}s elapsed)")
            with open(OUT_PATH, "w") as f:
                json.dump(state, f, indent=2)
    print(f"\n  ALL SUBMISSIONS DONE in {time.time()-t0:.1f}s. Saved -> {OUT_PATH}")
    return p, non_id_labels, kept, groups


def analyze(p=None, non_id_labels=None, kept=None, groups=None):
    if p is None:
        p, non_id_labels, kept, groups = setup_and_verify()
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

    print(f"\n  LiH exact energy: {p['exact_energy']:.6f} Ha  (K={K} classical truncation floor: "
          f"{p['truncation_err_vs_exact_kcal']:.4f} kcal/mol)")
    results = {}
    for backend_name in BACKENDS:
        blended = {name: {} for name in kept}
        accept_fracs = []
        for name in kept:
            entry = state["done"][f"{backend_name}|{name}"]
            for group, counts in zip(entry["groups"], entry["counts"]):
                total = sum(counts.values())
                filtered = {bs[1:]: c for bs, c in counts.items() if bs[0] == LEAKAGE_FREE}
                accept_fracs.append(sum(filtered.values()) / total if total else 0.0)
                for l in group:
                    blended[name][l] = expectation_from_counts(filtered, l)

        if backend_name == "ideal":
            data_for_fit = blended
        else:
            data_for_fit = apply_correction_n(build_lih_ansatz, blended, GPI2_SELECTED[backend_name],
                                               kept, non_id_labels, p["targets"], N_REGISTER, LEAKAGE_FREE)

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
    ap.add_argument("--classical-only", action="store_true", help="setup+verify only, no network calls")
    args = ap.parse_args()

    if args.classical_only:
        setup_and_verify()
    elif args.analyze:
        analyze()
    else:
        p, non_id_labels, kept, groups = submit()
        analyze(p, non_id_labels, kept, groups)
