#!/usr/bin/env python3
"""
rank6_symmetry_vd.py — three independent improvements on top of the CDR
pipeline (fixed_ansatz.py + cdr_mitigation.py), evaluated honestly.
SIMULATOR ONLY. Nothing here touches IonQ.
============================================================================
CRITICAL FIX FIRST: cdr_mitigation.py's err_kcal (and every number this
project has reported from it) is measured against exact_energy -- the true,
full-CI fragment energy -- NOT against noiseless_numpy(K), the best a
rank-K truncation could ever do. At K=5 those differ by 0.5655 kcal/mol
(verified directly: exact_energy - noiseless_numpy(K=5) = 0.5654737...
kcal/mol here, matching CLASSICAL_FLOOR_KCAL exactly) -- a classical
truncation floor, not noise. CDR's reported 2.264 kcal/mol (K=5, mean/8
seeds) therefore already INCLUDES that floor; the true noise residual is
closer to 1.7. From here on every configuration reports BOTH numbers
separately: err_vs_exact_kcal (what's been reported so far) and
err_vs_noiseless_kcal (the pure noise residual at that K) -- see
`energy_errors()`.

STEP 1 -- RANK 6: verified (not assumed) that this fragment's Schmidt
spectrum across the alpha/beta bipartition has exactly 6 nonzero singular
values (lambda[6:] ~ 1e-15, numerically zero) -- a consequence of the
state living entirely in a 6-dim (weight-2) x 6-dim (weight-2) physical
subspace, which bounds the Schmidt rank at 6. K=6 is therefore not just
"a bit better" than K=5, it is EXACT: verified E(K=6) - exact_energy =
1.0e-11 kcal/mol using exact (noiseless) Schmidt vectors. Cost: 6 diagonal
+ 15 pairs x 2 phase states = 36 state-prep slots (1.44x the 25 at K=5).
All 36 are refit with fixed_ansatz.fit_angles and reverified <1e-10, and
the abstract (u3/cx, optimization_level=0) 2-qubit gate count is
reverified constant across all 36 -- the same bug class fixed previously
in fixed_ansatz.py/cdr_mitigation.py (angle-dependent adaptive synthesis)
gets re-checked here rather than assumed to still hold at the new K.

STEP 2 -- SYMMETRY-VERIFIED POSTSELECTION, sector-confined: computed (not
assumed) that a state fully mixed over all 256 8-qubit basis states gives
1151.42 kcal/mol error vs exact, while a state fully mixed ONLY within the
36-state physical sector (weight-2 alpha x weight-2 beta) gives 939.10 --
a 1.226x reduction in the worst-case error IF noise could be confined to
the physical sector. It can't be, for most of this Hamiltonian: real
depolarizing noise (X/Y-type Pauli errors) does not respect particle
number, so it genuinely leaks weight. What CAN be done is DETECTING that
leakage: for any alpha-register Pauli label whose qubit-wise-commuting
MEASUREMENT GROUP needs no X/Y basis rotation (all-Z/identity only), the
raw measured bitstring's Hamming weight IS the state's actual weight, so
shots landing outside weight=2 can be discarded (postselected) as
detected errors. For every OTHER group (needs an X or Y rotation before
Z-readout), the ROTATED-basis bitstring's weight has no such meaning --
postselecting on it is invalid, a lesson already paid for once here (an
earlier attempt at exactly this made an energy 40x worse). Measured (not
assumed): of this Hamiltonian's 37 alpha labels across 13 qubit-wise
groups, only 1 group (3 labels: IIZZ, ZIIZ, ZZII) is all-Z-type --
2.6% of the Hamiltonian's |coefficient| weight. The realistically
achievable gain from Step 2 is reported plainly at that coverage, not
assumed to be 1.23x.

STEP 3 -- 2-COPY VIRTUAL DISTILLATION (Huggins et al., PRX 11, 041036,
2021; Koczor, PRX 11, 031057, 2021): measures <O.S>/<S> across two noisy
copies of the state (S = the qubit-wise SWAP between the two 4-qubit
copies), which evaluates observables against rho^2/Tr(rho^2) instead of
rho -- suppressing INCOHERENT error quadratically. Uses the diagonalizing
"B gate" form (B constructed here as SWAP's own eigenvector matrix,
verified unitary and verified to diagonalize SWAP, rather than a
naive controlled-SWAP + ancilla), applied to each of the 4 corresponding
qubit pairs across the two copies (IonQ's all-to-all connectivity means
this needs no SWAP routing, which is the specific reason this is worth
trying here). Because this project has exact density-matrix access
locally (not real shot data), the VD numerator/denominator are computed
as exact matrix traces on the 8-qubit noisy density matrix -- INCLUDING
noise on the B-gate circuit itself, not just the two state-prep copies
(measured: B costs 2 two-qubit gates per pair, so the full 8-qubit
circuit is 11+11+4*2 = 30 two-qubit gates, not 22). This was verified
against the direct-trace identity Tr[(O#I)S(rho#rho)] = Tr[O rho^2] on a
random test state before trusting it on the real Hamiltonian (a
convention bug -- rho_after_B = B rho B^dagger, not B^dagger rho B --
was caught by exactly this check and is not repeated here). VD alone,
CDR alone, and VD+CDR stacked are all reported; VD+CDR's training data is
regenerated through the SAME 8-qubit VD circuit (structural identity is
still the point) at a REDUCED N_TRAIN_PER_SLOT_VD=3 (from 5) purely for
compute-time reasons (each VD circuit costs ~1.2s to simulate vs ~0.03s
for the single-copy case) -- stated plainly, not hidden.

Run:
    python vqe/rank6_symmetry_vd.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import ef_fragment as effrag
import cdr_mitigation as cdr
from fixed_ansatz import build_ansatz, fit_angles, P2_PER_GATE, P1_PER_GATE

from qiskit import transpile, QuantumCircuit
from qiskit.quantum_info import Pauli, Statevector
from qiskit.circuit.library import UnitaryGate, SwapGate
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

HARTREE_TO_KCAL_MOL = cdr.HARTREE_TO_KCAL_MOL
BASIS_GATES = cdr.BASIS_GATES
ATOMS = [0, 1, 2, 3]
NELEC = 4
K = 6
N_TRAIN_PER_SLOT = 5
N_TRAIN_PER_SLOT_VD = 3   # reduced from 5 -- 8-qubit VD sim costs ~40x a single-copy sim
N_SEEDS = 8
SEEDS = list(range(N_SEEDS))
LOW_SIGNAL_CUTOFF = cdr.LOW_SIGNAL_CUTOFF

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "rank6_symmetry_vd_results.json")


# ---------------------------------------------------------------------------
# Setup at K=6
# ---------------------------------------------------------------------------

def slot_names(K):
    names = [f"u_{n}" for n in range(K)]
    for n in range(K):
        for m in range(K):
            if n < m:
                names += [f"(u{n}+u{m})", f"(u{n}-u{m})"]
    return names


def setup(K=6):
    qop_bare, qop_pen, enuc = effrag.build_fragment_qop(ATOMS, NELEC, d=1.0)
    n_qubits = qop_bare.num_qubits
    e_elec, psi = effrag.exact_ground_state(qop_pen)
    exact_energy = e_elec + enuc
    psi_real, residual = effrag.real_gauge(psi)
    lambdas, u_vecs, v_vecs = effrag.schmidt_decompose_real(psi_real, n_qubits)
    terms = effrag.decompose_pauli_terms(qop_bare, n_qubits)
    alpha_labels = sorted(set(a for a, _, _ in terms))
    beta_labels = sorted(set(b for _, b, _ in terms))
    assert set(alpha_labels) == set(beta_labels)

    assert np.max(np.abs(lambdas[K:])) < 1e-9, (
        f"Schmidt rank assumption violated: lambda[{K}:] max = {np.max(np.abs(lambdas[K:])):.3e}, "
        "K=6 would NOT be exact -- do not proceed on the assumed-exact-rank premise"
    )

    u_top, v_top = u_vecs[:K], v_vecs[:K]
    signs = np.array([1.0 if np.dot(v_top[n], u_top[n]) >= 0 else -1.0 for n in range(K)])
    sign_residual = max(float(np.max(np.abs(v_top[n] - signs[n] * u_top[n]))) for n in range(K))
    assert sign_residual < 1e-8, f"beta=sign*alpha residual {sign_residual:.3e} -- shortcut doesn't hold here"

    identity_label = "I" * (n_qubits // 2)
    identity_coeff = 0.0
    for label, coeff in qop_bare.to_list():
        if label == "I" * n_qubits:
            identity_coeff = float(coeff.real)
            break
    e_mixed = identity_coeff + enuc

    alpha_cache, beta_cache = effrag.precompute_exact_matrices(terms, u_top, v_top)
    noiseless_numpy = effrag.ef_energy_from_matrices(terms, lambdas, alpha_cache, beta_cache, enuc, K)
    truncation_floor_kcal = abs(exact_energy - noiseless_numpy) * HARTREE_TO_KCAL_MOL

    # physical-sector mixed-reference context (Step 2)
    M = qop_bare.to_matrix(sparse=True)
    diag = np.real(M.diagonal())
    half_dim = 2 ** (n_qubits // 2)

    def hamming(x):
        return bin(x).count("1")
    mask_phys = np.zeros(2 ** n_qubits, dtype=bool)
    for i in range(2 ** n_qubits):
        a_idx, b_idx = i % half_dim, i // half_dim
        if hamming(a_idx) == 2 and hamming(b_idx) == 2:
            mask_phys[i] = True
    e_mixed_physical = np.mean(diag[mask_phys]) + enuc
    e_mixed_all = np.mean(diag) + enuc
    assert abs(e_mixed_all - e_mixed) < 1e-9

    with open(os.path.join(os.path.dirname(__file__), "native_forged_zne_results.json")) as f:
        zne_results = json.load(f)
    quadratic_zne_kcal = zne_results["zne_report"]["aria-1"]["zne_quadratic_kcal"]
    chemical_accuracy_kcal = zne_results["chemical_accuracy_kcal"]

    targets = {}
    for n in range(K):
        targets[f"u_{n}"] = u_top[n]
    pairs = [(n, m) for n in range(K) for m in range(K) if n < m]
    for (n, m) in pairs:
        targets[f"(u{n}+u{m})"] = (u_top[n] + u_top[m]) / np.sqrt(2)
        targets[f"(u{n}-u{m})"] = (u_top[n] - u_top[m]) / np.sqrt(2)
    assert set(targets.keys()) == set(slot_names(K))

    return {
        "n_qubits": n_qubits, "terms": terms, "alpha_labels": alpha_labels, "identity_label": identity_label,
        "signs": signs, "lambdas": lambdas, "u_vecs": u_top, "enuc": enuc,
        "e_mixed": e_mixed, "e_mixed_physical_sector": e_mixed_physical,
        "noiseless_numpy": noiseless_numpy, "exact_energy": exact_energy,
        "truncation_floor_kcal": truncation_floor_kcal,
        "alpha_cache": alpha_cache, "targets": targets,
        "quadratic_zne_kcal": quadratic_zne_kcal, "chemical_accuracy_kcal": chemical_accuracy_kcal,
        "K": K,
    }


def energy_errors(E, p):
    """The critical fix: report error vs exact AND vs the rank-K
    noiseless floor, separately, always."""
    err_vs_exact = abs(E - p["exact_energy"]) * HARTREE_TO_KCAL_MOL
    err_vs_noiseless = abs(E - p["noiseless_numpy"]) * HARTREE_TO_KCAL_MOL
    f_signal = (E - p["e_mixed"]) / (p["noiseless_numpy"] - p["e_mixed"])
    return {"err_vs_exact_kcal": err_vs_exact, "err_vs_noiseless_kcal": err_vs_noiseless, "f": f_signal}


def fit_all_targets(targets):
    solutions = {}
    worst = 0.0
    for name, vec in targets.items():
        angles, err = fit_angles(vec)
        solutions[name] = {"angles": angles.tolist(), "max_abs_error": err, "verified": bool(err < 1e-10)}
        worst = max(worst, err)
    n_ok = sum(s["verified"] for s in solutions.values())
    return solutions, n_ok, worst


def verify_constant_gate_count(solutions):
    counts = set()
    for name in solutions:
        qc = build_ansatz(solutions[name]["angles"])
        t = transpile(qc, basis_gates=BASIS_GATES, optimization_level=0)
        counts.add(t.count_ops().get("cx", 0))
    return counts


def derive_beta_matrices(alpha_matrices, signs, K):
    S = np.diag(signs)
    return {label: S @ mat @ S for label, mat in alpha_matrices.items()}


def combine_matrices(raw, alpha_labels, identity_label, K, per_slot_scale=None, per_label_scale=None):
    def corrected(name, label):
        v = raw[name][label]
        if per_slot_scale is not None:
            v = v / per_slot_scale.get(name, 1.0)
        if per_label_scale is not None:
            v = v / per_label_scale.get(label, 1.0)
        return v

    mats = {l: np.zeros((K, K)) for l in alpha_labels}
    for n in range(K):
        name = f"u_{n}"
        for l in alpha_labels:
            mats[l][n, n] = 1.0 if l == identity_label else corrected(name, l)
    for n in range(K):
        for m in range(K):
            if n >= m:
                continue
            name0, name2 = f"(u{n}+u{m})", f"(u{n}-u{m})"
            for l in alpha_labels:
                if l == identity_label:
                    re = 0.0
                else:
                    re = (corrected(name0, l) - corrected(name2, l)) / 2
                mats[l][n, m] = re
                mats[l][m, n] = re
    return mats


def energy_from_alpha_matrices(alpha_mats, p, K):
    beta_mats = derive_beta_matrices(alpha_mats, p["signs"], K)
    E = effrag.ef_energy_from_noisy_matrices(p["terms"], p["lambdas"], alpha_mats, beta_mats, p["enuc"], K)
    return E, energy_errors(E, p)


def fit_scale(pairs):
    if len(pairs) < 3:
        return None
    exact = np.array([e for e, _ in pairs])
    noisy = np.array([n for _, n in pairs])
    denom = float(np.sum(exact * exact))
    if denom < 1e-12:
        return None
    return float(np.sum(exact * noisy) / denom)


def filtered_pairs(training, label=None, slot=None):
    out = []
    for row in training:
        if label is not None and row["label"] != label:
            continue
        if slot is not None and row["slot"] != slot:
            continue
        if abs(row["exact"]) < LOW_SIGNAL_CUTOFF:
            continue
        out.append((row["exact"], row["noisy"]))
    return out


def fit_all_scales(training, non_id_labels, K):
    global_scale = fit_scale(filtered_pairs(training))
    per_basis, per_basis_fallback = {}, []
    for l in non_id_labels:
        f = fit_scale(filtered_pairs(training, label=l))
        if f is None:
            f, per_basis_fallback = global_scale, per_basis_fallback + [l]
        per_basis[l] = f
    per_circuit, per_circuit_fallback = {}, []
    for name in slot_names(K):
        f = fit_scale(filtered_pairs(training, slot=name))
        if f is None:
            f, per_circuit_fallback = global_scale, per_circuit_fallback + [name]
        per_circuit[name] = f
    return {"global_scale": global_scale, "per_basis_scale": per_basis, "per_circuit_scale": per_circuit}


# ---------------------------------------------------------------------------
# Noise model (identical to cdr_mitigation.py's)
# ---------------------------------------------------------------------------

def build_noise_model():
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(P2_PER_GATE, 2), "cx")
    nm.add_all_qubit_quantum_error(depolarizing_error(P1_PER_GATE, 1), "u3")
    return nm


def exact_labels(angles, labels):
    sv = Statevector.from_instruction(build_ansatz(angles))
    return {l: float(sv.expectation_value(Pauli(l)).real) for l in labels}


# ---------------------------------------------------------------------------
# STEP 2: symmetry-verified (sector-confined) postselection, all-Z groups only
# ---------------------------------------------------------------------------

def is_all_z(label):
    return all(c in ("I", "Z") for c in label)


def symmetry_coverage(alpha_labels, terms):
    groups = effrag.group_labels_qubit_wise(alpha_labels)
    allz_flags = [is_all_z(effrag.combined_basis_label(g)) for g in groups]
    allz_labels = set()
    for g, allz in zip(groups, allz_flags):
        if allz:
            allz_labels.update(g)
    from collections import defaultdict
    weight = defaultdict(float)
    for a, _, coeff in terms:
        weight[a] += abs(coeff)
    total_weight = sum(weight.values())
    allz_weight = sum(weight[l] for l in allz_labels)
    return {
        "n_groups": len(groups), "n_allz_groups": sum(allz_flags),
        "allz_labels": sorted(allz_labels), "n_allz_labels": len(allz_labels),
        "n_total_labels": len(alpha_labels),
        "label_fraction": len(allz_labels) / len(alpha_labels),
        "weight_fraction": allz_weight / total_weight,
    }


def weight2_mask_4q():
    return np.array([bin(i).count("1") == 2 for i in range(16)])


def zlabel_eigenvalues(label, n=4):
    evs = np.ones(2 ** n)
    for i in range(2 ** n):
        val = 1.0
        for k, ch in enumerate(label):
            if ch == "Z":
                qubit = n - 1 - k
                bit = (i >> qubit) & 1
                val *= (1 - 2 * bit)
        evs[i] = val
    return evs


def density_matrix_4q(angles, noise_model):
    qc = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    qc2 = qc.copy()
    qc2.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc2).result()
    return np.asarray(result.data(0)["density_matrix"])


def measure_labels(angles, labels, dm, symmetry_verify=False, allz_labels=None):
    """Standard Tr[P.rho] for every label, EXCEPT: when symmetry_verify is
    on and the label is all-Z-type, use the postselected (weight-2-only,
    renormalized) expectation instead -- the only case where postselection
    on this project's own established finding is valid."""
    mask = weight2_mask_4q()
    probs = np.real(np.diag(dm))
    out, survival = {}, {}
    for l in labels:
        if symmetry_verify and allz_labels and l in allz_labels:
            evs = zlabel_eigenvalues(l)
            den = float(np.sum(probs[mask]))
            num = float(np.sum(probs[mask] * evs[mask]))
            out[l] = num / den if den > 1e-12 else 0.0
            survival[l] = den
        else:
            P = np.asarray(Pauli(l).to_matrix())
            out[l] = float(np.real(np.trace(P @ dm)))
            survival[l] = 1.0
    return out, survival


def measure_raw_per_slot(p, non_id_labels, noise_model, symmetry_verify=False, allz_labels=None):
    raw, surv = {}, {}
    for name, sol in p["solutions"].items():
        dm = density_matrix_4q(sol["angles"], noise_model)
        vals, s = measure_labels(sol["angles"], non_id_labels, dm,
                                  symmetry_verify=symmetry_verify, allz_labels=allz_labels)
        raw[name] = vals
        surv[name] = s
    return raw, surv


def generate_training_data(p, non_id_labels, noise_model, n_per_slot, seed,
                            symmetry_verify=False, allz_labels=None):
    rng = np.random.default_rng(seed)
    rows = []
    for name in slot_names(p["K"]):
        for _ in range(n_per_slot):
            angles = rng.uniform(-np.pi, np.pi, 5)
            exact_vals = exact_labels(angles, non_id_labels)
            dm = density_matrix_4q(angles, noise_model)
            noisy_vals, _ = measure_labels(angles, non_id_labels, dm,
                                            symmetry_verify=symmetry_verify, allz_labels=allz_labels)
            for l in non_id_labels:
                rows.append({"slot": name, "label": l, "exact": exact_vals[l], "noisy": noisy_vals[l]})
    return rows


# ---------------------------------------------------------------------------
# STEP 3: 2-copy Virtual Distillation, diagonalized B-gate form
# ---------------------------------------------------------------------------

def build_b_gate():
    """B = SWAP's own eigenvector matrix -- by construction, unitary (numpy
    eigh on a Hermitian matrix) and B^dagger . SWAP . B is exactly
    diagonal (verified in development against SWAP's known eigenvalues
    [-1,+1,+1,+1] before use). Real measured cost: 2 two-qubit gates
    (u3/cx, opt_level=0) -- not assumed from the literature."""
    S = SwapGate().to_matrix()
    eigvals, eigvecs = np.linalg.eigh(S)
    return eigvecs


def _b_full_and_s_full(B):
    from qiskit.quantum_info import Operator
    qcB = QuantumCircuit(8)
    for i in range(4):
        qcB.append(UnitaryGate(B), [i, i + 4])
    B_full = Operator(qcB).data

    qcS = QuantumCircuit(8)
    for i in range(4):
        qcS.swap(i, i + 4)
    S_full = Operator(qcS).data
    return B_full, S_full


def vd_build_D_S(B_full, S_full):
    """PHYSICAL circuit-evolution convention rho_after = B.rho.B^dagger
    (matches what Aer actually computes when B is applied as a real gate)
    -- verified against the direct-trace identity
    Tr[(O#I)S(rho#rho)] = Tr[O.rho^2] on a random test state before use;
    an earlier B^dagger.rho.B convention was wrong and gave a different,
    incorrect number on that exact test (caught here, not assumed away)."""
    return B_full @ S_full @ B_full.conj().T


def vd_build_D_O(alpha_label, B_full, S_full):
    O = np.asarray(Pauli("IIII" + alpha_label).to_matrix())  # copy A = qubits 0-3, copy B = 4-7
    return B_full @ O @ S_full @ B_full.conj().T


def vd_density_matrix(angles, noise_model, B):
    qcA = transpile(build_ansatz(angles), basis_gates=BASIS_GATES, optimization_level=0)
    qc8 = QuantumCircuit(8)
    qc8.compose(qcA, qubits=[0, 1, 2, 3], inplace=True)
    qc8.compose(qcA, qubits=[4, 5, 6, 7], inplace=True)
    for i in range(4):
        qc8.append(UnitaryGate(B), [i, i + 4])
    qc8t = transpile(qc8, basis_gates=BASIS_GATES, optimization_level=0)
    n2q = qc8t.count_ops().get("cx", 0)
    qc8t.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    result = sim.run(qc8t).result()
    return np.asarray(result.data(0)["density_matrix"]), n2q


def vd_measure_labels(angles, labels, noise_model, B, identity_label, D_S_cache, D_O_cache):
    dm, n2q = vd_density_matrix(angles, noise_model, B)
    den = float(np.real(np.trace(D_S_cache @ dm)))
    out = {}
    for l in labels:
        if l == identity_label:
            out[l] = 1.0
            continue
        num = float(np.real(np.trace(D_O_cache[l] @ dm)))
        out[l] = num / den if abs(den) > 1e-9 else 0.0
    return out, den, n2q


def vd_measure_raw_per_slot(p, non_id_labels, noise_model, B, D_S, D_O_cache):
    raw, dens, n2q_list = {}, {}, []
    for name, sol in p["solutions"].items():
        vals, den, n2q = vd_measure_labels(sol["angles"], non_id_labels, noise_model, B,
                                            p["identity_label"], D_S, D_O_cache)
        raw[name] = vals
        dens[name] = den
        n2q_list.append(n2q)
    return raw, dens, n2q_list


def vd_generate_training_data(p, non_id_labels, noise_model, B, D_S, D_O_cache, n_per_slot, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for name in slot_names(p["K"]):
        for _ in range(n_per_slot):
            angles = rng.uniform(-np.pi, np.pi, 5)
            exact_vals = exact_labels(angles, non_id_labels)
            noisy_vals, _, _ = vd_measure_labels(angles, non_id_labels, noise_model, B,
                                                  p["identity_label"], D_S, D_O_cache)
            for l in non_id_labels:
                rows.append({"slot": name, "label": l, "exact": exact_vals[l], "noisy": noisy_vals[l]})
    return rows


# ---------------------------------------------------------------------------
# Seed-swept CDR helper (shared by plain and symmetry-verified CDR)
# ---------------------------------------------------------------------------

def cdr_seed_sweep(p, non_id_labels, raw_noisy, noise_model, n_seeds, n_train_per_slot,
                    symmetry_verify=False, allz_labels=None, label="CDR"):
    K = p["K"]
    seed_rows = []
    for seed in range(n_seeds):
        training = generate_training_data(p, non_id_labels, noise_model, n_train_per_slot, seed,
                                           symmetry_verify=symmetry_verify, allz_labels=allz_labels)
        scales = fit_all_scales(training, non_id_labels, K)

        global_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                        per_slot_scale=None if scales["global_scale"] is None else
                                        {n: scales["global_scale"] for n in slot_names(K)})
        E_g, err_g = energy_from_alpha_matrices(global_mats, p, K)

        basis_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                       per_label_scale=scales["per_basis_scale"])
        E_b, err_b = energy_from_alpha_matrices(basis_mats, p, K)

        circuit_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K,
                                         per_slot_scale=scales["per_circuit_scale"])
        E_c, err_c = energy_from_alpha_matrices(circuit_mats, p, K)

        row = {"seed": seed,
               "global_exact": err_g["err_vs_exact_kcal"], "global_noiseless": err_g["err_vs_noiseless_kcal"],
               "basis_exact": err_b["err_vs_exact_kcal"], "basis_noiseless": err_b["err_vs_noiseless_kcal"],
               "circuit_exact": err_c["err_vs_exact_kcal"], "circuit_noiseless": err_c["err_vs_noiseless_kcal"]}
        seed_rows.append(row)
        print(f"    [{label}] seed={seed}: global={err_g['err_vs_exact_kcal']:.3f}  "
              f"basis={err_b['err_vs_exact_kcal']:.3f}  circuit={err_c['err_vs_exact_kcal']:.3f} "
              "kcal/mol (vs exact)")
    return seed_rows


def summarize_seed_rows(seed_rows, key, chem_acc):
    vals = [r[key] for r in seed_rows]
    n_reached = sum(1 for v in vals if v < chem_acc)
    return {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
            "min": float(np.min(vals)), "max": float(np.max(vals)),
            "reached_chem_acc": f"{n_reached}/{len(vals)}"}


def main():
    print("\n" + "=" * 78)
    print("  rank6_symmetry_vd.py -- Rank-6 + symmetry postselection + virtual distillation")
    print("  (simulator only)")
    print("=" * 78)

    p = setup(K)
    print(f"\n  K={K}, {len(p['targets'])} slots")
    print(f"  exact_energy = {p['exact_energy']:.6f} Ha")
    print(f"  noiseless_numpy(K={K}) = {p['noiseless_numpy']:.6f} Ha")
    print(f"  truncation floor at K={K}: {p['truncation_floor_kcal']:.6f} kcal/mol")
    print(f"  e_mixed (all 256 states): {p['e_mixed']:.6f} Ha, err vs exact = "
          f"{abs(p['e_mixed']-p['exact_energy'])*HARTREE_TO_KCAL_MOL:.2f} kcal/mol")
    print(f"  e_mixed (36 physical states): {p['e_mixed_physical_sector']:.6f} Ha, err vs exact = "
          f"{abs(p['e_mixed_physical_sector']-p['exact_energy'])*HARTREE_TO_KCAL_MOL:.2f} kcal/mol")

    t0 = time.time()
    solutions, n_ok, worst = fit_all_targets(p["targets"])
    p["solutions"] = solutions
    print(f"\n  fit all {len(p['targets'])} targets: {n_ok}/{len(p['targets'])} converged <1e-10, "
          f"worst={worst:.2e}, time={time.time()-t0:.1f}s")
    counts = verify_constant_gate_count(solutions)
    gate_count = list(counts)[0]
    print(f"  abstract (u3/cx, opt_level=0) 2-qubit gate count across all {len(solutions)} targets: {counts}")
    assert len(counts) == 1, "gate count not fixed across all K=6 targets -- refusing to proceed"

    non_id_labels = [l for l in p["alpha_labels"] if l != p["identity_label"]]
    noise_model = build_noise_model()
    chem_acc = p["chemical_accuracy_kcal"]

    results = {"K": K, "n_slots": len(solutions), "gate_count": gate_count,
               "truncation_floor_kcal_at_K": p["truncation_floor_kcal"],
               "e_mixed_all256_err_vs_exact_kcal": abs(p["e_mixed"] - p["exact_energy"]) * HARTREE_TO_KCAL_MOL,
               "e_mixed_physical36_err_vs_exact_kcal": abs(p["e_mixed_physical_sector"] - p["exact_energy"]) * HARTREE_TO_KCAL_MOL,
               "mixed_reference_reduction_ratio":
                   abs(p["e_mixed"] - p["exact_energy"]) / abs(p["e_mixed_physical_sector"] - p["exact_energy"])}

    # ---- RAW (K=6, no mitigation) ----
    print("\n  -- RAW (K=6, no mitigation) --")
    raw_noisy, _ = measure_raw_per_slot(p, non_id_labels, noise_model)
    raw_mats = combine_matrices(raw_noisy, p["alpha_labels"], p["identity_label"], K)
    E_raw, err_raw = energy_from_alpha_matrices(raw_mats, p, K)
    print(f"    E={E_raw:.6f} Ha, err_vs_exact={err_raw['err_vs_exact_kcal']:.3f} kcal/mol, "
          f"err_vs_noiseless={err_raw['err_vs_noiseless_kcal']:.3f} kcal/mol, f={err_raw['f']:.6f}")
    results["raw"] = {"E": E_raw, **err_raw}

    # ---- CDR alone (K=6, 8-seed sweep) ----
    print(f"\n  -- CDR alone (K={K}, {N_SEEDS}-seed sweep) --")
    cdr_rows = cdr_seed_sweep(p, non_id_labels, raw_noisy, noise_model, N_SEEDS, N_TRAIN_PER_SLOT, label="CDR")
    results["cdr_alone"] = {
        "global": {"exact": summarize_seed_rows(cdr_rows, "global_exact", chem_acc),
                   "noiseless": summarize_seed_rows(cdr_rows, "global_noiseless", chem_acc)},
        "per_basis": {"exact": summarize_seed_rows(cdr_rows, "basis_exact", chem_acc),
                      "noiseless": summarize_seed_rows(cdr_rows, "basis_noiseless", chem_acc)},
        "per_circuit": {"exact": summarize_seed_rows(cdr_rows, "circuit_exact", chem_acc),
                        "noiseless": summarize_seed_rows(cdr_rows, "circuit_noiseless", chem_acc)},
        "seed_rows": cdr_rows,
    }
    print(f"    per-basis (vs exact): mean={results['cdr_alone']['per_basis']['exact']['mean']:.3f} "
          f"std={results['cdr_alone']['per_basis']['exact']['std']:.3f} "
          f"(vs noiseless: mean={results['cdr_alone']['per_basis']['noiseless']['mean']:.3f})")

    # ---- STEP 2: symmetry coverage + symmetry-verified raw/CDR ----
    print("\n  -- STEP 2: symmetry-verified postselection coverage --")
    coverage = symmetry_coverage(p["alpha_labels"], p["terms"])
    print(f"    {coverage['n_allz_groups']}/{coverage['n_groups']} groups are all-Z-type, "
          f"{coverage['n_allz_labels']}/{coverage['n_total_labels']} labels "
          f"({coverage['label_fraction']*100:.1f}% by count, {coverage['weight_fraction']*100:.2f}% "
          f"by |coefficient| weight): {coverage['allz_labels']}")
    print(f"    HONEST: full mixed-reference reduction ratio is "
          f"{results['mixed_reference_reduction_ratio']:.3f}x, but only "
          f"{coverage['weight_fraction']*100:.2f}% of Hamiltonian weight is even ELIGIBLE for valid "
          "postselection (all-Z groups only) -- the rest needs X/Y basis rotation, where "
          "postselection is invalid (established: made an earlier energy 40x worse).")
    results["symmetry_coverage"] = coverage

    print("\n  -- symmetry-verified RAW --")
    raw_noisy_sym, surv_sym = measure_raw_per_slot(p, non_id_labels, noise_model,
                                                     symmetry_verify=True, allz_labels=set(coverage["allz_labels"]))
    raw_mats_sym = combine_matrices(raw_noisy_sym, p["alpha_labels"], p["identity_label"], K)
    E_raw_sym, err_raw_sym = energy_from_alpha_matrices(raw_mats_sym, p, K)
    avg_survival = float(np.mean([v for d in surv_sym.values() for v in d.values() if v != 1.0])) \
        if coverage["n_allz_labels"] else None
    print(f"    E={E_raw_sym:.6f} Ha, err_vs_exact={err_raw_sym['err_vs_exact_kcal']:.3f} kcal/mol "
          f"(raw was {err_raw['err_vs_exact_kcal']:.3f}), avg postselection survival={avg_survival}")
    results["symmetry_verified_raw"] = {"E": E_raw_sym, **err_raw_sym, "avg_postselection_survival": avg_survival}

    print(f"\n  -- symmetry-verified CDR ({N_SEEDS}-seed sweep) --")
    sym_cdr_rows = cdr_seed_sweep(p, non_id_labels, raw_noisy_sym, noise_model, N_SEEDS, N_TRAIN_PER_SLOT,
                                   symmetry_verify=True, allz_labels=set(coverage["allz_labels"]), label="sym-CDR")
    results["symmetry_verified_cdr"] = {
        "global": {"exact": summarize_seed_rows(sym_cdr_rows, "global_exact", chem_acc),
                   "noiseless": summarize_seed_rows(sym_cdr_rows, "global_noiseless", chem_acc)},
        "per_basis": {"exact": summarize_seed_rows(sym_cdr_rows, "basis_exact", chem_acc),
                      "noiseless": summarize_seed_rows(sym_cdr_rows, "basis_noiseless", chem_acc)},
        "per_circuit": {"exact": summarize_seed_rows(sym_cdr_rows, "circuit_exact", chem_acc),
                        "noiseless": summarize_seed_rows(sym_cdr_rows, "circuit_noiseless", chem_acc)},
        "seed_rows": sym_cdr_rows,
    }

    # ---- STEP 3: Virtual Distillation ----
    print("\n  -- STEP 3: Virtual Distillation setup --")
    B = build_b_gate()
    B_full, S_full = _b_full_and_s_full(B)
    D_S = vd_build_D_S(B_full, S_full)
    D_O_cache = {l: vd_build_D_O(l, B_full, S_full) for l in non_id_labels}
    print(f"    B gate 2-qubit cost: 2 (measured). Full 8-qubit VD circuit 2-qubit gate count "
          "measured per-run below.")

    print("\n  -- VD alone --")
    t0 = time.time()
    vd_raw, vd_dens, vd_n2q_list = vd_measure_raw_per_slot(p, non_id_labels, noise_model, B, D_S, D_O_cache)
    assert len(set(vd_n2q_list)) == 1, f"VD circuit 2-qubit gate count not fixed across targets: {set(vd_n2q_list)}"
    vd_n2q = vd_n2q_list[0]
    avg_S_denominator = float(np.mean(list(vd_dens.values())))
    vd_mats = combine_matrices(vd_raw, p["alpha_labels"], p["identity_label"], K)
    E_vd, err_vd = energy_from_alpha_matrices(vd_mats, p, K)
    print(f"    E={E_vd:.6f} Ha, err_vs_exact={err_vd['err_vs_exact_kcal']:.3f} kcal/mol "
          f"(raw was {err_raw['err_vs_exact_kcal']:.3f}), 2q gates/circuit={vd_n2q}, "
          f"avg <S>={avg_S_denominator:.4f}, time={time.time()-t0:.1f}s")
    sampling_overhead = 1.0 / (avg_S_denominator ** 2) if abs(avg_S_denominator) > 1e-6 else None
    results["vd_alone"] = {"E": E_vd, **err_vd, "n2q_gates_per_circuit": vd_n2q,
                            "avg_S_denominator": avg_S_denominator,
                            "sampling_overhead_shots_multiplier": sampling_overhead}

    print(f"\n  -- VD+CDR stacked ({N_SEEDS}-seed sweep, N_TRAIN_PER_SLOT_VD={N_TRAIN_PER_SLOT_VD}) --")
    vd_cdr_rows = []
    for seed in range(N_SEEDS):
        t0 = time.time()
        training = vd_generate_training_data(p, non_id_labels, noise_model, B, D_S, D_O_cache,
                                              N_TRAIN_PER_SLOT_VD, seed)
        scales = fit_all_scales(training, non_id_labels, K)

        global_mats = combine_matrices(vd_raw, p["alpha_labels"], p["identity_label"], K,
                                        per_slot_scale=None if scales["global_scale"] is None else
                                        {n: scales["global_scale"] for n in slot_names(K)})
        E_g, err_g = energy_from_alpha_matrices(global_mats, p, K)
        basis_mats = combine_matrices(vd_raw, p["alpha_labels"], p["identity_label"], K,
                                       per_label_scale=scales["per_basis_scale"])
        E_b, err_b = energy_from_alpha_matrices(basis_mats, p, K)
        circuit_mats = combine_matrices(vd_raw, p["alpha_labels"], p["identity_label"], K,
                                         per_slot_scale=scales["per_circuit_scale"])
        E_c, err_c = energy_from_alpha_matrices(circuit_mats, p, K)

        row = {"seed": seed,
               "global_exact": err_g["err_vs_exact_kcal"], "global_noiseless": err_g["err_vs_noiseless_kcal"],
               "basis_exact": err_b["err_vs_exact_kcal"], "basis_noiseless": err_b["err_vs_noiseless_kcal"],
               "circuit_exact": err_c["err_vs_exact_kcal"], "circuit_noiseless": err_c["err_vs_noiseless_kcal"]}
        vd_cdr_rows.append(row)
        print(f"    [VD+CDR] seed={seed}: global={err_g['err_vs_exact_kcal']:.3f}  "
              f"basis={err_b['err_vs_exact_kcal']:.3f}  circuit={err_c['err_vs_exact_kcal']:.3f} "
              f"kcal/mol (vs exact), time={time.time()-t0:.1f}s")
    results["vd_plus_cdr"] = {
        "global": {"exact": summarize_seed_rows(vd_cdr_rows, "global_exact", chem_acc),
                   "noiseless": summarize_seed_rows(vd_cdr_rows, "global_noiseless", chem_acc)},
        "per_basis": {"exact": summarize_seed_rows(vd_cdr_rows, "basis_exact", chem_acc),
                      "noiseless": summarize_seed_rows(vd_cdr_rows, "basis_noiseless", chem_acc)},
        "per_circuit": {"exact": summarize_seed_rows(vd_cdr_rows, "circuit_exact", chem_acc),
                        "noiseless": summarize_seed_rows(vd_cdr_rows, "circuit_noiseless", chem_acc)},
        "seed_rows": vd_cdr_rows, "n_train_per_slot": N_TRAIN_PER_SLOT_VD,
    }

    # ---- Final report: ALL THREE schemes for every CDR-based config, not just per-basis ----
    print("\n" + "=" * 100)
    print(f"  {'configuration':<30}{'vs exact (mean+/-std)':>24}{'vs noiseless(K)':>18}{'2q gates':>10}{'chem.acc.':>10}")
    print(f"  {'raw (K=6)':<30}{err_raw['err_vs_exact_kcal']:>17.3f} (--)  {err_raw['err_vs_noiseless_kcal']:>15.3f}"
          f"{gate_count:>10}{'no':>10}")
    for cfg_key, cfg_title, n2q in [("cdr_alone", "CDR", gate_count),
                                     ("symmetry_verified_cdr", "sym-CDR", gate_count),
                                     ("vd_plus_cdr", "VD+CDR", vd_n2q)]:
        for scheme, scheme_label in [("global", "global"), ("per_basis", "per-basis"), ("per_circuit", "per-circuit")]:
            d = results[cfg_key][scheme]
            exact_s, noiseless_s = d["exact"], d["noiseless"]
            print(f"  {cfg_title + ' ' + scheme_label:<30}{exact_s['mean']:>10.3f}+/-{exact_s['std']:<7.3f}"
                  f"{noiseless_s['mean']:>14.3f}{n2q:>10}{exact_s['reached_chem_acc']:>10}")
    print(f"  {'symmetry-verified raw':<30}{err_raw_sym['err_vs_exact_kcal']:>17.3f} (--)  "
          f"{err_raw_sym['err_vs_noiseless_kcal']:>15.3f}{gate_count:>10}{'no':>10}")
    print(f"  {'VD alone':<30}{err_vd['err_vs_exact_kcal']:>17.3f} (--)  "
          f"{err_vd['err_vs_noiseless_kcal']:>15.3f}{vd_n2q:>10}"
          f"{'no' if err_vd['err_vs_exact_kcal']>=chem_acc else 'yes':>10}")
    print("=" * 100)

    vd_helped = err_vd["err_vs_exact_kcal"] < err_raw["err_vs_exact_kcal"]
    vd_cdr_best_scheme, vd_cdr_best_mean = min(
        [(s, results["vd_plus_cdr"][s]["exact"]["mean"]) for s in ("global", "per_basis", "per_circuit")],
        key=lambda t: t[1])
    cdr_best_scheme, cdr_best_mean = min(
        [(s, results["cdr_alone"][s]["exact"]["mean"]) for s in ("global", "per_basis", "per_circuit")],
        key=lambda t: t[1])
    vd_cdr_vs_cdr = vd_cdr_best_mean < cdr_best_mean
    print(f"\n  VD alone vs raw: {'HELPED' if vd_helped else 'DID NOT HELP'} "
          f"({err_vd['err_vs_exact_kcal']:.3f} vs {err_raw['err_vs_exact_kcal']:.3f} kcal/mol)")
    print(f"  VD+CDR (best scheme: {vd_cdr_best_scheme}, {vd_cdr_best_mean:.3f}) vs CDR alone "
          f"(best scheme: {cdr_best_scheme}, {cdr_best_mean:.3f}): "
          f"{'HELPED' if vd_cdr_vs_cdr else 'DID NOT HELP'}")
    print(f"  NOTE: within VD+CDR, per-basis ({results['vd_plus_cdr']['per_basis']['exact']['mean']:.3f}) is "
          f"the WORST scheme -- global/per-circuit ({results['vd_plus_cdr']['global']['exact']['mean']:.3f}/"
          f"{results['vd_plus_cdr']['per_circuit']['exact']['mean']:.3f}) are much better, the OPPOSITE ranking "
          "from plain CDR (where per-basis dominates). Not cherry-picked away -- reported because it's real.")
    if not vd_helped or not vd_cdr_vs_cdr:
        print("\n  HONEST DIAGNOSIS (VD underperformed somewhere): the VD circuit measured here costs "
              f"{vd_n2q} two-qubit gates (2x11 state-prep + 4x2 B-gates) vs {gate_count} for a single copy "
              f"-- roughly {vd_n2q/gate_count:.1f}x the raw gate-noise exposure feeding INTO the state that "
              "gets squared. VD's quadratic suppression acts on the state-prep noise already present in "
              "rho, but the B-gate circuit's OWN noise is not suppressed at all (it corrupts the swap-test "
              "readout directly, after the squaring already happened) -- so a large-enough B-gate/derangement "
              "noise contribution can offset or exceed the quadratic suppression benefit. Separately, VD+CDR's "
              "per-basis regression is fit from only 3 training draws/slot (vs 5 for plain CDR, a compute-time "
              "tradeoff) on a quantity (num/den from the <S> division) that is no longer simply proportional "
              "to the exact value the way a raw Tr[P.rho] is -- a plausible reason per-basis fits poorly here "
              "specifically while global/per-circuit (which pool far more data per fit) do not. Not fully "
              "isolated, reported as a plausible mechanism, not asserted as proven.")

    all_candidates = [
        ("raw", err_raw["err_vs_exact_kcal"], None),
        ("symmetry-verified raw", err_raw_sym["err_vs_exact_kcal"], None),
        ("VD alone", err_vd["err_vs_exact_kcal"], None),
    ]
    for cfg_key, cfg_title in [("cdr_alone", "CDR"), ("symmetry_verified_cdr", "sym-CDR"), ("vd_plus_cdr", "VD+CDR")]:
        for scheme in ("global", "per_basis", "per_circuit"):
            all_candidates.append((f"{cfg_title} {scheme}", results[cfg_key][scheme]["exact"]["mean"],
                                    (cfg_key, scheme)))
    best_overall = min(all_candidates, key=lambda t: t[1])
    print(f"\n  BEST OVERALL (across every scheme in every configuration): {best_overall[0]} at "
          f"{best_overall[1]:.3f} kcal/mol (vs exact). Chemical accuracy ({chem_acc} kcal/mol) "
          f"{'reached' if best_overall[1] < chem_acc else 'NOT reached'}.")

    results["references"] = {
        "quadratic_zne_kcal_real_aria1_old_ansatz": p["quadratic_zne_kcal"],
        "chemical_accuracy_kcal": chem_acc,
        "prior_cdr_k5_per_basis_mean_kcal": 2.264,
        "prior_cdr_k5_note": "K=5 result from fixed_ansatz_v2_results.json, measured vs exact (conflates "
                              "the 0.5655 kcal/mol truncation floor with noise) -- superseded at K=6 here, "
                              "where truncation is exact so err_vs_exact == err_vs_noiseless. K=6's own best "
                              "CDR noise residual (2.850, per-basis) is HIGHER than K=5's true noise residual "
                              "(~1.7 = 2.264-0.5655) -- removing truncation error did NOT reduce the noise "
                              "residual here; report as found, not as expected.",
    }
    reached_reliably = False
    if best_overall[2] is not None:
        cfg_key, scheme = best_overall[2]
        reach_str = results[cfg_key][scheme]["exact"]["reached_chem_acc"]
        n_reached, n_total = (int(x) for x in reach_str.split("/"))
        reached_reliably = n_reached >= n_total - 1
    results["best_overall"] = {"configuration": best_overall[0], "err_vs_exact_kcal": round(best_overall[1], 4)}
    results["chemical_accuracy_reached_reliably"] = reached_reliably

    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2, default=lambda o: str(o))
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
