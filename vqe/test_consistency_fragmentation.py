#!/usr/bin/env python3
"""
test_consistency_fragmentation.py — two denoising ideas for molecular
tailoring, tested on H8 (d=1.0 A). All local, all free.

PART A — overlap extrapolation:
  Molecular tailoring's error shrinks as fragment overlap grows (already
  observed in covalent_fragment.py: block-4/overlap-2 has error 0.0107 Ha,
  block-6/overlap-4 has error 0.0018 Ha). If that shrinkage follows a
  predictable trend, we can EXTRAPOLATE from cheap small-overlap
  calculations to predict the large-overlap (near-exact) answer, without
  ever paying for the bigger calculation. This uses Aitken's Delta-squared
  process (a standard convergence-acceleration technique for a
  geometrically-converging sequence) on 3 overlap levels: no-overlap
  (block=2), overlap-2 (block=4), overlap-4 (block=6).

PART B — redundancy averaging:
  Independent (incoherent) hardware shot noise averages down like any
  other statistical noise: running the SAME measurement N times and
  averaging should shrink the noise-induced error by roughly 1/sqrt(N).
  This tests that directly: the overlap-2 (block=4) H8 fragmentation
  scheme, each fragment's energy sampled 20 times under a depolarizing
  noise model (qiskit-aer, local), then averaged and reassembled.

Both parts are entirely local simulation (qiskit statevector + qiskit-aer)
-- no Azure, no IBM, no network calls, no cost.

Run:
    python vqe/test_consistency_fragmentation.py
"""
import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
from hardware_fragmentation import _build_fragment_qop, hchain
from fragment_ansatz_test import exact_energy_no_penalty, _hf_bits
from molecules_real import ground_state

from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import ExcitationPreserving
from qiskit.primitives import StatevectorEstimator
from qiskit_aer.primitives import EstimatorV2 as AerEstimatorV2
from qiskit_aer.noise import NoiseModel, depolarizing_error

HARTREE_TO_KCAL_MOL = 627.5094740631
BASIS_GATES = ["u3", "cx"]
SHOTS = 4096
TWO_Q_ERROR = 0.01
N_SAMPLES = 20

FULL_EXACT_H8_HA = -4.307572
D = 1.0  # Angstrom, same bonded H-chain convention as covalent_fragment.py


# ============================================================================
# PART A -- overlap extrapolation (exact fragment energies, no noise)
# ============================================================================

def tailor_energy_exact(n_atoms, block, d):
    """
    Exact-diagonalization molecular tailoring (same scheme as
    covalent_fragment.py's tailor_energy: fragments slide by 2 atoms each
    step, overlap = block - 2), generalized to correctly handle a
    ZERO-overlap scheme (block=2 -> non-overlapping fragments, nothing to
    subtract) -- the original only ever ran with block=4/6, where the
    overlap is never empty.
    """
    shift = 2
    starts = list(range(0, n_atoms - block + 1, shift))
    fragments = [list(range(s, s + block)) for s in starts]

    covered = set()
    for frag in fragments:
        covered.update(frag)
    assert covered == set(range(n_atoms)), f"fragments {fragments} don't cover all {n_atoms} atoms"

    reassembled = 0.0
    for frag in fragments:
        geom = hchain(frag, d)
        _ehf, exact, _nq = ground_state(geom, nelec=len(frag), active=None)
        reassembled += exact

    for i in range(len(fragments) - 1):
        overlap_atoms = sorted(set(fragments[i]) & set(fragments[i + 1]))
        if not overlap_atoms:
            continue  # no-overlap scheme -- nothing to double-count, nothing to subtract
        geom = hchain(overlap_atoms, d)
        _ehf, exact, _nq = ground_state(geom, nelec=len(overlap_atoms), active=None)
        reassembled -= exact

    return reassembled


def aitken_extrapolate(y1, y2, y3):
    """
    Aitken's Delta-squared process: given 3 terms of a sequence converging
    geometrically toward a limit, estimate that limit directly:
        y_inf = y1 - (y2 - y1)^2 / (y3 - 2*y2 + y1)
    Exactly solves for E_inf, A, r in the model y(overlap) = E_inf - A*r^n
    from 3 data points (a fully-determined fit, not a least-squares one).
    """
    denom = (y3 - 2 * y2 + y1)
    if denom == 0:
        return None
    return y1 - (y2 - y1) ** 2 / denom


def run_part_a():
    print("\n" + "=" * 80)
    print("  PART A -- overlap extrapolation (exact fragment energies)")
    print("=" * 80)

    schemes = [("no-overlap (block=2)", 2, 0), ("overlap-2 (block=4)", 4, 2), ("overlap-4 (block=6)", 6, 4)]
    energies = []
    for label, block, overlap in schemes:
        e = tailor_energy_exact(8, block, D)
        err = abs(e - FULL_EXACT_H8_HA)
        print(f"  {label:24s}: E = {e:.6f} Ha  (error {err:.6f} Ha)")
        energies.append(e)

    y1, y2, y3 = energies
    e_extrap = aitken_extrapolate(y1, y2, y3)
    err_extrap = abs(e_extrap - FULL_EXACT_H8_HA) if e_extrap is not None else None

    print(f"\n  Aitken extrapolation (from the 3 points above) -> E_inf = {e_extrap:.6f} Ha "
          f"(error {err_extrap:.6f} Ha)" if e_extrap is not None else "\n  Aitken extrapolation failed (degenerate sequence)")

    worst_individual_err = max(abs(e - FULL_EXACT_H8_HA) for e in energies)
    best_individual_err = min(abs(e - FULL_EXACT_H8_HA) for e in energies)
    beats_all = err_extrap is not None and err_extrap < best_individual_err

    print(f"\n  Best individual scheme error  = {best_individual_err:.6f} Ha")
    print(f"  Worst individual scheme error = {worst_individual_err:.6f} Ha")
    print(f"  Extrapolated error             = {err_extrap:.6f} Ha" if err_extrap is not None else "  n/a")
    print(f"  -> Extrapolation beats every individual scheme: {beats_all}")
    print("=" * 80)

    return {
        "full_exact_h8_ha": FULL_EXACT_H8_HA,
        "schemes": [
            {"label": label, "block": block, "overlap": overlap, "energy_ha": round(e, 6),
             "error_ha": round(abs(e - FULL_EXACT_H8_HA), 6)}
            for (label, block, overlap), e in zip(schemes, energies)
        ],
        "aitken_extrapolated_energy_ha": round(e_extrap, 6) if e_extrap is not None else None,
        "aitken_extrapolated_error_ha": round(err_extrap, 6) if err_extrap is not None else None,
        "best_individual_error_ha": round(best_individual_err, 6),
        "worst_individual_error_ha": round(worst_individual_err, 6),
        "extrapolation_beats_all_individual_schemes": bool(beats_all),
    }


# ============================================================================
# PART B -- redundancy averaging (noisy ansatz measurements, N=20 samples)
# ============================================================================

def build_noise_model():
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(TWO_Q_ERROR, 2), "cx")
    return nm


def optimize_fragment_locally(qop_bare, enuc, exact, n_qubits, nelec, reps=2, n_restarts=5, maxiter=500):
    """Same approach as hardware_covalent.py / fragment_ansatz_test.py / fragmentation_noise_prediction.py."""
    bits = _hf_bits(n_qubits, nelec)
    base_ansatz = ExcitationPreserving(n_qubits, reps=reps)
    n_params = base_ansatz.num_parameters
    estimator = StatevectorEstimator()

    def make_circuit(params):
        qc = QuantumCircuit(n_qubits)
        for i, b in enumerate(bits):
            if b:
                qc.x(i)
        qc.compose(base_ansatz.assign_parameters(params), inplace=True)
        return qc

    def objective(params):
        qc = make_circuit(params)
        result = estimator.run([(qc, qop_bare)]).result()
        return float(result[0].data.evs) + enuc

    rng = np.random.default_rng(42)
    best_params, best_energy = None, None
    for r in range(n_restarts):
        x0 = np.zeros(n_params) if r == 0 else rng.uniform(-np.pi, np.pi, size=n_params)
        res = minimize(
            objective, x0=x0, method="L-BFGS-B",
            options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-10},
        )
        if best_energy is None or res.fun < best_energy:
            best_energy, best_params = float(res.fun), res.x
        error = abs(best_energy - exact)
        print(f"        restart {r + 1}/{n_restarts}: E = {res.fun:.6f} Ha "
              f"(best {best_energy:.6f} Ha, error {error:.6f} Ha)")
        if error < 0.0005:
            break

    qc = make_circuit(best_params)
    return best_energy, abs(best_energy - exact), qc


def sample_noisy_energies(qc, qop_bare, enuc, noise_model, n_samples):
    """Run the SAME optimized circuit N independent times under shot noise -- each call is a fresh random draw."""
    qc_t = transpile(qc, basis_gates=BASIS_GATES, optimization_level=1)
    samples = []
    # IMPORTANT: AerEstimatorV2 computes the EXACT expectation value under
    # the noise channel by default (via Aer's save_expectation_value) --
    # `shots` in run_options does NOT drive Monte Carlo sampling for this
    # code path, so without default_precision every "sample" comes back
    # bit-for-bit identical (confirmed: averaging 20 identical numbers
    # changes nothing). Setting default_precision to the expected
    # standard error for SHOTS shots (~1/sqrt(shots)) makes the primitive
    # inject genuine call-to-call statistical noise, which is what we
    # actually need to test whether averaging shrinks it.
    precision = 1 / np.sqrt(SHOTS)
    for _ in range(n_samples):
        estimator = AerEstimatorV2(options={
            "backend_options": {"noise_model": noise_model},
            "default_precision": precision,
        })
        result = estimator.run([(qc_t, qop_bare)]).result()
        samples.append(float(result[0].data.evs) + enuc)
    return samples


def process_fragment_for_part_b(name, indices, d, noise_model):
    geom = hchain(indices, d)
    nelec = len(indices)
    print(f"\n  --- Fragment '{name}' (atoms {indices}, {nelec} electrons) ---")

    qop_bare, _qop_penalized, enuc = _build_fragment_qop(geom, nelec)
    n_qubits = qop_bare.num_qubits
    exact = exact_energy_no_penalty(qop_bare, enuc)
    print(f"    qubits = {n_qubits}, exact = {exact:.6f} Ha")

    print("    Optimizing ExcitationPreserving(reps=2) noiselessly (L-BFGS-B)...")
    noiseless_energy, noiseless_error, qc = optimize_fragment_locally(qop_bare, enuc, exact, n_qubits, nelec)
    print(f"    -> noiseless best: E = {noiseless_energy:.6f} Ha, error = {noiseless_error:.6f} Ha")

    print(f"    Sampling {N_SAMPLES} independent noisy measurements...")
    samples = sample_noisy_energies(qc, qop_bare, enuc, noise_model, N_SAMPLES)
    avg = float(np.mean(samples))
    print(f"    -> single sample (#1) = {samples[0]:.6f} Ha, "
          f"average of {N_SAMPLES} = {avg:.6f} Ha")

    return {
        "name": name, "indices": indices, "n_qubits": n_qubits, "exact_energy_ha": exact,
        "samples_ha": samples, "single_sample_ha": samples[0], "average_ha": avg,
    }


def run_part_b():
    print("\n" + "=" * 80)
    print("  PART B -- redundancy averaging (overlap-2 / block=4 scheme, N=20 samples)")
    print("=" * 80)

    noise_model = build_noise_model()
    print(f"  Noise model: depolarizing_error({TWO_Q_ERROR}, 2) on every CX, {SHOTS} shots/sample")

    # overlap-2 (block=4) scheme for H8: fragments [0-3],[2-5],[4-7], overlaps [2,3] and [4,5]
    blocks = [("block1 (0-3)", [0, 1, 2, 3]), ("block2 (2-5)", [2, 3, 4, 5]), ("block3 (4-7)", [4, 5, 6, 7])]
    overlaps = [("overlap1 (2-3)", [2, 3]), ("overlap2 (4-5)", [4, 5])]

    block_results = [process_fragment_for_part_b(name, idx, D, noise_model) for name, idx in blocks]
    overlap_results = [process_fragment_for_part_b(name, idx, D, noise_model) for name, idx in overlaps]

    exact_h8 = (sum(b["exact_energy_ha"] for b in block_results)
                - sum(o["exact_energy_ha"] for o in overlap_results))

    e_single = (sum(b["single_sample_ha"] for b in block_results)
                - sum(o["single_sample_ha"] for o in overlap_results))
    e_averaged = (sum(b["average_ha"] for b in block_results)
                  - sum(o["average_ha"] for o in overlap_results))

    err_single = abs(e_single - exact_h8)
    err_averaged = abs(e_averaged - exact_h8)
    shrink_factor = err_single / err_averaged if err_averaged > 0 else float("inf")

    print("\n" + "-" * 80)
    print("  REASSEMBLED H8 (exact fragment sum, for reference within this noisy scheme)")
    print(f"  Exact (fragment-sum) H8       = {exact_h8:.6f} Ha")
    print(f"  Single noisy sample H8        = {e_single:.6f} Ha  (error {err_single:.6f} Ha)")
    print(f"  Averaged ({N_SAMPLES}x) H8              = {e_averaged:.6f} Ha  (error {err_averaged:.6f} Ha)")
    print(f"  Error shrink factor (single/averaged) = {shrink_factor:.3f}x  "
          f"(expected ~sqrt({N_SAMPLES}) = {np.sqrt(N_SAMPLES):.3f}x)")
    print("=" * 80)

    return {
        "exact_fragment_sum_h8_ha": round(exact_h8, 6),
        "shots_per_sample": SHOTS,
        "n_samples": N_SAMPLES,
        "two_qubit_depolarizing_error": TWO_Q_ERROR,
        "fragments": {
            "blocks": [
                {"name": b["name"], "indices": b["indices"], "n_qubits": b["n_qubits"],
                 "exact_energy_ha": round(b["exact_energy_ha"], 6),
                 "single_sample_ha": round(b["single_sample_ha"], 6),
                 "average_ha": round(b["average_ha"], 6)}
                for b in block_results
            ],
            "overlaps": [
                {"name": o["name"], "indices": o["indices"], "n_qubits": o["n_qubits"],
                 "exact_energy_ha": round(o["exact_energy_ha"], 6),
                 "single_sample_ha": round(o["single_sample_ha"], 6),
                 "average_ha": round(o["average_ha"], 6)}
                for o in overlap_results
            ],
        },
        "single_sample_h8_ha": round(e_single, 6),
        "single_sample_error_ha": round(err_single, 6),
        "averaged_h8_ha": round(e_averaged, 6),
        "averaged_error_ha": round(err_averaged, 6),
        "error_shrink_factor": round(shrink_factor, 4),
        "expected_shrink_factor_sqrt_n": round(float(np.sqrt(N_SAMPLES)), 4),
    }


# ============================================================================
# Main
# ============================================================================

def main():
    print("\n" + "=" * 80)
    print("  Consistency / denoising tests for fragmentation -- H8 chain")
    print("  (LOCAL simulation only -- no Azure, no IBM, no cost)")
    print("=" * 80)

    part_a = run_part_a()
    part_b = run_part_b()

    results = {"part_a_overlap_extrapolation": part_a, "part_b_redundancy_averaging": part_b}
    out_path = os.path.join(os.path.dirname(__file__), "consistency_fragmentation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out_path}\n")


if __name__ == "__main__":
    main()
