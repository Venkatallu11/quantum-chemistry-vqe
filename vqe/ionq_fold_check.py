#!/usr/bin/env python3
"""
ionq_fold_check.py — does gate folding actually scale noise on IonQ's
cloud simulator at all? (free tier only: ionq_simulator, never ionq_qpu)
============================================================================
ionq_run.py's ZNE result was a genuine negative finding: the H4 EF energy
under noise_model="aria-1" was statistically flat across fold=1/3/5
(125.07 / 124.75 / 126.49 kcal/mol) -- no detectable noise-vs-fold trend,
so extrapolation to fold=0 gave no improvement. That result came from a
185-term Hamiltonian sum, where per-term shot noise (even at 100k shots)
could plausibly be masking a real but small underlying trend.

This script isolates the question cleanly: take the simplest possible
circuit with a KNOWN exact answer -- a 2-qubit Bell state, H+CX, measured
in the XX basis (exact expectation = +1.0 for any fold, since folding a
gate with its own inverse never changes the ideal unitary) -- and use
ionq_run.py's own fold_circuit() to fold the CX gate across a WIDE range
of fold factors (not just 1/3/5), with a much higher shot count than the
full EF run could afford per circuit (only 2 circuits total here, vs.
~3330 for the full run, so precision is cheap).

Two noise_model runs:
  "ideal"  -- a control. Must stay at exactly +1.0 for every fold factor;
              if it doesn't, something is wrong with the folding mechanism
              itself, not the noise model.
  "aria-1" -- the real question. If <XX> decays with fold (as it would
              under any simple per-gate depolarizing/coherent error
              model), gate folding IS a usable noise-scaling knob here and
              the flat EF result was a shot-noise/precision problem, not a
              fundamental limitation. If it stays flat even out to a much
              higher fold factor than the EF run tested, that's stronger
              evidence aria-1 is a fixed noise *profile* that genuinely
              doesn't respond to stretched gate count -- worth knowing
              BEFORE spending money on ionq_qpu.

Run:
    python vqe/ionq_fold_check.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qiskit.circuit import QuantumCircuit

from ionq_backend import connect_provider, get_simulator
from ionq_run import fold_circuit, basis_change, pauli_expectation

FOLD_FACTORS = [1, 3, 5, 7, 9, 11, 15, 21, 31, 41, 61, 81]
SHOTS = 1_000_000  # only 2 circuits total per noise_model here (vs ~3330 for
# the full EF run), so we can afford far more shots -- per-fold noise floor
# at 1e6 shots is ~0.001, an order of magnitude tighter than what the full
# run could afford per term.
LABEL = "XX"


def bell_circuit():
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    return qc


def folded_measurement_circuit(fold):
    qc = fold_circuit(bell_circuit(), fold)
    basis_change(qc, LABEL)
    qc.measure_all()
    return qc


def run_sweep(backend, noise_model):
    circuits = [folded_measurement_circuit(f) for f in FOLD_FACTORS]
    job = backend.run(circuits, noise_model=noise_model, shots=SHOTS)
    result = job.result()
    probs_list = result.get_probabilities()
    if not isinstance(probs_list, list):
        probs_list = [probs_list]
    energies = [pauli_expectation(dict(p), LABEL) for p in probs_list]
    return energies


def fit_exponential_decay(folds, energies):
    """Fit E(fold) = A * exp(-lambda * fold) + offset via a linear fit on
    log(|E - offset|) vs fold, offset pinned to 0 since the ideal/no-noise
    value is exactly +1.0 and a simple depolarizing-type channel decays
    toward 0, not toward the ideal value."""
    folds = np.asarray(folds, dtype=float)
    energies = np.asarray(energies, dtype=float)
    safe = np.clip(np.abs(energies), 1e-9, None)
    slope, intercept = np.polyfit(folds, np.log(safe), 1)
    return {
        "log_linear_slope": float(slope),
        "log_linear_intercept": float(intercept),
        "implied_per_application_retention": float(np.exp(slope)),
    }


def main():
    print("\n" + "=" * 70)
    print("  IonQ fold-check -- does folding scale noise on ionq_simulator?")
    print("=" * 70)

    provider = connect_provider()
    backend = get_simulator(provider)
    print(f"\n  Connected. Target backend: {backend.name}")
    print(f"  Circuit: 2-qubit Bell (H, CX), measured in XX basis, exact = +1.0")
    print(f"  Fold factors: {FOLD_FACTORS}")
    print(f"  Shots per circuit: {SHOTS}")

    print("\n  -- noise_model=ideal (control -- must stay at +1.0) --")
    ideal_energies = run_sweep(backend, "ideal")
    for f, e in zip(FOLD_FACTORS, ideal_energies):
        print(f"    fold={f:3d}: <XX> = {e:.6f}")

    print("\n  -- noise_model=aria-1 (the real question) --")
    aria_energies = run_sweep(backend, "aria-1")
    for f, e in zip(FOLD_FACTORS, aria_energies):
        print(f"    fold={f:3d}: <XX> = {e:.6f}")

    ideal_ok = all(abs(e - 1.0) < 1e-3 for e in ideal_energies)
    print(f"\n  Control check (ideal stays at +1.0 for every fold): "
          f"{'PASS' if ideal_ok else 'FAIL -- folding mechanism itself may be broken'}")

    aria_range = max(aria_energies) - min(aria_energies)
    aria_trend = aria_energies[0] - aria_energies[-1]  # positive = decaying as expected
    print(f"\n  aria-1 <XX> range across all folds: {aria_range:.6f}")
    print(f"  aria-1 <XX> drop from fold={FOLD_FACTORS[0]} to fold={FOLD_FACTORS[-1]}: "
          f"{aria_trend:.6f}")

    fit = None
    if aria_range > 5 * SHOTS ** -0.5:
        fit = fit_exponential_decay(FOLD_FACTORS, aria_energies)
        print(f"\n  Detected a real trend -- exponential fit:")
        print(f"    implied per-application retention = "
              f"{fit['implied_per_application_retention']:.6f}")
        print(f"    (i.e. each extra folded gate application multiplies |<XX>| by this factor)")
        verdict = ("Gate folding DOES scale noise on aria-1 -- the flat EF result was "
                   "likely a shot-noise/precision problem at the full-Hamiltonian scale, "
                   "not a fundamental limitation. Worth revisiting EF+ZNE with more shots "
                   "or a narrower fold range before considering real hardware.")
    else:
        verdict = ("No detectable trend even out to fold={} (range {:.6f}, noise floor "
                   "~{:.6f} at {} shots). Gate folding does not appear to be a usable "
                   "ZNE knob against IonQ's aria-1 simulated noise profile -- consistent "
                   "with aria-1 being a fixed noise profile rather than a per-gate "
                   "depolarizing model that responds to stretched circuit depth.").format(
                       FOLD_FACTORS[-1], aria_range, 3 * SHOTS ** -0.5, SHOTS)

    print(f"\n  VERDICT: {verdict}")

    results = {
        "circuit": "2-qubit Bell (H, CX), measured in XX basis",
        "exact_value": 1.0,
        "fold_factors": FOLD_FACTORS,
        "shots": SHOTS,
        "ideal_energies": [round(float(x), 6) for x in ideal_energies],
        "ideal_control_pass": ideal_ok,
        "aria1_energies": [round(float(x), 6) for x in aria_energies],
        "aria1_range": round(float(aria_range), 6),
        "aria1_fit": fit,
        "verdict": verdict,
    }
    out = os.path.join(os.path.dirname(__file__), "ionq_fold_check_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


if __name__ == "__main__":
    main()
