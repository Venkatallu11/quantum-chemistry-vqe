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
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qiskit.circuit import QuantumCircuit
from qiskit_ionq.ionq_gates import GPI2Gate, MSGate, ZZGate

from ionq_backend import connect_provider, get_simulator, get_native_simulator
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


def abstract_fold_check():
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


# ---------------------------------------------------------------------------
# Native-gate experiment: does the compiler cancel abstract-gate folds?
# ---------------------------------------------------------------------------
# Hypothesis (unverified attribution, treated as a hypothesis to test, not
# fact -- see conversation): IonQ's cloud-side compiler may be cancelling
# the folded G G^-1 pairs before execution when circuits are submitted as
# abstract gates. Confirmed separately that qiskit-ionq's OWN optimizer
# plugin (TrappedIonOptimizerPlugin) ships FuseConsecutiveZZ/FuseConsecutiveMS
# passes that would do exactly this if invoked -- plausible mechanism, even
# though this repo's own code never calls that plugin (no backend= passed to
# transpile(), nothing re-transpiles after folding). "native" gateset
# submission is meant to bypass generic circuit-level recompilation, so this
# tests whether that holds.
#
# NATIVE_2Q_BY_FAMILY (qiskit-ionq's own mapping): "aria" -> ms, "forte" -> zz.
# Two different probe circuits are used because they entangle differently:
#   zz-probe (GPI2(1/4) on both qubits + ZZ(1/4)) is a clean XX eigenstate (+1).
#   ms-probe (bare MS(0,0,1/4) on |00>) is NOT an XX eigenstate (it has an i
#     relative phase between |00> and |11>) but IS a clean ZZ eigenstate (+1)
#     -- verified via Statevector before use, not assumed. ZZ needs no
#     basis-change gates at all (direct Z-basis measurement).
# Both verified locally: folding preserves the exact ideal value (+1.0) at
# every fold factor tested, with the correct native-gate count present in
# the circuit, before any real submission.
#
# PASS/FAIL: ideal staying at 1.0 proves nothing either way (a cancelled
# fold also returns 1.0 under noiseless execution) -- only aria-1/forte-1
# decaying with fold is evidence gateset=native bypasses the cancellation.

NATIVE_FOLD_FACTORS = [1, 3, 5, 9, 21, 81]
NATIVE_SHOTS = 1_000_000


def native_probe_zz():
    qc = QuantumCircuit(2)
    qc.append(GPI2Gate(0.25), [0])
    qc.append(GPI2Gate(0.25), [1])
    qc.append(ZZGate(0.25), [0, 1])
    return qc


def native_probe_ms():
    qc = QuantumCircuit(2)
    qc.append(MSGate(0, 0, 0.25), [0, 1])
    return qc


def fold_native_2q(qc, fold, gate_name):
    """Fold ONLY the single 2-qubit native gate: G -> G (G^-1 G)^reps.
    EXPLICIT correct inverse -- verified that the generic Gate.inverse()
    on ZZGate/MSGate returns a mis-parametrized "zz_dg"/"ms_dg" gate with
    the SAME (not negated) params, not a valid native-gateset instruction
    and not actually the inverse. Real bug caught by testing, not assumed."""
    if fold == 1:
        return qc.copy()
    assert fold % 2 == 1, "fold factor must be odd"
    reps = (fold - 1) // 2
    folded = qc.copy_empty_like()
    for instr in qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        folded.append(op, qargs, cargs)
        if op.name == gate_name:
            inv = ZZGate(-op.params[0]) if gate_name == "zz" else \
                MSGate(op.params[0], op.params[1], -op.params[2])
            for _ in range(reps):
                folded.append(inv, qargs, cargs)
                folded.append(op, qargs, cargs)
    return folded


def native_measurement_circuit(base, gate_name):
    qc = base.copy()
    if gate_name == "zz":  # XX observable needs an X-basis rotation
        qc.append(GPI2Gate(0.25), [0])
        qc.append(GPI2Gate(0.25), [1])
    # ms-probe: ZZ observable, direct Z-basis measurement, no rotation needed
    qc.measure_all()
    return qc


def native_run_sweep(backend, noise_model, gate_name):
    probe = native_probe_zz() if gate_name == "zz" else native_probe_ms()
    circuits = [
        native_measurement_circuit(fold_native_2q(probe, f, gate_name), gate_name)
        for f in NATIVE_FOLD_FACTORS
    ]
    job = backend.run(circuits, noise_model=noise_model, shots=NATIVE_SHOTS)
    result = job.result()
    probs_list = result.get_probabilities()
    if not isinstance(probs_list, list):
        probs_list = [probs_list]
    label = "XX" if gate_name == "zz" else "ZZ"
    return [pauli_expectation(dict(p), label) for p in probs_list]


def surviving_fraction(values, ideal_value=1.0, mixed_value=0.0):
    """f = (measured - mixed) / (ideal - mixed). Both probes are traceless
    2-qubit Paulis (XX, ZZ) so the fully-decohered/mixed-state expectation
    is exactly 0 -- same formula convention as ionq_tailoring.py."""
    return [(v - mixed_value) / (ideal_value - mixed_value) for v in values]


def run_native_experiment():
    print("\n" + "=" * 70)
    print("  Native-gate fold check -- does gateset=native bypass fold cancellation?")
    print("=" * 70)

    provider = connect_provider()
    native_backend = get_native_simulator(provider)
    print(f"\n  Connected. Target backend: {native_backend.name} (gateset=native)")
    print(f"  Fold factors: {NATIVE_FOLD_FACTORS}, shots={NATIVE_SHOTS}")

    results = {"fold_factors": NATIVE_FOLD_FACTORS, "shots": NATIVE_SHOTS, "models": {}}

    for noise_model in ("ideal", "aria-1", "forte-1"):
        gate_name = "ms" if noise_model in ("ideal", "aria-1") else "zz"
        print(f"\n  -- noise_model={noise_model} (native {gate_name.upper()} gate) --")
        values = native_run_sweep(native_backend, noise_model, gate_name)
        f_values = surviving_fraction(values)
        for fold, v, f in zip(NATIVE_FOLD_FACTORS, values, f_values):
            print(f"    fold={fold:3d}: measured={v:.6f}  f={f:.6f}")
        results["models"][noise_model] = {
            "gate": gate_name,
            "values": [round(float(v), 6) for v in values],
            "f": [round(float(v), 6) for v in f_values],
        }

    ideal_vals = results["models"]["ideal"]["values"]
    ideal_ok = all(abs(v - 1.0) < 1e-3 for v in ideal_vals)
    print(f"\n  ideal control (must stay ~1.0 at every fold): "
          f"{'PASS' if ideal_ok else 'FAIL'} -- note: even a CANCELLED fold "
          "returns 1.0 under ideal, so this alone proves nothing either way.")

    verdicts = {}
    for noise_model in ("aria-1", "forte-1"):
        f_vals = results["models"][noise_model]["f"]
        frange = max(f_vals) - min(f_vals)
        noise_floor = 3 * NATIVE_SHOTS ** -0.5
        decaying = frange > noise_floor and f_vals[-1] < f_vals[0] - noise_floor
        verdicts[noise_model] = "PASS (f decays with fold)" if decaying else \
            "FAIL (f pinned near baseline, same flat pattern as before)"
        print(f"  {noise_model}: f range={frange:.6f}, f[0]={f_vals[0]:.6f}, "
              f"f[-1]={f_vals[-1]:.6f} -> {verdicts[noise_model]}")

    results["ideal_control_pass"] = ideal_ok
    results["verdicts"] = verdicts

    out = os.path.join(os.path.dirname(__file__), "ionq_native_fold_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {out}\n")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--native", action="store_true",
                         help="run the native-gate fold-cancellation experiment "
                              "(vqe/ionq_native_fold_results.json) instead of the "
                              "original abstract-gate one")
    args = parser.parse_args()
    if args.native:
        return run_native_experiment()
    return abstract_fold_check()


if __name__ == "__main__":
    main()
