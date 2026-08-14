#!/usr/bin/env python3
"""
task31b_iyyi_convergence.py -- iteration 31, Task B. BLOCKING, cheap.
For slot (u0+u1), label IYYI: analytic PEC gave ~-0.996/-0.962 while
literal quasi-probability twirling (Task 30B addendum, N_MC=8) gave
~+0.757/+1.000 -- nearly the entire physical range, on one observable.
IYYI carries real Hamiltonian weight (121.6 kcal/mol, rank 10/37) -- not
negligible.
============================================================================
DESIGN, cheap by construction: IYYI's measurement group is just
['XYYX','IYYI'] -- ONE basis rotation, so N_MC=128 real twirled circuits
for this ONE (slot, group) pair is enough to test N_MC in
{8,16,32,64,128} CUMULATIVELY (running estimate using the first N_MC of
128 total draws) without separate submissions per level -- 128 real
circuits total, not thousands.

Uses Task 31A's corrected, position-consensus p2(zz)=0.0146 (9/11
positions, real CI-backed) instead of Task 30B's single-point 0.014.

Run:
    python vqe/task31b_iyyi_convergence.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from loop_pec import pec_inverse_weights, gamma_factor
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts
from native_stateprep import to_native
from task28d_all_gate_zne import optimized_native_circuit
from task30b_pec_application import analytic_A_and_B
from qiskit.quantum_info import Pauli

K = 6
SHOTS = 100_000
N_MC_LEVELS = [8, 16, 32, 64, 128]
N_MC_MAX = 128
GATE_NAME = "zz"
SLOT = "(u0+u1)"
GROUP = ["XYYX", "IYYI"]
TARGET_LABEL = "IYYI"
P2_ZZ_CORRECTED = 0.0146  # Task 31A's 9/11-position consensus, CI-backed
GPI2_FALLBACK_NOTE = "GPi mean (Task 30B fallback), GPi2 itself only bounded |p|<0.21 per Task 31A"
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task31b_iyyi_convergence.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task31b_iyyi_convergence_results.json")


def sample_twirled_circuit(base_qc, p2, p1_gpi, p1_gpi2, rng):
    qc = base_qc.copy_empty_like()
    total_sign = 1
    for instr in base_qc.data:
        op, qargs, cargs = instr.operation, instr.qubits, instr.clbits
        qc.append(op, qargs, cargs)
        if op.name == GATE_NAME:
            weights = pec_inverse_weights(p2, 2)
        elif op.name == "gpi":
            weights = pec_inverse_weights(p1_gpi, 1)
        elif op.name == "gpi2":
            weights = pec_inverse_weights(p1_gpi2, 1)
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
        n = len(chosen_label)
        if chosen_label != "I" * n:
            qc.append(Pauli(chosen_label).to_instruction(), qargs)
    qc = to_native(qc, GATE_NAME)
    return qc, total_sign


def main():
    print("\n" + "=" * 96)
    print("  task31b_iyyi_convergence.py -- BLOCKING: does IYYI converge to analytic or literal?")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)

    with open(os.path.join(os.path.dirname(__file__), "task30b_pec_calibration_results.json")) as f:
        learned = json.load(f)
    gpi_bins = learned["forte-1"]["gpi"]
    p1_gpi = float(np.mean([v["p"] for v in gpi_bins.values()]))
    p1_gpi2 = p1_gpi  # Task 30B's disclosed fallback, unchanged (Task 31A only bounded it, didn't replace it)
    print(f"  using p2(zz)={P2_ZZ_CORRECTED} (Task 31A consensus), p1(gpi)={p1_gpi:.6f}, "
          f"p1(gpi2)={p1_gpi2:.6f} ({GPI2_FALLBACK_NOTE})")

    base = optimized_native_circuit(fixed_solutions[SLOT]["angles"], GATE_NAME)
    combined = effrag_mod.combined_basis_label(GROUP)
    basis_qc = native_basis_change(combined, GATE_NAME)
    full_base = base.compose(basis_qc)

    gamma_gate2 = gamma_factor(P2_ZZ_CORRECTED, 2)
    gamma_gate1 = gamma_factor(p1_gpi, 1)
    n2q = sum(1 for instr in full_base.data if instr.operation.name == GATE_NAME)
    n1q = sum(1 for instr in full_base.data if instr.operation.name in ("gpi", "gpi2"))
    gamma_total = (gamma_gate2 ** n2q) * (gamma_gate1 ** n1q)
    print(f"  circuit: {n2q} zz gates, {n1q} 1q gates, gamma_total={gamma_total:.4f}")

    if os.path.exists(CKPT_PATH):
        print("  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"  connected, backend={backend.name}")

        rng = np.random.default_rng(stable_seed("task31b", SLOT, TARGET_LABEL))
        circuits, signs = [], []
        for draw in range(N_MC_MAX):
            twirled, sign = sample_twirled_circuit(full_base, P2_ZZ_CORRECTED, p1_gpi, p1_gpi2, rng)
            twirled.measure_all()
            circuits.append(twirled)
            signs.append(sign)
        print(f"  submitting {len(circuits)} twirled circuits (forte-1 only, aria-1 retired)...")
        job = submit_job(circuits, backend, "forte-1", shots=SHOTS)
        print(f"    job_id={job.job_id()}")
        counts = get_counts_list(job)
        ck = {"signs": signs, "counts": counts, "gamma_total": gamma_total}
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # -- analytic PEC value for this exact (slot, label), using the SAME corrected p2 --
    A, B = analytic_A_and_B(fixed_solutions[SLOT]["angles"], GATE_NAME, P2_ZZ_CORRECTED, gpi_bins, gpi_bins,
                             [TARGET_LABEL])
    analytic_val = B[TARGET_LABEL] / A[TARGET_LABEL] if abs(A[TARGET_LABEL]) > 1e-6 else None
    print(f"\n  analytic ratio B/{TARGET_LABEL}/A = {analytic_val}")

    # -- get the RAW (uncorrected) measured value for reference too --
    signs = ck["signs"]
    counts_list = ck["counts"]
    gamma_total = ck["gamma_total"]

    print(f"\n  -- N_MC convergence, cumulative over the SAME 128 real draws --")
    convergence = {}
    for n_mc in N_MC_LEVELS:
        vals_signed = []
        vals_raw = []
        rng = np.random.default_rng(stable_seed("task31b_combine", n_mc))
        for i in range(n_mc):
            resampled = bootstrap_counts(counts_list[i], SHOTS, rng)
            m = expectation_from_counts(resampled, TARGET_LABEL)
            vals_signed.append(signs[i] * m)
            vals_raw.append(m)
        literal_est = gamma_total * float(np.mean(vals_signed))
        literal_est_clamped = max(-1.0, min(1.0, literal_est))
        raw_est = float(np.mean(vals_raw))
        convergence[n_mc] = {"literal": literal_est, "literal_clamped": literal_est_clamped, "raw_mean": raw_est}
        diff_from_analytic = abs(literal_est_clamped - analytic_val) if analytic_val is not None else None
        print(f"    N_MC={n_mc:>4}: literal_PEC={literal_est:+.4f} (clamped {literal_est_clamped:+.4f})  "
              f"raw(unsigned mean)={raw_est:+.4f}  |diff from analytic|={diff_from_analytic:.4f}" if diff_from_analytic is not None
              else f"    N_MC={n_mc:>4}: literal_PEC={literal_est:+.4f}")

    print(f"\n" + "=" * 96)
    print(f"  DIAGNOSIS:")
    diffs = [abs(convergence[n]["literal_clamped"] - analytic_val) for n in N_MC_LEVELS]
    trend = "CONVERGING toward analytic" if diffs[-1] < diffs[0] * 0.5 else (
        "STAYING FAR from analytic (bimodal / real discrepancy)" if diffs[-1] > diffs[0] * 0.8 else "unclear trend")
    print(f"    |diff| trajectory across N_MC={N_MC_LEVELS}: {[round(d,4) for d in diffs]}")
    print(f"    -> {trend}")
    passes_tolerance = diffs[-1] < 0.02
    print(f"    ACCEPTANCE (|Delta P|<0.02 at N_MC=128, energy-weighted label, rank 10/37): "
          f"{'PASS' if passes_tolerance else 'FAIL'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({
            "analytic_val": analytic_val, "convergence": convergence, "diffs": diffs,
            "trend": trend, "passes_tolerance": bool(passes_tolerance),
        }, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
