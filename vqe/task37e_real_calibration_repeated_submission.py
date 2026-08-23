#!/usr/bin/env python3
"""
task37e_real_calibration_repeated_submission.py -- iteration 37, Task E
real follow-up, repeated-draw version. The single-draw submission
(`task37e_real_calibration_submission.py`) found aria-1 delta_zz_hat=
-0.00114 (~2.3 SE from zero) but aria-1/forte-1 DISAGREED IN SIGN for
delta_zz -- exactly the regime this project's own established finding
(quantified from 8 real submissions: real drift-std aria-1=4.01,
forte-1=2.31 kcal/mol-equivalent on OTHER circuits) warns is subject to
real submission-to-submission variation, not just shot noise. This is
NOT bootstrap resampling of one submission's counts (that only captures
shot noise, not true submission drift) -- these are N_DRAWS=8 INDEPENDENT
REAL API submissions per backend, matching the exact protocol that
established the cross-submission-drift finding in the first place, so
the resulting mean/std is directly comparable to it.

Run:
    PYTHONHASHSEED=0 python vqe/task37e_real_calibration_repeated_submission.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from task37e_real_calibration_submission import build_zz_circuit, build_gpi2_circuit, SLOPE_ZZ, SLOPE_GPI2
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, expectation_from_counts

SHOTS = 100_000
BACKENDS = ["ideal", "aria-1", "forte-1"]
N_DRAWS = 8
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37e_real_calibration_repeated_submission_results.json")


def main():
    print("\n" + "=" * 96)
    print(f"  task37e_real_calibration_repeated_submission.py -- {N_DRAWS} independent REAL submissions/backend")
    print("=" * 96)

    zz_circuit = build_zz_circuit()
    gpi2_circuit = build_gpi2_circuit()

    provider = connect_provider()
    backend = get_native_simulator(provider)
    print(f"  connected, backend={backend.name}")

    m_xi_by_backend = {b: [] for b in BACKENDS}
    m_z_by_backend = {b: [] for b in BACKENDS}
    raw_counts_log = []

    for draw in range(N_DRAWS):
        for backend_name in BACKENDS:
            c1 = zz_circuit.copy(); c1.measure_all()
            c2 = gpi2_circuit.copy(); c2.measure_all()
            job = submit_job([c1, c2], backend, backend_name, shots=SHOTS)
            counts = get_counts_list(job)
            m_xi = expectation_from_counts(counts[0], "XI")
            m_z = expectation_from_counts(counts[1], "Z")
            m_xi_by_backend[backend_name].append(m_xi)
            m_z_by_backend[backend_name].append(m_z)
            raw_counts_log.append({"draw": draw, "backend": backend_name, "counts_zz": counts[0], "counts_gpi2": counts[1]})
        print(f"    draw {draw+1}/{N_DRAWS} done ({len(BACKENDS)} backends x 2 circuits)")

    with open(RESULTS_PATH.replace(".json", "_raw.json"), "w") as f:
        json.dump(raw_counts_log, f, indent=2)

    print(f"\n  -- RESULTS: mean +/- std ACROSS {N_DRAWS} INDEPENDENT REAL SUBMISSIONS (real drift, not shot noise) --")
    results = {}
    for backend_name in BACKENDS:
        xi_vals = np.array(m_xi_by_backend[backend_name])
        z_vals = np.array(m_z_by_backend[backend_name])
        xi_mean, xi_std = float(xi_vals.mean()), float(xi_vals.std(ddof=1))
        z_mean, z_std = float(z_vals.mean()), float(z_vals.std(ddof=1))
        xi_sem = xi_std / np.sqrt(N_DRAWS)
        z_sem = z_std / np.sqrt(N_DRAWS)
        delta_zz_hat = xi_mean / SLOPE_ZZ
        delta_gpi2_hat = z_mean / SLOPE_GPI2
        se_delta_zz = xi_sem / abs(SLOPE_ZZ)
        se_delta_gpi2 = z_sem / abs(SLOPE_GPI2)
        sig_zz = abs(delta_zz_hat) > 2 * se_delta_zz
        sig_gpi2 = abs(delta_gpi2_hat) > 2 * se_delta_gpi2
        print(f"\n  === {backend_name} ({N_DRAWS} draws) ===")
        print(f"    <XI> across draws: {np.array2string(xi_vals, precision=6)}")
        print(f"    <XI> mean={xi_mean:+.6f}  std_across_draws={xi_std:.6f}  SEM={xi_sem:.6f}")
        print(f"    <Z>  across draws: {np.array2string(z_vals, precision=6)}")
        print(f"    <Z>  mean={z_mean:+.6f}  std_across_draws={z_std:.6f}  SEM={z_sem:.6f}")
        print(f"    delta_zz_hat   = {delta_zz_hat:+.6f}  (SEM~{se_delta_zz:.2e})  "
              f"{'SIGNIFICANT (>2 SEM)' if sig_zz else 'not significant vs real cross-draw variation'}")
        print(f"    delta_gpi2_hat = {delta_gpi2_hat:+.6f}  (SEM~{se_delta_gpi2:.2e})  "
              f"{'SIGNIFICANT (>2 SEM)' if sig_gpi2 else 'not significant vs real cross-draw variation'}")
        results[backend_name] = {
            "xi_vals": xi_vals.tolist(), "z_vals": z_vals.tolist(),
            "xi_mean": xi_mean, "xi_std_across_draws": xi_std, "xi_sem": xi_sem,
            "z_mean": z_mean, "z_std_across_draws": z_std, "z_sem": z_sem,
            "delta_zz_hat": delta_zz_hat, "delta_gpi2_hat": delta_gpi2_hat,
            "se_delta_zz": se_delta_zz, "se_delta_gpi2": se_delta_gpi2,
            "significant_zz": bool(sig_zz), "significant_gpi2": bool(sig_gpi2),
        }

    ideal_ok = abs(results["ideal"]["delta_zz_hat"]) < 2 * results["ideal"]["se_delta_zz"] + 1e-3 and \
               abs(results["ideal"]["delta_gpi2_hat"]) < 2 * results["ideal"]["se_delta_gpi2"] + 1e-3
    print(f"\n  -- SANITY CHECK: ideal backend should read ~0 for both -- "
          f"{'PASS' if ideal_ok else 'FAIL -- do not trust aria-1/forte-1 below'}")

    print(f"\n  -- CROSS-BACKEND COMPARISON (aria-1 vs forte-1 delta_zz sign disagreement from the single-draw run) --")
    az = results["aria-1"]["delta_zz_hat"]
    fz = results["forte-1"]["delta_zz_hat"]
    print(f"    aria-1 delta_zz_hat={az:+.6f}   forte-1 delta_zz_hat={fz:+.6f}   "
          f"{'SAME SIGN, now consistent' if az*fz > 0 else 'STILL OPPOSITE SIGN -- real backend difference or both consistent with zero'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"ideal_ok": bool(ideal_ok), "n_draws": N_DRAWS, "results": results,
                    "slope_zz": SLOPE_ZZ, "slope_gpi2": SLOPE_GPI2, "shots": SHOTS}, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
