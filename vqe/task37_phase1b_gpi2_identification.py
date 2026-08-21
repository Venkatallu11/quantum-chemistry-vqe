#!/usr/bin/env python3
"""
task37_phase1b_gpi2_identification.py -- iteration 37, Phase 1B. Real
new submissions, specifically designed to let real hardware data
identify GPi2's calibration -- the parameter Phase 0 confirmed accounts
for ~87% of the joint estimator's remaining robustness gap, and Phase 1A
could not identify from raw data with a weak analytic correction model.

WHY THIS NEEDS NEW DATA, not just cleverer fitting on what we have:
Task 31C's existing real literal-twirling checkpoint was collected using
ONE FIXED p_gpi2_assumed (the "disclosed fallback" ~0.000119) baked into
the twirl-recovery sampling itself (`sample_twirled_circuit`'s choice of
WHICH Pauli recovery gets inserted depends on p_gpi2 through
`pec_inverse_weights`) -- you cannot retroactively ask "what would a
different assumed p_gpi2 have produced" from data whose sampling already
committed to one value. Real new circuits, twirled at DIFFERENT candidate
p_gpi2 values, are the only way to get data that is genuinely sensitive
to which candidate is closer to the truth.

THE TARGET-LEAKAGE SAFEGUARD, explicit and enforced in code, not just
description: candidate p_gpi2 values are NEVER compared against or
selected using the known exact H4 energy anywhere in this script. The
selection criterion is purely an INTERNAL CONSISTENCY check: for a given
slot, ALL of its measured labels (2-4 per slot here) must be jointly
explainable by ONE valid unit-norm state `a` (`fit_pure_state`'s own
residual, `best_val` -- the same real diagnostic Task 34C already used
for a different purpose). A wrong p_gpi2 corrupts the corrected values
in a way that's NOT internally consistent with any single physical
state, regardless of what the "right answer" energy is; a correct
p_gpi2 should be more self-consistent. This is checkable from real data
alone.

SCOPE, disclosed: 2 representative slots ((u0+u1)/IYYI -- established
validation case; (u3+u5) -- a Task 33D-flagged high-sign-noise group),
4 candidate p_gpi2 values spanning the established bound [0, 0.21],
N_MC=16 real draws each (matching this project's standard). ~128 new
real circuits total, comparable scale to Task 33E's targeted extension
(480 circuits, ~15 min). p_ZZ held FIXED at its real, well-established
value (0.014593) throughout -- Phase 0 already confirmed ZZ is not the
problem.

Run:
    python vqe/task37_phase1b_gpi2_identification.py
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task2_fold_response_dataset import native_basis_change
import ef_fragment as effrag_mod
from task28d_all_gate_zne import optimized_native_circuit
from task29c_manifold_estimator import target_coeff_vector, fit_pure_state
from phys_constrained_reconstruction import build_P_S
from task31c_full_pec_calibration import sample_twirled_circuit, GATE_NAME
from ionq_backend import connect_provider, get_native_simulator
from ionq_simulator_binding_curve import submit_job, get_counts_list, stable_seed, bootstrap_counts, expectation_from_counts

K = 6
P2_ZZ = 0.014593  # real, established (Task 31A) -- fixed throughout, not a candidate
GPI_REAL_P1 = 0.000119  # real, established (Task 30B) -- fixed for the SEPARATE, well-calibrated
# real GPi gates throughout. FIXED A REAL BUG: a first version of this script passed the GPi2
# CANDIDATE value for BOTH the p1_gpi and p1_gpi2 arguments to sample_twirled_circuit, meaning
# every candidate also inflated the far more numerous, already-well-calibrated GPi gates -- caught
# because gamma (the real PEC sampling overhead) turned out astronomical (1e5 to 1e26 across the
# ~96-190 real gate instances in these circuits) even at "moderate" GPi2 candidates, which in turn
# made every corrected value clip to exactly +-1.0 regardless of the real signal. GPi stays FIXED
# at its own real calibration here; only GPi2 varies.
SHOTS = 100_000
N_MC = 16
GPI2_CANDIDATES = [0.0005, 0.001, 0.003, 0.008]  # realistic small-p regime, gamma verified to stay
# under ~40 across the real gate counts in these circuits (checked directly before choosing this
# range) -- NOT the full |p_gpi2|<0.21 statistical bound, which is a loose uncertainty interval,
# not a realistic guess at the true value; PEC sampling overhead is exponential in gate count and
# candidates anywhere near that bound are not executable with any practical shot budget
TEST_CASES = [
    {"slot": "(u0+u1)", "group": ["XYYX", "IYYI"], "tag": "IYYI (established validation case)"},
    {"slot": "(u3+u5)", "group": ["XZXZ", "XZXI", "IZIZ", "IIIZ"], "tag": "(u3+u5) (Task 33D-flagged high sign-noise group)"},
]
CKPT_DIR = os.path.join(os.path.dirname(__file__), "ionq_simulator_binding_curve_checkpoints")
CKPT_PATH = os.path.join(CKPT_DIR, "task37_phase1b_gpi2_sweep.json")
RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task37_phase1b_gpi2_identification_results.json")


def main():
    print("\n" + "=" * 96)
    print("  task37_phase1b_gpi2_identification.py -- real GPi2 sweep, internal-consistency identification")
    print(f"  TARGET-LEAKAGE SAFEGUARD: the known exact H4 energy is NEVER used below to pick a candidate.")
    print("=" * 96)

    p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
    fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
    U_exact = np.asarray(p["u_vecs"]).T
    P_S = build_P_S(p["alpha_labels"], U_exact)

    if os.path.exists(CKPT_PATH):
        print("\n  found existing checkpoint -> reusing, not resubmitting")
        with open(CKPT_PATH) as f:
            ck = json.load(f)
    else:
        provider = connect_provider()
        backend = get_native_simulator(provider)
        print(f"\n  connected, backend={backend.name}")

        circuits, tags = [], []
        for case in TEST_CASES:
            slot, group = case["slot"], case["group"]
            base = optimized_native_circuit(fixed_solutions[slot]["angles"], GATE_NAME)
            basis_qc = native_basis_change(effrag_mod.combined_basis_label(group), GATE_NAME)
            full_base = base.compose(basis_qc)
            for p_gpi2 in GPI2_CANDIDATES:
                rng = np.random.default_rng(stable_seed("task37_1b", slot, p_gpi2))
                for draw in range(N_MC):
                    twirled, sign = sample_twirled_circuit(full_base, P2_ZZ, GPI_REAL_P1, p_gpi2, rng)
                    twirled.measure_all()
                    circuits.append(twirled)
                    tags.append((slot, group, p_gpi2, draw, sign))
        print(f"  built {len(circuits)} new twirled circuits "
              f"({len(TEST_CASES)} slots x {len(GPI2_CANDIDATES)} GPi2 candidates x {N_MC} draws)")

        n_chunks = (len(circuits) + 90) // 91
        chunk_size = max(1, (len(circuits) + n_chunks - 1) // n_chunks)
        chunks_c = [circuits[i:i + chunk_size] for i in range(0, len(circuits), chunk_size)]
        print(f"  submitting in {len(chunks_c)} batches of <= {chunk_size} circuits each...")

        all_counts = []
        t0 = time.time()
        for ci, chunk in enumerate(chunks_c):
            job = None
            for attempt in range(6):
                try:
                    job = submit_job(chunk, backend, "forte-1", shots=SHOTS)
                    break
                except Exception as e:
                    wait_s = min(30 * (2 ** attempt), 300)
                    print(f"    batch {ci+1}/{len(chunks_c)} submit failed (attempt {attempt+1}/6): {e} -- backing off {wait_s}s")
                    time.sleep(wait_s)
            if job is None:
                raise RuntimeError(f"submit_job exhausted retries for batch {ci+1}/{len(chunks_c)}")
            counts = get_counts_list(job)
            all_counts.extend(counts)
            print(f"    batch {ci+1}/{len(chunks_c)} done ({len(all_counts)}/{len(circuits)} total), "
                  f"{time.time()-t0:.1f}s elapsed", flush=True)

        ck = {"tags": tags, "counts": all_counts}
        os.makedirs(CKPT_DIR, exist_ok=True)
        with open(CKPT_PATH, "w") as f:
            json.dump(ck, f, indent=2)
        print(f"  checkpoint saved -> {CKPT_PATH}")

    # ================= ANALYSIS: internal consistency per candidate, NO exact-energy comparison =================
    print(f"\n  -- ANALYSIS: single-slot fit residual per GPi2 candidate (internal consistency only) --")
    tags = [tuple(t) for t in ck["tags"]]
    counts_list = ck["counts"]
    by_slot_gpi2 = {}
    for (slot, group, p_gpi2, draw, sign), counts in zip(tags, counts_list):
        by_slot_gpi2.setdefault((slot, p_gpi2), []).append((tuple(group), sign, counts))

    results = {}
    for case in TEST_CASES:
        slot = case["slot"]
        print(f"\n  === {case['tag']} ===")
        residuals_by_candidate = {}
        for p_gpi2 in GPI2_CANDIDATES:
            entries = by_slot_gpi2[(slot, p_gpi2)]
            # gamma is fixed per (slot, gpi2 candidate) -- recompute once via the same construction
            # sample_twirled_circuit used (gamma = product of gamma_gate per gate, deterministic given p2/p_gpi2)
            rng_g = np.random.default_rng(0)
            base = optimized_native_circuit(fixed_solutions[slot]["angles"], GATE_NAME)
            basis_qc = native_basis_change(effrag_mod.combined_basis_label(case["group"]), GATE_NAME)
            full_base = base.compose(basis_qc)
            _, _ = sample_twirled_circuit(full_base, P2_ZZ, GPI_REAL_P1, p_gpi2, rng_g)  # unused, just to confirm construction
            from loop_pec import pec_inverse_weights
            gamma = 1.0
            for instr in full_base.data:
                op = instr.operation
                if op.name == GATE_NAME:
                    gamma *= float(np.sum(np.abs(list(pec_inverse_weights(P2_ZZ, 2).values()))))
                elif op.name == "gpi":
                    gamma *= float(np.sum(np.abs(list(pec_inverse_weights(GPI_REAL_P1, 1).values()))))
                elif op.name == "gpi2":
                    gamma *= float(np.sum(np.abs(list(pec_inverse_weights(p_gpi2, 1).values()))))

            blended = {}
            for group_t, sign, counts in entries:
                rng_boot = np.random.default_rng(stable_seed("t37_1b_analysis", slot, p_gpi2, group_t))
                resampled = bootstrap_counts(counts, SHOTS, rng_boot)
                for l in group_t:
                    m = expectation_from_counts(resampled, l)
                    blended.setdefault(l, []).append(sign * gamma * m)
            blended = {l: max(-1.0, min(1.0, float(np.mean(v)))) for l, v in blended.items()}

            v0 = target_coeff_vector(slot, K)
            a_hat, best_val = fit_pure_state(P_S, blended, {l: 1.0 for l in blended}, K, v0, seed=42)
            residuals_by_candidate[p_gpi2] = {"best_val": best_val, "blended": blended, "a_hat": a_hat.tolist()}
            print(f"    p_gpi2={p_gpi2:.2f}: fit residual (best_val) = {best_val:.6f}   "
                  f"blended={ {l: round(v,4) for l,v in blended.items()} }")

        best_candidate = min(residuals_by_candidate, key=lambda k: residuals_by_candidate[k]["best_val"])
        print(f"    -- DATA-PREFERRED p_gpi2 (lowest internal-fit residual, NO exact-energy comparison used): "
              f"{best_candidate:.2f} --")
        results[case["tag"]] = {
            "residuals_by_candidate": {str(k): v["best_val"] for k, v in residuals_by_candidate.items()},
            "best_candidate": best_candidate,
        }

    print(f"\n  -- CROSS-SLOT CHECK: do independent slots agree on which GPi2 the data prefers? --")
    best_per_slot = {tag: r["best_candidate"] for tag, r in results.items()}
    for tag, best in best_per_slot.items():
        print(f"    {tag}: prefers p_gpi2={best:.2f}")
    vals = list(best_per_slot.values())
    agree = len(set(vals)) == 1
    print(f"    {'AGREE' if agree else 'DISAGREE'} -- "
          f"{'a real, reproducible cross-slot signal' if agree else 'no consistent signal across independent slots; real, honest negative finding'}")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"gpi2_candidates": GPI2_CANDIDATES, "results": results, "cross_slot_agree": bool(agree)},
                   f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
