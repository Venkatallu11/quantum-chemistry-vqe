#!/usr/bin/env python3
"""
task6_error_budget.py — iteration 26, Task 6. THE ERROR BUDGET, for
every configuration this iteration touched: exact classical reference,
method error, shot error, hardware bias, mitigation bias, drift, TOTAL.
Chemical accuracy PASSES only when |total| < 1.0 kcal/mol AND the
drift-aware uncertainty interval also supports the claim.
============================================================================
Every number here is either cited from an EARLIER real measurement in
this ledger (method error, hardware bias, mitigation bias) or from THIS
session's own Tasks 1/A (shot error, drift) -- nothing is re-derived or
guessed. Where a component was measured for one representative circuit
(shot error: fixed ansatz; drift: Z2-tapered raw) rather than
independently per scheme, that is stated explicitly, not silently
generalized.

Run:
    python vqe/task6_error_budget.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task6_error_budget_results.json")
CHEM_ACC = 1.0

# -- components established this session / cited from this ledger --
METHOD_ERROR = 1.06e-11 * 627.5094740631  # essentially 0; setup_fragment(strict=True) verifies Schmidt rank <= K=6 exactly
SHOT_ERROR_REAL = 1.845     # Task 1's REAL validation at 300,000 shots/circuit (NOT the optimistic 0.087 local prediction)
SHOT_ERROR_REAL_STD = 0.198
DRIFT = {"aria-1": 4.01, "forte-1": 2.31}   # Task A, iteration 25, drift-std across 8 independent real submissions

# hardware bias (raw real - ideal) and mitigation bias (mitigated real - ideal), per configuration
CONFIGS = {
    "fixed (11 CX)": {
        "raw": {"aria-1": 33.09, "forte-1": 42.59}, "mitigated": {"aria-1": 29.55, "forte-1": 30.29},
    },
    "ADAPT (8.53 CX mean)": {
        "raw": {"aria-1": 62.18, "forte-1": 80.45}, "mitigated": {"aria-1": 25.88, "forte-1": 34.80},
    },
    "variational (4.3 CX mean)": {
        "raw": {"aria-1": 62.97, "forte-1": 79.25}, "mitigated": {"aria-1": 27.45, "forte-1": 36.38},
    },
    "full stack (ADAPT, 21 circ)": {
        "raw": {"aria-1": 66.41, "forte-1": 70.29}, "mitigated": {"aria-1": 27.74, "forte-1": 29.52},
    },
}
IDEAL_BASELINE = 0.0  # exact classical reference IS zero by construction (err_vs_exact)


def main():
    print("\n" + "=" * 96)
    print("  task6_error_budget.py -- the full error budget, every configuration")
    print("=" * 96)
    print(f"  method error (K=6 EF vs exact, verified no truncation): {METHOD_ERROR:.2e} kcal/mol")
    print(f"  shot error (REAL, 300k shots/circuit, ideal control): {SHOT_ERROR_REAL:.3f}+/-{SHOT_ERROR_REAL_STD:.3f} kcal/mol")
    print(f"    (measured on the fixed-ansatz circuit only -- generalized here as a platform floor, not")
    print(f"     independently re-measured per configuration; a real, disclosed simplification)")
    print(f"  drift-std (REAL, 8 independent submissions, Z2-tapered raw circuit): "
          f"aria-1={DRIFT['aria-1']}, forte-1={DRIFT['forte-1']} kcal/mol")
    print(f"    (measured on the Z2-tapered circuit only -- generalized here as a platform floor, not")
    print(f"     independently re-measured per configuration; a real, disclosed simplification)")

    budgets = {}
    for cfg_name, c in CONFIGS.items():
        budgets[cfg_name] = {}
        print(f"\n  ===== {cfg_name} =====")
        for model in ["aria-1", "forte-1"]:
            raw_bias = c["raw"][model] - IDEAL_BASELINE
            mit_bias = c["mitigated"][model] - IDEAL_BASELINE
            total = mit_bias  # TOTAL = mitigated hardware - reference (reference=exact=0 in err_vs_exact terms)
            combined_uncertainty = (SHOT_ERROR_REAL ** 2 + DRIFT[model] ** 2) ** 0.5
            passes_central = abs(total) < CHEM_ACC
            passes_with_uncertainty = passes_central and (abs(total) + combined_uncertainty < CHEM_ACC * 5 or combined_uncertainty < CHEM_ACC)
            # honest PASS rule: central value below 1.0 AND the uncertainty interval's UPPER bound also plausibly near/under it
            # (explicit, conservative: requires total+combined_uncertainty to at least be < 2x chemical accuracy to call it "close";
            # true PASS requires the interval to support the claim, which basically never happens at these error scales)
            true_pass = (abs(total) < CHEM_ACC) and (abs(total) + combined_uncertainty < CHEM_ACC)
            row = {
                "method_error": METHOD_ERROR, "shot_error": SHOT_ERROR_REAL,
                "hardware_bias_raw": raw_bias, "mitigation_bias": mit_bias,
                "drift_std": DRIFT[model], "combined_uncertainty": combined_uncertainty,
                "total": total, "passes_central_value": bool(passes_central),
                "passes_with_uncertainty": bool(true_pass),
            }
            budgets[cfg_name][model] = row
            verdict = "PASS" if true_pass else ("CENTRAL VALUE UNDER 1.0 BUT UNCERTAINTY TOO WIDE -- NOT A PASS" if passes_central else "FAIL")
            print(f"    {model}: method={METHOD_ERROR:.1e}  shot={SHOT_ERROR_REAL:.2f}  "
                  f"hw_bias(raw)={raw_bias:.2f}  mitigation_bias={mit_bias:.2f}  drift_std={DRIFT[model]:.2f}  "
                  f"combined_uncertainty=+/-{combined_uncertainty:.2f}")
            print(f"      TOTAL (mitigated - reference) = {total:.2f} kcal/mol  |total|+uncertainty = "
                  f"{abs(total)+combined_uncertainty:.2f}  -- {verdict}")

    print(f"\n  -- SUMMARY: does ANY configuration pass chemical accuracy, central value AND uncertainty? --")
    any_pass = False
    for cfg_name, models in budgets.items():
        for model, row in models.items():
            if row["passes_with_uncertainty"]:
                any_pass = True
                print(f"    PASS: {cfg_name} / {model}")
    if not any_pass:
        print(f"    NO configuration passes both the central-value bar (<1.0 kcal/mol) AND the drift-aware")
        print(f"    uncertainty interval. The closest central values (fixed forte-1=30.29, full-stack")
        print(f"    forte-1=29.52) are ~29-30x the chemical accuracy bar even before adding uncertainty --")
        print(f"    this project has NOT reached chemical accuracy on real IonQ hardware at K=6, and stating")
        print(f"    that plainly is the entire point of this task, per explicit instruction: 'a central")
        print(f"    estimate of 0.8 with a +/-4 bar is NOT a pass. Say so explicitly.'")

    results = {"method_error": METHOD_ERROR, "shot_error_real": SHOT_ERROR_REAL, "drift": DRIFT,
               "budgets": budgets, "any_configuration_passes": any_pass}
    with open(RESULTS_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved -> {RESULTS_PATH}\n")
    return results


if __name__ == "__main__":
    main()
