#!/usr/bin/env python3
"""
task37b_h4_noise_model.py -- iteration 37, Task B. Defines the reduced,
H4-task-specific noise model theta = (p_ZZ, p_GPi2, delta_ZZ, delta_GPi2,
p_readout) -- 5 parameters, per the proposal's own "start small, don't
overfit" instruction, NOT the 9-parameter task31h robustness-envelope
model and NOT full GST. GPi is deliberately EXCLUDED from theta and held
FIXED at its own real calibration -- Task 37 Phase 0 already established
it's tightly, genuinely known (std=0.000011, 14x tighter than the old
assumed range), so there's nothing for the H4 data to usefully re-learn
there; spending an identifiability degree of freedom on an
already-solved parameter would only make the OTHER four harder to pin
down.

EVERY prior/bound below is sourced from real data where it exists, and
explicitly flagged as weak/uninformative where it does not -- this
module does not fabricate confidence into anything, matching this
project's own established discipline.

  p_ZZ:      REAL, tight (Task 31A, 9/11 real positions):
             Gaussian(mean=0.014593, std=0.000124). A real, informative
             prior -- the joint fit should stay near here unless H4 data
             gives strong evidence otherwise (checked via the z-score of
             the fitted value, not silently trusted).
  p_GPi2:    NO real posterior exists -- Task 31A's own recalibration
             attempt is flagged "physical": false with nonsensical
             values (negative "probabilities", std up to 0.19). Given a
             WEAK (uninformative-within-bound) prior, bounded to
             [0, 0.05] -- NOT the historical |p_gpi2|<0.21 statistical
             bound. That bound was never a realistic guess at the true
             value; it was the residual width of a failed measurement.
             The [0, 0.05] cap here is a real, checked constraint: Task
             37 Phase 1B directly verified that literal-twirling PEC's
             sampling overhead (gamma) becomes astronomically infeasible
             (1e5-1e26 across this circuit's real gate counts) well
             before p_gpi2 reaches 0.05-0.08 -- so candidates above this
             range could never be tested with any practical shot budget
             regardless of what the true value is, making them
             operationally meaningless to include as "identifiable."
  delta_ZZ, delta_GPi2: coherent per-gate-TYPE angle bias (matching
             task31h_robustness_envelope's own established convention:
             one systematic bias per gate type, not independent per
             instance -- that was tried once, iteration 31, and produced
             an unrealistic Q95~827 kcal/mol until fixed). NO real
             measurement exists for either. Weak Gaussian(0, 0.01) prior
             -- wide enough to not presuppose an answer, narrow enough
             to stay in a physically plausible small-coherent-error
             regime (thousands of gate instances would make even a
             0.01-radian-scale systematic bias clearly visible in H4
             data if present, per the same N*delta_p accumulation logic
             this project has used since iteration 31).
  p_readout: NO real dedicated calibration exists in this project (
             checked directly: no readout-specific calibration result
             file found anywhere in vqe/, confirmed before writing this
             module, not assumed). Kept as task31h's own original
             Uniform(0, 0.02) bound, disclosed as uninformative.

Run (self-test only -- prints the model definition and one sanity draw):
    python vqe/task37b_h4_noise_model.py
"""
import numpy as np

PARAM_NAMES = ["p_zz", "p_gpi2", "delta_zz", "delta_gpi2", "p_readout"]
N_PARAMS = len(PARAM_NAMES)

# real, fixed (NOT part of theta -- Task 37 Phase 0 established these don't need re-learning)
GPI_REAL_MEAN = 0.000119   # Task 30B
GPI_REAL_STD = 0.000011    # Task 30B

# theta priors -- (mean, std) for a Gaussian penalty, (lo, hi) hard bound, and a disclosed
# confidence tag so downstream code (and future readers) can never mistake a weak prior for a
# real measurement
PRIORS = {
    "p_zz":      {"mean": 0.014593, "std": 0.000124, "bounds": (1e-6, 0.05),  "source": "REAL (Task 31A)"},
    "p_gpi2":    {"mean": 0.001,    "std": 0.02,      "bounds": (0.0, 0.05),  "source": "WEAK (no real posterior -- Task 31A recalibration failed; bound set by PEC feasibility, Task 37 Phase 1B)"},
    "delta_zz":  {"mean": 0.0,      "std": 0.01,      "bounds": (-0.05, 0.05), "source": "WEAK (no real measurement)"},
    "delta_gpi2": {"mean": 0.0,     "std": 0.01,      "bounds": (-0.05, 0.05), "source": "WEAK (no real measurement)"},
    "p_readout": {"mean": 0.005,    "std": 0.01,      "bounds": (0.0, 0.02),  "source": "WEAK (no dedicated real calibration exists in this project, confirmed by search before writing this)"},
}


def theta_to_dict(theta):
    return {name: float(theta[i]) for i, name in enumerate(PARAM_NAMES)}


def prior_neg_log_penalty(theta, weight_by_param=None):
    """Gaussian penalty per parameter, weighted by 1/std^2 (a real
    calibration-informed prior for p_zz, deliberately weak/uninformative
    penalties elsewhere -- returns per-parameter penalty terms, NOT
    summed, so callers can inspect which parameter is doing how much
    work)."""
    weight_by_param = weight_by_param or {}
    out = {}
    for i, name in enumerate(PARAM_NAMES):
        pr = PRIORS[name]
        w = weight_by_param.get(name, 1.0)
        out[name] = w * ((theta[i] - pr["mean"]) / pr["std"]) ** 2
    return out


def clip_to_bounds(theta):
    out = np.array(theta, dtype=float)
    for i, name in enumerate(PARAM_NAMES):
        lo, hi = PRIORS[name]["bounds"]
        out[i] = min(hi, max(lo, out[i]))
    return out


def sample_prior(rng):
    """One real, honest prior draw -- for pilot/diagnostic use, NEVER
    for choosing a point estimate to report."""
    theta = np.zeros(N_PARAMS)
    for i, name in enumerate(PARAM_NAMES):
        pr = PRIORS[name]
        lo, hi = pr["bounds"]
        theta[i] = min(hi, max(lo, rng.normal(pr["mean"], pr["std"])))
    return theta


def _self_test():
    print("\n" + "=" * 96)
    print("  task37b_h4_noise_model.py -- 5-parameter H4-native noise model definition")
    print("=" * 96)
    print(f"\n  theta = {PARAM_NAMES}  ({N_PARAMS} parameters, GPi held fixed at "
          f"{GPI_REAL_MEAN} +/- {GPI_REAL_STD}, not part of theta)")
    print(f"\n  {'param':<12}{'mean':>10}{'std':>10}{'bounds':>18}   source")
    for name in PARAM_NAMES:
        pr = PRIORS[name]
        print(f"  {name:<12}{pr['mean']:>10.5f}{pr['std']:>10.5f}   [{pr['bounds'][0]:.4f},{pr['bounds'][1]:.4f}]   {pr['source']}")

    rng = np.random.default_rng(0)
    theta0 = sample_prior(rng)
    print(f"\n  one sample draw from the prior: {theta_to_dict(theta0)}")
    penalties = prior_neg_log_penalty(theta0)
    print(f"  prior penalty per param at that draw: { {k: round(v,3) for k,v in penalties.items()} }")


if __name__ == "__main__":
    _self_test()
