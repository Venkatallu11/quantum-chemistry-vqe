#!/usr/bin/env python3
"""
qforge.floor_test — the mandatory floor test as a reusable function, not
re-implemented ad hoc per iteration. Every free parameter a mitigation
method has (training radius, training set size, number of copies,
regularization strength, fit granularity, ...) must be swept through this
before a result is recorded in RESEARCH_LEDGER.md.

WHAT IT CATCHES: iteration 2 (loop_local_perturbation_cdr.py) reported
0.0636 kcal/mol as "target reached" before this function existed as
reusable code -- a by-hand sweep of PERTURB_RADIUS afterward found error
falling monotonically with NO plateau (0.60 rad -> 0.439, 0.30 -> 0.132,
0.15 -> 0.028, 0.05 -> 0.008, 0.01 -> 0.007 kcal/mol, see
RESEARCH_LEDGER.md iteration 2), meaning the method was interpolating
toward the classically-known target rather than measuring device noise:
as the radius shrinks, the "training circuit" converges to the target
circuit itself. That diagnosis is encoded here as floor_test's own
regression test (`_self_test()`) so the exact historical failure this
function exists to catch stays caught.

A method PASSES the floor test when error stops falling (plateaus) as
the free parameter is pushed toward its most-aggressive setting -- that
plateau IS the floor: real, device-measured noise residual that no
amount of further parameter tuning removes. A method with NO floor
across its full tested range is not proven safe by that alone (the sweep
might just not have gone far enough), but it must be reported as
INCONCLUSIVE, never as a passing result -- silently treating
"not-yet-disqualified" as "passing" is exactly how iteration 2 happened.
"""
import numpy as np


def floor_test(param_values, errors, floor_ratio_threshold=3.0, plateau_ratio_threshold=1.5,
                min_consecutive_flat_steps=2):
    """
    param_values: the free parameter's swept values, ordered from LEAST
      aggressive to MOST aggressive (aggressive = the setting that would
      make the classical-cheat failure mode, if present, most extreme --
      e.g. increasing training-set size, or SHRINKING a training-radius
      parameter toward the target).
    errors: the method's error (kcal/mol, or any consistent unit) at each
      param_values entry, same order, same length.

    A genuine floor means CONSECUTIVE steps near the aggressive end stop
    changing much -- checked as `min_consecutive_flat_steps` IN A ROW
    consecutive-step ratios all falling below plateau_ratio_threshold, not
    just the overall min/max ratio of a trailing window. That distinction
    is not cosmetic: iteration 2's own disqualifying sweep (0.439, 0.132,
    0.028, 0.008, 0.007 kcal/mol) has a single small LAST-step ratio
    (0.008/0.007=1.14x) purely because both values are already tiny, while
    every OTHER consecutive step in the same sweep is still a 3-5x jump --
    a single quiet step is not a floor, it is a sequence still heading to
    zero that happened to sample two nearby points. An earlier version of
    this function used a trailing-window min/max ratio and was WRONG on
    exactly this data (called it a pass) -- caught by _self_test() below,
    which is why that test exists and must keep passing.

    Returns a dict with `disqualified` (bool) and `verdict` (str). Never
    returns a bare pass/fail without the diagnostic ratios that justify
    it -- every floor-test verdict recorded in the ledger must be
    traceable back to real numbers, not an assertion.
    """
    assert len(param_values) == len(errors) >= 2, "need at least 2 swept points to test for a floor"
    errors = [float(e) for e in errors]
    overall_ratio = (max(errors) / min(errors)) if min(errors) > 0 else float("inf")

    step_ratios = []
    for a, b in zip(errors[:-1], errors[1:]):
        step_ratios.append((max(a, b) / min(a, b)) if min(a, b) > 0 else float("inf"))

    n_flat = min_consecutive_flat_steps
    has_floor = (len(step_ratios) >= n_flat and
                 all(r < plateau_ratio_threshold for r in step_ratios[-n_flat:]))
    tail_ratio = max(step_ratios[-n_flat:]) if len(step_ratios) >= n_flat else max(step_ratios)

    no_floor_at_all = overall_ratio > floor_ratio_threshold and not has_floor

    if no_floor_at_all:
        verdict = ("DISQUALIFIED -- error falls with no plateau across the full sweep "
                    f"({overall_ratio:.1f}x overall; the last {n_flat} consecutive step-ratios "
                    f"are {[round(r, 2) for r in step_ratios[-n_flat:]]}, not all below "
                    f"{plateau_ratio_threshold}x). This is the signature of interpolating toward "
                    "a classically-known answer, not measuring real device noise -- see "
                    "RESEARCH_LEDGER.md iteration 2.")
        disqualified = True
    elif has_floor:
        verdict = (f"PASS -- the last {n_flat} consecutive step-ratios "
                   f"({[round(r, 2) for r in step_ratios[-n_flat:]]}) are all below "
                   f"{plateau_ratio_threshold}x. This plateau is the method's real floor.")
        disqualified = False
    else:
        verdict = (f"INCONCLUSIVE -- {overall_ratio:.1f}x overall change, step-ratios "
                   f"{[round(r, 2) for r in step_ratios]} have not yet plateaued for "
                   f"{n_flat} consecutive steps, but overall change is below the "
                   f"{floor_ratio_threshold}x disqualification threshold. Extend the sweep "
                   "before recording a result; do not report this as a pass.")
        disqualified = False

    return {
        "param_values": list(param_values), "errors": errors,
        "overall_ratio": overall_ratio, "step_ratios": step_ratios, "tail_ratio": tail_ratio,
        "still_falling_at_most_aggressive_point": bool(errors[-1] < errors[-2]),
        "has_floor": bool(has_floor), "disqualified": bool(disqualified),
        "verdict": verdict,
    }


def _self_test():
    """Regression test: floor_test() must disqualify iteration 2's actual
    historical PERTURB_RADIUS sweep. Run automatically by
    qforge/tests/test_qforge.py; also runnable standalone."""
    # radius, ordered LEAST aggressive (large radius) -> MOST aggressive
    # (small radius, closest to "just re-evaluate the target exactly")
    radii = [0.60, 0.30, 0.15, 0.05, 0.01]
    errors_kcal = [0.439, 0.132, 0.028, 0.008, 0.007]
    result = floor_test(radii, errors_kcal)
    assert result["disqualified"], f"floor_test failed to catch iteration 2's known-bad sweep: {result}"

    # a method with a genuine floor: error drops then plateaus
    param = [1, 2, 4, 8, 16, 32]
    errs = [10.0, 5.0, 2.6, 2.1, 2.0, 2.0]
    result2 = floor_test(param, errs)
    assert result2["has_floor"] and not result2["disqualified"], f"floor_test wrongly flagged a real floor: {result2}"
    return True


if __name__ == "__main__":
    ok = _self_test()
    print("qforge.floor_test self-test:", "PASS" if ok else "FAIL")
