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
    # A SEPARATE, MANDATORY check on top of the step-ratio test, added after a real false
    # positive: a MONOTONICALLY changing tail (e.g. linearly increasing/decreasing) can have
    # consecutive-step RATIOS that drift toward 1 purely because the values themselves are
    # growing (or shrinking) roughly additively -- (a+(n+1)d)/(a+nd) -> 1 as n grows, for ANY
    # constant step d, convergent or not. Caught on iteration 13's ZNE range sweep: errors
    # [15.967, 22.809, 29.893, 37.847, 46.403] are STRICTLY, MONOTONICALLY INCREASING (getting
    # WORSE at every step) yet the last-2 step-ratios (1.27, 1.23) were both under the 1.5x
    # threshold, so the ratio-only check called it a "PASS" -- clearly wrong, this is a
    # diverging/drifting sequence, not a plateau. A genuine plateau's tail should NOT be
    # perfectly monotonic across n_flat+1 consecutive points (real measurement noise breaks
    # monotonicity; a systematic, still-changing trend preserves it) -- checked here and
    # required, not optional, for has_floor to be True.
    tail_window = errors[-(n_flat + 1):] if len(errors) >= n_flat + 1 else errors
    strictly_increasing = all(b > a for a, b in zip(tail_window[:-1], tail_window[1:]))
    strictly_decreasing = all(b < a for a, b in zip(tail_window[:-1], tail_window[1:]))
    tail_is_monotonic = (strictly_increasing or strictly_decreasing) and len(tail_window) >= 3

    has_floor = (len(step_ratios) >= n_flat and
                 all(r < plateau_ratio_threshold for r in step_ratios[-n_flat:]) and
                 not tail_is_monotonic)
    tail_ratio = max(step_ratios[-n_flat:]) if len(step_ratios) >= n_flat else max(step_ratios)

    no_floor_at_all = (overall_ratio > floor_ratio_threshold or tail_is_monotonic) and not has_floor

    if no_floor_at_all:
        if tail_is_monotonic and overall_ratio <= floor_ratio_threshold:
            verdict = (f"DISQUALIFIED -- the tail ({[round(v, 3) for v in tail_window]}) is strictly "
                       f"{'increasing' if strictly_increasing else 'decreasing'} across all "
                       f"{len(tail_window)} most-aggressive points; small step-ratios alone "
                       "(consecutive-ratio -> 1 as values grow/shrink roughly additively) do NOT "
                       "make this a plateau -- it is still drifting, just slowly. See "
                       "RESEARCH_LEDGER.md iteration 13 for the case that motivated this check.")
        else:
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
        "tail_is_monotonic": bool(tail_is_monotonic),
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

    # iteration 13's real false positive: ZNE-linear extrapolated error monotonically
    # INCREASING as the noise-scale range widens (15.967 -> 22.809 -> 29.893 -> 37.847 ->
    # 46.403 kcal/mol) -- an earlier version of this function called this a "PASS" because
    # consecutive-step ratios (1.27, 1.23) were both under the 1.5x threshold, purely an
    # artifact of the values growing roughly additively. Must be DISQUALIFIED, not passed.
    range_sizes = [3, 4, 5, 6, 7]
    zne_linear_errs = [15.967, 22.809, 29.893, 37.847, 46.403]
    result3 = floor_test(range_sizes, zne_linear_errs)
    assert result3["disqualified"], f"floor_test failed to catch the monotonically-diverging ZNE-linear range sweep: {result3}"
    assert result3["tail_is_monotonic"], "expected tail_is_monotonic=True for a strictly increasing sequence"

    return True


if __name__ == "__main__":
    ok = _self_test()
    print("qforge.floor_test self-test:", "PASS" if ok else "FAIL")
