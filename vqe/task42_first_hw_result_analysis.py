#!/usr/bin/env python3
"""
task42_first_hw_result_analysis.py -- iteration 42. Sanity-check analysis
of Task 41's first real hardware result (100 shots, qpu.forte-enterprise-1,
u_0 slot, group=[XYYX,IYYI], ancilla-augmented native circuit).

ONE circuit at 100-500 shots cannot reconstruct the full H4 energy (that
needs all 21 slots x 13 groups combined) -- this is NOT an energy
estimate. What it CAN do, honestly, is answer "does the real hardware
behave like our own already-validated model expects for this exact
circuit": compare the real QPU's postselected <XYYX>, <IYYI> and ancilla
retained-fraction directly against the SAME slot/group/postselection
already collected on the real `ionq_simulator` (draw 0, Task 39C), which
itself already matches the noise model this project's correction is
built on. Large, structural disagreement would mean something is wrong
with the circuit/basis-change/ancilla construction on real hardware.
Close agreement is real, non-fabricated evidence the circuit runs as
intended.

Run:
    PYTHONHASHSEED=0 python vqe/task42_first_hw_result_analysis.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from ionq_simulator_binding_curve import expectation_from_counts

HW_RESULTS = os.path.join(os.path.dirname(__file__), "task41_first_hw_submission_results.json")
SIM_CKPT = os.path.join(os.path.dirname(__file__), "task39c_ancilla_real_submission.partial.json")
SLOT = "u_0"
GROUP = ("XYYX", "IYYI")


def postselect(counts, keep_zero=True):
    filtered = {}
    for bs, c in counts.items():
        anc, reg = bs[0], bs[1:]
        if (anc == "0") == keep_zero:
            filtered[reg] = filtered.get(reg, 0) + c
    return filtered


def main():
    print("\n" + "=" * 96)
    print("  task42_first_hw_result_analysis.py -- real HW vs real ionq_simulator, same circuit")
    print("=" * 96)

    with open(HW_RESULTS) as f:
        hw = json.load(f)
    hw_counts = hw["counts"]
    hw_total = sum(hw_counts.values())
    hw_ps = postselect(hw_counts, keep_zero=True)
    hw_retained = sum(hw_ps.values())
    print(f"\n  REAL HARDWARE (qpu.forte-enterprise-1, job {hw['job_id']}, {hw_total} shots):")
    print(f"    ancilla=0 retained: {hw_retained}/{hw_total} = {100*hw_retained/hw_total:.1f}%")
    for label in GROUP:
        e_ps = expectation_from_counts(hw_ps, label)
        print(f"    <{label}> postselected(anc=0) = {e_ps:+.4f}")

    print(f"\n  REAL ionq_simulator (Task 39C draw 0, SAME slot={SLOT}, SAME group={GROUP}):")
    with open(SIM_CKPT) as f:
        sim = json.load(f)
    for backend_name in ["aria-1", "forte-1"]:
        entry = sim["done"][f"{backend_name}|{SLOT}"]
        idx = entry["groups"].index(list(GROUP))
        sim_counts = entry["counts"][idx]
        sim_total = sum(sim_counts.values())
        sim_ps = postselect(sim_counts, keep_zero=True)
        sim_retained = sum(sim_ps.values())
        print(f"    [{backend_name}] {sim_total} shots, ancilla=0 retained: "
              f"{sim_retained}/{sim_total} = {100*sim_retained/sim_total:.1f}%")
        for label in GROUP:
            e_ps = expectation_from_counts(sim_ps, label)
            print(f"      <{label}> postselected(anc=0) = {e_ps:+.4f}")

    print(f"\n  -- READ --")
    print(f"  Only {hw_total} real HW shots (vs {sim_total} sim shots) -- HW numbers carry large")
    print(f"  statistical noise (~1/sqrt(N) ~ {1/hw_total**0.5:.2f} scale) and should NOT be expected to")
    print(f"  match precisely. What matters: same SIGN, same rough MAGNITUDE, and a retained-fraction")
    print(f"  in the same ballpark -- that is what confirms the circuit/basis-change/ancilla construction")
    print(f"  is behaving physically on real hardware, not producing garbage/degenerate output.")


if __name__ == "__main__":
    main()
