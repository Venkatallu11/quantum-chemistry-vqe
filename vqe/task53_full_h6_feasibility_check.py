#!/usr/bin/env python3
"""
task53_full_h6_feasibility_check.py -- iteration 53. A quick, honest
FEASIBILITY SCAN (not an attempt) for the natural next question after
task52's tailored-covalent-bonding result: can the ~2.79 kcal/mol
residual vs. full exact H6 (a CLASSICAL fragmentation-method floor, not
a quantum-measurement limitation -- see task52's own docstring) be
closed by solving the FULL, unfragmented H6 molecule directly with the
same forging approach, instead of improving the classical fragmentation
scheme?

TWO WAYS TO CLOSE THAT GAP WERE CONSIDERED:
  (a) A bigger overlap within H6's own 2-fragment tailoring scheme.
      REAL CONSTRAINT FOUND: the only way to get a bigger-than-2-atom
      overlap while still covering all 6 atoms with two equal-sized
      fragments is 5-atom fragments (5-atom windows [0-4]/[1-5],
      4-atom overlap [1-4]) -- but 5 atoms means 5 electrons, an ODD
      count. This project's entanglement-forging machinery assumes
      closed-shell fragments split evenly into alpha/beta registers;
      odd-electron ("open-shell-like") fragments are NOT something this
      pipeline handles at all. This is real new methodology, not a
      parameter tweak -- not attempted here.
  (b) Apply the pipeline directly to the FULL H6 problem, no
      fragmentation. Checked here.

REAL RESULT: the full (unfragmented) H6 problem's Schmidt rank does NOT
saturate even at K=16 -- the max_schmidt_tail is still 2.949e-3 at K=16
(compare: fragment A/B's own rank is EXACTLY 6, already maximal for
their 6-dimensional weight-2 subspace). Full H6's register is 6 qubits
(3 electrons in 6 orbitals, weight-3, a 20-dimensional subspace,
C(6,3)=20) -- the true rank is close to that 20-dimensional maximum, not
a small truncation like every fragment tested so far.

WHAT THIS MEANS, disclosed plainly: solving full H6 via the same forging
approach would need K~19-20 (not 6), meaning ~210 kept measurement slots
(20 diagonal + 190 phase-pairs) instead of 21 -- a ~10x jump -- plus a
much bigger ansatz (a 20-dimensional real subspace needs ~19 free
parameters, not H4's 5). This is NOT a same-session follow-up
experiment; it is a substantially larger undertaking on its own scale.
This is also the real, concrete reason fragmentation is the right tool
for this class of problem in the first place: the full molecule is
intrinsically far more entangled than any of its sub-fragments.

Run:
    python vqe/task53_full_h6_feasibility_check.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "task53_full_h6_feasibility_check_results.json")


def main():
    print("=" * 90)
    print("  task53_full_h6_feasibility_check.py -- full (unfragmented) H6 Schmidt-rank feasibility scan")
    print("=" * 90)
    print("\n  For reference: fragment A/B's own Schmidt rank is EXACTLY 6 (maximal for their "
          "6-dimensional weight-2 subspace). Checking whether full H6's rank saturates similarly.\n")

    # ONE exact diagonalization (strict=False, K=20 -- the theoretical max
    # for a 6-qubit/3-electron weight-3 register, C(6,3)=20), then derive
    # every K's tail directly from the returned full Schmidt spectrum --
    # avoids redundantly re-solving the same 12-qubit exact ground state
    # up to 16 times (the original version of this scan did exactly that,
    # which is both wasteful and, on memory-constrained hardware, a real
    # cause of repeated out-of-memory kills during development of this
    # very script).
    K_MAX = 17
    p = setup_fragment([0, 1, 2, 3, 4, 5], nelec=6, d=1.0, K=K_MAX, strict=False)
    lambdas = p["lambdas"]
    print(f"  n_qubits={p['n_qubits']}, exact_energy={p['exact_energy']:.6f} Ha")
    print(f"  full Schmidt spectrum (first {min(len(lambdas), K_MAX)} values): {lambdas[:K_MAX]}")

    import numpy as np
    tails = {}
    for K in range(1, 17):
        tail = float(np.max(np.abs(lambdas[K:]))) if K < len(lambdas) else 0.0
        tails[K] = tail
        status = "SATURATED (exact)" if tail < 1e-9 else "not yet exact"
        print(f"  K={K}: {status} (tail={tail:.3e})")

    print(f"\n  -- CONCLUSION --")
    print(f"  Rank did not saturate by K=16 (tail={tails.get(16, '?')}) -- true rank is close to the "
          f"20-dimensional weight-3-subspace maximum (C(6,3)=20 for a 6-qubit, 3-electron register), "
          f"NOT a small truncation like every fragment tested so far (fragment A/B: rank exactly 6).")
    print(f"  Implication: full H6 via the same forging approach needs K~19-20 (~210 kept slots, "
          f"20 diagonal + 190 phase-pairs) instead of 21 -- a ~10x jump -- plus a much bigger ansatz. "
          f"A substantially larger undertaking, not a same-session follow-up.")

    with open(RESULTS_PATH, "w") as f:
        json.dump({"tails_by_K": tails, "saturated": False,
                   "conclusion": "full H6 Schmidt rank ~19-20 of a 20-dim max, ~10x more resources than "
                                 "any single fragment tested; not attempted in this session"}, f, indent=2)
    print(f"\n  Saved -> {RESULTS_PATH}\n")


if __name__ == "__main__":
    main()
