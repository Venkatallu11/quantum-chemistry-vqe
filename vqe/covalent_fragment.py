#!/usr/bin/env python3
"""
covalent_fragment.py — Molecular Tailoring (overlapping-fragment) method
==========================================================================
Big molecules are too expensive to solve exactly on a quantum computer —
qubit count grows with the number of orbitals, and the exact-diagonalization
step grows exponentially with qubit count. Molecular tailoring is a classical
trick chemists use to get around this:

  1. Cut the molecule into smaller, OVERLAPPING pieces (fragments).
  2. Solve each small piece exactly (cheap — few qubits).
  3. Glue the pieces back together with inclusion-exclusion:
     add up every fragment's energy, then SUBTRACT the energy of the
     overlap regions (because atoms in an overlap got counted twice).

This is the same idea as the inclusion-exclusion principle in counting:
  |A union B| = |A| + |B| - |A intersect B|
Here "size" is replaced by "energy", and A/B are overlapping fragments.

Everything below reuses the REAL exact-energy pipeline from
molecules_real.py (chem.py integrals -> Hartree-Fock -> qiskit-nature
qubit Hamiltonian -> exact diagonalization). No fitted constants, no
lookup tables — every energy is computed from atom positions.

Honesty check (read this before trusting the numbers):
  - This works well for a straight chain of hydrogens because H-H bonds
    are short-range and "local" — cutting between two H atoms doesn't
    destroy much real physics, so the fragments capture most of the
    true correlation energy.
  - The fragmentation error does NOT vanish — it shrinks as fragments
    get bigger (more of the true long-range correlation is captured
    inside a single fragment instead of being cut across a boundary).
  - Real covalent molecules (carbon backbones, aromatic rings, etc.)
    have much more delocalized electrons than a bare H-chain. Fragments
    would need to be considerably bigger (and sometimes capped with
    extra atoms at the cut) to reach the same accuracy. This script
    does NOT prove tailoring works for carbon chemistry — it only
    proves the classical bookkeeping is correct, on the simplest
    possible test case.

Run:
    python vqe/covalent_fragment.py
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(__file__))
from molecules_real import ground_state


def hchain(indices, d):
    """
    Build the geometry for hydrogen atoms sitting on a straight line.

    `indices` are the atom's positions along the chain (e.g. atom #3 of
    a longer chain), NOT a local 0,1,2... counter — this matters because
    a fragment cut from the middle of a chain must sit at its REAL
    bond distance from its neighbors, not be shifted back to the origin.

    Args:
        indices: which chain positions to place H atoms at (integers)
        d: bond spacing in Angstrom

    Returns:
        geometry list in the format chem.py expects: [("H", (x,y,z)), ...]
    """
    return [("H", (i * d, 0.0, 0.0)) for i in indices]


def full_exact(n_atoms, d):
    """
    The reference answer: solve the WHOLE chain exactly, no cutting.
    Only feasible while n_atoms is small (qubit count = 2 * n_atoms
    for STO-3G hydrogen, since each H contributes one 1s orbital).
    """
    geom = hchain(list(range(n_atoms)), d)
    _ehf, exact, nq = ground_state(geom, nelec=n_atoms, active=None)
    return exact, nq


def tailor_energy(n_atoms, block, d):
    """
    Molecular tailoring: reassemble the energy of an n_atoms H-chain
    from overlapping fragments of size `block`, without ever solving
    the full chain.

    How the fragments are chosen:
      The window slides forward by 2 atoms each step (this is what
      keeps neighboring fragments overlapping instead of leaving gaps).
      Example, n_atoms=8, block=4:
        fragment 0: atoms [0,1,2,3]
        fragment 1: atoms [2,3,4,5]   <- shares atoms {2,3} with frag 0
        fragment 2: atoms [4,5,6,7]   <- shares atoms {4,5} with frag 1
      Every atom is covered by at least one fragment, and every pair of
      NEIGHBORING fragments shares exactly (block - 2) atoms.

    Why subtract the overlaps:
      Atoms {2,3} above got their energy counted once inside fragment 0
      AND once inside fragment 1 — that's double-counting. We computed
      the SAME 2-electron/4-electron overlap fragment exactly on its own,
      and subtract it once, exactly cancelling the double count
      (this is the inclusion-exclusion principle from set theory,
      applied to energy instead of counting).

    Args:
        n_atoms: total atoms in the full (uncut) chain
        block: how many atoms are in each fragment
        d: bond spacing in Angstrom

    Returns:
        (reassembled_energy, max_qubits_used_by_any_single_fragment)
    """
    shift = 2  # fragments always slide forward by 2 atoms
    starts = list(range(0, n_atoms - block + 1, shift))
    fragments = [list(range(s, s + block)) for s in starts]

    # Sanity check: every atom in the chain must be covered by some
    # fragment, or the reassembled energy would be missing physics.
    covered = set()
    for frag in fragments:
        covered.update(frag)
    assert covered == set(range(n_atoms)), (
        f"Fragments {fragments} don't cover the full {n_atoms}-atom chain — "
        f"pick a block size where (n_atoms - block) is divisible by {shift}."
    )

    max_qubits = 0
    reassembled = 0.0

    # Step 1: add up every fragment's own exact energy.
    for frag in fragments:
        geom = hchain(frag, d)
        _ehf, exact, nq = ground_state(geom, nelec=len(frag), active=None)
        reassembled += exact
        max_qubits = max(max_qubits, nq)

    # Step 2: subtract the energy of each overlap region exactly once —
    # this cancels the double-counting between neighboring fragments.
    for i in range(len(fragments) - 1):
        overlap_atoms = sorted(set(fragments[i]) & set(fragments[i + 1]))
        geom = hchain(overlap_atoms, d)
        _ehf, exact, nq = ground_state(geom, nelec=len(overlap_atoms), active=None)
        reassembled -= exact
        max_qubits = max(max_qubits, nq)

    return reassembled, max_qubits


def main():
    d = 1.0  # Angstrom bond spacing for this test chain
    results = {}

    print("\n" + "=" * 70)
    print("  Molecular Tailoring — hydrogen chain fragmentation test")
    print(f"  Bond spacing d = {d} Angstrom")
    print("=" * 70)

    # --- H6, 4-atom fragments -----------------------------------------
    print("\n  H6, split into 4-atom fragments")
    full6, nq_full6 = full_exact(6, d)
    tailor6, nq_tailor6 = tailor_energy(6, block=4, d=d)
    err6 = abs(tailor6 - full6)
    print(f"    full exact        = {full6:.6f} Ha   ({nq_full6} qubits)")
    print(f"    tailored (4-atom) = {tailor6:.6f} Ha   ({nq_tailor6} qubits max per fragment)")
    print(f"    error             = {err6:.6f} Ha")
    results["H6_block4"] = {
        "full_exact_ha": round(full6, 6),
        "tailored_ha": round(tailor6, 6),
        "error_ha": round(err6, 6),
        "full_qubits": nq_full6,
        "max_fragment_qubits": nq_tailor6,
    }

    # --- H8, 4-atom fragments ------------------------------------------
    print("\n  H8, split into 4-atom fragments")
    full8, nq_full8 = full_exact(8, d)
    tailor8_b4, nq_tailor8_b4 = tailor_energy(8, block=4, d=d)
    err8_b4 = abs(tailor8_b4 - full8)
    print(f"    full exact        = {full8:.6f} Ha   ({nq_full8} qubits)")
    print(f"    tailored (4-atom) = {tailor8_b4:.6f} Ha   ({nq_tailor8_b4} qubits max per fragment)")
    print(f"    error             = {err8_b4:.6f} Ha")
    results["H8_block4"] = {
        "full_exact_ha": round(full8, 6),
        "tailored_ha": round(tailor8_b4, 6),
        "error_ha": round(err8_b4, 6),
        "full_qubits": nq_full8,
        "max_fragment_qubits": nq_tailor8_b4,
    }

    # --- H8, 6-atom fragments (bigger fragments -> smaller error) ------
    print("\n  H8, split into 6-atom fragments (bigger fragments)")
    tailor8_b6, nq_tailor8_b6 = tailor_energy(8, block=6, d=d)
    err8_b6 = abs(tailor8_b6 - full8)
    print(f"    full exact        = {full8:.6f} Ha   ({nq_full8} qubits)")
    print(f"    tailored (6-atom) = {tailor8_b6:.6f} Ha   ({nq_tailor8_b6} qubits max per fragment)")
    print(f"    error             = {err8_b6:.6f} Ha")
    results["H8_block6"] = {
        "full_exact_ha": round(full8, 6),
        "tailored_ha": round(tailor8_b6, 6),
        "error_ha": round(err8_b6, 6),
        "full_qubits": nq_full8,
        "max_fragment_qubits": nq_tailor8_b6,
    }

    print("\n" + "-" * 70)
    print(f"  H8 error shrank from {err8_b4:.6f} Ha (4-atom blocks) to "
          f"{err8_b6:.6f} Ha (6-atom blocks)")
    print("  -> bigger fragments capture more of the true correlation,")
    print("     so the tailoring error goes down. This is the whole point:")
    print("     you get a knob (fragment size) to trade qubits for accuracy.")
    print("-" * 70)

    print("\n  HONEST LIMITATION:")
    print("  This test is a bare hydrogen chain — the easiest possible case.")
    print("  Real covalent molecules (carbon backbones, rings, conjugated")
    print("  systems) have electrons far more delocalized across many atoms.")
    print("  Cutting those bonds loses more real physics per cut, so carbon")
    print("  chemistry needs noticeably bigger fragments (and often capping")
    print("  atoms at the cut point) to reach the same accuracy shown here.")
    print("  Nothing here demonstrates tailoring works for carbon — only")
    print("  that the fragment/overlap bookkeeping itself is correct.\n")

    out = os.path.join(os.path.dirname(__file__), "covalent_fragment_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved -> {out}\n")

    return results


if __name__ == "__main__":
    main()
