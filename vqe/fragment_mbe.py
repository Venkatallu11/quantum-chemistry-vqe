#!/usr/bin/env python3
"""
fragment_mbe.py — Many-Body Expansion (MBE) for molecular clusters.
=====================================================================
PLAIN ENGLISH EXPLANATION OF WHAT THIS FILE DOES:

Exact quantum chemistry (the "chem.py" pipeline used in molecules_real.py)
needs one qubit per spin-orbital. A cluster of many small molecules needs
too many qubits to simulate all at once on a laptop (e.g. 6 H2 molecules
= 24 qubits, way too slow).

The trick, called a "Many-Body Expansion" (MBE), is to break the big
cluster into small pieces ("fragments") and only ever solve SMALL exact
problems:
  1. Solve each fragment alone (a "monomer" energy).
  2. Solve every PAIR of fragments together (a "dimer" energy).
  3. The energy of the whole cluster is approximated as:

        E_MBE = (sum of monomer energies)
              + (sum over every pair of: pair_energy - the two monomer energies)

     The second term is exactly the 2-body "interaction energy" between
     each pair — how much the energy changes because those two fragments
     can see each other. Adding it back on top of the plain monomer sum
     is what captures most of the real physics (van der Waals attraction,
     hydrogen bonding, etc.) without ever building the full cluster's
     Hamiltonian.

  This is a 2-BODY truncation: it misses 3-body-and-higher effects (three
  or more fragments interacting all at once), but those are usually tiny,
  so 2-body MBE is a very good and very cheap approximation.

Every energy number in this file (monomer or pair) is computed by the
EXACT SAME pipeline as molecules_real.py:
    chem.integrals()  ->  chem.rhf()  ->  build qubit Hamiltonian
    (qiskit_nature ElectronicEnergy.from_raw_integrals + JordanWignerMapper)
    ->  exact ground state via scipy.sparse.linalg.eigsh (smallest eigenvalue)
    ->  + nuclear repulsion energy
No shortcuts, no fitted numbers — only the small pieces are separately
exact-diagonalized instead of the whole cluster at once.

Run:
    python vqe/fragment_mbe.py
"""
import os, sys
import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(__file__))
import chem
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper


def exact_energy(geom, nelec):
    """
    Compute the EXACT ground-state energy of a geometry, full space, no
    approximation beyond the STO-3G basis itself. This is the identical
    pipeline used for the "active=None" (full space) molecules in
    molecules_real.py:

      1. chem.integrals(geom)   -> one- and two-electron integrals (real, from geometry)
      2. chem.rhf(...)          -> Hartree-Fock mean-field solution (gives us the
                                    molecular-orbital coefficients C)
      3. Rotate the integrals into the molecular-orbital basis (h1, h2)
      4. Build the qubit Hamiltonian for those integrals with qiskit_nature
         (ElectronicEnergy.from_raw_integrals + Jordan-Wigner mapping)
      5. Find the lowest eigenvalue of that qubit Hamiltonian with a sparse
         eigensolver (scipy eigsh) -- this is the EXACT electronic energy
      6. Add the nuclear repulsion energy (the classical charge-charge
         repulsion between the atomic nuclei) to get the total energy.
    """
    # Steps 1-2: integrals + Hartree-Fock, straight from the atomic geometry
    S, T, V, eri, enuc = chem.integrals(geom)
    ehf, C, Hc = chem.rhf(S, T, V, eri, enuc, nelec=nelec)

    # Step 3: rotate the atomic-orbital integrals into the molecular-orbital
    # basis found by Hartree-Fock (this is what the qubit Hamiltonian needs)
    h1 = C.T @ Hc @ C
    h2 = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, eri, C, C)

    # Step 4: build the qubit Hamiltonian for the FULL orbital space
    ee = ElectronicEnergy.from_raw_integrals(h1, h2)
    qop = JordanWignerMapper().map(ee.second_q_op())
    M = qop.to_matrix(sparse=True)

    # Step 5: the exact ground state is just the smallest eigenvalue of
    # that (sparse) qubit Hamiltonian matrix
    e_elec = float(eigsh(M, k=1, which="SA", return_eigenvectors=False)[0])

    # Step 6: total energy = electronic energy + nuclear repulsion
    return e_elec + enuc


def mbe_2body(fragments):
    """
    2-body Many-Body Expansion.

    `fragments` is a list of (geometry, n_electrons) tuples -- one entry
    per fragment (e.g. one entry per H2 molecule in a cluster).

    Plain-English recipe:
      - First, solve each fragment completely on its own (the "monomer"
        energies). Add them all up -> this is the naive "fragments don't
        talk to each other" estimate.
      - Then, for every PAIR of fragments (i, j) with i < j: glue their
        two geometries together into one bigger geometry, and solve that
        combined system exactly. Whatever energy that pair has ABOVE what
        the two monomers had separately (pair_energy - monomer_i - monomer_j)
        is the "2-body interaction energy" for that pair -- i.e. how much
        those two fragments' presence near each other actually changes
        the energy.
      - Add up all of those pairwise interaction corrections and add them
        on top of the monomer sum. That total is the 2-body MBE estimate
        of the whole cluster's energy.

    Returns a dict with the pieces so callers can report on them:
      monomer_sum        - sum of all monomer (fragment-alone) energies
      pairwise_interaction- sum of all (E_pair - E_i - E_j) corrections
      E_mbe               - monomer_sum + pairwise_interaction (the answer)
      monomer_energies    - list of each individual monomer energy
    """
    n_fragments = len(fragments)

    # --- Step 1: solve every fragment by itself ---
    monomer_energies = [exact_energy(geom, nelec) for geom, nelec in fragments]
    monomer_sum = sum(monomer_energies)

    # --- Step 2: solve every PAIR of fragments together, and keep only
    #             the extra ("interaction") energy each pair brings ---
    pairwise_interaction = 0.0
    for i in range(n_fragments):
        for j in range(i + 1, n_fragments):
            geom_i, nelec_i = fragments[i]
            geom_j, nelec_j = fragments[j]
            # glue the two fragments' atoms into one combined geometry
            geom_pair = geom_i + geom_j
            nelec_pair = nelec_i + nelec_j
            e_pair = exact_energy(geom_pair, nelec_pair)
            # how much extra energy comes purely from i and j interacting
            pairwise_interaction += e_pair - monomer_energies[i] - monomer_energies[j]

    E_mbe = monomer_sum + pairwise_interaction

    return {
        "monomer_sum": monomer_sum,
        "pairwise_interaction": pairwise_interaction,
        "E_mbe": E_mbe,
        "monomer_energies": monomer_energies,
    }


def make_h2_cluster(n, internal=0.74, sep=2.0):
    """
    Build a simple test cluster: n separate H2 molecules, lined up along
    the x-axis. Each H2 becomes its own fragment (2 atoms, 2 electrons).

      internal - the H-H bond length INSIDE each H2 molecule (Angstrom,
                 0.74 A is the real equilibrium H2 bond length)
      sep      - the spacing (Angstrom) between the START of one H2
                 molecule and the START of the next one along x

    Example for n=3: three H2 molecules sitting one after another along
    the x-axis, each internally bonded, with 2.0 A between their
    starting atoms.
    """
    fragments = []
    for i in range(n):
        x0 = i * sep
        geom = [("H", (x0, 0.0, 0.0)), ("H", (x0 + internal, 0.0, 0.0))]
        fragments.append((geom, 2))
    return fragments


def run_verify(n=3):
    """
    Sanity check: for a SMALL cluster (small enough that solving the
    whole thing exactly is still possible), compute:
      (a) the naive monomer sum (fragments pretending not to interact)
      (b) the 2-body MBE estimate (this file's approximation)
      (c) the TRUE full exact energy (every atom solved together, no
          fragmentation at all)

    Comparing (b) to (c) tells us how good the 2-body MBE approximation
    really is -- if MBE is close to the truth, we can trust it on bigger
    clusters where computing the truth (c) is impossible.
    """
    fragments = make_h2_cluster(n)

    # (a) + (b): the 2-body MBE calculation (only ever solves
    # monomers and pairs -- small qubit counts)
    result = mbe_2body(fragments)
    monomer_sum = result["monomer_sum"]
    E_mbe = result["E_mbe"]

    # (c): solve the ENTIRE cluster as one single system, exactly.
    # This is only feasible because n is small here (n=3 -> 6 H atoms).
    full_geom = []
    full_nelec = 0
    for geom, nelec in fragments:
        full_geom += geom
        full_nelec += nelec
    E_full = exact_energy(full_geom, full_nelec)

    # How wrong is the naive "fragments don't interact" guess?
    monomer_error = E_full - monomer_sum
    # How wrong is the 2-body MBE guess? (should be MUCH smaller)
    mbe_error = E_full - E_mbe
    # Of all the interaction energy that exists (monomer_error), what
    # fraction did the cheap 2-body MBE correction actually capture?
    pct_recovered = (1.0 - abs(mbe_error) / abs(monomer_error)) * 100.0

    print(f"\n=== run_verify(n={n}): {n} H2 molecules, checked against full exact ===")
    print(f"  monomer sum (no interaction)   = {monomer_sum:.6f} Ha   "
          f"(error vs full = {monomer_error:+.6f} Ha)")
    print(f"  MBE 2-body estimate            = {E_mbe:.6f} Ha   "
          f"(error vs full = {mbe_error:+.6f} Ha)")
    print(f"  FULL exact (whole cluster)     = {E_full:.6f} Ha")
    print(f"  -> 2-body MBE recovers {pct_recovered:.1f}% of the total interaction energy\n")

    return monomer_sum, E_mbe, E_full, pct_recovered


def run_scale(n=6):
    """
    Show why MBE matters: for n=6 H2 molecules, the FULL cluster has
    6 * 4 = 24 qubits, which is far too slow to exactly diagonalize on a
    laptop (the qubit Hamiltonian matrix would be 2^24 x 2^24). So we
    NEVER build the full system here.

    Instead, MBE only ever needs to solve:
      - monomers: 4 qubits each (fast)
      - pairs: 8 qubits each (also fast)
    and combines those small solves into an estimate for the whole
    24-qubit system we could never solve directly.
    """
    fragments = make_h2_cluster(n)
    result = mbe_2body(fragments)

    total_qubits_if_full = n * 4  # each H2 monomer is 4 qubits (2 orbitals x 2 spins)
    print(f"\n=== run_scale(n={n}): {n} H2 molecules ===")
    print(f"  Full exact treatment would need {total_qubits_if_full} qubits -> infeasible, not attempted.")
    print(f"  MBE only ever solved fragments of <= 8 qubits (monomers=4, pairs=8).")
    print(f"  monomer sum          = {result['monomer_sum']:.6f} Ha")
    print(f"  pairwise interaction = {result['pairwise_interaction']:.6f} Ha")
    print(f"  MBE 2-body estimate  = {result['E_mbe']:.6f} Ha\n")

    return result["E_mbe"]


def main():
    # Step 1: small enough to double-check MBE against the real, full answer
    run_verify(3)
    # Step 2: too big to solve exactly -- MBE is the only option, run it anyway
    run_scale(6)


if __name__ == "__main__":
    main()
