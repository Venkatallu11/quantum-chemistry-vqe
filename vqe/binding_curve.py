#!/usr/bin/env python3
"""
binding_curve.py — Binding-energy curve of a 3-molecule H2 cluster.
=====================================================================
PLAIN ENGLISH EXPLANATION OF WHAT THIS FILE DOES:

Take three H2 molecules and line them up along the x-axis. If we push them
close together, they repel each other hard (their electron clouds and nuclei
don't want to overlap). If we pull them far apart, they stop noticing each
other at all, and the cluster's energy should just settle down to "three
separate, non-interacting H2 molecules."

A "binding-energy curve" is just: compute the cluster's total energy at a
bunch of different separations, and compare each one to the "three
completely isolated H2 molecules" baseline. That comparison is the
"binding energy" -- how much extra (or less) energy the cluster has because
the molecules are near each other, as a function of distance.

This script does that for real, using EXACT diagonalization every time (no
approximations, no fragmentation, no fitting) -- it solves the FULL 6-atom,
12-qubit system directly at each separation. That's only possible because 3
H2 molecules is still small enough for a laptop to diagonalize exactly; this
script is the "ground truth" that things like fragment_mbe.py are checked
against.

Every energy number here comes from the EXACT SAME pipeline used everywhere
else in this repo (see molecules_real.py, reused directly, not
reimplemented):
    chem.integrals()  ->  chem.rhf()  ->  build qubit Hamiltonian
    (qiskit_nature ElectronicEnergy.from_raw_integrals + JordanWignerMapper)
    ->  exact ground state via scipy.sparse.linalg.eigsh (smallest eigenvalue)
    ->  + nuclear repulsion energy
No PySCF anywhere in this compute path -- PySCF is only ever used (elsewhere
in this repo) to double check answers after the fact.

HONEST LIMITATION: this uses the STO-3G basis set, which is a small, "minimal"
basis. STO-3G captures the strong, short-range physics correctly (electron/
nuclear repulsion pushing the molecules apart when they're squeezed together),
but it does NOT capture long-range dispersion (van der Waals attraction --
the very weak "London forces" that make real H2 molecules attract each other
a tiny bit at large distances). So don't expect this curve to dip negative
at large separation the way a real experimental H2-H2 potential would -- it's
a known, honest limitation of the small basis set, not a bug.

Run:
    python vqe/binding_curve.py
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")  # write straight to a file, no GUI window needed
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import molecules_real as mr  # reuses ground_state() -- the exact pipeline, not reimplemented

# Converts Hartree (the atomic unit of energy chemistry uses internally)
# into kcal/mol (the unit chemists usually talk about reaction energies in).
HARTREE_TO_KCAL = 627.509

INTERNAL_BOND = 0.74  # H-H bond length INSIDE each H2 molecule, Angstrom (real equilibrium value)

# Center-to-center separations between neighboring H2 molecules to scan, Angstrom.
# Small values = molecules squeezed together (should repel hard).
# Large values = molecules far apart (should approach "no interaction").
SEPARATIONS = [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]


def isolated_h2_energy():
    """
    Exact energy of ONE isolated H2 molecule at its equilibrium bond length
    (0.74 A). This is our "zero interaction" reference point: if the three
    H2 molecules in the cluster genuinely couldn't see each other at all,
    the cluster's total energy would be exactly 3x this number.
    """
    geom = [("H", (0.0, 0.0, 0.0)), ("H", (INTERNAL_BOND, 0.0, 0.0))]
    # active=None -> full space, no approximation (same pipeline as molecules_real.py)
    _, exact, _ = mr.ground_state(geom, nelec=2, active=None)
    return exact


def cluster_energy(sep):
    """
    Build THREE H2 molecules in a line along x, each internally bonded at
    0.74 A, with sep Angstrom between each molecule's starting atom and the
    next molecule's starting atom. Then solve the WHOLE 6-atom, 12-qubit
    system EXACTLY -- no fragmentation, no approximation. This is the
    ground-truth energy of the cluster at this particular separation.
    """
    geom = []
    for i in range(3):
        x0 = i * sep
        geom.append(("H", (x0, 0.0, 0.0)))
        geom.append(("H", (x0 + INTERNAL_BOND, 0.0, 0.0)))
    ehf, exact, nq = mr.ground_state(geom, nelec=6, active=None)
    return exact, nq


def main():
    # Step 1: the reference point -- one H2 molecule, all by itself.
    e_h2 = isolated_h2_energy()
    e_h2_x3 = 3.0 * e_h2
    print(f"\nIsolated H2 energy  = {e_h2:.6f} Ha  (reference, one molecule)")
    print(f"3 x isolated H2     = {e_h2_x3:.6f} Ha  ('no interaction' baseline)\n")

    results = {
        "isolated_h2_ha": round(e_h2, 6),
        "three_times_isolated_ha": round(e_h2_x3, 6),
        "points": [],
    }

    # Step 2: for every separation, solve the FULL 6-atom cluster exactly,
    # and see how far its energy is from the "3 separate molecules" baseline.
    print(f"{'sep(A)':>7} | {'E_cluster(Ha)':>14} | {'binding(Ha)':>12} | {'binding(kcal/mol)':>18}")
    for sep in SEPARATIONS:
        e_cluster, nq = cluster_energy(sep)

        # Binding energy = how much MORE (or less) energy the cluster has
        # compared to three molecules that can't see each other at all.
        # Positive = the cluster is higher energy than isolated (net repulsion
        # at this separation); negative would mean net attraction/binding.
        binding_ha = e_cluster - e_h2_x3
        binding_kcal = binding_ha * HARTREE_TO_KCAL

        print(f"{sep:7.1f} | {e_cluster:14.6f} | {binding_ha:+12.6f} | {binding_kcal:+18.3f}")

        results["points"].append({
            "sep_angstrom": sep,
            "E_cluster_ha": round(e_cluster, 6),
            "binding_ha": round(binding_ha, 6),
            "binding_kcal_mol": round(binding_kcal, 3),
            "qubits": nq,
        })

    # Step 3: save the raw numbers so this run is reproducible / inspectable
    # without re-running the (slow) exact diagonalizations.
    out_json = os.path.join(os.path.dirname(__file__), "binding_curve_results.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved -> {out_json}")

    # Step 4: plot binding energy (kcal/mol) vs separation (Angstrom).
    # Expect a steep positive wall at small separation (hard repulsion) that
    # decays toward ~0 as separation grows (molecules stop noticing each
    # other) -- but see the HONEST LIMITATION note at the top of this file:
    # STO-3G won't show the small negative (attractive) dip a real dispersion
    # interaction would produce at large separation.
    seps = [p["sep_angstrom"] for p in results["points"]]
    bindings = [p["binding_kcal_mol"] for p in results["points"]]

    plt.figure(figsize=(7, 5))
    plt.plot(seps, bindings, marker="o")
    plt.axhline(0, color="gray", linewidth=0.8)
    plt.xlabel("Separation between H2 molecules (Angstrom)")
    plt.ylabel("Binding energy (kcal/mol)")
    plt.title("Binding-energy curve: 3 H2 molecules in a line (exact, STO-3G)")
    plt.grid(True, alpha=0.3)

    out_png = os.path.join(os.path.dirname(__file__), "binding_curve.png")
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Plot saved -> {out_png}\n")


if __name__ == "__main__":
    main()
