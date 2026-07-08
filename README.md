# Quantum Chemistry VQE

*A from-scratch quantum chemistry engine — exact molecular ground-state energies computed from geometry, verified against PySCF, and run on real IBM Quantum hardware.*

Every simulator result in this repo has been independently checked against
[PySCF](https://pyscf.org/), the standard professional reference for
quantum chemistry. PySCF is **not** required to run any of the code here —
it's only used, off to the side, to double-check that the numbers this repo
produces are correct. Beyond simulation, this repo has also taken real
circuits to real IBM Quantum hardware, and reports what happened honestly —
including where today's hardware falls short.

---

## Why it matters

Quantum chemistry is the bottleneck behind some of the most important stuck
problems in science:

- **Nitrogen fixation.** The enzyme FeMoco (nitrogenase's active site) turns
  atmospheric nitrogen into ammonia at room temperature and pressure — a
  reaction we still don't fully understand. Instead, industry uses the
  Haber-Bosch process, which has burned roughly **2% of the world's energy
  supply every year since 1909**, because we can't classically simulate what
  FeMoco is actually doing.
- **Drug discovery** depends on accurately predicting how molecules bind —
  which is a quantum chemistry problem at its core.
- **Superconductors** and other advanced materials are governed by electron
  correlations that classical computers can only approximate.

The common thread: these are all systems where electrons are strongly
correlated, and classical computers can't simulate that exactly once the
system gets big enough — the cost grows exponentially. Quantum computers,
in principle, don't have that problem.

This project builds and verifies the lower rungs of that ladder — small
molecules where an exact classical answer is still possible, so every
result can be checked — as a path toward the systems that currently can't
be checked at all:

```
H2 → LiH → H2O → NH3 → N2 → CH4 → ... → FeMoco (54 qubits)
```

---

## Verified simulator results

Every number below was computed from geometry using `chem.py`'s integral
engine → Hartree-Fock → qubit Hamiltonian → exact diagonalization, and
matches PySCF exactly.

| Molecule | Qubits | Exact ground state (Ha) | Reference |
|---|---|---|---|
| H2  | 4  | -1.137284   | PySCF FCI |
| LiH | 12 | -7.882324   | PySCF FCI |
| H2O | 14 | -75.012647  | PySCF FCI |
| NH3 | 16 | -55.511677  | PySCF FCI |
| N2  | 12 | -107.538421 | PySCF CASCI(10e,6o) — frozen-core active space |
| CH4 | 16 | -39.805581  | PySCF CASCI(8e,8o) — frozen-core active space |

**Ions and atoms:**

| Species | Qubits | Exact ground state (Ha) |
|---|---|---|
| He   | 2 | -2.807784 |
| HeH+ | 4 | -2.851024 |
| H3+  | 6 | -1.261200 |

**Dissociation / binding curve** (`binding_curve.py`): a 3-molecule H2
cluster solved *exactly* (12 qubits, no approximation) at 7 separations,
compared against 3× isolated H2 (-3.411852 Ha):

| Separation (Å) | Binding energy (kcal/mol) |
|---|---|
| 1.2 | +508.648 |
| 1.5 | +157.377 |
| 2.0 | +28.747  |
| 2.5 | +5.143   |
| 3.0 | +0.793   |
| 4.0 | +0.010   |
| 5.0 | +0.002   |

Steep repulsion when squeezed together, decaying toward zero as the
molecules separate — exactly the expected short-range behavior (STO-3G
doesn't capture the weak long-range attraction; see Honesty section).

---

## Fragmentation: making bigger molecules tractable

Exact diagonalization needs one qubit per spin-orbital, so anything much
past ~16 qubits is out of reach on a laptop. Two fragmentation strategies
let this repo reach *much* larger systems using only ≤8-qubit pieces:

**Many-Body Expansion** (`fragment_mbe.py`) — for clusters of separate
molecules (e.g. several H2 molecules near each other): solve monomers,
then pairs, then triples, and sum the corrections. Verified against a full
exact calculation: 2-body MBE recovers **99.7%** of the interaction energy
for a 3-molecule cluster, and 2-body/3-body recovers **99.52%/99.96%** for
a 4-molecule cluster — evidence the expansion converges toward the truth,
used later (`run_scale`, 6 molecules / 24 qubits) where a full check is no
longer feasible.

**Covalent-bond fragmentation / molecular tailoring** (`covalent_fragment.py`)
— for a single *bonded* chain (can't just split into separate molecules;
overlapping fragments + inclusion-exclusion instead):

| Chain | Fragment size | Full exact (Ha) | Tailored (Ha) | Error (Ha) |
|---|---|---|---|---|
| H6 | 4-atom blocks | -3.236066 | -3.231625 | 0.004442 |
| H8 | 4-atom blocks | -4.307572 | -4.296862 | 0.010710 |
| H8 | 6-atom blocks | -4.307572 | -4.305745 | 0.001826 |

The error **shrinks systematically as fragments grow** (0.0107 → 0.0018 Ha
going from 4- to 6-atom blocks) — a real, controllable accuracy knob, not
a fixed approximation. This is the mechanism that, in principle, scales to
hundreds of qubits for favorable (weakly-coupled-fragment) systems while
every individual solve stays small enough to run today.

---

## Real hardware results

**H2, measured on `ibm_marrakesh`** — the verified H2 Hamiltonian
(exact = **-1.137284 Ha**), an `EfficientSU2` ansatz optimized on a
noiseless simulator, then measured ONCE on real hardware two independent
ways:

| Path | Hardware energy (Ha) | Error vs exact (Ha) | Error (kcal/mol) | Job ID |
|---|---|---|---|---|
| Direct (qiskit-ibm-runtime) | -1.116691 | 0.020593 | 12.922 | `d949hudgc6cc73fer3ig` |
| Via Lokesh's MCP server     | -1.117200 | 0.020084 | 12.603 | `d949bolgc6cc73feqt20` |

The two paths agree with each other to **0.000509 Ha** — strong evidence
the MCP layer adds no distortion of its own; it's a faithful wrapper
around the same `EstimatorV2` primitive. Both land ~12.6-12.9 kcal/mol
from exact: an honest noise floor for one un-mitigated 2-qubit hardware
measurement, not a tuned result.

**The 8-qubit fragment frontier** — pushing to `block1` of the H6 chain
(4 electrons, 8 qubits) exposes where today's hardware actually struggles.
A simulator-only ansatz comparison (`fragment_ansatz_test.py`) found:

| Ansatz | Error vs exact (Ha) | Converged (<0.002 Ha)? | Circuit depth | 2-qubit gates |
|---|---|---|---|---|
| UCCSD (+ Hartree-Fock init) | 0.0000794 | Yes | 2168 | **1440** |
| ExcitationPreserving(reps=2) | 0.004749  | No (close) | 232 | **112** |

UCCSD reaches the exact ground state *in simulation* — but needs 1440
two-qubit gates, an amount that would be overwhelmed by real gate error
before a single useful measurement. `ExcitationPreserving` is ~13x
shallower but can't quite close the last 0.0047 Ha gap. There is no free
lunch here yet: accuracy and hardware-feasibility trade directly against
each other at 8 qubits with today's ansätze.

An actual hardware attempt at this fragment (`hardware_mitigation_test.py`,
testing whether resilience_level 0/1/2 rescue it) instead surfaced a real
methodology bug — see Integration below — and its numbers are not a valid
noise-floor measurement as a result. A corrected re-run
(`hardware_covalent.py`, full H6 tailoring on real hardware with the bug
fixed) is the ongoing frontier experiment for this repo; results will be
added here once complete, reported exactly as measured, noisy or not.

---

## Noise-resilience experiments (simulator-only, negative results)

Before spending more real hardware time/credits on the 8-qubit fragment
frontier, three cheap ideas were tested locally (`qiskit-aer` noise
models only — no IBM or Azure calls) to see whether fragmentation error
could be predicted or reduced without paying for more hardware runs. All
three are reported honestly even though none of them worked.

**1. Does noise cancel in the difference structure?**
(`test_difference_cancellation.py`) — H6's tailoring sum
(block1 + block2 − overlap) shares real structural overlap (same H-H
bond, same physical gates) between fragments. Hypothesis: a systematic
(coherent) noise bias might partially cancel in that inclusion-exclusion
sum, unlike purely random (incoherent) noise. Tested two noise models
calibrated to the same 1% two-qubit gate infidelity — a coherent RZZ
over-rotation vs. incoherent depolarizing noise.

Result: **not confirmed.** Reassembled H6 error blew up to 1.15 Ha
(coherent) / 1.28 Ha (incoherent) — 721–806 kcal/mol — with no
cancellation advantage from the coherent case.

**2. Is Quantinuum-grade hardware worth it for this fragment size?**
(`fragmentation_noise_prediction.py`) — before spending real Azure
Quantum credit, a local depolarizing noise model compared Quantinuum-like
fidelity (0.2% 2-qubit error) against this repo's actually-observed
IBM-like fidelity (3% 2-qubit error) on the same H6 tailoring fragment.

Result: Quantinuum-level noise still gives 0.34 Ha reassembled error
(212 kcal/mol); IBM-level noise gives 2.45 Ha (1538 kcal/mol). Even the
better hardware profile lands nowhere near chemical accuracy
(~1.5 kcal/mol) at this fragment size — evidence that spending real
Azure credit here would not currently pay off.

**3. Can the error be extrapolated or averaged away cheaply?**
(`test_consistency_fragmentation.py`) — Part A: does fragment error
shrink predictably enough with overlap size to extrapolate (Aitken's
delta-squared process) toward the near-exact answer from cheap
small-overlap runs? Part B: does simple redundancy averaging (20 samples
under 1% depolarizing noise) shrink shot noise the expected 1/√N way?

Result: Part A — extrapolation (0.0025 Ha error) did not beat the best
individual scheme already measured (0.0018 Ha at block=6). Part B —
averaging only shrank error by 1.02x, far short of the expected
√20 ≈ 4.47x, meaning the noise floor here isn't simple shot noise that
averaging fixes.

**Bottom line:** none of these three cheap simulator-only tricks rescue
the 8-qubit fragment frontier. The honest conclusion from the hardware
results above stands — 8-qubit chemistry needs shallower ansätze or
better gate fidelities, not cleverer classical post-processing of noisy
fragment energies.

---

## Integration with Lokesh's Quantum Hardware MCP server

This repo's chemistry engine is fully independent, but it also connects to
the open-source **Quantum Hardware MCP server** by Lokesh Pullakandam
(credits below) for exactly one purpose: **device selection**
(`mcp_backend.best_device()` → `server.compare_devices`, read-only, no
job submitted). The actual chemistry measurement always runs directly
through this repo's own `qiskit-ibm-runtime` code — see `mcp_energy.py`
for the version that also routes the *measurement* through the MCP server
(`estimate_expectation`), included for comparison.

**Known issue, not yet fixed:** the MCP server's `best_qubits` tool ranks
physical qubits purely by individual readout/gate error, with no check for
whether they're actually connected on the chip. In testing, this picked 8
qubits with only 1 connected pair between them, forcing a flood of SWAP
gates during transpilation and producing hardware energies off by >1 Ha —
not real noise, a qubit-selection bug. The fix (this repo's own
`generate_preset_pass_manager` already does connectivity-aware auto-
placement correctly) is documented and used going forward; `best_qubits`
itself has not been patched.

**Fixed:** `best_device()` (`mcp_backend.py`) now also cross-checks
`list_devices()`'s status field and skips any device in MAINTENANCE.
`compare_devices()` alone can rank a device top by calibration/queue data
while it's actually in a maintenance window (a separate field it doesn't
expose) — discovered when a job submitted to a top-ranked device sat
QUEUED indefinitely during `hardware_covalent.py` testing.

---

## Honesty section

- Every simulator number in this repo is computed from geometry, with no
  hardcoded energies and no tuned constants (the one stated exception is
  `h2_vqe.py`'s original demo Hamiltonian, since replaced by the verified
  one — see `h2_vqe.py`'s own comments), and every number is independently
  verified against PySCF.
- N2 and CH4 both use clearly-labelled **frozen-core active spaces**
  (CASCI(10e,6o) and CASCI(8e,8o) respectively), because their full-space
  problems (20 and 18 qubits) are confirmed intractable to exactly
  diagonalize on a laptop — CH4's full-space attempt failed with a memory
  allocation error during this repo's own testing. This is standard
  practice in real quantum chemistry, not a shortcut, and both are
  verified against PySCF's CASCI result for the identical active space.
- **STO-3G is a minimal basis set** — it captures strong, short-range
  physics correctly but misses long-range dispersion (see the binding
  curve's honest limitation, documented in `binding_curve.py`).
- **Fragmentation methods (MBE, molecular tailoring) are established
  research techniques, not novel science.** The contribution here is
  making them accessible: a from-scratch, verified, openly-packaged
  implementation anyone can run and check against PySCF, rather than a
  new algorithm.
- **8-qubit chemistry is beyond current hardware**, documented rather than
  hidden: see the fragment frontier results above. Only very small
  circuits (H2, ≤4-qubit fragments) currently produce usable real-hardware
  numbers; 8-qubit fragments need either far shallower ansätze than exist
  today or meaningfully better gate fidelities.
- **Known limits:** closed-shell molecules only (no open-shell radicals);
  full exact diagonalization tops out around 16 qubits on a laptop —
  anything bigger needs an active space or fragmentation.

---

## Setup and usage

```bash
pip install -r requirements.txt
```

```bash
# Simulator — verified chemistry
python vqe/molecules_real.py          # all molecules/ions in the verified table
python vqe/fragment_mbe.py            # many-body expansion / fragmentation demo
python vqe/binding_curve.py           # H2-cluster binding-energy curve + plot
python vqe/covalent_fragment.py       # covalent-bond molecular tailoring (H6/H8)
python vqe/h2_vqe.py                  # VQE demo on H2 (local simulator, or --real for hardware)
python vqe/fragment_ansatz_test.py    # simulator ansatz comparison (ExcitationPreserving vs UCCSD)

# Real IBM Quantum hardware (needs IBM_QUANTUM_TOKEN in .env)
python vqe/h2_hardware_direct.py      # H2 on hardware, direct qiskit-ibm-runtime
python vqe/mcp_energy.py              # H2 on hardware, routed through the MCP server
python vqe/hardware_clean.py          # H2 on hardware, MCP used only for device selection
python vqe/hardware_covalent.py       # full H6 covalent fragmentation on real hardware
```

### Adding your own molecule

`vqe/molecules_real.py` defines a `MOLECULES` dict mapping a name to
`(geometry, n_electrons, active_space)`. To try a new molecule, add an entry
in the same shape:

```python
MOLECULES = {
    ...
    "MyMolecule": ([("O", (0, 0, 0)), ("H", (0, 0, 0.96))], 10, None),
}
```

- `geometry` is a list of `(element, (x, y, z))` tuples in Angstrom.
- `n_electrons` is the total electron count (account for charge if the
  species is an ion).
- `active_space` is `None` for a full-space exact calculation, or
  `(n_core_orbitals, n_active_orbitals)` for a frozen-core active space
  (used for N2 and CH4, since their full spaces are too large).

Then run `python vqe/molecules_real.py` — it computes, prints, and saves the
Hartree-Fock and exact energies for every entry in the dict.

---

## File structure

```text
vqe/
├── chem.py                            # Pure-Python integral engine + Hartree-Fock
├── molecules_real.py                  # Exact ground-state energies for the molecule/ion table
├── molecules_real_results.json        # Saved results from molecules_real.py
├── fragment_mbe.py                    # Many-body expansion (2-body and 3-body) for clusters
├── binding_curve.py                   # H2-cluster binding-energy curve (exact, no fragmentation)
├── binding_curve_results.json         # Saved results from binding_curve.py
├── binding_curve.png                  # Saved plot from binding_curve.py
├── covalent_fragment.py               # Covalent-bond molecular tailoring (H6/H8)
├── covalent_fragment_results.json     # Saved results from covalent_fragment.py
├── fragment_ansatz_test.py            # Simulator ansatz comparison (ExcitationPreserving vs UCCSD)
├── fragment_ansatz_test_results.json  # Saved results from fragment_ansatz_test.py
│
├── h2_vqe.py                          # VQE demo on H2: EfficientSU2 + COBYLA, local sim or real hardware
├── h2_vqe_results.json                # Saved results from h2_vqe.py
├── h2_hardware_full.py                # Full COBYLA optimization loop entirely on real hardware
├── h2_hardware_results.json           # H2 on hardware via the MCP server (mcp_energy.py's result)
├── h2_hardware_direct.py              # H2 on hardware, direct qiskit-ibm-runtime (no MCP)
├── h2_hardware_direct_results.json    # Saved results from h2_hardware_direct.py
│
├── mcp_backend.py                     # Read-only MCP connector: best_device(), best_qubits_for(), validate_circuit()
├── mcp_energy.py                      # H2 hardware measurement routed through the MCP server
├── hardware_clean.py                  # H2 on hardware; MCP used ONLY for device selection
├── hardware_fragmentation.py          # H6 fragmentation on hardware (early attempt, EfficientSU2+penalty)
├── hardware_fragmentation_mcp.py      # Same, routed through the MCP server
├── hardware_mitigation_test.py        # Error-mitigation test on one 8-qubit fragment (found the connectivity bug)
├── hardware_mitigation_results.json   # Saved results from hardware_mitigation_test.py
├── hardware_covalent.py               # Corrected full H6 fragmentation on real hardware (auto-placement, ZNE)
├── hardware_covalent_results.json     # Saved results from hardware_covalent.py (once complete)
│
├── test_difference_cancellation.py    # Simulator test: does noise cancel in the tailoring difference structure?
├── difference_cancellation_results.json # Saved results (hypothesis not confirmed)
├── fragmentation_noise_prediction.py  # Simulator test: Quantinuum-like vs IBM-like noise on H6 tailoring
├── fragmentation_noise_prediction.json # Saved results
├── test_consistency_fragmentation.py  # Simulator test: overlap extrapolation + redundancy averaging on H8
└── consistency_fragmentation_results.json # Saved results

requirements.txt
```

---

## Powered by / credits

Built alongside the open-source **Quantum Hardware MCP server** by
Lokesh Pullakandam, which connects AI assistants to live IBM Quantum
hardware — the bridge that lets these same molecules run on a real
quantum machine instead of only a classical simulation of one.

- MCP server: https://github.com/Lokesh-2025/quantum-hardware-mcp
- Lokesh: https://www.linkedin.com/in/lokesh-pullakandam/

---

## License

MIT — Venkata Rao Allu. See [LICENSE](LICENSE).
