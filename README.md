# Quantum Chemistry VQE

[![Qiskit Ecosystem](https://qisk.it/e-e52f2069)](https://qisk.it/e)

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

## Entanglement forging + zero-noise extrapolation

The noise-resilience experiments above ended on a negative note: even
Quantinuum-grade noise (0.2% two-qubit error) gave 212 kcal/mol of
reassembled error on the 8-qubit H6 tailoring fragment — nowhere close to
chemical accuracy. **Entanglement forging** (`entanglement_forging_h4.py`)
revisits that question with a different qubit-reduction trick: split an
8-qubit fragment (H4, atoms 0-3) into two independent 4-qubit registers
(alpha-spin / beta-spin orbitals), write the ground state as a Schmidt
decomposition `|psi> = sum_n lambda_n |u_n>_alpha |v_n>_beta`, and
reconstruct `<psi|H|psi>` from small circuits measured *separately* on
each 4-qubit register — the two registers never need to be entangled with
each other on real hardware.

**Stage 1 (noiseless, numpy)** confirms the bipartite bookkeeping itself is
correct — truncating to the top `K` Schmidt terms reproduces known
reference errors exactly:

| K | Error vs exact (Ha) |
|---|---|
| 1 | 0.064759 |
| 3 | 0.004396 |
| 5 | 0.000901 |
| 8 (full rank) | 0.000000 |

**Stage 2** replaces the exact matrix elements with ones measured from
actual 4-qubit circuits (`StatePreparation` + the standard EF
superposition trick for off-diagonal terms) run through a `qiskit-aer`
Quantinuum-like depolarizing noise model (0.2% two-qubit, 0.005%
one-qubit). At K=5, un-mitigated noise pushed the error to **20.15
kcal/mol** — still far from chemical accuracy (1 kcal/mol), confirming
raw fragmentation + EF doesn't survive noise on its own.

**Adding ZNE** (`entanglement_forging_zne.py`) closes that gap.
Zero-noise extrapolation reruns the *same* circuits at 1x, 2x, and 3x the
noise strength, then fits energy vs. scale and extrapolates back to the
noiseless limit (scale = 0) — without needing a noise-free device:

| Method | Error (kcal/mol) |
|---|---|
| No mitigation (scale=1) | 20.22 |
| ZNE, linear extrapolation | 1.12 |
| **ZNE, quadratic extrapolation** | **0.57** |

Quadratic ZNE lands *under* the 1 kcal/mol chemical-accuracy threshold —
a 35x reduction in error from the un-mitigated result, using only
classical post-processing of circuits already being run.

**Honest limitation:** both stages use the *exact* Schmidt vectors from
full diagonalization of the real H4 Hamiltonian, not a variationally
trained ansatz — this isolates the noise/mitigation question from the
optimization question. A full EF-VQE (variationally optimizing the two
4-qubit circuits instead of reading them off exact diagonalization) is
the natural next step and is **not** done here. The noise model is also
still a local `qiskit-aer` approximation of Quantinuum hardware, not a
measurement on an actual Quantinuum device — see the Azure Quantum
section below for the real-emulator connection that would let this run
be redone against genuine device noise.

---

## Azure Quantum / Quantinuum connection

To eventually replace the local Quantinuum-*like* noise model above with
real Quantinuum device/emulator noise, this repo now connects to a real
Azure Quantum workspace (`vqe/azure_backend.py`), authenticated via
device-code login (no secrets in code — workspace resource ID, location,
and tenant ID all come from `.env`, see `.env.example`).

`compare_targets()` mirrors the IBM side's `compare_devices` pattern:
never hardcode a target, always ask the workspace what's actually live.
As of this connection, the workspace has 4 simulator targets and no
provisioned QPU hardware yet:

| Provider | Target | Kind |
|---|---|---|
| quantinuum | `quantinuum.sim.h2-1sc` | syntax checker (structural only) |
| quantinuum | `quantinuum.sim.h2-1e` | emulator (physically realistic noise) |
| pasqal | `pasqal.sim.emu-free` | simulator |
| rigetti | `rigetti.sim.qvm` | simulator |

**Connectivity verified with two real jobs** on `quantinuum.sim.h2-1e`
(`vqe/emulator_smoke_test.py`) — a Bell-state circuit (`H` + `CX`), which
should split ~50/50 between `|00>` and `|11>` with small noise leakage:

| Job ID | `00` | `01` | `10` | `11` |
|---|---|---|---|---|
| `44469383-7d7e-11f1-8269-24ee9a60c281` | 49 | 0 | 0 | 51 |

A clean split confirms the pipeline (auth → connect → submit → real
result) works end to end. **Known gap:** `estimate_cost()` is not
implemented on this backend class (`QuantinuumEmulatorQirBackend`), so
there's currently no pre-submission cost preview before running a job —
worth fixing before submitting the full ~270-circuit EF+ZNE batch to this
emulator.

### Noise-model validation: local model vs. real emulator

The entanglement-forging + ZNE result above (20 → 0.57 kcal/mol) was
built entirely on a *local* `qiskit-aer` depolarizing noise model
standing in for Quantinuum-grade hardware. Once real Quantinuum emulator
access existed, the natural question became: **is that local model
actually a good approximation of real Quantinuum noise?** Four real
validation runs answered this — including one that had to be corrected
after a larger follow-up run reversed it. Every number below is from an
actual submitted job, reported as measured, corrections included.

**1. Diagonal terms, single term (`ZIII`, 500 shots)**
(`vqe/emulator_ef_validate.py`, job `835b03e6-7d81-11f1-b0be-24ee9a60c281`) —
first check, one diagonal Hamiltonian term measured on the alpha
register's leading Schmidt vector:

| | Value | Error vs. exact |
|---|---|---|
| Exact | 0.997900 | — |
| Local Aer (Quantinuum-like) | 0.991875 | 0.006025 |
| Real `quantinuum.sim.h2-1e` | 0.980000 | 0.017900 |

Initial read: the local model *underestimates* real noise by ~3x on this
term. **This did not hold up** — see next result.

**2. Diagonal terms, 10 terms, one job (4000 shots)**
(`vqe/emulator_ef_validate.py`, job `d1564bf2-7d83-11f1-9564-24ee9a60c281`) —
a single Z-basis measurement of the same Schmidt vector at higher shot
count contains the full 4-qubit bitstring distribution, so all 10
largest-weight diagonal terms were derived from *one* real job instead of
paying for 10 separate ones:

| Term | Exact | Local Aer | Real emulator | Local err | Real err |
|---|---|---|---|---|---|
| ZIII | 0.9979 | 0.9919 | 0.9975 | 0.0060 | 0.0004 |
| IIIZ | -0.9979 | -0.9734 | -0.9865 | 0.0244 | 0.0114 |
| IIZI | -0.9979 | -0.9818 | -0.9905 | 0.0161 | 0.0074 |
| IZII | 0.9979 | 0.9878 | 0.9915 | 0.0101 | 0.0064 |
| ZIIZ | -0.9999 | -0.9755 | -0.9860 | 0.0244 | 0.0139 |
| ZIZI | -1.0000 | -0.9839 | -0.9910 | 0.0161 | 0.0090 |
| IZIZ | -1.0000 | -0.9816 | -0.9840 | 0.0184 | 0.0160 |
| ZZII | 0.9999 | 0.9859 | 0.9900 | 0.0140 | 0.0099 |
| IZZI | -0.9999 | -0.9798 | -0.9860 | 0.0201 | 0.0139 |
| IIZZ | 0.9999 | 0.9796 | 0.9850 | 0.0203 | 0.0149 |

Mean local err = 0.0170, mean real err = 0.0103 → **real/local error
ratio = 0.61x**. With a real sample size, the conclusion *reverses*: the
local model actually **overestimates** real Quantinuum-grade noise on
diagonal terms. The single-term "~3x underestimate" from run 1 does not
survive a larger sample and should be read as statistical noise from one
low-shot measurement, not a real effect — left in this README as a record
of the correction, not scrubbed out.

**3. Off-diagonal terms, first attempt (`XXII`, 4 jobs × 1000 shots)**
(`vqe/emulator_offdiag_validate.py`) — off-diagonal (X/Y-containing)
Hamiltonian terms appear in the EF energy as *cross terms*
`<u_n|P|u_m>` between different Schmidt vectors, measured via the
superposition trick from `entanglement_forging_h4.py` (4 phase circuits,
`Re = (E0-E2)/2`, `Im = (E3-E1)/2`). `XXII` was picked as the
largest-*Hamiltonian-coefficient* off-diagonal term — but its actual
cross-term magnitude between the two leading Schmidt vectors turned out
tiny: `|<u_0|XXII|u_1>| ≈ 0.0029`. Real result: `+0.048 +0.004i`, error
0.048961 — **larger than the true signal itself**. Hamiltonian-coefficient
size and cross-term magnitude are not the same thing; this run was
statistically inconclusive, not evidence of anything, and is reported as
such rather than spun into a finding.

**4. Off-diagonal terms, corrected (`IXXI`, 4 jobs × 1000 shots)**
(`vqe/emulator_offdiag_validate.py`, jobs `8f38e187-7d88-11f1-abf7-24ee9a60c281`,
`d9a01da4-7d88-11f1-b4b5-24ee9a60c281`, `238a6d8a-7d89-11f1-a3cd-24ee9a60c281`,
`6df889ca-7d89-11f1-8eaa-24ee9a60c281`) — checked directly which
off-diagonal label actually has the largest cross-term magnitude between
Schmidt vectors 0 and 1: `IXXI`, at `0.961330`, comfortably above the
shot-noise floor:

| | Value | Error vs. exact |
|---|---|---|
| Exact | -0.287634 -0.917290i | — |
| Local Aer (Quantinuum-like) | -0.281079 -0.897009i | 0.021314 |
| Real `quantinuum.sim.h2-1e` | -0.254000 -0.913000i | 0.033906 |

This time the signal is real and resolvable. **Real error is ~1.6x the
local model's error** — the local model *underestimates* real noise here,
the opposite direction from the diagonal-terms finding (0.61x). Plausible
explanation: cross-term circuits are deeper (superposition state prep)
and combine two independently-noisy measurements via subtraction, which a
flat depolarizing model doesn't fully capture — but this is one label,
one pair, not a swept study, so held as a lead rather than a conclusion.

**Bottom line:** the local Quantinuum-like noise model's accuracy is
*not uniform* — it overestimates noise on simple diagonal measurements
(0.61x) and underestimates it on deeper cross-term circuits (~1.6x), on
the evidence gathered so far. The 0.57 kcal/mol EF+ZNE headline result
should be read with that in mind: it's a real, verified simulator result,
but not yet a claim about what real Quantinuum hardware would give for
the *full* EF energy (which mixes many diagonal and cross terms together)
— that full rerun against `quantinuum.sim.h2-1e` is the natural next step
and has not been done.

---

## IonQ Cloud connection

Alongside the IBM and Azure/Quantinuum connections above, this repo also
connects to IonQ's cloud API (`vqe/ionq_backend.py`), authenticated via an
API key from `.env`. It only ever targets `ionq_simulator` — free,
unmetered against QPU credits — never `ionq_qpu` (real paid trapped-ion
hardware), which this repo deliberately does not touch.

`list_backends()` mirrors the same pattern as `compare_devices`/
`compare_targets` on the other two providers: connect, list what's
actually live, never hardcode a target.

### Entanglement forging + ZNE on real IonQ circuits

The entanglement-forging + ZNE result above (20 → 0.57 kcal/mol) was built
on a *local* `qiskit-aer` noise model. `vqe/ionq_run.py` repeats the same
H4-fragment (K=5, all 185 Pauli terms) EF measurement, but via **real
circuits submitted to IonQ's cloud simulator** instead of a local
approximation — every diagonal and cross-term matrix element is measured
by an actual job round-trip, not estimated.

**Sanity check (`noise_model="ideal"`)**: a real cloud round-trip
reconstructing the full 185-term EF energy from scratch lands within
0.63 kcal/mol of the exact numpy answer (0.57 kcal/mol) — confirming the
circuit construction and measurement math are correct at full scale, not
just in the small-molecule cases already verified above.

**Noisy (`noise_model="aria-1"`, IonQ's own emulation of their real Aria
device) + Zero-Noise Extrapolation via gate folding**: IonQ's simulator
doesn't expose a continuous noise-scale knob the way a local Aer noise
model does, so ZNE here uses the standard real-hardware technique
instead — local unitary gate folding (every 2-qubit gate `G` replaced
with `G (G⁻¹G)ᵏ`, physically stretching circuit depth to 1x/3x/5x without
changing the ideal unitary), applied after transpilation so IonQ's own
circuit converter (a straight instruction-by-instruction translation, no
re-optimization) submits it exactly as folded:

| Fold | Hardware energy (Ha) | Error (kcal/mol) |
|---|---|---|
| 1 (no mitigation) | -1.967080 | 125.07 |
| 3 | -1.967587 | 124.75 |
| 5 | -1.964808 | 126.49 |
| ZNE, linear extrapolation | -1.968196 | 124.37 |
| ZNE, quadratic extrapolation | -1.965594 | 126.00 |

**Result: gate-folding ZNE did not work here.** The error is flat across
fold=1/3/5 (125.07 → 124.75 → 126.49 kcal/mol, well within run-to-run
noise), so both linear and quadratic extrapolation land within noise of
the unmitigated value — no improvement, unlike the 35x reduction the same
technique gave against a local depolarizing noise model above.

**Update, later in this repo's history:** this turned out to have a
specific, fixable cause — see "Native-gate folding" further down. Folds
submitted as *abstract* gates were being cancelled by IonQ's compiler
before execution; submitting the same folds as *native* gates
(`GPi`/`GPi2`/`MS`/`ZZ`) recovers a real, monotonic noise-vs-fold
response on both `aria-1` and `forte-1`. The flat result above is real
and correctly reported for what was actually run (abstract gates) — it
just wasn't the final word on whether folding works on this hardware at
all.

### Isolating why: does folding scale noise on IonQ's simulator at all?

The 185-term EF result above is noisy enough (100k shots per matrix
element, summed over many terms) that a real-but-small folding trend could
plausibly hide in shot noise. `vqe/ionq_fold_check.py` isolates the
question cleanly: a 2-qubit Bell state (H, CX), measured in the XX basis
(exact = +1.0 for any fold, since folding a gate with its own inverse
never changes the ideal unitary), folded from 1x up to **81x**, with
1,000,000 shots per circuit — two orders of magnitude more precision per
point than the full run could afford, on a circuit simple enough that a
real trend of any size should be visible.

**Control (`noise_model="ideal"`)**: exactly `+1.0` for all 12 fold
factors (1 through 81) — the folding mechanism itself is confirmed
artifact-free.

**`noise_model="aria-1"`**: `<XX> = 0.985492`, bit-for-bit identical, for
every single fold factor from 1 to 81. Not "flat within noise" — exactly
unchanged to 6 decimal places across a 81x range of physically stretched
circuit depth.

**Conclusion: `aria-1` (and, per the `forte-1` run below, IonQ's other
named simulator noise profile too) applies noise that does not scale with
gate count at all** — consistent with these being fixed, per-shot or
readout-style error profiles rather than literal per-gate
depolarizing-error-rate models. Gate folding is not a usable ZNE knob
against either, for a structural reason rather than a precision problem.
This was checked, not assumed.

### Cross-checking against a second device profile (`forte-1`)

The same K=5, 185-term EF measurement was repeated against
`noise_model="forte-1"` (IonQ's emulation of their Forte device family) to
see whether the flat-with-fold result was specific to `aria-1` or general:

| Fold | Hardware energy (Ha) | Error (kcal/mol) |
|---|---|---|
| 1 (no mitigation) | -1.951852 | 134.62 |
| 3 | -1.951542 | 134.82 |
| 5 | -1.951831 | 134.64 |
| ZNE, linear extrapolation | -1.951757 | 134.68 |
| ZNE, quadratic extrapolation | -1.952232 | 134.38 |

Same pattern: flat across folds, no ZNE improvement. `forte-1`'s baseline
error (134.6 kcal/mol) is noticeably higher than `aria-1`'s (125.1
kcal/mol) — a genuinely different, noisier fixed profile — but it shares
the same non-response to gate folding. The `ideal` control was re-run at
fold=3 and fold=5 too (0.38 and 0.49 kcal/mol, consistent with fold=1's
0.50-0.63 kcal/mol across runs), confirming folding stays artifact-free at
full EF scale, not just in the toy Bell-state check.

**Bottom line:** gate-folding ZNE against IonQ's cloud simulator doesn't
work, confirmed two independent ways (a 12-point, 81x-fold isolation test,
and a second full EF run under a different named noise profile) — not a
tuned result, and not chased into a fix, since the likely cause (fixed
noise profiles instead of per-gate error rates) isn't something a
different circuit or fold scheme would get around.

Run:
```bash
python vqe/ionq_run.py --config <ideal|fold1|fold3|fold5>       # aria-1 (default)
python vqe/ionq_run.py --config fold3 --noise-model forte-1      # any other IonQ noise_model
python vqe/ionq_run.py --assemble --noise-model forte-1          # fit + print ZNE (no network calls)
python vqe/ionq_fold_check.py                                    # the 81x-fold isolation check
```
Each `--config` run is independent and resumable, merging into a
noise_model-namespaced `vqe/ionq_run_results*.json` so different profiles
never overwrite each other.

---

## Multi-fragment molecular tailoring on real IonQ circuits

Everything above runs a *single* fragment (H4) on IonQ. `vqe/ionq_tailoring.py`
extends this to **multi-fragment molecular tailoring** — reconstructing the
H6 chain's ground state (`vqe/covalent_fragment.py`'s 4-atom-block layout)
from three overlapping fragments, each solved by entanglement forging on
real IonQ circuits, then combined via inclusion-exclusion:

```
fragment A = atoms [0,1,2,3]   (8 qubits -> two 4-qubit EF registers)
fragment B = atoms [2,3,4,5]   (8 qubits -> two 4-qubit EF registers)
overlap    = atoms [2,3]       (4 qubits -> two 2-qubit EF registers)
E_tailored = E(A) + E(B) - E(overlap)
```

Classical references (already in the repo, `covalent_fragment_results.json`):
H6 full exact = -3.236066 Ha, H6 tailored (4-atom blocks) = -3.231625 Ha
(2.79 kcal/mol method floor).

### Three circuit-reduction techniques, each independently verified before use

Naively, this would cost ~10,000+ circuits (the full H4 approach, doubled
for two 8-qubit fragments plus the overlap). Three techniques cut that by
roughly 10x — **each checked against exact local math before being
trusted on real hardware**, not assumed correct because they sound
plausible:

1. **Real gauge.** `eigsh`'s ground state is exact only up to an arbitrary
   global phase; rotating the largest-magnitude amplitude onto the real
   axis makes the *whole* state real (residual imaginary part ~1e-14,
   verified every run). Once real, every Pauli label is provably either
   purely real (even Y-count) or purely imaginary (odd Y-count) — never
   both — so each label needs only 2 of the 4 phase-circuit measurements
   (`k=0,2` for this Hamiltonian, where every label happened to be
   even-Y). Verified against the original 4-circuit reconstruction:
   matched to 3e-16 Ha on H4, and to 3e-16 across 222 sampled
   (pair, label) combinations for the phase-circuit shortcut itself.

2. **Beta-register sign reuse — verified per fragment, not assumed.**
   For fragment A (= H4) and fragment B, the beta Schmidt vectors turned
   out to be exactly `v_n = s_n * u_n` for a classically-computable sign
   `s_n` (residual ~1e-13), so `B_nm = s_n*s_m*A_nm` and only the alpha
   register needs real measurement — beta comes free. **This does NOT
   hold for the overlap fragment** (a different, smaller geometry):
   checked, found false (residual ~0.5-1.3, nowhere near zero), and the
   code falls back to measuring both registers independently there
   instead of silently assuming the shortcut generalizes. A real,
   fragment-specific finding, not a bug — exactly the kind of thing
   worth checking rather than trusting.

3. **Qubit-wise commuting Pauli grouping.** Labels that agree on every
   qubit's measurement axis share one circuit (37 alpha labels -> 13
   groups on the 8-qubit fragments, using qiskit's own
   `group_commuting(qubit_wise=True)` — the same grouping
   `BackendEstimatorV2` uses internally, not the riskier general-commuting
   + Clifford-diagonalization grouping, which qiskit has no built-in
   solver for and which would have meant implementing stabilizer-tableau
   math from scratch with no easy way to verify a subtle sign/frame bug).
   Verified: max error 3e-14 across 222 sampled (pair, group-member,
   phase) combinations vs exact Statevector math.

Gate folding is skipped entirely — already proven to produce zero noise
response on this simulator (`ionq_fold_check.py`: a Bell state folded 1x
to 81x gave a bit-identical result every time). Every real run below is
`fold=1`.

### Ideal sanity check: two different, both-honest comparisons

`ideal` was run and confirmed *before* any noisy config, per the same
discipline used throughout this repo — but the "does it match the
classical value" check needed two separate comparisons, not one:

| K | Real vs exact-numpy EF-K (pipeline correctness) | EF-K vs classical tailored (method truncation) |
|---|---|---|
| 3 | 1.03 kcal/mol | 4.49 kcal/mol |
| 5 | 1.30 kcal/mol | **0.17 kcal/mol** |

At K=3, the *pipeline* was already correct (real hardware matched the
same K=3-truncated method's exact-numpy prediction to ~1 kcal/mol, shot-noise
level) — but K=3's Schmidt-rank truncation itself has real, expected error
(4.49 kcal/mol) vs the fully-exact classical tailored value, since fragments
A and B's individual truncation errors don't cancel against the overlap's
(near-zero) truncation error in the inclusion-exclusion sum. K=5 fixes
this on both counts: pipeline-correctness gap stays at shot-noise level
(1.30 kcal/mol) and the method itself lands within 0.17 kcal/mol of the
classical reference — under chemical accuracy, on real hardware, at zero
noise. All noisy runs below use K=5.

### Noisy results (`aria-1`, `forte-1`) and the key question

| Noise model | E(A) err (kcal/mol) | E(B) err (kcal/mol) | overlap err (kcal/mol) | E_tailored err vs classical (kcal/mol) | circuits |
|---|---|---|---|---|---|
| ideal | 0.19 | 0.02 | 0.00 | 0.17 | 746 |
| aria-1 | 84.73 | 84.37 | 1.31 | 167.79 | 746 |
| forte-1 | 88.61 | 88.04 | 1.29 | 175.37 | 746 |

Surviving signal fraction `f` (1.0 = fully preserved, 0.0 = fully
decohered to the mixed-state value) at K=5: ideal ~1.0006 for A/B, ~1.0
for overlap; aria-1 f=0.927 (A), 0.927 (B), 0.997 (overlap); forte-1
f=0.923 (A), 0.924 (B), 0.997 (overlap) — the small overlap fragment
(4 qubits, 96 circuits) is barely touched by noise, while the two
8-qubit fragments (325 circuits each) lose ~7-8% of their signal.

**Answer to the key question: fragmentation error ADDS, it does not
cancel.** If A's and B's errors simply summed (both large, same sign,
since A and B are the identical physical system under the same noise
model) with the overlap's small error subtracted per inclusion-exclusion:
84.73 + 84.37 - (-1.31 contribution accounted for) ≈ 170.4 kcal/mol
predicted — closely matching the observed 167.79 kcal/mol. The overlap
subtraction, being both small and itself just as positively-signed an
error as the others, provides essentially no cancellation benefit here.
Multi-fragment tailoring was **not** more noise-robust than a single big
fragment in this run — if anything, spreading the same noise level across
two independently-measured large fragments instead of one means twice
the large-fragment noise contribution, only partially offset by a small
and equally-noisy overlap term. This is a real, checked negative finding,
consistent with this repo's practice of reporting what didn't work
alongside what did.

Run: `python vqe/ionq_tailoring.py --fragment <A|B|overlap> --noise-model
<ideal|aria-1|forte-1> --K <3|5>` (each fragment is independent and
resumable, refuses to run a noisy config until a matching `ideal` result
exists); `python vqe/ionq_tailoring.py --assemble --noise-model X --K Y`
assembles E_tailored from saved fragments (no network calls). Results in
`vqe/ionq_tailoring_results.json` / `..._aria-1.json` / `..._forte-1.json`.

---

## Native-gate folding: recovering ZNE from a compiler-cancellation bug

The "Gate-folding ZNE did not work here" finding above (`aria-1`/`forte-1`
flat across fold=1/3/5, `<XX>` bit-for-bit identical from fold=1 to
fold=81 in `ionq_fold_check.py`) turned out to have a specific,
addressable cause: circuits submitted as **abstract gates** go through
IonQ's compiler, and qiskit-ionq ships its own optimizer plugin
(`TrappedIonOptimizerPlugin`, with `FuseConsecutiveZZ`/`FuseConsecutiveMS`
passes) that would fuse a folded gate's `G G⁻¹ G` triple straight back
down to `G` if invoked — a real, working mechanism, even though this
repo's own code never calls that plugin locally. **Native-gate
submission** (`gateset="native"`, circuits built directly from IonQ's
`GPi`/`GPi2`/`MS`/`ZZ` instructions) is meant to bypass that kind of
circuit-level recompilation, so `ionq_fold_check.py --native` tests
whether it actually does.

### Result: it does. Native folding recovers real per-gate noise response

| Fold | `aria-1` f | `forte-1` f |
|---|---|---|
| 1 | 0.986496 | 0.985268 |
| 3 | 0.960874 | 0.957220 |
| 5 | 0.935376 | 0.930132 |
| 9 | 0.886640 | 0.878748 |
| 21 | 0.755160 | 0.740006 |
| 81 | **0.338992** | **0.312134** |

Both noise models show a clean, monotonic decay from ~0.985 down to
~0.31-0.34 as fold goes 1→81 — a sharp contrast to the abstract-gate
result's flat 0.985492 at *every single fold*. The `ideal` control stayed
at exactly 1.0 throughout (folding mechanism itself intact — but per the
experiment design, this alone proves nothing either way, since a
cancelled fold also returns 1.0 under noiseless execution; only the
noisy decay is evidence). Two probe circuits were used since they
entangle differently: a `GPI2`+`ZZ`-built probe is a clean `XX`
eigenstate, a bare `MS` gate on `|00⟩` is a clean `ZZ` eigenstate
(verified via Statevector before use, not assumed) — `aria-1` uses the
`MS` probe (its native family), `forte-1` the `ZZ` probe.

**Two real API constraints were discovered by actual submission
failures, not anticipated in advance:**
1. `Gate.inverse()` on `ZZGate`/`MSGate` returns a mis-parametrized
   `"zz_dg"`/`"ms_dg"` gate with the **same, not negated**, parameters —
   not a valid native-gateset instruction and not actually the inverse.
   Folding uses an explicit, numerically-verified inverse construction
   instead.
2. IonQ's real API rejects `ms`/`zz` angle outside `[0, 0.25]` (a genuine
   native-pulse-angle limit, not just a sign convention) — caught on a
   real rejected submission (`"angle" must be less than or equal to
   0.25`), not predicted. `MS` has a phase parameter to exploit instead
   (`MS(φ0+0.5, φ1, θ)` is matrix-identical to `MS(φ0, φ1, -θ)`, verified
   numerically, θ stays in-range). `ZZ` has no such parameter, so its
   inverse uses the conjugation identity `X · ZZ(θ) · X = ZZ(-θ)`
   (verified numerically; `GPI(0)` is exactly the Pauli `X` matrix) — a
   3-gate sandwich at the *same* in-range θ standing in for the single
   out-of-range gate, still exactly `fold` total 2-qubit gates.

### Resource estimate: genuine folding costs 3x the old assumption

Because the old abstract-gate folds were being silently cancelled, any
circuit-count-based resource estimate from that era was really only
paying for fold=1's worth of real gate work, no matter what fold factor
was requested. `vqe/ionq_resource_estimate.py` recomputes this using
genuine (uncancelled) native folding, from numbers already measured
elsewhere in this repo (325 circuits/register at K=5, from
`ionq_tailoring_results.json`'s fragment A; 14 native 2-qubit gates per
state-prep circuit, from `native_stateprep_results.json`, itself verified
on real IonQ, not just locally):

| Fold | 2q gates/circuit | Total 2q gates (1 register) |
|---|---|---|
| 1 | 14 | 4,550 |
| 3 | 42 | 13,650 |
| 5 | 70 | 22,750 |

Summed across fold={1,3,5}: **40,950** real two-qubit gate operations for
one register at K=5 under one noise model — **3.0x** the old
cancelled-fold assumption (13,650, since every fold was secretly costing
the same as fold=1). Both registers, one noise model: 81,900. Both
registers, `aria-1`+`forte-1`: 163,800. This is a resource estimate for
`ionq_simulator` (free); it is not a cost estimate for `ionq_qpu` (real
paid trapped-ion hardware), which this project does not touch and has
not priced out — the point is that any future real-hardware decision
needs to budget for the *real* total, not the old (cancelled-fold) one.

### Native state prep, re-verified: still 14, and it works on real IonQ

The Experiment B result stands (see below) — `native_stateprep.py`'s
hand-derived real-only state prep costs 14 two-qubit native gates, not
fewer than the 11-CX generic baseline. What's new: those circuits were
only checked against local Statevector math before; `python vqe/native_stateprep.py
--real-check` now actually submits a representative sample to the real
`ionq_simulator` (`noise_model="ideal"`) and confirms every emitted
`ms`/`zz` angle is within IonQ's real `[0, 0.25]` limit *and* that the
measured probabilities match `|target|²` to shot-noise-level precision
(8.6e-5 to 5.1e-4 across the sample) — closing the same
"local-math-isn't-enough" gap the fold-check experiment exposed.

Run: `python vqe/ionq_fold_check.py --native` (results in
`vqe/ionq_native_fold_results.json`); `python vqe/ionq_resource_estimate.py`
(results in `vqe/ionq_resource_estimate_results.json`); `python
vqe/native_stateprep.py --real-check` (results merged into
`vqe/native_stateprep_results.json`).

### The real test: does native folding help the FULL forged energy, not just a probe

The probe above is one Schmidt vector through one Z-type observable — a
clean, minimal test of whether folding survives at all. `vqe/ionq_native_forged_energy.py`
measures the real thing: the complete K=5, 185-term H4 forged energy
(same physical system as fragment A above), built entirely from native
gates (state prep via `native_stateprep.py`'s hand-derived real-only
tree, folded via `ionq_fold_check.fold_native_2q`'s verified inverses,
basis-change built abstractly then translated to native via IonQ's own
equivalence library and *appended* — never re-transpiled together with
the already-folded state prep, so the fold is never at risk of being
re-optimized away). `ideal` was run and confirmed (0.222 kcal/mol from
exact) before any noisy config, per this repo's standing rule.

**The per-fold effective error rate is NOT constant here** — unlike the
simple probe's tight ~0.0133–0.0134 consistency across fold 1 through 81:

| | fold=1 | fold=3 | fold=5 | spread |
|---|---|---|---|---|
| `aria-1`, 1−f^(1/L) | 0.157 | 0.142 | 0.137 | 0.0205 |
| `forte-1`, 1−f^(1/L) | 0.161 | 0.149 | 0.145 | 0.0159 |

This is checked, not assumed, and the honest conclusion is exactly what
the check is for: **the extrapolation below should not be trusted as a
clean per-gate-depolarizing fit** the way the probe's was. The full
185-term aggregate — many different circuits, many different basis
rotations, many different 4-qubit states pushed through the 14-gate
native tree — doesn't decay as simply as one Bell-like pair does.

| | `aria-1` | `forte-1` |
|---|---|---|
| Raw (fold=1) | 181.47 kcal/mol | 185.83 kcal/mol |
| ZNE linear | 88.31 kcal/mol | 87.58 kcal/mol |
| ZNE quadratic | 34.25 kcal/mol | 31.82 kcal/mol |
| ZNE exponential | n/a — energies hit the fit's log(0) edge case, reported as unavailable rather than faked | n/a, same reason |

Mitigation genuinely reduces the error (181 → 34 kcal/mol at best,
`aria-1` quadratic) — real, computed from real data — but even the best
case lands **nowhere near** chemical accuracy (1 kcal/mol) or the
classical EF-K=5 method floor (0.5655 kcal/mol), and comes with the
caveat above that the fit's own assumption isn't cleanly satisfied.
Also notable: native fold=1's raw error (181–186 kcal/mol) is
meaningfully *higher* than the earlier abstract-gate fold=1 result
(125.07/134.62 kcal/mol) — the hand-derived 14-gate native state-prep
circuit, while mathematically correct, is less gate-efficient than
whatever qiskit's abstract-to-qis-gateset path produced, so it starts
from a noisier baseline even before folding.

**Resource cost (actual, from this run, not a prediction):** 325
circuits per configuration, 14 native two-qubit gates per fold=1
state-prep circuit, scaling to 14×fold per circuit under folding —
**86,450 real two-qubit gate operations** across all 7 configurations
run (`ideal` + `aria-1`/`forte-1` × fold {1,3,5}), 2,275 real circuits
submitted in total. This project has no verified IonQ `ionq_qpu` pricing
data and does not estimate a dollar figure — that would be fabricating a
number this repo has no basis for, exactly what the honesty rules here
exist to prevent. The real, checkable number is the gate count above.

**Plain English:** no, H4 does not reach chemical accuracy on Forte (or
Aria) with this approach, even after ZNE — the best mitigated result is
roughly 30-35 kcal/mol from exact, about 30-35x too large. It costs
86,450 real two-qubit gate operations to find that out on the free
simulator. Native gates fixed the *compiler-cancellation* problem (the
probe proves gate folding is a real, working noise-scaling knob on this
hardware) — but they exposed a *second*, different problem: the
full forged energy's noise response isn't the simple single-parameter
decay the probe suggested, so folding's benefit here is real but limited,
not the clean fix the probe result on its own would have suggested.

Run: `python vqe/ionq_native_forged_energy.py --noise-model <ideal|aria-1|forte-1>
--fold <1|3|5>` (ideal/fold=1 first, required before any noisy config);
`python vqe/ionq_native_forged_energy.py --assemble` (no network calls —
fits ZNE and prints the consistency check). Results in
`vqe/native_forged_zne_results.json`.

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
- **The entanglement-forging + ZNE result (0.57 kcal/mol) uses exact
  Schmidt vectors, not a variationally trained ansatz**, and the noise
  model is a local `qiskit-aer` approximation of Quantinuum hardware, not
  a real device measurement. Direct validation against the real emulator
  (see Noise-model validation above) found that approximation is *not*
  uniformly accurate — it overestimates noise on diagonal terms (0.61x)
  and underestimates it on off-diagonal cross terms (~1.6x) — so 0.57
  kcal/mol should be read as a real, verified simulator result, not yet a
  claim about real hardware performance on the full EF energy.
- **One validation finding was wrong and got corrected in place, not
  hidden:** an early single-term, 500-shot check suggested the local
  noise model underestimates real noise ~3x; a follow-up with 10 terms
  and 4000 shots reversed that conclusion. Both runs and the correction
  are documented above rather than only keeping the final answer.
- **Gate-folding ZNE on real IonQ circuits did not reduce error**, unlike
  the 35x reduction the same conceptual technique gave against a local
  Aer noise model — confirmed two independent ways, not just one run:
  the full 185-term EF energy was flat across fold=1/3/5 under both
  `aria-1` (125.07 → 124.75 → 126.49 kcal/mol) and `forte-1` (134.62 →
  134.82 → 134.64 kcal/mol), and a dedicated isolation check
  (`ionq_fold_check.py`, a 2-qubit Bell state folded 1x-81x at 1M shots)
  found `aria-1`'s effect on `<XX>` bit-for-bit identical at every single
  fold factor. Likely cause: these are fixed noise profiles, not per-gate
  depolarizing models — reported as measured, not tuned, not chased into
  a fix.
- **Multi-fragment tailoring: fragmentation error adds across fragments,
  it does not cancel**, on real IonQ circuits. The overlap fragment's
  small error (~1.3 kcal/mol) doesn't meaningfully offset the two large
  fragments' errors (~85-89 kcal/mol each, same sign, same physical
  system under the same noise model) — measured, not assumed either way
  going in. Also: the beta-register sign-reuse shortcut (verified for
  the two 8-qubit fragments) was independently re-checked for the smaller
  overlap fragment and found NOT to hold there (residual ~0.5-1.3, not
  ~0) — the code falls back to measuring both registers instead of
  silently trusting a shortcut verified on a different-shaped fragment.

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

# Entanglement forging + ZNE (simulator, Quantinuum-like noise model)
python vqe/entanglement_forging_h4.py # EF on H4 fragment: noiseless + Quantinuum-noise stages
python vqe/entanglement_forging_zne.py # ZNE on top of EF: 20 -> 0.57 kcal/mol

# Real Azure Quantum / Quantinuum (needs AZURE_QUANTUM_* vars in .env)
python vqe/azure_backend.py           # connect + compare live targets (no job submitted)
python vqe/emulator_smoke_test.py     # 1 real job on quantinuum.sim.h2-1e (Bell state)
python vqe/check_job_status.py        # list recent jobs + status in the workspace
python vqe/emulator_ef_validate.py    # local noise model vs real emulator, 10 diagonal terms (1 job)
python vqe/emulator_offdiag_validate.py # local noise model vs real emulator, 1 cross term (4 jobs)

# Real IonQ Cloud (needs IONQ_API_KEY in .env; ionq_simulator only, free)
python vqe/ionq_backend.py            # connect + list live backends (no job submitted)
python vqe/ionq_run.py --config ideal # full 185-term EF energy, real circuits, sanity check
python vqe/ionq_run.py --config fold1 # aria-1 noise, no mitigation
python vqe/ionq_run.py --config fold3 # aria-1 noise, gate-folded 3x
python vqe/ionq_run.py --config fold5 # aria-1 noise, gate-folded 5x
python vqe/ionq_run.py --assemble     # fit + print ZNE extrapolation (no network calls)
python vqe/ionq_run.py --config fold3 --noise-model forte-1  # same, under a different noise profile
python vqe/ionq_fold_check.py         # isolates whether folding scales noise at all (1x-81x)

# Multi-fragment molecular tailoring on real IonQ (H6 chain, 3 fragments)
python vqe/ionq_tailoring.py --fragment A --noise-model ideal --K 5
python vqe/ionq_tailoring.py --fragment B --noise-model ideal --K 5
python vqe/ionq_tailoring.py --fragment overlap --noise-model ideal --K 5
python vqe/ionq_tailoring.py --assemble --noise-model ideal --K 5      # MUST pass before noisy runs
python vqe/ionq_tailoring.py --fragment A --noise-model aria-1 --K 5   # then repeat --fragment for B, overlap

# Native-gate folding: does it fix ZNE? (fold-cancellation isolation, then the real forged energy)
python vqe/ionq_fold_check.py --native        # 2-qubit probe, fold 1-81, native gates
python vqe/native_stateprep.py --real-check   # native state-prep, verified on real IonQ
python vqe/ionq_resource_estimate.py          # gate-count cost of genuine (uncancelled) folding
python vqe/ionq_native_forged_energy.py --noise-model ideal --fold 1     # required first
python vqe/ionq_native_forged_energy.py --noise-model aria-1 --fold 1    # then fold 3, 5, and forte-1
python vqe/ionq_native_forged_energy.py --assemble    # ZNE fits + consistency check, no network calls
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
├── consistency_fragmentation_results.json # Saved results
│
├── entanglement_forging_h4.py         # EF on H4 fragment: noiseless Stage 1 + Quantinuum-noise Stage 2
├── entanglement_forging_h4_results.json # Saved results from entanglement_forging_h4.py
├── entanglement_forging_zne.py        # ZNE on top of EF: 20 -> 0.57 kcal/mol (chemical accuracy)
├── entanglement_forging_zne_results.json # Saved results from entanglement_forging_zne.py
│
├── azure_backend.py                   # Real Azure Quantum workspace connection + compare_targets()
├── emulator_smoke_test.py             # Real job on quantinuum.sim.h2-1e (Bell-state connectivity check)
├── check_job_status.py                # List recent jobs + status in the Azure Quantum workspace
├── emulator_ef_validate.py            # Local noise model vs real emulator: 10 diagonal terms, 1 job
├── emulator_ef_validate_results.json  # Saved results from emulator_ef_validate.py
├── emulator_offdiag_validate.py       # Local noise model vs real emulator: 1 cross term, 4 jobs
├── emulator_offdiag_validate_results.json # Saved results from emulator_offdiag_validate.py
│
├── ionq_backend.py                    # Real IonQ Cloud connection + list_backends() (ionq_simulator only)
├── ionq_run.py                        # Full 185-term H4 EF energy measured on real IonQ circuits + gate-folded ZNE
├── ionq_run_results.json              # Saved results from ionq_run.py, noise_model=aria-1 (per-config, resumable)
├── ionq_run_results_forte-1.json      # Same, noise_model=forte-1 (namespaced so profiles don't overwrite)
├── ionq_run_results_ideal.json        # ideal-noise control at fold=3/5 (folding-artifact check)
├── ionq_fold_check.py                 # Isolates whether gate folding scales noise at all (Bell state, fold 1x-81x)
├── ionq_fold_check_results.json       # Saved results from ionq_fold_check.py (abstract gates)
├── ionq_native_fold_results.json      # Saved results from ionq_fold_check.py --native
│
├── ef_fragment.py                     # EF core generalized to ARBITRARY fragments (real gauge, sign-reuse, grouping)
├── ionq_tailoring.py                  # Multi-fragment H6 molecular tailoring measured on real IonQ circuits
├── ionq_tailoring_results.json        # Saved results, noise_model=ideal
├── ionq_tailoring_results_aria-1.json # Saved results, noise_model=aria-1
├── ionq_tailoring_results_forte-1.json # Saved results, noise_model=forte-1
│
├── native_stateprep.py                # Hand-derived real-only state prep in native gates (14 vs 11-CX baseline)
├── native_stateprep_results.json      # Saved results, incl. real-IonQ verification sample
├── ionq_resource_estimate.py          # Real 2q gate cost of genuine (uncancelled) native folding
├── ionq_resource_estimate_results.json
├── ionq_native_forged_energy.py       # Full K=5 H4 forged energy on real IonQ, native gates, folds 1/3/5
└── native_forged_zne_results.json     # Saved results: energies, f, ZNE fits, consistency check, resource cost

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
