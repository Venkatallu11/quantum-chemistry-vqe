# Quantum Chemistry VQE

[![Qiskit Ecosystem](https://qisk.it/e-e52f2069)](https://qisk.it/e)

**A from-scratch quantum chemistry engine, and an honest account of
where it breaks.** No PySCF, no Psi4, no lookup tables for the core math
— Gaussian integrals, Hartree-Fock, and the qubit Hamiltonian are all
built from atomic geometry alone, then verified against PySCF and run on
IonQ's Cloud service. Every number below is a real, reproducible
computation. Nothing is hardcoded, and nothing that failed has been
removed.

### The short version

In October 2026 we put an 8-qubit forged H4 calculation on a live
trapped-ion machine — a ~$3,882 job on IonQ Research Credits, with the
analysis frozen and committed before the data existed. One submission,
no retries.

**Our noise model predicted the raw hardware error to within 1%** —
23.84 kcal/mol measured against 24.06 predicted. **And then our error
correction missed by 4.16×**, returning 3.46 kcal/mol where it had
promised 0.83. Chemical accuracy is 1.0.

The fit quality tells you why: χ²/dof was ~2 on every simulator and
**5.93 on hardware**, in all 20 bootstrap draws. So we tagged all 732
residuals by measurement group and Pauli weight and found the two things
the model gets wrong — a correction whose error grows with its own
magnitude, and a systematic sign offset in one measurement group.

A number that passes a threshold is worth less than knowing exactly which
term in your error model is lying to you. [Jump to the hardware
result](#the-h4-real-hardware-result-job-01a138e7-16-oct-2026).

> If you are here to check whether the headline claims survive scrutiny:
> §1 opens with a caveat an IonQ reviewer raised about this repository's
> own code, which we adopted and which demoted our best number. The
> project's record of every negative result is in `vqe/RESEARCH_LEDGER.md`.

**Terminology, stated precisely up front:** almost every result in this
README runs on IonQ's free `ionq_simulator`, using noise models *named*
`aria-1` and `forte-1` — historical device noise profiles available for
simulation. **These are not live hardware submissions.** Aria-1 hardware
itself has been retired; it is not possible to submit to it, real or
otherwise. The only genuine real-hardware runs in this entire project
targeted `qpu.forte-enterprise-1` (the live Forte system) directly, are
explicitly labeled as such with real job IDs, and are called out
individually below — everywhere else, "aria-1"/"forte-1" means the
simulator's noise model of that name, not a physical machine.

---

## The chain

```
H4 (chemical accuracy)  →  H6 covalent chain (167.79 → 0.0094 kcal/mol)
        →  LiH  →  formaldehyde  →  hydrogen peroxide  →  acetaldehyde
                (four real carcinogens/ROS molecules, none built for)
```

Each link re-uses the *exact same* validated measurement method on a
harder or more different molecule than the last — real IonQ free
simulators for the chain, one funded live-hardware run at the H4 link,
every limitation disclosed below.

**Project status: the H4 real-hardware run is done, and it did not
reach chemical accuracy on the headline estimator.** On 16 Oct 2026 the
full 84-circuit H4 design ran on live `qpu.forte-enterprise-1` (job
`01a138e7`, estimated $3,882.21). The frame-free number — the honest one — came in
at **3.46 kcal/mol**, against 0.83 predicted by the forte-1 noise model.
The full result, including what went wrong and where, is in §1.
Everything past H4 is still simulator-only.

---

## 1. Chemical accuracy, on IonQ cloud simulators (noise-modeled, not hardware)

H4 (4 electrons, 8 qubits) split via **entanglement forging** into two
independently-measured 4-qubit registers. Four techniques stacked to get
there:

- **Ancilla-parity leakage detection** — postselect out shots where noise
  pushed the state out of its symmetry sector.
- **Conditioned PEC + GPi2 correction** — a noise-inversion model fit
  *without ever looking at the exact energy* (selected by training
  chi²/dof only), then validated out-of-sample.
- **Joint Schmidt-frame fit** — one shared 6×6 orthogonal frame fit
  across all 21 measurement slots at once, instead of 21 independent fits.
- **General-commuting measurement grouping** — exact minimum graph
  coloring cuts 273 circuits to 84 (69% fewer), verified to 6.66e-16
  before trusting it; confirmed on real data paired with the frame fit
  at 0.0228/0.0204 kcal/mol (see §1 caveat section below).

**Result — real submissions to IonQ's `ionq_simulator`, using the
`aria-1` and `forte-1` noise models (simulator-only, not live hardware
— see terminology note above), replicated across 4 independent
submissions each (8 numbers total):**

**0.0105 – 0.0192 kcal/mol**, every single run.

*(chemical accuracy = 1 kcal/mol — this is 50-95x under the bar.)*

**The full real-hardware run has since been done, and the simulator
numbers above did not transfer.** See
"[The H4 real-hardware result](#the-h4-real-hardware-result-job-01a138e7-16-oct-2026)"
below for the measured hardware numbers. Two earlier real-hardware jobs
on `qpu.forte-enterprise-1` preceded it: a circuit-correctness spot-check
(job `01a08910-7a2b-762b-b0ad-6207191241b6`, real billed cost matching a
from-scratch cost model to within $0.05) and the 10 Sep GC3 job whose
three circuits are reused in the full design.

### An important caveat, found by external review — read this before trusting the number above

An IonQ reviewer read this repository's own code and raised a real,
correct point: the 0.0105-0.0192 kcal/mol number above
comes from a **joint Schmidt-frame fit** whose reconstruction basis is
built directly from the exact, classically pre-computed FCI Schmidt
vectors — and the circuits' own state-prep angles are separately fit to
reproduce those same exact vectors. In plain terms: **this is not yet a
blind variational calculation.** The target answer is supplied to the
circuit construction and to the reconstruction basis; what's being
measured is whether real, noisy hardware (plus error mitigation) can
reproduce a state whose identity is already known — a genuine and
useful test of the noise-mitigation stack, but a different and weaker
claim than "found an unknown ground state from scratch."

**Removing the frame fit** — keeping only the physics-based correction
stack (ancilla-parity postselection + conditioned PEC + GPi2 correction),
with no reference anywhere to the exact Schmidt vectors in reconstruction
— gives a real, honest, harder number. On any single real 20,000-shot
draw this is **inconsistent**: aria-1 ranged 1.50-2.24 kcal/mol across 3
independent real draws (never under the 1 kcal/mol bar), forte-1 ranged
0.07-1.86 kcal/mol (passed on 1 of 3 draws). Diagnosed, not guessed: this
turned out to be mostly ordinary shot noise, not a systematic bias.

Two independent fixes were tried, both real, both using only the 3
already-collected real draws (zero new spending):

| Approach | aria-1 | forte-1 | Method |
|---|---|---|---|
| Simple pooling | 0.3630 kcal/mol | 0.8645 kcal/mol | Combine raw counts across draws, same reconstruction as always (`task60_h4_no_frame_fit_multidraw.py`) |
| **Covariance-aware GLS** | **0.2007 kcal/mol** | **0.4640 kcal/mol** | Per-slot physical-state fit via generalized least squares, properly accounting for real correlations between Pauli expectations measured from the same accepted bitstrings (`task61_h4_no_frame_gls.py`) |

**GLS is the better result** — a genuinely more correct statistical
treatment of the same real data (simple pooling wrongly treats
correlated measurements as independent) — and it's the one we're
reporting. **One honest caveat, not hidden**: the GLS fit's own
chi²/dof is 28-150 (a well-calibrated fit should be near 1) — a real,
quantified sign that the assumed noise-correction parameters don't
fully describe the real noise. The number is real and the method is
statistically sound; the underlying correction model is disclosed as
imperfect, exactly the same limitation the frame-fit result was already
resting on, now visible because nothing is smoothing over it. A
follow-up attempt to recalibrate those parameters via held-out
cross-validation (never touching the exact energy) completed but came
out worse (aria-1=0.4263, forte-1=0.4103, self-reported as failing both
the 0.25 and 0.02 kcal/mol bars) — a real negative result, not adopted.

Neither approach reaches 0.25 kcal/mol on both backends at once — stated
plainly, not stretched. What both do show: real chemical accuracy
(<1 kcal/mol) is achievable without any oracle information in
reconstruction, at a real, honest, disclosed cost in tightness compared
to the frame-fit number.

A separate attempt to also reduce circuit count via general-commuting
(GC) grouping (4 groups instead of 13, 84 circuits instead of 273) was
tried *without* the frame fit and made things *worse* (ideal-backend
error jumped to 0.88 kcal/mol) — a real, understood, disclosed finding:
denser measurement groups dilute per-label shot statistics, trading
circuit-count efficiency for statistical precision. See
`vqe/task59_h4_gc_no_frame_fit.py`.

**But pairing that same low-circuit GC design *with* the frame fit works
— and it's a real result, not a projection.** `vqe/task47_gc_noisy_champion_comparison.py`
first showed this analytically (GC+frame-fit at ~2x QWC+frame-fit's
error, both far under the 0.25 kcal/mol bar). `vqe/task66_h4_gc_frame_fit.py`
then confirmed it with real finite-shot data — reusing the exact same
already-collected real GC circuit counts from task59 above (zero new
submissions, zero new cost), fed through the joint Schmidt-frame fit
instead of the no-frame reconstruction:

| Backend | GC (84 circuits) + frame-fit | Same data, no frame fit (task59) |
|---|---|---|
| ideal | 0.0039 kcal/mol | 0.8784 kcal/mol |
| aria-1 | **0.0228 kcal/mol** | 9.5660 kcal/mol |
| forte-1 | **0.0204 kcal/mol** | 8.2435 kcal/mol |

Real chemical accuracy, in the same range as the original 273-circuit
result (0.0105-0.0192 kcal/mol), using **3.25x fewer circuits** — the
joint frame fit's cross-slot pooling is what absorbs the shot-noise
dilution that broke the no-frame-fit reconstruction on the same data.
This still carries the same oracle-dependency caveat as every frame-fit
result above (the reconstruction basis is built from the exact FCI
Schmidt vectors) — it is a real circuit-efficiency result, not a fix for
the no-frame-fit question.

### Pre-run checks for the H4 hardware benchmark (tasks 67–78)

**Status: these were the pre-run checks; the hardware job has since run —
see the next subsection for what actually happened.** Everything in *this*
subsection ran on IonQ's free `ionq_simulator` (noise models named
`aria-1`/`forte-1`). Commit messages and file names for tasks 72–77
sometimes say "real trials"/"real data": there, "real" means *fresh
simulator submissions* (not reused or bootstrapped data), never hardware.

**Plan, as submitted:** the full 84-circuit design (4 GC groups × 21
slots) at 1,100 shots on `qpu.forte-enterprise-1`; 81 new circuits plus 3
reused from the 10 Sep hardware job. Cost $3,875 (IonQ's own estimate
endpoint, read around an SDK field-name bug) to $3,882 (rate-card
formula), within the $5,352.16 approved. IonQ's predicted queue time was
~14.5 days; the job actually sat ~46 hours.

**Fresh simulator submissions at 1,100 shots, shared-frame error vs exact
(kcal/mol):**

| Design | ideal | aria-1 | forte-1 | Under 1.0 |
|---|---|---|---|---|
| 54 circuits, 4 trials | 0.21–0.61 | 0.18–0.67 | 0.17–1.75 | 11/12 |
| 84 circuits, 2 trials | 0.18, 0.61 | 0.49, 0.45 | 0.22, 0.64 | 6/6 |

Without the shared frame (corrected per slot only), the 84-circuit trials
give 0.90–2.05 kcal/mol on aria-1/forte-1 (1 of 4 under 1.0): improved,
but not chemical accuracy.

**Fit quality (chi²/dof, shot-noise weighted):** ideal 0.9–1.2; aria-1 and
forte-1 1.6–2.3 per trial, 2.4–2.5 pooled. On the noisy models, the
residual scatter is about twice what shot noise alone predicts, so the
corrected data carries systematic error the noise model does not
describe, and the shared frame is absorbing it. The bootstrap rehearsals
on older data (1.1–1.5) understated this. Before task78, these result
files stored a unit-weight mean squared residual (~0.0002–0.002) under
the name `chi2_dof`; that value is now `unit_weight_msr`.

**Two sets of task numbers 72–74 exist.** Both are in the repo and both
are current:

| Hardware-run pipeline | Fresh simulator validation |
|---|---|
| `task72_rehearsal_sign_and_chi2_check.py`: Schmidt-sign pin and real chi² | `task72_h4_full_experiment.py`: 54-circuit plan, 4 fresh simulator trials |
| `task73_hw_analysis.py`: locked analysis for the hardware counts | `task73_h4_pooled_experiment.py`: the 4 trials pooled (4,400 shots) |
| `task74_hw_submission.py`: builds, checks and submits the hardware job (dry run by default) | `task74_real_hardware_cost_check.py`: IonQ cost estimates |
| | `task75`–`task77`: shot sweep, shot reallocation, 84-circuit fresh trials |
| `task78_chi2_relabel_backfill.py`: chi² relabel for the right-hand column | |

### The H4 real-hardware result (job `01a138e7`, 16 Oct 2026)

The 84-circuit design ran on live `qpu.forte-enterprise-1`. One
submission, no retries, analysis frozen and committed before the data
existed. Here is what came back.

| | raw | **frame-free** | shared frame | χ²/dof |
|---|---|---|---|---|
| ideal noise model | 1.96 | 0.048 | 0.117 | ~1 |
| aria-1 noise model | 22.53 | 1.041 | 0.707 | ~2 |
| forte-1 noise model | 24.06 | 0.832 | 0.394 | ~2 |
| **live hardware** | **23.84** | **3.46** | **0.86** | **5.93** |

*(kcal/mol error against exact FCI. Chemical accuracy = 1.0.)*

**The noise model predicted the uncorrected hardware error to within 1%.**
23.84 measured against 24.06 predicted — a 0.99× ratio. Shot retention
after ancilla-parity postselection came in at 0.894 against 0.898
predicted. For an 8-qubit forged calculation on a live trapped-ion
machine, two months of noise-model work landed almost exactly on the
nose.

**And that is precisely what makes the next line credible: the correction
stack failed.** The frame-free estimator — the honest one, with no exact
Schmidt vectors anywhere in reconstruction — returned **3.46 kcal/mol**
against 0.83 predicted. A 4.16× miss. It still removed 85% of the raw
error (23.84 → 3.46), so the physics is not wrong; it is incomplete.

**The shared-frame number is 0.86, and we are not calling that a pass.**
Its bootstrap distribution (N=20, every circuit resampled at its own shot
total, whole pipeline re-run) centres at **1.055 ± 0.419**, with 13 of 20
draws above 1.0 and q95 at 1.866. The point estimate landing under the
bar is the luck of one draw, not a result that would replicate. Reported
because it was pre-registered, not because it passed.

**χ²/dof tripled, and that is the actual finding.** ~2 on every noise
model, **5.93** on hardware, and every one of the 20 bootstrap draws
agrees (6.41–7.22). A tripled χ² alongside a 1%-accurate raw prediction
says something specific: the hardware carries error structure the
ZZ/GPi/GPi2 correction does not describe at all. Not more noise —
*different* noise.

#### Where the missing physics is

Rather than leave that as a shrug, `task82_hw_chi2_residual_breakdown.py`
re-derives the identical 732-term residual vector (it asserts χ²/dof
matches the locked 5.9286 exactly) and tags every term with its slot, GC
group and Pauli weight. Mean r² per cell:

| | weight 1 | weight 2 | weight 3 | weight 4 |
|---|---|---|---|---|
| GC group 0 | 3.59 | **10.14** | — | — |
| GC group 1 | — | 2.80 | — | 4.70 |
| GC group 2 | — | — | 3.98 | **13.58** |
| GC group 3 | — | 3.35 | — | 4.77 |

Two distinct failure modes, not one:

1. **The error grows with the size of the correction applied.** Binned by
   correction ratio, labels the model barely touches (ratio 1.00–1.05)
   give mean r² = 3.69; labels it scales hardest (1.15–1.25) give **8.51,
   and carry 48% of the total χ²** on 33% of the terms. The correction's
   *magnitude* is miscalibrated, and the miscalibration compounds exactly
   where the model works hardest. The twelve worst residuals in the run
   are all weight-4 XY-type strings — `YYXX`, `YXXY`, `XYYX`, `XXYY`,
   `YYYY`, `XXXX` — every one with a correction ratio near 1.19.
2. **GC group 0 carries a systematic offset, not scatter.** Its mean
   *signed* residual is −1.32 (weight 1) and −1.73 (weight 2). Noise
   averages to zero; this does not. Something in group 0's correction is
   biased in one direction.

That is a concrete, falsifiable diagnosis of a specific error model on a
specific machine — which is worth more than a number that happened to
land under a threshold.

#### Conditions, as agreed with IonQ in advance

`debiasing=False` passed explicitly on every circuit (with debiasing on,
the identical job costs ~$13,600); 1,100 shots per circuit; a single
submission with no resubmission permitted; Job ID posted to the IonQ
ticket on submission. Submitted 2026-10-14 05:34 UTC, completed
2026-10-16 03:36 UTC. Cost $3,882.21 by IonQ's own estimate endpoint and
our independent rate-card model, against the $5,352.16 approved under the
IonQ Research Credits Program; the final billed figure is not yet
reported back through the API (`billed_cost_usd` is still null) and will
be updated here when it is.

**What this run is, stated plainly:** a hardware-versus-noise-model
comparison in which the model won on the raw observable and lost on the
corrected one. It is not a demonstration of chemical accuracy on real
hardware. The frame-free estimator missed the bar by 3.5×, the number is
reported anyway, and the honest conclusion is that the correction model
needs the two repairs the residual breakdown just localised.

Reproduce it: `python vqe/task73_hw_analysis.py --counts vqe/task73_hw_counts.json --design full`
(raw counts for all 81 circuits are committed), then
`python vqe/task82_hw_chi2_residual_breakdown.py --counts vqe/task73_hw_counts.json`.

---

## 2. A real negative result, resolved

Chaining three overlapping fragments to reconstruct a full H6 covalent
molecule (`E_tailored = E(A) + E(B) - E(overlap)`) using the *old*
entanglement-forging + ZNE method gave a genuine, honestly-reported
failure: **fragmentation error adds, it does not cancel.**

| | Real error, old method |
|---|---|
| aria-1 | 167.79 kcal/mol |
| forte-1 | 175.37 kcal/mol |

Re-running the *same three fragments* through the new ancilla-parity +
conditioned-PEC + joint-frame pipeline instead:

| | Real error, new method |
|---|---|
| aria-1 | **0.0094 kcal/mol** |
| forte-1 | **0.0025 kcal/mol** |

A **~20,000-70,000x improvement** over the old method's real result.

**Honest scope:** the remaining ~2.79 kcal/mol vs. the *fully exact* H6
energy is the classical tailoring method's own truncation floor — not a
quantum measurement limitation. Closing it needs bigger classical
fragments, not more quantum work (this repo's own H8 data already shows
the direction: 4-atom blocks floor at 6.72 kcal/mol, 6-atom blocks at
1.15 kcal/mol).

---

## 3. Does it generalize? Four molecules it was never built for

The validated pipeline was never designed around H4 specifically. Four
genuinely different real molecules tested it:

| Molecule | What it is | Real error (ideal / aria-1 / forte-1) |
|---|---|---|
| **LiH** | Heteronuclear, real p-orbitals — not an H-chain | 0.018 / 0.053 / 0.092 kcal/mol |
| **Formaldehyde (CH2O)** | A real IARC Group 1 confirmed human carcinogen | 0.0036 / 0.0246 / 0.0303 kcal/mol |
| **Hydrogen peroxide (H2O2)** | The real reactive-oxygen-species molecule behind oxidative DNA damage | 0.0022 / 0.1351 / 0.0975 kcal/mol |
| **Acetaldehyde (CH3CHO)** | A real IARC Group 1 carcinogen, ethanol's toxic metabolite; a genuinely harder (6e,6o), weight-3 register shape, generic (non-optimized) ansatz | 0.0125 / 0.2122 / 0.1230 kcal/mol |

All twelve numbers land **inside chemical accuracy** — on real
molecules, independently cross-checked against a separate exact reference
computation before any quantum measurement (zero deviation in all four
cases).

**Said plainly, no overclaiming:** this is an electronic-structure
computation, not drug discovery or a claim about biological activity.
Real drug molecules (10-30+ heavy atoms) are categorically out of reach —
exact diagonalization needs a matrix of size 2^(2×orbitals), and no
computer, this one included, can hold that for a real chemotherapy
molecule. LiH and formaldehyde used real, standard, disclosed
approximations (a deliberate Schmidt-rank truncation for LiH; a real
literature-precedented (4e,4o) active space for formaldehyde and H2O2) —
never hidden, always stated.

**One open, unexplained finding, left in rather than smoothed over:**
H2O2's real noisy-backend error (0.098-0.135 kcal/mol) is meaningfully
worse than formaldehyde's (0.025-0.030 kcal/mol) — despite matching gate
count and nearly identical retention fractions. Something about H2O2's
specific Hamiltonian is more noise-sensitive. Not yet understood.

---

## How it works

1. **`chem.py`** — a pure numpy/scipy integral engine (McMurchie-Davidson),
   s and p orbitals, STO-3G. No external quantum chemistry package.
2. **Hartree-Fock → Jordan-Wigner → exact diagonalization**, verified
   against PySCF at every molecule.
3. **Entanglement forging** — split the register in half, Schmidt-decompose
   the exact ground state, measure each half's Pauli expectations
   separately, reconstruct the energy classically.
4. **Real IonQ Cloud** — free `ionq_simulator` (with `ideal`/`aria-1`/
   `forte-1` noise models — simulator-only) for every result above;
   real `qpu.forte-enterprise-1` trapped-ion hardware for the few
   explicitly-labeled spot-checks; every real dollar spent disclosed
   with cost estimates given in advance.

---

## Honest limits

- Closed-shell molecules only — an open-shell (odd-electron) fragment was
  tested and found to hit a genuine ground-state degeneracy the current
  pipeline doesn't resolve (real, verified, not yet fixed).
- Exact diagonalization tops out around 16 qubits on a laptop; bigger
  systems need fragmentation or an active-space approximation, both used
  and disclosed above.
- LiH did **not** compress the way H4 does — its true Schmidt rank is 15
  (the full theoretical maximum for its register), not a small number.
  Bigger molecules don't automatically get cheaper.
- The H4 headline result was replicated across 4 independent real
  submissions per backend (8 numbers total) specifically because a
  single lucky run isn't evidence — every one of those 8 real numbers is
  reported above, not just the best one.
- **The correction stack does not reach chemical accuracy on live
  hardware.** Measured, not estimated: 3.46 kcal/mol frame-free on
  `qpu.forte-enterprise-1`, against 0.83 predicted by the forte-1 noise
  model. The two localised causes are in §1; neither is fixed yet.
- **Simulator noise models are optimistic about corrected observables,
  even when they are accurate about raw ones.** The same model that
  predicted the raw hardware error to 1% was wrong by 4.16× after
  correction. Treat any corrected-observable projection from a noise
  model as an upper bound on performance, not an estimate.
- **Bootstrap resampling runs optimistic relative to fresh submissions.**
  Bootstrap predicted a frame-free maximum of 0.90 kcal/mol and a 100%
  pass rate; fresh simulator submissions of the same design returned up
  to 2.46. Cross-submission drift exceeds pure shot noise, so bootstrap
  error bars understate reality.
- **At 1,100 shots per circuit the frame-free estimator is shot-noise
  limited, not noise limited.** The zero-noise ideal simulator itself
  returned 1.225 kcal/mol on one fresh 1,100-shot draw. No amount of
  error mitigation fixes that; only more shots or a lower-variance
  estimator does, and 6× the shots costs ~$23,400.

---

## Run it yourself

```bash
# --- the real hardware run (raw counts for all 81 circuits are committed) ---
python vqe/task73_hw_analysis.py --counts vqe/task73_hw_counts.json --design full
                                                                # job 01a138e7: raw 23.84, frame-free 3.46, shared-frame 0.86, chi2/dof 5.93
python vqe/task82_hw_chi2_residual_breakdown.py --counts vqe/task73_hw_counts.json
                                                                # localises the missing physics: group 2 weight-4, and a signed offset in group 0

# --- simulator results ---
python vqe/task40_robustness_envelope_new_pipeline.py         # H4 headline: Q95 ~0.003-0.05 kcal/mol
python vqe/task40_certification_ablation_adversarial.py       # H4 WITHOUT the frame fit, single draw
python vqe/task60_h4_no_frame_fit_multidraw.py                # H4 without the frame fit, simple pooling -> 0.36/0.86 kcal/mol
python vqe/task61_h4_no_frame_gls.py --backend aria-1 --p-gpi2 0.0006 --checkpoint vqe/task39c_ancilla_real_submission_draw1.partial.json,vqe/task39c_ancilla_real_submission_draw2.partial.json,vqe/task39c_ancilla_real_submission_draw3.partial.json
                                                                # H4 without the frame fit, covariance-aware GLS (better) -> 0.20/0.46 kcal/mol
python vqe/task66_h4_gc_frame_fit.py                            # H4, low-circuit GC design (84) + frame fit, real data -> 0.02 kcal/mol
python vqe/task52_covalent_tailoring_result.py                # the real tailored H6 result, no network calls
python vqe/task54_lih_new_molecule.py --analyze                # LiH, already-collected real data
python vqe/task55_formaldehyde_carcinogen.py --analyze         # formaldehyde, already-collected real data
python vqe/task56_hydrogen_peroxide.py --analyze                # H2O2, already-collected real data
python vqe/task57_acetaldehyde_carcinogen.py --analyze          # acetaldehyde, already-collected real data
```

Full iteration-by-iteration history (60 iterations, every bug found and
fixed, every negative result kept in) is preserved in git history and in
`vqe/RESEARCH_LEDGER.md` for anyone who wants the whole story.

---

## Credits

Built alongside the open-source **Quantum Hardware MCP server** by
Lokesh Pullakandam (https://github.com/Lokesh-2025/quantum-hardware-mcp),
which connects AI assistants to live IBM Quantum hardware.

## License

MIT — Venkata Rao Allu. See [LICENSE](LICENSE).
