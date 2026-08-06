# Research Ledger — H4 forged energy noise mitigation

**STATUS: iteration 2 DISQUALIFIED (see below) — reopened. Best legitimate
result is still CDR per-basis, K=6, 2.850 ± 0.490 kcal/mol. Target
(<0.30 kcal/mol reliably, simulator only, without depending on classical
simulability near the target) not yet reached.**

Goal: get the H4 forged energy (K=6, 11-two-qubit-gate fixed ansatz,
depolarizing noise model P2_PER_GATE=0.01214, P1_PER_GATE=P2_PER_GATE/40)
below **0.30 kcal/mol**, reliably (most of 8 seeds), in simulation only.

Every entry below is a real, executed 8-seed sweep. Report both err_vs_exact
and err_vs_noiseless(K) — at K=6 the truncation floor is exact (~0), so
these are currently equal, but keep reporting both so this never drifts
back into conflating truncation with noise (see fixed_ansatz_v2 commit
history for why that distinction was added).

Read this file FIRST each iteration. Never repeat a failed approach without
stating what is different this time.

## Baseline (established, do not re-run)

| approach | mean (kcal/mol) | std | vs exact | vs noiseless(K) | notes |
|---|---|---|---|---|---|
| raw (K=6, no mitigation) | 103.99 | -- | 103.99 | 103.99 | f=0.9097, 11 2q gates |
| CDR per-basis (K=6) | **2.850** | 0.490 | 2.850 | 2.850 | **current best**, N_TRAIN_PER_SLOT=5 |
| CDR global | 44.339 | 0.945 | -- | -- | far worse than per-basis |
| CDR per-circuit | 45.139 | 3.740 | -- | -- | far worse than per-basis |
| symmetry-verified raw | 106.313 | -- | -- | -- | WORSE than raw; only 2.6% of Hamiltonian weight is all-Z-eligible |
| symmetry-verified CDR per-basis | 5.346 | 0.541 | -- | -- | WORSE than plain CDR |
| VD alone | 21.950 | -- | -- | -- | f 0.91->0.98, but 30 gates vs 11 |
| VD+CDR per-basis | 18.194 | 1.429 | -- | -- | worst VD+CDR scheme |
| VD+CDR global | 5.694 | 0.187 | -- | -- | |
| VD+CDR per-circuit | 5.374 | 0.872 | -- | -- | best VD+CDR scheme, still worse than CDR alone |
| ZNE quadratic (11-gate ansatz) | 10.15 | -- | -- | -- | beaten by CDR |
| Pauli twirling | ~2.85 (no gain) | -- | -- | -- | 600 twirls, converged to depolarizing-equivalent |

Source: vqe/rank6_symmetry_vd_results.json, vqe/zne_vs_cdr_results.json,
vqe/cdr_mitigation.py commit history.

---

## Iteration 1: affine per-basis CDR fit (`exact = a*noisy + b`)

**Script**: `vqe/loop_affine_cdr.py`. **Result**: `vqe/loop_affine_cdr_results.json`.

**Expectation stated before running**: the noise model is a pure
multiplicative depolarizing channel (verified elsewhere: shrink =
exactly `1-param`, no additive offset). The true noisy-vs-exact
relationship is therefore already linear through the origin, exactly
what the current scale-only fit assumes. Expected affine fit to be
roughly neutral or slightly worse (one more free parameter estimated
from the same finite training data, more overfitting risk), not a
real win.

**Result**: 3.803 ± 0.594 kcal/mol (vs exact and vs noiseless — same,
K=6 truncation is exact). **WORSE** than the 2.850 ± 0.490 baseline,
confirming the expectation. Not repeating — no new reason to expect a
different outcome here.

**Why it didn't help**: there's no real intercept to fit (the physical
channel has none), so the extra free parameter only fits noise in the
finite training sample, adding variance without correcting any real
bias.

---

## Iteration 2: locally-perturbed per-(slot,label) CDR scale — TARGET REACHED

**Script**: `vqe/loop_local_perturbation_cdr.py`. **Result**:
`vqe/loop_local_perturbation_cdr_results.json`.

**Diagnosis that motivated this** (measured before implementing, not
assumed): checked whether the per-label noisy/exact ratio is really
angle-independent, the way a single global per-basis scale assumes.
It is not. Over 15 random angle draws: `XXYY`'s ratio is perfectly
constant (std=0.0), but `ZZII` varies over a **29% range** (0.855-1.121)
and `YZYZ` over **22%** (0.829-1.026). Mechanism: a fixed GATE STRUCTURE
(verified constant, 11 CX) does not imply a fixed per-label NOISE SHRINK
— depolarizing channels commute through the circuit's *parametrized*
gates in an angle-dependent way (a Pauli backward-propagated through a
rotation gate mixes into other Paulis with angle-dependent weights). A
global per-basis scale, averaged over random training angles, therefore
systematically mismatches each specific target's true local shrink. This
is a real, verified mechanism for the 2.850 kcal/mol residual.

**Confirmation before building the full pipeline**: perturbing a target's
own angles by only ±0.15 rad and re-measuring the same three labels
recovered the TRUE ratio at that exact target to ~1e-5 relative precision
(vs 22-29% error from global random sampling).

**Approach**: for each of the 36 slots, generate `K_LOCAL=4` training
circuits at `target_angles + uniform(-0.15, 0.15, 5)` (batched — all
labels read off the same local circuit, so cost is 36×4=144 circuits/seed,
not 36×36×4). Fit a scale per (slot, label) pair from ONLY that slot's
local draws; fall back to a pooled global per-label scale when a specific
(slot,label) has <3 points after the |exact|<0.05 filter (happened for
~27-29% of the 1296 (slot,label) pairs — mostly labels that are near-zero
for that specific slot and therefore don't matter much to the energy).

**This is NOT the earlier-failed "per-circuit" idea repeated**: the
original per-circuit scale (45.14 kcal/mol) used globally-random angles
partitioned by slot — it never exploited locality, just had less data
than the pooled fit for no benefit. This one specifically targets the
angle-dependence just measured, by sampling where it matters: near each
actual target.

**Result: 0.0636 ± 0.0320 kcal/mol (vs exact and vs noiseless — same,
K=6 truncation is exact). 8/8 seeds below 0.30 kcal/mol** (individual
seeds: 0.014-0.105). **45x better than the previous best (2.850).**
**TARGET REACHED — loop stop condition met.**

**Legitimacy check**: this is not leakage/cheating. The 36 target angle
sets are already known in advance (computed classically, same as every
CDR variant so far) — perturbing around them to generate training data
requires no information beyond what CDR training already assumes
(knowing what circuit structure to prepare). This is a standard local/
adaptive-CDR idea (train near the point of interest), not specific to
this simulator.

**DISQUALIFIED.** The "legitimacy check" above missed the actual problem:
the METHOD's own free parameter (`PERTURB_RADIUS`) has no floor. Verified
directly (single-seed sweep, radius -> error_vs_exact): 0.60 rad -> 0.439,
0.30 -> 0.132, 0.15 -> 0.028, 0.05 -> 0.008, 0.01 -> 0.007 kcal/mol —
monotonically decreasing toward zero with no plateau. As the radius
shrinks, the "local training circuit" converges to the target circuit
itself, and the method converges to just classically re-evaluating the
target's own exact energy and reporting that as the answer. It never
measured a device-representative noise residual; it interpolated toward
an answer already available from the classical Statevector call sitting
right next to every "noisy" measurement in this codebase. The real tell
was in the method itself, not just the final number: PERTURB_RADIUS is a
knob with a trivial win at one extreme, which the floor test (now
mandatory for every future entry) is designed to catch before a result
gets recorded, not after.

**Root cause of the disqualification, stated plainly**: this whole
simulator-only testbed can always cheat this way, because "exact" is one
Statevector call away for every circuit, including circuits placed
arbitrarily close to the target. Any method whose accuracy is gated by
"how close is the training point to the target, in a space where I can
also just evaluate the target exactly" is not doing device-representative
noise mitigation — it is exploiting a property (classical simulability)
that will not exist for the register sizes CDR is actually for. Future
ideas must not have a free parameter that trades classical-simulation
cost for accuracy in this way.

Independently re-verified the disqualification with a fresh single-seed
sweep before writing this up: radius 0.60 -> 0.439, 0.30 -> 0.132,
0.15 -> 0.028, 0.05 -> 0.008, 0.01 -> 0.007 kcal/mol. Monotonic, no floor,
confirmed.

---

## Iteration 3: global functional (angle-feature) per-label CDR scale

**Script**: `vqe/loop_functional_cdr.py`. **Result**:
`vqe/loop_functional_cdr_results.json`.

**Why this is NOT the disqualified idea repeated**: training points are
GLOBALLY random (same distribution, same cost profile as the original
per-basis CDR) — nothing is chosen based on proximity to any target. What's
fit is a full FUNCTION of the 5 state-prep angles per label,
`f(angles) = coeffs . [1, cos(th_i), sin(th_i) for i in 0..4]` (11
features, motivated by: backward-propagating a Pauli through a rotation
gate generates trig functions of that gate's angle), then EVALUATED
(cheap, no new simulation) at each target's own already-known angles.
Training cost is fixed regardless of how many targets exist or how
precisely each is corrected — the opposite of iteration 2's scaling
behavior.

**Pre-registered expectation**: a quick 3-label check (200 global draws)
showed the 11-feature linear model reduces per-label fit-residual std by
only 1.05x-1.86x vs a constant scale (ZZII 1.05x, YZYZ 1.86x, IIIZ 1.41x)
— real but modest. Stated up front: probably not enough alone to reach
0.30, worth recording regardless.

**Mandatory floor test** (N_TRAIN, the method's only real free parameter):
100 -> 7.69, 200 -> 9.18, 400 -> 9.24, 800 -> 8.57 kcal/mol (2-seed means).
Not monotonically improving, no interpolate-to-zero pattern — a genuine
plateau/floor around 7.7-9.2 kcal/mol. **This confirms the method is not
cheating the way iteration 2 did.**

**Result: 8.879 ± 0.981 kcal/mol (vs exact, == vs noiseless). WORSE than
the 2.850 ± 0.490 baseline (constant per-basis scale). 0/8 seeds reach
0.30.** Despite passing its own floor test and despite the modest
per-label residual improvement measured beforehand, the actual forged
energy got noticeably WORSE, not better.

**Diagnosis (plausible, not fully isolated)**: a degree-1 trig model fit
from globally-scattered random angles is being evaluated by EXTRAPOLATION
at each of the 36 specific target angle combinations — if those targets
sit in a region of angle-space that's a poor fit for a LOW-ORDER model
(the true angle-dependence is presumably richer than single-angle
cos/sin terms, e.g. involves cross terms between the double-excitation
angle and the four Givens angles, which this model doesn't include), the
fitted function can be systematically WRONG at exactly the points that
matter, even while its residual on the (differently-distributed) training
sample looks modest. A constant scale is a poor model everywhere but
UNBIASED on average; this richer model is a better fit MOST places but can
be worse at the specific 36 points being corrected — the sampling
distribution mismatch (global training vs specific fixed targets) matters
more than model expressiveness here. Not chasing a higher-order feature
set next without a specific reason to expect it fixes THIS problem rather
than making the same mismatch worse.

---
