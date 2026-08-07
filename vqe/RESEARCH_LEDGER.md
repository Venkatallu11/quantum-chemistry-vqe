# Research Ledger — H4 forged energy noise mitigation

**STATUS UPDATE (iteration 9, real hardware): every finding below iteration
9 used this project's OWN synthetic depolarizing noise model. Iteration 9
ran the same raw/CDR/PEC comparison for real, on IonQ's free
`ionq_simulator`, against real `aria-1`/`forte-1` noise, concurrently with
an `ideal` correctness control. Result: a clean, diagnosable negative.
Real target-circuit error (35-43 kcal/mol) is 15-25x larger than a
Clifford-learned channel (from a REDUCED, 1-pair/1-qubit calibration probe)
predicted — the probe under-samples the real noise, not a PEC failure.
CDR, which helped 36x locally, makes things 2.1-2.6x WORSE on real noise
and its binding curve becomes unstable (not just biased). PEC (run at
d=1.0 only) still clearly beats CDR (2.8x) and modestly beats raw
(8-25%), so its correction logic still works — it is just correcting a
channel that was under-characterized for real hardware. See iteration 9
below for the full table and the two honest candidate causes.

**STATUS / GOAL REFRAMED AGAIN (iteration 8): chemistry needs energy
DIFFERENCES (reaction energies, binding curves, barrier heights), not
absolute energies, and chemical accuracy is DEFINED on differences.
Iteration 8 tested whether CDR's/raw's systematic bias (established in
iteration 6) cancels between neighboring geometries of the SAME molecule
under the SAME fixed circuit — it does, for raw (~5x, consistently across
shot levels) and for CDR (grows from ~0.8x at 1e3 shots to ~4.8x at 1e7,
as the shrinking statistical component stops swamping the flat bias) —
**CDR's energy DIFFERENCE crosses chemical accuracy at ~1e6 shots/setting,
even though its ABSOLUTE energy never does at any shot count tested.**
PEC shows the opposite, equally honest pattern: cancellation factor stays
near 1.0 (0.8-1.0x, occasionally slight ANTI-cancellation) at every shot
level, because PEC has little bias left to cancel — its residual is
dominated by independent statistical noise, which does not cancel in a
difference (variances add). PEC's difference error is still the smallest
in absolute terms at every shot level tested, just not because of
cancellation. The binding curve confirms this at the shape level: at 1e5
shots/setting, CDR recovers the equilibrium bond length to 0.004 Å and
the well depth to 1.5 kcal/mol — both well inside chemical accuracy —
despite an absolute per-point error of ~3.3 kcal/mol. See iteration 8
below for the full table.

Prior status (iteration 6, superseded in emphasis but not contradicted --
absolute-energy statements below remain accurate): every result through iteration 5
was shot-noise-free (`density_matrix` estimator, exact expectation
values) — omitting the dominant real-hardware error source entirely.
"Reach 0.30 kcal/mol" was the wrong objective while that omission stood:
at the noiseless-CDR level, shot noise alone would cost 33 million
shots/setting (10.8 billion total) for 0.30 kcal/mol, and IonQ's ~10,000-
shot/job cap makes that unreachable at any budget. Iteration 6 adds a
verified shot-noise model and reframes the deliverable as **the shots-
vs-accuracy trade-off curve for raw / CDR / PEC**, not a single target.
Headline finding: raw and CDR are both **bias-limited** — neither crosses
chemical accuracy (1.0 kcal/mol) at ANY shot count tested up to 10⁷/setting,
because their residual error is systematic (gate-noise bias for raw;
angle-dependent noise-model mismatch for CDR, per iteration 2's diagnosis),
not statistical. **PEC, using iteration 5's Clifford-learned channel with
BOTH calibration and correction shot-limited (not just correction), is
unbiased by construction and DOES converge**: chemical accuracy at
~10⁵ shots/setting, the (now secondary) 0.30 kcal/mol figure at ~10⁶
shots/setting — both real, achievable shot budgets, unlike the naive
10.8-billion-shot estimate for CDR alone. Orbital rotation (iteration 7)
was tried as a classical lever on the shot-noise-driving Hamiltonian norm;
for this specific symmetric, minimal-basis H4 chain it gave only a
marginal reduction (L1 1.005x), reported honestly rather than oversold.

Best surviving-scrutiny **noiseless-estimator** result remains CDR
per-basis, K=6, 2.850 ± 0.490 kcal/mol (bias floor, shot-noise-free).
Best surviving-scrutiny **shot-noise-included, hardware-representative**
result is PEC on the honestly shot-limited learned channel, reaching
chemical accuracy at ~10⁵ shots/setting.

Original goal (kept for context, superseded by iteration 6): get the H4
forged energy (K=6, 11-two-qubit-gate fixed ansatz, depolarizing noise
model P2_PER_GATE=0.01214, P1_PER_GATE=P2_PER_GATE/40) below **0.30
kcal/mol**, reliably (most of 8 seeds), in simulation only.

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

## Iteration 4: gate-by-gate probabilistic error cancellation (PEC) — TARGET REACHED, legitimately

**Script**: `vqe/loop_pec.py`. **Result**: `vqe/loop_pec_results.json`.

**Why this is structurally different from every prior attempt (including
the disqualified iteration 2)**: every CDR variant and iteration 3 corrects
the FINAL, aggregate measured expectation value, using some model of how
noise degrades it — which is exactly what iteration 2's diagnosis showed
is fundamentally limited (a fixed gate structure does not give a fixed
per-label shrink, because of angle-dependent Heisenberg backpropagation
through the circuit's parametrized gates). PEC instead corrects the noise
WHERE IT HAPPENS — gate by gate, during the circuit — using the EXACTLY
KNOWN noise channel (P2_PER_GATE, P1_PER_GATE — known throughout this
project the same way they were used to BUILD the noise model in every
prior experiment, not something new assumed here). This needs **no
training data, no random angles, and no proximity to any target** — there
is no "radius" or "N_TRAIN" parameter to sweep the way iterations 2-3
needed, because there's nothing to fit at all; the correction is derived
analytically from the channel definition.

**Correctness verified BEFORE measuring performance** (two direct tests,
both to machine precision):
1. Reproducing qiskit-aer's own `depolarizing_error` output gate-by-gate
   via an explicit Pauli-mixture formula (`q_I = 1-p(d²-1)/d²`,
   `q_P = p/d²` for `P≠I`) matches Aer exactly: max error `1.67e-16`.
2. Applying that forward channel then its analytically-derived
   quasi-probability inverse (`η_I = 1+p(d²-1)/((1-p)d²)`,
   `η_P = -p/((1-p)d²)`) to a random test density matrix recovers the
   exact input: max error `1.67e-16`.

**Result: E = -2.16638745 Ha, err_vs_exact = 0.000000 kcal/mol (== err_vs_
noiseless, K=6 truncation is exact).** Deterministic (independent re-run
gives 0.00 kcal/mol difference — this method has no randomness, so an
"8-seed sweep" doesn't apply; verified determinism instead of skipping
the requirement). **Both the 0.30 kcal/mol loop target and full chemical
accuracy (1.0 kcal/mol) are reached.**

**Honest sampling-overhead accounting** (the real, non-hidden cost of
PEC on actual shot-based hardware, computed from this circuit's own real
gate counts — 11 CX + 51 single-qubit gates on the u3/cx-transpiled
circuit): per-gate quasi-probability L1 cost γ₂q=1.023 (cx), γ₁q=1.0005
(u3); over the full circuit γ_total=1.315; shot-count multiplier =
γ_total² ≈ **1.73x** vs a noiseless circuit's own shot budget. Modest here
specifically because the per-gate error rate (1.214%) and gate count (62)
are both small — PEC's well-known exponential-in-total-error-rate cost is
real but not yet punishing at this noise level/circuit depth.

**Mandatory floor-test analog** (this method's real free parameter is not
a training knob but the PRECISION of channel characterization — PEC's
exactness assumed EXACT channel knowledge, which real deployments only
approximate): swept the noise-model parameter PEC's inverse ASSUMES,
away from the TRUE injected value, by 0/1/2/5/10/20% relative error.
Result: err_vs_exact = 0.00 / 1.10 / 2.21 / 5.54 / 11.10 / 22.34 kcal/mol
— **error/relative-channel-error ratio is 110.4-111.7 across the whole
sweep, i.e. tightly PROPORTIONAL, not exploding or interpolating toward
zero.** This is the legitimate-degradation signature the disqualified
iteration 2 lacked: accuracy here is bounded by an INDEPENDENT, real
precision requirement (device noise characterization, e.g. the sparse
Pauli-Lindblad learning protocol from the literature), not by how close a
classically-simulated training point is to the target.

**Honest scope limitation, stated plainly**: this iteration used the
noise channel's parameters directly (as GIVEN, same as how they were used
to build the injected noise model itself throughout this whole project) —
it did not implement a separate Clifford-circuit noise-LEARNING step. The
mis-characterization sweep above is the stand-in for that: it shows
PRECISELY how much characterization precision would be needed on real
hardware (e.g., 1% relative error costs ~1.1 kcal/mol) to still clear
chemical accuracy, which is the honest way to report this without
overclaiming a full learn-then-cancel pipeline that wasn't actually built.

---

## Iteration 4b: correction — "TARGET REACHED, legitimately" was overstated

Iteration 4's headline (0.000000 kcal/mol) is real, but the status line
above claimed this legitimately clears the loop target, and that
overstates what was shown. **The 0.000 kcal/mol result holds ONLY at 0%
channel mis-characterization — exact channel knowledge — which no real
device provides.** That condition was documented in the entry (the
mis-characterization sweep exists precisely because of it), but the
STATUS banner didn't carry the condition with the number, which is exactly
the kind of framing the honesty rules now in force (see iteration 5)
exist to prevent.

**Correct framing**: iteration 4 is not an achieved error, it's a
**specification**, read directly off the mis-characterization sweep's own
proportionality (ratio 110.4-111.7 kcal/mol per unit relative channel
error, call it 110.9 as the working constant): to reach 1.0 kcal/mol
(chemical accuracy) the noise channel must be known to **0.90% relative
error**; to reach the loop's 0.30 kcal/mol target it must be known to
**0.27% relative error**. Whether that precision is achievable with a
REAL (Clifford-circuit-only, no classical-simulation-of-general-states)
learning protocol was NOT tested in iteration 4 and is an open question —
answered in iteration 5.

Best surviving-scrutiny result remains **CDR per-basis, K=6, 2.850 ± 0.490
kcal/mol** until iteration 5's answer is in.

---

## Iteration 5: sparse Pauli-Lindblad noise learning (Clifford circuits only) then PEC

**Script**: `vqe/loop_pauli_lindblad_pec.py`. **Result**:
`vqe/loop_pauli_lindblad_pec_results.json`.

**The question, precisely**: iteration 4b reframed PEC's 0.000 kcal/mol as
a specification — chemical accuracy needs the noise channel known to
0.90% relative error, the loop's 0.30 kcal/mol target needs 0.27%. Can a
REAL, Clifford-only learning protocol (van den Berg, Minev, Kandala,
Temme, *Nature Physics* **19**, 1116 (2023)) reach that?

**Protocol implemented**: learn each gate's depolarizing rate from
repeated-application decay on Clifford circuits only — CX (already exactly
Clifford) applied N times (odd N only) to `|+0⟩`, tracking `⟨XX⟩` (ideal
value exactly 1 for every odd N, verified before use, so the fit is a pure
exponential with no oscillation to disentangle); U3 calibrated the same
way using `U3Gate(0,0,0)` (verified its instruction name is literally
`"u3"`, so the noise model attaches to it exactly like the real circuit's
own U3 gates) repeated N times on `|0⟩`, tracking `⟨Z⟩`. **No target
circuit, no target angles, nothing requiring classical simulation of a
generic state was used anywhere in this calibration** — the entire point,
and the thing iteration 2 violated.

**Sparsity-assumption check** (mandatory floor test): the real circuit's
11 CX gates sit on 7 distinct qubit pairs. Learned the rate on each pair
independently rather than assuming uniformity: all 7 gave
`p2_learned=0.01214000`, spread `0.00e+00` — single global rate is
justified BY MEASUREMENT, not by assumption. Same check across all 4
qubits for U3: spread `0.00e+00`.

**Repetition-depth check** (mandatory floor test, and a real correction to
how such checks were framed in iterations 2-3): every depth tested (2 to
13 points) gave IDENTICAL error (`~3e-15`, floating-point level). This is
NOT the "diminishing returns as depth grows" pattern floor tests usually
show — with **exact, noise-free calibration data**, 2 points already fit
an exponential decay exactly, so there is nothing for more depth to
improve. This is the correct, expected signature of noiseless data, not a
red flag — but it also means this particular sweep cannot answer the real
question (how much depth does a SHOT-LIMITED fit need), which is why the
analytic shot-budget calculation below exists.

**Learned-channel result (best case)**: `p2_learned` and `p1_learned`
match the true injected values to `3.3e-15` / `4.0e-14` relative error.
Running PEC with these learned (not true) values on all 36 real targets:
**err_vs_exact = 0.000000 kcal/mol (== err_vs_noiseless)**, deterministic
(verified via independent re-run, diff = 0.00). Matches iteration 4's
110.9-kcal/mol-per-unit-error prediction exactly (predicted ≈ measured ≈
0 at this tiny relative error).

**THIS NUMBER MUST NOT BE READ AS "TARGET REACHED."** It is conditioned on
information unavailable on real hardware: this project's simulator has
**no shot-noise model anywhere**, so this calibration is exact in the same
way every other "noisy" measurement in this whole project has been exact.
The 0.000 kcal/mol here is a best case bounded only by numerical fit
precision, not evidence about what a real, shot-limited device could
achieve.

**The actual, answerable question — analytic shot-noise budget** (standard
error propagation: known variance of a ±1-eigenvalue projective
measurement, propagated through the weighted-least-squares decay fit —
the textbook way such budgets are planned for real experiments; explicitly
NOT a Monte Carlo simulation, since this project has no shot-sampling
machinery to run one):

| shots/circuit | 1-σ relative error |
|---|---|
| 1e2 | 23.70% |
| 1e3 | 7.49% |
| 1e4 | 2.37% |
| 1e5 | 0.75% |
| 1e6 | 0.24% |
| 1e7 | 0.075% |

Solving for the precision bars: **0.90% (chemical accuracy) needs ≈6.9×10⁴
shots per calibration circuit** (≈6.9×10⁶ total across the 100 calibration
circuits used); **0.27% (the 0.30 kcal/mol loop target) needs ≈7.7×10⁵
shots per circuit** (≈7.7×10⁷ total).

**ANSWER to the well-posed question**: both budgets (10⁴-10⁶ shots per
circuit, 10⁷-10⁸ total) sit squarely within the range of real, published
Pauli-Lindblad characterization campaigns (e.g. the original paper ran
comparable or larger budgets on 100+ qubit devices) — **a real
Clifford-only learn-then-cancel pipeline plausibly CAN reach the precision
PEC needs here, at a realistic, not exotic, shot cost.** This is the
honest form of "yes": an analytic estimate with a stated method and a
number, not a simulated proof, and not the exact-simulator's 0.000 kcal/mol
figure misread as a real-hardware result.

**Cost, honestly, not hidden**: 100 distinct calibration circuits (7 CX
pairs × 8 depths + 4 qubits × 11 depths); PEC's own cancellation overhead
on the learned channel is essentially unchanged from iteration 4 (γ_total²
≈ 1.73x), since the learned parameters match the true ones to the
precision this simulator can produce.

**What this iteration validated vs. what it could not**: validated — the
Clifford-only, target-independent PROTOCOL correctly recovers the channel
in structure (sparsity and depth checks both pass legitimately); the
analytic shot-budget calculation gives a real, actionable, favorable
answer. NOT validated — an actual end-to-end run with simulated shot noise
(this project has never built a shot-noise model, in any experiment, so
this is a pre-existing scope limit, not one specific to this iteration).

---

## Iteration 6: shot noise added to the simulator — the goal reframed

**Script**: `vqe/shot_noise_study.py`. **Result**: `vqe/shot_noise_study_results.json`.

**Task 1 — the shot-noise model.** Every prior measurement in this
project used `AerSimulator(method="density_matrix")` with no shot count:
exact Born-rule expectation values, not sampled ones. Replaced with the
exact Binomial sample-mean estimator (`n_plus ~ Binomial(N, (1+e)/2)`,
estimator `2*n_plus/N - 1`) — the true distribution real shot-based
execution produces for a given (possibly noisy) state, applied on top of
the already-exact density-matrix values (computed once, cached, cheaply
re-sampled per shot level/seed — not re-simulating the circuit per trial).

**Mandatory verification, run before Task 2, per instruction:**
- *Convergence*: shots swept 1e3→1e8; error vs the exact result fell from
  3.39 kcal/mol (1e3) to 0.03 kcal/mol (1e8), non-monotonically at
  intermediate points (expected single-trial statistical fluctuation) but
  clearly trending to zero. **PASS.**
- *1/√N scaling*: 100x more shots (1e4→1e6) gave a 10.86x reduction in
  std(E) (200-trial empirical std), vs the 10.00x the √N law predicts.
  **PASS.**
- *L2/√N absolute match*: measured std(E) at 1e5 shots = 0.7375 kcal/mol;
  the naive prediction L2/√N (L2=2.7555 Ha, verified directly from the
  Hamiltonian's own Pauli coefficients, matching the given value exactly)
  gives 5.4679 kcal/mol — **a 7.4x MISMATCH, measured smaller than
  predicted.** Reported honestly, not forced: the standard L2/√N result
  assumes E is a linear combination of independently-measured Pauli terms;
  this project's forged-energy estimator is NOT that — each term's
  contribution is bilinear (`coeff*(diag+cross)/norm2`, with
  `Bmat = S.Amat.S` DERIVED from the SAME measured alpha matrix via the
  beta-reuse shortcut, never independently measured) — a genuinely
  different, more favorable variance structure than the textbook linear
  case. The functional form (1/√N) still holds; the absolute constant does
  not match the simple formula, and that mismatch is itself an honest,
  reportable structural finding, not a bug.

Both mandatory pass conditions (convergence, scaling) passed, so Task 2
proceeded, with the L2 mismatch carried forward as a caveat rather than
gating.

**Task 2 — error vs shots, 8-seed sweep, three methods:**

| n_shots/setting | raw | CDR per-basis | PEC (optimistic cal.) | PEC (honest cal.) |
|---|---|---|---|---|
| 1e3 | 106.90 ± 6.44 | 4.96 ± 2.73 | 4.35 ± 2.43 | 3.84 ± 2.65 |
| 1e4 | 104.18 ± 0.95 | 2.78 ± 1.55 | 1.19 ± 0.89 | 1.47 ± 0.71 |
| 1e5 | 104.06 ± 0.32 | 2.92 ± 0.79 | 0.27 ± 0.13 | 0.48 ± 0.25 |
| 1e6 | 104.05 ± 0.18 | 2.70 ± 0.66 | 0.040 ± 0.028 | 0.132 ± 0.097 |
| 1e7 | 103.99 ± 0.07 | 2.81 ± 0.49 | 0.029 ± 0.020 | 0.082 ± 0.042 |

(kcal/mol vs exact; vs-noiseless is identical throughout, K=6 truncation
is exact.) "Optimistic" PEC reuses iteration 5's exact, infinite-shot
calibration and only shot-limits the target correction — an overstatement
of real performance. **"Honest" PEC shot-limits the CALIBRATION too, at
the same shot count as the target correction, re-learning p2/p1 from
noisy calibration data at every shot level** — this is where iteration
5's `3.3e-15` learned-channel error becomes a real, shot-count-dependent
number: measured calibration `p2` relative error was 3.08% at 1e3 shots,
0.91% at 1e4, 0.35% at 1e5, 0.043% at 1e6 — matching iteration 5's
*analytic* prediction (~0.75% at 1e5, ~6.9e4 shots/circuit needed for
0.90%) closely, from an actual (not analytic) shot-sampled re-run.

**Crossing points** (mean over 8 seeds first drops below the bar):

| method | chemical accuracy (1.0) | 0.30 kcal/mol |
|---|---|---|
| raw | never (tested to 1e7) | never |
| CDR per-basis | never (tested to 1e7) | never |
| PEC, optimistic calibration | 1e5 shots/setting | 1e5 shots/setting |
| PEC, honest calibration | **1e5 shots/setting** | **1e6 shots/setting** |

**The actual finding**: raw and CDR are **bias-limited**, not
statistics-limited — their error is flat (raw: pinned at ~104 kcal/mol;
CDR: pinned at ~2.7-3.0 kcal/mol) across four orders of magnitude of
shots, because the residual is systematic (gate-noise bias for raw;
iteration 2's angle-dependent noise-model mismatch for CDR) — more shots
cannot fix a bias. **PEC is unbiased by construction** (iteration 4/5),
so it genuinely converges with more shots, crossing chemical accuracy at
a real, achievable ~1e5 shots/setting even with fully honest (shot-
limited, Clifford-only) calibration. This directly answers why "reach
0.30 kcal/mol" was the wrong framing for CDR alone (no amount of shots
gets there) while giving PEC a concrete, favorable, hardware-realistic
number instead.

---

## Iteration 7: orbital rotation to shrink the Hamiltonian coefficient norm

**Script**: `vqe/orbital_rotation_study.py`. **Result**:
`vqe/orbital_rotation_study_results.json`.

**Approach**: parametrized an orthogonal 4x4 rotation (6 independent
Givens angles, via matrix exponential of an antisymmetric generator —
guarantees orthogonality by construction) applied to the RHF MO
coefficients post hoc, recomputed the one/two-electron integrals and the
mapped qubit Hamiltonian in the rotated basis, and minimized the
resulting L1 norm with `scipy.optimize.minimize` (Nelder-Mead; budget
capped at 3 restarts x 150 evaluations, ~0.45s/evaluation, ~200s/restart
— a LIMITED search, stated plainly, not an exhaustive global optimization).

**Physics-invariance check** (mandatory before trusting any rotated-basis
number): recomputed the exact ground-state energy in the optimized
rotated basis and compared to the untouched-basis value — **diff =
9.47e-12 kcal/mol**, i.e. exact to the solver's own numerical precision,
confirming the rotation is a pure basis change with zero physics impact,
as any orthogonal orbital rotation must be.

**Result: L1 = 9.7175 Ha (from 9.7694, a 1.005x reduction), L2 = 2.7555 Ha
(unchanged, 1.000x).** A genuinely modest, close-to-negative finding,
reported as measured rather than reframed as a win. Plausible reason: RHF
orbitals for this specific highly-symmetric H4 chain in a minimal STO-3G
basis (only 4 spatial orbitals, no room for the kind of localization gains
seen in larger/less-symmetric systems in the literature) are already close
to whatever basis a simple rotation search finds — combined with the
limited optimization budget (150 evals/restart is not exhaustive over a
6-parameter nonlinear objective), this result should be read as "orbital
rotation did not help much HERE," not "orbital rotation cannot help" in
general.

**Implied shot-budget impact** (variance ~ norm²/N, so shots for fixed
precision ~ norm²): 1.011x fewer shots via L1, ~1.000x via L2 — negligible.
Explicitly not re-verified against the actual forged-energy pipeline (that
would need re-deriving the whole Schmidt decomposition, fixed-ansatz
angle-fits, and CDR/PEC pipeline in the rotated basis, a substantial
follow-up not attempted here) — reported as a specification derived from
the norm reduction alone, consistently with the honesty rules, not as a
re-measured shot count.

---

## Iteration 8: energy DIFFERENCES, not absolute energies — does the bias cancel?

**Script**: `vqe/energy_difference_study.py`. **Result**:
`vqe/energy_difference_study_results.json` (plus per-shot-level
checkpoints `energy_difference_study_partial_*.json` — the full 5-level x
8-seed x 7-geometry sweep does not fit in one Bash-tool command even
backgrounded, 10-minute hard cap, so it runs as 5 separate `--shots N`
invocations checkpointed to disk, combined by `--assemble`).

**The insight tested** (established in iteration 6, not re-derived): CDR's
residual is BIAS, flat across shots (2.70-2.92 kcal/mol, 1e4-1e7). A bias
similar at two nearby geometries should cancel in their difference — and
chemistry runs on differences (reaction energies, binding curves, barrier
heights), which is where chemical accuracy is actually defined.

**Setup**: H4 chain at d = 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0 Å, the SAME
fixed 11-gate ansatz, K=6, the SAME noise model, the SAME 5 shot levels as
iteration 6. Re-verified (not assumed) at every geometry: Schmidt rank
stays ≤6 (exact K=6 truncation holds everywhere tested) and the 11-gate
count stays fixed. Verified the alpha-label set is IDENTICAL across all 7
geometries before relying on it to reuse CDR training and PEC calibration
across geometries (both are properties of the fixed CIRCUIT/gate noise,
not the target Hamiltonian) — a real efficiency win, not assumed.

**Cancellation factor (mean|abs error| / mean|diff error|, d_ref=1.0 Å),
8-seed means, every shot level:**

| n_shots/setting | raw | CDR | PEC (honest) |
|---|---|---|---|
| 1e3 | 4.52x | 0.80x | 1.01x |
| 1e4 | 5.13x | 1.12x | 0.82x |
| 1e5 | 5.17x | 2.99x | 1.02x |
| 1e6 | 5.12x | 4.42x | 0.81x |
| 1e7 | 5.12x | 4.81x | 0.85x |

**Absolute vs difference error (kcal/mol, mean over 8 seeds x 6 non-ref
geometries), at 1e6 shots/setting**: raw 98.9 abs / 19.3 diff; CDR 3.3 abs
/ **0.75 diff**; PEC 0.11 abs / 0.14 diff.

**Raw**: real, consistent cancellation (~5x) at every shot level — its
error is entirely a large, shot-noise-independent bias, so the bias
dominates the total error at any shot count tested, giving stable
cancellation.

**CDR: cancellation GROWS with shots (0.80x → 4.81x)** — a real, physically
sensible pattern, not noise: at low shots, CDR's error is a MIX of
(non-cancelling) statistical noise and (cancelling) bias, with statistics
dominating; as shots grow, the statistical part shrinks as 1/√N while the
bias stays flat, so bias comes to dominate and cancellation strengthens.
**Consequence: CDR's energy DIFFERENCE crosses chemical accuracy (1.0
kcal/mol) at ~1e6 shots/setting (0.747 kcal/mol) — real, even though CDR's
ABSOLUTE energy never crosses chemical accuracy at any shot count tested
in iteration 6 or here.** This is the reframing working exactly as
hypothesized, for CDR specifically.

**PEC: no reliable cancellation (0.8-1.0x, sometimes just below 1 —
mild ANTI-cancellation)**, and this is equally honest, not a failure to
find something that should be there: PEC is close to unbiased by
construction (iteration 4/5), so there is little systematic bias left TO
cancel — its residual is dominated by independent statistical noise at
each geometry, and differencing two INDEPENDENT noisy quantities of
similar size increases the combined variance rather than cancelling it
(variances add for independent measurements). **PEC's difference error is
still the smallest of the three at every shot level tested (e.g. 0.058
kcal/mol at 1e7, vs CDR's 0.682) — just not because of cancellation.
Different mechanism, still the best method.**

**Comparison to `vqe/difference_cancellation_results.json`** (pre-dates
this ledger, found fragment errors ADD across different molecules with
different circuits, no cancellation): that was the least favorable case
for cancellation (different molecules, different circuit structures).
This is the most favorable case in principle (same molecule, same fixed
11-gate circuit, only target angles differ) — and the cancellation factor
here came out real and substantial for the bias-dominated methods (raw,
high-shot CDR), confirming the mechanism the earlier study's negative
result did not rule out. It does NOT hold for PEC, and that is reported
plainly too, not glossed over.

**Binding curve shape, exact vs noisy (1e5 shots/setting, local quadratic
fit around the true minimum at d=0.9 Å — a real methodological fix made
here: an all-7-point fit spanning the anharmonic dissociation tail out to
2.0 Å gave a nonsensical d_eq near -4.4 Å on the first attempt; the fit is
correctly restricted to the 4 points [0.8, 0.9, 1.0, 1.1] Å bracketing the
actual minimum, the standard way to extract equilibrium geometry from a
sampled curve):**

| method | d_eq (Å) | d_eq error (Å) | well depth (kcal/mol) | well depth error (kcal/mol) |
|---|---|---|---|---|
| exact | 0.9001 | — | 176.146 | — |
| raw | 0.9875 | 0.087 | 147.841 | 28.31 |
| CDR | 0.9043 | **0.004** | 174.669 | **1.48** |
| PEC (honest) | 0.8989 | 0.001 | 176.131 | 0.02 |

CDR recovers the equilibrium bond length to 0.004 Å and the well depth to
1.5 kcal/mol — both comfortably inside chemical accuracy — despite a
~3.3 kcal/mol absolute error at every individual point. The curve SHAPE
survives even where the absolute energies do not, exactly matching the
difference-cancellation finding at the level of a full property (not just
one geometry pair).

**Mandatory floor-test note**: the free parameters here (shot level,
which geometry pair) were SWEPT, not tuned to a favorable outcome — the
full 5-level table is reported for both raw and CDR, including the low-
shot regime where CDR's cancellation factor is BELOW 1 (0.80x at 1e3
shots) and PEC's is also below 1 at several levels. No cherry-picking:
the pattern (CDR cancellation growing with shots, PEC staying near 1) is
consistent and monotonic-in-shots for CDR, which is itself evidence this
is a real effect and not noise in a single measurement.

## Iteration 9, Task 1: real IonQ QPU pricing — is the $25.79 floor per-circuit or per-job?

**THE ANSWER FIRST, per the request that prompted this task: real QPU
hardware is not an option at any bundling strategy, and this makes the
per-job-vs-per-circuit question moot rather than decisive.** Queried
IonQ's own real, free, read-only `GET /jobs/estimate` endpoint (no
hardware touched or reserved) for `qpu.forte-1` — the only backend name
of the ones tried (`qpu.forte-1`, `qpu.aria-1`, `aria-1`, `qpu.aria-2`,
`qpu.harmony`) that returned a quote; aria-1 pricing was simply
unavailable through this account/endpoint, not fabricated.

Real rate card: `job_cost_minimum=$25.7899`, `cost_1q_gate=$0.000164`,
`cost_2q_gate=$0.001121` (unit: gates). Three real quotes:

| job | gates (1q/2q) | shots | real quoted price |
|---|---|---|---|
| 1 circuit | 51/11 | 1 | $25.79 (floor dominates) |
| 1 circuit | 51/11 | 10,000 | **$206.95** (8.02x above the floor) |
| 125 circuits' worth of gates, merged into 1 job | 6375/1375 | 10,000 | $25,868.75 (**exactly 125.0x** the 1-circuit price) |

The 125x scaling test is exact and decisive on its own narrow question:
**the floor is charged once per JOB**, confirmed by direct measurement,
not assumed. But that finding is secondary here, because gate-execution
cost already exceeds the floor by 8x for a SINGLE circuit at this
project's actual 10,000-shots/setting — bundling more circuits into one
job never brings the floor back into play; it was never binding to begin
with at this shot count.

**The real number that matters**: this project's actual planned
real-hardware workload (36 K=6 targets x 13 qubit-wise-commuting groups x
7 geometries = 3,276 circuits/noise-model, at $206.95/circuit) costs
**$677,968 for ONE noise model, $1,355,936 for aria-1+forte-1 together**
— **452x the $3,000 award** — and that is a LOWER bound (excludes CDR
training and PEC calibration circuits entirely). No bundling strategy
changes this conclusion; gate-execution cost, not the per-job floor, is
what makes real hardware unaffordable here. Confirms the free
`ionq_simulator` (real submission, zero cost, per IonQ's own docs) is
the only viable path for Task 2 below — exactly what was already
specified, now backed by a real, queried number rather than an assumption.

Full data: `vqe/ionq_resource_estimate_results.json` (`real_pricing_check` key).
Code: `vqe/ionq_resource_estimate.py::real_pricing_check()`.

## Iteration 9, Task 2: does PEC's advantage survive noise it did not design? Real submission to IonQ's free `ionq_simulator`, concurrent ideal/aria-1/forte-1

**Every result in this ledger through iteration 8 used this project's OWN
depolarizing noise model** (`P2_PER_GATE=0.01214`, `P1_PER_GATE`
=`P2_PER_GATE/40`) — PEC's near-exactness (iterations 4-5) is close to
tautological against a channel built to be exactly the kind of channel
PEC inverts. This task ran the same raw/CDR/PEC comparison for real,
submitted to IonQ's free cloud simulator (`ionq_simulator`, zero cost —
never `ionq_qpu`), against `aria-1` and `forte-1`'s own real noise
models, with `ideal` as a third, concurrently-submitted correctness
control. **Headline: it is a clean negative, exactly the kind the
question anticipated — not because PEC breaks, but because the specific
noise-LEARNING probe used here badly under-estimates the real error on
the actual target circuits, and CDR turns out to make things
substantially worse than doing nothing.**

**Scope, reduced from iteration 8's design and disclosed here, not
hidden** (real network round-trips, not local computation, are now the
bottleneck): 3 geometries (0.9/1.0/1.1 Å, bracketing iteration 8's own
d_eq≈0.90 Å), not 7. PEC's own randomized-circuit protocol (see below) is
run only at d=1.0, not all 3 — its difference-error/cancellation-
factor/binding-curve fields are therefore correctly N/A, not missing
data. CDR training used 8 seeds × 5 random-angle draws/seed (genuinely
independent real submissions), not iteration 8's local per-slot scheme.
PEC calibration used 1 representative CX pair + 1 qubit, not every
distinct pair/qubit as iteration 5 verified locally — **this specific
reduction turns out to be the main story below, not a footnote**. The
raw/CDR "8 seeds" at the target-measurement step are bootstrap resamples
(multinomial resampling of the real integer counts from ONE real
10,000-shot execution per circuit) — stated once, applies throughout;
CDR training and PEC's quasi-probability circuit draws are genuinely
independent real executions, not resamples.

**Concurrency, verified by wall-clock time, not asserted**: every phase
submitted every job (`backend.run()`, non-blocking) before calling
`.result()` on any of them.

| phase | jobs | circuits | submit time | retrieve time |
|---|---|---|---|---|
| calibrate | 4 | 1,040 | 22.2s | 447.6s |
| targets d=0.9 | 3 | 1,404 | 12.6s | 522.1s |
| targets d=1.0 | 3 | 1,404 | 12.3s | 542.1s |
| targets d=1.1 | 3 | 1,404 | 12.7s | 512.4s |
| pec | 16 | 7,488 | 99.6s | 813.5s |

The `pec` phase submitted 16 jobs (7,488 circuits, 5.3x the circuit count
of one `targets` phase) but its retrieval time was only 1.5x longer
(813.5s vs ~525s) — sub-linear scaling in circuit count is the expected
signature of genuine concurrent server-side execution, not proof by
itself, but consistent with it and inconsistent with the jobs having run
one at a time.

**Ideal is a correctness control, checked immediately, not glossed
over**: raw energy on `noise_model="ideal"` must reproduce the noiseless
numpy energy within real 10,000-shot statistical noise, or the script
raises and stops (a pipeline bug, not a noise finding). All 3 geometries
passed: d=0.9 err=0.056 kcal/mol, d=1.0 err=2.689 kcal/mol, d=1.1
err=1.862 kcal/mol vs the noiseless energy — all consistent with real
shot noise at this shot count, none indicating a bug.

**Learned channels (real, Clifford-only, per noise model, never shared,
never from the local model)**:

| model | p2 (2-qubit) | p1 (1-qubit) | γ_total | fit residual (CX / ry) |
|---|---|---|---|---|
| aria-1 | 0.000173 | 0.000000 | 1.0036 | 0.0010 / 0.0000 |
| forte-1 | 0.000284 | 0.000000 | 1.0059 | 0.0013 / 0.0000 |

Both are **40-70x smaller** than this project's own local model
(`P2_PER_GATE=0.01214`, γ_total=1.315) and the fit residuals are tiny —
the exponential-decay fit itself is clean, not noisy or curved. Read
naively, this predicts almost no correction is needed on real IonQ
noise. **That prediction is wrong**, and the reason why is the real
finding here.

**The actual result table** (8-seed mean, `chemical_accuracy_kcal=1.0`):

| model | raw abs err (kcal/mol) | CDR abs err (kcal/mol) | PEC abs err (kcal/mol) |
|---|---|---|---|
| ideal (control) | 1.704 | — | — |
| aria-1 | 34.983 | **89.751** | 32.132 (d=1.0 only) |
| forte-1 | 43.026 | **90.257** | 32.422 (d=1.0 only) |

**Finding 1 — the calibration/target mismatch is the headline number**:
raw error on real noise (35-43 kcal/mol) is **15-25x larger** than both
the ideal-control baseline (1.7 kcal/mol, pure shot noise) and what the
tiny learned γ_total≈1.004-1.006 would predict. The Clifford CX/ry-decay
probe — deliberately reduced here to 1 representative pair and 1
representative qubit, unlike iteration 5's local verification (spread
<1e-6 across all distinct pairs, justifying a single global rate) — does
not generalize to the real 11-CX/51-1q target circuits. Two honest
candidate causes, not adjudicated between here: (a) the reduction itself
was unjustified for real hardware — other qubit pairs/qubits may carry
real error the single-pair probe never sampled, unlike the local
synthetic model where uniformity was independently verified; (b) IonQ's
real per-gate error is genuinely context-dependent (crosstalk, connectivity,
coherent/non-Pauli effects) in a way an isolated two-qubit Bell-decay
circuit cannot see, even if that one pair's own isolated error truly is
tiny. Both are real possibilities; distinguishing them needs the
all-pairs/all-qubits calibration iteration 5 ran locally, not done here
for real-network-cost reasons — the honest scope limit of this run, not
a claim resolved by it.

**Finding 2 — CDR makes it WORSE, a genuine reversal from every prior
iteration in this ledger**: CDR's abs error (89.8-90.3 kcal/mol) is
**2.1-2.6x raw**, not an improvement. Locally, CDR helped by 36x (2.850
vs raw's 103.99). On real IonQ noise it actively hurts. This is
consistent with — and sharpens — iteration 2/3's own original diagnosis:
CDR's per-basis linear scale assumes a fixed per-label attenuation, but a
Pauli's noisy attenuation depends on how it Heisenberg-propagates
backward through the circuit's PARAMETRIZED gates, which differs between
CDR's random training angles and the real target angles. That mismatch
was already known to cap CDR's local performance; on real hardware noise
that is evidently less uniform than this project's synthetic depolarizing
channel, the same mismatch is bad enough to overshoot in the wrong
direction rather than merely under-correct.

**Finding 3 — PEC gives a real but modest edge, consistent with
Finding 1's diagnosis**: PEC (d=1.0 only) reaches 32.1-32.4 kcal/mol,
beating raw by 8-25% and CDR by ~2.8x — a genuine, not cherry-picked,
improvement, but nowhere close to iterations 4-5's near-zero local
result. This is exactly what Finding 1 predicts: PEC is correcting for
the LEARNED channel (tiny, from the 1-pair/1-qubit probe), and if the
real noise affecting the full circuit is larger or differently
structured than that channel, PEC under-corrects rather than failing
outright — an honest partial result, not a null one.

**Finding 4 — CDR's binding curve is unstable, not just biased**: local
quadratic fit (3-point window, the reduced geometry set):

| model / method | d_eq (Å) | error vs exact (0.8539 Å) |
|---|---|---|
| ideal / raw | 0.9025 | 0.049 |
| aria-1 / raw | 0.8680 | 0.014 |
| forte-1 / raw | 0.8912 | 0.037 |
| aria-1 / CDR | 0.6394 | 0.215 |
| forte-1 / CDR | 0.0465 | 0.807 |

Raw's binding-curve shape survives reasonably (d_eq errors 0.01-0.05 Å,
comparable to the ideal control's own 0.05 Å shot-noise floor) even
though its absolute-energy error is large — echoing iteration 8's
cancellation finding, now confirmed on real hardware noise for the
UNCORRECTED signal. CDR's binding curve does NOT survive (errors 0.2-0.8
Å) — its per-geometry correction is erratic enough, not just biased
enough, that the 3-point quadratic fit is unstable. This is the opposite
of iteration 8's local finding (CDR's absolute bias was large but SMOOTH
across geometries, so it canceled in differences and gave a clean
binding curve) — on real IonQ noise, CDR's bias is not smooth enough
across geometries for that cancellation to hold.

**A real reproducibility bug caught and fixed during this run**: the
bootstrap-resample RNG seeds initially used Python's built-in `hash()` on
`(model, d, seed)` tuples — `hash()` on tuples containing strings is
randomized per-process (`PYTHONHASHSEED`) in Python 3, so re-running
`--assemble` on the SAME real checkpointed data gave different numbers
each time (caught by literally running `--assemble` twice and diffing).
Fixed with a `zlib.crc32`-based deterministic seed (`stable_seed()` in
`vqe/ionq_simulator_binding_curve.py`); confirmed identical output across
repeated `--assemble` runs before reporting the numbers above.

**Answering the task's question directly: does PEC keep its advantage on
real IonQ noise?** Partially, and for a diagnosable reason, not a mysterious
one. PEC still clearly beats CDR (2.8x) and modestly beats raw (8-25%),
so its DIRECTION of advantage over CDR survives intact — CDR's collapse is the
sharper story here. But PEC's MAGNITUDE of advantage over raw shrinks from
"eliminates the error" (iterations 4-5, exact/near-exact locally) to "a
modest dent" (this run) — consistent with the channel-learning probe,
not PEC's correction logic itself, being the bottleneck: gate-by-gate PEC
is only as good as the channel it inverts, and this run's deliberately
reduced 1-pair/1-qubit Clifford probe evidently does not capture the real
noise affecting the full 62-gate target circuit. **A full-coverage
Clifford calibration (every distinct CX pair, every qubit — iteration 5's
local protocol, not yet run for real) is the natural next real-hardware
experiment this result points to, not a re-run of what was done here.**

Code: `vqe/ionq_simulator_binding_curve.py` (phases: `--calibrate`,
`--targets --d D`, `--pec`, `--assemble`, each independently checkpointed
under `vqe/ionq_simulator_binding_curve_checkpoints/` since real network
round-trips exceed the 10-minute-per-command budget this project has
worked within since iteration 8). Full data:
`vqe/ionq_simulator_binding_curve_results.json`.

---
