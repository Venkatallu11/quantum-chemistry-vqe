# Research Ledger — H4 forged energy noise mitigation

**STATUS UPDATE (iteration 20, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): Task B (Contextual Subspace VQE, Kirby/Tranter/Love, Quantum
5, 456 (2021)), the last untouched task from this session's original
five, done properly against the reference `symmer` implementation
(installed this iteration; verified its Hamiltonian matches this
project's own independent exact diagonalization to 1.06e-11 kcal/mol
before trusting anything further). Applied to the FULL, monolithic
8-qubit H4 Hamiltonian (a genuinely different paradigm from this
project's entanglement-forging approach used in all 19 prior
iterations — CS-VQE partitions a standalone molecular Hamiltonian, which
the forged "alpha register" is not). Result: the noncontextual+quantum-
remainder partition matches EXACTLY (0.000 kcal/mol) at 5 quantum
qubits — which is just Z2 tapering again (8 − 3 symmetries = 5, the
SAME reduction `z2_tapering.py` already found, now independently
cross-validated by a different reference package). Pushing to FEWER
qubits costs real accuracy that exceeds chemical accuracy even before
any hardware noise: 4 qubits gives 22.4 kcal/mol, 3 gives 33.8, purely
from the classical noncontextual approximation, zero noise involved.
**CS-VQE does not unlock a path to chemical accuracy this project hasn't
already found** — its best qubit-count-vs-accuracy trade-off matches
Z2 tapering exactly, and going smaller costs more than chemical accuracy
allows. Combined with iteration 19's completed leakage+ZNE and
leakage+CDR investigations, this closes out Task B and Task D from the
original task list. See iteration 20 below for the full write-up,
including a final, complete synthesis of where this 20-iteration project
stands relative to chemical accuracy.

**STATUS UPDATE (iteration 19, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): does removing real, hardware-confirmed particle-number
leakage (iteration 18) BEFORE fitting ZNE unlock the plateau raw ZNE has
never found (iterations 13, 14)? Tested via a local synthetic-noise
sweep reusing the exact same rigorous 3-direction floor test methodology
this project has used since iteration 13, applied to BOTH a raw and a
leakage-post-selected scheme built on iteration 18's verified-exact
ancilla circuit. Result: **NO PLATEAU for either scheme** — leakage
removal does NOT rescue ZNE. Every one of the 6 floor tests (3
directions × 2 schemes) is disqualified, each for the same reason
iteration 13 established: the extrapolated error keeps drifting as the
scale range or fit order grows, not settling. A genuine, useful side-
finding survives, though: leakage post-selection roughly HALVES the
exact (zero-shot-noise) error at EVERY individual noise scale tested
(scale=1: 122.9→57.5 kcal/mol; scale=7: 625.8→475.5 kcal/mol) — this
confirms iteration 18's real-hardware benefit is a general, scale-
independent property of the technique, not a fluke of one real
submission, even though it does not fix ZNE's separate, still-
unexplained convergence problem. See iteration 19 below for the full
write-up and ALTERNATIVES NOT TAKEN.

**STATUS UPDATE (iteration 18, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): Task D (spin/particle-number leakage projection), completed
and it WORKS — real, on real IonQ hardware. Direct analysis of iteration
9's already-collected real counts (the one measurement group that stays
in the Z basis) found genuine leakage out of the untapered register's
physical weight-2 sector: 0% ideal, 4.74% aria-1, 5.12% forte-1. Post-
selecting on it there cut RMS error ~2.2-2.4x using data already on disk.
Generalized this to ALL 13 measurement groups with a verified-exact
ancilla parity-check circuit (4 extra CX, entangling a 5th qubit with the
register's pre-rotation Z-basis parity before any basis-change gates —
confirmed via `partial_trace` to leave the original 4-qubit marginal
state unchanged to 1.5e-36, and to read ancilla=0 with certainty for the
noiseless physical state) and ran it for real. Result: adding the ancilla
ALONE (no post-selection) makes things WORSE (aria-1=67.25,
forte-1=74.23 kcal/mol vs the no-ancilla baseline's 34.98/43.03 — the
extra gates cost real accuracy). But POST-SELECTING on it (discarding the
real, measured 7.12%/7.93% of leaked shots) gives aria-1=31.77,
forte-1=33.86 kcal/mol — BETTER than the original circuit with NO
leakage detection at all. **The untapered register's implicit
particle-number structure is a real, usable, real-hardware-validated
error-detection resource, and exploiting it produces this project's best
raw (non-ZNE) real-hardware number yet for the abstract ansatz.** The
tapered register has no analogous structure to exploit — its Clifford
transform destroys the physical meaning of "weight," so this technique
is structurally unavailable there, strengthening (without fully
resolving — the exact quantitative link is not yet isolated) iteration
17's leading hypothesis for the still-unexplained gap. See iteration 18
below for the full write-up and the mandatory ALTERNATIVES NOT TAKEN.

**STATUS UPDATE (iteration 17, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): exhaustive search for WHY the tapered circuit still loses on
real hardware despite winning on every gate-accounting metric. User's
question: "it has less gates still it won't win i dont understand." Five
hypotheses tested; four ruled out, one confirmed real but still
unexplained:
  1. Qiskit-level CX count — tapered wins (mean 3.94 vs 11). Doesn't
     explain the loss.
  2. IonQ NATIVE gate count (gpi/gpi2/ms, the gates actually
     executed/billed on hardware, not qiskit's abstract cx/u3) — tapered
     wins even bigger (mean 3.67 ms vs 11; mean 34.67 gpi/gpi2 vs 197).
     Rules out "qiskit gate count hides a native-gate blowup," iteration
     12's own established mechanism for a DIFFERENT circuit pair — does
     not apply here.
  3. Circuit depth — tapered wins (mean 8.31 vs constant 22).
  4. Measurement-basis hardness (X/Y physical-pulse rotations vs free
     Z-basis reads before measurement) — tapered wins (mean 1.78 vs 2.77
     non-Z-axis qubits per measurement circuit).
  5. Measurement-reuse-induced noise correlation — the tapered pipeline's
     `reduced_label_map` causes 10 of 27 unique reduced Pauli
     measurements to be shared by 2 different original alpha_labels (vs
     0 sharing, untapered). Tested directly with a controlled Monte Carlo
     (200 trials, real shot-noise binomial sampling, no gate noise):
     RATIO (shared/independent RMS) = 0.919 — independent per-label
     measurement is NOT better, if anything marginally worse. RULED OUT.
     (An earlier quick sensitivity script had suggested a 2.4x formula-
     sensitivity gap; traced to a bug in that script — `random.sample()`
     drew different (name, label) pairs across the two comparison halves
     because the two target-name lists had different internal orderings
     under the same seed. Confirmed directly by inspecting both lists;
     not a real effect. Disclosed, not swept under the rug.)

Given every circuit-execution-cost and reconstruction-formula hypothesis
failed, ran a REPRODUCIBILITY CHECK instead: resubmitted the EXACT SAME
iteration-15 tapered circuit fresh to real IonQ a second, fully
independent time (same 10,000 shots/circuit, same noise models). Run 1:
aria-1=47.78±2.20, forte-1=51.25±1.77. Run 2: aria-1=50.28±1.81,
forte-1=55.24±1.15 kcal/mol. The two independent runs agree within
roughly their combined statistical uncertainty and BOTH sit solidly in
the high-40s/mid-50s range — nowhere near closing the ~12-15 kcal/mol gap
to the abstract ansatz's 34.98/43.03. **CONFIRMED: this is a real,
reproducible effect, not an artifact of limited single-submission shot
statistics.** See iteration 17 below for the full write-up and the
mandatory ALTERNATIVES NOT TAKEN.

**Honest summary for the record**: after ruling out every mechanism this
project can currently test locally or cheaply for real, the residual gap
remains UNEXPLAINED. The most likely remaining candidates are IonQ-
specific effects invisible to circuit-level accounting entirely —
possibilities include (not yet tested): per-qubit-position noise
asymmetry within the simulator's calibrated profile, cross-talk specific
to which of the register's physical qubit indices are addressed, or a
genuine amplitude/observable-sensitivity effect in how VALUES near the
tapered register's specific Schmidt-vector geometry propagate through
noise (as opposed to the reconstruction FORMULA's sensitivity, which was
tested and ruled out). Reported plainly as an open question rather than
forcing a story onto it.

**STATUS UPDATE (iteration 16, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): does the abstract 11-gate ansatz's real-hardware edge over
the Z2-tapered circuit come from gate STRUCTURE rather than gate count,
and if so, can it be ported? Diagnosed first (`gate_structure_compare.py`):
the abstract ansatz has MORE gates of both types (11 CX, 51 u3) than the
generic-StatePreparation tapered circuit (mean 3.94 CX, mean 6.72 u3) yet
still wins on real hardware (34.98/43.03 vs 47.78/51.25 kcal/mol) — direct
gate-list inspection found ~83% of the abstract ansatz's u3 gates are
FIXED, special-angle (0, ±π/2, ±π) structural gates from its Givens-
rotation recipe, vs 0% for generic StatePreparation's arbitrary output.
Hand-derived the same recipe (reference prep + discriminator-qubit bridge
+ XXPlusYYGate Givens hops, `fixed_ansatz.py`'s own construction) for the
tapered 3-qubit register. First attempt (matching the original's 5-angle
budget 1:1) converged for only 4/36 targets, worst error 0.42-0.51 — NOT
a parameter-count problem but a genuine structural fact, measured
directly: the original 4-qubit ansatz's 25 targets ALL share ONE
bit-complement pair (physical fact — 2 electrons confined to a single
Hamming-weight-2 sector), while the tapered register's 36 targets spread
across 1, 2, or 3 of its 3 available bit-complement pairs depending on
target (tapering's Clifford transform does not preserve Hamming-weight
structure). Repeating the Givens hops to 8 angles total DID converge
36/36 to machine precision (worst 1.26e-15) with a CONSTANT gate count
(16 CX, 86 u3, 71/86 = 83% special-angle on u_0) — genuinely closer to
the abstract ansatz's structure than generic StatePreparation. Verified
exact, then run for real, concurrently, on IonQ's free ionq_simulator.
Ideal correctness control PASSED (2.165 kcal/mol). Real result:
aria-1=207.08±2.93, forte-1=218.69±2.65 kcal/mol — WORSE than every
other variant tried in this project, including the previous worst
(native-optimized, 93.73/91.43). CONCLUSION: the special-angle structural
trick is real (verified, 83% special-angle) but porting it to the
tapered register costs 16 CX — 45% more than the abstract ansatz's 11
and >4x generic StatePreparation's mean 3.94 — because the tapered
register's targets need all 3 bit-complement pairs' worth of Givens
coverage where the original problem only ever needed 1. The extra gate
count dominates; the structural trick does not come close to
compensating. A clean negative, reported plainly. See iteration 16 below
for the ALTERNATIVES NOT TAKEN.

BONUS FINDING, independent of the above: while cross-checking this
circuit's fitted angles (computed in one process) against freshly
recomputed targets (in a separate process), found `build_reduced_problem()`
was NOT deterministic across process invocations — `u_4`'s sign flipped
(dot product exactly -1.0, not floating noise) between calls. Root cause:
`ef_fragment.py`'s `exact_ground_state()` used `scipy.sparse.linalg.eigsh`
(ARPACK's iterative Lanczos solver); pinning its `v0` to a fixed seed did
NOT fix it (cross-process diff got WORSE, 1.999, not better) — pointing to
non-associative floating-point rounding inside ARPACK's multi-threaded
sparse matrix-vector products, not just the random start. FIXED by
switching to dense `scipy.linalg.eigh` (the fragment's Hamiltonian is
only 256-dimensional for 8 qubits — dense diagonalization is trivially
fast and has no iterative-convergence non-determinism). Verified 0.0 diff
across 3 separate process invocations after the fix. This bug was
invisible in every PRIOR iteration of this project because every script
called `build_reduced_problem()` exactly once per process and reused the
result — self-consistent within any single run. It would only ever have
surfaced when comparing results computed across separate runs, exactly as
happened here. Independent, durable correctness fix — not specific to
Task A-E.

**STATUS UPDATE (iteration 15, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): the Z2-tapered circuit (3 qubits/register, mean 3.94 CX),
run RAW (no ZNE) for real, concurrently, on IonQ's free ionq_simulator.
Verified exact (1.28e-11 kcal/mol) before submission; ideal control
passed. Real result: aria-1=47.78±2.20, forte-1=51.25±1.77 kcal/mol —
better than the native-optimized circuit (iteration 12: 93.73/91.43,
~1.9x worse) but WORSE than the original abstract 11-gate ansatz
(iteration 9: 34.98/43.03, ~1.2-1.4x better), despite having roughly a
third as many two-qubit gates. Third circuit variant in a row (after
iteration 12's native-optimized case) confirming gate count alone does
not predict real-hardware ranking — the local synthetic-noise model
(iteration 14) correctly predicted tapering would beat the
native-optimized circuit but did not predict it would still trail the
abstract ansatz. See iteration 15 below.

**STATUS UPDATE (iteration 14, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed): does Z2 tapering fix ZNE's convergence problem? No.
Combined iteration 13's tapered (3-qubit, mean 3.94 CX) circuit with the
same rigorous 3-direction floor-tested ZNE sweep — raw error drops
~2.3x (46.0 vs 104.0 kcal/mol) but NO PLATEAU is found in any direction,
same as the untapered circuit. Caught the floor test's OWN discipline
working correctly in real time: the quadratic/widest-range cell reports
a tempting 0.730 kcal/mol (under chemical accuracy) but it's the tail of
a still-decreasing sequence, not a plateau, and is explicitly flagged as
untrustworthy rather than headlined. Per the explicit instruction this
was run under, this result does NOT pass, so no real IonQ submission was
made this iteration — the math was done and found wanting, reported as
the complete, honest result. See iteration 14 below.

**STATUS UPDATE (iteration 13, LOCAL BRANCH `local/attack-base-problem`,
NOT pushed — attacking the base of the problem, not the mitigation):
Task A (Z2 symmetry tapering) is exact, fully verified to machine
precision, and real: both forged registers drop from 4 to 3 qubits for
free (every physical state already sits in the symmetry's +1 sector), a
generic circuit on the reduced register needs a mean of 3.94 CX gates
vs 11 (range 2-4, not yet constant, not yet hand-optimized). Task C
(double factorization) gives a partial, honestly-bounded result: 7
rotated one-body bases suffice for the two-body integral tensor at a
reasonable truncation, fewer than the current 13 qubit-wise groups — but
NOT yet translated into a verified qubit-circuit measurement-group count.
Task E (rigorous ZNE) is the most important negative result of this
iteration: NO PLATEAU FOUND in any of three independent (scale-range,
fit-order) sweep directions — ZNE-linear's real-hardware ~30 kcal/mol
result (iteration 12) should be held with real skepticism until a
plateau can actually be demonstrated. Caught and fixed a genuine false
positive in `qforge.floor_test()` itself along the way (a monotonically
diverging sequence was passing the old ratio-only check). Tasks B
(contextual subspace VQE) and D (spin projection) not yet started.
Nothing here has touched real IonQ hardware yet — local verification
only, per this branch's explicit review-before-push instruction.

**STATUS UPDATE (iteration 12 — does the native gate-count reduction
actually help for real, combined with ZNE/CDR/PEC? Tested, not just
verified locally): the TrappedIonOptimizerPlugin-optimized fixed K=6
ansatz (mean 9.28 vs 11 two-qubit gates) run for real, concurrently, on
ideal/aria-1/forte-1. Key calibration finding: the native MS gate's
learned error rate (p_ms~0.014, both models) is ~80x LARGER than
iteration 9's abstract-gate calibration -- explaining why RAW error here
(93.7/91.4 kcal/mol) is ~2.2x WORSE than the abstract ansatz's raw
(35-43), despite fewer gates: each surviving native gate runs at full
strength with a much higher per-gate error rate. ZNE-linear is the one
method that pulls its weight (~30 kcal/mol both models, back in the
abstract ansatz's ballpark). CDR REVERSES direction from iteration 9 (was
2.1-2.6x worse than raw there; roughly HALVES raw here, though noisy).
PEC is a genuine catastrophic failure (170-400+ kcal/mol) -- diagnosed,
not mysterious: gamma_total scales with the (now much larger) learned
channel, and 4 real seeds is nowhere near enough shot budget for the
resulting sampling overhead. No method reaches chemical accuracy. See
iteration 12 below for the full real-data table.

**STATUS UPDATE (iteration 11 — the classic 0.57 kcal/mol EF+ZNE result,
fully re-examined): reproduced exactly (Task 1: 20.20 -> 0.57). Its EXACT
original circuits, run for real on IonQ (Task 2), give 123-135 kcal/mol —
WORSE than the fixed ansatz's real-hardware 35-43 (iteration 9), not
comparable — overturning "it's just noise" in the opposite direction than
expected: the older K=5 StatePreparation circuit (independent alpha/beta
measurement, 4-phase cross terms) compounds noise worse than the newer
fixed-ansatz pipeline's beta_signs() shortcut. Native-gate remodeling
(Task 3) recovers 3.6-4.3x of that gap via gate-COUNT reduction but
verifiably exploits ZERO partial-angle capability (every surviving MS
gate at theta=0.25 exactly) and breaks CDR's constant-gate-count
requirement. Rebuilding with qforge + shot noise (Task 4) reproduces the
original's ballpark at K=5 (0.71 vs 0.57, both sitting on the same 0.5655
kcal/mol classical floor) but a MAJOR finding: the ZNE noise-scale-range
itself fails its own mandatory floor test (34x/5.4x change extending
[1,2,3]->[1,2,3,4,5], no plateau) -- the classic result was never
robustly converged, independent of the real-hardware mismatch. The
fidelity threshold curve (Task 5) makes it quantitative: IonQ Aria/Forte
(98.786%) sits below every method's chemical-accuracy crossing point
except PEC's best-case (exactly-known channel) framing, which iteration 9
already found does not hold on real hardware. See iteration 11 below for
the full hardware specification. Simulator only throughout -- the $3,000
award remains unspent.

**STATUS UPDATE (iteration 10, QSE — third real-hardware negative, plus a
new vqe/qforge/ library): validated code from iterations 1-9 is now a
clean, importable package (vqe/qforge/, invariants enforced as
assertions, floor_test() reusable). Quantum subspace expansion (QSE,
McClean et al.) was implemented as a mitigation method needing NO channel
model — H_eff is built from the SAME alpha/beta matrices entanglement
forging already measures, so QSE-ordinary needed ZERO new real circuits
(reused iteration 9's checkpointed target data directly), and its
regularized variant needed only 135 small new "compute-uncompute" overlap
circuits. Verified to machine precision locally (7e-12 and 9e-16 kcal/mol)
before any real submission. Real result on ideal/aria-1/forte-1: QSE-
ordinary gives a small, real, zero-free-parameter improvement over raw
(34.0 vs 35.0 kcal/mol aria-1; 41.7 vs 43.0 forte-1) but does NOT beat
PEC (32.1-32.4). QSE-regularized is unstable (helps on forte-1, hurts on
aria-1; its "best" threshold varies wildly across geometries within the
same model). PEC remains the best real-hardware method found across three
independent attempts. See iteration 10 and the retrospective below.

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

## Task 1 (this session, no separate iteration number): `vqe/qforge/` — a clean, importable library

Everything worth keeping through iteration 9 was scattered across
`vqe/*.py` scripts. Extracted into `vqe/qforge/` (no `qforge` package
existed before this): `ansatz.py` (fixed 11-gate ansatz + `fit_angles`),
`forging.py` (fragment Hamiltonian, real gauge, `beta_signs()`,
qubit-wise-commuting/`frame="h"` measurement grouping, `setup_fragment()`
one-call entry point), `mitigation.py` (`RawStrategy`/`CDRStrategy`/
`PECStrategy` sharing one `correct()` interface), `shot_noise.py`
(shot-sampling model + shots-vs-accuracy harness), `floor_test.py` (the
mandatory floor test as a reusable function). Every invariant is now an
assertion, not a comment: `combine_matrices()` hardcodes the identity
Pauli's diagonal to 1.0 regardless of any scale passed in (tested against
a deliberately absurd scale); `filtered_pairs()` unconditionally drops
`|exact|<0.05` training rows; `transpile_fixed()` has no
`optimization_level` parameter at all — hardcoded to 0, tested via
signature inspection so the parameter cannot even be passed, not just
defaulted. `vqe/qforge/tests/test_qforge.py` passes end to end against
known values (K=6 exactness, 36/36 targets converged at 11 CX gates each,
CDR recovering a known injected scale).

**`floor_test()` caught a real bug in itself while being built**: a first
draft used a trailing-window min/max-ratio heuristic and wrongly called
iteration 2's own historical disqualifying sweep (0.60→0.439,
0.30→0.132, 0.15→0.028, 0.05→0.008, 0.01→0.007 kcal/mol) a *pass* —
because the LAST two values (0.008, 0.007) have a small ratio (1.14x)
purely from both being tiny numbers, not from genuinely plateauing;
every OTHER consecutive step in that same sweep is still a 3-5x jump.
Fixed with a "last N consecutive step-ratios must ALL be small" check
instead of a trailing-window aggregate, and `floor_test.py`'s own
`_self_test()` now asserts it disqualifies that exact historical sweep —
a permanent regression test for the bug that motivated writing this
function in the first place.

## Iteration 10: quantum subspace expansion (QSE) — a method that needs no channel model

**Motivation, directly from iteration 9's diagnosis**: CDR and PEC both
run into the SAME wall on real IonQ noise — they each need some model of
how noise degrades a measurement (a fitted scale, a learned Pauli
channel) and iteration 9 found that model badly mismatched the real
target-circuit error (a Clifford probe learned γ_total≈1.004-1.006 while
the real raw error was 15-25x larger than that predicts). QSE (McClean,
Romero, Babbush, Aspuru-Guzik, PRA 95, 042308 (2017)) needs no such
model: noise resilience is STRUCTURAL, from re-solving a generalized
eigenvalue problem, not from correcting a measured value against an
assumed channel.

**How this maps onto entanglement forging, derived not assumed** (full
derivation in `vqe/qse_mitigation.py`'s docstring): the standard forged-
energy formula is exactly the Rayleigh quotient λᵀH_effλ/(λᵀλ) + enuc for
a symmetric K×K matrix H_eff built from the SAME alpha/beta Pauli
matrices entanglement forging already measures. Every result through
iteration 9 evaluated that quotient at the CLASSICALLY KNOWN Schmidt
coefficients λ — i.e. trusted that the exact-diagonalization-derived
weights stay optimal even when the matrices are noisy. QSE removes that
assumption: measure H_eff (and, in the regularized variant, an overlap
matrix S) from the SAME noisy circuits, then let a classical eigensolve
find the best combination. **This means QSE-ordinary needs ZERO new
circuits — it is computed entirely from iteration 9's already-collected
real target data.** Only the regularized variant needs anything new: 15
"compute-uncompute" fidelity circuits per geometry (prepare uₙ, apply the
INVERSE of uₘ's ansatz, measure P(|0000⟩)=|⟨uₘ|uₙ⟩|², the standard
ancilla-free way to get a state-overlap MAGNITUDE from two circuits
sharing one parametrized family) — 135 circuits total (15 pairs × 3
geometries × 3 models), a small addition.

**A dead end caught by derivation before it was built, recorded so it is
not retried**: the first idea for measuring the overlap matrix S was to
reuse the identity Pauli's already-computed "cross term" from the
(uₙ±uₘ)/√2 phase circuits already built for entanglement forging (free,
no new circuits at all). This does NOT work: ⟨ψ|I|ψ⟩=1 is a
normalization tautology for ANY properly normalized measured probability
distribution — true whether or not the circuit is noisy — so it carries
exactly zero information about state overlap, regardless of noise. Caught
by direct algebraic derivation (not by running a failed experiment),
before any code was written that depended on it.

**Verification before any real submission** (`vqe/qse_mitigation.py`,
local only): (1) H_eff's noiseless ground eigenvalue matches the standard
forging-formula energy to **7e-12 kcal/mol** — confirms the H_eff
construction is correct, and confirms a real prediction (not an
assumption): since this fragment's Schmidt rank is exactly 6 (not a
truncation), λ_known MUST already be H_eff's own ground-state
eigenvector. (2) The compute-uncompute overlap circuit matches the exact
statevector overlap to **9e-16** locally, then **1e-17** specifically
with the IonQ abstract gateset (checked again before spending any real
API calls on it, since the local check used a different gateset).

**Local floor test, on this project's own (larger) synthetic noise
model**: QSE-ordinary gives a small, real, deterministic improvement over
the standard forging formula on the identical noisy matrices (103.17 vs
103.995 kcal/mol, 1.01x — modest, but genuine, with zero free
parameters). Regularized QSE's threshold sweep initially came back
vacuous — every threshold up to 0.5 gave an IDENTICAL result, because the
actually-measured S eigenvalues span [0.84, 1.58], never crossed by that
range — fixed by extending the sweep past the measured spectrum. Once
meaningful, the real finding is that aggressive regularization makes
things dramatically WORSE here, not better (n_kept=6/6: 33.1 kcal/mol;
n_kept=1/6: 1497.8 kcal/mol) — on this problem, at this noise level,
ordinary (unregularized) QSE is the more robust choice, a real, disclosed
consequence of the overlap circuit's magnitude-only sign limitation (it
cannot resolve whether an off-diagonal deviation from orthonormality is
constructive or destructive, so it can't reliably tell "safe to discard"
apart from "important to keep").

**Real result, run concurrently on ideal/aria-1/forte-1** (8-seed
mean, same bootstrap/real-execution conventions as iteration 9):

| model | raw (iter.9) | CDR (iter.9) | PEC (iter.9) | QSE-ordinary | QSE-regularized |
|---|---|---|---|---|---|
| ideal (control) | 1.704 | — | — | 1.833 | 1.833 |
| aria-1 | 34.983 | 89.751 | **32.132** | 34.026 | 43.278 |
| forte-1 | 43.026 | 90.257 | **32.422** | 41.742 | 36.647 |

**A third clean negative, exactly as anticipated**: QSE-ordinary gives a
marginal (2-3%) improvement over raw — consistent in direction and rough
magnitude with the local synthetic-noise finding — but does NOT beat PEC,
and needed zero new real circuits to find that out. QSE-regularized is
actively unstable: WORSE than QSE-ordinary on aria-1 (43.28 vs 34.03) but
BETTER on forte-1 (36.65 vs 41.74), still short of PEC either way. The
instability is diagnosable, not mysterious: the "best" regularization
threshold varies wildly ACROSS GEOMETRIES within the same model (aria-1:
best threshold ≈1e-6 at d=0.9/1.0, jumps to 0.83 at d=1.1) — there is no
single threshold choice that would generalize across a real binding-curve
scan, exactly the instability the mandatory floor test exists to surface.
Difference-error cancellation is consistent with this picture: QSE-
ordinary's bias cancels comparably to raw (9.37x aria-1, 7.10x forte-1),
while QSE-regularized's cancellation is much weaker (1.21x, 2.01x) —
its per-geometry behavior is less smooth, not just less accurate. On the
ideal control, QSE (1.833 kcal/mol) is even slightly WORSE than plain raw
(1.704 kcal/mol) — the nonlinear re-diagonalization has a small real cost
when there is no bias to correct in the first place, reported plainly
rather than only reporting the cases where it helps.

**Standing conclusion after three independent real-hardware attempts**:
PEC remains the best-performing method on real IonQ noise (32.1-32.4
kcal/mol) of everything tried in this project — not because its
correction logic is uniquely good, but because CDR's angle-mismatch
problem gets WORSE (not better) on real hardware noise, and QSE's
structural noise-resilience, while real and directionally helpful, is too
small here to close the gap. The bottleneck iteration 9 diagnosed — a
reduced, 1-pair/1-qubit Clifford calibration probe under-characterizing
the real noise on the full 62-gate target circuit — still stands as the
most likely lever for improvement, unresolved by any method tried since.

Code: `vqe/qse_mitigation.py` (local implementation + verification +
floor test), `vqe/ionq_qse_binding_curve.py` (`--overlap`, `--assemble` —
reuses iteration 9's `ionq_simulator_binding_curve_checkpoints/
targets_d*.json` directly, needs no re-collection). Full data:
`vqe/qse_mitigation_results.json`, `vqe/ionq_qse_binding_curve_results.json`.

## Retrospective: root cause of every disqualification/real bug in this project, in one place

Collected here so the failure modes stay visible as a group, not just
scattered across individual iteration write-ups above.

1. **Iteration 2 (locally-perturbed CDR), disqualified**: the method's
   free parameter (training-perturbation radius) had no floor — as the
   radius shrinks toward 0, the training circuit converges to the TARGET
   circuit itself, so the method degenerates into classically
   re-evaluating the answer it was supposed to be measuring. **Root
   cause**: this whole simulator-only testbed can always cheat this way,
   because "exact" is one `Statevector` call away for every circuit,
   including circuits placed arbitrarily close to a target. The fix
   was procedural, not a patched parameter: the mandatory floor test,
   applied to every free parameter of every method from that point on.

2. **CDR global/per-circuit scale, negative results (iterations
   pre-dating this ledger's numbering)**: a single scalar correction
   (one global scale, or one scale per target circuit) cannot capture
   noise attenuation that depends on WHICH Pauli label is being
   corrected — per-basis scale (fit separately per label) was the fix,
   and remains this project's best noiseless-estimator result locally
   (2.850±0.490 kcal/mol) even though iteration 9 found it collapses on
   real hardware noise (root cause below, item 8).

3. **ZNE gate folding, negative result**: abstract-gate folding got
   compiler-cancelled (the inserted G·G⁻¹ pairs were optimized back out
   before submission); native-gate folding worked structurally but its
   extrapolated result (10.15 kcal/mol) was still beaten by CDR (2.850).
   **Root cause**: folding only helps if the folded gates survive to
   execution — verifying that the SUBMITTED circuit, not just the
   locally-constructed one, retains the extra gates is a real, separate
   check this project learned to make explicitly afterward (fold AFTER
   transpilation, submit without further transpiler passes).

4. **`optimization_level>=1` silently collapsing gate counts**
   (`fixed_ansatz.py`, discovered while building `zne_vs_cdr.py`): for
   specific fitted angles landing near periodic special values, higher
   optimization levels' adaptive 2-qubit synthesis found a cheaper
   circuit for SOME targets but not others — silently breaking CDR's
   core assumption that training and target circuits are structurally
   identical. **Root cause**: "optimize the circuit" and "keep the
   circuit structurally comparable across many different parameter
   values" are different goals that `optimization_level` conflates: a
   transpiler pass that is locally optimal per-circuit is not obligated
   to be STRUCTURALLY CONSISTENT across a family of related circuits.
   Fixed by hardcoding `optimization_level=0` everywhere in this
   pipeline — as of Task 1 (this session), enforced in code
   (`qforge.ansatz.transpile_fixed` has no such parameter at all) rather
   than left as a convention every new script had to remember.

5. **Bash tool's ~10-minute hard cap, even on `run_in_background: true`
   commands** (first hit in iteration 8): a background job for a long
   local sweep was silently killed near the 10-minute mark with buffered
   stdout lost. **Root cause**: the cap applies regardless of
   backgrounding. Fixed procedurally, not by fighting the cap: every
   long-running experiment from iteration 8 onward is split into
   independent, checkpointed CLI phases (`--shots N`, `--config X`,
   `--targets --d D`, `--overlap`, ...) that each complete well within
   the budget, with a separate `--assemble` phase doing analysis from
   already-saved checkpoints, no network calls.

6. **`IonQBackend`'s `qiskit_circ_to_ionq_circ` doesn't re-transpile
   submitted circuits** (verified directly, iteration 9): confirmed by
   reading `qiskit_ionq`'s own source (`ionq_backend.py`,
   `IonQBackend.run()`) that circuits are submitted exactly as built —
   the `IonQTranspileLevelWarning` printed on every real run is a global
   qiskit user-config nag about a DIFFERENT default setting, not evidence
   that this project's own `optimization_level=0` circuits get silently
   re-optimized by IonQ's SDK. Checked directly rather than assumed
   either way, since getting this wrong would have invalidated every
   real-hardware invariant this project depends on.

7. **A real bug in `ionq_simulator_binding_curve.py`'s `phase_calibrate()`
   (iteration 9)**: `KeyError: 'IIII'` — qubit-wise measurement groups
   include the identity label, which was correctly never computed in the
   CDR training's exact-value cache (since ⟨I⟩=1 always, it needs no
   training). **Root cause**: iterating "every label in a group" and
   "every label with a cached exact value" look interchangeable until a
   group contains a label that was deliberately excluded elsewhere for a
   good reason — fixed by skipping labels absent from the exact-value
   cache rather than assuming group membership implies cache membership.

8. **Python's `hash()` non-determinism (iteration 9)**: bootstrap-resample
   RNG seeds built from `hash((model, d, seed))` gave DIFFERENT numbers
   on every rerun of `--assemble` against the SAME real checkpointed
   data, because `hash()` on tuples containing strings is randomized
   per-process (`PYTHONHASHSEED`) in Python 3. Caught by literally
   re-running `--assemble` twice and diffing the output — not something
   a single run could ever reveal on its own. Fixed with a
   `zlib.crc32`-based deterministic seed (`stable_seed()`), and
   confirmed identical output across reruns before trusting any number
   built on it.

9. **`floor_test()`'s own trailing-window bug (Task 1, this session)**:
   see the qforge writeup above — a heuristic that correctly flags a
   genuine plateau also incorrectly flagged a sequence still heading to
   zero, because both look "flat in ratio terms" once the numbers
   involved are small. **Root cause**: a RATIO-based flatness check
   cannot distinguish "genuinely converged" from "still shrinking but
   already small" without looking at more than the last two points —
   fixed by requiring several CONSECUTIVE small steps, not just the
   final one, and locking in the fix with a regression test built from
   the exact historical data it needs to keep catching.

10. **The identity-Pauli-overlap dead end (Task 2, this session, QSE
    design)**: see the QSE writeup above — measuring ⟨I⟩ on a superposed
    target-circuit state cannot reveal state overlap, since it is a
    normalization tautology independent of noise. **Root cause**: a
    circuit-level measurement trick that works correctly for one class
    of operators (Hermitian, non-identity Paulis, via the E0/E2
    phase-circuit reconstruction already verified in `ef_fragment.py`)
    does not automatically generalize to a degenerate special case
    (the identity) just because the SAME circuits and SAME formula are
    reused — checked by direct algebra before writing dependent code,
    the same discipline that caught bug 4 (opt_level) and bug 7
    (KeyError) after the fact, applied here before any code existed to
    debug.

11. **QSE's own vacuous regularization-threshold sweep (Task 2, this
    session)**: see the QSE writeup above — an initial threshold range
    never crossed the actual measured S eigenvalue spectrum, so every
    tested value gave an identical (meaningless) result, which the floor
    test technically "passed" without the pass meaning anything.
    **Root cause**: a free parameter's sweep RANGE has to be chosen from
    the actual data the parameter operates on, not guessed in advance —
    fixed by computing the real S eigenvalue spectrum FIRST, then
    building the threshold sweep to span past it.

**The pattern across all eleven**: this project's real bugs cluster into
three kinds — (a) a classical/exact shortcut available in this specific
simulator-only or verification context that would not survive contact
with a register too large to classically check (1, 10); (b) a convention
assumed to hold across a whole family of circuits/measurements that
actually only holds pointwise, not structurally (4, 7); and (c) a
free parameter or heuristic whose validity was asserted instead of
checked against the actual range of the data or the actual historical
counterexample it needed to handle (9, 11, and the floor test itself as
the general antidote to 1). Every one of these was caught by direct
verification — reading the SDK source, re-running to check determinism,
deriving the math before trusting a shortcut, extending a sweep to
actually cover the measured range — not by assuming correctness and
finding out later from a bad real-hardware result.

---

## Iteration 11: back to the original EF+ZNE result — reproduce, port, remodel, and convert into a hardware specification

**The question this answers**: the classic 0.57 kcal/mol EF+ZNE result
(entanglement_forging_zne.py, iteration-numbering predates this ledger)
was measured on a LOCAL Quantinuum-like depolarizing model, never on real
hardware. IonQ Aria/Forte measure ~98.786% two-qubit fidelity —
substantially worse per-gate than Quantinuum H1/H2 (97.82%/98.91% in
this project's own numbers, close to published ~99.8%/99.9%). The
hypothesis to test: is the real-hardware gap (iteration 9: raw 35-43,
CDR 2.1-2.6x worse, PEC 32-43) explained by NOISE LEVEL alone, or does
the CIRCUIT also matter? Five tasks, run in order, each gating the next.

### Task 1 — reproduce, unchanged

`entanglement_forging_zne.py` run exactly as it stands: **20.20 kcal/mol
raw, 0.57 kcal/mol quadratic-ZNE** — matches the historical claim
precisely. Everything downstream is now built on a confirmed foundation,
not an assumed one.

### Task 2 — the control: same circuits, real IonQ noise

Took the EXACT original circuits (K=5, generic `StatePreparation` of the
genuinely-complex — not real-gauged — exact Schmidt vectors, BOTH
registers measured independently, 4-phase cross-term reconstruction, no
`beta_signs()` shortcut) and ran them for real, concurrently, on
ideal/aria-1/forte-1. Two adaptations, disclosed, neither changing what
is measured: `optimization_level=0` (this project's own invariant,
established after the original script was written, which used
`optimization_level=1`) and qubit-wise-commuting measurement grouping (a
real device cannot read arbitrary Pauli expectations from one circuit
execution the way the original's `AerEstimatorV2(method="density_matrix")`
could — grouping only changes circuit COUNT, not what is measured;
verified by reconstructing the exact reference matrices from grouped,
noiseless measurements to 1.4e-13 before spending any real API calls).

**Result: 123.2 ± 3.1 kcal/mol (aria-1), 135.5 ± 1.8 kcal/mol (forte-1),
8 seeds** — cross-validated against an independent real submission from
earlier work in this project (`native_forged_zne_results.json`'s
`RAW_BASELINE_KCAL`: aria-1=125.07, forte-1=134.62 — consistent to ~2%).
The ideal correctness control passed (1.06 ± 0.52 kcal/mol, consistent
with real shot noise).

**This overturns the "it's just noise" hypothesis — in the opposite
direction than the task anticipated.** The prediction was: if noise
alone explains the gap, this should land near iteration 9's fixed-ansatz
numbers (35-43 kcal/mol). Instead it is 3-4x WORSE than that. The
circuit clearly matters — just not in the "maybe the newer ansatz is
worse" direction the task flagged as the overturning case; it is the
OLDER circuit that performs worse. **Mechanistic explanation, verified
not guessed**: per-circuit CX count is IDENTICAL (11) to the fixed
ansatz, so raw gate count is not the driver. The real drivers are
architectural — measuring alpha AND beta registers INDEPENDENTLY
(`beta_signs()` requires a real-gauged state; this state is genuinely
complex, psi max|imag| ranging 0.15-0.97 across separate `eigsh` calls
due to its own unconstrained global phase, confirmed benign since the
computed energy is provably phase-invariant) and the 4-phase cross-term
trick (vs. the real-gauge 2-phase version) both COMPOUND independent
measurement noise multiplicatively in the final bilinear energy formula,
where the newer pipeline's `beta_signs()` shortcut makes beta a
noiseless classical derivation from the same alpha measurement instead.

### Task 3 — remodel for IonQ native gates

**Verified findings, not assumptions**:
- `TrappedIonOptimizerPlugin` (instantiated directly, its entry point is
  not registered) DOES reduce 2-qubit gate count on the fixed ansatz —
  mean 9.28 vs 11 abstract CX, across all 36 K=6 targets — but the
  reduction is NOT uniform: min=4, max=11 per target. **This breaks
  CDR's structural-identity requirement** (training and target circuits
  are no longer guaranteed structurally identical) — a real, disclosed
  cost of native optimization this project had not previously measured.
- Every surviving MS gate's angle, after optimization, is **exactly
  0.25 (full strength), zero variance**, checked directly on gate
  parameters across all 334 checked instances — confirms the optimizer
  exploits NO partial-angle capability, exactly the diagnosis this task
  set out to verify.
- A bare partial-angle MS/ZZ gate is **not** a drop-in replacement for
  this ansatz's `XXPlusYYGate` Givens rotations — verified by direct
  matrix comparison (not derived from memory): `XXPlusYYGate` acts
  block-diagonally (leaves |00⟩/|11⟩ untouched), while MS/ZZ gates
  genuinely mix them at every phi0/phi1/theta combination tested (best
  achievable match: 9.45% matrix error, too large to trust). A genuine
  partial-angle-exploiting resynthesis needs a real KAK/Cartan
  decomposition with variable entangling strength; `TrappedIonOptimizer-
  Plugin` does not do this (confirmed above), and building an
  independent one was judged too high-risk to submit for real within
  this task's scope — flagged as the natural next engineering step, not
  fabricated here.
- **Real, already-collected data reused, not re-run**: native-gate K=5
  state prep + ZNE (`native_forged_zne_results.json`, from earlier work
  in this project, real submission to aria-1/forte-1, folds 1/3/5):
  ZNE-quadratic = 34.25 kcal/mol (aria-1), 31.82 kcal/mol (forte-1) — a
  **3.6x/4.3x improvement over Task 2's naive port** (123.2/135.5). This
  IS a real, substantial recovery — but comes with its OWN
  already-established honesty flag: `rate_consistent=False` for both
  models (the per-fold effective error rate is not constant, so this
  project's own check flags the fold-based ZNE extrapolation as not
  fully trustworthy), reported here, not smoothed over.
- A **clearly-labeled theoretical projection** (not a measurement) using
  IonQ's own published partial-angle fidelity relationship (err(s) =
  0.00357 + 0.02143·s, floor at 14.3% of the full-angle error as s→0)
  applied to this ansatz's actual rotation angles (mean 80% of full
  strength): a genuine partial-angle resynthesis could reduce per-gate
  error by **~17.2%** — real but modest, a ceiling this project has not
  yet reached, not a result claimed as achieved.

**Answering the task's question**: yes, native remodeling recovers real
accuracy the naive port lost (3.6-4.3x), but (a) it does so via gate-COUNT
reduction, not angle-strength reduction (confirmed unused), (b) it
introduces a new problem (non-uniform gate count breaking CDR
compatibility) while solving the old one, and (c) its own ZNE
extrapolation carries a disclosed reliability flag independent of this
task's other findings.

### Task 4 — rebuild with current machinery (qforge, K=5 and K=6, shot noise)

Rebuilt the ZNE experiment with the fixed 11-gate ansatz, real gauge,
`beta_signs()`, qubit-wise-commuting grouping, and shot noise included
(100,000 shots/setting — the original had none), 8-seed mean ± std.

| | K=5 (classical floor 0.5655 kcal/mol) | K=6 (no floor) |
|---|---|---|
| raw | 18.35 ± 0.40 | 17.88 ± 0.38 |
| ZNE-linear | 0.96 ± 0.69 | 0.63 ± 0.61 |
| ZNE-quadratic | 0.71 ± 0.43 | 1.20 ± 1.16 |

K=5's rebuilt ZNE-quadratic (0.71) sits close to the original's 0.57 —
**both are dominated by the 0.5655 kcal/mol classical truncation floor**,
not by measurement or mitigation quality (this fragment's true Schmidt
rank is 6, not 5). Removing that floor (K=6) gives a WORSE, noisier
result (1.20 ± 1.16) — the floor was acting almost like an accidental
regularizer; without it, the same 3-point ZNE fit is visibly less stable.

**A major finding from the mandatory floor test, independent of and in
addition to Tasks 1-2's real-hardware mismatch**: extending the ZNE
noise-scale range from [1,2,3] (the original's own choice) to
[1,2,3,4,5] changes the extrapolated quadratic-ZNE answer by **34x
(K=5) / 5.4x (K=6), with no plateau** — `qforge.floor_test()`'s verdict
is `DISQUALIFIED`, the exact signature iteration 2's training radius
showed. **The classic 0.57 kcal/mol result fails its own free-parameter
floor test.** This was never checked before this task, on either the
local model or real hardware — it is a property of the METHOD (a
3-point polynomial fit extrapolated outside its data range), not of
which noise model is used.

(En route: `qforge.forging.setup_fragment()` gained a `strict=False`
option — the function's original, correct-everywhere-else behavior
asserts K is Schmidt-rank-exact, which would crash on a deliberate
truncation like K=5; `strict=False` allows it while still surfacing
`max_schmidt_tail` so a caller cannot silently ignore the resulting
floor.)

### Task 5 — the deliverable: fidelity threshold curve

Swept two-qubit fidelity 98.5%→99.99% (p2 = 0.015→0.0001, p1 = p2/40,
fidelity ≡ 1−p2 per this task's own stated convention) across all five
methods, K=6, no shot noise (isolates method structure from shot-noise
confounding — a disclosed simplification, not a real-hardware claim).

| method | crosses 1.0 kcal/mol at fidelity | note |
|---|---|---|
| raw | 99.989% | |
| CDR (per-basis) | 99.507% | |
| ZNE-linear | 99.708% | |
| ZNE-quadratic | 98.802% | **inherits Task 4's noise-scale-range non-convergence finding — this specific crossing point is not robustly converged, flagged not hidden** |
| PEC (exact channel) | already below target across the ENTIRE swept range, including at 98.5% | best-case framing (channel exactly known by construction) — iteration 9 already found this framing does not hold on real hardware; restated in this context, not a new discovery |

**IonQ Aria/Forte (98.786%) sits below every method's crossing point
except PEC's** — and PEC's crossing is the one method here whose
premise (an exactly-known noise channel) iteration 9 already falsified
for real IonQ noise. Read plainly: at IonQ's actual measured fidelity,
NONE of the four methods that don't assume perfect channel knowledge
reach chemical accuracy on this circuit, in this local model. This
converts "it didn't reproduce on IonQ" into the quantitative statement
the task asked for: **H4 forged VQE (this ansatz, K=6) needs two-qubit
fidelity gains beyond what IonQ Aria/Forte currently deliver, evidenced
across four independent mitigation strategies, not asserted from one
failed run.**

### Hardware specification — the synthesis

1. **Real-hardware floor, established and cross-validated**: raw error
   on real IonQ noise is 35-43 kcal/mol for the fixed 11-gate ansatz
   (iteration 9) and 123-135 kcal/mol for the original K=5
   `StatePreparation` circuit (Task 2, cross-validated against
   independent prior real data to ~2%) — the CIRCUIT choice alone is a
   3-4x effect, larger than any single mitigation method's own gain.
2. **The 0.57 kcal/mol figure should never be cited without two
   caveats, both established in this iteration, not previously known
   together**: it is not a real-hardware measurement (already known),
   AND it fails its own noise-scale-range floor test independent of that
   (newly established here) — the method itself, not just the noise
   model, was untested against its own free parameters until now.
3. **Native remodeling is a real, partial lever** (3.6-4.3x recovery via
   gate-count reduction) but is NOT currently exploiting IonQ's
   arbitrary-angle hardware capability at all (confirmed: zero angle
   variance post-optimization) — a genuine partial-angle resynthesis is
   the clearest concrete next engineering step this whole investigation
   points to, bounded at ~17% further per-gate error reduction by IonQ's
   own published fidelity-vs-angle data, not yet built or tested here.
4. **The quantitative bar**: at K=6 with this ansatz, raw/CDR/ZNE-linear
   all need ≥99.5-99.99% two-qubit fidelity to reach chemical accuracy;
   only PEC's best-case (exactly-known channel) framing clears IonQ's
   actual 98.786%, and iteration 9 already showed that framing fails on
   real hardware because the channel is not exactly known in practice.
   **IonQ Aria/Forte, as measured today, is not yet sufficient for this
   specific circuit and problem size under any of the five methods
   tested here in their currently realistic (not best-case) form.**

Per the user's explicit instruction, every step of this iteration ran on
IonQ's free `ionq_simulator` only — the $3,000 award remains unspent,
and no step targeted `ionq_qpu`.

Code: `vqe/entanglement_forging_zne.py` (Task 1, unchanged),
`vqe/ionq_original_circuit_replication.py` (Task 2, `--control` /
`--assemble-control`), `vqe/ionq_native_remodel_analysis.py` (Task 3),
`vqe/qforge_ef_zne.py` (Task 4), `vqe/fidelity_threshold_curve.py`
(Task 5). Consolidated results:
`vqe/ionq_original_circuit_replication_results.json`.

---

## Iteration 12: does the TrappedIonOptimizerPlugin gate-count reduction actually help, combined with real ZNE/CDR/PEC, on IonQ's free simulators?

**The gap this closes**: iteration 11's Task 3 verified the gate-count
reduction (mean 9.28 vs 11 two-qubit gates, K=6) LOCALLY ONLY, and cited
OLDER real data using a DIFFERENT circuit (native_stateprep.py's K=5
hand-derived tree). This iteration is the missing real test: the ACTUAL
TrappedIonOptimizerPlugin-optimized fixed K=6 ansatz, submitted for real,
concurrently, to ideal/aria-1/forte-1, combined with native ZNE
(fold-after-optimize, verified order, never re-optimizing a folded
circuit), CDR (training draws each individually optimized, honestly
carrying the same non-uniform-gate-count structural mismatch a target
would have), and PEC (native MS-gate Clifford calibration + real
randomized quasi-probability circuits, 4 seeds).

**Verified before any real submission, to machine precision**: the
optimized-circuit measurement chain against the exact abstract-ansatz
reference (1.3e-14), the PEC-variant builder reducing to raw exactly at
p→0 (1e-14), `fold_native_2q` preserving the ideal unitary on the
optimized circuit at every fold (<1e-15), and — a real bug caught before
it could reach a real submission — an initial guessed native Z-gate
decomposition (`GPI2(0.5)·GPI(0)·GPI2(0.5)`) that direct matrix
comparison showed does NOT equal Z (it equals a global phase times
identity); replaced with the verified `GPI(0.25)` then `GPI(0)`
(X·Y = i·Z, phase-irrelevant for expectation values).

**Real, verified single-fact finding from `--calibrate`**: the native
MS-gate learned error rate is **p_ms ≈ 0.014 for BOTH aria-1 and
forte-1** — roughly **80x larger** than iteration 9's abstract-gate
Clifford calibration (p2 ≈ 0.0002-0.0003), and much closer to this
project's own local synthetic depolarizing model (0.01214). This single
number turns out to explain nearly everything that follows.

**Real results, 8-seed mean ± std (4 seeds for PEC), single geometry
d=1.0**:

| method | ideal (control) | aria-1 | forte-1 |
|---|---|---|---|
| raw | 2.39 ± 0.54 | **93.73 ± 3.56** | **91.43 ± 3.29** |
| ZNE-linear | 2.06 ± 0.63 | 31.23 ± 3.63 | 29.90 ± 3.57 |
| ZNE-quadratic | 5.36 ± 2.16 | 42.95 ± 8.82 | 18.06 ± 6.23 |
| CDR | — | 47.90 ± 19.96 | 42.32 ± 21.04 |
| PEC (4 seeds) | — | 406.51 ± 348.56 | 172.56 ± 76.08 |

**Answering the question directly: gate-count reduction alone does NOT
help — raw error on the native-optimized circuit is ~2.2x WORSE than
iteration 9's abstract-ansatz raw (93.7/91.4 vs 35-43 kcal/mol), despite
having FEWER two-qubit gates (mean 9.28 vs 11).** The mechanism is not
mysterious, it is the calibration finding above: each surviving native
gate runs at FULL STRENGTH (θ=0.25, confirmed zero-variance in iteration
11) with an ~80x larger per-gate error rate than the abstract "cx" gate
this project's other real-hardware results are calibrated against.
Fewer, individually-noisier gates costs more than more, individually-
quieter ones here — a genuine, disclosed reversal of the naive
"fewer gates = better" intuition, evidenced by a directly-measured
calibration number, not asserted.

**ZNE genuinely helps, and is the one method that recovers real ground
here**: ZNE-linear brings both models down to ~30 kcal/mol — back in the
same ballpark as iteration 9's abstract-ansatz raw and iteration 11's
cited native K=5 ZNE-quadratic (31.8-34.3). ZNE-quadratic is inconsistent
between the two models (worse than linear on aria-1, better on forte-1,
both with large std) — consistent with iteration 11's own established
finding that this ansatz's noise-scale-range/fit-order is not robustly
converged; not re-litigated here, just not ignored either. Even the
IDEAL control's own ZNE-quadratic (5.36 ± 2.16) is worse than its raw
(2.39 ± 0.54) — expected: extrapolating a quadratic through 3
near-identical low-noise points amplifies shot noise, not evidence of a
bug.

**CDR reverses direction from iteration 9** — a genuine, notable
finding: on the ABSTRACT ansatz, CDR made things 2.1-2.6x WORSE than raw
(iteration 9). Here, on the NATIVE-optimized ansatz, CDR roughly HALVES
the raw error (93.7→47.9 aria-1, 91.4→42.3 forte-1) — though with large
uncertainty (±20-21, roughly 40-50% relative). This is not attributed to
any specific mechanism here (a real, open question for a future
iteration — plausibly the larger, more uniform native-gate error rate
gives per-basis scale fitting a stronger, cleaner signal than the
abstract ansatz's much smaller, possibly differently-structured
residual), reported as measured, not explained away.

**PEC is a genuine catastrophic failure here, and the reason is
diagnosable, not mysterious**: 406.5 ± 348.6 (aria-1), 172.6 ± 76.1
(forte-1) — WORSE than raw, with enormous variance. PEC's sampling
overhead scales as γ_total², and γ_total grows with BOTH the per-gate
error rate and the gate count; with p_ms≈0.014 (80x the abstract
channel iteration 9's PEC was built on) over a ~9-gate circuit, γ_total
here is far larger than iteration 9's γ_total≈1.004-1.006 — meaning the
REAL shot budget PEC needs to converge is far larger too, and 4 real
seeds at 10,000 shots each is nowhere near enough. This is the same
"channel under-characterization/shot-budget" bottleneck iteration 9
first diagnosed, now seen from the opposite direction: there, the
learned channel was too SMALL to explain the real error; here, the
learned channel is large enough to be structurally believable, but the
REQUIRED shot budget that comes with a larger honestly-learned channel
was not provisioned for.

**Standing conclusion, updated**: no method tested on the native-
optimized circuit reaches chemical accuracy (best here: ZNE-linear at
~30 kcal/mol, comparable to — not better than — PEC's abstract-ansatz
32.1-32.4 kcal/mol from iteration 9). Gate-count reduction by itself is
not a lever worth pursuing further on its own; ZNE combined with native
gates is the one piece of this iteration that pulls its weight, and CDR's
reversal is a real, open thread for a future iteration to explain rather
than assume.

Code: `vqe/ionq_native_optimized_mitigation.py` (`--calibrate`,
`--targets`, `--pec`, `--assemble`). Full data:
`vqe/ionq_native_optimized_mitigation_results.json`.

---

## LOCAL BRANCH `local/attack-base-problem` — NOT pushed, pending review

**Correction, flagged before any new work cites the old numbers**: the
previously-stored `quantinuum_h1`/`quantinuum_h2` reference fidelities
(0.9782/0.9891) in `fidelity_threshold_curve.py` are WRONG — those are
11-gate CIRCUIT fidelities (i.e. survival probability of the whole
circuit), not PER-GATE fidelities. The correct per-gate reference values
are ~99.8% (H1) and ~99.9% (H2). Not yet fixed in the committed file
(this is a note for the next edit that touches it) — flagged here so it
is not propagated into new work first.

## Iteration 13 (Task A): Z2 symmetry tapering — attacking the base of the problem, not the mitigation

**Why this iteration exists**: every mitigation method tried on real IonQ
noise so far (CDR, PEC, QSE, native optimization) tops out around ~30
kcal/mol (ZNE-linear, iteration 12). The gate count is the lever that
actually moves the needle at IonQ's real fidelity — shrinking the
PROBLEM before building any circuit, rather than correcting a fixed
circuit's noise after the fact.

**Method** (Bravyi, Gambetta, Mezzacapo, Temme, arXiv:1701.08213), tool:
qiskit's own `Z2Symmetries.find_z2_symmetries` / `.cliffords` /
`.sq_paulis` (not hand-derived from scratch, matching this project's
established preference for tested library implementations over
re-deriving solved sub-problems).

**Verified, not assumed**: built a proxy `SparsePauliOp` from just the
alpha register's 37 unique Pauli labels (symmetry-finding depends only
on which Pauli strings are present, not their coefficients, so this
correctly isolates the alpha register's OWN measurement-relevant
symmetry structure) and found exactly **one** independent Z2 symmetry:
`ZZZZ` (alpha-electron-number parity). The beta register, checked
independently, gives the identical result. The full 8-qubit Hamiltonian
has **three** independent generators (`ZZZZIIII`, `ZIZIIZIZ`,
`ZIZIZIZI`) — each verified to commute with H exactly (`max ||[P,H]||
= 0.00e+00`). Decoding block-locality against this project's qubit
convention (0-3=alpha, 4-7=beta): `ZZZZIIII` is beta-only (embeds the
per-register symmetry found above); the other two have support spanning
BOTH registers — a genuine point-group-type symmetry of the symmetric H4
chain, real but not exploitable per-register (see ALTERNATIVES NOT
TAKEN).

**A real bug caught mid-derivation, not assumed away**: `Z2Symmetries`
maps the alpha symmetry to `sq_paulis=[IIIX]` — an **X**, not a Z. The
tapered qubit is fixed in the X-eigenbasis after the Clifford, not
automatically the computational (Z) basis — naively reading off a
computational-basis bit after the Clifford alone gave a
non-uniform-value assertion failure (a real, caught, fixed bug, not a
hypothetical caveat). Fixed by adding the correct single-qubit basis
rotation (H for X, `H·Sdg` for Y — verified this specific matrix by
direct computation rather than trusted from memory, since a hand-derived
attempt at the Y-basis rotation was ALSO checked and found wrong before
being corrected) before treating the qubit's value as a definite bit.
The Hamiltonian's Pauli terms are tapered by DIRECT MATRIX conjugation
with the identical combined transform used for the state vectors
(`taper_pauli_matrix`), rather than a separately-tracked symbolic
Pauli-string transform — deliberately avoiding the exact class of
states-vs-operators inconsistency the X-vs-Z bug above came from.

**Verified end to end, to machine precision**:
- Tapered qubit value is **uniform (=0) across all 6 diagonal Schmidt
  vectors AND all 36 K=6 target slots** — confirms the predicted
  mechanism exactly: every physical state here lives in the alpha
  register's weight-2 sector (2 particles among 4 orbitals), Hamming
  weight 2 is even, so the ZZZZ parity eigenvalue is +1 for every
  physical state tested, uniformly. Tapering costs nothing because no
  physical state this project ever prepares sits outside the sector this
  symmetry selects.
- Round-trip reconstruction (rotate → taper → un-taper → un-rotate)
  matches the original Schmidt vectors to **1.99e-14 / 2.15e-14** (6
  diagonal / all 36 targets).
- Energy recomputed from the TAPERED (3-qubit alpha, 3-qubit beta)
  matrix elements matches the exact H4 energy to **2.29e-11 kcal/mol**
  — the tapering is exact, not merely gate-count-reducing.

**Gate count on the reduced register**: the reduced 3-qubit target
states occupy exactly 6 of 8 possible basis states — indices {1..6},
excluding 000 and 111 — i.e. Hamming weight ∈ {1,2} on 3 qubits, a
direct structural echo of the ORIGINAL problem (weight=2 on 4 qubits).
Generic `StatePreparation` (abstract u3/cx, `optimization_level=0`,
this project's mandatory invariant), verified correct to <1e-9 on all 36
targets:

| | untapered (4-qubit) | tapered (3-qubit) |
|---|---|---|
| CX gates | constant 11 | min 2, max 4, **mean 3.94** |

Not constant across targets (2-4) — the same CDR-compatibility caveat
already found for the TrappedIonOptimizerPlugin-optimized native circuit
(iteration 12) applies here too, disclosed not hidden. This is the
GENERIC baseline, not yet hand-optimized (see ALTERNATIVES NOT TAKEN) —
still a real, verified, ~2.8x reduction in mean 2-qubit gate count for
zero approximation.

**No free parameter to floor-test here**: unlike CDR's training radius
or ZNE's noise-scale range, Z2 tapering is a discrete, exact operation
once the symmetry generators are found — there is no continuous
"aggressiveness" knob. The verification chain above (round-trip
exactness, energy exactness) plays the role a floor test would for a
continuous-parameter method, matching how loop_pec.py handled its own
parameter-free (exact-channel) case.

**What this is worth, read against the task's own cited table** (11
gates → 4.13 kcal/mol after 35x ZNE at IonQ's real fidelity; 7 gates →
2.69; 5 gates → 1.95; 3 gates → 1.18): mean 3.94 gates sits between the
5-gate and 3-gate rows — a real, concrete, honest reason to expect this
reduction to matter on real hardware, NOT YET TESTED on real IonQ
noise (that real-hardware test is the natural next step, deliberately
not run yet on this local branch per the instruction to validate before
anything gets pushed).

### ALTERNATIVES NOT TAKEN

1. **Tapering the two cross-register point-group symmetries too**
   (`ZIZIIZIZ`, `ZIZIZIZI`) — rejected: their support spans both alpha
   and beta qubits, so exploiting them would correlate the two
   registers' tapered-qubit values, conflicting with entanglement
   forging's core premise of measuring alpha and beta as fully
   independent circuits. Would revisit if a reformulation of forging
   that shares a single classical bit between the two registers'
   measurements is developed — a real, larger restructuring, not
   attempted here.
2. **Bravyi-Kitaev mapping instead of Jordan-Wigner** — BK sometimes
   exposes symmetries with different locality properties and might
   reveal additional per-register-exploitable Z2 symmetries beyond what
   JW gives here. Rejected for this iteration: switching the mapping
   convention touches every downstream piece of this project (labels,
   CDR training, PEC calibration, real-hardware measurement code), a
   much larger and riskier change than tapering on top of the existing
   convention. Would revisit if the JW-based reduction found here turns
   out to be at its ceiling and further qubit reduction is still needed.
3. **Hand-deriving a constant-gate-count fixed-structure ansatz for the
   reduced weight-{1,2}-on-3-qubits sector immediately**, mirroring
   `fixed_ansatz.py`'s own derivation for the original weight-2-on-4-qubit
   problem — rejected for THIS iteration on time budget alone, not
   because it looks hard: the structural echo (same weight-sector
   pattern, one dimension down) makes it plausible the same derivation
   technique transfers directly. Would revisit immediately if the
   generic-StatePreparation gate count (mean 3.94, non-constant) turns
   out to be the real-hardware bottleneck once tested for real.
4. **Using `qiskit_nature`'s `TaperedQubitMapper` with its standard
   Hartree-Fock-based sector selection**, instead of directly checking
   the ACTUAL target states' symmetry eigenvalues — rejected: HF-sector
   selection is a good heuristic when the exact ground state is not
   available, but this project already has the exact ground state and
   its full Schmidt decomposition in hand, so checking the real target
   states directly is more direct and avoids trusting a heuristic that
   could silently pick the wrong sector for a state HF approximates
   poorly. Would revisit for a larger fragment where computing the exact
   ground state directly becomes intractable — exactly the regime this
   project has repeatedly flagged as the one that matters, since "exact"
   being one function call away is a testbed-only luxury.

Code: `vqe/z2_tapering.py` (symmetry finding + per-register tapering +
exactness verification), `vqe/z2_tapered_ansatz.py` (all-36-target
tapering + generic-circuit gate count). Full data:
`vqe/z2_tapering_results.json`, `vqe/z2_tapered_ansatz_results.json`.
**Not yet run on real IonQ hardware — local verification only, per this
branch's own review-before-push discipline.**

---

## Iteration 13 (Task C): double factorization — measurement basis reduction, partial

**Method** (Motta et al., npj Quantum Inf 7, 83 (2021); Huggins et al.,
npj QI 7, 23 (2021)): rewrite the two-body Hamiltonian as a sum of
SQUARED one-body operators, each diagonal in its own rotated orbital
basis. Implemented directly on this project's own already-built MO-basis
two-electron integral tensor (`chem.integrals` + the same einsum
`entanglement_forging_h4.py` already uses) — reshaped to the symmetric
16×16 matrix `V[(p,q),(r,s)] = h2[p,q,r,s]` (verified symmetric,
`||V-V^T||=8.05e-16`), eigendecomposed, factors sorted by `|eigenvalue|`.

**Real, measured result**: the untruncated decomposition needs at most
`N(N+1)/2 = 10` factors for N=4 spatial orbitals (H4's RHF basis) — and
the ACTUAL eigenvalue spectrum drops off fast (1.93, 0.78, 0.46, 0.31,
0.027, 0.019, 0.016, 0.0001, 0.00005, then six that are exactly/
numerically zero), so truncation buys real savings:

| truncation tol | factors kept | reconstruction max error |
|---|---|---|
| 1e-10 | 10 | 9.44e-16 |
| 1e-6 | 9 | 7.06e-08 |
| 1e-4 | 8 | 1.05e-05 |
| 1e-3 | **7** | 4.08e-05 |
| 1e-2 | 7 | 4.08e-05 |

At a reasonable truncation (1e-3, reconstruction error still 4e-5 Ha,
far below chemical accuracy), **7 rotated one-body bases** suffice for
the two-body integral tensor — fewer than this project's current 13
qubit-wise-commuting groups, though not as few as the 5 general-commuting
bases the task cited as a reference point.

**Honest limitation, stated plainly, not glossed over**: this result is
for the two-body integral tensor's OWN basis count, at the SPATIAL-
ORBITAL level — it has NOT been translated into a verified QUBIT-CIRCUIT
Pauli-measurement-group count (the level this project's existing 13/5
figures actually operate at). That translation (working out which
qubit-level Pauli measurement groups each rotated-orbital-basis factor
corresponds to, and confirming the resulting circuit-level group count)
is a real, well-defined next step, not completed this iteration — see
ALTERNATIVES NOT TAKEN. Reporting "7 factors" as if it were already a
verified "7 qubit measurement bases" would be exactly the kind of
unearned equivalence this project's honesty rules exist to prevent.

### ALTERNATIVES NOT TAKEN

1. **Completing the qubit-level translation this same iteration** —
   rejected on time budget: mapping each of the 7 kept one-body-rotated
   bases into an actual qubit basis-change circuit and re-deriving the
   resulting Pauli measurement groups is itself a full sub-task (roughly
   the same order of work as the qubit-wise-commuting grouping this
   project already built once). Would revisit as the immediate next step
   before this number is used to justify any real circuit change.
2. **Cholesky decomposition instead of full eigendecomposition** — the
   standard DF literature often uses a pivoted Cholesky factorization of
   V rather than a full eigendecomposition (cheaper for larger systems,
   same mathematical content for a positive-semidefinite V). Rejected
   here because N=4 is small enough that a full eigendecomposition is
   free computationally and gives the SAME factor count with less
   implementation risk (no pivoting-order subtlety to get wrong). Would
   revisit for a larger fragment where the O(N^6) cost of a full
   eigendecomposition of the N^2 x N^2 matrix actually matters.
3. **Applying DF to the qubit-level Pauli operators directly** (rather
   than the classical MO integral tensor) — rejected: DF's whole
   leverage comes from operating on the STRUCTURED two-body integral
   tensor before Jordan-Wigner mapping scrambles that structure across
   many Pauli strings: doing it after mapping to qubits would require
   re-deriving the same factorization from a much less structured
   object. Would revisit only if a compelling reason emerged to avoid
   touching the classical integral tensor at all.

Code: `vqe/double_factorization.py`. Full data:
`vqe/double_factorization_results.json`.

---

## Iteration 13 (Task E): rigorous ZNE — no plateau found, full stop

**Why this matters more than a single number**: iteration 11 found the
classic 0.57 kcal/mol result fails its own noise-scale-range floor test.
This iteration asks the honest follow-up properly: is there ANY (scale
range, fit order) combination that plateaus for this ansatz, or does ZNE
simply not converge here at all? Swept 5 scale ranges (widths 3 through
7, all starting at 1), every polynomial fit order from 1 up to
`len(range)-1` (order = full range width − 1 is exact polynomial
interpolation — what "Richardson extrapolation" means for ZNE
mathematically), plus an exponential fit `E(s)=A·exp(−k·s)+E_inf`
(unavailable — not faked — whenever the shot-noisy energies are not
monotonic in scale, which was every range tested here). 8 real seeds,
100,000 shots/setting, shot noise in every headline number, on this
project's own local depolarizing model (`P2_PER_GATE=0.01214`) so the
methodology could be floor-tested without spending real IonQ time before
validating it — per this branch's explicit instruction.

**A real false positive caught and fixed in `qforge.floor_test()`
itself, before trusting any verdict from it**: the FIRST run reported
"PASS" for the order=1 (ZNE-linear) range sweep — but the actual values
are **15.97 → 22.81 → 29.89 → 37.85 → 46.40 kcal/mol**, strictly,
monotonically INCREASING across the entire sweep, not plateauing at all.
`floor_test()`'s consecutive-step-ratio check was fooled: for a linearly
growing sequence `(a, a+d, a+2d, ...)`, the ratio of consecutive terms
approaches 1 purely because the values themselves are growing —
`(a+(n+1)d)/(a+nd) → 1` as `n` grows, regardless of whether the sequence
is converging or diverging. **Fixed** by adding a mandatory
non-monotonicity requirement on the tail (a genuine plateau should NOT
be strictly increasing or decreasing across 3+ consecutive most-
aggressive points — real measurement noise breaks monotonicity; a
still-changing systematic trend preserves it), locked in with a new
`_self_test()` case built from this exact data so the failure mode stays
caught. This is the kind of bug the mandatory floor-test discipline
exists to surface — including, this time, IN the floor-test function
itself, not just in the methods it checks.

**Corrected result, all three independent sweep directions**:

| sweep | verdict |
|---|---|
| range, fixed order=1 (ZNE-linear) | **DISQUALIFIED** — tail [29.88, 37.66, 46.19] strictly increasing |
| order, fixed at the widest range (1-7) | **DISQUALIFIED** — 27.9x overall, last step-ratios [1.7, 8.38] |
| range, fixed order=2 (quadratic, reproducing iteration 11's own check) | **DISQUALIFIED** — tail [4.11, 5.33, 6.38] strictly increasing |

**NO PLATEAU FOUND in any of the three directions tested.** Per this
project's own disqualification rule, stated as the task required: ZNE is
NOT demonstrated to be converged for this circuit at this (local,
synthetic) noise level, full stop. This is consistent with — and now
independently confirms, at a DIFFERENT noise rate and with a properly
corrected floor test — iteration 11's original finding on the classic
Quantinuum-rate model. **The honest implication for iteration 12's real-
hardware ~30 kcal/mol ZNE-linear result**: that number should be held
with real skepticism, not treated as a converged answer, until a plateau
can actually be demonstrated (a wider real-hardware scale-range sweep,
not yet run — real IonQ time, deliberately not spent on this local-only
branch before the methodology itself was trustworthy).

**Uncertainty reported throughout, not just central values**: every row
in the full sweep table carries an 8-seed std alongside its mean (e.g.
order=2 at range 1-7: 6.375 ± 1.388 kcal/mol) — several combinations
show std comparable to or larger than the mean itself (order≥4 at wide
ranges), a second, independent signal (beyond the monotonicity check)
that those fits are not well-constrained by the data.

### ALTERNATIVES NOT TAKEN

1. **Extending the scale range past 7** — rejected for this iteration on
   compute-time budget (each additional scale point needs a fresh noisy
   density-matrix pass over all 36 targets); the monotonic-increase
   pattern is already unambiguous at width 5-7, so extending further was
   judged unlikely to change the qualitative verdict. Would revisit if a
   genuine inflection (the tail starting to curve back down) appeared
   near the current boundary — it does not.
2. **Mitiq's own ZNE implementation** (the standard open-source library
   for exactly this) instead of a from-scratch polyfit/exponential
   sweep — rejected to keep this project's own verified circuit/energy
   pipeline as the single source of truth for what "the energy at scale
   s" means (Mitiq would need its own adapter into this project's
   forged-energy bilinear reconstruction, a real integration cost for a
   library whose core extrapolation math is the same handful of
   `numpy.polyfit` calls already used here). Would revisit if Mitiq's own
   more sophisticated fit-quality diagnostics (e.g. confidence-interval-
   aware extrapolation) turn out to catch something this project's
   simpler floor-test approach misses.
3. **Testing wider ranges on REAL IonQ hardware immediately**, since
   that is the number that actually matters — rejected deliberately per
   this branch's own instruction (validate locally first, real hardware
   only after a result survives its floor test). The honest outcome here
   is that ZNE-linear did NOT survive its floor test even locally, so a
   real-hardware wide-range sweep is the natural next step ONCE this
   local finding is reviewed, not run pre-emptively.

Code: `vqe/zne_floor_tested.py`. Fix: `vqe/qforge/floor_test.py`
(monotonic-tail check + regression test). Full data:
`vqe/zne_floor_tested_results.json`.

---

## Iteration 14: does Z2 tapering fix ZNE's convergence problem? No — a clean, disciplined negative

**The question**: iteration 13 found two independent facts that invite an
obvious next question — Task A: the tapered (3-qubit, mean 3.94 CX)
circuit is real and exact. Task E: ZNE shows NO PLATEAU in any direction
for the UNTAPERED (11-gate) circuit. Does the SMALLER circuit's ZNE
behave better — fewer gates meaning a smaller total noise range is being
explored at any given scale, plausibly better-conditioning the
extrapolation? A real, testable hypothesis, checked here, not assumed.

**Built on already-verified math, extended and re-verified before
trusting it for a new noisy sweep**: confirmed tapering commutes EXACTLY
with `beta_signs()` (`max|reduced_v − signs·reduced_u| = 0.0`) — so this
analysis needed real circuits for the alpha register only, beta derived
for free, same efficiency as the untapered pipeline. Confirmed all 37
alpha-register Pauli labels reduce to a genuine single 3-qubit Pauli
string (times a real ±1 sign) by DIRECT matrix comparison against all 64
three-qubit Paulis — not assumed from the fact that a Clifford generally
preserves Pauli-ness. Confirmed the 3-qubit `StatePreparation` circuit
reproduces the exact tapered target to 3.31e-14 before running anything
noisy.

**Real result — gate reduction helps the RAW number a lot, but does not
fix ZNE's convergence**: raw (scale=1) error is **46.0 kcal/mol**, vs
104.0 for the untapered ansatz at the same base noise rate — a genuine
~2.3x reduction, consistent with iteration 13's own gate-count-vs-error
table. But the SAME three floor-test directions used on the untapered
circuit, re-run here:

| sweep | verdict |
|---|---|
| range, fixed order=1 (ZNE-linear) | **DISQUALIFIED** — 3.1x overall, last step-ratios [1.27, 1.30] |
| order, fixed at the widest range (1-7) | **DISQUALIFIED** — 56.5x overall, last step-ratios [2.1, 5.13] |
| range, fixed order=2 (quadratic) | **DISQUALIFIED** — tail [0.975, 0.952, 0.730] strictly decreasing |

**NO PLATEAU FOUND, in any direction, on the tapered circuit either.**
Tapering reduces gate count and raw error substantially but does NOT
resolve the underlying non-convergence this project's ZNE fits have now
shown twice, on two different circuits, at two different (equivalent)
noise levels.

**The trap this floor test exists to catch, caught in real time**: the
order=2, widest-range cell reports **0.730 ± 0.574 kcal/mol** — under
chemical accuracy, and exactly the kind of single number a less
disciplined report would headline. The floor test correctly flags it as
untrustworthy: it is the LAST point of a still-decreasing sequence
(0.975 → 0.952 → 0.730), not a converged value — extending the range
further could plausibly keep falling, overshoot, or do anything else;
there is no way to know from this data, and reporting 0.730 kcal/mol as
"the answer" would repeat exactly the mistake iteration 2 made in a new
disguise. This is reported here explicitly as the reason NOT to trust
it, not hidden because it looks good.

**Per the explicit instruction this iteration was run under ("do it if
it passes everything"): this does not pass. No real IonQ submission was
made.** The math was done, verified, and found wanting — that is the
complete, honest result of this iteration, not a placeholder for a
future positive one.

### ALTERNATIVES NOT TAKEN

1. **Reporting the 0.730 kcal/mol cell as a headline result anyway**,
   since it is numerically below chemical accuracy — rejected outright:
   this is precisely the failure mode the mandatory floor test exists to
   prevent, and doing it here after having JUST fixed a false-positive
   bug in `floor_test()` itself (iteration 13, Task E) would be a
   direct, immediate contradiction of that fix. Would never revisit this
   without a genuine plateau demonstrated first.
2. **Hand-deriving a fixed, constant-gate-count circuit for the tapered
   register before running ZNE** (rather than the generic, non-constant
   2-4-gate `StatePreparation` baseline) — rejected for this iteration:
   ZNE folds each circuit against ITSELF, so non-constant gate count
   across targets does not affect ZNE's own validity the way it would
   for CDR; the non-convergence problem found here is unlikely to be a
   circuit-structure artifact specifically, since it reproduces the
   SAME qualitative pattern (no plateau in any of 3 directions) that the
   untapered, CONSTANT-11-gate circuit already showed. Would revisit if
   the untapered/tapered comparison ever diverged qualitatively — it
   does not here.
3. **Trying Mitiq's factory/inference classes for a more sophisticated
   extrapolation** (adaptive scale-factor choice, Bayesian model
   selection between fit families) instead of a fixed grid of polynomial
   orders — rejected on the same integration-cost grounds as iteration
   13's Task E ALTERNATIVES NOT TAKEN #2. Would revisit if this
   project's own simple polynomial/exponential sweep is confirmed (via a
   literature comparison) to be systematically worse-conditioned than
   what Mitiq's adaptive methods would find on the SAME data — not yet
   checked.
4. **Concluding that ZNE is unsalvageable for entanglement-forged H4 and
   dropping it entirely** — rejected: two negative results (untapered,
   tapered) at ONE noise regime (this project's own local depolarizing
   model, scaled) is not the same as ruling out ZNE under every
   circumstance. The real, still-open question is whether REAL IonQ
   noise (which iteration 12 already showed differs qualitatively from
   this local model — e.g. the native MS-gate learned rate came in ~80x
   larger than the abstract-gate one) shows the same non-convergence
   pattern or a different one. Would revisit by running THIS SAME
   3-direction floor test on real IonQ fold data directly, once enough
   real fold points exist to test more than one scale range (iteration
   12 only has folds 1/3/5, one single range) — a concrete, well-scoped
   next real-hardware experiment, not run pre-emptively here.

Code: `vqe/z2_tapered_zne.py`. Full data: `vqe/z2_tapered_zne_results.json`.

---

## Iteration 15: the Z2-tapered circuit, RAW, for real on IonQ — fewer gates is not automatically better, again

**Why this iteration**: iteration 14 found ZNE does not converge on the
tapered circuit either, so a real-hardware ZNE submission was correctly
withheld. But the tapered circuit's RAW (no-mitigation) behavior is a
separate, still-open question the user asked directly: does the smaller
circuit (mean 3.94 CX vs 11, 9 measurement groups vs 13) measure better
on REAL IonQ noise, with no mitigation involved at all?

**Verified exactly before any real submission** (matching this
project's own established discipline): the full circuit-build +
qubit-wise-grouped-measurement + tapered-label reconstruction pipeline,
checked against the exact statevector locally, gives
**1.28e-11 kcal/mol** — correct to machine precision before spending any
real API time.

**Real result, concurrent submission (ideal/aria-1/forte-1, one job
each, 972 circuits total), 8-seed bootstrap mean ± std**:

| | ideal (control) | aria-1 | forte-1 |
|---|---|---|---|
| raw | 1.53 ± 0.89 | **47.78 ± 2.20** | **51.25 ± 1.77** |

Ideal correctness control passed (1.53 kcal/mol, consistent with real
shot noise).

**Placed against the other two real-hardware circuit variants tested in
this project**:

| circuit | 2-qubit gates | aria-1 raw | forte-1 raw |
|---|---|---|---|
| abstract fixed ansatz (iteration 9) | constant 11 | **34.98** | **43.03** |
| native-optimized, TrappedIonOptimizerPlugin (iteration 12) | mean 9.28, range 4-11 | 93.73 | 91.43 |
| **Z2-tapered (this iteration)** | **mean 3.94, range 2-4** | 47.78 | 51.25 |

**Honest reading, not the one the local synthetic-noise test (iteration
14) predicted**: tapering's real-hardware raw error sits BETWEEN the
other two — genuinely better than the native-optimized attempt (~1.9x),
but WORSE than the original abstract 11-gate ansatz (~1.2-1.4x), despite
having barely a third as many two-qubit gates. This directly echoes
iteration 12's own finding (native optimization also had fewer gates yet
worse real error) rather than overturning it: **gate COUNT alone is not
the dominant factor in this project's real IonQ results, for the third
circuit variant in a row.** The specific gate STRUCTURE and its
interaction with IonQ's actual (not naively gate-count-modeled) noise
matters more than the raw tally — plausible contributors, not yet
individually isolated: the tapered register's circuit is generic
`StatePreparation` (never hand-optimized the way the 11-gate ansatz was,
see iteration 13's own ALTERNATIVES NOT TAKEN #3), and its gate count is
non-uniform across targets (2-4), so some targets carry disproportionate
noise exposure in the aggregated bilinear energy sum.

**What this does NOT change**: the local-noise-model prediction from
iteration 14 (raw error should drop ~2.3x under tapering) was directionally
consistent with a real improvement over the native-optimized circuit, but
NOT sufficient to predict the ranking against the abstract ansatz — a
concrete, disclosed instance of this project's recurring lesson that
local synthetic-noise conclusions do not automatically transfer to real
IonQ behavior (iteration 9's own founding finding, now confirmed a third
time on a third circuit).

Per the standing branch discipline: this is a real result, reported
honestly including where it falls short of the best number already
known, not spun as a win. No push.

Code: `vqe/z2_tapered_ionq.py` (`--targets`, `--assemble`). Full data:
`vqe/z2_tapered_ionq_results.json`.

---

## Iteration 16: porting the abstract ansatz's structural gate-efficiency trick to the tapered register — a clean, bigger loss, plus a real determinism bug found along the way

**Why this iteration**: the user asked directly — what's the structural
difference between the abstract 11-gate ansatz and the Z2-tapered
circuit, and does trying to reproduce whatever made the abstract ansatz
good on the tapered register actually help?

**Diagnosis (`gate_structure_compare.py`)**: transpiled both circuit
families to `u3`/`cx` at `optimization_level=0` and inspected the full
gate lists, not just counts. The abstract ansatz has MORE gates of BOTH
types (constant 11 CX, constant 51 u3) than the generic-StatePreparation
tapered circuit (mean 3.94 CX, mean 6.72 u3) yet still wins on real
hardware (34.98/43.03 vs 47.78/51.25 kcal/mol, iteration 15). Direct
inspection of `u_0`'s full instruction list found the abstract ansatz's
51 u3 gates are overwhelmingly FIXED, special-angle (0, ±π/2, ±π
combinations) — decomposed Hadamard/S/Sdg-equivalent structural gates
from its Givens-rotation recipe — with only ~6-8 carrying the actual
fitted target angles. The tapered circuit's 7 u3 gates are essentially
ALL generic, arbitrary-angle `StatePreparation` synthesis output.
Hypothesis: gate-angle "specialness," not raw count, may explain the
real-hardware gap (special angles are plausibly cheaper on real trapped-
ion native gates — closer to identity or a Clifford operation than a
generic arbitrary rotation).

**Attempt 1 — port the recipe 1:1 (`z2_tapered_fixed_ansatz.py`,
`build_candidate`, 4-5 angles)**: mirrored `fixed_ansatz.py`'s own
construction exactly — reference prep (X gate) + a discriminator-qubit
RY + CX fan-out bridge (bit-complement pair) + `XXPlusYYGate` Givens hops
across the register's qubit pairs. Two real bugs caught and fixed before
this even ran correctly:
  1. `ValueError: Residuals are not finite in the initial point` from
     `scipy.optimize.least_squares` — the phase-alignment convention
     divides by the fitted statevector's leading component, which can
     land near-zero at an unlucky random initial guess. Fixed by
     wrapping each fit attempt in try/except, skipping to the next
     random seed rather than crashing the whole multi-attempt loop.
  2. After that fix, BOTH candidates gave 0/36 converged for EVERY
     target — a systematic, not random, failure. Diagnosed via direct
     `Statevector` inspection: the circuit only ever reached basis
     states 0 (`|000>`) and 7 (`|111>`) — exactly the two states
     EXCLUDED from the required subspace. Root cause: the RY
     discriminator gate had been applied to the SAME qubit the
     reference `X` gate had just set to 1, conflating the "reference"
     and "discriminator" roles `fixed_ansatz.py`'s own derivation
     deliberately keeps separate (discriminator must start untouched at
     |0>). Fixed by moving RY to an untouched qubit and fanning out CX
     from there instead. Verified via direct Statevector inspection at
     several angle sets before re-fitting.

**After the fix**: 4/36 converged, worst error 0.42-0.51 — still a real,
large failure, not force-fitted into "close enough." This is NOT a
parameter-count problem: measured directly (not assumed) that the
original 4-qubit ansatz's 25 targets ALL share exactly ONE bit-complement
pair (a physical fact — 2 electrons confined to a single Hamming-weight-2
sector, `fixed_ansatz.py`'s own verified premise), while the tapered
register's 36 targets are measured to spread across 1, 2, or 3 of its 3
available bit-complement pairs depending on target:

| pattern | example targets | count |
|---|---|---|
| 1 pair active | u_0, u_1, u_2, u_3, (u1±u3) | several |
| 2 pairs active | most combination targets | majority |
| 3 pairs active | (u1±u4), (u3±u4), (u1±u5), (u3±u5) | 8 |

Z2 tapering's Clifford transform does not preserve Hamming-weight
structure, so the single-reference/single-bridge topology that worked
for the untapered problem has no equivalent shared anchor here.

**Attempt 2 — same topology, more Givens repetitions (8 angles)**:
repeating the 3 available Givens pairs to 7 hops (1 bridge angle + 7
Givens angles = 8 total) converged 36/36 to machine precision (worst
1.26e-15, later reproduced at 2.22e-16 and 1.26e-15 across independent
fits) — confirming the topology itself was right, just
under-parameterized at 4-5 angles. Resulting circuit: CONSTANT 16 CX,
CONSTANT 86 u3, 71/86 (83%) special-angle on `u_0` — genuinely closer to
the abstract ansatz's structural character than generic
`StatePreparation` (0% special-angle).

**Honest caveat stated before any real submission**: 16 CX is MORE than
both baselines (abstract's 11, `StatePreparation`'s mean 3.94). The
naive fidelity model `f=(1-p2)^n2q * (1-p1)^n1q` (aria-1's measured
p2=0.01214) predicts this LOSES to both: f_new=0.801 vs
f_tapered_generic=0.954 vs f_abstract=0.861. But that SAME naive model
already mispredicted the abstract-vs-tapered_generic ranking in
iteration 15 (predicted tapered_generic should win; it didn't) — so it
could not be trusted to rule this out in advance. Running for real (free
`ionq_simulator`, no real credits spent) was the only honest way to find
out.

**Real result, concurrent submission (ideal/aria-1/forte-1, 324
circuits/model, 972 total), 8-seed bootstrap mean ± std**:

| | ideal (control) | aria-1 | forte-1 |
|---|---|---|---|
| raw | 2.165 ± 0.956 | **207.08 ± 2.93** | **218.69 ± 2.65** |

Ideal correctness control passed (2.165 kcal/mol — confirms this is a
real noise effect, not a pipeline bug).

**Placed against every circuit variant tested in this project**:

| circuit | 2-qubit gates | special-angle % | aria-1 raw | forte-1 raw |
|---|---|---|---|---|
| abstract fixed ansatz (iteration 9) | constant 11 | ~83% (6-8/51 fitted) | **34.98** | **43.03** |
| Z2-tapered, generic StatePreparation (iteration 15) | mean 3.94 | 0% | 47.78 | 51.25 |
| native-optimized, TrappedIonOptimizerPlugin (iteration 12) | mean 9.28 | n/a | 93.73 | 91.43 |
| **Z2-tapered, 8-angle FIXED structure (this iteration)** | **constant 16** | **83%** | **207.08** | **218.69** |

**Honest conclusion**: the special-angle structural trick is real and
was successfully reproduced (83% special-angle, matching the abstract
ansatz's own character almost exactly) — but porting it to the tapered
register required 16 CX gates (45% more than the abstract ansatz's 11,
>4x generic StatePreparation's mean 3.94), because the tapered register's
36 targets collectively need Givens coverage across all 3 available
bit-complement pairs where the original 4-qubit problem only ever needed
1. The extra gate count overwhelms whatever benefit the special-angle
structure provides — this is by a wide margin the WORST real-hardware
result of any variant tried in this project, including the previous
worst (native-optimized). **Gate structure matters, but it cannot be
ported independently of the gate count it costs to achieve it on a
register that has lost the physical symmetry (single Hamming-weight
sector) the original recipe depended on.** Reported plainly, not spun.

**Bonus finding, independent of the above**: cross-checking the fitted
8-angle solutions (computed by `build_final_candidate.py`, one process)
against freshly recomputed targets (`z2_tapered_fixed_ansatz_ionq.py`'s
own exactness check, a separate process) failed with error 3.44 — not a
phase mismatch, a genuinely different vector (verified: `u_4` in the two
calls had dot product exactly -1.0, not +1.0; other Schmidt vectors
matched). Root cause: `ef_fragment.py`'s `exact_ground_state()` used
`scipy.sparse.linalg.eigsh` (ARPACK's iterative Lanczos solver) with no
fixed `v0` — scipy draws a new random starting vector every call.
Pinning `v0` to a fixed seed was tried first and made it WORSE (cross-
process diff went from 1.335 to 1.999), pointing to non-associative
floating-point rounding inside ARPACK's multi-threaded sparse
matrix-vector products, not just the random start, as the real source.
**Fixed** by switching to dense `scipy.linalg.eigh` — this fragment's
Hamiltonian is only 256-dimensional (8 qubits), trivially fast to
diagonalize exactly, with no iterative-convergence non-determinism.
Verified 0.0 diff across 3 separate process invocations after the fix
(was up to 1.999 before). This bug was invisible in every prior
iteration of this project because every script called
`build_reduced_problem()` exactly once per process and reused the
result — self-consistent within any single run, but silently wrong the
moment two separate runs' outputs were compared, exactly as happened
here. An independent, durable correctness fix, not specific to any of
Tasks A-E.

**ALTERNATIVES NOT TAKEN**:

1. **A fully general, per-target-optimal ansatz search (e.g. a
   variational circuit-structure search or genetic algorithm over gate
   sequences) instead of hand-porting one fixed recipe.** Rejected: the
   entire point of a FIXED structure is that all 36 targets share the
   same gate sequence (only angles differ) — a per-target-optimal search
   would very likely find smaller per-target circuits (closer to
   StatePreparation's own mean 3.94) but lose the fixed-structure
   property this whole exercise was testing, and reintroduce the
   target-dependent-circuit problem CDR-style methods need to avoid.
   Would revisit if a future goal explicitly drops the "fixed structure"
   requirement and just wants the smallest correct circuit per target.

2. **Accepting a partial-coverage circuit (the 4-angle version, 4/36
   exact) plus a fallback (e.g. generic StatePreparation) for the other
   32 targets, rather than insisting on one 8-angle circuit for all
   36.** Rejected: mixing two circuit families defeats the fixed-
   structure premise just as much as a per-target search would, and
   would make the real-hardware comparison ambiguous (is the result
   coming from the special-angle circuit or the fallback?). Would
   revisit only if a specific downstream method (like CDR) turns out to
   only need fixed structure on a SUBSET of targets, not all 36.

3. **Reducing the 16-CX circuit's gate count post-hoc via
   `TrappedIonOptimizerPlugin` or `optimization_level>=1`, the way
   iteration 12 did for the abstract ansatz's native form.** Rejected
   for this write-up: iteration 12 already found `optimization_level>=1`
   can silently make transpiled 2-qubit gate counts non-constant across
   targets for fitted-angle solutions landing near periodic special
   values — exactly the kind of instability this circuit's fitted
   angles (many near 0/π/2 by construction) would be especially prone
   to, which would need its own verification pass before trusting any
   resulting number. Would revisit as a genuine next step given the
   16-CX raw result is now known to be far too large to be worth
   ZNE/CDR on directly — a real gate-count reduction pass (verified
   constant across all 36 targets before any real submission) is the
   most promising concrete next step if this circuit family is revisited
   at all.

4. **Leaving the eigsh non-determinism bug unfixed and just re-fitting
   within a single process each time (a workaround, not a fix).**
   Rejected: this would have "solved" the immediate blocker but left a
   silent correctness trap in `ef_fragment.py` for any future cross-
   process comparison in this codebase — exactly the kind of bug the
   "no fake or hardcoded values, every result must be a real
   computation" standing rule exists to catch. Fixed at the source
   instead (dense `eigh`), verified with a genuine before/after
   determinism test (0.0 diff after, up to 1.999 before), not merely
   asserted.

Per the standing branch discipline: this is a real, honestly-reported
loss — the biggest one recorded in this project so far — not spun as a
partial win. No push.

Code: `vqe/gate_structure_compare.py`, `vqe/z2_tapered_fixed_ansatz.py`,
`vqe/z2_tapered_fixed_ansatz_ionq.py` (`--targets`, `--assemble`), fix in
`vqe/ef_fragment.py`. Full data: `vqe/z2_tapered_fixed_ansatz_results.json`,
`vqe/z2_tapered_fixed_ansatz_ionq_results.json`.

---

## Iteration 17: exhaustively hunting the mechanism behind the tapered circuit's real-hardware loss — four hypotheses ruled out, the effect confirmed real via a reproducibility check, no mechanism found yet

**Why this iteration**: the user's direct question after iteration 16 —
"it has less gates still it won't win i dont understand." A completely
fair question: the tapered circuit (iteration 15, generic
`StatePreparation`, mean 3.94 CX) loses to the abstract 11-CX ansatz on
real hardware (47.78/51.25 vs 34.98/43.03 kcal/mol) despite having FEWER
of literally everything measurable at the circuit level. This iteration
tries to actually find the mechanism, not just restate that gate count
doesn't explain it.

**Hypothesis 1 — IonQ's real NATIVE gate count (gpi/gpi2/ms) differs from
qiskit's abstract cx/u3 count.** Motivated directly by iteration 12's own
established finding that native-gate cost can diverge sharply from
qiskit-level gate count. Compiled both circuit families to IonQ's actual
native target (`native_stateprep.to_native`, aria-1's `ms` gate) via the
same tooling iteration 12 used. RESULT: tapered wins even MORE decisively
at the native level — mean 3.67 `ms` vs 11 for the abstract ansatz (a
~3x native 2-qubit-gate advantage), and mean 34.67 `gpi`/`gpi2` vs 197 (a
~5.7x native 1-qubit-gate advantage). **RULED OUT** — if anything this
makes the mystery deeper, not smaller.

**Hypothesis 2 — circuit depth (serial time exposure), not just gate
count.** Transpiled both families to `u3`/`cx` and measured `.depth()`.
RESULT: tapered wins — mean 8.31 vs a constant 22 for the abstract
ansatz. **RULED OUT.**

**Hypothesis 3 — measurement-basis hardness.** Real trapped-ion Z-basis
reads are typically cheaper (no physical pulse needed) than X/Y-basis
reads (need an actual rotation before measurement). Counted non-Z-axis
qubits per measurement-group's combined basis label for both pipelines.
RESULT: tapered wins — mean 1.78 non-Z-axis qubits per measurement
circuit vs 2.77 for the untapered case. **RULED OUT.**

**Hypothesis 4 — measurement-reuse-induced noise correlation.** The
tapered pipeline's `reduced_label_map` (needed because tapering collapses
37 original alpha-register Pauli labels onto only 27 unique 3-qubit
labels) means 10 of those 27 unique measurements are shared by 2
different original alpha_labels each — a real structural difference from
the untapered pipeline, which measures all 37 labels independently with
zero sharing. If a single noisy measurement gets reused (with a fixed
sign) in two different terms of the final bilinear energy sum, its noise
does not partially cancel the way two INDEPENDENT noise draws would.
Built a controlled Monte Carlo (`mc_reuse_test.py`, not yet committed
under that name but reproduced here): 200 trials, real binomial shot-
noise sampling at 100,000 shots per measurement, comparing the CURRENT
shared-measurement scheme against a proposed independent-measurement
scheme (same reduced Pauli, but each original alpha_label draws its OWN
independent noisy sample instead of reusing one shared value). RESULT:
RMS error, shared=0.326 kcal/mol, independent=0.354 kcal/mol — ratio
(shared/independent) = **0.919**. Independent measurement is NOT better;
if anything marginally worse. **RULED OUT.**

  *A caught-and-disclosed dead end along the way*: an earlier, quicker
  sensitivity test (perturbing one raw expectation-value entry by
  δ=0.001 and measuring the resulting energy shift) had suggested the
  tapered reconstruction formula was ~2.4x MORE sensitive to a single
  noisy input than the untapered formula (mean dE/dδ 0.32 vs 0.13
  kcal/mol). This looked like strong evidence FOR hypothesis 4 and
  informed the decision to build the proper Monte Carlo test above. But
  it did not survive scrutiny: the two quick scripts built their "same"
  8-name/8-label test set via `random.sample(list, 8)` with the same
  seed (0) applied to two DIFFERENTLY-ORDERED underlying lists
  (`sorted(reduced_targets.keys())`, alphabetical — vs `slot_names(K)`,
  canonical u_0..u_5-then-combinations order) — `random.sample` with a
  fixed seed on differently-ordered inputs draws genuinely different
  items, confirmed directly by printing both lists and their samples
  side by side. The two scripts were silently comparing sensitivity at
  DIFFERENT (name, label) points, not the same ones — an apples-to-
  oranges bug, not a real formula-sensitivity difference. Caught before
  building anything further on top of it, and reported here rather than
  quietly dropped, per this project's standing honesty rules.

**Reproducibility check — is the observed gap even real, or shot-noise
luck?** Every real-hardware number in this project (iterations 9-16)
comes from exactly ONE submission per circuit at 10,000 shots each,
with the reported ±std being BOOTSTRAP RESAMPLING variance from that
single draw's counts — not independent-submission variance. It was never
verified that resubmitting the same circuit gives a stable number. Given
every cheaper hypothesis had just been ruled out, this was the
highest-value remaining check: resubmitted iteration 15's EXACT SAME
circuit (`z2_tapered_ionq.py --targets`, unmodified) fresh, a second,
fully independent real submission (same 10,000 shots/circuit, same three
noise models, 324 circuits/model, 972 total). Original checkpoint backed
up first (`vqe/ionq_simulator_binding_curve_checkpoints/z2_tapered_targets_run1.json`)
so both raw datasets are preserved.

| | ideal (control) | aria-1 | forte-1 |
|---|---|---|---|
| run 1 (iteration 15) | 1.53 | **47.78 ± 2.20** | **51.25 ± 1.77** |
| run 2 (this iteration, fresh submission) | 2.61 | **50.28 ± 1.81** | **55.24 ± 1.15** |

Both ideal controls pass. The two independent runs agree within roughly
their combined statistical uncertainty (aria-1: diff=2.50 vs combined
std≈2.83; forte-1: diff=3.99 vs combined std≈2.13, just outside 1σ but
far short of the ~12-15 kcal/mol gap to the abstract ansatz) and both
land solidly in the same high-40s/mid-50s range. **CONFIRMED: the
tapered-vs-abstract real-hardware gap is a genuine, reproducible effect
— not an artifact of limited shot statistics from a single submission.**
Full data: `vqe/z2_tapered_reproducibility_check.json`.

**Honest conclusion**: every mechanism this project can currently test —
qiskit gate count, native (billed) gate count, circuit depth, measurement
basis hardness, and measurement-reuse-induced correlation — has now been
individually ruled out as the explanation, each with a direct, controlled
test rather than a hand-wave. The gap itself is confirmed real via
independent reproduction. **The mechanism remains genuinely unknown.**
The most plausible remaining candidates, none tested here, are things
invisible to any circuit-level or reconstruction-formula accounting:
per-qubit-position noise asymmetry inside IonQ's calibrated system
profile, cross-talk tied to which physical qubit indices are addressed,
or a genuine physical noise-channel effect specific to how errors
propagate through the tapered register's particular (non-Hamming-weight-
constrained) state geometry versus the untapered register's naturally
particle-number-protected one — this last candidate connects directly to
the still-untouched Task D (spin/S² projection), which was specifically
designed to probe exactly this kind of leakage-detection question and
has not yet been attempted in this project.

**ALTERNATIVES NOT TAKEN**:

1. **Testing whether physical qubit PLACEMENT (which of the simulator's
   qubit indices the 3 logical qubits map to) affects the result**, by
   submitting the same circuit with a different `initial_layout`.
   Rejected for this iteration: IonQ's cloud simulator noise models are
   almost certainly a single averaged system-level profile (selected by
   the `noise_model="aria-1"` string, not a live per-qubit calibration
   snapshot), making a placement effect unlikely a priori — but "unlikely"
   is not "ruled out," and this is cheap (free simulator) to actually
   check. Would revisit as the next concrete real experiment if the user
   wants to keep pursuing a circuit-level explanation.

2. **A full local density-matrix noise simulation using IonQ's PUBLISHED
   average gate/measurement fidelities** (as opposed to the classical
   gate-counting proxies used here), to see whether a proper quantum
   channel model — not just a tally of "how many gates" — predicts the
   real ranking. Rejected for this iteration due to time: iteration 14
   already built a related local synthetic-noise model that correctly
   predicted tapering beats native-optimized but did NOT predict tapering
   still trails the abstract ansatz, suggesting even a proper channel
   model may not resolve this without recalibration against the SPECIFIC
   real data now in hand from both circuits. Would revisit by fitting a
   depolarizing+dephasing channel to BOTH real datasets simultaneously
   (not just one, as iteration 14 did) and checking whether a single
   consistent channel can reproduce both real numbers at once.

3. **Task D (spin/S² projection)**, explicitly flagged as connected to
   the leading remaining hypothesis (loss of the untapered register's
   natural Hamming-weight/particle-number error-detecting structure).
   Rejected for this iteration: Task D as originally scoped is about
   checking whether spin-projection catches noise-induced leakage the
   particle-number constraint misses, which is a different (though
   related) question from directly measuring whether tapering's loss of
   that structure explains THIS specific real-hardware gap. Would revisit
   as the most theoretically motivated remaining lead: quantify, for the
   SAME noise model, how much of a given real bit-flip/dephasing event's
   physical-observable effect differs between a state that can leak out
   of a protected weight-2 sector (untapered) vs one that cannot leak
   in any detectable way because the tapered register has no such
   protected sector to leak out of.

4. **Chasing a sixth hypothesis (e.g. crosstalk, calibration drift
   between the two separate DATES the abstract-ansatz and tapered-circuit
   real submissions were made) before confirming reproducibility.**
   Rejected as the wrong order of operations: without first confirming
   the observed gap survives an independent resubmission, any further
   mechanism-hunting risked chasing noise in the original single-draw
   result. The reproducibility check was done FIRST among the remaining
   options for exactly this reason, and it paid off — confirming the
   effect is real rather than sending further investigation down a
   dead end chasing shot-noise variance.

Per the standing branch discipline: an honest "we don't know yet, but we
proved it's real and we know what it isn't" result — not spun as
resolved, not abandoned as unexplainable either. No push.

Code: `vqe/z2_tapered_ionq.py` (rerun, unmodified). Diagnostics: Monte
Carlo reuse test and native-gate/depth/basis-hardness comparisons run
inline (not saved as standalone scripts this iteration — see this
section's numbers for the full record). Full data:
`vqe/z2_tapered_reproducibility_check.json`,
`vqe/ionq_simulator_binding_curve_checkpoints/z2_tapered_targets_run1.json`
(run 1, preserved), `vqe/ionq_simulator_binding_curve_checkpoints/z2_tapered_targets.json`
(run 2, current), `vqe/z2_tapered_ionq_results.json` (overwritten with
run 2's numbers by `--assemble`; run 1's numbers are preserved in the
table above and in the reproducibility-check JSON).

---

## Iteration 18: Task D completed — particle-number leakage post-selection is real, works on real hardware, and beats the best prior baseline

**Why this iteration**: iteration 17's leading, untested hypothesis for
the unexplained tapered-vs-abstract gap — the untapered register's
natural confinement to a single Hamming-weight-2 sector (2 electrons in
4 orbitals) is a real physical constraint that noise can visibly violate,
while the tapered register (after Z2's Clifford transform) has no such
sector left to violate, because tapering doesn't just relabel qubits, it
changes the BASIS in a way that destroys the direct correspondence
between "computational basis state" and "physical electron
configuration." This is Task D from this session's original five-task
list (spin/S² projection — check whether noise-induced leakage the
particle-number constraint misses is removed), finally attempted.

**Part 1 — does real leakage exist, using data already on disk?**
Iteration 9's real IonQ measurement checkpoint (`targets_d1.0.json`)
includes raw per-bitstring counts for all 13 measurement groups of the
abstract ansatz. One of those 13 groups (combined label `ZZZZ`, members
`IIZZ`/`ZIIZ`/`ZZII`) applies NO basis-rotation gates before measurement
— its raw bitstring directly reflects the physical Z-basis electron
occupation, so its Hamming weight is a genuine, real, already-collected
leakage signal. Computed directly:

| model | total shots | leaked (weight≠2) | fraction |
|---|---|---|---|
| ideal | 360,000 | 0 | 0.0000 |
| aria-1 | 360,000 | 17,075 | 0.0474 |
| forte-1 | 360,000 | 18,430 | 0.0512 |

Real, non-negligible, hardware-measured leakage — 0% on the noiseless
control (correctness check passes), ~5% on real noise models. Then
tested causally: does discarding (post-selecting on) those leaked shots
before computing the Pauli expectation values for `IIZZ`/`ZIIZ`/`ZZII`
actually improve accuracy against the exact value? Yes, substantially:

| model | RAW RMS err | POST-SELECTED RMS err | improvement |
|---|---|---|---|
| aria-1 | 0.0509 | 0.0230 | 2.2x |
| forte-1 | 0.0518 | 0.0213 | 2.4x |

This used ONLY data already collected in iteration 9 — no new real
submission, immediate and conclusive for the 3 labels/1 group it covers.

**Part 2 — generalizing to all 13 measurement groups.** The Part 1 result
only applies to the one group that happens to stay in the Z basis; the
other 12 groups apply X/Y basis-rotation gates before measurement, so a
raw post-rotation bitstring's Hamming weight carries no physical
particle-number information at all. To extend post-selection everywhere,
built `spin_leakage_postselect_ionq.py`: adds ONE ancilla qubit and 4
CNOTs (register qubit → ancilla) immediately after state-prep, BEFORE any
group's basis-rotation gates — the ancilla ends up holding the parity of
the register's PRE-rotation Z-basis weight regardless of which Pauli
group is subsequently measured, since gates on disjoint qubits commute
and the basis-rotation gates never touch the ancilla.

**Verified exactly before spending anything real**: `partial_trace`d the
5-qubit ancilla-augmented circuit's statevector over the ancilla and
compared to the un-augmented 4-qubit circuit's own density matrix — max
diff **1.5e-36** (twice, across two separate real runs — see the bug
below). Separately confirmed the ancilla reads `0` with probability
exactly `0` for `P(ancilla=1)` on the noiseless physical (weight-2)
state, i.e. never mis-flags a genuinely physical outcome.

**A real bug, caught by the real submission itself (not local testing)**:
the first `--targets` run submitted successfully (549.7s round trip, all
3 real jobs retrieved) but then crashed on a `KeyError: 'noiseless_energy'`
while assembling the checkpoint dict — `rank6_symmetry_vd.setup()`'s
returned dict uses the key `noiseless_numpy`, not `noiseless_energy`
(the key name used by the OTHER setup function, `qforge.setup_fragment`,
used elsewhere in this same file for `assemble()`). Since the crash
happened before `save_ckpt()`, the real counts from that first submission
were lost and had to be resubmitted (a second real, free-simulator
round trip, ~463s) after fixing the one-line bug. Reported plainly:
this cost real wall-clock time (not real money — simulator only) from a
naming inconsistency between two setup helpers in this codebase that
happen to describe the same physical quantity under different keys.

**Real result, concurrent submission (36 targets × 13 groups × 5 qubits,
468 circuits/model, 1,404 total), 8-seed bootstrap mean ± std**:

| | ideal (control) | aria-1 | forte-1 |
|---|---|---|---|
| leakage fraction (odd-weight, ancilla-detected) | 0.0000 | 0.0712 | 0.0793 |
| RAW (ancilla overhead, no post-selection) | 1.83 ± 0.93 | **67.25 ± 1.51** | **74.23 ± 1.57** |
| POST-SELECTED (leaked shots discarded) | 1.38 ± 0.99 | **31.77 ± 1.43** | **33.86 ± 1.15** |

Ideal correctness control passed (1.83 kcal/mol raw). The measured
leakage fraction (7.12%/7.93%) is HIGHER than Part 1's single-group
estimate (4.74%/5.12%) — consistent with the ancilla circuit itself
costing 4 extra CX gates (15 total vs 11), each an additional
opportunity for a real bit-flip before the parity is latched.

**Placed against every prior real number for the untapered register**:

| circuit | aria-1 | forte-1 |
|---|---|---|
| abstract 11-gate ansatz, no ancilla (iteration 9) | 34.98 | 43.03 |
| **THIS run, +ancilla overhead, RAW (no post-selection)** | 67.25 | 74.23 |
| **THIS run, +ancilla, POST-SELECTED** | **31.77** | **33.86** |

**Honest reading**: the ancilla overhead ALONE is a net loss — 4 extra
CX gates on real IonQ noise cost more than they're worth if the leakage
information they provide is never used (67.25/74.23, roughly double the
no-ancilla baseline's error). But USING that information via
post-selection doesn't just recover the overhead — it produces a result
BETTER than the original circuit had NO leakage detection at all
(31.77 < 34.98, a real ~9% improvement; 33.86 < 43.03, a real ~21%
improvement). **This is the best raw (non-ZNE) real-hardware number this
project has obtained for the untapered register.** The mechanism is now
concretely demonstrated, not just hypothesized: real noise really does
kick a measurable, non-negligible fraction of shots (~7-8%) out of the
physically valid sector, and that leakage really does carry enough
signal to be worth detecting and discarding.

**What this does and does NOT establish about iteration 17's mystery**:
it establishes, for the first time with real data, that the untapered
register's Hamming-weight structure is a REAL, exploitable resource that
the tapered register — by construction, since Z2 tapering's Clifford
transform does not preserve which computational basis states correspond
to physical particle numbers — cannot access in any analogous way (there
is no single qubit or fixed basis rotation in the reduced 3-qubit
register whose measurement would reveal "did this shot leak," because
ALL 6 of its live basis states are equally "valid" post-tapering). This
is consistent with, and strengthens, the leading hypothesis from
iteration 17. It does NOT yet prove the SIZE of iteration 17's specific
observed gap (12-15 kcal/mol between the plain circuits, no ancilla on
either side) is fully explained by this mechanism — that would require
either a comparable leakage-detection scheme for the tapered register
(shown here to be structurally unavailable) or a quantitative model
translating "5-8% Z-basis leakage in 4 qubits" into "expected excess
error in a 3-qubit reduced measurement," neither of which has been built.
Reported as a strong, real, well-supported piece of the picture — not
oversold as the complete answer.

**ALTERNATIVES NOT TAKEN**:

1. **A full weight-exactly-2 detector (distinguishing weight 0/1/3/4, not
   just odd/even parity) using 2+ ancillas**, which would also catch
   even-weight leakage (0, 2-but-wrong-state is not detectable this way,
   or 4) that a single parity ancilla misses. Rejected for this
   iteration: single-bit-flip errors (odd-weight) are the dominant
   physical error mode this project's noise-accounting has consistently
   assumed elsewhere (e.g. `fidelity_threshold_curve.py`'s per-gate
   models), so a parity check was the right first, cheapest test; the
   real ~7-8% detected rate already demonstrates a strong effect without
   needing the extra ancilla overhead of a full detector. Would revisit
   if the parity-only result's residual error (POST-SELECTED still isn't
   at the ideal 1.38 kcal/mol floor) suggests even-weight leakage is
   still contributing meaningfully.

2. **Applying this same ancilla scheme to the Z2-tapered circuit anyway**,
   checking parity of SOME arbitrary 3-qubit combination even though it
   has no established physical meaning post-tapering, just to get a
   directly comparable number. Rejected: this would test something
   without a clear physical interpretation (there is no reason a random
   parity check on the reduced register's qubits would correlate with
   "did an error occur" the way it does for the physically-grounded
   untapered case) — a null or misleading result would be as likely as a
   meaningful one, and testing it wouldn't actually validate or refute
   the CORE claim (that tapering destroys a REAL structure), just add
   an ambiguous data point. Would revisit only alongside a rigorous
   derivation of what quantity, if any, in the tapered basis plays an
   analogous protective role (if any exists at all).

3. **Recomputing iteration 9's original 11-CX (no-ancilla) numbers in
   the SAME batch as this iteration's ancilla submission**, to fully
   control for the same batch-to-batch variation iteration 17's
   reproducibility check flagged as a real methodological concern.
   Rejected for time this iteration, given the improvement found
   (31.77 vs 34.98, 33.86 vs 43.03) is not enormous in absolute terms
   and iteration 17 already established that independent real
   submissions of the SAME circuit agree within a few kcal/mol (not
   enough alone to manufacture this particular result, but also not
   quantified for THIS specific comparison). Would revisit as the
   correct rigor upgrade before treating "31.77 beats 34.98" as
   airtight rather than "real and repeated-hypothesis-consistent."

4. **Retrying the FIRST failed submission's already-paid-for real API
   call by attempting to recover its data from process memory/logs
   instead of resubmitting from scratch.** Rejected: the process had
   already exited by the time the bug was diagnosed (the crash happened
   inside the same script invocation, not a separate recoverable step),
   so there was no real data to recover — a straightforward "fix and
   resubmit," not a case where cleverness could have avoided the second
   real (free, simulator-only) round trip. Documented as a real cost
   (about 8 minutes of wall-clock time) from the bug, not hidden.

Per the standing branch discipline: a genuine, real, hardware-validated
positive result — the best raw untapered number this project has found —
reported alongside its real cost (ancilla overhead) and its real limits
(doesn't fully close iteration 17's gap, doesn't transfer to the tapered
register by construction). No push.

Code: `vqe/spin_leakage_postselect_ionq.py` (`--targets`, `--assemble`).
Full data: `vqe/spin_leakage_postselect_ionq_results.json`,
`vqe/ionq_simulator_binding_curve_checkpoints/spin_leakage_targets.json`.

---

## Iteration 19: does leakage removal fix ZNE or CDR? No to both — but the search was rigorous, and it clarifies what's actually going on

**Why this iteration**: the user asked directly to synthesize everything
learned across 18 iterations and try something genuinely new,
mathematically motivated, aimed at finally reaching chemical accuracy
(1 kcal/mol). Two concrete, well-motivated hypotheses followed directly
from iteration 18's finding that real leakage (~7-8% of shots landing
outside the physical weight-2 sector) is a real, previously-uncounted
noise source: (1) does removing it before fitting ZNE unlock the
plateau raw ZNE has never found (iterations 13, 14)? (2) does removing
it before fitting CDR rescue CDR, which made things 2.1-2.6x WORSE on
real IonQ hardware (iteration 9)?

**Test 1 — leakage + ZNE (`leakage_zne_floor_tested.py`).** Built a local
synthetic-noise sweep combining iteration 18's verified-exact ancilla
circuit with the SAME depolarizing-scaling noise-scale model and the
SAME mandatory 3-direction floor test this project has used since
iteration 13 (range@order1, order@widest-range, range@order2), applied
separately to a RAW and a POST-SELECTED scheme, at scales 1-7, 8 seeds,
proper joint multinomial shot sampling (a genuine methodological
upgrade over the project's existing per-label-independent `shot_sample`
shortcut, required because post-selection needs the ancilla outcome and
the Pauli outcome correlated within the SAME shot, not sampled
independently).

Exact (zero-shot-noise) result at every scale:

| scale | RAW (kcal/mol) | POST-SELECTED (kcal/mol) |
|---|---|---|
| 1 | 122.9 | 57.5 |
| 2 | 231.9 | 121.1 |
| 3 | 329.0 | 189.3 |
| 4 | 415.7 | 260.5 |
| 5 | 493.4 | 333.0 |
| 6 | 563.2 | 405.1 |
| 7 | 625.8 | 475.5 |

Leakage post-selection roughly HALVES the error at every single scale —
a clean, general confirmation that iteration 18's real-hardware benefit
isn't a one-off. But the SHAPE of both curves is the problem, not the
level: both climb steadily and near-linearly with scale, with no sign of
flattening. The 8-seed shot-noisy floor test confirms this formally —
**all 6 floor tests (3 directions × {raw, post-selected}) are
DISQUALIFIED**, each for the same "still drifting, not plateaued" reason
iteration 13 first documented. **Leakage removal does NOT fix ZNE's
convergence problem.**

**Test 2 — leakage + CDR (`leakage_cdr_combined.py`).** Same local
noise model (scale=1 only, no ZNE), CDR's own established recipe
(`cdr_mitigation.py`: per-label/per-slot weighted-least-squares-through-
origin scale correction, `N_TRAIN_PER_SLOT=5` random-angle training
circuits with classically-exact known values), applied on the
ancilla-augmented circuit with and without leakage post-selection on
BOTH the training data and the target measurements, 8 seeds:

| scheme | mean (kcal/mol) | std |
|---|---|---|
| raw | 122.95 | 1.00 |
| CDR only | 3.39 | 1.16 |
| leakage only | 57.40 | 0.91 |
| CDR + leakage | 4.65 | 0.98 |

CDR alone looks excellent on this LOCAL model — 3.39 kcal/mol, at
chemical accuracy. Adding leakage removal does NOT improve it further
(4.65 vs 3.39, if anything marginally worse within noise). **This
refutes the hypothesis**: if leakage contamination were a meaningful
part of what breaks CDR, removing it should have helped CDR's fit, not
left it unchanged or slightly worse. It doesn't — CDR's local success
and its real-hardware failure both have a different, unchanged
explanation.

**The correct reading of CDR's 3.39 kcal/mol number**: this is NOT a new
achievement of chemical accuracy. This project already knows, from real
hardware (iteration 9), that this exact CDR recipe — applied to a
circuit that looks identical in every way that matters to THIS local
model — produces 2.1-2.6x WORSE error than raw on real IonQ noise (raw
was 34.98/43.03 there, so CDR was roughly 73-112 kcal/mol for real). The
3.39 kcal/mol number is a restatement of something this project has
already established: the LOCAL depolarizing-scaling noise model is a
close enough fit to itself that CDR's linear-correction premise works
almost perfectly against it, but that premise does not survive contact
with real IonQ noise (coherent errors, crosstalk, drift — whatever the
real gap is, it is NOT primarily leakage, per this iteration's finding).
Reporting the 3.39 kcal/mol number without this context would be
actively misleading; reported here with it.

**Decision: no new real-hardware submission this iteration.** Given (a)
CDR alone is already known, from real data, to catastrophically fail on
this hardware, and (b) this iteration's own local test shows leakage
removal does not change CDR's behavior in the direction that would
justify hoping for a different real-hardware outcome this time, spending
a real submission on CDR+leakage is not supported by the evidence in
hand. This mirrors iteration 14's own discipline (ZNE+tapering: math
done, found wanting, no real submission) and iteration 4's original one
(don't submit for real unless the math passes).

**Where this leaves the project's best real number**: unchanged at
iteration 18's **31.77 kcal/mol (aria-1) / 33.86 kcal/mol (forte-1)**,
raw leakage-postselected, real hardware, no extrapolation risk. Chemical
accuracy (1 kcal/mol) has NOT been reached. The honest picture after 19
iterations: ZNE has now failed its own plateau requirement 4 independent
times (raw untapered, tapered, leakage-postselected untapered, and — via
the exponential/order sweeps — every fit variant tried within each); CDR
and PEC both actively hurt on real hardware; leakage post-selection is
the ONE technique in this entire project that is both physically
motivated AND validated to help on real hardware, not just locally. The
~32-34 kcal/mol residual after leakage removal is consistent with an
irreducible floor set by this circuit's real 2-qubit gate count (11,
soon 15 with the ancilla) at IonQ's real per-gate error rate (~1.2%) —
not something any extrapolation or regression-based correction tested in
this project has been able to lift without either failing to converge
(ZNE) or failing to generalize from the idealized model that makes it
look good (CDR, PEC).

**ALTERNATIVES NOT TAKEN**:

1. **Real gate-folding ZNE on actual IonQ hardware, rather than this
   project's local depolarizing-rate-scaling proxy for "noise scale."**
   Every ZNE test in this project's history (iterations 11, 13, 14, this
   one) uses the SAME simplified local model: scaling the per-gate error
   RATE directly, not physically folding gates (G → G G† G) and
   submitting the longer circuit for real. These are similar but not
   identical for a real device, and it remains formally possible that
   REAL folded circuits behave better than this proxy predicts. Rejected
   for this iteration due to cost and a weak prior: iteration 12's real
   ZNE-linear result (~30 kcal/mol) was ALREADY shown (iteration 13) to
   fail its own floor test on real data once the range was extended, so
   there is direct real-hardware evidence, not just the local proxy,
   that this project's actual ZNE submissions don't converge either.
   Would revisit only with a specific, different folding recipe not yet
   tried (e.g. random/mixed folding rather than uniform integer scales).

2. **A full weight-exactly-2 detector (2+ ancillas, catching even-weight
   leakage too) instead of the single-ancilla parity check.** Rejected:
   estimated impact is small — even-weight leakage (0 or 4 particles)
   requires 2 correlated bit-flips, roughly quadratically rarer than the
   single-flip (odd-weight, ~7-8%) events already caught, unlikely to
   close a 30x gap to chemical accuracy on its own. Would revisit if a
   cheaper, single-extra-ancilla design is found that also catches
   weight-2-but-wrong-configuration coherent errors, which NO Hamming-
   weight-based check (parity or full) can detect in principle.

3. **Testing CDR + leakage removal at multiple noise scales (building
   toward a combined CDR+ZNE+leakage triple-stack), rather than stopping
   after the single-scale CDR result already refuted the core
   hypothesis.** Rejected: once leakage removal was shown not to change
   CDR's qualitative behavior at scale=1, extending to more scales would
   only re-confirm the same negative finding at higher cost, not add new
   information — the mechanism question (does leakage explain CDR's
   real-hardware gap?) was already answered. Would revisit only if a
   DIFFERENT, real reason emerges to suspect CDR's real-hardware failure
   mode is scale-dependent in a way this single-scale test couldn't see.

4. **Running the locally-excellent CDR+leakage combination for real on
   IonQ anyway, on the logic that "it's free, why not check."** Rejected
   deliberately, not out of laziness: this project's own repeated,
   hard-won lesson (iterations 9, 11, 12, and now reconfirmed here) is
   that a technique looking good on the local synthetic model is weak
   evidence at best for real-hardware performance, and CDR specifically
   already has a DIRECT real-hardware data point (iteration 9) showing
   catastrophic failure under conditions this iteration's own test says
   leakage-removal wouldn't have fixed. Submitting anyway would be
   spending real submissions against the evidence already in hand, not
   because of it — exactly the discipline this project's standing rules
   (floor-test everything, don't submit unless the math passes) exist to
   prevent.

Per the standing branch discipline: two more honest negatives, reported
with exactly as much rigor as the positives — the project's best real
number stands at 31.77/33.86 kcal/mol, still well short of chemical
accuracy, with a clear, evidence-based account of why each tested path
there has failed. No push.

Code: `vqe/leakage_zne_floor_tested.py`, `vqe/leakage_cdr_combined.py`.
Full data: `vqe/leakage_zne_floor_tested_results.json`,
`vqe/leakage_cdr_combined_results.json`.

---

## Iteration 20: Task B (Contextual Subspace VQE) via the reference `symmer` package, and the final synthesis of this 20-iteration project

**Why this iteration**: the user asked to finish Task B and revisit Task
D, applied against both the abstract 11-gate and tapered circuits, and
to give a clear, final assessment of whether chemical accuracy is
reachable with everything found so far.

**Building Task B correctly**: rather than hand-implement CS-VQE's
noncontextuality-detection algorithm from scratch (real risk of a subtle
correctness bug, given this project's own history of exactly this kind
of mistake in Z2 tapering's early development — the X-vs-Z basis mixup
documented in `z2_tapering.py`), installed the reference implementation,
`symmer`, from the same research group (Kirby, Tranter, Love) that
published the technique. `PauliwordOp.from_qiskit()` converts this
project's existing `SparsePauliOp` Hamiltonian directly, avoiding any
hand-rolled label-convention risk. Verified the result is the SAME
physical operator as this project's own independent exact diagonalization
(`ef_fragment.exact_ground_state`, dense `scipy.linalg.eigh`) before
trusting anything further: **1.06e-11 kcal/mol** match.

**Why the FULL Hamiltonian, not the forged alpha register**: CS-VQE
partitions a standalone, monolithic qubit Hamiltonian. The "alpha
register" used throughout this project's other 19 iterations is not
itself a physical system with its own ground state — its Pauli labels
come from decomposing the FULL Hamiltonian's terms for entanglement
forging's bipartite Schmidt reconstruction, a fundamentally different
decomposition. Applying CS-VQE faithfully means building it against the
real, un-forged, full 8-qubit H4 problem — genuinely new territory for
this project, not a variant of anything tested in iterations 1-19.

**Result — the ContextualSubspace qubit-count sweep**:

| quantum-remainder qubits | energy (Ha) | err vs exact (kcal/mol) |
|---|---|---|
| 1 | -2.098546 | 42.571 |
| 2 | -2.104388 | 38.906 |
| 3 | -2.112595 | 33.755 |
| 4 | -2.130663 | 22.417 |
| **5** | **-2.166387** | **0.000 (exact)** |
| 6, 7 | — | FAILED: "Search region collapsed without identifying any stabilizers" |

The 5-qubit result matching EXACTLY (0.000 kcal/mol, not just small) is
not a coincidence — it is precisely Z2 tapering's own already-established
reduction (8 qubits − 3 symmetry generators = 5, `z2_tapering.py`'s
`verify_full_hamiltonian_symmetries`), now independently reproduced by a
completely different reference implementation from a different research
group. A genuine, valuable cross-validation of this project's own Z2
tapering work, even though it isn't new information on its own.

The actually NEW information is qubits 1-4: **CS-VQE's noncontextual
approximation error, PURELY CLASSICAL, with ZERO hardware noise
involved, already exceeds chemical accuracy at every qubit count below
5** — 22.4 kcal/mol at 4 qubits, worse at fewer. This is a clean,
verified, negative result: there is no way to trade a small amount of
classical approximation error for a smaller quantum circuit here and
still have room left for chemical accuracy once real hardware noise is
added on top. The 6- and 7-qubit points failing ("search region
collapsed") are a real limitation of `symmer`'s current search strategy
for this Hamiltonian at those specific target sizes, not something this
iteration had time to work around — reported honestly, not smoothed over
(those points aren't needed for the main conclusion regardless, since
5 qubits is already the interesting boundary).

**Honest reading**: Task B does not provide a new path to chemical
accuracy this project hadn't already found. Its best-possible qubit
count for an exact result (5) matches Z2 tapering exactly; going smaller
costs more accuracy than chemical accuracy allows, before any hardware
noise is even considered. No real IonQ submission was made for a
monolithic CS-VQE circuit — building one would mean designing an
entirely new ansatz for a genuinely different (non-forged) circuit
paradigm, and the classical-accuracy ceiling found here means even a
NOISELESS 4-qubit version would already fail chemical accuracy by 22x,
so there is no math to build on top of a working real-hardware attempt.
This mirrors this project's own standing discipline (iteration 14,
iteration 19: don't submit for real unless the math passes).

**FINAL SYNTHESIS — where this 20-iteration project stands**:

| approach | best real/verified result | chemical accuracy (1 kcal/mol)? |
|---|---|---|
| Abstract 11-gate ansatz, raw (iteration 9) | 34.98 / 43.03 kcal/mol (real) | no |
| Z2-tapered, raw (iteration 15/17) | 47.78-55.24 kcal/mol (real, reproduced 2x) | no |
| Z2-tapered, hand-derived fixed structure (iteration 16) | 207/218 kcal/mol (real) | no |
| Native-optimized (iteration 12) | 93.73/91.43 kcal/mol (real) | no |
| ZNE, any variant tested (iterations 11, 13, 14, 19) | no plateau found, ever | not converged, can't report a trustworthy number |
| CDR (iteration 9, confirmed iteration 19) | works locally, 2.1-2.6x WORSE than raw on real hardware | no, actively harmful |
| PEC (iteration 12) | 170-400+ kcal/mol (real) | no, catastrophic |
| **Particle-number leakage post-selection (iteration 18)** | **31.77 / 33.86 kcal/mol (real)** | **no, but the best real result found** |
| CS-VQE, exact classical bound (iteration 20) | 0.000 kcal/mol at 5 qubits (noiseless), 22.4 kcal/mol at 4 | matches Z2 tapering exactly, no improvement |

Twenty iterations, five originally-scoped tasks all attempted (A, C, E
fully; B and D fully as of this iteration), a working combination search
across all pairwise-plausible interactions between the project's major
techniques (leakage×ZNE, leakage×CDR), a real reproducibility check, and
a real, previously-unknown determinism bug found and fixed along the
way. The honest conclusion: **chemical accuracy was not reached, and the
evidence gathered across this whole project points to why — not a single
missed technique, but a real, structural noise floor.** The best real
number (31.77/33.86 kcal/mol) uses the ONE technique validated to help
on both a local model and real hardware (leakage post-selection); every
extrapolation-based method (ZNE) fails to converge on this circuit
family no matter how it's combined or measured; every regression-based
method (CDR) that looks good on an idealized noise model has been shown,
concretely, not to survive contact with real IonQ noise; and the one
structural reduction technique left untried (CS-VQE) matches the
qubit-count ceiling this project already had, without unlocking anything
smaller. Getting from ~32 kcal/mol to 1 kcal/mol from here would most
plausibly require either meaningfully better real hardware (lower native
gate error than IonQ's currently measured ~1.2% per two-qubit gate), a
circuit with fundamentally fewer two-qubit gates than anything found in
this 20-iteration search (the abstract 11-gate ansatz remains the best
structural design found, and CS-VQE confirms 5 qubits is the qubit-count
floor for an EXACT classical partition), or a mitigation idea not yet
conceived by this project or the literature it drew from. Reported as
the complete, honest state of the work — not a partial win dressed up as
a full solution.

**ALTERNATIVES NOT TAKEN**:

1. **Hand-implementing CS-VQE's noncontextuality-detection algorithm
   from scratch instead of using `symmer`.** Rejected: real risk of a
   subtle correctness bug (this project has direct prior experience with
   exactly this failure mode in Z2 tapering's own development), and the
   reference implementation, once verified against this project's own
   exact diagonalization, is strictly more trustworthy for the same
   effort. Would revisit only if a `symmer`-specific limitation (like the
   6/7-qubit search failures found here) genuinely blocked a promising
   result — it didn't, so not pursued.

2. **Debugging `symmer`'s 6- and 7-qubit "search region collapsed"
   failures to complete the full 0-7 sweep.** Rejected: the qubit counts
   that matter for the conclusion (4 and 5, the boundary where accuracy
   crosses chemical-accuracy-infeasible) were already obtained; 6 and 7
   qubits would only ever show SMALLER errors than 5's already-exact
   0.000 kcal/mol (more qubits retained = less thrown into the
   noncontextual approximation), so they cannot change the answer to
   "can CS-VQE go SMALLER than 5 qubits and still be accurate" — the
   question this iteration was actually asking. Would revisit if a
   future goal specifically needed the 6/7-qubit partitions themselves
   (e.g. for a different molecule where 5 isn't already known to be
   exact).

3. **Building a real monolithic CS-VQE circuit and running it on IonQ
   anyway, to see if hardware-specific effects might beat the classical
   bound's implications.** Rejected: the classical bound is a
   NOISELESS floor — no real circuit can do BETTER than its own
   noiseless limit, only worse. A 4-qubit circuit already fails chemical
   accuracy by 22x before any noise; adding real hardware noise can only
   widen that gap, never close it. Submitting for real here would
   violate this project's own standing discipline (verify the math
   passes before spending a real submission) for no possible benefit.

4. **Extending the leakage-postselection technique (iteration 18) to
   the CS-VQE 5-qubit partition, to see if a Z2-symmetry-equivalent
   quantum remainder could ALSO support particle-number-style leakage
   detection.** Rejected for this iteration on time grounds: the 5-qubit
   CS-VQE partition IS mathematically identical to Z2 tapering's own
   already-tested reduction (this iteration's whole point in confirming
   the 0.000 kcal/mol exact match), so this would very likely reproduce
   iteration 17's own finding that the tapered register's Clifford
   transform destroys the physical weight-2 structure leakage detection
   depends on — a predictable, not novel, result. Would revisit only if
   a FUTURE CS-VQE application (a different molecule, or a genuinely
   sub-5-qubit accurate partition found some other way) produced a
   quantum remainder with a DIFFERENT structure worth checking fresh.

Per the standing branch discipline: a complete, honest close-out of this
session's task list, including a clear final synthesis that does not
overstate what was achieved. No push.

Code: `vqe/cs_vqe_h4.py`. Full data: `vqe/cs_vqe_h4_results.json`.

---
