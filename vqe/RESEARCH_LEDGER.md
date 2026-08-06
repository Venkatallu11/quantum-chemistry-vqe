# Research Ledger — H4 forged energy noise mitigation

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
