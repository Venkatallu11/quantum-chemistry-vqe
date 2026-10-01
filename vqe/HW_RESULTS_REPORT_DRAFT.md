# H4 Hardware Benchmark — Results Report (DRAFT, not yet run)

**Status:** plan submitted for approval. No hardware data has been collected for
this run yet. Every value marked `[ ]` will be filled only from the real
hardware data, by the analysis script, after the run.

## 1. What was run

| Item | Value |
|---|---|
| Backend | Forte Enterprise 1 *(to confirm)* |
| Job ID | `[ ]` |
| Date submitted / completed | `[ ]` / `[ ]` |
| Circuits | 51 new + 3 reused from 10 Sep 2026 (job `01a08910-7a2b-762b-b0ad-6207191241b6`) = 54 |
| Design | GC0 + GC1 on all 21 slots (42); GC2 + GC3 on panel u_0, u_1, u_2, (u1+u3), (u2+u5), (u4+u5) (12, of which 3 reused) |
| Shots | 1,100 per new circuit (reused circuits: 2,000) |
| Cost | estimated $2,259.00 · actual billed `[ ]` · difference `[ ]` |
| Code version | analysis locked at commit `c95817c` or later pre-run commit `[ ]`; no changes after data arrives |

## 2. Data-quality checks (before any analysis)

| Check | Pass condition | Result |
|---|---|---|
| All 51 new circuits returned | 51/51 | `[ ]` |
| Postselection keep rate | simulator predicted ~89% | `[ ]` |
| Reused circuits still valid | gate counts match today's code | Pass (checked 25 Sep) |
| Schmidt sign-convention gate (ideal data at exact frame) | residual < 1e-3 | Pass (3.2e-5) |

## 3. Main result — energy error vs exact (kcal/mol; chemical accuracy < 1.0)

| Analysis | Hardware result ± uncertainty | Simulator predicted (forte-1) | Agrees? |
|---|---|---|---|
| Raw | `[ ] ± [ ]` | 26.3 ± 4.9 | `[ ]` |
| Corrected, no shared frame | `[ ] ± [ ]` | 5.6 ± 3.3 | `[ ]` |
| Corrected, shared frame | `[ ] ± [ ]` | 0.30 ± 0.22 | `[ ]` |

Uncertainty: bootstrap resampling of the measured counts. Simulator predictions
from task71 / task72 (noise-model data bootstrapped to 1,100 shots, 4 trials).

## 4. Fit quality

| Quantity | Hardware | Simulator predicted |
|---|---|---|
| Shared-frame chi²/dof | `[ ]` | 1.03 – 1.46 |
| Labels excluded (correction ratio \|B/A\| < 0.01) | `[ ]` of ~516 | 23 – 24 |

- chi²/dof ≈ 1: residuals consistent with shot noise; hardware noise is captured by the noise model.
- chi²/dof ≫ 1: hardware has error the model does not include — reported as such.

## 5. Where exact-state information enters each analysis

- **All three:** circuits prepare the FCI Schmidt states and pair
  combinations; each slot is fitted inside the span of the exact Schmidt
  vectors; the energy uses the exact Schmidt coefficients and beta-side signs.
  This benchmarks measurement and mitigation on a known state, not a VQE
  optimization.
- **No shared frame, additionally:** per-label correction factors B/A come from
  simulating the intended (exact-state) circuit under the fixed, pre-set noise
  model. For labels with |B/A| < 0.01 the correction forces the value to ≈ 0,
  its exact answer.
- **Shared frame, additionally:** assumes one orthonormal frame (15
  parameters), starts the fit at the exact frame, and depends on the Schmidt
  vector signs (pinned in `qforge.forging.SCHMIDT_SIGN_REFERENCE`).

## 6. Did the mitigation help on real hardware?

- Correction alone: raw `[ ]` → no shared frame `[ ]`
- Shared frame: `[ ]` → `[ ]`
- Chemical accuracy reached: `[ yes / no ]`
- Where hardware disagreed with the simulator, and likely reasons: `[ ]`

## 7. Conclusion

`[ Written only after the results are in — positive or negative. ]`

## 8. What the remaining 30 circuits would add

GC2 + GC3 on the other 15 slots, ≈ $1,623 at the same rate card (≈ $1,700 with
margin). This would:

- show directly whether using 54 instead of 84 circuits cost any accuracy;
- give the no-shared-frame analysis full label coverage on every slot.

## 9. Reproducing this report

```
PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py
```

*(analysis script to be written and committed before submission)*
