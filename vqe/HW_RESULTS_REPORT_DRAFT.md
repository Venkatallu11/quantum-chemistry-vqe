# H4 Hardware Benchmark — Results Report (pre-run template)

**Status:** Vadim approved the 81-circuit submission on 11 October 2026. This template is
pre-populated with the fixed design and the latest *simulator-only* references. It is
not evidence that the QPU job has been submitted or that chemical accuracy has been
achieved on hardware. Fill hardware-result fields only from the returned IonQ job data.

## 1. Approved run and accounting

| Item | Fixed value / result |
|---|---|
| Backend | `qpu.forte-enterprise-1` |
| Job ID | **[fill only after accepted submission]** |
| Submitted / completed | [timestamp] / [timestamp] |
| Logical design | 84 measurement settings = 21 slots × 4 general-commuting groups |
| New circuits in this job | 81, in one submission |
| Reused measurements | 3 GC3 circuits (u_0, u_1, u_2), from the 10 Sep 2026 job `01a08910-7a2b-762b-b0ad-6207191241b6` |
| Shots | 1,100 for each of the 81 new circuits; the 3 reused circuits retain their original 2,000-shot counts |
| Cost estimate | $3,875–$3,882, consistent with Vadim's model (~$3,875) |
| Actual billed cost | [fill from IonQ job metadata] |
| Debiasing | Must be `false` in the serialized API job settings; the run call explicitly passes `debiasing=False` and the SDK's `ErrorMitigation.NO_DEBIASING` setting |
| Code revision | Record submission commit SHA here: [SHA] |
| Fixed analysis | `vqe/task73_hw_analysis.py`; bootstrap uncertainty from the same counts |

**Submission protocol:** run the dry-run checks first. Submit exactly once. If submission
or retrieval fails, stop and write to Vadim before any retry. Post the accepted Job ID
to the IonQ thread immediately. Do not resubmit to obtain a more favorable outcome.

## 2. Pre-submission gates and data integrity

Record the actual stdout for each item; do not substitute a pass inferred from a later
energy result.

| Gate | Requirement | Result |
|---|---|---|
| Native ancilla circuit | Equivalent to the verified abstract version within tolerance | [stdout/result] |
| Schmidt-sign gate | Stored ideal data fit the exact frame below the locked threshold | [stdout/result] |
| Reused circuits | Current GC3 circuits reproduce the recorded Task 49 gate counts | [stdout/result] |
| Noiseless circuit check | Every circuit reproduces exact target Pauli expectations within (10^{-6}) | [stdout/result] |
| Cost gate | 81 new circuits × 1,100 shots; estimate in the approved $3,875–$3,882 range | [stdout/result] |
| Debiasing gate | Serialized payload has `settings.error_mitigation.debiasing = false` | [stdout/result] |
| Completeness | 81 new circuit outputs plus the 3 reused settings = 84 required settings | [fill after retrieval] |
| Shot counts | Every new circuit returns exactly 1,100 counts; reused settings use original 2,000-shot counts | [fill after retrieval] |
| Retention / leakage | Report mean, minimum, and any circuits below the preregistered 0.70 flag | [fill after retrieval] |

## 3. Hardware result — headline the frame-free comparison

Report **all three estimators from the same counts**, with the fixed noise correction
and bootstrap procedure. Vadim's instruction is that the frame-free value must be
reported either way; **only the shared-frame estimator may be called “passing.”** Show
signed energy error and absolute error so the sign cannot conceal a large discrepancy.

Exact reference used by the committed pipeline: (-2.166387) Ha.

| Estimator | Hardware energy (Ha) | Signed error (kcal/mol) | Absolute error (kcal/mol) | Bootstrap mean ± SD (kcal/mol) | Bootstrap 5–95% interval | Fit χ²/dof / notes |
|---|---:|---:|---:|---:|---:|---|
| **Corrected, no shared frame — required headline** | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |
| Raw, no correction | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |
| Corrected, shared frame | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |

For the shared-frame fit, also report the number of labels used in χ²/dof and how many
near-zero correction-ratio labels ((|B/A| < 0.01)) were excluded. State whether the
residuals look consistent with shot noise (χ²/dof near 1) or reveal model mismatch
(χ²/dof substantially above 1). Do not treat a large χ²/dof as merely a cosmetic
statistic; it indicates the assumed correction model does not explain all the residuals.

Use the same bootstrap resampling of counts for raw, frame-free, and shared-frame
results. Record the seed, number of resamples, bootstrap mean, SD, and 5th/95th
percentiles. The recommended analysis command is:

```bash
PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py \
  --counts vqe/task73_hw_counts.json --design full --n-boot 1000
```

## 4. Frozen noise-model comparison (simulator only, not hardware)

These are the three fresh Task 77 simulator submissions of the full 84-circuit design,
each at 1,100 shots/circuit. They are the pre-existing comparison set, not measured-QPU
results. With only three independent submissions, their sample ranges are descriptive;
they do **not** establish a reliable Q95.

| Noise model | Estimator | Trial 1 | Trial 2 | Trial 3 | Mean ± sample SD | Range |
|---|---|---:|---:|---:|---:|---:|
| aria-1 | Raw | 20.964 | 27.189 | 30.233 | 26.129 ± 4.725 | 20.964–30.233 |
| aria-1 | No shared frame | 0.899 | 2.055 | 0.504 | 1.153 ± 0.806 | 0.504–2.055 |
| aria-1 | Shared frame | 0.488 | 0.449 | 0.466 | 0.468 ± 0.020 | 0.449–0.488 |
| forte-1 | Raw | 32.025 | 24.780 | 20.374 | 25.726 ± 5.883 | 20.374–32.025 |
| forte-1 | No shared frame | 1.258 | 1.138 | 2.458 | 1.618 ± 0.730 | 1.138–2.458 |
| forte-1 | Shared frame | 0.225 | 0.640 | 0.204 | 0.356 ± 0.246 | 0.204–0.640 |

All values are errors in kcal/mol against the exact H4 reference. The no-frame Forte
trial reaching **2.458 kcal/mol** is a real warning from the simulator rehearsal, not an
outlier to omit. The measured hardware result should be compared directly with this
spread and with the corresponding fit-quality statistics. Record any disagreement rather
than retuning the fixed corrections to force agreement.

The hardware shot design is not perfectly identical to the simulator trials: three GC3
measurements are reused from the 10 September hardware job at 2,000 shots, while the
fresh simulator trials generated all 84 settings at 1,100 shots. Keep this difference
visible when interpreting variation.

## 5. What “passing” means — and what this experiment does not demonstrate

- Apply a pass/fail label **only to the shared-frame estimator** and only relative to
  the explicitly stated threshold. Report the conventional 1 kcal/mol benchmark and
  the project's stricter 0.25 kcal/mol target separately if both are discussed.
- Report the frame-free corrected value regardless of whether it is below 1 kcal/mol.
  It is the required headline for evaluating the correction without joint-frame fitting.
- This is an **exact-state-targeted entanglement-forging measurement/mitigation benchmark**,
  not a blind VQE ground-state search. There is no variational optimization in this run.
- For all three estimators, state preparation targets the classically computed exact FCI
  Schmidt vectors; each fitted slot is restricted to their span; the energy reconstruction
  uses exact Schmidt coefficients and beta-side signs.
- The no-frame estimator additionally uses per-label correction factors simulated from
  the intended exact-state circuits; labels with nearly zero correction ratios can be
  forced toward their exact near-zero value.
- The shared-frame estimator additionally fits a common Schmidt frame, starts at the
  exact frame, and depends on Schmidt-vector sign conventions. Do not present its low
  error as independent of this exact-state information.

## 6. Hardware versus simulator and conclusions

**Headline:** [summarize the frame-free hardware error next to the frozen three-trial
simulator range; state whether hardware agrees with the rehearsal].

**Raw-to-corrected change:** [fill from Section 3].

**Shared-frame estimate:** [fill from Section 3; label passing only if the stated
threshold is met and report the uncertainty and χ²/dof alongside it].

**Main discrepancy / negative finding:** [fill in; include poor fits, anomalous retention,
missing counts, forced-zero exclusions, or noise-model disagreement].

**Actual spend / balance after job:** [fill from IonQ metadata/dashboard].

Conclusion must say plainly whether hardware measurements support, weaken, or fail to
support the simulator prediction. A negative result is still a valid funded benchmark;
do not tune correction parameters after seeing hardware results or imply that the
simulator values were collected on the QPU.

## 7. Publication and attribution

Vadim approved a preprint without separate approval but requested that the draft be sent
to him **one week before posting** so IonQ-specific statements can be checked. Credit
IonQ Research Credits in the acknowledgement (not co-authorship). Vadim is the media
promotion contact. IonQ said it does not post externally written blog content, so do not
plan on an IonQ blog post. The first two SDK bugs (multi-circuit result parsing and
histogram scaling) were fixed in qiskit-ionq PR #273 on `main`, merged 2 September
2026, but were not in a tagged release as of Vadim's 12 October email. The third issue
is the estimate cost-field mismatch (`estimated_cost` versus API `estimated_total_cost`);
the open PR in progress is [PR #280](https://github.com/qiskit-community/qiskit-ionq/pull/280).
Confirm with Vadim before describing final release status in the preprint.

## 8. Reproduction

```bash
# Pre-submit validation (no QPU spend):
PYTHONHASHSEED=0 python vqe/task74_hw_submission.py

# Submit only after all gates pass; use the exact rounded cost printed by the dry run.
# Vadim requires one attempt only. If anything fails, STOP and write to him.
PYTHONHASHSEED=0 python vqe/task74_hw_submission.py --submit --confirm-cost 3882

# After the one accepted job completes:
PYTHONHASHSEED=0 python vqe/task74_hw_submission.py --fetch <JOB_ID>
PYTHONHASHSEED=0 python vqe/task73_hw_analysis.py \
  --counts vqe/task73_hw_counts.json --design full --n-boot 1000
```

If the dry run prints a rounded estimate other than `3882`, use the value it prints
for `--confirm-cost` only after checking that the estimate remains within Vadim's
approved $3,875–$3,882 range. Do not guess or override a failing pre-submission gate.
