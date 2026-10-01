# IonQ Real-Hardware Results: H4 Measurement Checks (Forte Enterprise 1)

**Six real hardware jobs run on `qpu.forte-enterprise-1`, 25 Aug – 10 Sep 2026.** All values below are computed directly from the stored hardware counts.

**Important:** these jobs were circuit-correctness checks on 1–3 of the 21 slots. A molecular energy in kcal/mol needs all 21 slots measured, so **no real-hardware kcal/mol energy exists yet**. That is what the planned 54-circuit run produces. For reference, the simulator rehearsal (noise-model data, not hardware) predicts 26.3 raw, 5.6 corrected, and 0.30 kcal/mol with the shared frame on the forte-1 noise model.

## 1. Summary by job

Values are raw (no noise correction), postselected on the parity ancilla. Deviation = measured − exact expectation value of the prepared state. z = deviation / shot-noise standard error.

| Task | Date | What was run | Slots | Shots | Kept after postselection | Values measured | Mean \|deviation\| | Mean z² | Max \|z\| |
|---|---|---|---|---|---|---|---|---|---|
| 41 | 2026-08-25 | 1 circuit, first test (100 shots) | u_0 | 100 | 92.0% | 2 | 0.088 | 0.72 | 0.9 |
| 42 | 2026-08-25 | same circuit, 500 shots | u_0 | 500 | 90.0% | 2 | 0.019 | 0.29 | 0.8 |
| 43 | 2026-08-25 | 3 circuits in one batch | u_0 | 2,000 | 90.5% | 7 | 0.023 | 1.66 | 2.0 |
| 44A | 2026-09-09 | baseline circuit | u_0 | 2,000 | 90.8% | 2 | 0.043 | 3.59 | 2.2 |
| 44B | 2026-09-09 | same circuit, 1q gates reduced 120→81 | u_0 | 2,000 | 90.1% | 2 | 0.029 | 1.86 | 1.8 |
| 49 | 2026-09-10 | 3 slots, hardest measurement group (GC3) | u_0, u_1, u_2 | 2,000 | 87.5% | 24 | 0.037 | 3.41 | 3.5 |

**All jobs together:** 39 measured values, mean z² = 2.73, 29/39 within 2σ of exact, 37/39 within 3σ.

## 2. What this shows

- **The circuits run correctly on real hardware.** Measured values sit close to the exact values, with no structural errors (wrong signs, wrong basis, broken ancilla).
- **The leakage check works.** The parity ancilla keeps about 85–92% of shots, in line with the simulator's prediction of about 89%.
- **Real noise is present, as expected.** Mean z² ≈ 2.7 means the scatter is somewhat larger than shot noise alone. That is the uncorrected hardware error the mitigation is designed to remove.
- **The hardest measurement group (GC3, Task 49) works on hardware.** It runs with its 3–5 extra two-qubit gates and gives physically sensible values on all three slots.
- **Cost model confirmed.** Real billed costs matched our cost formula (Task 49: $288.85 billed, formula $288.95).

## 3. Every measured value

| Task | Slot | Pauli label | Measured ± std. error | Exact | Deviation | z |
|---|---|---|---|---|---|---|
| 41 | u_0 | XYYX | -0.152 ± 0.103 | -0.064 | -0.088 | -0.9 |
| 41 | u_0 | IYYI | -0.087 ± 0.104 | -0.000 | -0.087 | -0.8 |
| 42 | u_0 | XYYX | -0.067 ± 0.047 | -0.064 | -0.003 | -0.1 |
| 42 | u_0 | IYYI | +0.036 ± 0.047 | -0.000 | +0.036 | +0.8 |
| 43 | u_0 | XYYX | -0.018 ± 0.023 | -0.064 | +0.046 | +2.0 |
| 43 | u_0 | IYYI | +0.002 ± 0.023 | -0.000 | +0.002 | +0.1 |
| 43 | u_0 | YXXY | -0.016 ± 0.024 | -0.064 | +0.048 | +2.0 |
| 43 | u_0 | IXXI | +0.037 ± 0.024 | -0.000 | +0.037 | +1.6 |
| 43 | u_0 | XXYY | +0.063 ± 0.024 | +0.064 | -0.001 | -0.0 |
| 43 | u_0 | IIYY | +0.024 ± 0.024 | -0.000 | +0.024 | +1.0 |
| 43 | u_0 | XXII | +0.000 ± 0.024 | +0.000 | -0.000 | -0.0 |
| 44A | u_0 | XYYX | -0.116 ± 0.023 | -0.064 | -0.052 | -2.2 |
| 44A | u_0 | IYYI | -0.035 ± 0.023 | -0.000 | -0.035 | -1.5 |
| 44B | u_0 | XYYX | -0.022 ± 0.024 | -0.064 | +0.042 | +1.8 |
| 44B | u_0 | IYYI | +0.016 ± 0.024 | -0.000 | +0.016 | +0.7 |
| 49 | u_0 | IIYY | -0.006 ± 0.024 | -0.000 | -0.006 | -0.3 |
| 49 | u_0 | IYIY | -0.054 ± 0.024 | -0.012 | -0.042 | -1.8 |
| 49 | u_0 | IYYI | +0.001 ± 0.024 | -0.000 | +0.001 | +0.0 |
| 49 | u_0 | XZXZ | -0.018 ± 0.024 | +0.010 | -0.028 | -1.2 |
| 49 | u_0 | XZZX | -0.037 ± 0.024 | -0.000 | -0.037 | -1.6 |
| 49 | u_0 | YIYI | -0.002 ± 0.024 | -0.010 | +0.008 | +0.3 |
| 49 | u_0 | YYII | +0.013 ± 0.024 | +0.000 | +0.013 | +0.5 |
| 49 | u_0 | ZXZX | +0.063 ± 0.024 | +0.012 | +0.051 | +2.2 |
| 49 | u_1 | IIYY | -0.030 ± 0.024 | +0.000 | -0.030 | -1.3 |
| 49 | u_1 | IYIY | -0.045 ± 0.024 | -0.000 | -0.045 | -1.9 |
| 49 | u_1 | IYYI | +0.067 ± 0.024 | +0.000 | +0.067 | +2.8 |
| 49 | u_1 | XZXZ | -0.018 ± 0.024 | -0.000 | -0.018 | -0.7 |
| 49 | u_1 | XZZX | +0.064 ± 0.024 | +0.000 | +0.064 | +2.7 |
| 49 | u_1 | YIYI | +0.012 ± 0.024 | +0.000 | +0.012 | +0.5 |
| 49 | u_1 | YYII | +0.006 ± 0.024 | -0.000 | +0.006 | +0.3 |
| 49 | u_1 | ZXZX | +0.060 ± 0.024 | +0.000 | +0.060 | +2.5 |
| 49 | u_2 | IIYY | -0.042 ± 0.024 | -0.000 | -0.042 | -1.7 |
| 49 | u_2 | IYIY | -0.045 ± 0.024 | -0.002 | -0.043 | -1.8 |
| 49 | u_2 | IYYI | -0.006 ± 0.024 | -0.000 | -0.006 | -0.3 |
| 49 | u_2 | XZXZ | -0.062 ± 0.024 | -0.003 | -0.059 | -2.4 |
| 49 | u_2 | XZZX | -0.084 ± 0.024 | -0.000 | -0.084 | -3.5 |
| 49 | u_2 | YIYI | +0.076 ± 0.024 | +0.003 | +0.073 | +3.0 |
| 49 | u_2 | YYII | +0.063 ± 0.024 | +0.000 | +0.063 | +2.6 |
| 49 | u_2 | ZXZX | +0.036 ± 0.024 | +0.002 | +0.034 | +1.4 |

## 4. Job IDs

- Task 41: `01a03ad5-0399-7684-ada3-6d47d3cd00f4`
- Task 42: `01a03b00-2b0b-751c-8f55-834368f0ec47`
- Task 43: `01a03b1a-82c4-7525-9a3d-dff1f392deb2`
- Task 44A: `01a087d8-66df-717e-be93-6ac376d5c62e`
- Task 44B: `01a087e0-877a-75db-9a86-c8f21cc74b97`
- Task 49: `01a08910-7a2b-762b-b0ad-6207191241b6`

## 5. Next step

The planned 54-circuit run (≈ $2,259, within the $2,352.16 balance) measures all 21 slots, giving the first real-hardware H4 energy in kcal/mol, reported three ways: raw, corrected without the shared frame, and shared frame.
