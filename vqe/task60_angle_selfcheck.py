#!/usr/bin/env python3
"""Unit test for Task 60's five-angle S^5 parameterization."""
import numpy as np

from task60_vadim_no_frame_h4 import angles_to_vector, vector_to_angles, K

rng = np.random.default_rng(60060)
worst = 0.0
for _ in range(2000):
    v = rng.normal(size=K)
    v /= np.linalg.norm(v)
    a = angles_to_vector(vector_to_angles(v))
    err = min(np.linalg.norm(a - v), np.linalg.norm(a + v))
    worst = max(worst, float(err))

print(f"Task60 angle round-trip worst error: {worst:.3e}")
assert worst < 1e-10
