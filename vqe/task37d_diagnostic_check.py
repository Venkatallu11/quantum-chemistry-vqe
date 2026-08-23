#!/usr/bin/env python3
"""Ad-hoc diagnostic (not part of the committed pipeline): is Task 37D's
near-zero Fisher info for p_readout/delta_zz/delta_gpi2/p_gpi2 a real
finding or an artifact (boundary clipping for p_readout sitting exactly
at its clamp lower bound; finite-difference underflow for the others)?
Directly perturb each by a real, non-infinitesimal amount and see if
predictions actually move."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qforge import setup_fragment, fit_all_targets
from task27c_full_h4_folds import kept_slots_for_K
from phys_constrained_reconstruction import build_P_S
from task37c_pec_corrected_joint_fit import predict_corrected, load_real_blended, split_train_val, VAL_FRACTION

K = 6
p = setup_fragment([0, 1, 2, 3], nelec=4, d=1.0, K=K, strict=True)
non_id_labels = sorted(l for l in p["alpha_labels"] if l != p["identity_label"])
fixed_solutions, n_ok, worst = fit_all_targets(p["targets"], tol=1e-10)
diag, plus, kept = kept_slots_for_K(K)
U_exact = np.asarray(p["u_vecs"]).T
P_S = build_P_S(p["alpha_labels"], U_exact)
real_blended, _ = load_real_blended(kept, non_id_labels)
train_blended, _ = split_train_val(real_blended, kept, VAL_FRACTION, seed=37)

name, l = "u_0", non_id_labels[0]
for lbl in non_id_labels:
    if lbl in train_blended[name]:
        l = lbl
        break
m_pec = train_blended[name][l]
print(f"testing slot={name} label={l} m_pec={m_pec}")

base = (0.0, 0.001, 0.0, 0.0, 0.0)
pred_base = predict_corrected(name, l, m_pec, base, fixed_solutions, non_id_labels)
print(f"base (all-zero-ish) prediction: {pred_base}")

for i, param_name in enumerate(["p_zz", "p_gpi2", "delta_zz", "delta_gpi2", "p_readout"]):
    perturbed = list(base)
    perturbed[i] += 0.01
    pred_pert = predict_corrected(name, l, m_pec, tuple(perturbed), fixed_solutions, non_id_labels)
    print(f"  perturb {param_name} by +0.01: pred={pred_pert}  delta={pred_pert - pred_base:+.8f}")

print("\n-- p_readout specifically, perturbing POSITIVE only (since 0.0 sits at its clamp lower bound) --")
for dp in [0.001, 0.005, 0.01, 0.02]:
    perturbed = list(base)
    perturbed[4] = dp
    pred_pert = predict_corrected(name, l, m_pec, tuple(perturbed), fixed_solutions, non_id_labels)
    print(f"  p_readout={dp}: pred={pred_pert}  delta={pred_pert - pred_base:+.8f}")

print("\n-- delta_zz specifically, larger perturbations --")
for dz in [0.001, 0.01, 0.05, -0.01]:
    perturbed = list(base)
    perturbed[2] = dz
    pred_pert = predict_corrected(name, l, m_pec, tuple(perturbed), fixed_solutions, non_id_labels)
    print(f"  delta_zz={dz}: pred={pred_pert}  delta={pred_pert - pred_base:+.8f}")
