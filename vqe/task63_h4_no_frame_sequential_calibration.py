#!/usr/bin/env python3
"""
Task 63 -- fast sequential no-frame calibration.

Scientific rules:
  * no Schmidt-frame fit
  * no exact-energy selection
  * no new circuits
  * same conditioned PEC + GPi2 correction
  * per-slot independent physical 5-angle fit
  * p_ZZ and p_GPi2 selected only by held-out measurement residuals

Stage 1 scans p_ZZ at the previously selected measurement-only p_GPi2.
Stage 2 scans p_GPi2 at the selected p_ZZ.
Candidate fits use the deterministic spectral initializer only; the final
selected model is refit on all labels with four restarts.
"""

from __future__ import annotations
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37b_h4_noise_model import GPI_REAL_MEAN
from phys_constrained_reconstruction import build_P_S
from task60_vadim_no_frame_h4 import (
    K, GATE_NAME, ZZ_ASSUMED, fit_slot_data_only, split_labels,
    load_pooled_postselected, variance_weights,
)

ZZ_MEAN, ZZ_STD = 0.014593, 0.000124
PZZ_GRID = np.array([ZZ_MEAN + ZZ_STD*s for s in (-3,-2,-1,0,1,2,3)], dtype=float)
GPI2_GRID = np.arange(0.00030, 0.000701, 0.000050, dtype=float)
VAL_FRACTION = 0.30

CACHE = {}


def corrected(postselected, pzz, pg2, kept, labels, fixed):
    out = {n:{} for n in kept}
    for name in kept:
        key=(name,round(float(pzz),8),round(float(pg2),8))
        if key not in CACHE:
            A,B,_,_=analytic_A_and_B_conditioned(
                fixed[name]["angles"], GATE_NAME, float(pzz), GPI_REAL_MEAN,
                float(pg2), 0.0, 0.0, labels)
            CACHE[key]=(A,B)
        A,B=CACHE[key]
        for l in labels:
            den=float(A[l]); ratio=float(B[l]/den) if abs(den)>1e-6 else 1.0
            out[name][l]=float(np.clip(postselected[name][l]*ratio,-1.0,1.0))
    return out


def score(corrected_data, weights, P_S, kept, train, val):
    total=0.0; n=0
    for i,name in enumerate(kept):
        tm={l:corrected_data[name][l] for l in train[name]}
        tw={l:weights[name][l] for l in train[name]}
        fit=fit_slot_data_only(P_S,tm,tw,seed=63000+i,n_restarts=0)
        a=fit.vector
        for l in val[name]:
            P=np.real_if_close(np.asarray(P_S[l])).astype(float)
            r=float(a@P@a)-float(corrected_data[name][l])
            total+=float(weights[name][l])*r*r; n+=1
    return total/max(1,n-len(kept)*(K-1))


def build_full(fits,diag,labels,P_S):
    full={n:{} for n in diag}
    for n in diag:
        for l in labels:
            P=np.real_if_close(np.asarray(P_S[l])).astype(float)
            full[n][l]=float(fits[n].vector@P@fits[n].vector)
    for n in range(K):
        for m in range(n+1,K):
            plus=f"(u{n}+u{m})"; minus=f"(u{n}-u{m})"
            full[plus]={}; full[minus]={}
            for l in labels:
                P=np.real_if_close(np.asarray(P_S[l])).astype(float)
                full[plus][l]=float(fits[plus].vector@P@fits[plus].vector)
                full[minus][l]=full[f"u_{n}"][l]+full[f"u_{m}"][l]-full[plus][l]
    return full


def run(backend, checkpoints):
    p=setup_fragment([0,1,2,3],nelec=4,d=1.0,K=K,strict=True)
    labels=sorted(l for l in p["alpha_labels"] if l!=p["identity_label"])
    fixed,n_ok,_=fit_all_targets(p["targets"],tol=1e-10); assert n_ok==len(p["targets"])
    diag,_,kept=kept_slots_for_K(K)
    P_S=build_P_S(p["alpha_labels"],np.asarray(p["u_vecs"]).T)
    train,val=split_labels(kept,labels,seed=63,val_fraction=VAL_FRACTION)
    post,kept_shots=load_pooled_postselected(checkpoints,kept,backend)
    weights=variance_weights(post,kept_shots)

    historical_pg2={"aria-1":0.0006,"forte-1":0.0004}[backend]

    print(f"\n[{backend}] STAGE 1: p_ZZ scan at p_GPi2={historical_pg2:.7f}")
    stage1=[]
    for pzz in PZZ_GRID:
        c=corrected(post,pzz,historical_pg2,kept,labels,fixed)
        s=score(c,weights,P_S,kept,train,val)
        stage1.append((float(pzz),float(s)))
        print(f"  pZZ={pzz:.8f} heldout_chi2/dof={s:.6g}",flush=True)
    best_pzz=min(stage1,key=lambda x:x[1])[0]

    print(f"\n[{backend}] STAGE 2: p_GPi2 scan at p_ZZ={best_pzz:.8f}")
    stage2=[]
    for pg2 in GPI2_GRID:
        c=corrected(post,best_pzz,pg2,kept,labels,fixed)
        s=score(c,weights,P_S,kept,train,val)
        stage2.append((float(pg2),float(s)))
        print(f"  pGPi2={pg2:.7f} heldout_chi2/dof={s:.6g}",flush=True)
    best_pg2=min(stage2,key=lambda x:x[1])[0]

    final_c=corrected(post,best_pzz,best_pg2,kept,labels,fixed)
    final={}
    for i,name in enumerate(kept):
        final[name]=fit_slot_data_only(P_S,final_c[name],weights[name],seed=64000+i,n_restarts=4)
    full=build_full(final,diag,labels,P_S)
    mats=combine_matrices(full,p["alpha_labels"],p["identity_label"],K)
    E,errs=energy_from_alpha_matrices(
        mats,p["terms"],p["lambdas"],p["enuc"],p["signs"],K,
        exact_energy=p["exact_energy"],noiseless_energy=p["noiseless_energy"])
    out={
        "backend":backend,"p_zz":best_pzz,"p_gpi2":best_pg2,
        "stage1":stage1,"stage2":stage2,"energy_ha":float(E),
        "err_vs_exact_kcal":float(errs["err_vs_exact_kcal"]),
        "strict_chemical_accuracy_0p25":bool(abs(errs["err_vs_exact_kcal"])<=0.25),
        "stretch_0p02":bool(abs(errs["err_vs_exact_kcal"])<=0.02),
        "uses_joint_schmidt_frame":False,"uses_exact_energy_for_selection":False,
        "uses_exact_schmidt_basis":True,"uses_exact_target_preparation":True,
    }
    outpath=os.path.join(os.path.dirname(__file__),f"task63_h4_no_frame_sequential_calibration_{backend}_result.json")
    with open(outpath,"w") as f: json.dump(out,f,indent=2)
    print("\nFINAL INFORMATIONAL ENERGY CHECK")
    print(json.dumps({k:out[k] for k in ("backend","p_zz","p_gpi2","energy_ha","err_vs_exact_kcal","strict_chemical_accuracy_0p25","stretch_0p02")},indent=2))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--backend",choices=["aria-1","forte-1"],required=True)
    ap.add_argument("--checkpoint",required=True)
    args=ap.parse_args()
    paths=[p.strip() for p in args.checkpoint.split(",") if p.strip()]
    for pth in paths:
        if not os.path.exists(pth): raise SystemExit(f"checkpoint not found: {pth}")
    run(args.backend,paths)

if __name__=="__main__":
    main()
