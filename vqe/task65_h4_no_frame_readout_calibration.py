#!/usr/bin/env python3
"""
Task 65 -- no-frame H4 readout-attenuation calibration.

Use the existing physical mitigation stack and the independently fitted
5-angle state per slot, but add only the project's already-documented
readout attenuation model:

    m_true ~= m_raw / (1 - 2 p_readout)^weight(label)

p_readout is NOT known from a dedicated hardware calibration in this project.
Therefore it is selected solely by held-out measurement consistency.

No exact energy is used for parameter selection. Exact energy is printed only
as the final informational validation target.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from qforge import setup_fragment, fit_all_targets, combine_matrices, energy_from_alpha_matrices
from task27c_full_h4_folds import kept_slots_for_K
from task39e_conditioned_correction import analytic_A_and_B_conditioned
from task37c_extended_forward_model import readout_attenuation
from task37b_h4_noise_model import GPI_REAL_MEAN
from phys_constrained_reconstruction import build_P_S
from task60_ionq_no_frame_h4 import (
    K, GATE_NAME, fit_slot_data_only, split_labels,
    load_pooled_postselected,
)

PZZ_BY_BACKEND={"aria-1":0.014593,"forte-1":0.014593}
PG2_BY_BACKEND={"aria-1":0.0006,"forte-1":0.0004}

# Project's weak disclosed prior spans 0--0.02. Focused grid includes the
# separately observed ~0.002 value while retaining the prior endpoints.
P_READOUT_GRID=np.array([0.0,0.001,0.002,0.003,0.005,0.0075,0.010,0.015,0.020],dtype=float)

CACHE={}

def corrected(post,p_ro,backend,kept,labels,fixed):
    pzz=PZZ_BY_BACKEND[backend]; pg2=PG2_BY_BACKEND[backend]
    out={n:{} for n in kept}
    for name in kept:
        key=(name,backend,round(float(p_ro),8))
        if key not in CACHE:
            A,B,_,_=analytic_A_and_B_conditioned(
                fixed[name]["angles"],GATE_NAME,pzz,GPI_REAL_MEAN,
                pg2,0.0,0.0,labels)
            CACHE[key]=(A,B)
        A,B=CACHE[key]
        for l in labels:
            ro=readout_attenuation(l,float(p_ro))
            ro=ro if abs(ro)>1e-8 else 1.0
            den=float(A[l])
            ratio=float(B[l]/den) if abs(den)>1e-6 else 1.0
            out[name][l]=float(np.clip((post[name][l]/ro)*ratio,-1.0,1.0))
    return out

def weights(corrected_data,kept_shots):
    W={}
    for name,vals in corrected_data.items():
        W[name]={}
        for l,m in vals.items():
            n=max(int(kept_shots[name][l]),1)
            var=max(1.0-float(m)**2,1e-4)/n
            W[name][l]=1.0/var
    return W

def heldout_score(data,W,P_S,kept,train,val):
    total=0.0;n=0
    for i,name in enumerate(kept):
        tm={l:data[name][l] for l in train[name]}
        tw={l:W[name][l] for l in train[name]}
        fit=fit_slot_data_only(P_S,tm,tw,seed=65000+i,n_restarts=1)
        a=fit.vector
        for l in val[name]:
            P=np.real_if_close(np.asarray(P_S[l])).astype(float)
            r=float(a@P@a)-float(data[name][l])
            total+=float(W[name][l])*r*r;n+=1
    return total/max(1,n-len(kept)*(K-1))

def build_full(fits,diag,labels,P_S):
    full={n:{} for n in diag}
    for n in diag:
        for l in labels:
            P=np.real_if_close(np.asarray(P_S[l])).astype(float)
            a=fits[n].vector
            full[n][l]=float(a@P@a)
    for i in range(K):
        for j in range(i+1,K):
            plus=f"(u{i}+u{j})"; minus=f"(u{i}-u{j})"
            full[plus]={};full[minus]={}
            for l in labels:
                P=np.real_if_close(np.asarray(P_S[l])).astype(float)
                ap=fits[plus].vector
                full[plus][l]=float(ap@P@ap)
                full[minus][l]=full[f"u_{i}"][l]+full[f"u_{j}"][l]-full[plus][l]
    return full

def run(backend,checkpoints):
    p=setup_fragment([0,1,2,3],nelec=4,d=1.0,K=K,strict=True)
    labels=sorted(l for l in p["alpha_labels"] if l!=p["identity_label"])
    fixed,n_ok,_=fit_all_targets(p["targets"],tol=1e-10);assert n_ok==len(p["targets"])
    diag,_,kept=kept_slots_for_K(K)
    P_S=build_P_S(p["alpha_labels"],np.asarray(p["u_vecs"]).T)
    train,val=split_labels(kept,labels,seed=65,val_fraction=0.30)
    post,kept_shots=load_pooled_postselected(checkpoints,kept,backend)

    rows=[]
    for p_ro in P_READOUT_GRID:
        d=corrected(post,p_ro,backend,kept,labels,fixed)
        W=weights(d,kept_shots)
        s=heldout_score(d,W,P_S,kept,train,val)
        rows.append((float(p_ro),float(s)))
        print(f"[{backend}] p_readout={p_ro:.4f} heldout_chi2/dof={s:.6g}",flush=True)

    best_ro=min(rows,key=lambda x:x[1])[0]
    final_d=corrected(post,best_ro,backend,kept,labels,fixed)
    W=weights(final_d,kept_shots)

    fits={}
    for i,name in enumerate(kept):
        fits[name]=fit_slot_data_only(P_S,final_d[name],W[name],seed=66000+i,n_restarts=4)

    full=build_full(fits,diag,labels,P_S)
    mats=combine_matrices(full,p["alpha_labels"],p["identity_label"],K)
    E,errs=energy_from_alpha_matrices(
        mats,p["terms"],p["lambdas"],p["enuc"],p["signs"],K,
        exact_energy=p["exact_energy"],noiseless_energy=p["noiseless_energy"])

    out={
        "backend":backend,"p_readout":best_ro,"grid":rows,
        "energy_ha":float(E),"err_vs_exact_kcal":float(errs["err_vs_exact_kcal"]),
        "strict_chemical_accuracy_0p25":bool(abs(errs["err_vs_exact_kcal"])<=0.25),
        "stretch_0p02":bool(abs(errs["err_vs_exact_kcal"])<=0.02),
        "p_zz_assumed":PZZ_BY_BACKEND[backend],"p_gpi2_assumed":PG2_BY_BACKEND[backend],
        "uses_joint_schmidt_frame":False,"uses_exact_energy_for_selection":False,
        "uses_exact_schmidt_basis":True,"uses_exact_target_preparation":True
    }
    outpath=os.path.join(os.path.dirname(__file__),f"task65_h4_no_frame_readout_calibration_{backend}_result.json")
    with open(outpath,"w") as f:json.dump(out,f,indent=2)
    print("\nFINAL",json.dumps({k:out[k] for k in ("backend","p_readout","energy_ha","err_vs_exact_kcal","strict_chemical_accuracy_0p25","stretch_0p02")},indent=2))
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--backend",choices=["aria-1","forte-1"],required=True)
    ap.add_argument("--checkpoint",required=True)
    args=ap.parse_args()
    paths=[x.strip() for x in args.checkpoint.split(",") if x.strip()]
    for pth in paths:
        if not os.path.exists(pth):raise SystemExit(f"checkpoint not found: {pth}")
    run(args.backend,paths)

if __name__=="__main__":main()
