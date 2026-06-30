#!/usr/bin/env python3
"""
chem.py — Real pure-Python quantum chemistry engine (numpy + scipy only).
Handles s and p orbitals via the McMurchie-Davidson method.
Everything is computed FROM GEOMETRY. The only tabulated inputs are STO-3G basis
parameters (orbital definitions used by every chemistry code), never answers.
"""
import numpy as np
from scipy.special import hyp1f1

ANG = 1.0/0.529177210903  # Angstrom -> Bohr

# STO-3G basis: element -> list of shells. Each shell = (l, [(exp, coef), ...])
# l=0 is s, l=1 is p. Values extracted from the standard STO-3G basis set.
STO3G = {
 'H':[(0,[(3.42525091,0.15432897),(0.62391373,0.53532814),(0.1688554,0.44463454)])],
 'Li':[(0,[(16.119575,0.15432897),(2.9362007,0.53532814),(0.7946505,0.44463454)]),
       (0,[(0.6362897,-0.09996723),(0.1478601,0.39951283),(0.0480887,0.70011547)]),
       (1,[(0.6362897,0.15591627),(0.1478601,0.60768372),(0.0480887,0.39195739)])],
 'C':[(0,[(71.616837,0.15432897),(13.045096,0.53532814),(3.5305122,0.44463454)]),
      (0,[(2.9412494,-0.09996723),(0.6834831,0.39951283),(0.2222899,0.70011547)]),
      (1,[(2.9412494,0.15591627),(0.6834831,0.60768372),(0.2222899,0.39195739)])],
 'N':[(0,[(99.106169,0.15432897),(18.052312,0.53532814),(4.8856602,0.44463454)]),
      (0,[(3.7804559,-0.09996723),(0.8784966,0.39951283),(0.2857144,0.70011547)]),
      (1,[(3.7804559,0.15591627),(0.8784966,0.60768372),(0.2857144,0.39195739)])],
 'O':[(0,[(130.70932,0.15432897),(23.808861,0.53532814),(6.4436083,0.44463454)]),
      (0,[(5.0331513,-0.09996723),(1.1695961,0.39951283),(0.380389,0.70011547)]),
      (1,[(5.0331513,0.15591627),(1.1695961,0.60768372),(0.380389,0.39195739)])],
}
ZTAB = {'H':1,'Li':3,'C':6,'N':7,'O':8}
# cartesian components for each l
COMPS = {0:[(0,0,0)], 1:[(1,0,0),(0,1,0),(0,0,1)]}

def _fact2(n):
    r=1
    while n>1: r*=n; n-=2
    return r

def _norm(a,l,m,n):
    return (2*a/np.pi)**0.75 * (4*a)**((l+m+n)/2) / np.sqrt(_fact2(2*l-1)*_fact2(2*m-1)*_fact2(2*n-1))

def _E(i,j,t,Qx,a,b):
    p=a+b; q=a*b/p
    if t<0 or t>i+j: return 0.0
    if i==j==t==0: return np.exp(-q*Qx*Qx)
    if j==0:
        return (1/(2*p))*_E(i-1,j,t-1,Qx,a,b)-(q*Qx/a)*_E(i-1,j,t,Qx,a,b)+(t+1)*_E(i-1,j,t+1,Qx,a,b)
    return (1/(2*p))*_E(i,j-1,t-1,Qx,a,b)+(q*Qx/b)*_E(i,j-1,t,Qx,a,b)+(t+1)*_E(i,j-1,t+1,Qx,a,b)

def _boys(n,t):
    return hyp1f1(n+0.5,n+1.5,-t)/(2*n+1)

def _R(t,u,v,n,p,PCx,PCy,PCz,RPC):
    val=0.0
    if t==u==v==0:
        return (-2*p)**n*_boys(n,p*RPC*RPC)
    if t>0:
        if t>1: val+=(t-1)*_R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
        val+=PCx*_R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
    elif u>0:
        if u>1: val+=(u-1)*_R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
        val+=PCy*_R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
    else:
        if v>1: val+=(v-1)*_R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
        val+=PCz*_R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
    return val

def _overlap_prim(a,la,A,b,lb,B):
    S=1.0
    for d in range(3):
        S*=_E(la[d],lb[d],0,A[d]-B[d],a,b)
    return S*(np.pi/(a+b))**1.5

def _kinetic_prim(a,la,A,b,lb,B):
    l,m,n=lb
    t0=b*(2*(l+m+n)+3)*_overlap_prim(a,la,A,b,(l,m,n),B)
    t1=-2*b**2*(_overlap_prim(a,la,A,b,(l+2,m,n),B)+_overlap_prim(a,la,A,b,(l,m+2,n),B)+_overlap_prim(a,la,A,b,(l,m,n+2),B))
    t2=-0.5*(l*(l-1)*_overlap_prim(a,la,A,b,(l-2,m,n),B)+m*(m-1)*_overlap_prim(a,la,A,b,(l,m-2,n),B)+n*(n-1)*_overlap_prim(a,la,A,b,(l,m,n-2),B))
    return t0+t1+t2

def _nuclear_prim(a,la,A,b,lb,B,C):
    p=a+b; P=(a*np.array(A)+b*np.array(B))/p; RPC=np.linalg.norm(P-C)
    val=0.0
    for t in range(la[0]+lb[0]+1):
        for u in range(la[1]+lb[1]+1):
            for v in range(la[2]+lb[2]+1):
                val+=_E(la[0],lb[0],t,A[0]-B[0],a,b)*_E(la[1],lb[1],u,A[1]-B[1],a,b)*_E(la[2],lb[2],v,A[2]-B[2],a,b)*_R(t,u,v,0,p,*(P-C),RPC)
    return val*2*np.pi/p

def _eri_prim(a,la,A,b,lb,B,c,lc,C,d,ld,D):
    p=a+b; q=c+d; alpha=p*q/(p+q)
    P=(a*np.array(A)+b*np.array(B))/p; Q=(c*np.array(C)+d*np.array(D))/q
    RPQ=np.linalg.norm(P-Q); val=0.0
    for t in range(la[0]+lb[0]+1):
     for u in range(la[1]+lb[1]+1):
      for v in range(la[2]+lb[2]+1):
       e1=_E(la[0],lb[0],t,A[0]-B[0],a,b)*_E(la[1],lb[1],u,A[1]-B[1],a,b)*_E(la[2],lb[2],v,A[2]-B[2],a,b)
       for tau in range(lc[0]+ld[0]+1):
        for nu in range(lc[1]+ld[1]+1):
         for phi in range(lc[2]+ld[2]+1):
          e2=_E(lc[0],ld[0],tau,C[0]-D[0],c,d)*_E(lc[1],ld[1],nu,C[1]-D[1],c,d)*_E(lc[2],ld[2],phi,C[2]-D[2],c,d)
          val+=e1*e2*(-1)**(tau+nu+phi)*_R(t+tau,u+nu,v+phi,0,alpha,*(P-Q),RPQ)
    return val*2*np.pi**2.5/(p*q*np.sqrt(p+q))

def build_basis(atoms):
    """atoms = [(element, np.array([x,y,z]) in bohr), ...]. Returns list of basis funcs."""
    bfs=[]
    for el,R in atoms:
        for l,prims in STO3G[el]:
            for comp in COMPS[l]:
                exps=np.array([e for e,_ in prims]); coefs=np.array([c for _,c in prims])
                norms=np.array([_norm(e,*comp) for e in exps])
                bfs.append({'R':R,'l':comp,'exps':exps,'coefs':coefs,'norms':norms})
    return bfs

def _contract(f,bf1,bf2,*rest):
    tot=0.0
    e1,c1,n1,l1,R1=bf1['exps'],bf1['coefs'],bf1['norms'],bf1['l'],bf1['R']
    e2,c2,n2,l2,R2=bf2['exps'],bf2['coefs'],bf2['norms'],bf2['l'],bf2['R']
    if not rest:
        for a,ca,na in zip(e1,c1,n1):
            for b,cb,nb in zip(e2,c2,n2):
                tot+=ca*cb*na*nb*f(a,l1,R1,b,l2,R2)
        return tot

def integrals(atoms_ang):
    atoms=[(el,np.array(xyz)*ANG) for el,xyz in atoms_ang]
    bfs=build_basis(atoms); n=len(bfs)
    S=np.zeros((n,n));T=np.zeros((n,n));V=np.zeros((n,n))
    nuclei=[(ZTAB[el],R) for el,R in atoms]
    for i in range(n):
        for j in range(n):
            S[i,j]=_contract(_overlap_prim,bfs[i],bfs[j])
            T[i,j]=_contract(_kinetic_prim,bfs[i],bfs[j])
            v=0.0
            e1,c1,nn1=bfs[i]['exps'],bfs[i]['coefs'],bfs[i]['norms']
            e2,c2,nn2=bfs[j]['exps'],bfs[j]['coefs'],bfs[j]['norms']
            for a,ca,na in zip(e1,c1,nn1):
                for b,cb,nb in zip(e2,c2,nn2):
                    for Z,C in nuclei:
                        v+=-Z*ca*cb*na*nb*_nuclear_prim(a,bfs[i]['l'],bfs[i]['R'],b,bfs[j]['l'],bfs[j]['R'],C)
            V[i,j]=v
    eri=np.zeros((n,n,n,n))
    for i in range(n):
     for j in range(n):
      for k in range(n):
       for l in range(n):
        s=0.0
        for a,ca,na in zip(bfs[i]['exps'],bfs[i]['coefs'],bfs[i]['norms']):
         for b,cb,nb in zip(bfs[j]['exps'],bfs[j]['coefs'],bfs[j]['norms']):
          for c,cc,nc in zip(bfs[k]['exps'],bfs[k]['coefs'],bfs[k]['norms']):
           for d,cd,nd in zip(bfs[l]['exps'],bfs[l]['coefs'],bfs[l]['norms']):
            s+=ca*cb*cc*cd*na*nb*nc*nd*_eri_prim(a,bfs[i]['l'],bfs[i]['R'],b,bfs[j]['l'],bfs[j]['R'],c,bfs[k]['l'],bfs[k]['R'],d,bfs[l]['l'],bfs[l]['R'])
        eri[i,j,k,l]=s
    enuc=0.0
    for x in range(len(nuclei)):
        for y in range(x+1,len(nuclei)):
            enuc+=nuclei[x][0]*nuclei[y][0]/np.linalg.norm(nuclei[x][1]-nuclei[y][1])
    return S,T,V,eri,enuc

def rhf(S,T,V,eri,enuc,nelec):
    Hc=T+V; ev,U=np.linalg.eigh(S); X=U@np.diag(ev**-0.5)@U.T
    P=np.zeros_like(S); occ=nelec//2; E=0; C=None
    for _ in range(300):
        J=np.einsum('ijkl,kl->ij',eri,P); K=np.einsum('ikjl,kl->ij',eri,P)
        F=Hc+J-0.5*K; Fp=X.T@F@X; e,Cp=np.linalg.eigh(Fp); C=X@Cp
        Pn=2*C[:,:occ]@C[:,:occ].T; En=0.5*np.sum((Hc+F)*Pn)+enuc
        if abs(En-E)<1e-10: E=En; break
        P,E=Pn,En
    return E,C,Hc

if __name__=="__main__":
    # Verify LiH against pyscf
    geom=[('Li',(0,0,0)),('H',(0,0,1.6))]
    S,T,V,eri,enuc=integrals(geom)
    E,C,Hc=rhf(S,T,V,eri,enuc,nelec=4)
    print(f"LiH pure-python RHF = {E:.6f} Ha  (nbasis={len(S)})")
    from pyscf import gto,scf
    m=gto.M(atom='Li 0 0 0; H 0 0 1.6',basis='sto-3g',verbose=0)
    print(f"LiH pyscf       RHF = {scf.RHF(m).kernel():.6f} Ha")
