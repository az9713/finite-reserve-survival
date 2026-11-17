"""Numerical weighted Fourier Galerkin and sampled energy-residual diagnostics.

This is an exploratory calculator. Trapezoidal quadrature of residual norms,
truncation of Fourier fluxes, and sampled check times have no enclosures.
"""
import sys,time
import numpy as np
from scipy.fft import dctn,idctn,idstn
from scipy.integrate import quad_vec
from scipy.special import ive
from scipy.linalg import eigvalsh, cholesky, solve_triangular
from itertools import product,permutations,combinations_with_replacement
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import argparse
import json
from pathlib import Path
from run_continuum_fourier import ROOT, Organization, chemical_generators, equilibrium, _chebyshev_degree

M=Organization()
D,a,offset,E,N=M.diffusivity,M.radius,M.offset,M.depth,M.fuel
states,damage,repair=chemical_generators(M)
sidx=[i for i,(h,b) in enumerate(states) if h]
damage=damage[np.ix_(sidx,sidx)]
repair=repair[np.ix_(sidx,sidx)]
states=[states[i] for i in sidx]
weights=equilibrium(M)[sidx]
launch=states.index((1,N))
weights=np.sqrt(weights/weights[launch])
damage=weights[:,None]*damage/weights[None,:]
repair=weights[:,None]*repair/weights[None,:]
accept=a*np.sqrt(1-np.sqrt(offset/E))

def assemble(K):
    axes=(0,1,2)
    size=2*K+1
    m2=np.arange(3*size*size+1)
    freq=2*np.pi*np.sqrt(m2)
    def fn(r):
        c=1-(r/a)**2
        U=offset-E*c*c
        ww=np.exp(-U)
        g=c**4
        return 4*np.pi*r*r*np.array([ww-np.exp(-offset),g,g*ww,g*np.exp(-max(0.,U))])[:,None]*np.sinc(freq*r/np.pi)
    coeff,err=quad_vec(fn,0,a,points=[accept],epsrel=1e-12,epsabs=1e-13,limit=300)
    coeff[0,0]+=np.exp(-offset)
    reps=np.array(list(combinations_with_replacement(range(K+1),3)),dtype=int)
    orbits=[]
    mult=[]
    ind=np.empty((K+1,K+1,K+1),dtype=int)
    for i,k in enumerate(reps):
        orbit=sorted(set(tuple(np.array(p)*s) for p in permutations(k) for s in product((-1,1),repeat=3)))
        orbits.append(np.array(orbit))
        mult.append(len(orbit))
        for p in set(permutations(k)):ind[p]=i
    mult=np.sqrt(np.array(mult))
    C=np.empty((5,len(reps),len(reps)))
    for j,orb in enumerate(orbits):
        diff2=np.sum((reps[:,None,:]-orb[None,:,:])**2,axis=-1)
        v=coeff[:,diff2]
        C[:4,:,j]=np.sum(v,axis=-1)*mult/mult[j]
        C[4,:,j]=np.sum(v[0]*(reps@orb.T),axis=-1)*mult/mult[j]*4*np.pi**2*D
    C=(C+C.swapaxes(1,2))/2
    return reps,mult,ind,C,err

def weighted(K,policy,times,residual_checks=True):
    started=time.perf_counter()
    reps,mult,ind,(mass,G0,G1,A,stiff),err=assemble(K)
    print('assembled',K,len(reps),err,time.perf_counter()-started,flush=True)
    chol=cholesky(mass,lower=True)
    V=solve_triangular(chol,np.eye(len(reps)),lower=True)
    # transformed coeff y by z=L^T y. M-inner products are Euclidean.
    # block bare diffusion, weighted contact, and paid/passive switch
    W1=V@stiff@V.T
    R1=V@G1@V.T
    A01=A@V.T
    A1=V@A@V.T
    Lambda=4*np.pi**2*D*np.sum(reps**2,axis=1)
    ix={s:i for i,s in enumerate(states)}
    C=[]
    for i,(h,b) in enumerate(states):
      if b:
        j=ix[(h,b-1)]
        s=(b-1)/(N-1) if policy=='auto' else policy
        f=M.paid*s*b;r=M.paid*s*(N-b+1)/8
        C.append((i,j,f,r,np.sqrt(f*r)))
    def H(u):
        out=np.empty_like(u)
        out[0]=damage@u[0]+(repair@u[0])@G0-u[0]*Lambda
        out[1]=damage@u[1]+(repair@u[1])@R1-u[1]@W1
        out[0]+=M.passive*(u[1]@A01.T-u[0]@A)
        out[1]+=M.passive*(u[0]@A01-u[1]@A1)
        for i,j,f,r,pair in C:
            out[0,i]+=pair*(A01@u[1,j])-f*(A@u[0,i])
            out[1,j]+=pair*(A01.T@u[0,i])-r*(A1@u[1,j])
        return out
    q=(max(np.max(Lambda),eigvalsh(W1,subset_by_index=[len(reps)-1,len(reps)-1])[0])+170)/2
    degree=_chebyshev_degree(q*max(times),1e-15)
    coeff=ive(np.arange(degree+1)[:,None],q*np.array(times)[None,:]);coeff[1:]*=2
    launchv=np.zeros((2,len(states),len(reps)));launchv[0,launch,0]=1
    output=np.zeros_like(launchv)
    output[0,:,0]=weights
    output[1]=weights[:,None]*chol.T[0:len(reps),0][None,:]
    def propagate(v):
      prev=v
      curr=v+H(v)/q
      res=coeff[0,:,None,None,None]*prev+coeff[1,:,None,None,None]*curr
      for k in range(2,degree+1):
        prev,curr=curr,2*(curr+H(curr)/q)-prev
        res+=coeff[k,:,None,None,None]*curr
      return res
    p=propagate(launchv)
    d=propagate(output)
    print('propagated',K,policy,degree,np.sum(p[:]*output[None],axis=(1,2,3)),time.perf_counter()-started,flush=True)
    values = np.sum(p*output[None],axis=(1,2,3))
    if not residual_checks:return values.tolist(),None,None,None,None
    # evaluate strong residual r(t)=L v(t)-v'(t) in Gibbs norm
    nx=max(64,3*K+4)
    x=np.arange(nx+1)/(2*nx)
    r2=x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2
    core=np.maximum(0,1-r2/a**2)
    U=offset-E*core**2
    ww=np.exp(-U)
    g=core**4
    grow=g*np.exp(-np.maximum(U,0));gshrink=g*np.exp(np.minimum(U,0))
    nweight=np.ones(nx+1)/nx;nweight[[0,-1]]*=0.5
    quad=nweight[:,None,None]*nweight[None,:,None]*nweight[None,None,:]
    def coeff_grid(z,lap=False):
      y=z.copy();y[1]=z[1]@V
      if lap:y=-y*Lambda
      fill=np.zeros((2,len(states),nx+1,nx+1,nx+1))
      fill[:,:,:K+1,:K+1,:K+1]=(y/mult)[...,ind]
      return fill
    def idir(c):return idctn(c,type=1,norm='forward',axes=(-3,-2,-1),workers=2)
    abs_potential_grad=[4*E*core/a**2 * xx for xx in (x[:,None,None],x[None,:,None],x[None,None,:])]
    def residual(z):
      cg=coeff_grid(z)
      u=idir(cg)
      r=-idir(coeff_grid(H(z)))
      r+=idir(coeff_grid(z,lap=True))
      r+=np.einsum('ij,qjxyz->qixyz',damage,u,optimize=True)
      r+=np.einsum('ij,qjxyz->qixyz',repair,u,optimize=True)*g
      r[0]+=M.passive*grow*(u[1]-u[0])
      r[1]+=M.passive*gshrink*(u[0]-u[1])
      for i,j,f,b,pair in C:
        r[0,i]+=grow*(pair*u[1,j]-f*u[0,i])
        r[1,j]+=gshrink*(pair*u[0,i]-b*u[1,j])
      for axis in range(3):
        cg1=np.moveaxis(cg[1],axis+1,-1)
        sine=np.zeros_like(cg1[...,1:-1])
        sine[...,:K]=cg1[...,1:K+1]*(-2*np.pi*np.arange(1,K+1))
        physical=idctn(sine,type=1,norm='forward',axes=(-3,-2),workers=2)
        physical=np.moveaxis(idstn(physical,type=1,norm='forward',axes=(-1,),workers=2),-1,axis+1)
        ss=[slice(None)]*4;ss[axis+1]=slice(1,-1)
        r1=r[1]
        grad=abs_potential_grad[axis]
        r1[tuple(ss)]-=D*physical*grad[tuple(ss[1:])]
      norm=np.sqrt(np.sum((r[0]**2+r[1]**2*ww)*quad))
      rz = r.copy()
      rz[1] *= ww
      zhat = dctn(rz,type=1,norm='forward',axes=(-3,-2,-1),workers=2)
      low_mode=np.max(np.abs(zhat[...,:K+1,:K+1,:K+1]))
      zhat[...,:K+1,:K+1,:K+1] = 0.
      modes=np.arange(nx+1)
      k2=modes[:,None,None]**2+modes[None,:,None]**2+modes[None,None,:]**2
      k2[0,0,0]=1.
      div_hat=zhat / (2*np.pi*k2)
      flux_norm_squared=0.
      for axis in range(3):
        moved=np.moveaxis(div_hat,axis+2,-1)
        sine=moved[...,1:-1]*modes[1:-1]
        physical=idctn(sine,type=1,norm='forward',axes=(-3,-2),workers=2)
        physical=np.moveaxis(idstn(physical,type=1,norm='forward',axes=(-1,),workers=2),-1,axis+2)
        ss=[slice(None)]*3;ss[axis]=slice(1,-1)
        flux_norm_squared+=np.sum((physical[0]**2+physical[1]**2/ww[tuple(ss)])*quad[tuple(ss)])/D
      return norm,float(np.sqrt(flux_norm_squared)),float(low_mode)
    pp=np.array([residual(z) for z in p])
    dd=np.array([residual(z) for z in d])
    print('primal',K,pp,flush=True)
    print('dual',K,dd,flush=True)
    print('seconds',time.perf_counter()-started,flush=True)
    return values.tolist(),pp.tolist(),dd.tolist(),float(err),degree

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoffs",nargs="+",type=int,default=[8,12])
    args=parser.parse_args()
    times=np.linspace(0,2,5)
    cases=[]
    for cutoff in args.cutoffs:
        for policy in ["auto",1.]:
            values,primal,dual,estimate,degree=weighted(cutoff,policy,times)
            cases.append(dict(cutoff=cutoff,policy=policy,survival=values,
                              sampled_primal=primal,sampled_dual=dual,
                              radial_quadrature_error_estimate=estimate,chebyshev_degree=degree))
    sources=["scripts/run_continuum_weighted.py","scripts/run_continuum_fourier.py",
             "src/purposeful/organization.py","src/purposeful/spatial_repair.py"]
    receipt=dict(
        generated_utc=datetime.now(timezone.utc).isoformat(), model=asdict(M), times=times.tolist(),
        method="Weighted Fourier forms in the Gibbs space; cube orbit reduction",
        scope="Exploratory, unenclosed quadrature/roundoff/truncation/time and intensity sampling",
        residual_columns=["Gibbs strong L2", "weighted norm of Fourier-truncated flux divided by sqrt(D)",
                          "maximum absolute raw low coefficient in weighted residual"],
        cases=cases,source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},
    )
    (ROOT/"results/continuum-weighted.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
