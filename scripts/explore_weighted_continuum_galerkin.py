"""Screen a weighted Fourier residual approach using unvalidated sampled quadrature."""
import argparse
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.fft import dctn, idctn, idst
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.special import ive

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from purposeful.organization import Organization
from purposeful.spatial_repair import chemical_generators, equilibrium, _chebyshev_degree


def measure(cutoff, nodes, policy, time_nodes):
    started = perf_counter()
    model = Organization()
    stock = model.fuel
    states, damage, repair = chemical_generators(model)
    live = [i for i, (h, _) in enumerate(states) if h]
    live_states = [(h,b,q) for i in live for h,b in [states[i]] for q in (0,1)]
    state_id = {key: i for i, key in enumerate(live_states)}
    chemistry_mass = np.repeat(equilibrium(model)[live], 2)
    nstate = len(live_states)
    modes = list(combinations_with_replacement(range(cutoff + 1), 3))
    ids = {key:i for i,key in enumerate(modes)}
    ind = np.indices([cutoff + 1]*3)
    orbit = np.array([ids[tuple(sorted(row))] for row in ind.reshape(3,-1).T]).reshape(ind.shape[1:])
    norm = np.sqrt(np.bincount(orbit.ravel()))
    low = tuple(slice(0,cutoff+1) for _ in range(3))
    dims = (nodes,)*3
    axes = (-3,-2,-1)
    grid_axes = [((np.arange(nodes)+.5)/(2*nodes)).reshape(tuple(nodes if i==axis else 1 for i in range(3))) for axis in range(3)]
    r2 = sum(x**2 for x in grid_axes)
    core = np.maximum(0., 1-r2/model.radius**2)
    g = core**4
    U = model.offset-model.depth*core**2
    spatial_mass = np.exp(-U)
    mass_factor = np.stack([spatial_mass if q else np.ones(dims) for h,b,q in live_states])
    grad_U = [4*model.depth*core*x/model.radius**2 for x in grid_axes]
    on,off = g*np.exp(-np.maximum(U,0)),g*np.exp(np.minimum(U,0))
    constant = np.kron(-damage[np.ix_(live,live)],np.eye(2))
    contact = np.kron(-repair[np.ix_(live,live)],np.eye(2))
    edges = []
    for i,(h,b,q) in enumerate(live_states):
        if q:
            continue
        edges.append((i,state_id[h,b,1],model.passive*on,model.passive*off))
        if b:
            s=(b-1)/(stock-1) if policy=='spectator' else float(policy)
            edges.append((i,state_id[h,b-1,1],model.paid*b*s*on,model.paid*(stock-b+1)*s*off/8))
    lap = model.diffusivity*(2*np.pi)**2*np.sum(np.square(modes),axis=1)

    def coeff(values):
        transformed=dctn(values,axes=axes,norm='ortho',workers=4)[(...,*low)]
        if values.ndim == 3:
            return np.bincount(orbit.ravel(), weights=transformed.ravel())/norm
        return np.stack([np.bincount(orbit.ravel(), weights=state.ravel()) for state in transformed])/norm

    def modes_grid(values):
        ans=np.zeros((*values.shape[:-1],*dims))
        ans[(...,*low)]=(values/norm)[...,orbit]
        return ans

    def samples(values):
        return idctn(modes_grid(values),axes=axes,norm='ortho',workers=4,overwrite_x=True)

    def gradient(values):
        padded=modes_grid(values)
        for i in range(3):
            axis=padded.ndim-3+i
            raw=np.moveaxis(padded,axis,-1)
            temp=np.zeros(raw.shape)
            temp[...,:cutoff]=-raw[...,1:cutoff+1]*(2*np.pi*np.arange(1,cutoff+1))
            shifted=idst(np.moveaxis(temp,-1,axis),axis=axis, type=2, norm='ortho', workers=4,overwrite_x=True)
            parallel=tuple(j for j in axes if j!=i-3)
            yield idctn(shifted, axes=parallel,norm='ortho',workers=4,overwrite_x=True)

    def chemistry(y):
        out=(constant @ y.reshape(nstate,-1)).reshape(y.shape)
        out+=(contact @ y.reshape(nstate,-1)).reshape(y.shape)*g
        for i,j,kf,kr in edges:
            diff=y[i]-y[j]
            out[i]+=kf*diff
            out[j]-=kr*diff
        return out

    # Galerkin sample mass and weighted diffusion via form, no second derivative of U.
    matrix_dim=len(modes)
    mass=np.zeros((matrix_dim,matrix_dim))
    stiff=np.zeros_like(mass)
    for col in range(matrix_dim):
        probe=np.zeros(matrix_dim);probe[col]=1.
        mass[:,col]=coeff(spatial_mass*samples(probe))
        drift=sum(f*z for f,z in zip(grad_U,gradient(probe)))
        stiff[:,col]=mass[:,col]*lap[col]+model.diffusivity*coeff(spatial_mass*drift)
    mass=(mass+mass.T)/2
    stiff=(stiff+stiff.T)/2
    chol=cho_factor(mass,check_finite=False)
    diffusion=cho_solve(chol,stiff,check_finite=False)

    def compressed(vectors):
        projected=coeff(mass_factor*chemistry(samples(vectors)))
        projected[::2]+=vectors[::2]*lap
        projected[1::2]=cho_solve(chol,projected[1::2].T,check_finite=False).T+vectors[1::2]@diffusion.T
        return projected

    eig_max=eigh(stiff,mass,subset_by_index=[matrix_dim-1,matrix_dim-1],eigvals_only=True)[0]
    rho=(max(lap.max(),eig_max)+180)/2
    times, weights=leggauss(time_nodes)
    ts=np.concatenate(([0.,2.],times+1))
    degree=_chebyshev_degree(rho*ts.max(),1.e-16)
    weights_t=ive(np.arange(degree+1)[:,None],rho*ts[None,:]);weights_t[1:]*=2

    def evolve(init):
        prev=init
        dest=weights_t[0,:,None,None]*prev
        cur=prev-compressed(prev)/rho
        for k in range(1,degree+1):
            dest+=weights_t[k,:,None,None]*cur
            if k!=degree:
                prev,cur=cur,2*(cur-compressed(cur)/rho)-prev
        return dest

    launch=state_id[1,stock,0]
    reward=np.zeros((nstate,matrix_dim)); reward[:,0]=1
    d_init=np.zeros_like(reward); d_init[launch,0]=1/chemistry_mass[launch]
    primal=evolve(reward);dual=evolve(d_init)
    survival=float(primal[1,launch,0])
    # Sampled evaluation of the negative norm of the residual via Fourier antidivergence.
    kgrid=2*np.pi*np.arange(nodes)
    k2=kgrid[:,None,None]**2+kgrid[None,:,None]**2+kgrid[None,None,:]**2

    def residual(vectors):
        av=chemistry(samples(vectors))+samples(vectors*lap)
        for g,qgrad in zip(grad_U,gradient(vectors[1::2])):
            av[1::2]+=model.diffusivity*g*qgrad
        weighted=mass_factor*(av-samples(compressed(vectors)))
        full=dctn(weighted,axes=axes,norm='ortho',workers=4)
        mean=full[:,0,0,0]
        full[:,0,0,0]=0
        k2[0,0,0]=1
        phi=full/k2
        norm2=0.
        for i in range(3):
            ax=i+1
            raw=np.moveaxis(phi,ax,-1)
            temp=np.zeros(raw.shape);temp[...,:-1]=-raw[...,1:]*kgrid[1:]
            shifted=idst(np.moveaxis(temp,-1,ax),axis=ax,type=2,norm='ortho',workers=4)
            z=idctn(shifted,axes=tuple(j for j in axes if j!=i-3),norm='ortho',workers=4)
            norm2+=np.einsum('i,ijxy,ijxy->',chemistry_mass,z,z/mass_factor)/model.diffusivity
        # A tiny residual mean introduced by sampled quadrature can be enclosed via A>=lambda0.
        mean_norm=np.sqrt(np.dot(mean**2,chemistry_mass)/((7-np.sqrt(17))*np.exp(-2)/4))
        return (np.sqrt(norm2)+mean_norm)**2

    pr=np.sqrt(weights @ np.array([residual(v) for v in primal[2:]]))
    dr=np.sqrt(weights @ np.array([residual(v) for v in dual[2:]]))
    return {'cutoff':cutoff,'spatial_quadrature_points_per_axis':nodes,'time_quadrature_points':time_nodes,
            'policy':policy,'survival':survival,'primal_negative_norm':pr,'dual_negative_norm':dr,
            'sampled_product':pr*dr,'label':'exploratory quadrature, not an enclosure',
            'seconds':perf_counter()-started}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cutoff',type=int,default=10)
    p.add_argument('--nodes',type=int,default=24)
    p.add_argument('--times',type=int,default=12)
    p.add_argument('--policy',default='spectator')
    a=p.parse_args()
    print(json.dumps(measure(a.cutoff,a.nodes,a.policy,a.times),indent=2))
