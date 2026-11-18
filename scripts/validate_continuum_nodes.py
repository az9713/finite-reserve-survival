"""Build exact dyadic continuum trial archives and verified interpolation node plans."""
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


def build_trial_operator(cutoff, nodes, policy):
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

    return locals()


def build_curve_archive(cutoff=12, nodes=32, policy='spectator', degree=12,
                        bits=40, output=None):
    """Generate, then quantize trials. This generation is NOT validation."""
    import hashlib
    from fractions import Fraction
    from purposeful.validated_time import DyadicTimeCurve
    policy_label=str(policy)
    if policy != 'spectator':
        policy=Fraction(policy)
    op = build_trial_operator(cutoff, nodes, policy)
    nstate, matrix_dim, rho = op['nstate'], op['matrix_dim'], op['rho']
    compressed, live_states = op['compressed'], op['live_states']
    knots = [Fraction(0)] + [Fraction(2)**j for j in range(-9, 2)]
    local_nodes = (1-np.cos(np.pi*np.arange(degree+1)/degree))/2
    times = np.concatenate([float(a)+float(b-a)*local_nodes
                            for a,b in zip(knots, knots[1:])])
    expansion_degree = _chebyshev_degree(rho*times.max(), 1.e-16)
    weights = ive(np.arange(expansion_degree+1)[:,None], rho*times[None,:])
    weights[1:] *= 2
    launch = op['state_id'][1, op['stock'], 0]
    reward = np.zeros((nstate, matrix_dim)); reward[:,0]=1
    density = np.zeros_like(reward); density[launch,0]=1/op['chemistry_mass'][launch]
    def evolve(init):
        previous = init
        dest = weights[0,:,None,None]*previous
        current = previous-compressed(previous)/rho
        for k in range(1, expansion_degree+1):
            dest += weights[k,:,None,None]*current
            if k != expansion_degree:
                previous,current = current,2*(current-compressed(current)/rho)-previous
        return dest
    # Raw orbit sums of cosine products have exact rational Fourier expansion.
    positive = np.count_nonzero(op['modes'],axis=1)
    basis_conversion = 2.**(positive/2) / op['norm']
    vandermonde = np.polynomial.polynomial.polyvander(local_nodes, degree)
    trial_arrays=[]
    for initial in (reward,density):
        samples = evolve(initial) * basis_conversion
        segments = samples.reshape(len(knots)-1,degree+1,nstate,matrix_dim)
        coefficients = np.linalg.solve(vandermonde,
            np.moveaxis(segments,1,0).reshape(degree+1,-1))
        coefficients = np.moveaxis(coefficients.reshape(degree+1,len(knots)-1,nstate,matrix_dim),0,1)
        initial_exact = np.zeros((nstate,matrix_dim),dtype=object)
        if initial is reward:
            initial_exact[:,0]=1
        else:
            initial_exact[launch,0]=Fraction(59049,4)
        trial = DyadicTimeCurve.quantize(knots,coefficients,bits,initial_exact)
        trial_arrays.append(trial)
    if output is None:
        output=ROOT/'results'/'continuum-trials'/f'K{cutoff}-{policy_label.replace("/","_")}.npz'
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    arrays = [np.asarray(curve.integers,dtype=np.int64) for curve in trial_arrays]
    if any(not np.array_equal(raw, curve.integers) for raw,curve in zip(arrays,trial_arrays)):
        raise OverflowError('trial quantization exceeds signed64 integer range')
    np.savez_compressed(output, primal=arrays[0],dual=arrays[1],bits=np.array(bits),
        knot_numerators=np.array([t.numerator for t in knots],dtype=np.int64),
        knot_denominators=np.array([t.denominator for t in knots],dtype=np.int64),
        modes=np.asarray(op['modes'],dtype=np.int64),states=np.asarray(live_states,dtype=np.int64))
    receipt={'schema':'dyadic-continuum-trial-v1','policy':policy_label,'cutoff':cutoff,
        'quadrature_nodes_per_folded_axis':nodes,'degree':degree,'bits':bits,
        'shape':list(arrays[0].shape),'basis':'sum of distinct permutations of raw products cos(2*pi*k*x)',
        'archive':str(output.relative_to(ROOT)).replace('\\','/'),
        'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
        'initial_primal_exact':True,'initial_dual_exact':True,'continuity_exact':True,
        'status':'exact trial archive; continuum residual NOT yet validated',
        'seconds':perf_counter()-op['started']}
    output.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def load_curve_archive(path):
    from fractions import Fraction
    from purposeful.validated_time import DyadicTimeCurve
    data=np.load(path,allow_pickle=False)
    knots=tuple(Fraction(int(n),int(d)) for n,d in zip(data['knot_numerators'],data['knot_denominators']))
    return (DyadicTimeCurve(knots,data['primal'],int(data['bits'])),
            DyadicTimeCurve(knots,data['dual'],int(data['bits'])),
            data['modes'],data['states'])


def verify_curve_archive(path):
    """Replay exact initial, continuity, mode/state, hash and shape assertions."""
    import hashlib
    from fractions import Fraction
    from purposeful.validated_time import fraction_record
    path=Path(path)
    manifest=json.loads(path.with_suffix('.json').read_text())
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != manifest['sha256']:
        raise ValueError('trial archive hash mismatch')
    primal,dual,modes,states=load_curve_archive(path)
    if primal.knots != dual.knots or primal.knots[0] != 0 or primal.knots[-1] != 2:
        raise ValueError('trial time horizons inconsistent')
    expected_primal=np.zeros(primal.integers.shape[2:],dtype=object)
    expected_primal[:,0]=1
    expected_dual=np.zeros_like(expected_primal)
    launch=np.where(np.all(states==[1,4,0],axis=1))[0]
    if len(launch)!=1:
        raise ValueError('missing unique launch state')
    expected_dual[int(launch[0]),0]=Fraction(59049,4)
    if not np.array_equal(primal.endpoint(),expected_primal) or not np.array_equal(dual.endpoint(),expected_dual):
        raise ValueError('incorrect exact initial data')
    cutoff=int(modes.max())
    expected_modes=np.asarray(list(combinations_with_replacement(range(cutoff+1),3)))
    if not np.array_equal(modes,expected_modes):
        raise ValueError('spatial mode convention mismatch')
    # Constructors above check all shared integer endpoint coefficients exactly.
    return {'schema':'dyadic-trial-replay-v1','archive':manifest['archive'],'sha256':digest,
        'checks':{'source_archive_hash':True,'exact_primal_initial_data':True,
                  'exact_dual_initial_data':True,'exact_time_continuity':True,
                  'raw_symmetric_cosine_modes':True,'horizon_two':True},
        'shape':list(primal.integers.shape),'integer_coefficients_checked':int(primal.integers.size+dual.integers.size),
        'knots':[fraction_record(t) for t in primal.knots],
        'status':'exact input replay passed; no continuum enclosure inferred'}


def build_root_node(index,cutoff=12,nodes=32,degree=12,bits=40):
    """Archive one dyadic model at an enclosed exact interpolation root."""
    from fractions import Fraction
    from purposeful.validated_time import chebyshev_node_plan
    plan=chebyshev_node_plan()
    if not 0 <= index < plan['count']:
        raise ValueError('Chebyshev root index out of range')
    node=plan['nodes'][index]
    policy=Fraction(int(node['alpha_numerator']),2**node['alpha_bits'])
    output=ROOT/'results'/'continuum-trials'/f'K{cutoff}-node{index:02d}.npz'
    receipt=build_curve_archive(cutoff,nodes,policy,degree,bits,output)
    receipt['chebyshev_root']=node
    receipt['interpolation_node_count']=plan['count']
    output.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cutoff',type=int,default=12)
    parser.add_argument('--nodes',type=int,default=32)
    parser.add_argument('--policy',default='spectator')
    parser.add_argument('--degree',type=int,default=12)
    parser.add_argument('--bits',type=int,default=40)
    parser.add_argument('--output')
    parser.add_argument('--node-plan',action='store_true')
    parser.add_argument('--verify-archive')
    parser.add_argument('--node-index',type=int)
    args=parser.parse_args()
    if args.node_index is not None:
        result=build_root_node(args.node_index,args.cutoff,args.nodes,args.degree,args.bits)
    elif args.verify_archive:
        result=verify_curve_archive(args.verify_archive)
    elif args.node_plan:
        from purposeful.validated_time import chebyshev_node_plan
        result=chebyshev_node_plan()
    else:
        result=build_curve_archive(args.cutoff,args.nodes,args.policy,args.degree,args.bits,args.output)
    print(json.dumps(result,indent=2))
