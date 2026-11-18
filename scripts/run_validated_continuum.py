"""Validate exact dyadic trial curves using continuum spatial residual forms."""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from purposeful.validated_linear import BallArray as BA, rational_ball
from purposeful.validated_residual import ResidualMatrices, ContinuumValidator, power_to_legendre
from purposeful.validated_time import arb_endpoint_fraction
from purposeful.validated_refinement import fast_refine

spec=importlib.util.spec_from_file_location('trial_nodes',ROOT/'scripts'/'validate_continuum_nodes.py')
trial_nodes=importlib.util.module_from_spec(spec);spec.loader.exec_module(trial_nodes)


def integer_ball(integers,bits):
    """Exact int64/object-integer dyadics enclosed after binary64 conversion."""
    centers=np.asarray(integers,dtype=float)*2.**(-bits)
    # The int->float conversion is correctly rounded; scaling by power2 is exact
    # here. A deliberately doubled epsilon covers conversion and radius rounding.
    radius=np.nextafter(2*np.finfo(float).eps*np.abs(centers)+np.finfo(float).tiny,np.inf)
    return BA(centers,radius)


def stack_balls(rows):
    return BA(np.stack([r.center for r in rows]),np.stack([r.radius for r in rows]))


def legendre_batch(curve):
    """Return exact-weight positive temporal decomposition in certified balls."""
    values=[]; derivatives=[]; weights=[]
    for j,(a,b) in enumerate(zip(curve.knots,curve.knots[1:])):
        coeff=integer_ball(curve.integers[j],curve.bits)
        deriv_rows=[coeff[k+1]*rational_ball(Fraction(k+1,b-a))
                    for k in range(coeff.center.shape[0]-1)]
        deriv_rows.append(BA.exact(np.zeros(coeff.center.shape[1:])))
        derivative=stack_balls(deriv_rows)
        value_leg=power_to_legendre(coeff)
        deriv_leg=power_to_legendre(derivative)
        for k in range(len(coeff.center)):
            values.append(value_leg[k]);derivatives.append(deriv_leg[k])
            weights.append((b-a)/ (2*k+1))
    return stack_balls(values),stack_balls(derivatives),weights


def correction_batches(primal,dual):
    reverse=primal.reversed()
    knots=tuple(sorted(set(reverse.knots+dual.knots)))
    reverse,dual=fast_refine(reverse,knots),fast_refine(dual,knots)
    u,_,weights=legendre_batch(reverse)
    v,dv,weights2=legendre_batch(dual)
    if weights!=weights2:
        raise ArithmeticError('correction temporal weights disagree')
    return u,v,dv,weights


def _arb_float(x):
    return arb(float(x))


def integrate_energy(validator,curve,policy):
    values,derivatives,weights=legendre_batch(curve)
    # Vectorized validator accepts all temporal Legendre coefficients in a batch.
    bounds,diagnostics=validator.residual_energy(values,derivatives,policy)
    bounds=np.asarray(bounds)
    if bounds.shape != (len(weights),):
        raise ArithmeticError('residual validator returned wrong temporal batch shape')
    with ctx.workprec(160):
        total=sum(arb(w.numerator)/w.denominator*_arb_float(v) for w,v in zip(weights,bounds))
        upper=arb_endpoint_fraction(total)
    diagnostics={k:np.asarray(v).tolist() for k,v in diagnostics.items()}
    return upper,{'legendre_terms':len(weights),'energy_upper':str(upper),'diagnostics':diagnostics}


def integrate_correction(validator,primal,dual,policy):
    u,v,dv,weights=correction_batches(primal,dual)
    pairings=validator.pairing(u,validator.sources(v,dv,policy))
    if pairings.center.shape != (len(weights),):
        raise ArithmeticError('correction validator returned wrong temporal batch shape')
    answer=BA.exact(0.)
    for j,w in enumerate(weights):
        answer=answer+pairings[j]*rational_ball(w)
    return answer,{'legendre_terms':len(weights),'center':float(answer.center),'radius':float(answer.radius)}


def validate_archive(path,flux_cutoff=24,cache=None,output=None):
    started=perf_counter();path=Path(path)
    source_files=['scripts/run_validated_continuum.py','scripts/validate_continuum_nodes.py',
                  'src/purposeful/validated_time.py','src/purposeful/validated_linear.py',
                  'src/purposeful/validated_spatial.py','src/purposeful/validated_residual.py',
                  'src/purposeful/validated_refinement.py']
    source_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in source_files}
    replay=trial_nodes.verify_curve_archive(path)
    primal,dual,modes,states=trial_nodes.load_curve_archive(path)
    manifest=json.loads(path.with_suffix('.json').read_text())
    policy=manifest['policy']
    matrices=ResidualMatrices(cutoff=int(modes.max()),flux_cutoff=flux_cutoff,cache=cache)
    validator=ContinuumValidator(matrices)
    if not np.array_equal(states,np.asarray(validator.states)):
        raise ValueError('chemical state convention mismatch')
    eu,ureceipt=integrate_energy(validator,primal,policy)
    ev,vreceipt=integrate_energy(validator,dual,policy)
    correction,creceipt=integrate_correction(validator,primal,dual,policy)
    terminal=integer_ball(dual.integers[-1].sum(axis=0),dual.bits)
    observable=validator.observable(terminal)
    corrected=observable-correction
    with ctx.workprec(160):
        radius=((arb(eu.numerator)/eu.denominator)*(arb(ev.numerator)/ev.denominator)).sqrt()
        lower=_arb_float(float(corrected.lower()))-radius
        upper=_arb_float(float(corrected.upper()))+radius
        lower_q=arb_endpoint_fraction(lower,False);upper_q=arb_endpoint_fraction(upper)
        radius_q=arb_endpoint_fraction(radius)
    end_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in source_files}
    if end_hashes != source_hashes:
        raise ArithmeticError('critical source changed during validation; rerun from a frozen source snapshot')
    cache_hashes={}
    if matrices.cache:
        for key in matrices._matrices:
            target=matrices.cache/f'{key}-K{matrices.K}-M{matrices.M}.npz'
            digest=hashlib.sha256(target.read_bytes()).hexdigest()
            if target.with_suffix('.sha256').read_text().strip() != digest:
                raise ArithmeticError('consumed matrix cache hash changed during validation')
            cache_hashes[target.name]=digest
    result={'schema':'validated-continuum-node-v1','policy':policy,'cutoff':int(modes.max()),
        'flux_cutoff':flux_cutoff,'horizon':2,'trial':replay,
        'primal':ureceipt,'dual':vreceipt,'correction':creceipt,
        'corrected_observable_center':float(corrected.center),'corrected_observable_radius':float(corrected.radius),
        'residual_product_upper':str(radius_q),'survival_lower':str(lower_q),'survival_upper':str(upper_q),
        'source_sha256':source_hashes,'source_unchanged_during_validation':True,
        'consumed_matrix_cache_sha256':cache_hashes,
        'matrix_cache':str(matrices.cache.resolve().relative_to(ROOT)).replace('\\','/') if matrices.cache else None,
        'seconds':perf_counter()-started,
        'arithmetic_model':'IEEE754 binary64 round-to-nearest with gradual underflow; Arb radial integrals and rational constants',
        'status':'validated continuum node enclosure'}
    if output:
        Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive')
    parser.add_argument('--flux-cutoff',type=int,default=24)
    parser.add_argument('--cache',default=str(ROOT/'results'/'continuum-residual-matrices'))
    parser.add_argument('--output')
    args=parser.parse_args()
    print(json.dumps(validate_archive(args.archive,args.flux_cutoff,args.cache,args.output),indent=2))
