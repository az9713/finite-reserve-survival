"""Run all continuum nodes and assemble a full-interval comparison receipt.

This driver can be restarted. It verifies input/source hashes before reusing
node enclosures; the separate checker independently replays parameter coverage.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
from purposeful.validated_time import chebyshev_node_plan,node_interpolant
from run_validated_continuum import validate_archive
from check_validated_continuum import check_parameter_receipt,certified_lebesgue_constant


def digest(p): return sha256(Path(p).read_bytes()).hexdigest()


def read_valid_node(path, archive):
    data=json.loads(path.read_text())
    if data['status']!='validated continuum node enclosure':
        raise ValueError('wrong node status')
    for name,expected in data['source_sha256'].items():
        if digest(ROOT/name)!=expected:
            raise ValueError('node source changed: '+name)
    if digest(archive)!=data['trial']['sha256']:
        raise ValueError('node trial changed')
    # Caches are rebuildable intermediates, excluded from Git. A fresh clone
    # must recompute nodes rather than silently reuse an unreplayable chain.
    cache=ROOT/data['matrix_cache']
    for name,expected in data['consumed_matrix_cache_sha256'].items():
        matrix=cache/name
        if not matrix.exists() or digest(matrix)!=expected:
            raise ValueError('node matrix cache missing or changed: '+name)
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flux-cutoff',type=int,default=24)
    parser.add_argument('--cache',default='results/continuum-residual-matrices-v1')
    parser.add_argument('--wait-for-trials',action='store_true')
    args=parser.parse_args()
    plan=chebyshev_node_plan()
    output=ROOT/'results/continuum-validated-nodes';output.mkdir(exist_ok=True)
    requests=[('spectator',None)]+[(f'node{j:02d}',row) for j,row in enumerate(plan['nodes'])]
    receipts=[];records=[]
    for name,row in requests:
        archive=ROOT/'results/continuum-trials'/f'K12-{name}.npz'
        while not archive.exists() or not archive.with_suffix('.json').exists():
            if not args.wait_for_trials:
                raise FileNotFoundError(archive)
            print('Waiting for exact trial',name,flush=True)
            time.sleep(10)
        receipt=output/f'{name}.json'
        reuse=False
        if receipt.exists():
            try:
                data=read_valid_node(receipt,archive);reuse=data['flux_cutoff']==args.flux_cutoff
            except (ValueError,KeyError):
                pass
        if not reuse:
            print('Validating',name,flush=True)
            data=validate_archive(archive,args.flux_cutoff,ROOT/args.cache,receipt)
        receipts.append({'path':receipt.relative_to(ROOT).as_posix(),'sha256':digest(receipt),
                         'trial':archive.relative_to(ROOT).as_posix(),'trial_sha256':digest(archive)})
        lo,hi=F(data['survival_lower']),F(data['survival_upper'])
        print(name,'survival',[float(lo),float(hi)],'radius',float(F(data['residual_product_upper'])),flush=True)
        if row is None:
            if data['policy']!='spectator':
                raise ValueError('spectator receipt does not identify spectator policy')
            spectator={'lower':str(lo),'upper':str(hi)}
        else:
            alpha=F(int(row['alpha_numerator']),2**row['alpha_bits'])
            if F(data['policy'])!=alpha:
                raise ValueError('node trial policy differs from exact dyadic parameter')
            center_bits=60
            center=round(((lo+hi)/2)*2**center_bits)
            centerq=F(center,2**center_bits)
            records.append({**row,'center_numerator':str(center),'center_bits':center_bits,
                            'radius_upper':str(max(centerq-lo,hi-centerq)),
                            'node_receipt':receipt.relative_to(ROOT).as_posix()})
    interpolant=node_interpolant(records,coefficient_bits=160)
    candidate={'schema':'validated-continuum-comparison-v1','horizon':2,
        'generated_utc':datetime.now(timezone.utc).isoformat(),'node_records':records,
        'power_coefficients':interpolant['power_coefficients'],
        'coefficient_error_upper':interpolant['coefficient_error_upper'],
        'spectator':spectator,'node_receipts':receipts,
        'scope':'Brownian continuum, original parameters and uniform launch, spectator versus every constant paid intensity in [0,1]',
        'arithmetic_model':'Arb analytic integrals plus explicit binary64 roundoff enclosures; exact rational time and parameter coverage',
        'source_sha256':{'scripts/complete_validated_continuum.py':digest(__file__)}}
    lebesgue=certified_lebesgue_constant()
    cover=check_parameter_receipt(candidate,lebesgue_certificate=lebesgue)
    candidate['parameter_cover']=cover
    candidate['certified_lebesgue_constant']=lebesgue
    delta=F(cover['spectator_minus_best_fixed_lower'])
    if delta>0:
        candidate['advantage_time_window']={
            'open_left':str(F(2)-delta),'open_right':str(F(2)+delta),
            'gap_lower_formula':'delta - abs(t-2)',
            'proof':'Only killing is h=1 to h=0 at rate 1; every survival derivative lies in [-1,0], so every arm difference, and their infimum, is 1-Lipschitz.'}
    candidate['status']='STRICT_CONTINUUM_ADVANTAGE_CERTIFIED' if cover['strict_gap_arithmetically_verified'] else 'CONTINUUM_COMPARISON_ENCLOSED_SIGN_UNRESOLVED'
    target=ROOT/'results/continuum-validated-comparison.json'
    target.write_text(json.dumps(candidate,indent=2)+'\n')
    print(candidate['status'],flush=True)
    print('gap lower',float(F(cover['spectator_minus_best_fixed_lower'])),flush=True)


if __name__=='__main__': main()
