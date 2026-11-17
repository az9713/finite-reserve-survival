"""Replay the full rational cover and cross-check every point by expm."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import hashlib
import json
import sys

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from purposeful.assembly import Assembly, success_interval, curvature_bounds
from purposeful.autonomous import spectral_interval


def main():
    path = ROOT / 'results/assembly.json'
    record = json.loads(path.read_text(encoding='utf-8'))
    model = Assembly(**{key: int(value) if key in ['pairs','fuel','fuel_factor','damage_factor'] else F(value)
                        for key,value in record['model'].items()})
    deadline,lifetime = record['deadline'],record['lifetime']
    cert = record['certificate']
    for name,digest in record['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    points = {}
    widest = 0.
    for i,entry in enumerate(cert['points']):
        alpha,z=F(entry['alpha']),F(entry['z'])
        lo,hi=F(entry['lower']),F(entry['upper'])
        assert (lo,hi)==success_interval(model,deadline,lifetime,alpha,z)
        independent=model.formation_and_maintenance(deadline,lifetime,alpha,z)
        assert float(lo)-1e-13 <= independent <= float(hi)+1e-13
        points[alpha,z]=(lo,hi)
        widest=max(widest,float(hi-lo))
        if (i+1)%5000==0:
            print(f'Replayed {i+1} rational point enclosures and expm checks',flush=True)
    nodes=cert['nodes']
    seen=set()
    leaves=[]
    def visit(identifier,expected):
        assert identifier not in seen
        seen.add(identifier)
        node=nodes[identifier]
        box=tuple(map(F,node['box']))
        assert box==expected
        a0,a1,z0,z1=box
        assert 0<=a0<a1<=1 and 0<=z0<z1<=1
        ca,cz=curvature_bounds(model,deadline,lifetime,box)
        upper=min(F(1),max(points[a,z][1] for a in [a0,a1] for z in [z0,z1])
                  +ca*(a1-a0)**2/8+cz*(z1-z0)**2/8)
        assert F(node['upper'])==upper
        children=node['children']
        if not children:
            leaves.append(node)
            return
        assert len(children)==2
        first=tuple(map(F,nodes[children[0]]['box']))
        if first[1]<a1:
            mid=(a0+a1)/2
            expected_children=[(a0,mid,z0,z1),(mid,a1,z0,z1)]
        else:
            mid=(z0+z1)/2
            expected_children=[(a0,a1,z0,mid),(a0,a1,mid,z1)]
        for child,childbox in zip(children,expected_children):
            visit(child,childbox)
    visit(0,(F(0),F(1),F(0),F(1)))
    assert len(seen)==len(nodes)
    lower=max(pair[0] for pair in points.values())
    upper=max(F(node['upper']) for node in leaves)
    assert lower==F(cert['lower']) and upper==F(cert['upper'])
    assert upper-lower<=F(cert['tolerance'])
    construction=tuple(map(F,cert['construction']))
    assert points[construction][0]==lower
    boundary=F(0)
    for node in leaves:
        a0,a1,z0,z1=map(F,node['box'])
        for axis,value in [('a',F(0)),('a',F(1)),('z',F(0)),('z',F(1))]:
            if axis=='a' and value in (a0,a1):
                _,cz=curvature_bounds(model,deadline,lifetime,(value,value,z0,z1))
                value_bound=max(points[value,z][1] for z in [z0,z1])+cz*(z1-z0)**2/8
            elif axis=='z' and value in (z0,z1):
                ca,_=curvature_bounds(model,deadline,lifetime,(a0,a1,value,value))
                value_bound=max(points[a,value][1] for a in [a0,a1])+ca*(a1-a0)**2/8
            else:
                continue
            boundary=max(boundary,min(F(node['upper']),value_bound))
    assert boundary==F(cert['boundary_upper'])
    assert boundary<lower
    states,matrix,den=model.killed_generator(*construction)
    initial=np.array([states.index((1,0,model.fuel-1))-1])
    spectral=spectral_interval(matrix[1:,1:],den,initial)
    live=states[1:]
    distribution=expm(lifetime*matrix[1:,1:]/den)[initial[0]]
    distribution/=distribution.sum()
    allocation={'assembled':sum(p*(h+d) for p,(h,d,b) in zip(distribution,live)),
                'remaining_fuel':sum(p*b for p,(h,d,b) in zip(distribution,live)),
                'net_repairs':sum(p*(model.fuel-b-h-d) for p,(h,d,b) in zip(distribution,live))}
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),
            'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'source_hashes_match':True,'replayed_points':len(points),'replayed_nodes':len(nodes),
            'leaf_rectangles':len(leaves),'full_square_covered':True,'all_point_enclosures_recomputed':True,
            'all_points_cross_checked_with_expm':True,'max_point_enclosure_width':widest,
            'global_lower':str(lower),'global_upper':str(upper),'global_gap':str(upper-lower),
            'boundary_upper':str(boundary),'strictly_interior_optimum':True,
            'attaining_lower_bound_design':list(map(str,construction)),
            'reachable_spectral_enclosure':list(map(str,spectral)),
            'conditional_allocation_at_certified_design':allocation,
            'verification_scope':'Replay of every arithmetic enclosure and tree partition. This is an additional local implementation, with mathematical proof obligations stated in the study.'}
    (ROOT/'results/assembly-certificate-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
