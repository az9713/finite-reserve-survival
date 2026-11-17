"""Retain exact analytic controls for the global constant-intensity optimizer."""
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from verify_organization import certificate, fraction_interval, exp_interval

ROOT=Path(__file__).resolve().parents[1]


def main():
    cases=[]
    for name,kind,alpha in [('interior','single',.5),('left_endpoint','decreasing',0.),
                            ('right_endpoint','increasing',1.)]:
        values=[]
        for m in range(81):
            value=(Fraction(int(m==1)) if kind=='single' else
                   Fraction(1,2**m) if kind=='decreasing' else 1-Fraction(1,2**m))
            values.append(fraction_interval(value))
        coefficients=tuple([v[d] for v in values] for d in (0,1))
        result=certificate(coefficients,Fraction(2),Decimal('1e-9'))
        el,eh=exp_interval(Decimal(-1),Decimal(-1))
        exact=(el,eh) if kind=='single' else (Decimal(1),Decimal(1)) if kind=='decreasing' else (1-eh,1-el)
        assert Decimal.from_float(result['lower'])<=exact[0]
        assert Decimal.from_float(result['upper'])>=exact[1]
        assert abs(result['attaining_sample_alpha']-alpha)<1e-6
        cases.append({'name':name,'known_optimum_alpha':alpha,'exact_optimum_interval':[str(v) for v in exact],
                      'certificate':result,'passed':True})
    names=['scripts/check_organization_certificate_controls.py','scripts/verify_organization.py','tests/test_organization.py']
    receipt={'status':'PASS','scope':'Analytic Poisson-series optimizer controls; no additional physical-model result.',
        'cases':cases,'source_sha256':{n:sha256((ROOT/n).read_bytes()).hexdigest() for n in names}}
    (ROOT/'results/organization-certificate-controls.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','cases':len(cases)}))


if __name__=='__main__':main()
