"""Certified n16 diagnostic removing the attractive well, not a retuned main case."""
from dataclasses import asdict, replace
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from purposeful.organization import Organization
from verify_organization import matrices, paid_coefficients, certificate, one_kernel_survival


def main():
    model=replace(Organization(),depth=0.)
    started=time.perf_counter();n=16;horizon=Fraction(2)
    kernels,rates,indices,probability=matrices(model,n)
    coefficients,tail=paid_coefficients(kernels['base'],kernels['paid'],
        [rates[name]*horizon for name in ['base','paid']],indices,probability)
    fixed=certificate(coefficients,rates['paid']*horizon,Decimal('2e-6'))
    auto,auto_tail=one_kernel_survival(kernels['spectator'],rates['spectator']*horizon,indices,probability)
    files=['src/purposeful/organization.py','scripts/verify_organization.py',
           'tests/test_organization.py','scripts/run_organization_ablation.py']
    record={'scope':'Exploratory removal of the attractive well: depth=0 changes mechanical force AND conformation free-energy differences. This is not a pure-force causal isolation.',
        'model':asdict(model),'n':n,'horizon':2,'fixed':fixed,'spectator_interval':auto,
        'gap_interval':[float(np.nextafter(auto[0]-fixed['upper'],-np.inf)),
                        float(np.nextafter(auto[1]-fixed['lower'],np.inf))],
        'base_poisson_tail_upper':tail,'spectator_poisson_tail_upper':auto_tail,
        'seconds':time.perf_counter()-started,
        'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in files}}
    (ROOT/'results/organization-force-ablation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
