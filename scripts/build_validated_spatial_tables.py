"""Build reusable, outward radial coefficient tables for the T=2 residual proof."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
import flint
from purposeful.validated_spatial import SpatialOracle, RESIDUAL_PROFILES


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff",type=int,default=12)
    parser.add_argument("--resolved",type=int,default=24)
    parser.add_argument("--bits",type=int,default=100)
    parser.add_argument("--target-bits",type=int,default=70)
    args=parser.parse_args()
    oracle=SpatialOracle(bits=args.bits,target_bits=args.target_bits)
    arrays={}; entries=[]; started=time.perf_counter()
    single_bound=3*(args.resolved+args.cutoff)**2
    product_bound=12*args.cutoff**2
    def add(key,p,bound,derivative=0):
        center,radius=oracle.table(p,bound,derivative)
        if not np.all(np.isfinite(center)) or not np.all(np.isfinite(radius)):
            raise ArithmeticError("nonfinite coefficient export: "+key)
        if np.max(radius)>1e-10:
            raise ArithmeticError("valid coefficient enclosure too wide for residual use: "+key)
        arrays[key+"__center"]=center;arrays[key+"__radius"]=radius
        entries.append(dict(key=key,derivative=derivative,max_n2=bound,
            maximum_radius=float(max(radius)),minimum_radius=float(min(radius)),
            inner=[(t.power,str(t.exponent),str(t.coefficient)) for t in p.inner],
            outer=[(t.power,str(t.exponent),str(t.coefficient)) for t in p.outer],
            exterior=[(t.power,str(t.exponent),str(t.coefficient)) for t in p.exterior]))
        print(key,bound,'maxradius',float(max(radius)),
              'seconds',round(time.perf_counter()-started,2),flush=True)
    for name,p in RESIDUAL_PROFILES.items():
        add(name,p,single_bound)
        if name=="drift":add(name+"__d1",p,single_bound,1)
    names=list(RESIDUAL_PROFILES)
    for i,name in enumerate(names):
        for other in names[i:]:
            p=RESIDUAL_PROFILES[name]*RESIDUAL_PROFILES[other]
            key=name+"__times__"+other
            add(key,p,product_bound)
            if "drift" in (name,other):add(key+"__d1",p,product_bound,1)
            if name==other=="drift":add(key+"__d2",p,product_bound,2)
    # Flux squared in the q1 block has the inverse Gibbs spatial weight.
    from purposeful.validated_spatial import INVERSE_WEIGHT
    add("inverse_weight",INVERSE_WEIGHT,12*args.resolved**2)
    out=ROOT/"results/continuum-validated-spatial-tables.npz"
    np.savez_compressed(out,**arrays)
    sources=["src/purposeful/validated_spatial.py","scripts/build_validated_spatial_tables.py",
             "tests/test_validated_spatial.py"]
    receipt=dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        method="Arb validated complex analytic quadrature; exact rational declarations; outward binary64 center-radius export",
        scope="Spatial coefficient enclosures only. A survival certificate additionally needs validated temporal residual contractions and tail bounds.",
        flint_version=flint.__version__,precision_bits=args.bits,target_bits=args.target_bits,
        cutoff=args.cutoff,resolved_cutoff=args.resolved,entries=entries,
        integration_calls=oracle.integral_calls,seconds=time.perf_counter()-started,
        array_sha256=sha256(out.read_bytes()).hexdigest(),
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})
    (ROOT/"results/continuum-validated-spatial-tables.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
