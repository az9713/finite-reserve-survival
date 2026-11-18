"""Validate one disjoint shard of the 38 constant-intensity nodes."""
import argparse
import json
import os
from pathlib import Path
import time
from run_validated_continuum import validate_archive
from complete_validated_continuum import read_valid_node,ROOT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--step',type=int,default=1)
    p.add_argument('--cache',default='results/continuum-residual-matrices-final')
    args=p.parse_args()
    target=ROOT/'results/continuum-validated-nodes';target.mkdir(exist_ok=True)
    for j in range(args.start,38,args.step):
        name=f'node{j:02d}'
        archive=ROOT/'results/continuum-trials'/f'K12-{name}.npz'
        while not archive.exists() or not archive.with_suffix('.json').exists(): time.sleep(5)
        receipt=target/f'{name}.json'
        if receipt.exists():
            try:
                read_valid_node(receipt,archive)
                print(name,'valid existing receipt',flush=True)
                continue
            except (ValueError,KeyError): pass
        print(name,'begin validation',flush=True)
        data=validate_archive(archive,24,ROOT/args.cache,receipt)
        print(name,'completed',data['seconds'],'seconds; radius',float(__import__('fractions').Fraction(data['residual_product_upper'])),flush=True)


if __name__=='__main__': main()
