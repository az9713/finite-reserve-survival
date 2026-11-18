"""Exact dyadic time restriction using integer array operations only."""
from fractions import Fraction
from math import comb

import numpy as np

from .validated_time import DyadicTimeCurve


def fast_refine(curve,knots):
    """Same polynomial as curve.refine, avoiding per-entry Fraction arithmetic.

    Affine restriction s_old=a+b*s_new gives the exact rational factors
    binom(k,j)*a**(k-j)*b**j. A common power-two denominator permits every
    coefficient-array operation below to be an exact Python bigint operation.
    """
    knots=tuple(Fraction(t) for t in knots)
    if knots[0]!=curve.knots[0] or knots[-1]!=curve.knots[-1] or not set(curve.knots).issubset(knots):
        raise ValueError('new partition must refine the existing partition')
    if any(b<=a for a,b in zip(knots,knots[1:])):
        raise ValueError('new knots must strictly increase')
    plans=[];max_exponent=0;source=0;degree=curve.integers.shape[1]-1
    for left,right in zip(knots,knots[1:]):
        while curve.knots[source+1]<=left:
            source+=1
        duration=curve.knots[source+1]-curve.knots[source]
        offset=(left-curve.knots[source])/duration
        scale=(right-left)/duration
        rows=[]
        for j in range(degree+1):
            row=[]
            for k in range(j,degree+1):
                factor=comb(k,j)*offset**(k-j)*scale**j
                d=factor.denominator
                if d & (d-1):
                    raise ValueError('restriction factors must be dyadic')
                exponent=d.bit_length()-1
                max_exponent=max(max_exponent,exponent)
                row.append((k,factor.numerator,exponent))
            rows.append(row)
        plans.append((source,rows))
    blocks=[]
    for source,rows in plans:
        block=[]
        for row in rows:
            result=np.zeros(curve.integers.shape[2:],dtype=object)
            for k,numerator,exponent in row:
                if numerator:
                    factor=numerator << (max_exponent-exponent)
                    result+=curve.integers[source,k]*factor
            block.append(result)
        blocks.append(np.stack(block))
    return DyadicTimeCurve(knots,np.stack(blocks),curve.bits+max_exponent)
