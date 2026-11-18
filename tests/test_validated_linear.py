"""Adversarial small checks against independent Arb exact-float arithmetic."""
import unittest
from fractions import Fraction
import numpy as np
from flint import arb, ctx
from purposeful.validated_linear import BallArray as BA, rational_ball,positive_matmul
from purposeful.validated_residual import ResidualMatrices
from purposeful.validated_spatial import SpatialOracle, SpatialTerm, ONE, DRIFT


class LinearTests(unittest.TestCase):
    def test_upward_tiny_floor_preserves_bound(self):
        a=np.array([[1.e-307,0.,2.]])
        b=np.array([1.e-20,1.e-307,3.e-307])
        got=positive_matmul(a,b)
        with ctx.workprec(150):
            exact=sum(arb(float(x))*arb(float(y)) for x,y in zip(a[0],b))
            self.assertTrue(arb(float(got[0]))>=exact)

    def test_cancelling_dot_against_arb(self):
        x=np.array([1.e12,1.,-1.e12,1.e-12]); y=np.array([1.,1.,1.,1.e12])
        ans=BA.exact(x).matmul(BA.exact(y))
        with ctx.workprec(150):
            exact=sum(arb(float(a))*arb(float(b)) for a,b in zip(x,y))
            self.assertTrue(arb(float(ans.center),float(ans.radius)).contains(exact))

    def test_ball_product_corner(self):
        a=BA(np.array([[1.,-3.],[2.,4.]]),np.full((2,2),.01))
        b=BA(np.array([5.,-7.]),np.full(2,.02))
        got=a.matmul(b)
        for sa in (-1,1):
            for sb in (-1,1):
                with ctx.workprec(150):
                    for i in range(2):
                        val=sum((arb(float(a.center[i,j]))+sa*arb(float(a.radius[i,j])))
                                *(arb(float(b.center[j]))+sb*arb(float(b.radius[j]))) for j in range(2))
                        self.assertTrue(arb(float(got.center[i]),float(got.radius[i])).contains(val))

    def test_raw_orbit_gram_and_drift_sign(self):
        m=ResidualMatrices(1,2); o=SpatialOracle()
        r=tuple(m.modes[-1]); s=tuple(m.modes[1])
        def terms(k,drift):
            from purposeful.validated_residual import orbit
            from flint import acb
            vals=[]; nz=sum(v!=0 for v in k)
            with ctx.workprec(120):
                for freq in orbit(k):
                    kk=tuple(int(v) for v in freq)
                    if drift:
                        for j in range(3):
                            if kk[j]:
                                vals.append(SpatialTerm(DRIFT,kk,(j,),acb(0,2*arb.pi()*kk[j]/2**nz)))
                    else:
                        vals.append(SpatialTerm(ONE,kk,(),acb(1)/2**nz))
            return vals
        i=len(m.modes)-1; j=1
        for dl,dr in ((0,0),(0,1),(1,0),(1,1)):
            name_l='drift' if dl else 'one'; name_r='drift' if dr else 'one'
            got=m.gram((name_l,dl),(name_r,dr))[i,j]
            exact=o.inner(terms(r,dl),terms(s,dr)).real
            self.assertTrue(arb(float(got.center),float(got.radius)).contains(exact))

    def test_rational_enclosure(self):
        x=rational_ball(Fraction(1,3))
        with ctx.workprec(150):
            self.assertTrue(arb(float(x.center),float(x.radius)).contains(arb(1)/3))


if __name__=='__main__': unittest.main()
