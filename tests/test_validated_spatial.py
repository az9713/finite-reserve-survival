"""Analytic checks of the outward radial/Fourier enclosure engine."""
import unittest
from fractions import Fraction

from flint import acb, arb, ctx

from purposeful.validated_spatial import (SpatialOracle, SpatialTerm, profile,
    ONE, CONTACT, WEIGHT, SWITCH, DRIFT, ball, midpoint_radius,
    longdouble_midpoint_radius, tail_flux_bound)


class ValidatedSpatialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.oracle=SpatialOracle(bits=120,target_bits=80)

    def test_exact_constant_modes_and_cube_moments(self):
        o=self.oracle
        with ctx.workprec(120):
            self.assertTrue(o.moment(ONE).contains(1))
            self.assertTrue(o.moment(ONE,(2,3,4)).contains(0))
            self.assertTrue(o.moment(ONE,(0,0,0),(0,0)).contains(arb(1)/12))
            self.assertTrue(o.moment(ONE,(1,0,0),(0,)).contains(acb(0,-1/(2*arb.pi()))))
            self.assertTrue(o.moment(ONE,(1,0,0),(0,0)).contains(-1/(2*arb.pi()**2)))

    def test_ball_indicator_closed_transform(self):
        o=self.oracle
        p=profile(compact=True)
        with ctx.workprec(120):
            a=ball(Fraction(8,25))
            self.assertTrue(o.moment(p).contains(4*arb.pi()*a**3/3))
            for k in ((1,0,0),(2,3,1),(12,12,12)):
                s=2*arb.pi()*arb(sum(n*n for n in k)).sqrt()
                exact=4*arb.pi()*((s*a).sin()-s*a*(s*a).cos())/s**3
                self.assertTrue(o.moment(p,k).contains(exact),str(k))

    def test_zero_mode_polynomial_integral(self):
        with ctx.workprec(120):
            a=ball(Fraction(8,25))
            # 4*pi*a^3*int_0^1 t²(1-t²)^4 dt.
            exact=512*arb.pi()*a**3/3465
            self.assertTrue(self.oracle.moment(CONTACT).contains(exact))

    def test_cartesian_moment_isotropic_at_zero(self):
        p=profile(compact=True)
        with ctx.workprec(120):
            a=ball(Fraction(8,25))
            exact=4*arb.pi()*a**5/15
            for i in range(3):
                self.assertTrue(self.oracle.moment(p,axes=(i,i)).contains(exact))
                self.assertTrue(self.oracle.moment(p,axes=(i,)).contains(0))
            self.assertTrue(self.oracle.moment(p,axes=(0,1)).contains(0))

    def test_cartesian_first_moment_closed_formula(self):
        p=profile(compact=True)
        with ctx.workprec(120):
            a=ball(Fraction(8,25)); pi=arb.pi()
            # derivative of 4pi[sin(sa)-sa cos(sa)]/s³ wrt s,
            # times i because x-transform = i d/ds at axial frequency.
            for n in (1,7):
                s=2*pi*n
                d=4*pi*(a*a*(s*a).sin()/s**2-
                    3*((s*a).sin()-s*a*(s*a).cos())/s**4)
                self.assertTrue(self.oracle.moment(p,(n,0,0),(0,)).contains(acb(0,d)))

    def test_cartesian_second_moments_nonzero_closed_formula(self):
        p=profile(compact=True)
        with ctx.workprec(120):
            a=ball(Fraction(8,25));pi=arb.pi()
            for k in ((1,0,0),(2,3,1),(7,0,0)):
                norm=arb(sum(n*n for n in k)).sqrt();s=2*pi*norm;t=s*a
                num=t.sin()-t*t.cos()
                hp=4*pi*(a*a*t.sin()/s**2-3*num/s**4)
                hpp=4*pi*(a**3*t.cos()/s**2-5*a*a*t.sin()/s**3+12*num/s**5)
                for i,j in ((0,0),(0,1),(1,1)):
                    direction=arb(k[i]*k[j])/(norm*norm)
                    exact=-(hpp*direction+hp/s*((i==j)-direction))
                    self.assertTrue(self.oracle.moment(p,k,(i,j)).contains(exact))

    def test_inner_product_includes_conjugation_and_weight(self):
        terms=[SpatialTerm(ONE,(1,2,3),coefficient=acb(1,2))]
        self.assertTrue(self.oracle.inner(terms,terms).contains(5))
        with ctx.workprec(120):
            self.assertTrue(self.oracle.inner(terms,terms,WEIGHT).overlaps(
                5*self.oracle.moment(WEIGHT)))

    def test_metropolis_branch_profile_product(self):
        # g*m(U), times exp(U), equals g inside U>0 and g*exp(U) inside U<0.
        from purposeful.validated_spatial import RadialProfile, Elementary
        expected=RadialProfile((Elementary(4,1),),(Elementary(4),))
        self.assertEqual(SWITCH*profile(exponent=1),expected)

    def test_serialization_float_conversions_enclose(self):
        with ctx.workprec(120):
            value=self.oracle.moment(WEIGHT).real
            for convert in (midpoint_radius,longdouble_midpoint_radius):
                m,r=convert(value)
                mn,md=m.as_integer_ratio();rn,rd=r.as_integer_ratio()
                exact_mid=arb(mn)/md; exact_rad=arb(rn)/rd
                self.assertGreaterEqual(exact_rad,(value-exact_mid).abs_upper())

    def test_tail_bound_rejects_inconsistent_parseval(self):
        with self.assertRaises(ArithmeticError):
            tail_flux_bound(arb(1),arb(2),arb(1),12)
        with ctx.workprec(120):
            value=tail_flux_bound(arb(3),arb(1),arb(2),4)
            self.assertTrue(value.contains(1/(10*arb.pi())**2))

    def test_float_parameter_rejected(self):
        with self.assertRaises(TypeError):
            SpatialOracle(radius=.32)

    def test_shell_table_includes_exterior_and_encloses(self):
        c,r=self.oracle.table(WEIGHT,3)
        with ctx.workprec(120):
            for n2 in range(4):
                value=self.oracle.transform(WEIGHT,n2)
                if n2==0:value+=self.oracle.exterior_value(WEIGHT)
                self.assertGreaterEqual(arb(float(r[n2])),
                    (value-arb(float(c[n2]))).abs_upper())
        with self.assertRaises(ValueError):
            self.oracle.table(WEIGHT,3,1)


if __name__=="__main__":
    unittest.main()
