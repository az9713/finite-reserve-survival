"""Exact temporal bookkeeping and full-interval node interpolation checks."""
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import importlib.util
import unittest

import numpy as np

from purposeful.validated_time import (DyadicTimeCurve, corrected_observable,
    integrate_curves, integrate_pair, chebyshev_node_plan, node_interpolant,
    arb_endpoint_fraction, shifted_legendre_coefficients, residual_legendre_terms)


class ExactTimeTests(unittest.TestCase):
    def test_fast_refinement_matches_fraction_reference(self):
        from purposeful.validated_refinement import fast_refine
        random=np.random.default_rng(4927)
        for degree in (1,3,7,12):
            raw=random.normal(size=(4,degree+1,2,3))
            curve=DyadicTimeCurve.quantize([0,F(1,8),F(1,4),F(1,2),1],raw,bits=20)
            knots=(F(0),F(1,32),F(1,16),F(1,8),F(3,16),F(1,4),F(3,8),F(1,2),F(3,4),F(1))
            old=curve.refine(knots)
            new=fast_refine(curve,knots)
            np.testing.assert_array_equal(old.integers*2**new.bits,new.integers*2**old.bits)
        root=Path(__file__).resolve().parents[1]
        archive=root/'results'/'continuum-trials'/'K12-spectator.npz'
        if archive.exists():
            data=np.load(archive)
            knots=tuple(F(int(n),int(d)) for n,d in zip(data['knot_numerators'],data['knot_denominators']))
            curve=DyadicTimeCurve(knots,data['primal'][:,:,:2,:6],int(data['bits'])).reversed()
            common=tuple(sorted(set(knots+curve.knots)))
            old=curve.refine(common);new=fast_refine(curve,common)
            np.testing.assert_array_equal(old.integers*2**new.bits,new.integers*2**old.bits)

    def test_validated_batch_encloses_exact_legendre_data(self):
        root=Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('validated_runner_test',root/'scripts'/'run_validated_continuum.py')
        runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
        integers=np.array([[[2**60+1],[-2**59+3],[2**58-7]]],dtype=object)
        curve=DyadicTimeCurve((F(0),F(1,8)),integers,40)
        values,derivatives,weights=runner.legendre_batch(curve)
        from flint import arb,ctx
        exact=list(residual_legendre_terms(curve))
        with ctx.workprec(160):
            for j,(weight,value,derivative) in enumerate(exact):
                self.assertEqual(weights[j],weight)
                for balls,expected in [(values,value),(derivatives,derivative)]:
                    q=F(expected[0]); exact_value=arb(q.numerator)/q.denominator
                    self.assertTrue(exact_value >= arb(float(balls.lower()[j,0])))
                    self.assertTrue(exact_value <= arb(float(balls.upper()[j,0])))

    def test_quantized_curve_continuity_initial_and_archive(self):
        curve = DyadicTimeCurve.quantize([0, F(1, 2), 1],
            np.array([[[1.01], [-.21], [.03]], [[.6], [-.1], [.04]]]),
            bits=20, initial=np.array([F(1)]))
        self.assertEqual(curve.endpoint()[0], 1)
        self.assertEqual(curve.integers[0].sum(), curve.integers[1, 0, 0])
        with TemporaryDirectory() as folder:
            path = Path(folder) / "trial.json"
            curve.write(path)
            other = DyadicTimeCurve.read(path)
            np.testing.assert_array_equal(other.integers, curve.integers)
            self.assertEqual(other.knots, curve.knots)

    def test_reject_discontinuity(self):
        with self.assertRaisesRegex(ValueError, "discontinuous"):
            DyadicTimeCurve((F(0), F(1, 2), F(1)), np.array([[[1], [2]], [[2], [0]]]), 0)

    def test_derivative_and_time_reversal_exact(self):
        curve = DyadicTimeCurve((F(0), F(1, 2), F(1)),
            np.array([[[1], [2], [3]], [[6], [1], [4]]]), 0)
        np.testing.assert_array_equal(curve.derivative_coefficients(0), [[F(4)], [F(12)]])
        rev = curve.reversed()
        np.testing.assert_array_equal(rev.integers, [[[11], [-9], [4]], [[6], [-8], [3]]])
        np.testing.assert_array_equal(rev.reversed().integers, curve.integers)

    def test_exact_time_pairing_and_correction(self):
        pair = lambda x, y: sum(a*b for a, b in zip(x.flat, y.flat))
        curve = DyadicTimeCurve((F(0), F(1)), np.array([[[1], [-1]]]), 0)
        self.assertEqual(integrate_curves(curve, curve, pair), F(1, 3))
        self.assertEqual(corrected_observable(curve, curve, np.array([F(1)]), pair, pair), F(1, 3))
        from flint import arb
        actual = (-arb(1)).exp()
        # r=-t, A=1 gives both squared negative norms 1/3.
        self.assertTrue(abs(actual-arb(1)/3) < arb(1)/3)

    def test_legendre_energy_and_exact_common_refinement(self):
        pair = lambda x,y: sum(a*b for a,b in zip(x.flat,y.flat))
        coeff=np.array([[F(1)], [F(2)], [F(3)]],dtype=object)
        leg=shifted_legendre_coefficients(coeff)
        direct=integrate_pair(coeff,coeff,F(1),pair)
        orthogonal=sum(pair(c,c)/F(2*k+1) for k,c in enumerate(leg))
        self.assertEqual(orthogonal,direct)
        curve=DyadicTimeCurve((F(0),F(1,4),F(1)),
            np.array([[[1],[2],[3]], [[6],[-1],[2]]]),0)
        # Here one segment has duration3/4, so use a refinement preserving it.
        refined=curve.refine((F(0),F(1,8),F(1,4),F(1)))
        self.assertEqual(integrate_curves(curve,curve,pair),integrate_curves(refined,refined,pair))
        value=sum(weight*pair(uprime+u,uprime+u) for weight,u,uprime in residual_legendre_terms(curve))
        direct=0
        for j,(a,b) in enumerate(zip(curve.knots,curve.knots[1:])):
            c=curve.coefficients(j).copy()
            c[:-1]+=curve.derivative_coefficients(j)
            direct+=integrate_pair(c,c,b-a,pair)
        self.assertEqual(value,direct)

    def test_node_plan_covers_actual_roots(self):
        from flint import arb, ctx
        plan = chebyshev_node_plan(38)
        with ctx.workprec(256):
            for record in plan["nodes"]:
                j = record["index"]
                true = (1+(arb(2*j+1)/76).cos_pi())/2
                lo = record["root_lower"]; hi = record["root_upper"]
                lower = arb(int(lo["numerator"]))/int(lo["denominator"])
                upper = arb(int(hi["numerator"]))/int(hi["denominator"])
                self.assertTrue(true >= lower)
                self.assertTrue(true <= upper)
        rem = plan["derivative_interpolation_remainder"]
        self.assertLess(F(int(rem["numerator"]), int(rem["denominator"])), F(18, 100_000_000))

    def test_interpolant_linear_data_and_error_separation(self):
        # Exact node centers are dyadic evaluations of 1+2a. Root displacement
        # belongs to the separate uncertainty term, not DCT conversion error.
        plan = chebyshev_node_plan(8)
        records = []
        for node in plan["nodes"]:
            n = int(node["alpha_numerator"])
            records.append({**node, "center_numerator": str(2**60+2*n),
                            "center_bits": 60, "radius_upper": "0"})
        result = node_interpolant(records, count=8)
        values = [F(int(c["numerator"]), 2**c["bits"]) for c in result["power_coefficients"]]
        for x in [F(0), F(1, 3), F(1)]:
            value = sum(c*x**k for k,c in enumerate(values))
            self.assertLess(abs(value - (1+2*x)), F(1, 10**12))
        self.assertLess(F(result["coefficient_error_upper"]), F(1, 10**40))


if __name__ == "__main__":
    unittest.main()
