"""Independent parameter-cover checks, including the sharper stability proof."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

_path = Path(__file__).resolve().parents[1] / 'scripts/check_validated_continuum.py'
_spec = importlib.util.spec_from_file_location('independent_continuum_checker', _path)
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)


class ParameterCoverTests(unittest.TestCase):
    def test_known_global_maxima_and_exact_constants(self):
        self.assertEqual(checker.self_check()['known_polynomial_maxima'], 5)
        constants = checker.analytic_constants()
        self.assertLess(F(constants['interpolation_remainder_upper']), F(18, 10**8))

    def test_interval_restriction_preserves_polynomial(self):
        p = [F((-1)**j, j+1) for j in range(38)]
        left, right = F(5, 32), F(7, 8)
        local = checker.restrict_power(p, left, right)
        for s in [F(0), F(1, 17), F(1, 2), F(1)]:
            self.assertEqual(checker.evaluate_power(local, s),
                             checker.evaluate_power(p, left+(right-left)*s))

    def test_synthetic_full_cover_and_incomplete_data_rejection(self):
        receipt = checker.synthetic_cover_check()
        self.assertTrue(receipt['full_parameter_interval_replayed'])
        self.assertTrue(receipt['incomplete_cover_rejected'])

    def test_actual_lebesgue_constant_bound(self):
        receipt = checker.certified_lebesgue_constant()
        bound = F(receipt['upper'])
        self.assertGreater(bound, F(3278, 1000))
        self.assertLess(bound, F(3279, 1000))
        self.assertEqual(len(receipt['intervals']), 39)
        self.assertTrue(receipt['full_unit_interval_covered'])
        # Each adjacent pair overlaps, retaining the algebraic boundary roots.
        intervals = [tuple(map(F, row['interval'])) for row in receipt['intervals']]
        self.assertEqual(intervals[0][1], 1)
        self.assertEqual(intervals[-1][0], 0)
        for higher, lower in zip(intervals, intervals[1:]):
            self.assertLessEqual(higher[0], lower[1])


if __name__ == '__main__':
    unittest.main()
