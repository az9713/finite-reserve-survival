from fractions import Fraction
import unittest

import numpy as np
from scipy.linalg import eigvals, expm

from purposeful.autonomous import (
    Circuit, Network, benchmark_kinetics, reactions, spectral_interval,
    survival_curves, survival_intervals, uniform_activity_enclosure,
)


class AutonomousTests(unittest.TestCase):
    def setUp(self):
        self.kinetics = benchmark_kinetics()
        self.network = Network(self.kinetics)
        self.circuits, self.matrices = self.network.digital_circuits()
        self.initial = self.network.initial(2)

    def test_energy_ledger_and_reverses(self):
        expected = {"environment": 0, "sense": -1, "control": -1,
                    "harvest": 1, "wear": 0, "failure": 0, "repair": -1}
        for informed in [True, False]:
            for edge in reactions(self.kinetics, informed):
                self.assertGreater(edge.reverse, 0)
                self.assertEqual(edge.target[0] - edge.source[0], expected[edge.name])
                self.assertEqual(edge.forward / edge.reverse,
                                 edge.work_factor * self.kinetics.energy_factor(edge.source)
                                 / self.kinetics.energy_factor(edge.target))
                if edge.name == "failure":
                    self.assertEqual(edge.target[-1], 2)
        expected_killing = [
            self.kinetics.failure if x != a and d == 1 else 0
            for _, x, _, a, d in self.network.states
        ]
        np.testing.assert_array_equal(
            self.matrices[0].sum(axis=1), -self.network.denominator * np.array(expected_killing)
        )

    def test_rational_bounds_enclose_expm_and_isolate_optimum(self):
        lo, hi, scale = survival_intervals(self.matrices, self.network.denominator, 4, self.initial)
        winner = int(lo.argmax())
        self.assertEqual(self.circuits[winner], Circuit((0, 1), (0, 1, -1, 1)))
        self.assertGreater(lo[winner], np.delete(hi, winner).max())
        for i in [0, winner, len(self.circuits) - 1]:
            exact = expm(self.matrices[i] * (4 / self.network.denominator))[self.initial].sum(axis=1).mean()
            self.assertLessEqual(lo[i] / scale, exact)
            self.assertLessEqual(exact, hi[i] / scale)
        zero = survival_curves(self.matrices[[winner]], self.network.denominator, [0, 4], self.initial)
        self.assertAlmostEqual(zero[0, 0], 1.0)
        self.assertAlmostEqual(zero[0, 1], 0.533766303938, places=10)

    def test_continuum_enclosure_against_an_interior_setting(self):
        lo, hi = uniform_activity_enclosure(self.network, Fraction(1, 4), self.initial, subdivisions=16)
        _, base = self.network.control_bank()
        sensor = self.network.components[1,] + self.network.components[2,]
        inside = survival_curves(7 * base + 4 * sensor, 7 * self.network.denominator, [.25], self.initial)
        self.assertLessEqual(inside.max(), float(hi))
        self.assertGreater(lo, Fraction(9, 10))

    def test_asymptotic_positive_vector_enclosure(self):
        candidate = self.circuits.index(Circuit((0, 0), (1, -1, 1, -1)))
        lo, hi = spectral_interval(self.matrices[candidate], self.network.denominator, self.initial)
        visited = {int(i) for i in self.initial}
        pending = list(visited)
        while pending:
            for j in np.flatnonzero(self.matrices[candidate, pending.pop()] > 0):
                if j not in visited:
                    visited.add(j)
                    pending.append(j)
        block = self.matrices[candidate][np.ix_(sorted(visited), sorted(visited))]
        abscissa = max(eigvals(block / self.network.denominator).real)
        self.assertLessEqual(float(lo), abscissa + 1e-14)
        self.assertLessEqual(abscissa, float(hi) + 1e-14)
        self.assertLess(float(hi - lo), 1e-8)


if __name__ == "__main__":
    unittest.main()
