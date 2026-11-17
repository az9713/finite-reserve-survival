"""Independent thermodynamic, first-passage, and enclosure checks."""
from fractions import Fraction as F
from math import exp, factorial, log
import unittest

import numpy as np
from scipy.linalg import expm

from purposeful.assembly import (Assembly, joint_interval, maintenance_interval,
                                one_fuel_joint, success_interval, curvature_bounds)


class AssemblyTests(unittest.TestCase):
    def test_stoichiometry_and_detailed_balance(self):
        model = Assembly(fuel=4)
        for name, s, t, forward, reverse in model.edges(F(1, 2), F(1, 4)):
            self.assertGreater(forward, 0)
            self.assertGreater(reverse, 0)
            self.assertEqual(model.equilibrium_weight(s, F(1, 4)) * forward,
                             model.equilibrium_weight(t, F(1, 4)) * reverse)
            def counts(state):
                h, d, b = state
                return np.array([model.pairs-h-d, model.pairs-h-d, h, d, b, model.fuel-b])
            change = counts(t)-counts(s)
            expected = {'assembly': [-1, -1, 1, 0, -1, 1],
                        'damage': [0, 0, -1, 1, 0, 0],
                        'repair': [0, 0, 1, -1, -1, 1]}[name]
            np.testing.assert_array_equal(change, expected)
            self.assertEqual(sum(change[[0, 2, 3]]), 0)
            self.assertEqual(sum(change[[1, 2, 3]]), 0)
            self.assertEqual(sum(change[[4, 5]]), 0)

    def test_full_ensemble_stationarity_and_entropy_decrease(self):
        model = Assembly(fuel=3)
        states, generator = model.full_generator(F(2, 3), F(1, 3))
        weights = np.array([float(model.equilibrium_weight(s, F(1, 3))) for s in states])
        equilibrium = weights / weights.sum()
        np.testing.assert_allclose(equilibrium @ generator, 0, atol=1e-15)
        initial = np.ones(len(states)) / len(states)
        divergences = []
        for time in [0, .1, 1, 5]:
            distribution = initial @ expm(time * generator)
            divergences.append(np.sum(distribution * np.log(distribution / equilibrium)))
        self.assertTrue(all(a >= b for a, b in zip(divergences, divergences[1:])))

    def test_sharp_resource_and_energy_constructions(self):
        model = Assembly(pairs=3, fuel=5)
        state = (0, 0, 5)
        ledger_a = ledger_r = 0
        sequence = ['assembly']*3 + ['damage','repair']*2 + ['damage']*3
        for reaction in sequence:
            candidates = [edge for edge in model.edges(F(1, 2), F(1, 4))
                          if edge[0] == reaction and edge[1] == state]
            self.assertEqual(len(candidates), 1)
            _, _, state, forward, _ = candidates[0]
            self.assertGreater(forward, 0)
            ledger_a += reaction == 'assembly'
            ledger_r += reaction == 'repair'
            h,d,b = state
            self.assertEqual(ledger_a, h+d)
            self.assertEqual(ledger_a+ledger_r, 5-b)
        self.assertEqual(state, (0,3,0))
        energy = -log(4)*3-log(2)*3
        self.assertAlmostEqual(5*log(8)-energy, 5*log(8)+3*(log(4)+log(2)))

    def test_launch_and_one_fuel_closed_form(self):
        for pairs in [1, 2, 3]:
            model = Assembly(pairs=pairs, fuel=1)
            for alpha in [F(0), F(1, 4), F(1)]:
                for z in [F(0), F(1, 2), F(1)]:
                    for time in [0, .3, 2]:
                        self.assertAlmostEqual(model.joint(time, alpha, z),
                                               one_fuel_joint(model, time, alpha, z), places=12)
        # Equal-rate limit of the formula.
        model = Assembly(pairs=1, fuel=1, damage_rate=F(7, 8))
        self.assertAlmostEqual(one_fuel_joint(model, 2, 1, 1), 2*exp(-2))

    def test_independent_block_convolution(self):
        model, alpha, z = Assembly(fuel=4), F(1, 2), F(1, 4)
        states, integers, denominator = model.killed_generator(alpha, z)
        generator = integers / denominator
        arrival = float(alpha * model.assembly_rate * model.pairs**2 * model.fuel)
        from scipy.integrate import quad
        start = states.index((1, 0, model.fuel - 1))-1
        integral = quad(lambda u: arrival*np.exp(-arrival*u)
                        * expm((2-u)*generator[1:, 1:])[start].sum(), 0, 2,
                        epsabs=1e-12)[0]
        self.assertAlmostEqual(integral, model.joint(2, alpha, z), places=11)

    def test_rational_enclosures_and_curvature(self):
        model = Assembly(fuel=4)
        for alpha, z in [(F(1, 2), F(1, 4)), (F(1), F(0)), (F(0), F(1))]:
            lo, hi = joint_interval(model, 2, alpha, z)
            value = model.joint(2, alpha, z)
            self.assertLessEqual(float(lo), value)
            self.assertGreaterEqual(float(hi), value)
            lo, hi = success_interval(model, 1, 4, alpha, z)
            value = model.formation_and_maintenance(1, 4, alpha, z)
            self.assertLessEqual(float(lo), value)
            self.assertGreaterEqual(float(hi), value)
        box = (F(1, 4), F(3, 4), F(1, 8), F(1, 2))
        ca, cz = curvature_bounds(model, 1, 4, box)
        f = lambda a,z: model.formation_and_maintenance(1, 4, a, z)
        step = .0001
        a, z = .5, .3
        self.assertLess(abs(f(a+step,z)-2*f(a,z)+f(a-step,z))/step**2, float(ca))
        self.assertLess(abs(f(a,z+step)-2*f(a,z)+f(a,z-step))/step**2, float(cz))

    def test_single_complex_sharp_allocation(self):
        model = Assembly(pairs=1, fuel=3)
        deadline, lifetime, z = 1, 16, F(1)
        optimum = log(1+8*model.fuel*deadline/(float(z)*lifetime)) / (model.fuel*deadline)
        best = model.formation_and_maintenance(deadline, lifetime, optimum, z)
        for alpha in np.linspace(0, 1, 101):
            self.assertLessEqual(model.formation_and_maintenance(deadline, lifetime, alpha, z), best+1e-12)
        for alpha in [F(1, 4), F(1)]:
            expected = (1-exp(-float(alpha)*3))*exp(-(1+float(alpha)*float(z)/8)*16)
            self.assertAlmostEqual(model.formation_and_maintenance(1, 16, alpha, z), expected)

    def test_zero_resources_and_invalid_parameters(self):
        self.assertEqual(Assembly(fuel=0).formation_and_maintenance(1, 4, .5, .5), 0)
        self.assertEqual(Assembly(pairs=0).joint(2, .5, .5), 0)
        for kwargs in [{'fuel': -1}, {'pairs': True}, {'assembly_rate': 0}, {'fuel_factor': 1}]:
            with self.assertRaises(ValueError):
                Assembly(**kwargs)
        with self.assertRaises(ValueError):
            Assembly().joint(-1, .5, .5)
        with self.assertRaises(ValueError):
            Assembly().joint(1, 2, .5)


if __name__ == '__main__':
    unittest.main()
