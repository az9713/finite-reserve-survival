import unittest
import numpy as np
from scipy.linalg import expm
from purposeful.model import (Sensor, sensor_for_accuracy, frontier_power,
    four_state_thermodynamics, hazard_generator, hazard_survival, hazard_decay,
    damage_survival, damage_generator, simulate_survival)


class ModelTests(unittest.TestCase):
    def test_local_detailed_balance_and_stationarity(self):
        for r in [.01, .1, 1., 10.]:
            for k in [.03, 1., 20.]:
                for affinity in [0., .2, 3., 12.]:
                    s = Sensor(r, k, affinity)
                    self.assertAlmostEqual(np.log(s.a/s.b), affinity)
                    np.testing.assert_allclose(s.stationary() @ s.generator(), 0, atol=2e-14)
                    np.testing.assert_allclose(s.generator().sum(axis=1), 0, atol=2e-14)

    def test_flux_accounting_independent_of_closed_form(self):
        rng = np.random.default_rng(81)
        for _ in range(60):
            s = Sensor(10**rng.uniform(-2, 1), 10**rng.uniform(-2, 2), rng.uniform(0, 9))
            f = four_state_thermodynamics(s)
            np.testing.assert_allclose(f['stationary'], s.stationary(), atol=2e-13)
            for key in ['heat', 'work', 'entropy_production']:
                self.assertAlmostEqual(f[key], s.power, places=10)
            self.assertAlmostEqual(f['learning_rate'], s.learning_rate, places=10)
            self.assertLessEqual(s.learning_rate, s.power+1e-12)

    def test_zero_affinity(self):
        s = Sensor(1, 3, 0)
        self.assertEqual(s.q, .5)
        self.assertEqual(s.power, 0)
        self.assertEqual(s.mutual_information, 0)

    def test_inverse_and_frontier(self):
        q = .2; r = .1; cap = 6
        best = frontier_power(r, cap, q)
        for k in np.geomspace(.31, cap, 50):
            s = sensor_for_accuracy(r, k, q)
            self.assertAlmostEqual(s.q, q)
            self.assertGreaterEqual(s.power+1e-12, best)

    def test_infeasible_target(self):
        with self.assertRaises(ValueError): sensor_for_accuracy(.1, .1, .2)

    def test_frontier_accuracy_monotonicity(self):
        costs = [frontier_power(.1, 10, q) for q in [.4, .3, .2, .1]]
        self.assertTrue(np.all(np.diff(costs)>0))

    def test_no_hazard_no_death(self):
        self.assertAlmostEqual(hazard_survival(Sensor(.1, 1, 2), 100, h=0), 1)

    def test_hazard_aggregation_matches_full_generator(self):
        s = Sensor(.2, 2, 3)
        g = s.generator()-np.diag([0, 1.3, 1.3, 0])
        for t in [0, .1, 4, 30]:
            full = s.stationary() @ expm(g*t) @ np.ones(4)
            self.assertAlmostEqual(full, hazard_survival(s, t, 1.3), places=12)

    def test_jensen_lower_bound(self):
        for k in [.4, 1, 5.8]:
            s = sensor_for_accuracy(.1, k, .2)
            for t in [.1, 1, 20]:
                self.assertGreaterEqual(hazard_survival(s,t)+1e-12, np.exp(-s.q*t))

    def test_decay_eigenvalue_and_monotonicity(self):
        last = -1
        for k in np.geomspace(.31, 30, 30):
            s = sensor_for_accuracy(.1, k, .2)
            eigen = -max(np.linalg.eigvals(hazard_generator(s)).real)
            self.assertAlmostEqual(eigen, hazard_decay(s), places=12)
            self.assertGreater(hazard_decay(s), last)
            last = hazard_decay(s)

    def test_ranking_reversal(self):
        slow = sensor_for_accuracy(.1, .4, .2)
        fast = sensor_for_accuracy(.1, 5.8, .2)
        self.assertLess(fast.power, slow.power)
        self.assertAlmostEqual(slow.mutual_information, fast.mutual_information)
        self.assertGreater(hazard_survival(slow,20,h=.5), hazard_survival(fast,20,h=.5))
        self.assertLess(damage_survival(slow,20), damage_survival(fast,20))

    def test_one_stage_damage_equals_hazard(self):
        s = Sensor(.3, 4, 1.5)
        for t in [.2, 3, 20]:
            self.assertAlmostEqual(damage_survival(s,t,stages=1,duration=2),
                                   hazard_survival(s,t,h=.5), places=12)

    def test_gillespie_against_matrix_exponential(self):
        s = sensor_for_accuracy(.1,.4,.2)
        g = hazard_generator(s)
        result = simulate_survival(g,[.8,.2],5,samples=10000,seed=123)
        exact = hazard_survival(s,5)
        self.assertLess(abs(result['estimate']-exact), 5*np.sqrt(exact*(1-exact)/10000))

    def test_damage_probability_and_initial_preparation(self):
        s = Sensor(.1,1,2)
        for q0 in [0., .2, 1.]:
            self.assertAlmostEqual(damage_survival(s,0,initial_q=q0), 1)
            self.assertTrue(0 <= damage_survival(s,20,initial_q=q0) <= 1)

    def test_invalid_parameters(self):
        for args in [(0,1,2), (1,0,2), (1,1,-1), (1,1,float('nan'))]:
            with self.assertRaises(ValueError): Sensor(*args)
        with self.assertRaises(ValueError): hazard_survival(Sensor(1,1,1),-1)
        with self.assertRaises(ValueError): damage_generator(Sensor(1,1,1), stages=0)


if __name__ == '__main__':
    unittest.main()
