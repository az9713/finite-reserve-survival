import itertools
import math
import unittest

import numpy as np
from scipy.linalg import expm

from purposeful.reserve import (
    RepairContact,
    minimum_energy_for_known_stress_count,
    minimum_tokens,
    stress_count_distribution,
    stress_transition,
    survival_frontier,
)


def full_checkpoint_survival(rounds, q, mixing_rate, interval, contact, reserve):
    """Independent propagation over stress, damage, and explicit battery states."""
    transition = stress_transition(q, mixing_rate, interval)
    states = list(itertools.product(range(2), range(2), range(reserve + 1)))
    index = {state: i for i, state in enumerate(states)}
    size = len(states)
    load = np.zeros((size, size))
    repair = np.zeros_like(load)
    pair = contact.generator()
    forward, reverse = pair[0, 1], pair[1, 0]
    live = np.zeros(size, dtype=bool)
    for i, (stress, damage, fuel) in enumerate(states):
        live[i] = damage == 0
        for next_stress in range(2):
            j = index[next_stress, damage ^ next_stress, fuel]
            load[i, j] = transition[stress, next_stress]
        if stress and damage and fuel > 0:
            repair[i, index[1, 0, fuel - 1]] = forward
        if stress and not damage and fuel < reserve:
            repair[i, index[1, 1, fuel + 1]] = reverse
        repair[i, i] = -repair[i].sum()
    step = (load @ expm(repair * contact.duration)) * live
    initial = np.zeros(size)
    initial[index[0, 0, reserve]], initial[index[1, 0, reserve]] = 1 - q, q
    return initial @ np.linalg.matrix_power(step, rounds) @ np.ones(size)


class ReserveTests(unittest.TestCase):
    def test_reversible_paid_repair_and_rate_cap(self):
        contact = RepairContact(token_energy=1.3, rate_ceiling=2.2, duration=0.75)
        generator = contact.generator()
        self.assertAlmostEqual(math.log(generator[0, 1] / generator[1, 0]), 1.3)
        self.assertAlmostEqual(generator[0, 1] + generator[1, 0], 2.2)
        self.assertAlmostEqual(expm(generator * contact.duration)[0, 1], contact.success)

        throttled = expm(contact.generator(0.7) * 0.3) @ expm(generator * 0.45)
        self.assertLess(throttled[0, 1], contact.success)

    def test_frontier_matches_full_state_propagation(self):
        contact = RepairContact(1.7, 3.1, 1.0)
        for q, mixing_rate in [(0.0, 1.0), (0.2, 0.6), (0.3, 6.0), (1.0, 2.0)]:
            for rounds in [0, 2, 5]:
                distribution = stress_count_distribution(rounds, q, mixing_rate, 0.9)
                frontier = survival_frontier(distribution, contact)
                exact = [
                    full_checkpoint_survival(rounds, q, mixing_rate, 0.9, contact, reserve)
                    for reserve in range(rounds + 1)
                ]
                np.testing.assert_allclose(frontier, exact, rtol=1e-13, atol=1e-14)

    def test_reserve_reverses_ranking(self):
        contact = RepairContact(math.log(9), math.log(10), 1.0)
        slow = survival_frontier(stress_count_distribution(2, 0.2, 0.6, 1.0), contact)
        fast = survival_frontier(stress_count_distribution(2, 0.2, 6.0, 1.0), contact)
        self.assertGreater(slow[0], fast[0])
        self.assertLess(slow[1], fast[1])
        self.assertGreater(slow[2], fast[2])

    def test_minimum_reserve_and_infeasibility(self):
        distribution = stress_count_distribution(2, 0.2, 6.0, 1.0)
        contact = RepairContact(math.log(9), math.log(10), 1.0)
        frontier = survival_frontier(distribution, contact)
        self.assertEqual(minimum_tokens(distribution, contact, frontier[0]), 0)
        self.assertEqual(minimum_tokens(distribution, contact, 0.88), 1)
        self.assertEqual(minimum_tokens(distribution, contact, 0.90), 2)
        self.assertIsNone(minimum_tokens(distribution, contact, 0.95))

    def test_known_count_energy_bound_and_construction(self):
        stresses, target, rate, duration = 3, 0.4, 2.0, 1.0
        total = minimum_energy_for_known_stress_count(stresses, target, rate, duration)
        equal = RepairContact(total / stresses, rate, duration)
        self.assertAlmostEqual(equal.success ** stresses, target)
        unequal = math.prod(
            RepairContact(total * fraction, rate, duration).success
            for fraction in [0.1, 0.2, 0.7]
        )
        self.assertLess(unequal, target)
        self.assertEqual(minimum_energy_for_known_stress_count(0, 1, rate, duration), 0)
        self.assertEqual(minimum_energy_for_known_stress_count(3, 0.01, rate, duration), 0)
        self.assertTrue(math.isinf(minimum_energy_for_known_stress_count(3, 0.9, rate, duration)))


if __name__ == "__main__":
    unittest.main()
