from dataclasses import replace
from functools import lru_cache
from itertools import product
from math import comb, log, sqrt
import unittest

import numpy as np
from scipy.linalg import expm

from purposeful.closed_loop import (
    Device,
    STATIONARY_PRIOR,
    best_fixed_frontier,
    block_failure_bound,
    evaluate_fixed_schedule,
    minimum_reserve,
    optimal_policy,
)
from purposeful.reserve import RepairContact


def tree_survival(device, rounds, reserve, previous, schedule=None, policy=None):
    """Enumerate executable decisions before observations and both live outcomes."""
    quantum = device.repair.token_energy
    max_repair = expm(device.repair.generator() * device.repair.duration)[0, 1]

    @lru_cache(None)
    def visit(left, battery, last):
        if left == 0:
            return 1.0
        p = device.opportunity
        if last != STATIONARY_PRIOR:
            p = (1 - device.persistence) * p + device.persistence * last
        affordable = battery >= device.sensing_cost
        options = [False, True] if affordable else [False]
        if schedule is not None:
            options = [schedule[rounds - left] and affordable]
        if policy is not None:
            options = [bool(policy.sense[left, last, battery])]
        totals = []
        for sensing in options:
            sensing_spend = device.sensing_cost if sensing else 0
            after_sensing = battery - sensing_spend
            acts = range(2) if after_sensing >= device.actuation_cost else range(1)
            maps = list(product(acts, repeat=2 if sensing else 1))
            if policy is not None:
                recorded = policy.action_after_readout[left, last, battery] if sensing else (
                    policy.blind_action[left, last, battery],
                )
                assert tuple(recorded) in maps
                maps = [tuple(recorded)]
            for actions in maps:
                total = 0.0
                for flag, chance in enumerate((1 - p, p)):
                    for readout, action in enumerate(actions):
                        likelihood = 1.0
                        if sensing:
                            likelihood = device.fidelity if readout == flag else 1 - device.fidelity
                        actuator_spend = int(action) * device.actuation_cost
                        after_action = after_sensing - actuator_spend
                        if flag == action:
                            harvest = min(device.harvest, device.capacity - after_action) if action else 0
                            repair_work, survival = 0, 1
                        elif after_action:
                            harvest, repair_work, survival = 0, 1, max_repair
                        else:
                            continue
                        next_battery = after_action + harvest - repair_work
                        accounted = (next_battery + sensing_spend + actuator_spend
                                     + repair_work - harvest) * quantum
                        assert abs(battery * quantum - accounted) < 1e-12
                        total += chance * likelihood * survival * visit(left - 1, next_battery, flag)
                totals.append(total)
        return max(totals)

    return visit(rounds, reserve, previous)


class ClosedLoopTests(unittest.TestCase):
    def setUp(self):
        self.device = Device(2, 1, 1, 2, 0.5, RepairContact(log(3), log(3), 1.0))

    def test_exact_two_round_gap(self):
        policy = optimal_policy(self.device, 2)
        self.assertAlmostEqual(policy.survival[2, STATIONARY_PRIOR, 2], 7 / 8)
        self.assertFalse(policy.sense[1, 0, 1])
        self.assertTrue(policy.sense[1, 1, 2])
        for schedule, correct in zip(product([False, True], repeat=2), [9/16, 5/8, 3/4, 3/4]):
            self.assertAlmostEqual(evaluate_fixed_schedule(self.device, schedule)[STATIONARY_PRIOR, 2], correct)
        self.assertEqual(minimum_reserve(policy, 0.8), 2)
        self.assertIsNone(minimum_reserve(policy, 0.9))

    def test_exhaustive_tree_and_attaining_controller(self):
        device = replace(self.device, capacity=3, harvest=1, fidelity=0.82, persistence=0.6)
        rounds = 3
        policy = optimal_policy(device, rounds)
        envelope, schedules, _ = best_fixed_frontier(device, rounds)
        for previous, reserve in np.ndindex(envelope.shape):
            with self.subTest(previous=previous, reserve=reserve):
                exact = tree_survival(device, rounds, reserve, previous)
                attaining = tree_survival(device, rounds, reserve, previous, policy=policy)
                best_fixed = max(
                    tree_survival(device, rounds, reserve, previous, schedule=schedule)
                    for schedule in product([False, True], repeat=rounds)
                )
                self.assertAlmostEqual(policy.survival[-1, previous, reserve], exact)
                self.assertAlmostEqual(attaining, exact)
                self.assertAlmostEqual(envelope[previous, reserve], best_fixed)
                self.assertAlmostEqual(
                    tree_survival(device, rounds, reserve, previous, schedule=schedules[previous, reserve]),
                    best_fixed,
                )

    def test_noise_threshold_and_positive_gap_with_persistence(self):
        for fidelity, strict in [(0.6, False), (2 / 3, False), (0.7, True), (0.9, True)]:
            device = replace(self.device, fidelity=fidelity)
            adaptive = optimal_policy(device, 2).survival[-1, STATIONARY_PRIOR, 2]
            fixed = best_fixed_frontier(device, 2)[0][STATIONARY_PRIOR, 2]
            self.assertEqual(adaptive > fixed + 1e-12, strict)
        device = replace(self.device, fidelity=0.9, persistence=0.5)
        self.assertAlmostEqual(optimal_policy(device, 2).survival[-1, STATIONARY_PRIOR, 2]
                               - best_fixed_frontier(device, 2)[0][STATIONARY_PRIOR, 2], 0.028125)

    def test_zero_reserve_monotonic_frontier_and_block_bound(self):
        device = replace(self.device, capacity=4, fidelity=0.87, persistence=0.4)
        length, failure = block_failure_bound(device)
        frontier = optimal_policy(device, 3 * length).survival
        self.assertTrue(np.all(np.diff(frontier, axis=-1) >= -1e-14))
        self.assertTrue(np.all(frontier[length] <= 1 - failure))
        self.assertTrue(np.all(frontier[-1] <= (1 - failure) ** 3))
        self.assertAlmostEqual(frontier[4, STATIONARY_PRIOR, 0], 0.5 * 0.7 ** 3)

    def test_long_horizon_construction_and_certificate(self):
        rounds = 20
        rate = (1 + sqrt(5)) / 4
        frontier = optimal_policy(self.device, rounds).survival[:, STATIONARY_PRIOR]
        self.assertTrue(np.all(frontier[-1] <= rate ** rounds * np.array([1, 1, 2 * rate])))
        self.assertGreaterEqual(frontier[-1, 2], rate ** rounds)
        repeating = evaluate_fixed_schedule(self.device, [True, False] * (rounds // 2))
        self.assertGreaterEqual(repeating[STATIONARY_PRIOR, 2], 2 ** (-rounds / 2))
        self.assertGreater(frontier[-1, 2], 10 * best_fixed_frontier(self.device, rounds)[0][STATIONARY_PRIOR, 2])

        strong = replace(self.device, repair=RepairContact(log(9), log(9), 1.0))
        kappa = strong.repair.success
        root_rate = (1 + sqrt(1 + 4 * kappa)) / 4
        self.assertGreater(root_rate, 1 / sqrt(2))
        blind = evaluate_fixed_schedule(strong, [False] * rounds)
        self.assertGreaterEqual(blind[STATIONARY_PRIOR, 2], root_rate ** rounds)

    def test_no_harvest_bounds(self):
        rounds, rho = 25, 0.3
        device = replace(self.device, harvest=0, persistence=rho, fidelity=0.9)
        p00 = (1 + rho) / 2
        lower = 0.5 * p00 ** (rounds - 1)
        upper = sum(comb(rounds, j) * p00 ** (rounds - 1 - 2 * j)
                    for j in range(device.capacity + 1))
        adaptive = optimal_policy(device, rounds).survival[-1, STATIONARY_PRIOR, 2]
        self.assertLessEqual(lower, adaptive)
        self.assertLessEqual(adaptive, upper)


if __name__ == "__main__":
    unittest.main()
