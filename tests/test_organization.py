"""Independent channel enumeration for the spatial organization model."""

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import comb, exp, log
import unittest

import numpy as np
from scipy.linalg import expm

from purposeful.organization import Organization, build_system


@dataclass(frozen=True)
class Channel:
    destination: tuple
    prefactor: Fraction
    exponent: Fraction
    event: int

    def rate(self):
        return float(self.prefactor) * exp(float(self.exponent))


# Opposite event indices denote opposite directions in the physical ledger.
DAMAGE, REPAIR, PASSIVE, PAID = 1, 2, 3, 4


def rational(value):
    """Treat declared decimal model parameters as exact rational numbers."""
    return Fraction(str(value))


def reference_graph(model, n, *, folded, killed):
    """Enumerate per-edge rates as an exact rational times an exponential."""
    canonical = lambda x: tuple(sorted(min(v, n - v) for v in x))
    positions = tuple(product(range(n), repeat=3))
    orbits = {}
    for position in positions:
        representative = canonical(position) if folded else position
        orbits[representative] = orbits.get(representative, 0) + 1
    sites = sorted(orbits)
    lengths = [(min(j, n - j) * rational(model.length) / n) ** 2 for j in range(n)]
    radius_squared = rational(model.radius) ** 2
    offset, depth = rational(model.offset), rational(model.depth)
    geometry = {}
    for position in sites:
        core = max(Fraction(0), 1 - sum(lengths[j] for j in position) / radius_squared)
        geometry[position] = (core**4, offset - depth * core**2)
    neighbors = {}
    for position in sites:
        values = []
        for axis, direction in product(range(3), (-1, 1)):
            result = list(position)
            result[axis] = (result[axis] + direction) % n
            values.append(canonical(result) if folded else tuple(result))
        neighbors[position] = values
    diffusional = rational(model.diffusivity) * (n / rational(model.length)) ** 2
    eta, kappa = rational(model.eta), rational(model.activity)
    passive, paid = rational(model.passive), rational(model.paid)
    inventory = model.fuel
    states = [(h, b, tethered, position)
              for h, b, tethered, position in product(range(1 if killed else 0, 3),
                                                       range(inventory + 1), (0, 1), sites)]
    start = [Fraction(orbits[position], n**3) if (h, b, q) == (1, inventory, 0) else Fraction(0)
             for h, b, q, position in states]
    out = []
    for h, b, q, position in states:
        gate, potential = geometry[position]
        branches = []

        def add(destination, prefactor, event, exponent=Fraction(0)):
            if prefactor:
                branches.append(Channel(destination, prefactor, exponent, event))

        for target in neighbors[position]:
            difference = q * (geometry[target][1] - potential)
            add((h, b, q, target), diffusional, 0, -difference / 2)
        for label in (0, 1):
            healthy = label < h
            other_healthy = 1 - label < h
            if healthy:
                add((h - 1, b, q, position), eta, DAMAGE)
                if other_healthy and b < inventory:
                    add((h - 1, b + 1, q, position), kappa * gate * (inventory - b) / 4, -REPAIR)
            else:
                add((h + 1, b, q, position), eta / 2, -DAMAGE)
                if other_healthy and b:
                    add((h + 1, b - 1, q, position), kappa * gate * b, REPAIR)
        toggle_exponent = min(Fraction(0), potential if q else -potential)
        add((h, b, 1 - q, position), passive * gate, -PASSIVE if q else PASSIVE, toggle_exponent)
        if q and b < inventory:
            add((h, b + 1, 0, position), paid * gate * (inventory - b) / 8, -PAID, toggle_exponent)
        if not q and b:
            add((h, b - 1, 1, position), paid * gate * b, PAID, toggle_exponent)
        out.append(branches)
    return states, out, start, geometry, orbits


def policy_multiplier(policy, fuel, state, event):
    if abs(event) != PAID:
        return Fraction(1)
    if policy != "spectator":
        return rational(policy)
    reactant_fuel = state[1] if event == PAID else state[1] + 1
    return Fraction(reactant_fuel - 1, max(1, fuel - 1))


def dense_reference(model, n, policy, *, folded=True, killed=True):
    states, channels, start, geometry, orbits = reference_graph(model, n, folded=folded, killed=killed)
    matrix = np.zeros((len(states), len(states)))
    index = {state: i for i, state in enumerate(states)}
    for i, (state, branches) in enumerate(zip(states, channels)):
        for branch in branches:
            rate = branch.rate() * float(policy_multiplier(policy, model.fuel, state, branch.event))
            matrix[i, i] -= rate
            if branch.destination in index:
                matrix[i, index[branch.destination]] += rate
    return states, matrix, np.asarray(start, dtype=float), geometry, orbits


class OrganizationTests(unittest.TestCase):
    def test_periodic_reference_survival(self):
        model = Organization(fuel=2)
        n = 4
        system = build_system(model, n)
        for policy in (0., .61, "spectator"):
            _, matrix, launch, _, _ = dense_reference(model, n, policy, folded=False)
            expected = launch @ expm(.33 * matrix) @ np.ones(len(launch))
            np.testing.assert_allclose(system.survival(policy, [0., .33]), [1, expected], atol=5e-13)

    def test_channel_matrix_and_thermal_equilibrium(self):
        model = Organization()
        n = 4
        system = build_system(model, n, killed=False)
        for policy in (0., .61, 1., "spectator"):
            states, matrix, _, geometry, orbits = dense_reference(model, n, policy, killed=False)
            expected = np.array([comb(2, h) * 2**(2 - h) * comb(model.fuel, b) * 8.**(-b)
                                 * orbits[pos] * exp(-q * float(geometry[pos][1]))
                                 for h, b, q, pos in states])
            expected /= expected.sum()
            actual = system.generator(policy).toarray()
            np.testing.assert_allclose(actual, matrix, atol=5e-13)
            np.testing.assert_allclose(system.equilibrium(), expected, atol=3e-15)
            flow = expected[:, None] * actual
            np.testing.assert_allclose(flow, flow.T, atol=5e-15)

    def test_first_collapse_preserves_killing_flux(self):
        model = Organization()
        full = build_system(model, 4, killed=False)
        killed = build_system(model, 4)
        site_count = len(full.grid.weight)
        rows = np.concatenate([i * site_count + np.arange(site_count)
                               for i, (h, _, _) in enumerate(full.states) if h])
        for policy in (0., .5, 1., "spectator"):
            whole = full.generator(policy).toarray()
            live = killed.generator(policy).toarray()
            np.testing.assert_allclose(live, whole[np.ix_(rows, rows)], atol=1e-13)
            np.testing.assert_allclose(live.sum(axis=1), np.repeat(
                [-model.eta if h == 1 else 0 for h, _, _ in killed.states], site_count), atol=1e-13)
        for system in (full, killed):
            start = np.zeros((len(system.states), site_count))
            start[system.states.index((1, model.fuel, 0))] = system.grid.weight
            np.testing.assert_array_equal(system.initial(), start.ravel())

    def test_ledger_and_reverse_gate(self):
        model = Organization()
        states, branches, _, geometry, _ = reference_graph(model, 4, folded=True, killed=False)
        index = dict(zip(states, branches))
        checked = set()
        for state, row in zip(states, branches):
            h, b, q, pos = state
            for channel in row:
                event = channel.event
                other_h, other_b, other_q, other_pos = channel.destination
                sign = 1 if event >= 0 else -1
                damage, repair, passive, paid = [sign * (abs(event) == kind)
                                                for kind in (DAMAGE, REPAIR, PASSIVE, PAID)]
                self.assertEqual(other_h - h, -damage + repair)
                self.assertEqual(other_b - b, -repair - paid)
                self.assertEqual(other_q - q, passive + paid)
                prior_energy = b * log(8) - (2 - h) * log(2) + q * float(geometry[pos][1])
                after_energy = (other_b * log(8) - (2 - other_h) * log(2)
                                + other_q * float(geometry[other_pos][1]))
                expected = (-damage * log(2) - repair * log(4) - paid * log(8)
                            + other_q * float(geometry[other_pos][1]) - q * float(geometry[pos][1]))
                self.assertAlmostEqual(after_energy - prior_energy, expected, places=13)
                reverse = [entry for entry in index[channel.destination]
                           if entry.destination == state and entry.event == -event]
                self.assertTrue(reverse)
                self.assertEqual(policy_multiplier("spectator", model.fuel, state, event),
                                 policy_multiplier("spectator", model.fuel, channel.destination, -event))
                if paid:
                    checked.add((b if paid > 0 else other_b, policy_multiplier("spectator", model.fuel,
                                                                               state, event)))
        self.assertEqual(checked, {(1, 0), (2, Fraction(1, 3)), (3, Fraction(2, 3)), (4, 1)})

    def test_zero_paid_cap_removes_policy_dependence(self):
        model = Organization(paid=0.)
        system = build_system(model, 4)
        t = np.array([0., .1, 1., 2.])
        expected = system.survival(0., t)
        for policy in (1., "spectator"):
            np.testing.assert_array_equal(system.generator(policy).toarray(), system.base.toarray())
            np.testing.assert_allclose(system.survival(policy, t), expected, atol=8e-13)

    def test_stopped_resource_expectations_match_integrated_gross_events(self):
        model = Organization()
        system = build_system(model, 6)
        times = [0., .2, .8, 2.]
        for policy in (0., .6, 'spectator'):
            o = system.observables(policy, times)
            paid = o['paid_forward'] - o['paid_reverse']
            passive = o['passive_forward'] - o['passive_reverse']
            repair = o['repair_forward'] - o['repair_reverse']
            np.testing.assert_allclose(o['stopped_fuel'], model.fuel - paid - repair, atol=5e-12)
            np.testing.assert_allclose(o['stopped_tether'], paid + passive, atol=5e-12)


if __name__ == "__main__":
    unittest.main()
