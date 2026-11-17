"""Closed, reversible assembly of catalysts that repair other assembled catalysts.

Counts are (healthy complexes, damaged complexes, fuel molecules). Free
complementary monomers and waste follow from the two conserved inventories.
All propensities use stochastic mass action at a fixed reference volume.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial, lcm

import numpy as np
from scipy.linalg import expm

from .autonomous import _poisson_weights, _ceil


@dataclass(frozen=True)
class Assembly:
    pairs: int = 2
    fuel: int = 3
    assembly_rate: Fraction = Fraction(1)
    repair_rate: Fraction = Fraction(4)
    damage_rate: Fraction = Fraction(1)
    fuel_factor: int = 8
    damage_factor: int = 2

    def __post_init__(self):
        for name in ('pairs', 'fuel'):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(name + ' must be a nonnegative integer')
        if not isinstance(self.fuel_factor, int) or not isinstance(self.damage_factor, int):
            raise ValueError('energy factors must be integers')
        if not self.fuel_factor > self.damage_factor >= 1:
            raise ValueError('energy factors must satisfy fuel > damage >= 1')
        for name in ('assembly_rate', 'repair_rate', 'damage_rate'):
            value = Fraction(getattr(self, name))
            if value <= 0:
                raise ValueError(name + ' must be positive')
            object.__setattr__(self, name, value)

    def states(self):
        return [(h, d, b) for h in range(self.pairs + 1)
                for d in range(self.pairs - h + 1) for b in range(self.fuel + 1)]

    def edges(self, alpha, inverse_binding):
        """Yield distinct reversible reaction pairs, including collapse returns."""
        alpha, z = Fraction(alpha), Fraction(inverse_binding)
        if not 0 <= alpha <= 1 or not 0 <= z <= 1:
            raise ValueError('alpha and inverse binding must be in [0,1]')
        a = self.assembly_rate * alpha
        r = self.repair_rate * (1 - alpha) * z
        for state in self.states():
            h, d, b = state
            free, waste = self.pairs - h - d, self.fuel - b
            if free and b:
                yield ('assembly', state, (h + 1, d, b - 1),
                       a * free ** 2 * b,
                       a * z * (h + 1) * (waste + 1) / self.fuel_factor)
            if h:
                yield ('damage', state, (h - 1, d + 1, b),
                       self.damage_rate * h,
                       self.damage_rate * (d + 1) / self.damage_factor)
            if h and d and b:
                yield ('repair', state, (h + 1, d - 1, b - 1),
                       r * h * d * b,
                       r * Fraction(self.damage_factor, self.fuel_factor) * h * (h + 1) * (waste + 1))

    def equilibrium_weight(self, state, inverse_binding):
        h, d, b = state
        free, waste = self.pairs - h - d, self.fuel - b
        z = Fraction(inverse_binding)
        if z <= 0:
            raise ValueError('finite binding requires positive inverse binding')
        return (z ** (-h - d) * self.damage_factor ** d * Fraction(1, self.fuel_factor ** b)
                / (factorial(free) ** 2 * factorial(h) * factorial(d) * factorial(b) * factorial(waste)))

    def full_generator(self, alpha, inverse_binding):
        states = self.states()
        index = {s: i for i, s in enumerate(states)}
        matrix = np.zeros((len(states), len(states)))
        for _, s, t, forward, reverse in self.edges(alpha, inverse_binding):
            for start, end, rate in [(s, t, forward), (t, s, reverse)]:
                i, j = index[start], index[end]
                matrix[i, j] += float(rate)
                matrix[i, i] -= float(rate)
        return states, matrix

    def killed_generator(self, alpha, inverse_binding):
        """Keep launch and post-assembly live states. A return to h=0 kills."""
        launch = (0, 0, self.fuel)
        live = [state for state in self.states() if state[0] > 0]
        states = [launch, *live]
        index = {state: i for i, state in enumerate(states)}
        matrix = [[Fraction(0) for _ in states] for _ in states]
        for _, s, t, forward, reverse in self.edges(alpha, inverse_binding):
            for start, end, rate in [(s, t, forward), (t, s, reverse)]:
                if start not in index:
                    continue
                i = index[start]
                matrix[i][i] -= rate
                # Launch is a history-labelled state. Post-launch paths may
                # never re-enter it, even if their physical counts coincide.
                if end in index and end != launch:
                    matrix[i][index[end]] += rate
        denominator = lcm(*(value.denominator for row in matrix for value in row))
        integers = [[int(value * denominator) for value in row] for row in matrix]
        maximum = max(sum(abs(value) for value in row) for row in integers)
        if maximum > np.iinfo(np.int64).max:
            raise ValueError('integer generator is too large')
        return states, np.array(integers, dtype=np.int64), denominator

    def joint(self, time, alpha, inverse_binding):
        if time < 0:
            raise ValueError('time must be nonnegative')
        states, full = self.full_generator(alpha, inverse_binding)
        launch = states.index((0, 0, self.fuel))
        indices = [launch] + [i for i, state in enumerate(states) if state[0] > 0]
        matrix = full[np.ix_(indices, indices)]
        matrix[1:, 0] = 0
        return float(expm(float(time) * matrix)[0, 1:].sum())

    def formation_and_maintenance(self, deadline, lifetime, alpha, inverse_binding):
        if deadline < 0 or lifetime < 0:
            raise ValueError('durations must be nonnegative')
        if not self.pairs or not self.fuel:
            return 0.0
        states, full = self.full_generator(alpha, inverse_binding)
        indices = [i for i, state in enumerate(states) if state[0] > 0]
        live = [states[i] for i in indices]
        initial = live.index((1, 0, self.fuel - 1))
        survival = expm(float(lifetime) * full[np.ix_(indices, indices)])[initial].sum()
        arrival = float(Fraction(alpha) * self.assembly_rate * self.pairs ** 2 * self.fuel)
        return float(-np.expm1(-arrival * deadline) * survival)


def joint_interval(model, time, alpha, inverse_binding, scale=2 ** 27):
    """Rational uniformization enclosure with terminal reward 0 at launch."""
    time = Fraction(time)
    if time < 0:
        raise ValueError('time must be nonnegative')
    _, matrix, denominator = model.killed_generator(alpha, inverse_binding)
    exit_rate = int(np.max(-matrix.diagonal()))
    if exit_rate == 0:
        return Fraction(0), Fraction(0)
    transition = matrix + exit_rate * np.eye(len(matrix), dtype=np.int64)
    assert np.all(transition >= 0)
    if (scale + 1) * exit_rate > np.iinfo(np.int64).max:
        raise ValueError('uniformization multiplication could overflow')
    weights_lo, weights_hi, tail = _poisson_weights(time * exit_rate / denominator, scale)
    if (len(weights_hi) + 2) * scale ** 2 > np.iinfo(np.int64).max:
        raise ValueError('uniformization accumulation could overflow')
    lower = np.full(len(matrix), scale, dtype=np.int64)
    upper = lower.copy()
    lower[0] = upper[0] = 0
    total_lower, total_upper = 0, tail * scale
    for k, (weight_lo, weight_hi) in enumerate(zip(weights_lo, weights_hi)):
        if k:
            lower = transition @ lower // exit_rate
            upper = (transition @ upper + exit_rate - 1) // exit_rate
        total_lower += weight_lo * int(lower[0])
        total_upper += weight_hi * int(upper[0])
    return Fraction(total_lower, scale ** 2), min(Fraction(1), Fraction(total_upper, scale ** 2))


def maintenance_interval(model, time, alpha, inverse_binding, scale=2 ** 27):
    """Use the existing nonnegative survival certifier on post-launch states."""
    from .autonomous import survival_intervals
    states, matrix, denominator = model.killed_generator(alpha, inverse_binding)
    if not model.pairs or not model.fuel:
        return Fraction(0), Fraction(0)
    initial = states.index((1, 0, model.fuel - 1)) - 1
    lower, upper, divisor = survival_intervals(
        matrix[1:, 1:][None], denominator, time, np.array([initial]), scale=scale)
    return Fraction(int(lower[0]), divisor), Fraction(int(upper[0]), divisor)


def success_interval(model, deadline, lifetime, alpha, inverse_binding, scale=2 ** 27):
    deadline = Fraction(deadline)
    if deadline < 0:
        raise ValueError('deadline must be nonnegative')
    arrival = Fraction(alpha) * model.assembly_rate * model.pairs ** 2 * model.fuel
    if not arrival:
        return Fraction(0), Fraction(0)
    lower, upper, _ = _poisson_weights(arrival * deadline, scale)
    formation_lo = 1 - Fraction(upper[0], scale)
    formation_hi = 1 - Fraction(lower[0], scale)
    survival_lo, survival_hi = maintenance_interval(model, lifetime, alpha, inverse_binding, scale)
    return formation_lo * survival_lo, formation_hi * survival_hi


@lru_cache(maxsize=16)
def _curvature_data(model):
    """Cache parameter-independent rational coefficients for rectangle bounds."""
    launch = model.assembly_rate * model.pairs ** 2 * model.fuel
    matrices = {}
    for a, z in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        _, integers, denominator = model.killed_generator(a, z)
        matrices[a, z] = [[Fraction(int(x), denominator) for x in row[1:]] for row in integers[1:]]
    def norm_difference(left, right):
        return max(sum(abs(x-y) for x, y in zip(row1, row2))
                   for row1, row2 in zip(left, right))
    na = max(norm_difference(matrices[1, z], matrices[0, z]) for z in [0, 1])
    nz = max(norm_difference(matrices[a, 1], matrices[a, 0]) for a in [0, 1])
    # Labelled-channel intensities, retaining rates that lead into failure.
    totals = {}
    for a, z in [(1, 0), (0, 1), (1, 1)]:
        by_state = {state: {'assembly': Fraction(0), 'repair': Fraction(0), 'z': Fraction(0)}
                    for state in model.states() if state[0] > 0}
        for name, s, t, forward, reverse in model.edges(a, z):
            for start, rate, backward in [(s, forward, False), (t, reverse, True)]:
                if start not in by_state or name == 'damage':
                    continue
                by_state[start][name] += rate
                if name == 'repair' or backward:
                    by_state[start]['z'] += rate
        totals[a, z] = by_state
    return launch, na, nz, totals


def curvature_bounds(model, deadline, lifetime, box):
    """Exact Duhamel and labelled-path Fisher pure-derivative bounds.

    For affine labelled intensities the second path likelihood derivative
    is score squared minus observed jump information. Its positive and
    negative parts each have expectation at most the Fisher information.
    """
    a0, a1, z0, z1 = map(Fraction, box)
    deadline, lifetime = Fraction(deadline), Fraction(lifetime)
    launch, na, nz, totals = _curvature_data(model)
    bound_a = (deadline * 2 * launch + lifetime * na) ** 2
    bound_z = (lifetime * nz) ** 2
    if a0 > 0 and a1 < 1:
        fisher = max(totals[1, 1][state]['assembly'] / a0
                     + totals[0, 1][state]['repair'] / (1-a1) for state in totals[1, 1])
        bound_a = min(bound_a, deadline * launch / a0 + lifetime * fisher)
    if z0 > 0:
        fisher = max(max(totals[1, 1][state]['z'], totals[0, 1][state]['z'])
                     for state in totals[1, 1]) / z0
        bound_z = min(bound_z, lifetime * fisher)
    return bound_a, bound_z


def one_fuel_joint(model, time, alpha, inverse_binding):
    """Closed form with one fuel molecule, including the equal-rate limit."""
    if model.fuel != 1:
        raise ValueError('closed form requires one initial fuel molecule')
    rate = float(model.assembly_rate * Fraction(alpha))
    arrival = rate * model.pairs ** 2
    loss = float(model.damage_rate + rate * Fraction(inverse_binding) / model.fuel_factor)
    if abs(arrival - loss) < 1e-12:
        return arrival * time * np.exp(-arrival * time)
    return arrival * (np.exp(-loss * time) - np.exp(-arrival * time)) / (arrival - loss)
