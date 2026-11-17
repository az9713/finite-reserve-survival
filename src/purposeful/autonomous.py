"""Autonomous shared-reserve reaction network and checked survival enclosures."""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import lcm

import numpy as np
from scipy.linalg import eig
from scipy.stats import poisson


State = tuple[int, int, int, int, int]  # battery, environment, sensor, actuator, damage
Gate = tuple[int, ...]


@dataclass(frozen=True)
class Kinetics:
    capacity: int
    token_factor: int
    damage_factor: int
    fuel_factor: int
    environment: Fraction
    sensing: Fraction
    control: Fraction
    harvest: Fraction
    wear: Fraction
    failure: Fraction
    repair: Fraction
    both_matches: bool = True

    def __post_init__(self):
        if not isinstance(self.capacity, int) or self.capacity < 1:
            raise ValueError("capacity must be a positive integer")
        factors = (self.damage_factor, self.token_factor, self.fuel_factor)
        if not all(isinstance(factor, int) for factor in factors) or not 1 <= factors[0] < factors[1] < factors[2]:
            raise ValueError("integer factors must satisfy 1 <= damage < token < fuel")
        for name in ("environment", "sensing", "control", "harvest", "wear", "failure", "repair"):
            rate = Fraction(getattr(self, name))
            if rate <= 0:
                raise ValueError("all reference rates must be positive")
            object.__setattr__(self, name, rate)

    def energy_factor(self, state):
        battery, _, _, _, damage = state
        return Fraction(self.token_factor ** battery, self.damage_factor ** damage)


@dataclass(frozen=True)
class Reaction:
    name: str
    source: State
    target: State
    forward: Fraction
    reverse: Fraction
    work_factor: int
    gate: Gate


@dataclass(frozen=True)
class Circuit:
    """Sensing activity at each funded level and actuator targets at each (b,y)."""

    sensing: tuple[int, ...]
    targets: tuple[int, ...]  # -1 means off. Otherwise: 0 or 1.


def benchmark_kinetics():
    return Kinetics(
        capacity=2, token_factor=8, damage_factor=2, fuel_factor=128,
        environment=Fraction(1, 2), sensing=11, control=Fraction(21, 2),
        harvest=Fraction(17, 4), wear=Fraction(5, 2), failure=Fraction(5, 2),
        repair=Fraction(21, 4), both_matches=False,
    )


def reactions(kinetics, informative=True):
    """Yield physical reversible channels, including return from the failed state."""
    C = kinetics.capacity
    for b, x, y, a, d in product(range(C + 1), range(2), range(2), range(2), range(3)):
        state = b, x, y, a, d
        channels = []
        if x == 0:
            channels.append(("environment", (b, 1, y, a, d), kinetics.environment, 1, ()))
        if b and (y != x or not informative):
            channels.append(("sense", (b - 1, x, 1 - y, a, d), kinetics.sensing, 1, (b,)))
        if b:
            channels.append(("control", (b - 1, x, y, 1 - a, d), kinetics.control, 1, (b, y, 1 - a)))
        if b < C and x == a and (x or kinetics.both_matches):
            channels.append(("harvest", (b + 1, x, y, a, d), kinetics.harvest, kinetics.fuel_factor, ()))
        if x != a and d < 2:
            name, rate = ("wear", kinetics.wear) if d == 0 else ("failure", kinetics.failure)
            channels.append((name, (b, x, y, a, d + 1), rate, 1, ()))
        if b and d == 1:
            channels.append(("repair", (b - 1, x, y, a, 0), kinetics.repair, 1, ()))
        for name, target, forward, work, gate in channels:
            ratio = work * kinetics.energy_factor(state) / kinetics.energy_factor(target)
            yield Reaction(name, state, target, forward, forward / ratio, work, gate)


class Network:
    """Compile common integer generators for the family of physical gate settings."""

    def __init__(self, kinetics, informative=True):
        self.kinetics = kinetics
        self.states = list(product(range(kinetics.capacity + 1), range(2), range(2), range(2), range(2)))
        self.index = {state: i for i, state in enumerate(self.states)}
        channels = list(reactions(kinetics, informative))
        self.denominator = lcm(*(rate.denominator for edge in channels for rate in (edge.forward, edge.reverse)))
        total = sum(int(rate * self.denominator) for edge in channels for rate in (edge.forward, edge.reverse))
        if total > np.iinfo(np.int64).max:
            raise ValueError("reference rates cannot be represented in the integer generator")
        size = len(self.states)
        self.components = {
            gate: np.zeros((size, size), dtype=np.int64)
            for gate in [(), *((b,) for b in range(1, kinetics.capacity + 1)),
                         *product(range(1, kinetics.capacity + 1), range(2), range(2))]
        }
        for edge in channels:
            part = self.components[edge.gate]
            for source, target, rate in [(edge.source, edge.target, edge.forward),
                                         (edge.target, edge.source, edge.reverse)]:
                i, j = self.index.get(source), self.index.get(target)
                if i is not None:
                    integer_rate = int(rate * self.denominator)
                    part[i, i] -= integer_rate
                    if j is not None:
                        part[i, j] += integer_rate

    def initial(self, battery):
        """Unknown stationary environment, with sensor/actuator at zero and no damage."""
        if battery not in range(self.kinetics.capacity + 1):
            raise ValueError("initial reserve must lie within capacity")
        return np.array([self.index[battery, x, 0, 0, 0] for x in range(2)])

    def control_bank(self):
        """All digital control truth tables, each matrix including ungated dynamics."""
        keys = list(product(range(1, self.kinetics.capacity + 1), range(2)))
        targets = list(product((-1, 0, 1), repeat=len(keys)))
        return targets, np.stack([
            sum((self.components[b, y, t] for (b, y), t in zip(keys, table) if t != -1),
                self.components[()].copy())
            for table in targets
        ])

    def digital_circuits(self):
        """All binary sensing gates paired with all permitted control truth tables."""
        targets, control = self.control_bank()
        masks = list(product((0, 1), repeat=self.kinetics.capacity))
        sense = np.stack([
            sum((self.components[b,] for b, flag in enumerate(mask, start=1) if flag),
                np.zeros_like(self.components[()]))
            for mask in masks
        ])
        return ([Circuit(mask, table) for mask, table in product(masks, targets)],
                (sense[:, None] + control).reshape(-1, *control.shape[-2:]))


def _uniformization(matrices, exit_rate=None):
    if exit_rate is None:
        exit_rate = int(np.max(-np.diagonal(matrices, axis1=-2, axis2=-1)))
    scaled = matrices + exit_rate * np.eye(matrices.shape[-1], dtype=np.int64)
    if np.any(scaled < 0) or np.any(scaled.sum(axis=-1) > exit_rate):
        raise ValueError("expected integer subgenerators")
    return scaled, exit_rate


def survival_curves(matrices, denominator, times, initial):
    """Floating-point curves for illustration. Use survival_intervals for proofs."""
    scaled, exit_rate = _uniformization(matrices)
    means = np.asarray(times) * (exit_rate / denominator)
    order = int(poisson.isf(1e-14, means.max()))
    probabilities = poisson.pmf(np.arange(order + 1)[:, None], means[None, :])
    transition = scaled / exit_rate
    value = np.ones(matrices.shape[:-1])
    powers = [value[..., initial].mean(axis=-1)]
    for _ in range(order):
        value = np.einsum("...ij,...j->...i", transition, value)
        powers.append(value[..., initial].mean(axis=-1))
    return np.stack(powers, axis=-1) @ probabilities


@lru_cache(maxsize=16)
def _poisson_weights(mean, scale):
    """Outward-rounded normalized weights and an upper tail bound."""
    raw, total = [Fraction(1)], Fraction(1)
    while True:
        next_weight = raw[-1] * mean / len(raw)
        ratio = mean / (len(raw) + 1)
        if ratio < 1:
            tail = next_weight / (1 - ratio)
            if tail * scale * 4 <= total:
                break
        raw.append(next_weight)
        total += next_weight
    lower = [int(weight * scale // (total + tail)) for weight in raw]
    upper = [_ceil(weight * scale / total) for weight in raw]
    return lower, upper, _ceil(tail * scale / total)


def _ceil(value):
    return -(-value.numerator // value.denominator)


def survival_intervals(matrices, denominator, time, initial, scale=2 ** 29, exit_rate=None):
    """Return certified lower/upper integer numerators with denominator scale squared."""
    transition, exit_rate = _uniformization(matrices, exit_rate)
    time = Fraction(time)
    if time < 0:
        raise ValueError("time must be nonnegative")
    w_lo, w_hi, tail = _poisson_weights(time * exit_rate / denominator, scale)
    bound = np.iinfo(np.int64).max
    if (scale + 1) * exit_rate > bound or (2 * scale + len(w_hi) + len(initial)) * scale > bound:
        raise ValueError("certification scale would overflow int64")
    lo = np.full(matrices.shape[:-1], scale, dtype=np.int64)
    hi = lo.copy()
    result_lo = np.zeros(matrices.shape[:-2], dtype=np.int64)
    result_hi = np.full(matrices.shape[:-2], tail * scale, dtype=np.int64)
    for k, (lower, upper) in enumerate(zip(w_lo, w_hi)):
        if k:
            lo = np.einsum("...ij,...j->...i", transition, lo) // exit_rate
            hi = (np.einsum("...ij,...j->...i", transition, hi) + exit_rate - 1) // exit_rate
        start_lo = lo[..., initial].sum(axis=-1) // len(initial)
        start_hi = (hi[..., initial].sum(axis=-1) + len(initial) - 1) // len(initial)
        result_lo += lower * start_lo
        result_hi += upper * start_hi
    return result_lo, np.minimum(result_hi, scale ** 2), scale ** 2


def uniform_activity_enclosure(network, time, initial, subdivisions=256):
    """Bound the best uniform sensing intensity in [0,1] with digital control gates."""
    if not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("subdivisions must be a positive integer")
    time = Fraction(time)
    _, base = network.control_bank()
    sensor = sum(network.components[b,] for b in range(1, network.kinetics.capacity + 1))
    max_exit = int(np.max(-np.diagonal(base + sensor, axis1=-2, axis2=-1)))
    norm = Fraction(int(np.max(np.abs(sensor).sum(axis=-1))), network.denominator)
    denominator = subdivisions * network.denominator
    if subdivisions * max_exit > np.iinfo(np.int64).max:
        raise ValueError("uniform activity grid would overflow the integer generator")
    fixed = subdivisions * base
    lower = upper = Fraction(0)
    for node in range(subdivisions + 1):
        lo, hi, scale = survival_intervals(
            fixed + node * sensor, denominator, time, initial, exit_rate=subdivisions * max_exit
        )
        lower = max(lower, Fraction(int(lo.max()), scale))
        upper = max(upper, Fraction(int(hi.max()), scale))
    interpolation = (time * norm) ** 2 / (8 * subdivisions ** 2)
    return lower, min(Fraction(1), upper + interpolation)


def spectral_interval(matrix, denominator, initial):
    """Certify the reachable spectral abscissa by an exact positive-vector inequality."""
    reachable = set(map(int, initial))
    edges = matrix > 0
    frontier = list(reachable)
    while frontier:
        neighbors = set(np.flatnonzero(edges[frontier.pop()])) - reachable
        reachable.update(neighbors)
        frontier.extend(neighbors)
    indices = sorted(reachable)
    sub = matrix[np.ix_(indices, indices)]
    eigenvalues, vectors = eig(sub / denominator)
    vector = np.abs(vectors[:, np.argmax(eigenvalues.real)].real)
    witness = [max(1, int(round(value * 2 ** 38 / vector.max()))) for value in vector]
    ratios = [
        Fraction(sum(int(coefficient) * x for coefficient, x in zip(row, witness)), denominator * x0)
        for row, x0 in zip(sub, witness)
    ]
    return min(ratios), max(ratios)
