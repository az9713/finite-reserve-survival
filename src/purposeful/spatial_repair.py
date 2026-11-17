"""Reversible finite-fuel maintenance for two Brownian complexes on a cubic torus.

The relative periodic grid approximates Brownian motion. Reflections, axis
permutations, and interchange of the complex labels are exact symmetries of
that grid process. Rows are backward generators and failure means first h=0.
"""
from dataclasses import dataclass
from itertools import combinations_with_replacement
from numbers import Integral

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.special import beta, comb, ive
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class SpatialRepair:
    fuel: int = 4
    eta: float = 1.
    diffusivity: float = .025
    length: float = 1.
    radius: float = .32
    activity: float = 16.
    local_width: float = .075

    def __post_init__(self):
        if isinstance(self.fuel, bool) or not isinstance(self.fuel, Integral) or self.fuel < 0:
            raise ValueError('fuel must be a nonnegative integer')
        for key in ('eta', 'diffusivity', 'length', 'radius', 'activity', 'local_width'):
            value = getattr(self, key)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(key + ' must be finite and positive')
        if self.radius > self.length / 2:
            raise ValueError('contact radius exceeds half the box length')


def gate_activity(model, r):
    """Dimensionless compact radial gate; intrinsic rate is model.activity."""
    return np.maximum(0., 1. - (np.asarray(r) / model.radius) ** 2) ** 4


def continuum_gate_mean(model):
    return 2 * np.pi * (model.radius / model.length) ** 3 * beta(1.5, 5.)


@dataclass(frozen=True)
class SpatialGrid:
    n: int
    points: np.ndarray
    r: np.ndarray
    weight: np.ndarray
    gate: np.ndarray
    laplacian: sparse.csr_matrix


def build_grid(model, n):
    """Lump the n³ periodic six-neighbor grid by its 48 cube symmetries."""
    if not isinstance(n, Integral) or isinstance(n, bool) or n < 2 or n % 2:
        raise ValueError('n must be an even integer of at least two')
    half = n // 2
    representatives = list(combinations_with_replacement(range(half + 1), 3))
    index = {point: i for i, point in enumerate(representatives)}
    points = np.asarray(representatives, dtype=int)
    rows, cols, multiplicities = [], [], []
    for i, point in enumerate(representatives):
        distinct = len(set(point))
        permutations = {1: 1, 2: 3, 3: 6}[distinct]
        signs = 2 ** sum(0 < x < half for x in point)
        multiplicities.append(permutations * signs)
        for axis in range(3):
            for step in (-1, 1):
                neighbor = list(point)
                moved = (neighbor[axis] + step) % n
                neighbor[axis] = min(moved, n - moved)
                rows.append(i)
                cols.append(index[tuple(sorted(neighbor))])
    spacing = model.length / n
    scale = model.diffusivity / spacing ** 2
    size = len(representatives)
    jumps = sparse.coo_matrix((np.full(len(rows), scale), (rows, cols)), shape=(size, size)).tocsr()
    diffusion = jumps - (6 * scale) * sparse.eye(size, format='csr')
    radius = spacing * np.linalg.norm(points, axis=1)
    weights = np.asarray(multiplicities, dtype=float) / n ** 3
    return SpatialGrid(n, points, radius, weights, gate_activity(model, radius), diffusion)


def preparation_weights(model, grid):
    """Uniform, centered truncated Gaussian, and equilibrium repair-event Palm laws."""
    localized = grid.weight * np.exp(-.5 * (grid.r / model.local_width) ** 2)
    reaction = grid.weight * grid.gate
    return {'uniform': grid.weight, 'localized': localized / localized.sum(),
            'reaction': reaction / reaction.sum()}


def chemical_generators(model, rates=None):
    """Damage and contact row generators after exact label exchange lumping.

    rates[b] multiplies BOTH directions of the edge HD,b <-> HH,b-1.
    HH has two directed reverse channels, one for each possible helper.
    """
    if rates is None:
        rates = np.full(model.fuel + 1, model.activity)
    rates = np.asarray(rates, dtype=float)
    if rates.shape != (model.fuel + 1,) or np.any(~np.isfinite(rates[1:])) or np.any(rates[1:] < 0):
        raise ValueError('supply one finite nonnegative edge rate per positive fuel count')
    states = [(h, b) for h in range(3) for b in range(model.fuel + 1)]
    index = {state: i for i, state in enumerate(states)}
    damage = np.zeros((len(states), len(states)))
    contact = np.zeros_like(damage)
    for i, (h, b) in enumerate(states):
        if h:
            damage[i, index[h - 1, b]] = h * model.eta
        if h < 2:
            damage[i, index[h + 1, b]] = (2 - h) * model.eta / 2
        if h == 1 and b:
            contact[i, index[2, b - 1]] = b * rates[b]
        if h == 2 and b < model.fuel:
            contact[i, index[1, b + 1]] = (model.fuel - b) * rates[b + 1] / 2
    for matrix in (damage, contact):
        np.fill_diagonal(matrix, -matrix.sum(axis=1))
    return states, damage, contact


def equilibrium(model):
    """Count weights for the uncensored chain; relative position is independent uniform."""
    states, _, _ = chemical_generators(model)
    weights = np.asarray([comb(2, h) * 2. ** (2 - h) * comb(model.fuel, b) * (1 / 8) ** b
                          for h, b in states])
    return weights / weights.sum()


def calibrate(model, grid, diagnostics=None):
    """Match each forward edge's uniform-start isolated mean encounter time.

    (-D Delta + b*kappa*g)T_b = 1; k_b = 1/(b E_uniform[T_b]). The reverse
    edge shares k_b to conserve the equilibrium ratio. b=0 has no edge.
    """
    coefficients = np.full(model.fuel + 1, np.nan)
    ones = np.ones(len(grid.r))
    residuals, means = [], []
    for b in range(1, model.fuel + 1):
        operator = sparse.diags(b * model.activity * grid.gate) - grid.laplacian
        times = spsolve(operator.tocsc(), ones)
        mean = np.dot(grid.weight, times)
        coefficients[b] = 1 / (b * mean)
        if diagnostics is not None:
            residuals.append(float(np.max(np.abs(operator @ times - ones))))
            means.append(float(mean))
    if diagnostics is not None:
        diagnostics.update(fuel_before_forward=list(range(1, model.fuel + 1)),
                           uniform_mean_reaction_time=means, infinity_residual=residuals)
    return coefficients


def spatial_generator(model, grid, killed=True):
    """Return (chemical states, CSR backward generator), spatial index fastest."""
    states, damage, contact = chemical_generators(model)
    if killed:
        live = [i for i, (h, _) in enumerate(states) if h]
        states = [states[i] for i in live]
        damage = damage[np.ix_(live, live)]
        contact = contact[np.ix_(live, live)]
    space_id = sparse.eye(len(grid.r))
    chemistry_id = sparse.eye(len(states))
    generator = (sparse.kron(chemistry_id, grid.laplacian, format='csr')
                 + sparse.kron(damage, space_id, format='csr')
                 + sparse.kron(contact, sparse.diags(grid.gate), format='csr'))
    return states, generator


def _times_array(times):
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or not len(times) or np.any(~np.isfinite(times)) or np.any(times < 0):
        raise ValueError('times must be a nonempty vector of finite nonnegative values')
    return times


def _initial_index(states, initial):
    try:
        return states.index(tuple(initial))
    except (TypeError, ValueError) as exc:
        raise ValueError('initial must be a live (h, b) state') from exc


def count_survival(model, times, rates, initial=None):
    times = _times_array(times)
    states, damage, contact = chemical_generators(model, rates)
    live = [i for i, (h, _) in enumerate(states) if h]
    i = _initial_index([states[k] for k in live], (1, model.fuel) if initial is None else initial)
    killed = (damage + contact)[np.ix_(live, live)]
    return np.asarray([expm(t * killed)[i].sum() for t in times])


def _chebyshev_tail_log(s, first_omitted):
    if s == 0:
        return -np.inf
    k = first_omitted
    theta = np.arcsinh(k / s)
    return np.log(2.) + k * (k / (np.hypot(k, s) + s) - theta)


def _chebyshev_degree(s, target):
    """Choose K such that the discarded scalar Bessel tail is below target.

    e^-s I_|k|(s) is the probability mass of a symmetric Skellam variate.
    Chernoff's bound optimized over its exponential moment gives the tail.
    """
    if s == 0:
        return 0
    target_log = np.log(target)
    k = 1
    while _chebyshev_tail_log(s, k) > target_log:
        k += 1
    return k - 1


def spatial_survivals(model, grid, times, initials, diagnostics=None):
    """Apply the finite-grid killed semigroup at all requested times and starts.

    Pointwise detailed balance gives a self-adjoint generator in the equilibrium
    inner product. Gershgorin puts its spectrum in [-2q,0], q=max exit rate.
    With B=I+Q/q and s=qt, exp(tQ)=exp(-s) sum (2-delta_k0) I_k(s) T_k(B).
    Only the projections of T_k(B)1 on the requested starts are retained.
    """
    times = _times_array(times)
    states, generator = spatial_generator(model, grid)
    initials = [tuple(initial) for initial in initials]
    if not initials:
        raise ValueError('at least one initial chemical state is required')
    indices = [_initial_index(states, initial) for initial in initials]
    laws = preparation_weights(model, grid)
    weights = np.asarray(list(laws.values()))
    q = max(-generator.diagonal())
    scaled_times = q * times
    # Cauchy-Schwarz in L²(pi): ||density of start|| * ||live payoff||.
    # The chemical probability is the full equilibrium mass on each state.
    live_pi = equilibrium(model)[model.fuel + 1:]
    spatial_norm_squared = (weights ** 2 / grid.weight).sum(axis=1)
    amplification = np.sqrt(live_pi.sum() * spatial_norm_squared.max() / live_pi[indices].min())
    target = min(2e-16, 2e-14 / amplification)
    maximum_time = scaled_times.max()
    degree = _chebyshev_degree(maximum_time, target)
    if diagnostics is not None:
        tail = np.exp(_chebyshev_tail_log(maximum_time, degree + 1))
        diagnostics.update(maximum_exit_rate=float(q), chebyshev_degree=int(degree),
                           scalar_tail_upper=float(tail), projection_norm_factor=float(amplification),
                           exact_arithmetic_observable_tail_upper=float(tail * amplification))
    coefficients = ive(np.arange(degree + 1)[:, None], scaled_times[None, :])
    coefficients[1:] *= 2
    scaled = sparse.eye(generator.shape[0], format='csr') + generator / q
    size = len(grid.r)
    projections = np.empty((degree + 1, len(initials), len(laws)))

    def project(k, vector):
        selected = vector.reshape(len(states), size)[indices]
        projections[k] = np.einsum('ij,pj->ip', selected, weights)

    previous = np.ones(generator.shape[0])
    project(0, previous)
    if degree:
        current = scaled @ previous
        project(1, current)
        for k in range(2, degree + 1):
            following = 2 * (scaled @ current) - previous
            project(k, following)
            previous, current = current, following
    values = np.einsum('kt,kip->ipt', coefficients, projections)
    return {initial: {name: values[i, p] for p, name in enumerate(laws)}
            for i, initial in enumerate(initials)}


def spatial_survival(model, grid, times, initial=None):
    initial = (1, model.fuel) if initial is None else tuple(initial)
    return spatial_survivals(model, grid, times, [initial])[initial]

