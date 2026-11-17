"""Closed finite-fuel spatial maintenance with a reversible contact-formed tether."""
from dataclasses import dataclass
from numbers import Real

import numpy as np
from scipy import sparse
from scipy.integrate import quad_vec
from scipy.special import ive

from .spatial_repair import (
    SpatialGrid, SpatialRepair, _chebyshev_degree, _chebyshev_tail_log,
    _times_array, build_grid, chemical_generators, equilibrium as chemical_equilibrium,
)


@dataclass(frozen=True)
class Organization(SpatialRepair):
    passive: float = 1.
    paid: float = 4.
    offset: float = 2.
    depth: float = 6.

    def __post_init__(self):
        super().__post_init__()
        if self.fuel < 2:
            raise ValueError('spectator comparison needs at least two fuel molecules')
        if not np.isfinite(self.offset):
            raise ValueError('tether offset must be finite')
        for name in ('depth', 'passive', 'paid'):
            value = getattr(self, name)
            if not np.isfinite(value) or value < 0:
                raise ValueError(name + ' must be finite and nonnegative')


def _tether_diffusion(grid, potential):
    """Tilt each directed nearest-neighbor rate by the half energy difference."""
    jumps = grid.laplacian.tocoo()
    off = jumps.row != jumps.col
    source, target = jumps.row[off], jumps.col[off]
    values = jumps.data[off] * np.exp((potential[source] - potential[target]) / 2)
    rates = sparse.coo_matrix((values, (source, target)), shape=jumps.shape).tocsr()
    return rates - sparse.diags(np.asarray(rates.sum(axis=1)).ravel())


def _conductance(policy, b, inventory):
    if isinstance(policy, str) and policy == 'spectator':
        return (b - 1) / (inventory - 1)
    if isinstance(policy, bool) or not isinstance(policy, Real) or not 0 <= policy <= 1:
        raise ValueError("policy must be a constant in [0,1] or 'spectator'")
    return float(policy)


@dataclass
class OrganizationSystem:
    model: Organization
    grid: SpatialGrid
    states: tuple
    base: sparse.csr_matrix
    paid: sparse.csr_matrix
    spectator: sparse.csr_matrix
    potential: np.ndarray
    killed: bool

    def generator(self, policy):
        if isinstance(policy, str) and policy == 'spectator':
            return self.base + self.spectator
        alpha = _conductance(policy, 1, self.model.fuel)
        return self.base + alpha * self.paid

    def initial(self, h=1, b=None, q=0):
        """Common uniform relative-position law in the selected chemical state."""
        state = (h, self.model.fuel if b is None else b, q)
        if state not in self.states:
            raise ValueError('initial state is outside the represented state space')
        out = np.zeros((len(self.states), len(self.grid.r)))
        out[self.states.index(state)] = self.grid.weight
        return out.ravel()

    def equilibrium(self):
        """Full equilibrium, or its normalized live restriction for Hilbert norms."""
        chemical = chemical_equilibrium(self.model).reshape(3, self.model.fuel + 1)
        spatial = [self.grid.weight, self.grid.weight * np.exp(-self.potential)]
        weights = np.concatenate([chemical[h, b] * spatial[q] for h, b, q in self.states])
        return weights / weights.sum()

    def _series(self, policy, times, payoff, initial, cumulative, diagnostics):
        """Project one forward Chebyshev recurrence onto all endpoint/reward columns."""
        times = _times_array(times)
        generator = self.generator(policy)
        rate = float(np.max(-generator.diagonal()))
        pi = self.equilibrium()
        norm = np.sqrt(np.sum(initial ** 2 / pi))
        factor = norm * np.max(np.sqrt(np.einsum('rs,rs,s->r', payoff, payoff, pi)))
        maximum_time = float(max(times))
        scaled_end = rate * maximum_time
        # Integrated observables acquire a further factor of at most T.
        scale = max(1., factor * max(1., maximum_time))
        degree = _chebyshev_degree(scaled_end, min(2e-16, 2e-14 / scale))
        orders = np.arange(degree + 1)
        coeff = ive(orders[:, None], rate * times[None, :])
        coeff[1:] *= 2
        basis = (sparse.eye(generator.shape[0], format='csr') + generator.T.tocsr() / rate)
        proj = np.empty((degree + 1, len(payoff)))

        def project(order, vector):
            proj[order] = np.einsum('rs,s->r', payoff, vector)

        prev = initial
        project(0, prev)
        if degree:
            cur = basis @ prev
            project(1, cur)
            for k in range(2, degree + 1):
                nxt = 2 * (basis @ cur) - prev
                project(k, nxt)
                prev, cur = cur, nxt
        values = np.einsum('kr,kt->rt', proj, coeff)
        if cumulative:
            # Quadrature is one-dimensional and independent of the spatial grid.
            for i, t in enumerate(times):
                if t == 0:
                    values[-cumulative:, i] = 0.
                else:
                    integ = quad_vec(lambda u: ive(orders, rate * u), 0, t,
                                     epsabs=1e-13, epsrel=1e-13)[0]
                    integ[1:] *= 2
                    values[-cumulative:, i] = proj[:, -cumulative:].T @ integ
        if diagnostics is not None:
            tail = float(np.exp(_chebyshev_tail_log(scaled_end, degree + 1)))
            diagnostics.update(maximum_exit_rate=rate, chebyshev_degree=int(degree),
                               scalar_tail_upper=tail, projection_norm_factor=float(factor),
                               exact_arithmetic_truncation_upper=tail * scale)
        return values

    def survival(self, policy, times, diagnostics=None):
        """Live probability from uniform h=1,b=N,q=0. Requires killed dynamics."""
        if not self.killed:
            raise ValueError('survival requires the killed generator')
        initial = self.initial()
        return self._series(policy, times, np.ones((1, len(initial))), initial,
                            0, diagnostics)[0]

    def observables(self, policy, times):
        """Unconditional endpoint masses and gross event expectations per launch.

        A killed trajectory has zero live endpoint mass and no further events.
        Stopped fuel/tether moments retain the values on entry to the failure set.
        For the full process, absorption rewards vanish and mass is always one.
        """
        model, grid = self.model, self.grid
        names = ('mass', 'tether', 'contact', 'gate', 'fuel',
                 'paid_forward', 'paid_reverse', 'repair_forward', 'repair_reverse',
                 'passive_forward', 'passive_reverse', 'absorbed_fuel', 'absorbed_tether')
        payoffs = np.zeros((len(names), len(self.states), len(grid.r)))
        gate_on = grid.gate * np.exp(-np.maximum(self.potential, 0.))
        gate_off = grid.gate * np.exp(np.minimum(self.potential, 0.))
        for i, (h, b, q) in enumerate(self.states):
            payoffs[0, i] = 1.
            payoffs[1, i] = q
            payoffs[2, i] = grid.r < model.radius
            payoffs[3, i] = grid.gate
            payoffs[4, i] = b
            if not q and b:
                payoffs[5, i] = model.paid * b * _conductance(policy, b, model.fuel) * gate_on
            if q and b < model.fuel:
                payoffs[6, i] = (model.paid * (model.fuel - b) / 8 *
                                  _conductance(policy, b + 1, model.fuel) * gate_off)
            if h == 1:
                payoffs[7, i] = model.activity * b * grid.gate
            if h == 2:
                payoffs[8, i] = model.activity * (model.fuel - b) / 2 * grid.gate
            payoffs[9, i] = model.passive * (1 - q) * gate_on
            payoffs[10, i] = model.passive * q * gate_off
            if self.killed and h == 1:
                # The sole killing jump damages the last healthy unit at rate eta.
                payoffs[11, i] = model.eta * b
                payoffs[12, i] = model.eta * q
        values = self._series(policy, times, payoffs.reshape(len(names), -1),
                              self.initial(), len(names) - 5, None)
        result = dict(zip(names, values))
        result['stopped_fuel'] = result['fuel'] + result['absorbed_fuel']
        result['stopped_tether'] = result['tether'] + result['absorbed_tether']
        return result


def build_system(model, n, killed=True):
    """Assemble the exactly cube-lumped finite grid; row/backward generators."""
    grid = build_grid(model, n)
    chem_states, damage, repair = chemical_generators(model)
    included = [i for i, (h, _) in enumerate(chem_states) if h or not killed]
    damage = damage[np.ix_(included, included)]
    repair = repair[np.ix_(included, included)]
    included_chemistry = [chem_states[i] for i in included]
    states = tuple((h, b, q) for h, b in included_chemistry for q in (0, 1))
    index = {state: i for i, state in enumerate(states)}
    potential = model.offset - model.depth * np.maximum(0, 1 - (grid.r / model.radius)**2)**2
    tether_diffusion = _tether_diffusion(grid, potential)
    size = len(grid.r)
    eye_q = sparse.eye(2, format='csr')
    base = (sparse.block_diag([tether_diffusion if q else grid.laplacian for _, _, q in states],
                              format='csr')
            + sparse.kron(sparse.kron(damage, eye_q), sparse.eye(size), format='csr')
            + sparse.kron(sparse.kron(repair, eye_q), sparse.diags(grid.gate), format='csr'))
    on = sparse.diags(grid.gate * np.exp(-np.maximum(potential, 0.)))
    off = sparse.diags(grid.gate * np.exp(np.minimum(potential, 0.)))
    # Each pair is chemical direction matrix x its spatial acceptance factor.
    def switch_pair(policy):
        forward = np.zeros((len(states), len(states)))
        reverse = np.zeros_like(forward)
        for i, (h, b, q) in enumerate(states):
            if q or (policy is not None and b == 0):
                continue
            if policy is None:
                j = index[h, b, 1]
                f = r = model.passive
            else:
                j = index[h, b - 1, 1]
                cap = model.paid * _conductance(policy, b, model.fuel)
                f = cap * b
                r = cap * (model.fuel - b + 1) / 8
            forward[i, j], forward[i, i] = f, -f
            reverse[j, i], reverse[j, j] = r, -r
        return sparse.kron(forward, on, format='csr') + sparse.kron(reverse, off, format='csr')

    return OrganizationSystem(model, grid, states, base + switch_pair(None),
                              switch_pair(1.), switch_pair('spectator'), potential, killed)
