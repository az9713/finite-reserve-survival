"""Exact finite-state models; energies are in k_B T, times in a chosen unit.

Sensor thermodynamics is explicit. Failure/repair layers are phenomenological;
their physical work costs are NOT included in the sensor power. No claim of a
complete autonomous organism or universal thermodynamic bound is made.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import expm
from scipy.special import expit, xlogy


@dataclass(frozen=True)
class Sensor:
    r: float  # symmetric external signal switching rate
    k: float  # sum of the two sensor transition rates (kinetic capacity)
    affinity: float  # energy gap / (k_B T), nonnegative

    def __post_init__(self):
        if not all(math.isfinite(x) for x in (self.r, self.k, self.affinity)):
            raise ValueError("parameters must be finite")
        if self.r <= 0 or self.k <= 0 or self.affinity < 0:
            raise ValueError("r,k > 0 and affinity >= 0 required")

    @property
    def a(self):
        """Mismatch -> match sensor transition rate."""
        return self.k * expit(self.affinity)

    @property
    def b(self):
        return self.k * expit(-self.affinity)

    @property
    def u(self):
        return self.r + self.b  # match -> mismatch, either channel

    @property
    def v(self):
        return self.r + self.a

    @property
    def s(self):
        return 2 * self.r + self.k

    @property
    def q(self):
        """Mismatch probability in the unfailed stationary ensemble."""
        return self.u / self.s

    @property
    def d(self):
        return self.k * math.tanh(self.affinity / 2) / self.s

    @property
    def power(self):
        """Net external work = bath heat at stationarity, in k_B T / time."""
        return self.affinity * self.r * self.d

    @property
    def mutual_information(self):
        q = self.q
        return float(math.log(2) + xlogy(q, q) + xlogy(1-q, 1-q))

    @property
    def learning_rate(self):
        return self.r * self.d * math.log((1-self.q) / self.q)

    def generator(self):
        """Row generator of four states (x,m) = 00,01,10,11.

        Keep environment and sensor channels separate for physical accounting.
        """
        g = np.zeros((4, 4))
        for x in range(2):
            for m in range(2):
                i = 2*x+m
                g[i, 2*(1-x)+m] = self.r
                g[i, 2*x+1-m] = self.b if x == m else self.a
                g[i, i] = -g[i].sum()
        return g

    def stationary(self):
        return np.array([(1-self.q)/2, self.q/2, self.q/2, (1-self.q)/2])


def sensor_for_accuracy(r, k, q):
    if not 0 < q <= .5:
        raise ValueError("q must be in (0, 1/2]")
    if r <= 0 or k <= 0:
        raise ValueError("r,k > 0 required")
    z = (1-2*q) * (1+2*r/k)
    if z >= 1:
        raise ValueError("accuracy infeasible at finite affinity for this capacity")
    return Sensor(r, k, 2*np.arctanh(z))


def frontier_power(r, k_max, q):
    """Exact minimum within this sensor family, at accuracy >= 1-q, k<=k_max."""
    return sensor_for_accuracy(r, k_max, q).power


def hazard_generator(sensor, h=1., baseline=0.):
    if h < 0 or baseline < 0:
        raise ValueError("hazards must be nonnegative")
    return np.array([[-sensor.u-baseline, sensor.u],
                     [sensor.v, -sensor.v-baseline-h]])


def survival_from_generator(g, initial, horizon):
    if horizon < 0:
        raise ValueError("horizon must be nonnegative")
    initial = np.asarray(initial)
    if len(initial) != len(g) or np.any(initial < 0) or not np.isclose(initial.sum(), 1):
        raise ValueError("initial must be a probability vector")
    return float(initial @ expm(g*horizon) @ np.ones(len(g)))


def hazard_survival(sensor, horizon, h=1., baseline=0., initial_q=None):
    q0 = sensor.q if initial_q is None else initial_q
    return survival_from_generator(hazard_generator(sensor, h, baseline),
                                   [1-q0, q0], horizon)


def hazard_decay(sensor, h=1., baseline=0.):
    """Asymptotic exponent, not equal to finite-time -log(S)/T in general."""
    if h < 0 or baseline < 0:
        raise ValueError("hazards must be nonnegative")
    disc = math.sqrt((sensor.s+h)**2-4*sensor.u*h)
    return baseline + (2*sensor.u*h)/(sensor.s+h+disc)


def damage_generator(sensor, stages=4, duration=2., baseline=0.):
    """Effective damage clock. Matching resets damage instantly.

    In mismatch, each of n stages advances at n/duration; last stage kills.
    Reset cost and molecular damage physics are not modeled.
    """
    if not isinstance(stages, int) or stages < 1 or duration <= 0 or baseline < 0:
        raise ValueError("positive integer stages, positive duration, baseline>=0")
    z = stages/duration
    g = np.zeros((stages+1, stages+1))
    g[0, 1] = sensor.u
    g[0, 0] = -sensor.u-baseline
    for j in range(1, stages+1):
        g[j, 0] = sensor.v
        if j < stages:
            g[j, j+1] = z
        g[j, j] = -sensor.v-z-baseline
    return g


def damage_survival(sensor, horizon, stages=4, duration=2., baseline=0., initial_q=None):
    q0 = sensor.q if initial_q is None else initial_q
    initial = np.zeros(stages+1)
    initial[0:2] = [1-q0, q0]  # fresh damage clock at t=0, documented preparation
    return survival_from_generator(damage_generator(sensor, stages, duration, baseline),
                                   initial, horizon)


def simulate_survival(g, initial, horizon, samples=10000, seed=0):
    """Independent Gillespie simulation of a killed continuous-time chain."""
    if horizon < 0 or samples < 1:
        raise ValueError("horizon>=0, samples>=1 required")
    rng = np.random.default_rng(seed)
    survivors = 0
    n = len(g)
    rates = -np.diag(g)
    channels = []
    for i in range(n):
        row = g[i].copy(); row[i] = 0
        death = max(0., -g[i].sum())
        probs = np.r_[row, death]
        channels.append(np.cumsum(probs/rates[i]) if rates[i] else None)
    for _ in range(samples):
        state = rng.choice(n, p=initial)
        time = 0.
        while True:
            rate = rates[state]
            if rate == 0:
                survivors += 1; break
            time += rng.exponential(1/rate)
            if time >= horizon:
                survivors += 1; break
            state = min(np.searchsorted(channels[state], rng.random(), side='right'), n)
            if state == n:
                break
    p = survivors/samples
    # Wilson interval is informative even for zero/all successes.
    z = 1.959963984540054
    center = (p+z*z/(2*samples))/(1+z*z/samples)
    half = z*math.sqrt(p*(1-p)/samples+z*z/(4*samples*samples))/(1+z*z/samples)
    return {"estimate": p, "ci95": [center-half, center+half], "samples": samples, "seed": seed}


def four_state_thermodynamics(sensor):
    """Independent flux accounting from the four-state stationary solve."""
    g = sensor.generator()
    matrix = g.T.copy(); matrix[-1] = 1
    rhs = np.array([0., 0., 0., 1.])
    p = np.linalg.solve(matrix, rhs)
    heat = work = sigma = learning = 0.
    for i in range(4):
        x, m = divmod(i, 2)
        ei = sensor.affinity * (x != m)
        for j in range(4):
            if i == j or g[i, j] == 0:
                continue
            xx, mm = divmod(j, 2)
            ej = sensor.affinity * (xx != mm)
            flow = p[i]*g[i, j]
            if x == xx:
                heat += flow*(ei-ej)
                learning += flow*math.log(p[j]/p[i])  # m marginal uniform
            else:
                work += flow*(ej-ei)
            sigma += flow*math.log(flow/(p[j]*g[j, i]))
    return {"stationary": p, "heat": heat, "work": work,
            "entropy_production": sigma, "learning_rate": learning}
