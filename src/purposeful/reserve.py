"""Checkpoint survival with finite reserves and reversible, fuel-paid repair."""
from dataclasses import dataclass
from math import exp, expm1, inf, isfinite, log
from numbers import Integral

import numpy as np
from scipy.special import expit


def _nonnegative_integer(value, name):
    if not isinstance(value, Integral) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")


def _probability(value, name):
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be a probability")


def _contact_fraction(rate_ceiling, duration):
    if not isfinite(rate_ceiling) or rate_ceiling <= 0:
        raise ValueError("rate_ceiling must be finite and positive")
    if not isfinite(duration) or duration <= 0:
        raise ValueError("duration must be finite and positive")
    return -expm1(-rate_ceiling * duration)


@dataclass(frozen=True)
class RepairContact:
    """Repair pair (damaged, b) <-> (healthy, b-1), with energies in k_B T."""

    token_energy: float
    rate_ceiling: float
    duration: float

    def __post_init__(self):
        if not isfinite(self.token_energy) or self.token_energy <= 0:
            raise ValueError("token_energy must be finite and positive")
        _contact_fraction(self.rate_ceiling, self.duration)

    def generator(self, rate=None):
        """Row generator in the order damaged, healthy, for a fueled pair."""
        rate = self.rate_ceiling if rate is None else rate
        if not 0 <= rate <= self.rate_ceiling:
            raise ValueError("rate must be between zero and rate_ceiling")
        forward = rate * expit(self.token_energy)
        reverse = rate * expit(-self.token_energy)
        return np.array([[-forward, forward], [reverse, -reverse]])

    @property
    def success(self):
        """Maximum probability of ending a contact healthy, starting damaged."""
        return float(expit(self.token_energy) * _contact_fraction(self.rate_ceiling, self.duration))


def stress_transition(q, mixing_rate, interval):
    """Sample a binary stationary stress chain after one exposure interval."""
    _probability(q, "q")
    if not isfinite(mixing_rate) or mixing_rate <= 0:
        raise ValueError("mixing_rate must be finite and positive")
    if not isfinite(interval) or interval <= 0:
        raise ValueError("interval must be finite and positive")
    clock = mixing_rate * interval
    stay = exp(-clock)
    refresh = -expm1(-clock)
    return stay * np.eye(2) + refresh * np.array([[1 - q, q], [1 - q, q]])


def stress_count_distribution(rounds, q, mixing_rate, interval):
    """P(K=m) for the count of stressed checkpoints under stationary preparation."""
    _nonnegative_integer(rounds, "rounds")
    transition = stress_transition(q, mixing_rate, interval)
    # Mass by accumulated count and latest stress bit, before the first checkpoint.
    mass = np.zeros((rounds + 1, 2))
    mass[0] = [1 - q, q]
    for completed in range(1, rounds + 1):
        arriving = mass[:completed] @ transition
        mass[:completed, 0] = arriving[:, 0]
        mass[1:completed + 1, 1] = arriving[:, 1]
        mass[0, 1] = 0
    return mass.sum(axis=1)


def survival_frontier(count_distribution, contact):
    """Optimal survival at reserves 0 through N, for N checkpoint rounds."""
    mass = np.asarray(count_distribution, dtype=float)
    if (mass.ndim != 1 or not mass.size or np.any(mass < 0)
            or not np.isclose(mass.sum(), 1, rtol=0, atol=1e-12)):
        raise ValueError("count_distribution must be a probability vector")
    weighted = mass * contact.success ** np.arange(mass.size)
    return np.minimum(1.0, np.cumsum(weighted))


def minimum_tokens(count_distribution, contact, target):
    """First reserve attaining the target, or None if the target is infeasible."""
    _probability(target, "target")
    frontier = survival_frontier(count_distribution, contact)
    tokens = int(np.searchsorted(frontier, target))
    return tokens if tokens < len(frontier) else None


def minimum_energy_for_known_stress_count(stresses, target, rate_ceiling, duration):
    """Sharp energy infimum over token ladders, conditional on exactly m stresses.

    Zero means an infimum when stresses > 0 and arbitrarily small positive
    quanta are required. Infinity means no finite reserve attains the target.
    """
    _nonnegative_integer(stresses, "stresses")
    _probability(target, "target")
    fraction = _contact_fraction(rate_ceiling, duration)
    if not stresses or target == 0:
        return 0.0
    if fraction == 0:
        return inf
    log_required_bias = log(target) / stresses - log(fraction)
    if log_required_bias >= 0:
        return inf
    return stresses * max(0.0, log_required_bias - log(-expm1(log_required_bias)))
