"""Sharp checkpoint survival with one reserve for sensing, action, and repair."""
from dataclasses import dataclass

import numpy as np

from purposeful.reserve import RepairContact, _nonnegative_integer, _probability

IDLE, ENGAGE = 0, 1
STATIONARY_PRIOR = 2


@dataclass(frozen=True)
class Device:
    capacity: int
    sensing_cost: int
    actuation_cost: int
    harvest: int
    opportunity: float
    repair: RepairContact
    fidelity: float = 1.0
    persistence: float = 0.0

    def __post_init__(self):
        for name in ("capacity", "sensing_cost", "actuation_cost", "harvest"):
            _nonnegative_integer(getattr(self, name), name)
        if not self.sensing_cost or not self.actuation_cost:
            raise ValueError("sensing and actuation must each cost at least one unit")
        _probability(self.opportunity, "opportunity")
        _probability(self.persistence, "persistence")
        if not 0.5 <= self.fidelity <= 1:
            raise ValueError("fidelity must be in [1/2, 1]")

    @property
    def opportunity_by_prior(self):
        q = self.opportunity
        return (q * (1 - self.persistence), q + self.persistence * (1 - q), q)


@dataclass(frozen=True)
class SurvivalPolicy:
    """Tables indexed by remaining rounds, previous environment, and reserve.

    Use STATIONARY_PRIOR when the preceding environment is unknown at launch.
    action_after_readout has a final index (zero or one) for the observation.
    """

    survival: np.ndarray
    sense: np.ndarray
    blind_action: np.ndarray
    action_after_readout: np.ndarray


class _Checkpoint:
    """Expected continuation for the two observations and the blind decision."""

    def __init__(self, device):
        self.capacity = device.capacity
        self.cost = device.sensing_cost
        reserve = np.arange(device.capacity + 1)
        after_action = reserve[None, :] - np.array([[0], [device.actuation_cost]])
        # Axes are action, current environment, reserve. A surviving round
        # discloses the previous environment through its reserve change.
        self.next_reserve = np.stack([
            np.stack([reserve, reserve - 1]),
            np.stack([after_action[ENGAGE] - 1, after_action[ENGAGE] + device.harvest]),
        ]).clip(0, device.capacity)
        can_repair = after_action > 0
        repair_success = device.repair.success
        self.live = np.stack([
            np.stack([np.ones_like(reserve), repair_success * can_repair[IDLE]]),
            np.stack([repair_success * can_repair[ENGAGE], after_action[ENGAGE] >= 0]),
        ])
        opportunity = np.array(device.opportunity_by_prior)
        self.prior = np.stack([1 - opportunity, opportunity], axis=1)  # prior, current
        self.readout = np.array([
            [device.fidelity, 1 - device.fidelity],
            [1 - device.fidelity, device.fidelity],
        ])  # observation, current

    def options(self, future):
        next_environment = np.array([0, 1])[None, :, None]
        payoff = self.live * future[next_environment, self.next_reserve]
        blind = np.einsum("px,axb->pab", self.prior, payoff)
        informed = np.einsum("px,yx,axb->pyab", self.prior, self.readout, payoff)
        # Infeasible actions have zero payoff. Break zero ties in favor of idle.
        best_blind = np.argmax(blind, axis=1)
        values_blind = np.max(blind, axis=1)
        values_sense = np.full_like(values_blind, -np.inf)
        actions_sense = np.zeros((*values_blind.shape, 2), dtype=np.int8)
        if self.cost <= self.capacity:
            funded = self.capacity + 1 - self.cost
            values_sense[:, self.cost:] = np.max(informed[:, :, :, :funded], axis=2).sum(axis=1)
            best_informed = np.argmax(informed[:, :, :, :funded], axis=2)
            actions_sense[:, self.cost:] = np.moveaxis(best_informed, 1, -1)
        return values_blind, values_sense, best_blind, actions_sense

    def fixed_options(self, future):
        blind, sense, _, _ = self.options(future)
        # A precommitted probe falls back to blind action when unaffordable.
        return blind, np.where(np.isfinite(sense), sense, blind)


def optimal_policy(device, rounds):
    """Exact Bellman frontier and an attaining policy, computed in O(N*C)."""
    _nonnegative_integer(rounds, "rounds")
    checkpoint = _Checkpoint(device)
    shape = (rounds + 1, 3, device.capacity + 1)
    survival = np.ones(shape)
    sense = np.zeros(shape, dtype=bool)
    blind_action = np.zeros(shape, dtype=np.int8)
    action_after_readout = np.zeros((*shape, 2), dtype=np.int8)
    for remaining in range(1, rounds + 1):
        blind, informed, blind_action[remaining], action_after_readout[remaining] = (
            checkpoint.options(survival[remaining - 1])
        )
        sense[remaining] = informed > blind
        survival[remaining] = np.maximum(blind, informed)
    return SurvivalPolicy(survival, sense, blind_action, action_after_readout)


def evaluate_fixed_schedule(device, schedule):
    """Evaluate an ex ante probe schedule, with optimal conditional actuation."""
    checkpoint = _Checkpoint(device)
    future = np.ones((3, device.capacity + 1))
    for sense in reversed(tuple(schedule)):
        future = checkpoint.fixed_options(future)[int(bool(sense))]
    return future


def best_fixed_frontier(device, rounds):
    """Exact best fixed schedule per launch state, pruning dominated suffixes.

    At worst the number of retained vectors is exponential in rounds.
    Returns the survival envelope, schedule attaining it per state, and
    the number of nondominated suffixes retained at each horizon.
    """
    _nonnegative_integer(rounds, "rounds")
    checkpoint = _Checkpoint(device)
    candidates = [(np.ones((3, device.capacity + 1)), ())]
    counts = [len(candidates)]
    for _ in range(rounds):
        kept = []
        for future, suffix in candidates:
            for mode, values in enumerate(checkpoint.fixed_options(future)):
                if any(np.all(other >= values) for other, _ in kept):
                    continue
                kept = [(other, schedule) for other, schedule in kept if not np.all(values >= other)]
                kept.append((values, (bool(mode), *suffix)))
        candidates = kept
        counts.append(len(candidates))
    vectors = np.stack([values for values, _ in candidates])
    winner = np.argmax(vectors, axis=0)
    schedules = {
        (prior, reserve): candidates[winner[prior, reserve]][1]
        for prior, reserve in np.ndindex(winner.shape)
    }
    return vectors.max(axis=0), schedules, counts


def minimum_reserve(policy, target, prior=STATIONARY_PRIOR):
    """Least initial reserve meeting the target, or None if it cannot be met."""
    _probability(target, "target")
    feasible = np.flatnonzero(policy.survival[-1, prior] >= target)
    return int(feasible[0]) if feasible.size else None


def block_failure_bound(device):
    """Return (block length, uniform failure lower bound) for fixed device parameters."""
    p_min = min(q for p in device.opportunity_by_prior for q in (p, 1 - p))
    length = device.capacity // device.sensing_cost + 1
    return length, p_min ** length * (1 - device.repair.success)
