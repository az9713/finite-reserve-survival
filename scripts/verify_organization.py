"""Independent interval/SSA verification of the finite-grid organization experiment."""

import argparse
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import heapq
import json
from math import ceil, log, sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "tests")]
from purposeful.organization import Organization, build_system
from test_organization import DAMAGE, PAID, PASSIVE, REPAIR, policy_multiplier, rational, reference_graph


LOW = Context(prec=45, rounding=ROUND_FLOOR)
HIGH = Context(prec=45, rounding=ROUND_CEILING)
ZERO, ONE = Decimal(0), Decimal(1)
FLOAT_EPS = np.finfo(float).eps
UNDERFLOW_PAD = 1e-300


def fraction_interval(value):
    value = Fraction(value)
    numerator, denominator = Decimal(value.numerator), Decimal(value.denominator)
    return LOW.divide(numerator, denominator), HIGH.divide(numerator, denominator)


def exp_interval(lo, hi):
    # Decimal.exp is correctly rounded half-even. One adjacent Decimal is outward
    # regardless of the configured context rounding; zero is evaluated exactly.
    return (ONE if not lo else LOW.next_minus(LOW.exp(lo)),
            ONE if not hi else HIGH.next_plus(HIGH.exp(hi)))


@lru_cache(maxsize=100_000)
def rate_interval(prefactor, exponent):
    lo, hi = fraction_interval(prefactor)
    exp_lo, exp_hi = exp_interval(*fraction_interval(exponent))
    return LOW.multiply(lo, exp_lo), HIGH.multiply(hi, exp_hi)


def to_float(value, upper):
    rounded = float(value)
    if upper and Decimal.from_float(rounded) < value:
        return np.nextafter(rounded, np.inf)
    if not upper and Decimal.from_float(rounded) > value:
        return np.nextafter(rounded, -np.inf)
    return rounded


def matrices(model, n):
    """Directed interval arithmetic for all three independently enumerated kernels."""
    states, rows, launch, _, _ = reference_graph(model, n, folded=True, killed=True)
    index = {state: i for i, state in enumerate(states)}
    size = len(states)
    labels = ('base', 'paid', 'spectator')
    blocks = {name: [] for name in labels}
    rates = {name: 1 for name in labels}

    for row, (state, branches) in enumerate(zip(states, rows)):
        split = {name: {} for name in labels}
        for branch in branches:
            destination = index.get(branch.destination)
            mapping = ({'paid': Fraction(1),
                        'spectator': policy_multiplier('spectator', model.fuel, state, branch.event)}
                       if abs(branch.event) == PAID else {'base': Fraction(1), 'spectator': Fraction(1)})
            for name, multiplier in mapping.items():
                if not multiplier:
                    continue
                lo, hi = rate_interval(branch.prefactor * multiplier, branch.exponent)
                old_lo, old_hi = split[name].get(destination, (ZERO, ZERO))
                split[name][destination] = LOW.add(old_lo, lo), HIGH.add(old_hi, hi)
        for name in labels:
            exit_lo = exit_hi = ZERO
            for lo, hi in split[name].values():
                exit_lo = LOW.add(exit_lo, lo)
                exit_hi = HIGH.add(exit_hi, hi)
            rates[name] = max(rates[name], int(exit_hi.to_integral_value(rounding=ROUND_CEILING)))
            blocks[name].append((row, split[name], exit_lo, exit_hi))

    output = {}
    for name in labels:
        rate = Decimal(rates[name])
        values_lo, values_hi, columns, index_ptr = [], [], [], [0]
        for row, branches, exit_lo, exit_hi in blocks[name]:
            branches[row] = LOW.subtract(rate, exit_hi), HIGH.subtract(rate, exit_lo)
            for col, (lo, hi) in sorted((col, value) for col, value in branches.items() if col is not None):
                columns.append(col)
                values_lo.append(to_float(LOW.divide(lo, rate), upper=False))
                values_hi.append(to_float(HIGH.divide(hi, rate), upper=True))
            index_ptr.append(len(columns))
        output[name] = (csr_matrix((values_lo, columns, index_ptr), shape=(size, size)),
                        csr_matrix((values_hi, columns, index_ptr), shape=(size, size)))
    launch_indices = [i for i, value in enumerate(launch) if value]
    launch_probability = tuple(np.array([to_float(fraction_interval(launch[i])[d], upper=d == 1)
                                         for i in launch_indices]) for d in (0, 1))
    return output, rates, launch_indices, launch_probability


def poisson_intervals(mean, order):
    """Mass intervals through order and a geometric bound for the discarded tail."""
    x_lo, x_hi = fraction_interval(mean)
    lo, hi = exp_interval(x_hi.copy_negate(), x_lo.copy_negate())
    lower, upper = [lo], [hi]
    for k in range(1, order + 2):
        lo = LOW.divide(LOW.multiply(lo, x_lo), Decimal(k))
        hi = HIGH.divide(HIGH.multiply(hi, x_hi), Decimal(k))
        lower.append(lo)
        upper.append(hi)
    ratio = HIGH.divide(x_hi, Decimal(order + 2))
    tail = HIGH.divide(upper.pop(), LOW.subtract(ONE, ratio))
    lower.pop()
    return lower, upper, tail


def poisson_order(mean):
    return ceil(float(mean) + 12 * sqrt(float(mean) + 1) + 20)


def rounded_positive(value, term_count, upper):
    """Enclose a nonnegative float expression using a conservative gamma bound.

    At most 2*term_count+12 operations affect any sum before this padding.
    32*epsilon per counted term dominates the relative gamma bound, including
    float ratio formation and the padding arithmetic. UNDERFLOW_PAD also
    covers underflow of intermediate products and the final update.
    """
    if term_count > 1_000_000:
        raise ValueError('float roundoff bound requires at most one million terms')
    relative = 32 * FLOAT_EPS * (term_count + 2)
    if upper:
        return np.minimum(1., value * (1 + relative) + UNDERFLOW_PAD)
    return np.maximum(0., value * (1 - relative) - UNDERFLOW_PAD)


def poisson_arrays(mean):
    lo, hi, tail = poisson_intervals(mean, poisson_order(mean))
    return (np.array([to_float(value, upper=False) for value in lo]),
            np.array([to_float(value, upper=True) for value in hi]), to_float(tail, upper=True))


def project(vector, indices, probability, upper):
    return rounded_positive(probability @ vector[indices], len(indices), upper)


def paid_coefficients(base, paid, means, indices, probability):
    base_weights = poisson_arrays(means[0])
    m_max = poisson_order(means[1])
    n_max = len(base_weights[0]) - 1
    size = base[0].shape[0]
    coeff = []
    for d, upper in enumerate((False, True)):
        v = np.empty((n_max + 1, size))
        u = np.empty(m_max + 1)
        row_terms = max(np.diff(base[d].indptr).max(), np.diff(paid[d].indptr).max())
        for m in range(m_max + 1):
            terms = np.empty(n_max + 1)
            for n in range(n_max + 1):
                if n + m == 0:
                    v[n] = 1.
                else:
                    if m == 0:
                        value = base[d] @ v[n - 1]
                    elif n == 0:
                        value = paid[d] @ v[n]
                    else:
                        value = (n / (n + m)) * (base[d] @ v[n - 1]) + (m / (n + m)) * (paid[d] @ v[n])
                    v[n] = rounded_positive(value, row_terms, upper)
                terms[n] = project(v[n], indices, probability[d], upper)
            expectation = base_weights[d] @ terms + (base_weights[2] if upper else 0.)
            u[m] = rounded_positive(expectation, n_max + 2, upper)
        coeff.append([Decimal.from_float(value) for value in u])
    return coeff, base_weights[2]


def one_kernel_survival(kernel, mean, indices, probability):
    lo, hi, tail = poisson_arrays(mean)
    bounds = []
    for d, weights in enumerate((lo, hi)):
        upper = d == 1
        terms = np.empty(len(weights))
        row_terms = np.diff(kernel[d].indptr).max()
        vector = np.ones(kernel[d].shape[0])
        for k in range(len(terms)):
            terms[k] = project(vector, indices, probability[d], upper)
            vector = rounded_positive(kernel[d] @ vector, row_terms, upper)
        bounds.append(float(rounded_positive(weights @ terms + (tail if upper else 0), len(terms), upper)))
    return bounds, tail


def signed_expectation(weights, values):
    lo = hi = ZERO
    for p_lower, p_upper, value_lower, value_upper in zip(*weights, *values):
        left = LOW.multiply(p_upper if value_lower < 0 else p_lower, value_lower)
        right = HIGH.multiply(p_lower if value_upper < 0 else p_upper, value_upper)
        lo, hi = LOW.add(lo, left), HIGH.add(hi, right)
    return lo, hi


def certificate(coeff, paid_mean, tolerance):
    """Cover [0,1] by Taylor boxes using paid-Poisson second differences."""
    lower, upper = coeff
    means_lo, means_hi = fraction_interval(paid_mean)
    differences = ([LOW.subtract(low, high) for low, high in zip(lower[1:], upper[:-1])],
                   [HIGH.subtract(high, low) for high, low in zip(upper[1:], lower[:-1])])
    absolute_second = []
    for i in range(len(lower) - 2):
        low = LOW.subtract(differences[0][i + 1], differences[1][i])
        high = HIGH.subtract(differences[1][i + 1], differences[0][i])
        absolute_second.append(max(low.copy_negate(), high))
    _, _, curvature_tail = poisson_intervals(paid_mean, len(lower) - 3)
    curvature = HIGH.multiply(HIGH.multiply(means_hi, means_hi), HIGH.add(max(absolute_second),
                                                                          HIGH.multiply(Decimal(2), curvature_tail)))

    @lru_cache(maxsize=None)
    def point(alpha):
        lo, hi, tail = poisson_intervals(alpha * paid_mean, len(lower) - 1)
        s_lower, s_upper = signed_expectation((lo, hi), (lower, upper))
        p_lower, p_upper = signed_expectation((lo[:-1], hi[:-1]), differences)
        derivative_tail = HIGH.add(hi[-1], tail)
        p_lower = LOW.subtract(p_lower, derivative_tail)
        p_upper = HIGH.add(p_upper, derivative_tail)
        derivative_lower = LOW.multiply(p_lower, means_hi if p_lower < 0 else means_lo)
        derivative_upper = HIGH.multiply(p_upper, means_lo if p_upper < 0 else means_hi)
        return s_lower, HIGH.add(s_upper, tail), derivative_lower, derivative_upper

    best = max((point(alpha)[0], alpha) for alpha in (Fraction(0), Fraction(1)))

    def upper_bound(left, right):
        nonlocal best
        center = (left + right) / 2
        _, r = fraction_interval((right - left) / 2)
        s_lower, s_upper, derivative_lower, derivative_upper = point(center)
        derivative = max(derivative_lower.copy_negate(), derivative_upper)
        best = max(best, (s_lower, center))
        curvature_term = HIGH.divide(HIGH.multiply(curvature, HIGH.multiply(r, r)), Decimal(2))
        bound = HIGH.add(HIGH.add(s_upper, HIGH.multiply(derivative, r)), curvature_term)
        return bound.copy_negate(), left, right

    pending = [upper_bound(Fraction(0), Fraction(1))]
    while HIGH.subtract(pending[0][0].copy_negate(), best[0]) > tolerance:
        if len(pending) >= 10_000:
            raise RuntimeError('cannot reach the requested tolerance within 10,000 cover intervals')
        _, left, right = heapq.heappop(pending)
        midpoint = (left + right) / 2
        heapq.heappush(pending, upper_bound(left, midpoint))
        heapq.heappush(pending, upper_bound(midpoint, right))
    global_upper = pending[0][0].copy_negate()
    endpoint_margin = None
    endpoint_cells = None
    if best[1] in (0, 1):
        # Prove monotonicity if a modest uniform signed-derivative cover suffices.
        for cells in (8, 32, 128, 512):
            cell_margin = ONE
            _, radius = fraction_interval(Fraction(1, 2 * cells))
            variation = HIGH.multiply(curvature, radius)
            for j in range(cells):
                _, _, d_lower, d_upper = point(Fraction(2 * j + 1, 2 * cells))
                aligned_derivative = d_lower if best[1] == 1 else d_upper.copy_negate()
                cell_margin = min(cell_margin, LOW.subtract(aligned_derivative, variation))
                if cell_margin <= 0:
                    break
            if cell_margin > 0:
                global_upper = min(global_upper, point(best[1])[1])
                endpoint_margin, endpoint_cells = cell_margin, cells
                break
    _, _, maximum_paid_tail = poisson_intervals(paid_mean, len(lower) - 1)
    result = {
        'lower': to_float(best[0], upper=False), 'upper': to_float(global_upper, upper=True),
        'attaining_sample_alpha': float(best[1]), 'sample_points': point.cache_info().currsize,
        'cover_intervals': len(pending), 'absolute_second_derivative_upper': to_float(curvature, upper=True),
        'maximum_paid_count_tail_upper': to_float(maximum_paid_tail, upper=True),
        'narrowest_cover_interval': float(min(right - left for _, left, right in pending)),
    }
    if endpoint_margin is not None:
        result.update(certified_attaining_alpha=float(best[1]), derivative_cover_intervals=endpoint_cells,
                      minimum_signed_derivative_lower=to_float(endpoint_margin, upper=False))
    return result


def exact_lattice_paths(model, n, policy, times, size, seed):
    """Independent direct jump paths, with the first failure absorbing for observation."""
    states, branches, launch, geometry, _ = reference_graph(model, n, folded=True, killed=False)
    index = {state: i for i, state in enumerate(states)}
    width = max(map(len, branches))
    rates = np.zeros((len(states), width))
    destinations = np.tile(np.arange(len(states))[:, None], (1, width))
    events = np.zeros((len(states), width), dtype=np.int8)
    for i, (state, choices) in enumerate(zip(states, branches)):
        for j, choice in enumerate(choices):
            rates[i, j] = choice.rate() * float(policy_multiplier(policy, model.fuel, state, choice.event))
            destinations[i, j] = index[choice.destination]
            events[i, j] = choice.event
    cutoffs = np.cumsum(rates, axis=1)
    state_data = np.asarray([state[:3] for state in states])
    local_energy = np.array([q * float(geometry[position][1]) for _, _, q, position in states])
    chemical_energy = state_data[:, 1] * log(8) - (2 - state_data[:, 0]) * log(2)
    total_energy = chemical_energy + local_energy
    rng = np.random.default_rng(seed)
    current = rng.choice(len(states), p=np.asarray(launch, dtype=float), size=size)
    elapsed = np.zeros(size)
    failure = np.full(size, np.inf)
    ledger = np.zeros((size, 2 * PAID + 1), dtype=np.int32)
    path_energy = np.zeros(size)
    live = np.arange(size)
    while len(live):
        exit_rate = cutoffs[current[live], -1]
        waiting = rng.exponential(1 / exit_rate)
        happens = waiting + elapsed[live] <= times[-1]
        if not np.any(happens):
            break
        live, waiting = live[happens], waiting[happens]
        elapsed[live] += waiting
        indices = current[live]
        total = cutoffs[indices, -1]
        marks = np.minimum(rng.random(len(live)) * total, np.nextafter(total, 0.))
        selected = (marks[:, None] >= cutoffs[indices]).sum(axis=1)
        event = events[indices, selected]
        new_indices = destinations[indices, selected]
        ledger[live, event + PAID] += 1
        path_energy[live] += total_energy[new_indices] - total_energy[indices]
        current[live] = new_indices
        collapsed = state_data[new_indices, 0] == 0
        failure[live[collapsed]] = elapsed[live[collapsed]]
        live = live[~collapsed]
    net = ledger[:, PAID + 1:] - ledger[:, PAID - 1::-1]
    h, b, q = state_data[current].T
    if not (np.array_equal(b, model.fuel - net[:, REPAIR - 1] - net[:, PAID - 1]) and
            np.array_equal(h, 1 - net[:, DAMAGE - 1] + net[:, REPAIR - 1]) and
            np.array_equal(q, net[:, PASSIVE - 1] + net[:, PAID - 1])):
        raise RuntimeError('path ledger did not reproduce final states')
    expected_energy = (-log(2) * net[:, DAMAGE - 1] - log(4) * net[:, REPAIR - 1]
                       - log(8) * net[:, PAID - 1] + local_energy[current])
    discrepancy = float(np.max(np.abs(path_energy - expected_energy)))
    if discrepancy > 1e-10:
        raise RuntimeError('path ledger did not reproduce the net energy change')
    sorted_failure = np.sort(failure)
    result = 1 - np.searchsorted(sorted_failure, times, side='right') / size
    return result, dict(seed=seed, paths=size, channel_counts_by_signed_event=ledger.sum(axis=0).tolist(),
                        signed_event_order=list(range(-PAID, PAID + 1)),
                        paths_with_collapse=int(np.isfinite(failure).sum()),
                        max_energy_ledger_error=discrepancy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--levels', nargs='+', type=int, default=[4, 8, 16, 24, 32, 48])
    parser.add_argument('--horizon', type=Fraction, default=Fraction(2))
    parser.add_argument('--tolerance', type=Decimal, default=Decimal('2e-6'))
    parser.add_argument('--paths', type=int, default=20000)
    parser.add_argument('--path-level', type=int, help='Run path verification at this certified level')
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 'organization-verification.json')
    args = parser.parse_args()
    if args.horizon <= 0 or args.tolerance <= 0 or args.paths < 1 or any(n < 2 or n % 2 for n in args.levels):
        parser.error('positive time, tolerance and path count, and even grids of at least two are required')
    if args.path_level is not None and args.path_level not in args.levels:
        parser.error('--path-level must be one of --levels')
    source_paths = [Path('src/purposeful/organization.py'), Path('src/purposeful/spatial_repair.py'),
                    Path('tests/test_organization.py'), Path('scripts/verify_organization.py')]
    source_fingerprints = {path.as_posix(): sha256((ROOT / path).read_bytes()).hexdigest() for path in source_paths}
    model = Organization()
    runs = []
    for n in args.levels:
        print(f'Certifying n={n} ...', flush=True)
        kernels, rates, indices, probability = matrices(model, n)
        coeff, base_tail = paid_coefficients(kernels['base'], kernels['paid'],
                                             [rates[name] * args.horizon for name in ('base', 'paid')],
                                             indices, probability)
        fixed = certificate(coeff, rates['paid'] * args.horizon, args.tolerance)
        auto, auto_tail = one_kernel_survival(kernels['spectator'], rates['spectator'] * args.horizon,
                                              indices, probability)
        system = build_system(model, n)
        time = [float(args.horizon)]
        alpha = fixed['attaining_sample_alpha']
        check_auto = float(system.survival('spectator', time)[0])
        check_fixed = float(system.survival(alpha, time)[0])
        if not (auto[0] <= check_auto <= auto[1] and abs(check_fixed - fixed['lower']) < 2e-9):
            raise RuntimeError('independent and production survival values disagree')
        row = dict(n=n, horizon=time[0], states=len(system.initial()), uniformization_rates=rates,
                   fixed=fixed, spectator_interval=auto, chebyshev_spectator=check_auto,
                   chebyshev_fixed_at_attaining_sample=check_fixed,
                   base_poisson_tail_upper=base_tail, spectator_poisson_tail_upper=auto_tail,
                   max_paid_coefficient_interval_width=max(to_float(HIGH.subtract(hi, lo), upper=True)
                                                           for lo, hi in zip(*coeff)),
                   spectator_minus_best_fixed_interval=[float(np.nextafter(auto[0] - fixed['upper'], -np.inf)),
                                                         float(np.nextafter(auto[1] - fixed['lower'], np.inf))])
        runs.append(row)
        print(f"n={n} gap enclosure {row['spectator_minus_best_fixed_interval']}", flush=True)

    n = args.path_level if args.path_level is not None else (32 if 32 in args.levels else args.levels[-1])
    system = build_system(model, n)
    times = np.linspace(0, float(args.horizon), 201)
    fixed_policy = next(record['fixed']['attaining_sample_alpha'] for record in runs if record['n'] == n)
    policies = (0., fixed_policy, 'spectator')
    envelope = sqrt(log(2 * len(policies) / .01) / (2 * args.paths))
    path_checks = []
    for arm, policy in enumerate(policies):
        seed = 19891101 + 1729 * arm
        actual, ledger = exact_lattice_paths(model, n, policy, times, args.paths, seed)
        expected = system.survival(policy, times)
        maximum_error = float(np.max(np.abs(actual - expected)))
        if maximum_error > envelope:
            raise RuntimeError('path and production survival disagree beyond the prescribed band')
        path_checks.append(dict(policy=policy, maximum_grid_difference=maximum_error,
                                endpoint_path_survival=float(actual[-1]), endpoint_solver_survival=float(expected[-1]),
                                passed_dkw_band=True, **ledger))
    document = dict(model=vars(model), design='alpha in [0,1] with common hardware and launch',
                    arithmetic='directed Decimal rate and tail intervals plus padded nonnegative binary64 sums',
                    continuum_discretization_error_enclosed=False, source_sha256=source_fingerprints, certificates=runs,
                    paths=dict(n=n, observations=len(times), dkw_familywise_confidence=.99,
                               dkw_probability_half_width=envelope, checks=path_checks))
    if any(sha256((ROOT / name).read_bytes()).hexdigest() != expected
           for name, expected in source_fingerprints.items()):
        raise RuntimeError('source files changed during verification')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')
    print(args.output)


if __name__ == '__main__':
    main()
