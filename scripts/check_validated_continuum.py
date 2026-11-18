"""Independent rational parameter-cover replay for the continuum certificate.

This module deliberately separates parameter coverage from PDE validation.
An interpolant and 38 numerical node values alone are not a continuum proof.
The final receipt must also identify validated spatial/time residual bounds.
The Bernstein maximum calculation uses only exact Fraction arithmetic.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
import heapq
import json
from math import comb, factorial
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def fraction(value):
    """Read an exact dyadic, rational, or decimal string; reject binary floats."""
    if isinstance(value, dict):
        if 'denominator' in value:
            return F(int(value['numerator']), int(value['denominator']))
        return F(int(value['numerator']), 2 ** int(value['bits']))
    if isinstance(value, float):
        raise TypeError('certificate scalars must be exact strings, integers, or dyadics')
    return F(value)


def exp_upper(x, degree=64):
    """Exact Taylor upper bound, using the geometric bound for the positive tail."""
    x = F(x)
    if x < 0 or degree + 2 <= x:
        raise ValueError('Taylor tail requires 0 <= x < degree + 2')
    term = F(1)
    value = term
    for k in range(1, degree + 1):
        term *= x / k
        value += term
    next_term = term * x / (degree + 1)
    return value + next_term / (1 - x / (degree + 2))


def analytic_constants():
    """Replay all rational constants needed by the 38-root cover theorem."""
    assert exp_upper(F(4)) < 55
    # ||f||^2 <= (5/9)(1+e^4), ||ell||^2 = 59049/4.
    squared_norm_product_upper = F(5, 9) * 56 * F(59049, 4)
    assert squared_norm_product_upper < 700 ** 2
    count = 38
    operator_upper = F(33, 2)
    horizon = F(2)
    # For roots of T_n(2 alpha - 1), max |product(alpha-alpha_j)| = 2^(1-2n).
    interpolation_error = F(1400) * F(33, 4) ** count / factorial(count)
    assert interpolation_error < F(18, 100_000_000)
    # A rational upper enclosure for sqrt(2n-1)=sqrt(75).
    amplification = F(8_660_254_038, 1_000_000_000)
    assert amplification ** 2 > 2 * count - 1
    return {
        'node_count': count,
        'horizon': str(horizon),
        'norm_product_upper': '700',
        'norm_product_squared_upper': str(squared_norm_product_upper),
        'paid_operator_norm_upper': str(operator_upper),
        'parameter_lipschitz_upper': str(700 * horizon * operator_upper),
        'interpolation_remainder_upper': str(interpolation_error),
        'node_uncertainty_amplification_upper': str(amplification),
    }


def evaluate_power(coefficients, x):
    value = F(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def restrict_power(coefficients, left, right):
    """Exact power coefficients of p(left+(right-left)*s), 0<=s<=1."""
    left, right = fraction(left), fraction(right)
    width = right - left
    return tuple(width**j * sum((a * comb(k, j) * left**(k-j)
                                for k, a in enumerate(coefficients) if k >= j), F(0))
                 for j in range(len(coefficients)))


def power_to_bernstein(coefficients):
    """Convert p(x)=sum a_k x^k to Bernstein coefficients on [0,1]."""
    n = len(coefficients) - 1
    if n < 0:
        raise ValueError('empty polynomial')
    return tuple(sum((coefficients[k] * F(comb(i, k), comb(n, k))
                      for k in range(i + 1)), F(0)) for i in range(n + 1))


def split_bernstein(coefficients):
    """Exact de Casteljau midpoint restriction, preserving both end values."""
    row = tuple(coefficients)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = tuple((a + b) / 2 for a, b in zip(row, row[1:]))
        left.append(row[0])
        right.append(row[-1])
    return tuple(left), tuple(reversed(right))


def polynomial_maximum(coefficients, tolerance=F(1, 10**12), max_splits=100_000):
    """Certify max p on [0,1] between exact rationals by convex hulls.

    Every live heap entry is a dyadic interval. Replacing one entry by its
    two children preserves their union, hence covers the full domain.
    A Bernstein polynomial is a convex combination of its coefficients.
    Endpoints and split midpoints give attained lower bounds.
    """
    coefficients = tuple(map(fraction, coefficients))
    tolerance = fraction(tolerance)
    if tolerance <= 0:
        raise ValueError('positive tolerance required')
    bernstein = power_to_bernstein(coefficients)
    lower = max(bernstein[0], bernstein[-1])
    maximizing_sample = F(0) if bernstein[0] >= bernstein[-1] else F(1)
    queue = [(-max(bernstein), 0, F(0), F(1), bernstein)]
    serial, splits, pruned = 0, 0, 0
    while -queue[0][0] - lower > tolerance and splits < max_splits:
        negative_upper, _, left, right, row = heapq.heappop(queue)
        if -negative_upper <= lower:
            pruned += 1
            continue
        a, b = split_bernstein(row)
        middle = (left + right) / 2
        if a[-1] > lower:
            lower, maximizing_sample = a[-1], middle
        for lo, hi, child in ((left, middle, a), (middle, right, b)):
            upper = max(child)
            # Keep a degenerate attained lower leaf too, so the heap is never empty.
            if upper < lower:
                pruned += 1
                continue
            serial += 1
            heapq.heappush(queue, (-upper, serial, lo, hi, child))
        splits += 1
    upper = max(lower, -queue[0][0])
    return {
        'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
        'attaining_sample': str(maximizing_sample),
        'requested_tolerance': str(tolerance), 'tolerance_met': upper - lower <= tolerance,
        'splits': splits, 'active_intervals': len(queue), 'pruned_intervals': pruned,
        'full_unit_interval_covered': True,
        'arithmetic': 'exact rational Bernstein conversion and de Casteljau bisection',
    }


def self_check():
    """Adversarial examples for endpoints, interior maxima, and cancellation."""
    examples = [([F(3)], F(3)), ([F(0), F(1)], F(1)),
                ([F(0), F(4), F(-4)], F(1)),
                ([F(-1, 9), F(2, 3), F(-1)], F(0)),
                ([F(-1), F(32), F(-160), F(256), F(-128)], F(1))]
    for coefficients, exact in examples:
        result = polynomial_maximum(coefficients, F(1, 10**14))
        assert F(result['lower']) <= exact <= F(result['upper'])
        assert result['tolerance_met']
    try:
        fraction(0.1)
    except TypeError:
        pass
    else:
        raise AssertionError('binary float must not silently enter proof arithmetic')
    return {'known_polynomial_maxima': len(examples), 'binary_float_rejected': True}


def arb_endpoint_fraction(value, upper=True):
    mantissa, exponent = (value.upper() if upper else value.lower()).man_exp()
    return F(int(mantissa)) * F(2) ** int(exponent)


def independent_interpolant_check(record, precision=320):
    """Recompute actual-root shifts and interpolation conversion in Arb.

    The supplied coefficient error is not trusted. We recompute a sufficient
    error, using integer Chebyshev basis coefficients and a separate DCT.
    This is independent of the trial/PDE implementation.
    """
    from flint import arb, ctx
    rows = record['node_records']
    count = len(rows)
    if count != 38:
        raise ValueError('the analytic constants are for 38 nodes')
    with ctx.workprec(precision):
        def ball(q):
            q = fraction(q)
            return arb(q.numerator) / q.denominator
        shifts, centers = [], []
        for j, row in enumerate(rows):
            alpha = F(int(row['alpha_numerator']), 2**int(row['alpha_bits']))
            if not F(0) <= alpha <= F(1):
                raise ValueError('dyadic node outside the fixed comparison family')
            true_root = (arb(1) + (arb(2*j+1) / (2*count)).cos_pi()) / 2
            shifts.append(arb_endpoint_fraction(abs(true_root - ball(alpha))))
            centers.append(arb(int(row['center_numerator'])) / 2**int(row['center_bits']))
        # Integer polynomial coefficients of T_k(2 alpha - 1).
        bases = [[1], [-1, 2]]
        for k in range(2, count):
            previous, before = bases[-1], bases[-2]
            current = [0] * (k + 1)
            for p, value in enumerate(previous):
                current[p] -= 2*value
                current[p+1] += 4*value
            for p, value in enumerate(before):
                current[p] -= value
            bases.append(current)
        power = [arb(0) for _ in range(count)]
        for k, basis in enumerate(bases):
            coefficient = sum(value * (arb(k*(2*j+1))/(2*count)).cos_pi()
                              for j, value in enumerate(centers))
            coefficient *= arb(1 if k == 0 else 2) / count
            for p, integer in enumerate(basis):
                power[p] += integer * coefficient
        supplied = list(map(fraction, record['power_coefficients']))
        if len(supplied) > count:
            raise ValueError('degree exceeds 37')
        supplied += [F(0)] * (count - len(supplied))
        conversion = sum(abs(c-ball(q)) for c, q in zip(power, supplied))
        return {
            'root_distances_upper': shifts,
            'coefficient_error_upper': arb_endpoint_fraction(conversion),
            'arb_precision_bits': precision,
        }


def check_parameter_receipt(record, tolerance=F(1, 10**12), lebesgue_certificate=None):
    """Cover all alpha from dyadic interpolant coefficients and node radii.

    Required schema: power_coefficients, coefficient_error_upper,
    node_records[38] with exact dyadic alpha, center, radius_upper, and spectator
    with lower/upper. Radius includes validated PDE error, while root_distance
    is independently recomputed between the dyadic alpha and the exact root.
    This function verifies arithmetic, not provenance of supplied PDE bounds.
    """
    constants = analytic_constants()
    nodes = record['node_records']
    if len(nodes) != constants['node_count']:
        raise ValueError('exactly 38 Chebyshev-root node records required')
    coefficients = list(map(fraction, record['power_coefficients']))
    if len(coefficients) > 38:
        raise ValueError('interpolant degree exceeds 37')
    recomputed = independent_interpolant_check(record)
    conversion_error = recomputed['coefficient_error_upper']
    if conversion_error < 0:
        raise ValueError('coefficient error cannot be negative')
    radii = []
    for row, distance in zip(nodes, recomputed['root_distances_upper']):
        radius = fraction(row['radius_upper'])
        if radius < 0 or distance < 0:
            raise ValueError('node errors cannot be negative')
        radii.append(radius + fraction(constants['parameter_lipschitz_upper']) * distance)
    maximization = polynomial_maximum(coefficients, tolerance)
    amplification = fraction(constants['node_uncertainty_amplification_upper'])
    if lebesgue_certificate is not None:
        if lebesgue_certificate['count'] != 38 or not lebesgue_certificate['full_unit_interval_covered']:
            raise ValueError('incorrect Lebesgue cover')
        amplification = fraction(lebesgue_certificate['upper'])
    node_error = amplification * max(radii)
    interpolation_error = fraction(constants['interpolation_remainder_upper'])
    fixed_upper = F(maximization['upper']) + conversion_error + node_error + interpolation_error
    spectator_lower = fraction(record['spectator']['lower'])
    spectator_upper = fraction(record['spectator']['upper'])
    if not F(0) <= spectator_lower <= spectator_upper <= F(1):
        raise ValueError('invalid spectator survival interval')
    gap_lower = spectator_lower - fixed_upper
    return {
        'analytic_constants': constants,
        'polynomial_maximum': maximization,
        'coefficient_error_upper': str(conversion_error),
        'arb_precision_bits': recomputed['arb_precision_bits'],
        'actual_root_shifts_recomputed': True,
        'node_center_interpolant_recomputed': True,
        'maximum_root_adjusted_node_radius': str(max(radii)),
        'node_uncertainty_amplification_upper': str(amplification),
        'node_interpolation_uncertainty_upper': str(node_error),
        'fixed_family_global_upper': str(fixed_upper),
        'spectator_lower': str(spectator_lower),
        'spectator_minus_best_fixed_lower': str(gap_lower),
        'strict_gap_arithmetically_verified': gap_lower > 0,
        'advantage_time_interval_open_conditional_on_pde_bounds':
            [str(F(2)-gap_lower), str(F(2)+gap_lower)] if gap_lower > 0 else None,
        'time_window_justification': 'Only killing has rate1 in h1; every survival decrement lies in [-dt,0], so each policy difference is1-Lipschitz.',
        'pde_validation_replayed': False,
        'scope': 'Independent exact parameter-cover arithmetic; supplied PDE node enclosures need a separate validated residual replay.',
    }


def synthetic_cover_check():
    """Exercise the complete cover on a known constant, without research data."""
    from flint import arb, ctx
    rows = []
    with ctx.workprec(192):
        for j in range(38):
            root = (1 + (arb(2*j+1)/76).cos_pi())/2
            midpoint = arb_endpoint_fraction(root.mid()) * 2**60
            numerator = midpoint.numerator // midpoint.denominator
            rows.append({'alpha_numerator': str(numerator), 'alpha_bits': 60,
                         'center_numerator': '1', 'center_bits': 2,
                         'radius_upper': '0'})
    data = {'node_records': rows, 'power_coefficients': ['1/4'],
            'coefficient_error_upper': '0',
            'spectator': {'lower': '3/10', 'upper': '3/10'}}
    result = check_parameter_receipt(data)
    assert F(result['polynomial_maximum']['upper']) == F(1, 4)
    assert F(result['fixed_family_global_upper']) < F(25_000_018, 100_000_000)
    assert result['strict_gap_arithmetically_verified']
    bad = dict(data, node_records=rows[:-1])
    try:
        check_parameter_receipt(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('incomplete node cover must be rejected')
    return {'known_constant': '1/4', 'node_count': 38,
            'root_shifts_and_center_polynomial_recomputed': True,
            'full_parameter_interval_replayed': True,
            'incomplete_cover_rejected': True,
            'scope': 'Synthetic checker regression, not a PDE survival result'}


def certified_lebesgue_constant(count=38, precision=320, bits=160, endpoint_bits=40):
    """Compute a stronger nodal uncertainty factor by exact interval coverage.

    On every interval between Chebyshev roots, all Lagrange basis signs are
    known from their roots. The Lebesgue function there is one ordinary
    signed-sum interpolation polynomial. We enclose that polynomial by Arb,
    quantize it to dyadics, and maximize it on an outward dyadic interval by
    exact Bernstein arithmetic. This is a finite proof, without assuming the
    location of the Lebesgue maximum from an external theorem.
    """
    from flint import arb, ctx
    with ctx.workprec(precision):
        roots = [(1+(arb(2*j+1)/(2*count)).cos_pi())/2 for j in range(count)]
        cosines = [[(arb(k*(2*j+1))/(2*count)).cos_pi()
                    for j in range(count)] for k in range(count)]
        bases = [[1], [-1, 2]]
        for k in range(2, count):
            current = [0]*(k+1)
            for j, value in enumerate(bases[-1]):
                current[j] -= 2*value
                current[j+1] += 4*value
            for j, value in enumerate(bases[-2]):
                current[j] -= value
            bases.append(current)
        upper_root, lower_root = [], []
        for root in roots:
            low = arb_endpoint_fraction(root, False)*2**endpoint_bits
            high = arb_endpoint_fraction(root, True)*2**endpoint_bits
            lower_root.append(F(low.numerator//low.denominator, 2**endpoint_bits))
            upper_root.append(F(-((-high.numerator)//high.denominator), 2**endpoint_bits))
        intervals = []
        for m in range(count+1):
            # roots are in descending order; m roots lie to the right.
            signs = [(-1)**(m-(1 if i < m else 0)-i) for i in range(count)]
            power = [arb(0) for _ in range(count)]
            for k in range(count):
                coefficient = sum(sign*cosine for sign, cosine in zip(signs, cosines[k]))
                coefficient *= arb(1 if k == 0 else 2)/count
                for j, value in enumerate(bases[k]):
                    power[j] += value*coefficient
            exact, error = [], arb(0)
            for value in power:
                midpoint = arb_endpoint_fraction(value.mid())*2**bits
                integer = midpoint.numerator//midpoint.denominator
                q = F(integer, 2**bits)
                exact.append(q)
                error += abs(value-arb(integer)/2**bits)
            low = F(0) if m == count else lower_root[m]
            high = F(1) if m == 0 else upper_root[m-1]
            maximum = polynomial_maximum(restrict_power(exact, low, high), F(1, 10**10))
            bound = F(maximum['upper']) + arb_endpoint_fraction(error)
            intervals.append({'interval': [str(low), str(high)],
                              'roots_to_right': m, 'upper': str(bound),
                              'bernstein_splits': maximum['splits'],
                              'polynomial_conversion_error_upper': str(arb_endpoint_fraction(error))})
        global_upper = max(F(row['upper']) for row in intervals)
        return {'count': count, 'upper': str(global_upper),
                'upper_decimal_for_display': float(global_upper),
                'intervals': intervals, 'arb_precision_bits': precision,
                'arithmetic': 'Known Lagrange signs, Arb interpolation, exact rational interval restriction and Bernstein maximum',
                'full_unit_interval_covered': True}


def replay_receipt_chain(record):
    """Independently check archived inputs and scalar enclosure assembly.

    This verifies exact trial endpoints/initial values, all recorded source and
    artifact hashes, and Eu*Ev/radius/observable interval arithmetic. It does
    not recompute spatial matrices or their residual contractions.
    """
    import numpy as np
    from itertools import combinations_with_replacement
    links = record['node_receipts']
    if len(links) != 39 or len({row['path'] for row in links}) != 39:
        raise ValueError('39 unique PDE node receipts required')
    required_sources = {'scripts/run_validated_continuum.py','scripts/validate_continuum_nodes.py',
                        'src/purposeful/validated_time.py','src/purposeful/validated_linear.py',
                        'src/purposeful/validated_spatial.py','src/purposeful/validated_residual.py',
                        'src/purposeful/validated_refinement.py'}
    checked_sources, checked_coefficients, checked_caches = set(), 0, {}
    for name, expected in record.get('source_sha256',{}).items():
        if sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('comparison source hash changed: '+name)
    for index, link in enumerate(links):
        path, trial = ROOT/link['path'], ROOT/link['trial']
        if sha256(path.read_bytes()).hexdigest() != link['sha256']:
            raise ValueError('PDE receipt hash changed: '+link['path'])
        if sha256(trial.read_bytes()).hexdigest() != link['trial_sha256']:
            raise ValueError('PDE trial hash changed: '+link['trial'])
        node = json.loads(path.read_text(encoding='utf-8'))
        if node['status'] != 'validated continuum node enclosure' or node['horizon'] != 2:
            raise ValueError('incorrect node enclosure scope')
        if node['trial']['sha256'] != link['trial_sha256']:
            raise ValueError('PDE receipt refers to another trial')
        if not required_sources.issubset(node['source_sha256']):
            raise ValueError('PDE receipt omits a critical source dependency')
        if node.get('source_unchanged_during_validation') is not True:
            raise ValueError('PDE receipt lacks the source-freeze check')
        for name, expected in node['source_sha256'].items():
            if sha256((ROOT/name).read_bytes()).hexdigest() != expected:
                raise ValueError('source hash changed: '+name)
            checked_sources.add(name)
        cache_directory = node.get('matrix_cache')
        cache_records = node.get('consumed_matrix_cache_sha256',{})
        if not cache_directory or not cache_records:
            raise ValueError('PDE receipt omits consumed matrix cache provenance')
        cache_root = ROOT/cache_directory
        identity = sha256((ROOT/'src/purposeful/validated_residual.py').read_bytes()
                          +(ROOT/'src/purposeful/validated_linear.py').read_bytes()
                          +(ROOT/'src/purposeful/validated_spatial.py').read_bytes()).hexdigest()
        if (cache_root/'SOURCE-SHA256.txt').read_text().strip() != identity:
            raise ValueError('matrix cache source identity changed')
        for name, expected in cache_records.items():
            cache_file = cache_root/name
            key = cache_file.relative_to(ROOT).as_posix()
            if key not in checked_caches:
                checked_caches[key] = sha256(cache_file.read_bytes()).hexdigest()
            if checked_caches[key] != expected or cache_file.with_suffix('.sha256').read_text().strip() != expected:
                raise ValueError('consumed matrix cache changed: '+key)
        lo, hi = fraction(node['survival_lower']), fraction(node['survival_upper'])
        eu, ev = fraction(node['primal']['energy_upper']), fraction(node['dual']['energy_upper'])
        product = fraction(node['residual_product_upper'])
        if eu < 0 or ev < 0 or product < 0 or product*product < eu*ev:
            raise ValueError('invalid residual product arithmetic')
        # These fields are exact binary64 data, unlike certificate decimal text.
        center = F.from_float(node['corrected_observable_center'])
        rounding = F.from_float(node['corrected_observable_radius'])
        if rounding < 0 or lo > center-rounding-product or hi < center+rounding+product:
            raise ValueError('observable interval omits residual or rounding error')
        if index == 0:
            if node['policy'] != 'spectator':
                raise ValueError('candidate receipt is not the spectator policy')
            if lo != fraction(record['spectator']['lower']) or hi != fraction(record['spectator']['upper']):
                raise ValueError('spectator interval differs from PDE node')
        else:
            interpolation = record['node_records'][index-1]
            alpha = F(int(interpolation['alpha_numerator']), 2**int(interpolation['alpha_bits']))
            if fraction(node['policy']) != alpha:
                raise ValueError('fixed PDE receipt has the wrong exact alpha')
            value = F(int(interpolation['center_numerator']), 2**int(interpolation['center_bits']))
            radius = fraction(interpolation['radius_upper'])
            if value-radius > lo or value+radius < hi:
                raise ValueError('interpolation node omits part of PDE enclosure')
        with np.load(trial, allow_pickle=False) as archive:
            bits = int(archive['bits'])
            if bits < 2 or bits > 100:
                raise ValueError('unexpected exact trial precision')
            knots = tuple(F(int(a), int(b)) for a, b in zip(archive['knot_numerators'], archive['knot_denominators']))
            if knots[0] != 0 or knots[-1] != 2 or any(b <= a for a,b in zip(knots,knots[1:])):
                raise ValueError('invalid trial time partition')
            modes, states = archive['modes'], archive['states']
            expected_modes = np.asarray(list(combinations_with_replacement(range(int(modes.max())+1),3)))
            expected_states = np.asarray([(h,b,q) for h in (1,2) for b in range(5) for q in (0,1)])
            if not np.array_equal(modes,expected_modes) or not np.array_equal(states,expected_states):
                raise ValueError('trial mode or chemical convention mismatch')
            launch = list(map(tuple,states)).index((1,4,0))
            for label in ('primal','dual'):
                integers = archive[label]
                if integers.dtype.kind not in ('i','u'):
                    raise ValueError('archive coefficients must be exact integers')
                if integers.shape[0] != len(knots)-1 or integers.shape[2:] != (20,len(modes)):
                    raise ValueError('inconsistent trial dimensions')
                maximum = max(int(integers.max()),-int(integers.min()))
                # Use machine-integer sums only after proving no prefix can
                # overflow; otherwise fall back to exact Python bigints.
                safe_int64 = (integers.dtype.kind == 'i' and integers.dtype.itemsize <= 8
                              and maximum*integers.shape[1] <= np.iinfo(np.int64).max)
                if not safe_int64:
                    integers = np.asarray(integers,dtype=object)
                initial = np.zeros((20,len(modes)),dtype=object)
                if label == 'primal':
                    initial[:,0] = 2**bits
                else:
                    initial[launch,0] = 59049*2**(bits-2)
                if not np.array_equal(integers[0,0],initial):
                    raise ValueError('incorrect exact trial initial condition')
                endpoints = integers[:-1].sum(axis=1,dtype=np.int64 if safe_int64 else object)
                if not np.array_equal(endpoints,integers[1:,0]):
                    raise ValueError('discontinuous archived time curve')
                checked_coefficients += integers.size
    return {'node_receipts_checked':39,'sources_checked':sorted(checked_sources),
            'consumed_matrix_cache_files_checked':len(checked_caches),
            'exact_integer_coefficients_checked':checked_coefficients,
            'exact_launch_and_continuity_replayed':True,
            'residual_product_and_observable_interval_assembly_replayed':True,
            'all_artifact_and_source_hashes_match':True,
            'spatial_residual_contractions_recomputed':False,
            'scope':'Independent input/provenance and final scalar assembly replay; the spatial/time node validator supplies the PDE residual energy premise.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/continuum-independent-check.json')
    parser.add_argument('--tolerance', default='1/1000000000000')
    parser.add_argument('--lebesgue', action='store_true', help='also certify a sharper 38-node interpolation stability factor')
    args = parser.parse_args()
    result = {'generated_utc': datetime.now(timezone.utc).isoformat(),
              'status': 'ANALYTIC_CONSTANTS_AND_CHECKER_VERIFIED',
              'self_check': self_check(), 'analytic_constants': analytic_constants(),
              'synthetic_cover_check': synthetic_cover_check(),
              'source_sha256': {name: sha256((ROOT/name).read_bytes()).hexdigest()
                                for name in ['scripts/check_validated_continuum.py',
                                             'research/validated-continuum-review.md',
                                             'tests/test_validated_parameter_cover.py']},
              'continuum_sign_proved_by_this_receipt': False}
    if args.lebesgue:
        result['certified_lebesgue_constant'] = certified_lebesgue_constant()
    if args.receipt:
        supplied = json.loads(args.receipt.read_text(encoding='utf-8'))
        result['parameter_cover'] = check_parameter_receipt(supplied, F(args.tolerance), result.get('certified_lebesgue_constant'))
        if 'node_receipts' in supplied:
            result['receipt_chain'] = replay_receipt_chain(supplied)
            result['positive_primary_continuum_certificate_supported_by_replayed_chain'] = bool(
                result['parameter_cover']['strict_gap_arithmetically_verified'])
        result['input_receipt_sha256'] = sha256(args.receipt.read_bytes()).hexdigest()
        result['status'] = ('PARAMETER_COVER_AND_RECEIPT_CHAIN_REPLAYED' if 'receipt_chain' in result
                            else 'PARAMETER_ARITHMETIC_REPLAYED')
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    summary = dict(result)
    if 'certified_lebesgue_constant' in summary:
        summary['certified_lebesgue_constant'] = {
            k: v for k, v in summary['certified_lebesgue_constant'].items() if k != 'intervals'}
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
