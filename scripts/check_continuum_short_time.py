"""Replay exact rational bounds for short-horizon Brownian survival."""
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def constants():
    """Evaluate the inequalities used in the Duhamel remainder proof."""
    radius = Fraction(8, 25)
    diffusivity = Fraction(1, 40)
    depth = 6
    fuel = 4
    repair = 16
    paid = 4
    gate_mean = Fraction(2, 125)
    mean_exit = Fraction(14, 5)
    maximum_drift = Fraction(3, 4)
    generator_bound = 9000
    rate_cap = fuel * paid
    gate_gradient = 8 / radius
    gate_laplacian = 24 / radius**2

    # The rational bounds 3 < pi < 22/7 enclose the radial integrals.
    assert Fraction(512, 3465) * Fraction(22, 7) * radius**3 < gate_mean
    assert Fraction(3, 2) + (repair * fuel + 1 + rate_cap) * gate_mean < mean_exit
    assert (4 * diffusivity * depth / radius)**2 * Fraction(4, 27) < maximum_drift**2
    assert 2 - depth * Fraction(3, 4)**2 < 0
    f2_low = -2 - Fraction(repair * fuel, 2)
    f2_high = Fraction(3, 2) + repair * fuel
    base_exit = Fraction(3, 2) + repair * fuel + 1
    generator_bound_check = repair * fuel * (diffusivity * gate_laplacian + maximum_drift * gate_gradient)
    generator_bound_check += base_exit * (f2_high - f2_low)
    assert generator_bound_check < generator_bound

    contact_integral = 12 * radius**3 * sum(
        Fraction((-1)**k * comb(8, k), 2**(2*k+3) * (2*k+3)) for k in range(9)
    )
    cubic_lower = rate_cap * repair * contact_integral
    linear_remainder = (2 * rate_cap * repair * mean_exit + 2 * rate_cap * generator_bound * gate_mean) / 4
    quadratic_remainder = 2 * rate_cap * generator_bound * mean_exit / 20
    horizon = Fraction(1, 1250)
    margin = cubic_lower - linear_remainder * horizon - quadratic_remainder * horizon**2
    assert margin > Fraction(1, 10)
    return {
        'radius': radius, 'diffusivity': diffusivity, 'uniform_gate_mean_upper': gate_mean,
        'first_jump_mean_rate_upper': mean_exit, 'closed_drift_upper': maximum_drift,
        'base_generator_f2_absolute_upper': generator_bound_check,
        'restricted_contact_integral_lower': contact_integral,
        'initial_cubic_loss_coefficient_lower': cubic_lower,
        'normalized_remainder_linear': linear_remainder,
        'normalized_remainder_quadratic': quadratic_remainder,
        'certified_short_horizon': horizon, 'remaining_margin_at_horizon': margin,
    }


def main():
    values = constants()
    inputs = ['scripts/check_continuum_short_time.py', 'src/purposeful/organization.py',
              'research/organization-preregistration.md']
    record = {
        'status': 'PASS',
        'scope': 'Exact rational inequalities for the Brownian short-horizon proof. The sign at time two is a separate question.',
        'result': 'For 0<t<=1/1250, spectator loses to alpha=0 by at least t^3/60. Each fixed alpha>0 loses to alpha=0 by at least alpha*t^3/60.',
        'rational': {key: str(value) for key, value in values.items()},
        'decimal_for_display': {key: float(value) for key, value in values.items()},
        'source_sha256': {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
    }
    (ROOT/'results/continuum-short-time.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
