"""Build the continuum claims, numerical comparison and verification receipt."""
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from playwright.sync_api import sync_playwright

from build_research_docs import page
from verify_research_docs import Links

ROOT = Path(__file__).resolve().parents[1]


def read(relative):
    return json.loads((ROOT/relative).read_bytes())


def digest(relative):
    return sha256((ROOT/relative).read_bytes()).hexdigest()


def euler_constants():
    gamma_upper = Fraction(2587, 10000)
    rate_upper = 17 * Fraction(512, 3465) * Fraction(22, 7) * Fraction(8, 25)**3
    assert rate_upper < gamma_upper
    occupation_polynomial = sum((-1)**(k+1) * gamma_upper**k * 2**(k+1)
                                / factorial(k+1) for k in range(1, 6))
    occupation_upper = Fraction(439, 1000)
    assert occupation_polynomial < occupation_upper
    squared_drift = Fraction(25, 48)
    lipschitz = Fraction(375, 32)
    diffusivity = Fraction(1, 40)
    rows = []
    for step, cap in [(Fraction(1, 100), Fraction(481, 1000)),
                      (Fraction(1, 200), Fraction(339, 1000)),
                      (Fraction(1, 400), Fraction(239, 1000))]:
        squared_bound = lipschitz**2 * occupation_upper * (
            3*step/8 + squared_drift*step**2/(24*diffusivity))
        assert squared_bound <= cap**2
        rows.append({'maximum_step': str(step), 'squared_bound': str(squared_bound),
                     'round_up': str(cap)})
    return {'gamma_upper': str(gamma_upper), 'occupation_polynomial': str(occupation_polynomial),
            'occupation_upper': str(occupation_upper), 'rows': rows}


def large_time_constants():
    """Check the rational relaxations in the ground-state coercivity proof."""
    # 4.1 < sqrt(17) < 4.15 puts r² between 2.4 and 2.5.
    assert Fraction(41, 10)**2 < 17 < Fraction(83, 20)**2
    assert Fraction(12, 5) < (Fraction(31, 10)/2)**2
    assert (Fraction(63, 20)/2)**2 < Fraction(5, 2)
    lower = Fraction(144, 25)*Fraction(3, 10)**4/8
    upper = 10*Fraction(5, 4)*81
    assert lower > Fraction(1, 200) and upper < 1100
    assert all(comb(4,b)*Fraction(5,16)**b <= Fraction(5,4) for b in range(5))
    # e < 11/4, by summing through order 4 and geometrically bounding k>=5.
    e_upper = sum(Fraction(1, factorial(k)) for k in range(5)) + Fraction(1,120)/(1-Fraction(1,6))
    assert e_upper < Fraction(11,4) and e_upper**2 < 8 and e_upper**4 < 81
    assert 4*Fraction(4,25)**3 > Fraction(1,64)
    assert Fraction(1,200)*Fraction(3,4)**4 > Fraction(1,800)
    assert Fraction(16,1024)*Fraction(3,4)**4 > Fraction(1,256)
    assert Fraction(12,5)**2/2 > 2 and Fraction(5,2)**4 < 40
    local_path_bound = 51200
    anchor_passive = 20*3*local_path_bound
    all_passive = 2*anchor_passive + 2*19*local_path_bound
    h_coefficient = 1100*(2*Fraction(65*8000,36) + 40*all_passive)
    p_coefficient = 1100*40*16384
    scale = 400_000_000_000
    assert max(h_coefficient,p_coefficient) < scale
    prefactor = 32*81*2
    assert 20*324 < 81**2
    assert sum(Fraction(10**k,factorial(k)) for k in range(11)) > prefactor
    return {'spectator_excess_exponent_lower': str(Fraction(1,scale)),
            'passive_form_coefficient': str(h_coefficient), 'paid_form_coefficient': str(p_coefficient),
            'ratio_prefactor_upper': prefactor, 'sufficient_long_horizon': 10*scale}


def parameter_budget():
    """Sufficient target for a future validated 38-node comparator calculation."""
    nodes = 38
    remainder = 1400*Fraction(33,4)**nodes/factorial(nodes)
    assert remainder < Fraction(18,100_000_000)
    # sqrt(2n-1) is less than 9. Allow the same radius at the spectator point.
    assert 2*nodes-1 < 9**2
    radius = Fraction(1,100_000)
    total = 10*radius+remainder
    assert total < Fraction(101,1_000_000)
    return {'chebyshev_nodes': nodes, 'analytic_remainder_upper': str(remainder),
            'proposed_node_radius': str(radius), 'combined_error_budget_upper': str(total),
            'validated_node_enclosures_available': False}


def figures(short, fourier, weighted, lattice):
    background, ink, muted = '#0d141d', '#e6edf7', '#adc0d6'
    colors = ['#72dbc1', '#edbd75', '#d899e8']
    plt.rcParams.update({'figure.facecolor': background, 'axes.facecolor': background,
        'text.color': ink, 'axes.labelcolor': ink, 'axes.edgecolor': muted,
        'xtick.color': muted, 'ytick.color': muted, 'axes.spines.top': False,
        'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.6), layout='constrained')
    bounds = short['decimal_for_display']
    limit = bounds['certified_short_horizon']
    times = np.linspace(0, .001, 301)
    curve = (-bounds['initial_cubic_loss_coefficient_lower']
             + bounds['normalized_remainder_linear']*times
             + bounds['normalized_remainder_quadratic']*times**2)/6
    axes[0].plot(times*1e3, curve, color=colors[0], label='Rigorous upper bound')
    axes[0].axhline(0, color=muted, linewidth=1)
    axes[0].axvspan(0, limit*1e3, color=colors[0], alpha=.13, label='Proved loss interval')
    axes[0].set(xlabel='Horizon × 1,000', ylabel='Spectator minus passive, divided by $t^3$',
                title='Brownian short-horizon theorem')
    axes[0].legend(facecolor=background, labelcolor=ink, fontsize=8)
    calculations = fourier['rows']
    axes[1].plot([row['cutoff'] for row in calculations],
                 [100*row['paired_gap'] for row in calculations], 'o-', color=colors[0],
                 label='Fourier, versus fixed α = 1')
    pairs = {}
    for case in weighted['cases']:
        pairs.setdefault(case['cutoff'], {})[case['policy']] = case['survival'][-1]
    cutoffs = sorted(pairs)
    axes[1].plot(cutoffs, [100*(pairs[k]['auto']-pairs[k][1.0]) for k in cutoffs],
                 's--', color=colors[1], label='Weighted Fourier, versus α = 1')
    lower, upper = lattice['spectator_minus_best_fixed_interval']
    axes[1].axhspan(100*lower, 100*upper, color=colors[2], alpha=.4,
                    label='48³ lattice, versus best fixed')
    axes[1].set(xlabel='Fourier cutoff per axis', ylabel='Gap (percentage points)',
                title='Independent spatial checks at T = 2')
    axes[1].legend(facecolor=background, labelcolor=ink, fontsize=8)
    steps = np.geomspace(1e-11, .02, 301)
    euler = (375/32)*np.sqrt(.439 * (3*steps/8 + (25/48)*steps**2/(24*.025)))
    axes[2].loglog(steps, euler, color=colors[1], label='Per-policy Euler bias bound')
    axes[2].axhline(lower/2, color=colors[2], linestyle='--', label='Half the lattice margin')
    axes[2].axvline(.0025, color=muted, linestyle=':', label='Finest retained Euler step')
    axes[2].set(xlabel='Maximum Euler step', ylabel='Absolute probability',
                title='Bias bound versus target precision')
    axes[2].legend(facecolor=background, labelcolor=ink, fontsize=8)
    for extension in ['png', 'pdf']:
        fig.savefig(ROOT/f'results/figures/continuum-organization-bounds.{extension}', dpi=180)
    plt.close(fig)


def verify_page(target, pending_receipt):
    links = Links()
    links.feed(target.read_text(encoding='utf-8'))
    for url in links.urls:
        parts = urlsplit(url)
        if parts.scheme or parts.netloc:
            continue
        path = target.parent/unquote(parts.path) if parts.path else target
        assert path.exists() or path.resolve() == pending_receipt, url
        if path.resolve() == target and parts.fragment:
            assert unquote(parts.fragment) in links.ids, url
    records = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel='chrome', headless=True)
        for width in [1440, 390]:
            tab = browser.new_page(viewport={'width': width, 'height': 1000})
            errors = []
            tab.on('pageerror', lambda error: errors.append(str(error)))
            tab.on('console', lambda message: errors.append(message.text) if message.type == 'error' else None)
            tab.goto(target.as_uri(), wait_until='networkidle')
            tab.evaluate('document.fonts.ready')
            result = tab.evaluate('''() => ({equations:document.querySelectorAll('.math').length,
                rendered:document.querySelectorAll('.math .katex').length,
                mathErrors:document.querySelectorAll('.math[data-error]').length,
                brokenImages:[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length,
                width:innerWidth,scrollWidth:document.documentElement.scrollWidth})''')
            assert result['equations'] == result['rendered'] and result['equations'] > 0, (errors, result)
            assert not errors and not result['mathErrors'] and not result['brokenImages'], (errors, result)
            assert result['scrollWidth'] <= width, result
            result['console_errors'] = errors
            tab.screenshot(path=str(ROOT/f'results/document-previews/continuum-organization-{width}.png'))
            records.append(result)
            tab.close()
        browser.close()
    return records


def main():
    short = read('results/continuum-short-time.json')
    fourier = read('results/continuum-fourier.json')
    weighted = read('results/continuum-weighted.json')
    residual_screen = read('results/continuum-residual-screen.json')
    lattice = read('results/organization-verification.json')['certificates'][-1]
    assert short['status'] == 'PASS'
    for record in [short, fourier, weighted, residual_screen]:
        for name, value in record['source_sha256'].items():
            assert digest(name) == value, name
    euler = euler_constants()
    long_time = large_time_constants()
    budget = parameter_budget()
    figures(short, fourier, weighted, lattice)
    source = (ROOT/'docs/continuum-organization.md').read_text(encoding='utf-8')
    fine = fourier['rows'][-1]
    formula = (f"At Fourier cutoff {fine['cutoff']} per axis, it gives spectator survival "
               f"{fine['values']['spectator']:.10f} and constant α=1 survival {fine['values']['1.0']:.10f}. "
               f"Their difference is {fine['paired_gap']:.10f}. The number is a comparison to α=1. "
               "Its equality to the gap against the optimized constant family requires a separate upper bound on that family. "
               "The Fourier integrals use ordinary numerical quadrature, so refinement cannot certify their error. "
               "[Fourier derivation](../research/continuum-spectral.md) · [Numerical receipt](../results/continuum-fourier.json)")
    fine_weighted = {case['policy']:case['survival'][-1] for case in weighted['cases']
                     if case['cutoff'] == max(c['cutoff'] for c in weighted['cases'])}
    residual = (f"A second Fourier calculation works directly in the Gibbs inner product. It gives a paired difference "
                f"of {fine_weighted['auto']-fine_weighted[1.0]:.10f} at its finest cutoff. "
                "The exact error identity uses the product of the forward and backward residuals in a norm adapted to diffusion. "
                "An independent numerical screen estimates those products at about 2.5e−6 at cutoff 12 for both arms, "
                "before accounting for unresolved frequencies, time quadrature, and the correction for inexact orthogonality. "
                "[Weighted calculation](../results/continuum-weighted.json) · [Residual screen](../results/continuum-residual-screen.json)\n\n"
                "A rigorous parameter interpolation estimate reduces a future global comparator calculation to 38 Chebyshev node values. "
                "For example, if each node and the spectator value has a validated radius of 1e−5, the combined value and interpolation budget is below 1.01e−4. "
                "A verified maximum of the resulting interpolant must still clear that margin. The current calculations do not supply those nodal enclosures. "
                "[Residual identity and interval-integration recipe](../research/continuum-residual.md) · "
                "[Derivation of the 38-node bound](../research/continuum-spectral.md)")
    source = source.replace('@@FOURIER@@', formula).replace('@@RESIDUAL@@', residual)
    assert '@@' not in source
    target = ROOT/'docs/CONTINUUM-ORGANIZATION.html'
    receipt_path = ROOT/'results/continuum-organization-verification.json'
    page(target.name, source, 'Brownian continuum survival and its error bounds')
    browser = verify_page(target, receipt_path)
    env = {**os.environ, 'PYTHONPATH': str(ROOT/'src')}
    test = subprocess.run([sys.executable,'-B','-m','unittest','discover','-s','tests','-v'],
                          cwd=ROOT, env=env, capture_output=True, text=True, check=True)
    test_output = test.stdout+test.stderr
    files = ['docs/continuum-organization.md', 'docs/CONTINUUM-ORGANIZATION.html',
             'research/continuum-coupling.md', 'research/continuum-residual.md',
             'research/continuum-spectral.md', 'results/continuum-short-time.json',
             'results/continuum-fourier.json', 'results/continuum-residual-screen.json', 'results/organization-verification.json',
             'results/organization-continuum-check.json', 'scripts/build_continuum_organization.py',
             'results/figures/continuum-organization-bounds.png',
             'results/figures/continuum-organization-bounds.pdf',
             *[p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'results').glob('continuum-*weighted*.json'))],
             *[p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'scripts').glob('*continuum*.py'))]]
    receipt = {'generated_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Brownian short- and long-horizon results, and numerics at T=2.',
        'time_two_status': 'INCONCLUSIVE_CONTINUUM_SIGN',
        'short_time_exact_rational_checks': 'PASS', 'euler_bias_exact_rational_checks': euler,
        'long_time_exact_rational_checks': long_time,
        'parameter_interpolation_budget': budget,
        'test_count': int(re.search(r'Ran (\d+) tests?', test_output)[1]),
        'test_output': test_output, 'browser_checks': browser,
        'sha256': {name: digest(name) for name in sorted(set(files))}}
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: receipt[key] for key in ['time_two_status','test_count','browser_checks']}, indent=2))


if __name__ == '__main__':
    main()
