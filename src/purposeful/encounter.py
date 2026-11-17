"""Classical diffusion-to-capture benchmarks, in micrometres and seconds.

The exact geometry is one mobile particle outside a fixed central target,
with a reflecting outer spherical boundary. It is a search benchmark,
separate from the six-species chemical network.
"""
from math import exp, isfinite, pi, sqrt

from scipy.integrate import quad


def shell_radius(volume, capture_radius):
    """Outer radius for accessible shell volume, excluding the target."""
    if not all(isfinite(x) and x > 0 for x in (volume, capture_radius)):
        raise ValueError('volume and capture radius must be finite and positive')
    outer = (capture_radius**3 + 3*volume/(4*pi))**(1/3)
    if outer <= capture_radius:
        raise ValueError('the shell is too thin for floating-point resolution')
    return outer


def mean_capture_time(volume, diffusion, capture_radius):
    """Exact MFPT for a uniform initial point in a reflecting spherical shell."""
    if not isfinite(diffusion) or diffusion <= 0:
        raise ValueError('diffusion must be finite and positive')
    a, outer = capture_radius, shell_radius(volume, capture_radius)
    polynomial = 5*outer**3 + 6*outer**2*a + 3*outer*a*a + a**3
    return (outer-a)**2*polynomial/(15*diffusion*a*(outer**2+outer*a+a*a))


def capture_before_closure(radius, volume, diffusion, capture_radius, closure_rate):
    """Laplace transform of the capture time at an individual starting radius.

    The observation window closes at an independent exponential time.
    The analytic function solves D laplacian(q) - gamma q = 0, with
    q(a)=1 and q'(R)=0. Small arguments use a cancellation-free series.
    """
    if not isfinite(diffusion) or diffusion <= 0:
        raise ValueError('diffusion must be finite and positive')
    if not isfinite(closure_rate) or closure_rate < 0:
        raise ValueError('closure rate must be finite and nonnegative')
    a, outer = capture_radius, shell_radius(volume, capture_radius)
    if not a <= radius <= outer:
        raise ValueError('starting radius must lie in the accessible shell')
    if closure_rate == 0:
        return 1.0
    k = sqrt(closure_rate/diffusion)
    if k*outer < .05:
        def v(r):
            s2 = (k*(outer-r))**2
            sinhc = 1+s2/6+s2**2/120+s2**3/5040+s2**4/362880+s2**5/39916800
            difference = s2/3+s2**2/30+s2**3/840+s2**4/45360+s2**5/3991680
            return r/outer*sinhc+difference
    else:
        def v(r):
            return ((1-1/(k*outer))*exp(-k*(r-a))
                    +(1+1/(k*outer))*exp(-k*(2*outer-r-a)))
    return a/radius*v(radius)/v(a)


def uniform_capture_probability(volume, diffusion, capture_radius, closure_rate):
    """Average the capture-before-closure probability over the uniform shell."""
    a, outer = capture_radius, shell_radius(volume, capture_radius)
    integral = quad(lambda r: 4*pi*r*r*capture_before_closure(
        r, volume, diffusion, a, closure_rate), a, outer, epsabs=1e-12, epsrel=1e-11)[0]
    return integral/volume


def partial_reactivity_mean(volume, diffusion, capture_radius, surface_reactivity):
    """Exact shell MFPT with a radiation boundary of reactivity in um/s."""
    if not isfinite(surface_reactivity) or surface_reactivity <= 0:
        raise ValueError('surface reactivity must be finite and positive')
    return (mean_capture_time(volume, diffusion, capture_radius)
            + volume/(4*pi*capture_radius**2*surface_reactivity))
