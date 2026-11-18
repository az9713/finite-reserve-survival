"""Arb enclosures of radial Fourier integrals for the continuum model.

All declared parameters and radial-profile coefficients are rational; floats
are rejected there. SpatialTerm coefficients may be floats: these declare the
exact represented binary trial coefficients. Numerical integration uses Arb's *validated* complex
integrator.  Its tolerance is a work target, never an asserted error bound:
the returned ball is the enclosure, and nonfinite results raise an exception.

Requires python-flint >= 0.9.  The global FLINT context is process-global;
use separate processes, not threads, for concurrent calls at different precision.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import nextafter, inf
from typing import Iterable

from flint import acb, arb, ctx


def rational(x) -> Fraction:
    if isinstance(x, float):
        raise TypeError("Use Fraction, int, or a decimal string, not a float")
    return Fraction(x)


def ball(x) -> arb:
    """Enclose a declared rational exactly, with directed Arb arithmetic."""
    x = rational(x)
    return arb(x.numerator) / x.denominator


@dataclass(frozen=True)
class Elementary:
    """coefficient * core**power * exp(exponent * U), inside the ball."""
    power: int = 0
    exponent: Fraction = Fraction(0)
    coefficient: Fraction = Fraction(1)

    def __post_init__(self):
        if self.power < 0 or not isinstance(self.power, int):
            raise ValueError("core power must be a nonnegative integer")
        object.__setattr__(self, "exponent", rational(self.exponent))
        object.__setattr__(self, "coefficient", rational(self.coefficient))


def _canonical(terms):
    values = {}
    for t in terms:
        key = t.power, t.exponent
        values[key] = values.get(key, Fraction(0)) + t.coefficient
    return tuple(Elementary(p, e, c) for (p, e), c in sorted(values.items()) if c)


@dataclass(frozen=True)
class RadialProfile:
    """Piecewise analytic radial field, with the split at U(r)=0.

The exterior is a sum of exp(exponent * offset), encoded as Elementary
terms with power=0.  This retains exact parameter declarations even for
non-rational exterior constants such as exp(-2).
"""
    inner: tuple[Elementary, ...]
    outer: tuple[Elementary, ...]
    exterior: tuple[Elementary, ...] = ()

    def __post_init__(self):
        for name in ("inner", "outer", "exterior"):
            object.__setattr__(self, name, _canonical(getattr(self, name)))
        if any(t.power for t in self.exterior):
            raise ValueError("Exterior terms have no core power")

    def __mul__(self, other):
        if not isinstance(other, RadialProfile):
            q = rational(other)
            return RadialProfile(*(tuple(Elementary(t.power, t.exponent,
                q*t.coefficient) for t in terms) for terms in
                (self.inner, self.outer, self.exterior)))
        def mul(a, b):
            return tuple(Elementary(x.power+y.power, x.exponent+y.exponent,
                x.coefficient*y.coefficient) for x, y in product(a, b))
        return RadialProfile(mul(self.inner, other.inner),
            mul(self.outer, other.outer), mul(self.exterior, other.exterior))

    __rmul__ = __mul__

    def __add__(self, other):
        return RadialProfile(self.inner+other.inner, self.outer+other.outer,
                             self.exterior+other.exterior)


def profile(power=0, exponent=0, *, compact=False) -> RadialProfile:
    """core**power exp(exponent U), optionally zero outside contact."""
    t = Elementary(power, rational(exponent))
    return RadialProfile((t,), (t,), () if compact or power else (t,))


ONE = profile()
CONTACT = profile(4)
WEIGHT = profile(exponent=-1)
INVERSE_WEIGHT = profile(exponent=1)
SWITCH = RadialProfile((Elementary(4),), (Elementary(4, Fraction(-1)),))
DRIFT = profile(1, -1, compact=True)
RESIDUAL_PROFILES = {
    "one": ONE,
    "contact": CONTACT,
    "switch": SWITCH,
    "weight": WEIGHT,
    "weighted_contact": WEIGHT*CONTACT,
    "drift": DRIFT,
}


def residual_product_profiles():
    """Distinct profile products needed for squared strong sources, with aliases."""
    result = {}
    names=list(RESIDUAL_PROFILES)
    for i,name in enumerate(names):
        for other in names[i:]:
            p=RESIDUAL_PROFILES[name]*RESIDUAL_PROFILES[other]
            result.setdefault(p,[]).append((name,other))
    return result


@dataclass(frozen=True)
class SpatialTerm:
    """coefficient * x[axes] * profile(|x|) * exp(2 pi i k.x)."""
    profile: RadialProfile
    frequency: tuple[int, int, int]
    axes: tuple[int, ...] = ()
    coefficient: object = 1

    def __post_init__(self):
        if len(self.frequency) != 3 or any(not isinstance(k, int) for k in self.frequency):
            raise ValueError("frequency must contain three integers")
        if len(self.axes)>2 or any(j not in (0, 1, 2) for j in self.axes):
            raise ValueError("Cartesian moments of degree at most two are supported")


class SpatialOracle:
    """Validated radial transform/moment oracle, with per-instance caching."""
    def __init__(self, *, radius=Fraction(8,25), offset=2, depth=6,
                 bits=100, target_bits=70):
        self.radius, self.offset, self.depth = map(rational, (radius, offset, depth))
        if not 0 < self.radius < Fraction(1,2):
            raise ValueError("contact ball must lie inside the fundamental cube")
        if not 0 < self.offset < self.depth:
            raise ValueError("requires an interior U=0 branch split")
        if target_bits < 10 or bits < target_bits+16:
            raise ValueError("working precision needs at least 16 guard bits")
        self.bits, self.target_bits = bits, target_bits
        self.integral_calls = 0
        self.transform = lru_cache(maxsize=None)(self._transform)

    def exterior_value(self, p):
        with ctx.workprec(self.bits):
            return sum((ball(t.coefficient)*(ball(t.exponent)*ball(self.offset)).exp()
                        for t in p.exterior), arb(0))

    def _transform(self, p: RadialProfile, n2: int, derivative: int = 0):
        """d**derivative/dn2**derivative of the compact radial transform.

        The constant exterior is excluded.  n2 denotes |k|**2, not |2pi k|**2.
        Derivatives use entire trigonometric kernels, avoiding r=0 division.
        """
        if n2 < 0 or not isinstance(n2, int) or derivative not in (0, 1, 2):
            raise ValueError("n2 is nonnegative integer and derivative is 0, 1, or 2")
        with ctx.workprec(self.bits):
            a, delta, E = map(ball, (self.radius, self.offset, self.depth))
            split = a*(1-(delta/E).sqrt()).sqrt()
            pi = arb.pi()
            ext = self.exterior_value(p)
            answer = acb(0)
            for left, right, terms in ((arb(0), split, p.inner), (split, a, p.outer)):
                def integrand(t, analytic):
                    # Polynomial, exponential and trigonometric functions are
                    # entire in r. Thus the analytic flag needs no branch test.
                    # Integrate on exact endpoints 0,1; the algebraic split is
                    # an enclosed parameter, never an implicitly exact endpoint.
                    r = left+(right-left)*t
                    core = 1-(r/a)**2
                    U = delta-E*core**2
                    field = sum((ball(t.coefficient)*core**t.power*
                        (ball(t.exponent)*U).exp() for t in terms), acb(0))-ext
                    t = 2*pi*arb(n2).sqrt()*r
                    if derivative == 0:
                        kernel = t.sinc()
                    elif n2 == 0:
                        kernel = (-2*(pi*r)**2/3 if derivative==1
                                  else 4*(pi*r)**4/15)
                    elif derivative == 1:
                        kernel = (t.cos()-t.sinc())/(2*n2)
                    else:
                        kernel = (3*t.sinc()-3*t.cos()-t*t.sin())/(4*n2*n2)
                    return (right-left)*4*pi*r*r*field*kernel
                value = acb.integral(integrand, 0, 1,
                    rel_tol=arb(2)**(-self.target_bits),
                    abs_tol=arb(2)**(-self.target_bits), eval_limit=100000,
                    depth_limit=40, use_heap=True)
                self.integral_calls += 1
                if not value.is_finite():
                    raise ArithmeticError("Arb failed to enclose the radial integral")
                answer += value
            if not answer.imag.contains(0):
                raise ArithmeticError("radial transform unexpectedly nonreal")
            return answer.real

    def moment(self, p: RadialProfile, k=(0,0,0), axes=()):
        """Enclose integral x[axes] p(|x|) exp(-2 pi i k.x) over the torus."""
        if len(k)!=3 or any(not isinstance(n,int) for n in k):
            raise ValueError("three integer frequencies required")
        if len(axes)>2 or any(j not in (0,1,2) for j in axes):
            raise ValueError("only degree<=2 moments are supported")
        with ctx.workprec(self.bits):
            n2 = sum(n*n for n in k)
            pi = arb.pi()
            if not axes:
                compact = acb(self.transform(p,n2))
            elif len(axes)==1:
                compact = acb(0, arb(k[axes[0]])/pi*self.transform(p,n2,1))
            else:
                i,j=axes
                compact = acb(-(2*(i==j)*self.transform(p,n2,1)+
                    4*k[i]*k[j]*self.transform(p,n2,2))/(4*pi*pi))
            cube = acb(1)
            for j,n in enumerate(k):
                power = axes.count(j)
                if power == 0:
                    piece=acb(int(n==0))
                elif power == 1:
                    piece=acb(0) if n==0 else acb(0,arb((-1)**n)/(2*pi*n))
                else:
                    piece=acb(arb(1)/12 if n==0 else arb((-1)**n)/(2*pi*pi*n*n))
                cube *= piece
            return compact+self.exterior_value(p)*cube

    def integral(self, terms: Iterable[SpatialTerm]):
        """Enclose a finite complex Fourier/profile sum, without truncation."""
        with ctx.workprec(self.bits):
            return sum((acb(t.coefficient)*self.moment(t.profile,
                tuple(-k for k in t.frequency),t.axes) for t in terms),acb(0))

    def inner(self, left: Iterable[SpatialTerm], right: Iterable[SpatialTerm],
              weight=ONE):
        """Enclose integral conjugate(left)*right*weight (finite exact sums)."""
        with ctx.workprec(self.bits):
            answer=acb(0)
            for x,y in product(tuple(left),tuple(right)):
                k=tuple(a-b for a,b in zip(x.frequency,y.frequency))
                answer += acb(x.coefficient).conjugate()*acb(y.coefficient)*self.moment(
                    x.profile*y.profile*weight,k,x.axes+y.axes)
            return answer

    def shell_table(self, p: RadialProfile, max_coordinate: int, derivatives=(0,)):
        """Certified compact transforms on all attained squared-frequency shells.

        Rows retain the enclosing ball text as well as outward float midpoint /
        radius pairs for downstream coefficient perturbation estimates. The
        constant exterior must still be added for the zero Fourier mode.
        """
        shells=sorted({i*i+j*j+k*k for i in range(max_coordinate+1)
            for j in range(i,max_coordinate+1) for k in range(j,max_coordinate+1)})
        return {str(n2): {str(d): ball_receipt(self.transform(p,n2,d))
                         for d in derivatives} for n2 in shells}

    def table(self, p: RadialProfile, max_n2: int, derivative=0):
        """Float64 midpoint/radius arrays of certified radial shell coefficients.

        At derivative=0 the zero-mode exterior constant is included. At positive
        derivative p must be compact, as Cartesian derivatives of a torus-wide
        constant have cube rather than radial-ball moments. Entries for every
        integer 0<=n²<=max_n2 are included, even unattained three-square shells.
        """
        import numpy as np
        if derivative and p.exterior:
            raise ValueError("Cartesian derivative tables require compact profiles")
        center=np.empty(max_n2+1); radius=np.empty(max_n2+1)
        with ctx.workprec(self.bits):
            for n2 in range(max_n2+1):
                value=self.transform(p,n2,derivative)
                if n2==0 and derivative==0:
                    value+=self.exterior_value(p)
                center[n2],radius[n2]=midpoint_radius(value)
        return center,radius


def nonnegative_upper(x: arb) -> arb:
    """Return a valid nonnegative upper bound; reject a negative upper endpoint."""
    upper=x.upper()
    if upper < 0:
        raise ArithmeticError("negative certified upper endpoint for squared norm")
    return upper


def tail_flux_bound(total_norm_squared: arb, resolved_energy: arb,
                    minimum_density: arb, cutoff: int) -> arb:
    """Upper bound for squared unweighted Fourier-antidivergence tail /rho.

    Parseval: total upper minus resolved LOWER; no unproved numerical clipping.
    """
    if cutoff < 0 or not minimum_density > 0:
        raise ValueError("nonnegative cutoff and positive density bound required")
    remainder = total_norm_squared.upper()-resolved_energy.lower()
    upper = nonnegative_upper(remainder)
    return upper/(minimum_density.lower()*(2*arb.pi()*(cutoff+1))**2)


def ball_receipt(value):
    """Round-trippable decimal ball text; never substitute midpoint output."""
    if isinstance(value, acb):
        return {"real": value.real.str(40, more=True),
                "imag": value.imag.str(40, more=True)}
    return value.str(40, more=True)


def midpoint_radius(value: arb):
    """Enclose an Arb ball by an exact float64 midpoint and upward float radius.

    The returned pair is meaningful as [mid-radius,mid+radius] in exact real
    arithmetic; a consumer must additionally round any endpoint arithmetic out.
    """
    midpoint=float(value.mid())
    # arb(float) represents the exact binary64 value, not its printed decimal.
    radius=nextafter(float((value-arb(midpoint)).abs_upper()),inf)
    return midpoint,radius


def longdouble_midpoint_radius(value: arb):
    """Outward NumPy long-double center/radius, accounting for decimal conversion."""
    import numpy as np
    with ctx.workprec(max(ctx.prec, 128)):
        midpoint=np.longdouble(value.mid().str(50,radius=False,more=True))
        numerator,denominator=midpoint.as_integer_ratio()
        exact_mid=arb(numerator)/denominator
        bound=(value-exact_mid).abs_upper()
        radius=np.longdouble(bound.str(50,radius=False,more=True))
        radius=np.nextafter(radius,np.longdouble("inf"))
        # Verify conversion itself by exact integer-ratio arithmetic. More than
        # one step is exceedingly rare but this loop makes that irrelevant.
        while True:
            rn,rd=radius.as_integer_ratio()
            if arb(rn)/rd >= bound:
                break
            radius=np.nextafter(radius,np.longdouble("inf"))
        return midpoint,radius
