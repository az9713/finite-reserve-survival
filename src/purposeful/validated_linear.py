"""Binary64 midpoint/radius linear algebra with explicit rounding budgets.

The arithmetic model is IEEE 754 round-to-nearest, gradual underflow, and
finite operands. BLAS dot products may reorder additions or use FMA: the
deliberately conservative gamma(4*n) bound covers either implementation.
This module never treats an unadorned numerical dot product as an enclosure.
"""
from dataclasses import dataclass
import numpy as np

U = np.finfo(float).eps / 2
TINY = np.finfo(float).tiny


def gamma(n):
    nu = n * U
    if nu >= .01:
        raise ArithmeticError("rounding budget too large")
    return np.nextafter(nu / (1-nu), np.inf)


def positive_sum(a, axis=None):
    a = np.asarray(a)
    if np.any(a < 0):
        raise ValueError("nonnegative summands required")
    n = a.size if axis is None else a.shape[axis]
    ans = np.sum(a, axis=axis)
    return np.nextafter((ans + n*TINY) / (1-gamma(4*n+4)), np.inf)


def positive_matmul(a, b):
    if np.any(a<0) or np.any(b<0):
        raise ValueError('nonnegative factors required')
    # This routine computes an UPPER bound only. Raising tiny positive
    # operands is exact monotonicity, not a flush-to-zero assumption. It
    # avoids extremely slow subnormal products in radius propagation.
    a=np.where((a>0)&(a<1.e-100),1.e-100,a)
    b=np.where((b>0)&(b<1.e-100),1.e-100,b)
    n = a.shape[-1]
    ans = a @ b
    return np.nextafter((ans+n*TINY)/(1-gamma(4*n+4)), np.inf)


@dataclass
class BallArray:
    center: np.ndarray
    radius: np.ndarray

    def __post_init__(self):
        self.center = np.asarray(self.center, dtype=float)
        self.radius = np.broadcast_to(np.asarray(self.radius, dtype=float), self.center.shape)
        if not np.all(np.isfinite(self.center)) or not np.all(np.isfinite(self.radius)) or np.any(self.radius < 0):
            raise ArithmeticError("nonfinite or invalid ball")

    @classmethod
    def exact(cls, a):
        """Declare *binary64 values* exact, not their intended decimal values."""
        return cls(np.asarray(a, dtype=float), 0.)

    def __getitem__(self, key):
        return BallArray(self.center[key], self.radius[key])

    def __add__(self, other):
        if not isinstance(other, BallArray):
            other = BallArray.exact(other)
        c = self.center + other.center
        r = self.radius + other.radius + gamma(8)*(np.abs(self.center)+np.abs(other.center))+8*TINY
        return BallArray(c, np.nextafter(r/(1-gamma(8)), np.inf))

    __radd__ = __add__

    def __neg__(self):
        return BallArray(-self.center, self.radius)

    def __sub__(self, other):
        return self + (-other if isinstance(other, BallArray) else -np.asarray(other))

    def __mul__(self, other):
        if not isinstance(other, BallArray):
            other = BallArray.exact(other)
        c = self.center * other.center
        r = (np.abs(self.center)*other.radius + self.radius*(np.abs(other.center)+other.radius)
             + gamma(12)*np.abs(c)+12*TINY)
        return BallArray(c, np.nextafter(r/(1-gamma(12)), np.inf))

    __rmul__ = __mul__

    def sum(self, axis=None):
        c = self.center.sum(axis=axis)
        n = self.center.size if axis is None else self.center.shape[axis]
        r = positive_sum(self.radius + gamma(4*n+4)*np.abs(self.center), axis=axis)
        return BallArray(c, r)

    def matmul(self, other):
        if not isinstance(other, BallArray):
            other = BallArray.exact(other)
        c = self.center @ other.center
        a, b = np.abs(self.center), np.abs(other.center)
        rounding = positive_matmul(a, b)*gamma(4*self.center.shape[-1]+4)
        r = (positive_matmul(a, other.radius)
             + positive_matmul(self.radius, np.nextafter(b+other.radius, np.inf)) + rounding)
        return BallArray(c, np.nextafter(r/(1-gamma(8)), np.inf))

    @property
    def T(self):
        return BallArray(self.center.T, self.radius.T)

    def lower(self):
        return np.nextafter(self.center-self.radius, -np.inf)

    def upper(self):
        return np.nextafter(self.center+self.radius, np.inf)

    def squared(self):
        c = self.center*self.center
        r = 2*np.abs(self.center)*self.radius + self.radius*self.radius + gamma(16)*np.abs(c)+16*TINY
        # For positive sums of squares, a nonnegative lower endpoint is valid.
        return BallArray(c, np.nextafter(r/(1-gamma(16)), np.inf))


def rational_ball(q):
    from fractions import Fraction
    from flint import arb, ctx
    from .validated_spatial import midpoint_radius
    q = Fraction(q)
    with ctx.workprec(120):
        c,r = midpoint_radius(arb(q.numerator)/q.denominator)
    return BallArray(c,r)


def quadratic(x, matrix, y=None):
    """Enclose x M y, optionally for a batch of row vectors."""
    if y is None:
        y = x
    return (x.matmul(matrix)*y).sum(axis=-1)
