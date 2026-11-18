"""Exact time trials and validated constant-intensity interpolation.

No ODE accuracy claim is needed for a stored trial. Its exact polynomial
derivative, including coefficient quantization, enters the continuum residual.
Arb is imported only by the node-cover routines (optional python-flint package).
"""
from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial
import json

import numpy as np


def _arb_fraction(x):
    from flint import arb
    q = Fraction(x)
    return arb(q.numerator) / q.denominator


def arb_endpoint_fraction(x, upper=True):
    """Exact dyadic outward endpoint, without a binary64 conversion."""
    mantissa, exponent = (x.upper() if upper else x.lower()).man_exp()
    return Fraction(int(mantissa)) * Fraction(2) ** int(exponent)


def fraction_record(x):
    q = Fraction(x)
    return {"numerator": str(q.numerator), "denominator": str(q.denominator)}


@dataclass
class DyadicTimeCurve:
    """Integer power coefficients on local s=(t-t_j)/(t_{j+1}-t_j).

    Coefficient values are integers times 2**(-bits). Knots are exact rational
    dyadics. The first two axes are segment and power; remaining axes describe
    arbitrary spatial coefficient arrays. Object integer arrays avoid overflow.
    """
    knots: tuple
    integers: np.ndarray
    bits: int

    def __post_init__(self):
        self.knots = tuple(Fraction(t) for t in self.knots)
        self.integers = np.asarray(self.integers, dtype=object)
        if self.bits < 0 or len(self.knots) != self.integers.shape[0] + 1:
            raise ValueError("inconsistent curve dimensions")
        if any(t.denominator & (t.denominator - 1) for t in self.knots):
            raise ValueError("time knots must be dyadic")
        if any(b <= a for a, b in zip(self.knots, self.knots[1:])):
            raise ValueError("knots must strictly increase")
        if any(not isinstance(v, (int, np.integer)) for v in self.integers.flat):
            raise ValueError("coefficients must be exact integers")
        for j in range(len(self.knots) - 2):
            if not np.array_equal(self.integers[j].sum(axis=0), self.integers[j+1, 0]):
                raise ValueError("trial is discontinuous at a time knot")

    @classmethod
    def quantize(cls, knots, coefficients, bits=40, initial=None):
        """Quantize approximate local polynomials and enforce exact continuity.

        Each segment's constant is copied from the previous exact endpoint.
        The highest coefficient is adjusted so the independently rounded target
        endpoint is retained. This is a definition of a trial, not an error bound.
        """
        raw = np.asarray(coefficients, dtype=float)
        if raw.shape[1] < 2:
            raise ValueError("at least linear local polynomials required")
        scale = 2 ** bits
        ints = np.array([int(round(v * scale)) for v in raw.flat], dtype=object).reshape(raw.shape)
        targets = np.array([int(round(v * scale)) for v in raw.sum(axis=1).flat], dtype=object).reshape(raw[:, 0].shape)
        if initial is not None:
            init = np.asarray(initial, dtype=object)
            if init.shape != raw.shape[2:]:
                raise ValueError("initial value has wrong shape")
            exact = [Fraction(v) * scale for v in init.flat]
            if any(v.denominator != 1 for v in exact):
                raise ValueError("initial value is not represented at chosen dyadic precision")
            ints[0, 0] = np.array([int(v) for v in exact], dtype=object).reshape(init.shape)
        for j in range(raw.shape[0]):
            if j:
                ints[j, 0] = targets[j-1]
            ints[j, -1] = targets[j] - ints[j, :-1].sum(axis=0)
        return cls(tuple(knots), ints, bits)

    def coefficients(self, segment):
        return self.integers[segment] * Fraction(1, 2**self.bits)

    def derivative_coefficients(self, segment):
        h = self.knots[segment+1] - self.knots[segment]
        a = self.coefficients(segment)
        return np.stack([(j / h) * a[j] for j in range(1, len(a))])

    def endpoint(self, final=False):
        a = self.integers[-1].sum(axis=0) if final else self.integers[0, 0]
        return a * Fraction(1, 2**self.bits)

    def reversed(self):
        """Exact U(T-t), requiring the initial knot to be zero."""
        if self.knots[0] != 0:
            raise ValueError("time reversal requires initial time zero")
        a = self.integers[::-1]
        rev = np.empty_like(a)
        for j in range(a.shape[1]):
            rev[:, j] = (-1)**j * sum(comb(k, j)*a[:, k] for k in range(j, a.shape[1]))
        end = self.knots[-1]
        return DyadicTimeCurve(tuple(end-t for t in self.knots[::-1]), rev, self.bits)

    def refine(self, knots):
        """Exact restriction to a finer dyadic partition; no resampling."""
        knots = tuple(Fraction(t) for t in knots)
        if knots[0] != self.knots[0] or knots[-1] != self.knots[-1] or not set(self.knots).issubset(knots):
            raise ValueError("new partition must refine the existing partition")
        blocks=[]
        source=0
        for left,right in zip(knots,knots[1:]):
            while self.knots[source+1] <= left:
                source+=1
            duration=self.knots[source+1]-self.knots[source]
            offset=(left-self.knots[source])/duration
            scale=(right-left)/duration
            old=self.coefficients(source)
            block=np.stack([sum(comb(k,j)*offset**(k-j)*scale**j*old[k]
                                for k in range(j,len(old))) for j in range(len(old))])
            blocks.append(block)
        values=np.stack(blocks)
        denominators=[Fraction(v).denominator for v in values.flat]
        if any(d & (d-1) for d in denominators):
            raise ValueError("refinement ratios must preserve dyadic coefficients")
        bits=max(d.bit_length()-1 for d in denominators)
        integers=np.array([int(v*2**bits) for v in values.flat],dtype=object).reshape(values.shape)
        return DyadicTimeCurve(knots,integers,bits)

    def write(self, path):
        """Portable, exact archive; integer strings do not lose precision."""
        payload = {"schema": "dyadic-time-curve-v1", "bits": self.bits,
                   "knots": [fraction_record(t) for t in self.knots],
                   "shape": list(self.integers.shape),
                   "integers": [str(int(v)) for v in self.integers.flat]}
        with open(path, "w", encoding="utf-8") as out:
            json.dump(payload, out, separators=(",", ":"))

    @classmethod
    def read(cls, path):
        with open(path, encoding="utf-8") as source:
            data = json.load(source)
        if data["schema"] != "dyadic-time-curve-v1":
            raise ValueError("unsupported trial format")
        knots = [Fraction(int(t["numerator"]), int(t["denominator"])) for t in data["knots"]]
        integers = np.array([int(v) for v in data["integers"]], dtype=object).reshape(data["shape"])
        return cls(tuple(knots), integers, int(data["bits"]))


def integrate_pair(left, right, duration, pairing):
    """Integrate local polynomial bilinear pairing, with exact time weights.

    `pairing` receives two exact coefficient arrays and must return a rigorous
    Arb enclosure (or an exact rational). Spatial validation is its responsibility.
    """
    value = 0
    for p, u in enumerate(left):
        for q, v in enumerate(right):
            weight = Fraction(duration) / (p + q + 1)
            term = pairing(u, v)
            value += term * (_arb_fraction(weight) if hasattr(term, "rad") else weight)
    return value


def integrate_curves(left, right, pairing, right_derivative=False):
    if left.knots != right.knots:
        knots=tuple(sorted(set(left.knots+right.knots)))
        left,right=left.refine(knots),right.refine(knots)
    return sum(integrate_pair(left.coefficients(j),
                              right.derivative_coefficients(j) if right_derivative else right.coefficients(j),
                              b-a, pairing)
               for j, (a, b) in enumerate(zip(left.knots, left.knots[1:])))


def shifted_legendre_coefficients(power):
    """Exact conversion: int_0^1 P_j(2s-1)P_k(2s-1)=delta_jk/(2k+1).

    This eliminates temporal off-diagonal cancellation when computing residual
    quadratic forms. The source operator commutes with this finite transform.
    """
    power=np.asarray(power,dtype=object)
    n=len(power)
    return np.stack([sum(Fraction((2*k+1)*factorial(p)**2,
                                factorial(p-k)*factorial(p+k+1))*power[p]
                         for p in range(k,n)) for k in range(n)])


def residual_legendre_terms(curve):
    """Yield (exact time weight, U_k, U'_k) on each panel.

    Any squared spatial Hilbert norm of U'+AU integrates to the sum of these
    positive weights times the squared norm of U'_k+AU_k. No time quadrature
    or interval evaluation of a high-degree monomial is involved.
    """
    for panel,(left,right) in enumerate(zip(curve.knots,curve.knots[1:])):
        coefficients=curve.coefficients(panel)
        derivative=np.zeros_like(coefficients)
        derivative[:-1]=curve.derivative_coefficients(panel)
        values=shifted_legendre_coefficients(coefficients)
        derivatives=shifted_legendre_coefficients(derivative)
        for k,(value,deriv) in enumerate(zip(values,derivatives)):
            yield (right-left)/ (2*k+1), value, deriv


def corrected_observable(primal, dual, initial_payoff, mass_pairing, operator_pairing):
    """C=<f,V(T)>-int <U(T-t),V'+AV>; all rounding is enclosed."""
    rev = primal.reversed()
    endpoint = mass_pairing(initial_payoff, dual.endpoint(final=True))
    return endpoint - integrate_curves(rev, dual, mass_pairing, True) - integrate_curves(rev, dual, operator_pairing)


def chebyshev_node_plan(count=38, bits=60, precision=192):
    """Enclose exact roots and choose exact dyadic model parameters.

    An actual node survival enclosure must be enlarged by `root_shift_bound`:
    ||f||||ell|| T||B|| <=700*2*16.5=23100. This handles the distinction
    between transcendental interpolation nodes and dyadic model parameters.
    """
    from flint import arb, ctx
    with ctx.workprec(precision):
        nodes = []
        for j in range(count):
            alpha = (1 + (arb(2*j+1)/(2*count)).cos_pi()) / 2
            numerator = int(round(float(alpha) * 2**bits))
            dyadic = arb(numerator) / 2**bits
            distance = abs(alpha-dyadic)
            nodes.append({"index": j, "alpha_numerator": str(numerator), "alpha_bits": bits,
                          "root_lower": fraction_record(arb_endpoint_fraction(alpha, False)),
                          "root_upper": fraction_record(arb_endpoint_fraction(alpha)),
                          "root_shift_bound": fraction_record(arb_endpoint_fraction(23100*distance))})
        remainder = Fraction(1400) * Fraction(33, 4)**count / factorial(count)
        return {"schema": "chebyshev-root-cover-plan-v1", "count": count,
                "precision_bits": precision, "nodes": nodes,
                "derivative_interpolation_remainder": fraction_record(remainder),
                "node_error_amplification_squared": 2*count-1}


def node_interpolant(records, count=38, precision=256, coefficient_bits=160):
    """Convert dyadic node centers to an exact dyadic power polynomial.

    Returned conversion error bounds the sup norm on [0,1]. Data uncertainties
    are kept separate and amplified only by sqrt(2n-1), avoiding cancellation
    blow-up from uncertain nodal values in the monomial basis.
    """
    from flint import arb, ctx
    if len(records) != count:
        raise ValueError("all interpolation nodes are required")
    plan = chebyshev_node_plan(count, precision=precision)
    with ctx.workprec(precision):
        ys = [arb(int(v["center_numerator"])) / 2**int(v["center_bits"]) for v in records]
        errors = []
        for record, node in zip(records, plan["nodes"]):
            if (str(record["alpha_numerator"]), int(record["alpha_bits"])) != (node["alpha_numerator"], node["alpha_bits"]):
                raise ValueError("node parameter does not match the verified node plan")
            err = Fraction(str(record["radius_upper"]))
            if err < 0:
                raise ValueError("negative node radius")
            shift = node["root_shift_bound"]
            errors.append(_arb_fraction(err + Fraction(int(shift["numerator"]), int(shift["denominator"]))))
        coeffs = []
        for k in range(count):
            c = sum(y*(arb(k*(2*j+1))/(2*count)).cos_pi() for j, y in enumerate(ys)) * (arb(1 if k == 0 else 2)/count)
            coeffs.append(c)
        # Chebyshev recurrence T_k(2a-1), all polynomial arithmetic in Arb.
        basis0, basis1 = [arb(1)], [arb(-1), arb(2)]
        power = [arb(0) for _ in range(count)]
        for k, c in enumerate(coeffs):
            basis = basis0 if k == 0 else basis1
            for p, value in enumerate(basis):
                power[p] += c*value
            if k >= 1 and k+1 < count:
                nxt = [arb(0) for _ in range(len(basis1)+1)]
                for p, value in enumerate(basis1):
                    nxt[p] -= 2*value
                    nxt[p+1] += 4*value
                for p, value in enumerate(basis0):
                    nxt[p] -= value
                basis0, basis1 = basis1, nxt
        midpoint_records, conversion_error = [], arb(0)
        for c in power:
            # Round an exact dyadic midpoint with integer arithmetic.
            midpoint = arb_endpoint_fraction(c.mid()) * 2**coefficient_bits
            numerator = (2*midpoint.numerator + midpoint.denominator)//(2*midpoint.denominator)
            q = arb(numerator)/2**coefficient_bits
            conversion_error += abs(c-q)
            midpoint_records.append({"numerator": str(numerator), "bits": coefficient_bits})
        max_error = max(errors, key=lambda x: arb_endpoint_fraction(x))
        amplification = arb(2*count-1).sqrt()
        remainder = Fraction(1400)*Fraction(33,4)**count/factorial(count)
        return {"schema": "chebyshev-root-cover-polynomial-v1", "count": count,
                "power_coefficients": midpoint_records,
                "coefficient_error_upper": str(arb_endpoint_fraction(conversion_error)),
                "node_uncertainty_upper": str(arb_endpoint_fraction(amplification*max_error)),
                "interpolation_remainder_upper": str(remainder),
                "precision_bits": precision}
