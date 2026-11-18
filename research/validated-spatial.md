# Validated spatial integrals for the Brownian residual certificate

This note describes the executable spatial integration oracle in
`src/purposeful/validated_spatial.py`. It encloses integrals for the actual
continuum coefficients. It is one component of the residual proof; generating
these tables alone is not a survival certificate.

## Exact model declarations and analytic branches

The fundamental cube is \([-1/2,1/2]^3\), the contact radius is the exact rational
\(a=8/25\), and inside the contact ball

\[
c(r)=1-r^2/a^2,\qquad U(r)=2-6c(r)^2,\qquad g(r)=c(r)^4.
\]

Outside the ball, \(U=2\) and \(g=0\). The Metropolis acceptance changes its
analytic expression at

\[
r_*=a\sqrt{1-\sqrt{1/3}}.
\]

A `RadialProfile` stores separate finite sums of
\(C c^p\exp(eU)\), with rational \(C,e\) and nonnegative integer \(p\), on the
two intervals \([0,r_*]\) and \([r_*,a]\). It also stores a constant exterior
as a sum of \(C\exp(2e)\). Products are exact symbolic products followed by
collection of equal \((p,e)\) pairs. Input floats are rejected for these
physical declarations. A finite trial function's floating coefficients are
instead interpreted as the exact dyadic numbers they represent.

The six source multipliers are

| Name | Interior multiplier | Exterior |
| --- | --- | --- |
| `one` | \(1\) | \(1\) |
| `contact` | \(g\) | \(0\) |
| `switch` | \(g\min(1,e^{-U})\) | \(0\) |
| `weight` | \(e^{-U}\) | \(e^{-2}\) |
| `weighted_contact` | \(g e^{-U}\) | \(0\) |
| `drift` | \(c e^{-U}\) | \(0\) |

The coefficient \(4DE/a^2\), the chemical equilibrium constants, and the
chemical rates multiply these profiles outside the oracle. Thus the raw
weighted residual has a sum of profile times Fourier polynomial terms, and
one term of the form \(c e^{-U}x\cdot\nabla P\). Squared residuals require
Cartesian moments of degree at most two.

## Radial Fourier reduction

Let \(p_\infty\) denote the exterior constant. Its compact correction is
\(f(r)=p(r)-p_\infty\) in the ball and zero outside. Spherical integration gives

\[
\widehat p(k)=p_\infty\mathbf1_{k=0}+F(|k|^2),\qquad
F(n)=4\pi\int_0^a r^2 f(r)
\operatorname{sinc}(2\pi r\sqrt n)\,dr.
\]

This identity is exact because the ball lies entirely within the fundamental
cube. The oracle integrates on each analytic branch separately. In code each
branch is reparameterized to exact endpoints \(t=0,1\); the algebraic endpoints
are themselves Arb balls. This avoids treating an approximate decimal branch
location as exact.

For Cartesian moments the integrator uses entire kernels rather than division
by powers of \(r\) near the origin. With \(z=-\pi^2r^2\),

\[
\frac{d^j}{dn^j}\operatorname{sinc}(2\pi r\sqrt n)
=\frac{z^j}{(3/2)_j}\,{}_0F_1(3/2+j;nz),\qquad j=0,1,2.
\]

The hypergeometric identity proves the derivative formulas, but direct
hypergeometric ball evaluation at large frequency produces needlessly wide
enclosures. The implementation uses the following equivalent entire
trigonometric kernels, with \(t=2\pi r\sqrt n\):

\[
K_1=\frac{\cos t-\operatorname{sinc}t}{2n},\qquad
K_2=\frac{3\operatorname{sinc}t-3\cos t-t\sin t}{4n^2},\quad n>0.
\]

At \(n=0\), their exact limits are \(-2\pi^2r^2/3\) and
\(4\pi^4r^4/15\). There is no division by the integration variable, so no
analyticity check is lost at \(r=0\). For \(j=0\), Arb's entire `sinc`
evaluation is used directly.
For compact profiles, differentiation under the finite integral gives

\[
\int x_i p(|x|)e^{-2\pi i k\cdot x}\,dx
=\frac{i k_i}{\pi}F'(|k|^2),
\]

\[
\int x_i x_j p(|x|)e^{-2\pi i k\cdot x}\,dx
=-\frac{2\delta_{ij}F'(|k|^2)+4k_i k_jF''(|k|^2)}{4\pi^2}.
\]

The formulas include \(k=0\) without a limiting division. For a nonzero
exterior constant the general `moment` API additionally supplies the exact
separable cube moments. Positive-derivative shell tables reject noncompact
profiles, preventing accidental use of a ball-moment formula for a cube-wide
constant.

## What the numerical enclosure actually guarantees

The engine is `python-flint` 0.9.0, which exposes Arb's validated complex
integration. Each branch integrand is a polynomial, exponential, or entire
trigonometric expression in the complex integration variable. Consequently
it is analytic on every evaluation ball; there is no unchecked square-root,
absolute-value, or `min` branch in the callback. The acceptance kink has already
been split symbolically. A nonfinite integral enclosure raises an exception.

The requested integration tolerance is only a work target. The mathematical
output is the returned ball, including its radius. If the integrator returns
a wider finite ball, that wider radius is retained rather than replaced by
the requested tolerance.

For efficient downstream matrix contractions, `SpatialOracle.table` exports
binary64 midpoint/radius arrays. If the Arb result is \(I\), the exported
midpoint \(m\) is treated as an exact dyadic number and the radius is rounded
up from

\[
\sup\{|x-m|:x\in I\}.
\]

Thus \(I\subseteq[m-r,m+r]\) in exact real arithmetic. Consumers must still
enclose their own summation and multiplication errors. A BLAS or FFT result
without an accompanying rounding bound does not become certified merely
because its coefficient table is enclosed.

## Fourier tails

For a real source \(R\), let \(E\) enclose \(\int R^2\), and let \(E_M\) enclose
the resolved Fourier energy through the full cube \(\|k\|_\infty\le M\), with
all multiplicities included. Parseval and \(|k|\ge M+1\) outside the cube imply

\[
\int\frac{|Z-Z^{(M)}|^2}{\rho}
\le\frac{E^{\rm upper}-E_M^{\rm lower}}
{\rho_{\min}^{\rm lower}[2\pi(M+1)]^2}.
\]

`tail_flux_bound` implements this upper-minus-lower formula and rejects an
inconsistent negative upper endpoint. It does not silently clip an
underestimated squared norm to zero. The resolved and tail fluxes must then
be combined with the triangle inequality in the weighted norm.

## Reproducible checks and receipts

Run the analytic checks with

```powershell
$env:PYTHONPATH='src'
python -m unittest discover -s tests -p test_validated_spatial.py -v
```

The checks compare the enclosures with independently written exact formulas:
the Fourier transform of a ball indicator, its first Cartesian moment, the
zero-mode beta integral of \(c^4\), isotropic second moments, and the cube
moments of a constant. Further checks cover conjugation, Metropolis profile
products, outward export, and rejection of inconsistent Parseval data. These
checks exercise implementation details; the derivation above justifies the
continuum identities.

Build the reusable coefficient package with

```powershell
python scripts/build_validated_spatial_tables.py --cutoff 12 --resolved 24
```

The NPZ retains every midpoint and radius. The adjacent JSON records all
profile declarations, precision, actual maximum radii, source SHA-256 values,
and the NPZ SHA-256. The retained scope explicitly says that temporal
residual contractions and tail bounds are additional requirements.

The retained \(K=12,M=24\) package contains 36 tables, 82,548 enclosed
coefficients and 137,432 validated branch integrals. Its largest exported
radius is \(4.4053599786\times10^{-16}\). These are observed certificate
properties, not a replacement for using each individual stored radius.
The first implementation's hypergeometric derivative evaluation returned
valid but unusably wide balls at high frequency; the trigonometric rewrite
above resolved that numerical conditioning problem. The builder now rejects
any coefficient radius exceeding \(10^{-10}\), preventing such a package
from being silently used in a residual calculation.

`tests/test_validated_residual.py` independently expands the symmetry orbits
into complex exponentials and evaluates their Cartesian moments with this
oracle. It checks the reduced transform matrices, mixed and drift-squared
Gram matrices, a dyadic contraction, and a separately enumerated reaction
edge ledger for both spectator control and fixed intensity \(3/7\). The
retained test transcripts are
`results/continuum-validated-spatial-tests.txt` and
`results/continuum-validated-residual-tests.txt`.
