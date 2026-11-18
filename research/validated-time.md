# Exact time reconstruction and the complete fixed-intensity comparison

This note records the temporal and parameter parts of the continuum certificate.
The spatial residual theorem is in [continuum-residual.md](continuum-residual.md).
An exact trial archive alone is not a continuum survival certificate: the
continuum residual and the observable correction must still be enclosed.

## The stored trial is the mathematical object being validated

For each dyadic interval \([t_j,t_{j+1}]\), set
\(s=(t-t_j)/h_j\), \(h_j=t_{j+1}-t_j\), and store

\[
U(t,x)=\sum_{p=0}^{d}2^{-B}a_{jp}(x)s^p,
\qquad
V(t,x)=\sum_{p=0}^{d}2^{-B}b_{jp}(x)s^p.
\]

Every scalar spatial coefficient of \(a,b\) is an integer. The spatial basis is
the sum over distinct coordinate permutations of products
\(\prod_{i=1}^3\cos(2\pi k_i x_i)\), with sorted nonnegative frequency triples.
This basis has a rational finite complex Fourier expansion; there are no hidden
square-root orbit normalizations in the stored definition.

The generator uses an approximate weighted Galerkin evolution only to choose
these integers. It rounds the initial constant coefficients to the exact prescribed
values, copies each previous exact endpoint into the next constant coefficient,
and adjusts the last coefficient to retain a rounded target endpoint. Hence
\(U(0)=f\), \(V(0)=\ell\), and both curves are continuous by exact integer equality.
For this model the unique nonzero launch coefficient is \(59049/4\), which is
dyadic. The derivative of the stored trial is exactly

\[
U'(t,x)=h_j^{-1}\sum_{p=1}^{d}p2^{-B}a_{jp}(x)s^{p-1}.
\]

No assumption that the numerical ODE solver is exact enters the proof. Numerical
propagation, fitting and coefficient rounding errors are part of this exact
trial's residual \(r_U=U'+AU\), and likewise for \(V\).

The current archives use degree 12, \(B=40\), and the eleven panels whose knots
are \(0,2^{-9},2^{-8},\ldots,2^1\). These panels resolve the initial transient
of the trial; they add no external clock or controller to the physical model.
Each K=12 trial has shape \(11\times13\times20\times455\), and its pair contains
2,602,600 exact integers. The exact initial values, continuity and hashes are
replayed by `verify_curve_archive`, not merely asserted by a label.

## Exact time integration without interval cancellation

Let \(L_k(s)=P_k(2s-1)\) be the shifted Legendre polynomial. For \(p\ge k\),

\[
\int_0^1s^pL_k(s)\,ds
=\frac{(p!)^2}{(p-k)!(p+k+1)!},
\qquad
\int_0^1 L_kL_l=\frac{\delta_{kl}}{2k+1}.
\]

The first identity follows from the Rodrigues formula
\(L_k=(k!)^{-1}d^k/ds^k[s^k(s-1)^k]\), integrating by parts k times and evaluating
the beta integral. Boundary terms vanish because \(s^k(s-1)^k\) has zeros of
order k at both endpoints. Thus a power polynomial \(c(s)=\sum_pc_ps^p\) has
exact Legendre coefficients

\[
\widetilde c_k=(2k+1)\sum_{p=k}^d
\frac{(p!)^2}{(p-k)!(p+k+1)!}c_p.
\]

Because the spatial operator is linear, the residual's coefficients are
\(\widetilde r_k=\widetilde{U'}_k+A\widetilde U_k\). In any spatial Hilbert norm,

\[
\int_{t_j}^{t_{j+1}}\|r(t)\|^2dt
=h_j\sum_{k=0}^d\frac{\|\widetilde r_k\|^2}{2k+1}.
\]

This identity avoids sampling time and avoids cancellation between large power
coefficients in a double sum. Every temporal weight is an exact positive
rational. The implementation encloses the rational transformation and spatial
quadratic forms using midpoint/radius arithmetic, then sums their outward upper
endpoints using Arb and exact rational weights.

## The correction is retained

The corrected scalar is

\[
C=(f,V(T))_\rho-
\int_0^T(U(T-t),V'(t)+AV(t))_\rho\,dt.
\]

Time reversal is exact: the coefficient of \(s^j\) in a reversed local power
polynomial is \((-1)^j\sum_{k\ge j}\binom{k}{j}c_k\). The reversed primal and
dual partitions generally differ. Their union supplies a common refinement;
restriction uses the exact binomial substitution \(s_{old}=a+bs_{new}\), not
interpolation or resampling. The correction is then integrated by the same
Legendre orthogonality identity, this time as a bilinear form. The nonzero
computed correction is subtracted before adding the residual-product error.
It is never discarded on an assumption of exact Galerkin orthogonality.

For practical replay, `validated_refinement.fast_refine` first computes the
dyadic binomial factors as scalar rationals, selects a common power-two
denominator, and performs all array arithmetic with exact Python integers.
It is algebraically identical to the slower Fraction-per-entry reference.
Tests compare both definitions by exact cross multiplication, including
degree-12 restrictions of the actual K12 archive. A full reflected K12 curve
was restricted to 20 panels in under one second; its coefficients use 148
dyadic bits. This is an implementation improvement, not an additional
mathematical approximation.

## Covering every constant paid intensity

Write \(A_\alpha=A_0+\alpha B\), \(0\le\alpha\le1\), with
\(\|B\|\le33/2\), and \(S(\alpha)=(f,e^{-2A_\alpha}\ell)_\rho\).
The bounded-perturbation Dyson formula has n insertions of B, integrated over
an ordered n-simplex of volume \(2^n/n!\). Differentiation contributes n!;
contractivity of each real-parameter semigroup and \(\|f\|\|\ell\|<700\) give

\[
|S^{(n)}(\alpha)|\le700\,33^n.
\]

For the n mapped Chebyshev roots \(\alpha_j=(1+\cos((2j+1)\pi/(2n)))/2\),
the monic nodal product has uniform norm \(2^{1-2n}\). The ordinary Lagrange
remainder therefore bounds degree n-1 interpolation by

\[
\|S-I_{n-1}S\|_\infty
\le \frac{1400(33/4)^n}{n!}.
\]

At n=38 this is below \(1.8\times10^{-7}\). This is a full interval argument,
not a claim that the best fixed intensity occurs at an endpoint.

The exact Chebyshev roots are transcendental. The model is validated at stored
dyadic parameters \(\widehat\alpha_j\). Arb encloses each true root, and the
node radius is enlarged by

\[
|S(\alpha_j)-S(\widehat\alpha_j)|
\le23100|\alpha_j-\widehat\alpha_j|.
\]

The constant is \(700\cdot2\cdot(33/2)\). This small correction is required even
when binary64 rounding makes the displayed parameters look identical.

The interpolant of exact stored dyadic node centers is constructed by a DCT
with Arb-enclosed cosine coefficients, converted to a power polynomial and
rounded to exact dyadic coefficients. The sum of conversion-error coefficient
radii bounds its uniform conversion error on [0,1]. True node uncertainty is
kept separate. A simple independently proved amplification bound is
\(\sqrt{2n-1}=\sqrt{75}\): discrete cosine orthogonality bounds the interpolant's
coefficient Euclidean norm, then Cauchy-Schwarz bounds its value. The independent
checker also computes a sharper certified Lebesgue constant by exact rational
Bernstein subdivision. The final comparison must use one amplification bound,
not count both or absorb uncertainty twice.

## Receipts and limitations

- `scripts/validate_continuum_nodes.py`: generates and replays exact trial archives;
  `--node-plan` emits rigorous Chebyshev roots and dyadic shift budgets.
- `scripts/run_validated_continuum.py`: encloses residual energies, the complete
  correction and the survival interval for an actual archive.
- `results/continuum-time-verification.json`: exact K12 input replays and hashes.
- `tests/test_validated_time.py`: exact polynomial identities, initial values,
  continuity rejection, reversal, refinement, archived round trips, root enclosure
  and interpolation error separation. The integer refinement is checked against
  the independent Fraction implementation.

These temporal and parameter statements do not by themselves establish the
positive survival gap. That conclusion requires successfully executed spatial
enclosures for the spectator and all required fixed nodes, followed by a rigorous
maximum of the resulting fixed-family interpolant and a strict final subtraction.
