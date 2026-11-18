# Executed residual certificate: derivation and arithmetic contract

This note specifies the validation performed by `validated_residual.py` and
`run_validated_continuum.py`. A final positive comparison additionally requires
all node receipts and the independent parameter-interval replay. A software
implementation alone is not that comparison.

## Exact trial, rather than an allegedly exact numerical solution

Let \(\phi_k\) be the sum, over distinct coordinate permutations of a sorted
triple \(k\), of \(\prod_{j=1}^3\cos(2\pi k_jx_j)\). Its complex Fourier
coefficients are exactly \(2^{-\#\{j:k_j\ne0\}}\) on the signed permutation
orbit, and zero elsewhere. Thus storing dyadic coefficients in this basis
defines an exact smooth periodic spatial function. The trial on every time
panel is a polynomial with integer coefficients divided by a power of two.
Initial values and every junction are checked by integer arithmetic.

The approximate Galerkin solve merely proposes these coefficients. Its
floating-point integration, quadrature, and coefficient fitting do not need
separate accuracy assumptions: the residual of the actual stored trial
includes all of them. The survival enclosure uses the corrected observable
and dual-residual theorem in [continuum-residual.md](continuum-residual.md).

## The actual strong residual

Write \(w=e^{-U}\), \(g=(1-r^2/a^2)_+^4\),
\(j=g\min(1,e^{-U})\), and \(d=(1-r^2/a^2)_+e^{-U}\).
For \(R_i=\rho_i(V_i'+AV_i)\), remove the constant chemical mass \(c_i\).
Then, with \(\mathcal D=x\cdot\nabla\),

\[
R_{h,b,0}/c_{h,b}=P_{h,b,0}+gQ_{h,b,0}+jJ_{h,b,0},
\]
\[
R_{h,b,1}/c_{h,b}=wP_{h,b,1}+wgQ_{h,b,1}+jJ_{h,b,1}
 +\frac{4DE}{a^2}d\,\mathcal D V_{h,b,1}.
\]

Here \(P=V'-D\Delta V+A_{\rm thermal}V\),
\(Q=A_{\rm repair}/g\,V\), and \(J\) is the passive and paid switch
difference with its chemical rate factor. The same radial profile \(j\)
occurs in both conformation equations because
\(e^{-U}g\min(1,e^U)=g\min(1,e^{-U})\). The killed thermal matrix is
\(\left(\begin{smallmatrix}3/2&-1/2\\-2&2\end{smallmatrix}\right)\).
The implementation explicitly retains the exit-to-failure diagonal.

## Rigorous, deliberately conservative negative norm

Put \(F_i=R_i/c_i\). Let \(\widehat F_i(k)\) denote its Fourier coefficient
and let \(M\) be the complete coordinate cutoff. Parseval gives

\[
\sum_{\|k\|_\infty>M}\frac{|\widehat F_i(k)|^2}{4\pi^2|k|^2}
\le\frac{\|F_i\|_2^2-
 \sum_{\|k\|_\infty\le M}|\widehat F_i(k)|^2}
 {4\pi^2(M+1)^2}.
\]

The numerator uses an **upper** bound for the total squared norm and a
**lower** bound for every resolved square. A negative upper bound triggers
failure; it is not clipped into a certificate. With \(m_0=1,m_1=e^2\),
\(1/e^{-U_q}\le m_q\), and thermal coercivity
\(\lambda_0=(7-\sqrt{17})/4\), the implemented bound is

\[
\|A^{-1/2}r\|_\rho^2\le\sum_i c_i m_{q_i}
\left[\frac1D\sum_{0<\|k\|_\infty\le M}
\frac{|\widehat F_i(k)|^2}{4\pi^2|k|^2}
+\frac{\|F_i\|_2^2-\sum_{\|k\|_\infty\le M}|\widehat F_i(k)|^2}
{4D\pi^2(M+1)^2}+\frac{|\widehat F_i(0)|^2}{\lambda_0}\right].
\]

This follows from the gradient/scalar direct-sum Cauchy–Schwarz argument in
the earlier residual proof. No unobserved Fourier mode is discarded.

## Exact time integration without sampling between time nodes

For shifted Legendre polynomials \(L_n(2s-1)\),

\[
\int_0^1s^kL_n(2s-1)\,ds=
\frac{(k!)^2}{(k-n)!(k+n+1)!}\quad(k\ge n),
\]

and the integral is zero for \(k<n\). Transform the time-polynomial residual
to this orthogonal basis using the exact rational weights. On a panel of
length \(h\), its integrated squared dual norm is the sum of the modal squared
norms with positive weights \(h/(2n+1)\). Applying the spatial upper bound to
each mode therefore bounds the entire time integral. For the correction,
common-refine the knots of \(U(T-t)\) and \(V(t)\) exactly, transform both,
and sum their matching modal pairings with the same weights.

## Spatial integrals and roundoff

The radial-transform proof and the split at \(U=0\) are documented in
[validated-spatial.md](validated-spatial.md). If
\(H(s)=\int p(|x|)e^{-2\pi i k\cdot x}dx\), \(s=|k|^2\),
the orbit sums for ordinary/ordinary, ordinary/\(\mathcal D\), and
\(\mathcal D/\mathcal D\) Gram entries use, respectively,
\(H(|k+l|^2)\), \(2l\cdot(k+l)H'\), and
\(2(k\cdot l)H'+4[k\cdot(k+l)][l\cdot(k+l)]H''\).
The Fourier coefficient of the drift term uses
\(-2k\cdot(n-k)H'(|n-k|^2)\). Independent direct complex-moment checks
test these signs and orbit multiplicities.

Arb integrals are converted to binary64 midpoint/radius pairs with the
conversion error included. `validated_linear.py` encloses every subsequent
sum, product, and BLAS contraction. Its contract assumes IEEE 754 binary64,
round-to-nearest, gradual underflow, and finite inputs. The standard bound
\(\gamma_n=nu/(1-nu)\), \(u=2^{-53}\), is used with \(4n+4\) in dot
products, allowing reassociation and fused multiply-add. All nonnegative
error sums are inflated again for their own rounding. Overflow/nonfinite
balls abort. In nonnegative matrix products used **only for error upper
bounds**, a positive operand smaller than the exact binary64 threshold
`1e-100` is replaced by that threshold. Monotonicity of the nonnegative dot
product proves that this only raises its upper bound; no tiny contribution
is discarded. This avoids slow gradual-underflow products and does not rely
on flushing subnormals to zero. Tiny absolute budgets cover underflow. This is an explicit
hardware arithmetic assumption, not a claim of a machine-independent formal
proof. Matrix caches are source-bound and checksummed, and can be rebuilt.

## Final comparison and replay

The final scalar enclosure is \(C\pm\sqrt{E_UE_V}\), with the correction
also enclosed. Perform this for the spectator and all 38 fixed-intensity
nodes, include exact-node shift errors, and run the independent rational
Bernstein checker over the full unit interval. All raw receipts, not the
sampled diagnostic values, determine whether a strict sign has been proved.

The completed run is retained in `results/continuum-validated-comparison.json`.
The primary certificate encloses all 39 continuum observables and proves a
positive difference over the full constant-intensity interval. The independent
checker replays the parameter arithmetic, exact input archives, source/cache
hash chain and scalar interval assembly. It does not independently recompute
every spatial contraction. Floating component diagnostics in individual node
receipts aid inspection; the exact stored energy upper endpoints, corrected
observable balls and residual-product bounds carry the certificate.
