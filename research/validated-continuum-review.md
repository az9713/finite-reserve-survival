# Independent review of the validated continuum route

This note records an implementation review, explicit mathematical checks, and the scope of the independent replay. It does not treat a small sampled residual, a matching plot, or passing unit tests as a continuum error enclosure.

## The observable and the residual identity

In the common Gibbs space \(H=L^2(\rho)\), the killed operator \(A=-\mathcal L\) is self-adjoint and

\[
\mathfrak a(w,w)\ge \lambda_0\|w\|_\rho^2
 +D\sum_i\int\rho_i|\nabla w_i|^2,
\qquad \lambda_0=\frac{7-\sqrt{17}}4.
\]

The launch is uniform in position, in the single chemical state \((h,b,q)=(1,4,0)\). Hence \(f=1\) and \(\ell=\mathbf 1_{(1,4,0)}/c_{1,4}\), with \(c_{1,4}=4/59049\). The desired probability is \(S(T)=(f,e^{-TA}\ell)_\rho\), with \(T=2\). A numerical matrix calculation using another launch or a different normalization would not certify this probability.

For trial curves with **exact** initial values, let \(r_U=U'+AU\), \(r_V=V'+AV\), and

\[
C=(f,V(T))_\rho-\int_0^T(U(T-t),r_V(t))_\rho\,dt.
\]

The dual error \(e=e^{-tA}\ell-V\) satisfies \(e'+Ae=-r_V\). Direct differentiation gives

\[
\frac d{dt}(U(T-t),e(t))_\rho
=-(r_U(T-t),e(t))_\rho-(U(T-t),r_V(t))_\rho.
\]

Consequently \(S(T)-C=-\int(r_U(T-t),e(t))_\rho\,dt\). Testing the error equation with \(e\), applying \(2ab\le a^2+b^2\), and integrating yields

\[
\int_0^T\|A^{1/2}e\|_\rho^2\,dt
\le\int_0^T\|A^{-1/2}r_V\|_\rho^2\,dt.
\]

Thus the corrected center has error at most the product of the primal and dual residual norms. The sign of the correction in the implementation agrees with this independent derivation.

Continuous piecewise polynomial trial curves are admissible. Their derivatives are taken on each open panel; continuity makes the integrated identity telescope without jump terms. Discontinuous curves require jump corrections. The stored polynomial coefficients, rather than an unrounded ODE solution, define the trial being certified. The source implementation enforces exact initial values and exact continuity; unit tests must also exercise rejection of discontinuous input.

## Spatial and roundoff checks

For \(R_i=\rho_i r_i\), set \(Z_i=\nabla(-\Delta)^{-1}(R_i-\bar R_i)\). Since \(R_i=\bar R_i-\nabla\cdot Z_i\), a single Cauchy–Schwarz inequality on the scalar and gradient direct sum gives

\[
\|A^{-1/2}r\|_\rho^2\le
\sum_i\left(D^{-1}\int|Z_i|^2/\rho_i
+\lambda_0^{-1}|\bar R_i|^2\int1/\rho_i\right).
\]

For a Fourier cube with cutoff \(M\), Parseval gives the tail bound

\[
\|(Z_i-Z_i^{(M)})/\sqrt{\rho_i}\|_2^2
\le\frac{\int|R_i|^2-\sum_{\|k\|_\infty\le M}|\widehat R_i(k)|^2}
{\inf\rho_i\,[2\pi(M+1)]^2}.
\]

The upper endpoint of the total norm must be used and the lower endpoint of resolved energy must be subtracted. A negative upper endpoint is an error; setting it to zero would hide an invalid enclosure. For a real cosine basis, multiplicities must agree with Parseval normalization. Merely increasing a sampled FFT grid does not verify this tail.

The radial coefficient oracle was inspected for the sign and normalization of its Cartesian moments. For a compact radial transform \(F(|k|^2)\), it uses

\[
\int x_jF_xe^{-2\pi ik\cdot x}\,dx
=\frac{ik_j}{\pi}F'(|k|^2),\quad
\int x_ix_jF_xe^{-2\pi ik\cdot x}\,dx
=-\frac{2\delta_{ij}F'+4k_ik_jF''}{4\pi^2}.
\]

These formulas and the conjugate pairing convention are correct. Splitting the radial integral at the contact radius and the algebraic root \(U(r)=0\) avoids integrating across the Metropolis branch change. The Arb integration output is the enclosure; its requested tolerance is only a work target. Ball endpoints must be retained when crossing into rational arithmetic. A float midpoint alone carries no enclosure claim.

For the raw symmetric basis \(\phi_k=2^{-\#\{k_j\ne0\}}\sum_{v\in\operatorname{orbit}(k)}e^{2\pi iv\cdot x}\), radial invariance reduces one orbit sum to its cardinality. With \(p=k+l\), the reviewed Gram factors are \(2l\cdot p\,F'\) for one \(x\cdot\nabla\), and \(2k\cdot l\,F'+4(k\cdot p)(l\cdot p)F''\) for two. The Fourier transform of \(F_x x\cdot\nabla\phi_k\) has factor \(-2k\cdot(n-k)F'(|n-k|^2)\). All three signs and factors agree with direct differentiation.

The weighted switch profile is shared between the two mode equations because \(e^{-U}m(-U)=m(U)\). The thermal, repair, chemical-equilibrium, and drift prefactors in `validated_residual.py` were inspected against the declared generator. The reciprocal density bound is \(1\) for \(q=0\) and \(e^2\) for \(q=1\). This is conservative but valid; it need not approximate the actual weighted flux norm closely to give a proof.

The binary64 midpoint/radius contractions assume finite IEEE 754 round-to-nearest arithmetic, allowing reordered BLAS sums and fused multiply-add. Their \(\gamma(4n+4)\) budgets exceed the standard sequential dot-product budget, and the explicit normal-underflow terms conservatively dominate subnormal errors. The review requested rejection of a negative tail upper endpoint, an Arb enclosure for the reciprocal Fourier factor, and retention of source/model hashes for cached matrices. These are proof-integrity requirements: silently reusing a cache produced by different code would not be a replay of the stated calculation.

A safe performance improvement for the **positive** products used to propagate radii is to replace each operand \(0<x<c\) by the positive binary64 constant \(c=\operatorname{float64}(10^{-100})\), retaining exact zeros. Coordinatewise monotonicity of nonnegative matrix multiplication proves that the inflated exact product upper-bounds the original product. Applying the same outward rounding budget therefore remains valid. Every nonzero elementary product is then at least \(c^2\), which is normal in binary64 and avoids slow subnormal BLAS arithmetic. The added conservatism is propagated in the computed radius; no universal absolute estimate for its size is assumed. This improvement does not alter center values or the physical model.

The spatial implementation also replaced high-frequency hypergeometric derivative evaluations, whose balls were valid but excessively wide, with the equivalent entire-in-\(r\) kernels. Writing \(t=2\pi r\sqrt n\),

\[
K_1=\frac{\cos t-\operatorname{sinc}t}{2n},\qquad
K_2=\frac{3\operatorname{sinc}t-3\cos t-t\sin t}{4n^2}.
\]

Their limits at \(n=0\) are \(-2\pi^2r^2/3\) and \(4\pi^4r^4/15\), respectively. Independent differentiation confirms these formulas. The variables in the denominators are fixed nonzero frequency indices, so they introduce no singularity in the integration variable.

## A full constant-intensity cover with 38 nodes

The independent checker proves the following parameter coverage theorem. It applies to actual continuum node enclosures when they are available; it cannot turn unvalidated sampled node values into such enclosures.

For \(A_\alpha=A_0+\alpha B\), \(\alpha\in[0,1]\), the paid channel is a bounded nonnegative self-adjoint operator in the same Gibbs space, with \(\|B\|\le33/2\). Duhamel's ordered simplex formula and contraction give

\[
|S^{(n)}(\alpha)|\le\|f\|\|\ell\|(T\|B\|)^n
<700\,33^n.
\]

The norm estimate follows without rounded transcendental constants: \(e^4<55\) is checked by a rational Taylor tail, and

\[
\|f\|^2\|\ell\|^2
\le\frac59(1+e^4)\frac{59049}4
<459270<700^2.
\]

Take \(n=38\) exact Chebyshev roots
\(\alpha_j=(1+\cos((2j+1)\pi/(2n)))/2\).
The nodal polynomial has sup norm \(2^{1-2n}\), so the degree-37 interpolation error is at most

\[
E_{\rm int}=\frac{1400(33/4)^{38}}{38!}<1.8\times10^{-7}.
\]

Suppose each true node value has radius at most \(\epsilon\). Discrete cosine orthogonality gives, for the Chebyshev coefficients \(a_k\) of its error interpolant,

\[
a_0^2+\tfrac12\sum_{k=1}^{37}a_k^2
=\frac1{38}\sum_{j=0}^{37}e_j^2\le\epsilon^2.
\]

Since \(|T_k(x)|\le1\), weighted Cauchy–Schwarz implies
\(\|p_{\rm error}\|_\infty\le\sqrt{75}\epsilon\).
This avoids amplifying uncertain node data by the enormous intermediate monomial coefficients.

A stronger stability factor was subsequently **computed and certified**, rather than assumed. Every cardinal Lagrange polynomial has roots at all the other interpolation nodes. Therefore its sign is known on each of the 39 intervals between adjacent nodes and the domain endpoints. On each interval the Lebesgue function \(\sum_j|\ell_j(\alpha)|\) is one signed-sum degree-37 polynomial. The checker encloses its coefficients in Arb, enlarges the interval to outward dyadic endpoints, and applies the same exact rational Bernstein maximum bound. All 39 intervals were covered, giving

\[
\Lambda_{38}\le 3.278312312437324.
\]

The `--lebesgue` option recomputes and retains this certificate. It reduces the amplification of a common node radius from approximately 8.66 to approximately 3.28. This is useful when the PDE residual radii are around \(1.7\times10^{-5}\), for which the coarser bound would consume the observed gap. The computation does not presume that the maximum occurs at a domain endpoint.

If a node is validated at a dyadic parameter \(\widetilde\alpha_j\), rather than the exact transcendental root, enlarge its radius by
\(23100|\widetilde\alpha_j-\alpha_j|\). The independent checker encloses each root afresh using Arb and recomputes this distance. It also recomputes the DCT interpolation polynomial of the exact stored dyadic node centers. It does not trust an asserted conversion error in an input receipt.

Finally, the exact dyadic midpoint power polynomial is converted to Bernstein form on \([0,1]\) using only rational arithmetic. Its range lies in the convex hull of its Bernstein coefficients. Exact midpoint de Casteljau subdivision supplies a global upper bound, while endpoints and split midpoints supply attained lower bounds. Every split preserves the interval cover. The independently replayed fixed-family bound is

\[
\sup_{\alpha\in[0,1]}S_\alpha(2)
\le P_{\max}^{\rm upper}+E_{\rm conversion}
+\sqrt{75}\max_j(\epsilon_j+23100|\widetilde\alpha_j-\alpha_j|)
+E_{\rm int}.
\]

## What the independent receipt establishes

Run `python scripts/check_validated_continuum.py` for exact-constant and known-polynomial self-checks. `--receipt PATH` additionally replays supplied node-center interpolation and full-interval parameter arithmetic. Its output states `pde_validation_replayed: false`: this checker intentionally has no ability to certify a supplied residual radius merely because it appears in JSON. The separate spatial/time validation receipt must establish that premise.

The polynomial tests include an endpoint maximum, a rational interior maximum, a constant polynomial, and a quartic with oscillating coefficients. Exact binary float input is rejected at the certificate arithmetic boundary. A synthetic 38-node constant-function pipeline was executed to test the independent conversion and coverage path; this is a checker test, not data from the research model.

No claim of a strict continuum survival advantage follows until all 39 PDE observable enclosures, the full fixed-family cover, and the final positive gap have actually been computed with their errors retained. Any receipt reporting only the verified machinery must continue to distinguish that from a completed sign certificate.

## A time-window corollary, conditional on a certified positive gap

Suppose the completed calculation gives \(\delta>0\) with
\(S_\star(2)-\sup_\alpha S_\alpha(2)\ge\delta\).
The only exit to failure is thermal damage from \(h=1\), at rate \(\eta=1\). The forward equation therefore gives, for every policy in the comparison,

\[
S_p'(t)=-\Pr_p\{h_t=1,\,\tau>t\}\in[-1,0].
\]

For every fixed \(\alpha\), the difference \(D_\alpha=S_\star-S_\alpha\) has derivative in \([-1,1]\). Integration, uniformly in \(\alpha\), proves

\[
S_\star(t)-\sup_{\alpha\in[0,1]}S_\alpha(t)
\ge\delta-|t-2|.
\]

Thus the strict advantage persists throughout the explicit interval \(|t-2|<\delta\). This is a conservative certified time interval and requires no new numerical solve.

There is also a qualitative statement about a larger component. Joint continuity in time and intensity, together with compactness of the intensity interval, makes the best-fixed survival continuous. The earlier strict short-time and eventual passive-comparator losses then imply that the connected positive component containing \(t=2\) has a finite left boundary above the short-time loss interval and a finite right boundary below the eventual-loss regime. Its boundary values are zero. This does not prove that there are exactly two crossings, nor that the entire positive set is connected, nor locate the actual transition times.
