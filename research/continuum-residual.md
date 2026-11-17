# A weighted dual-residual bound for continuum survival

The earlier lattice optimizer proves a gap near \(1.46\times10^{-4}\) at \(T=2\) for the finite chains. This note gives an exact continuum error identity and an exploratory assessment of the size of its error term. The sampled inner products used in the assessment have no interval enclosure. The strict sign at \(T=2\) remains unresolved by this method.

## Killed operator in the Gibbs space

Let \(i=(h,b,q)\) run over the 20 live chemical states. On the unit torus take \(\rho_i(x)=c_{h,b}e^{-U_q(x)}\), where \(c_{h,b}\) is the full chemical equilibrium mass and \(U_q\) is the published potential. The real inner product is \((u,v)_\rho=\sum_i\int\rho_i u_i v_i\). Use the unit survival payoff \(f_i(x)=1\) and the launch density \(\ell_i(x)=\mathbf1_{i=(1,4,0)}/c_{1,4}\). The desired scalar is

\[
S(T)=(f,e^{-TA}\ell)_\rho,\qquad A=-\mathcal L_{\rm kill}.
\]

The generator of the closed model is reversible. Its killed form consists of the reversible diffusion term \(D\sum_i\int\rho_i|\nabla v_i|^2\), nonnegative jump-energy terms across live edges, and the killing term for exits to failure. The term comprising thermal health jumps and killing has, at each \((x,b,q)\), the reversible two-state matrix \(\eta\left(\begin{smallmatrix}3/2 &-1/2\\-2&2\end{smallmatrix}\right)\). Its smaller eigenvalue is \(\lambda_0=\eta(7-\sqrt{17})/4>0\). All other chemical jump forms are nonnegative. It follows that

\[
(v,Av)_\rho\ge \lambda_0\|v\|_\rho^2+D\sum_i\int\rho_i|\nabla v_i|^2.
\]

The form is closed on the product periodic \(H^1\) space because every \(\rho_i\) is bounded above and below by positive constants and the rates are bounded. It defines a positive self-adjoint operator. The density is continuously differentiable with a Lipschitz gradient, and smooth periodic trial functions lie in the operator domain. This reasoning uses the diffusion form directly, so no derivatives of the Metropolis kink at \(U_1=0\) enter.

## Exact error identity

Choose sufficiently regular curves \(U,V\) with \(U(0)=f\), \(V(0)=\ell\). Let \(r_U=U'+AU\), \(r_V=V'+AV\). The result also applies to continuous time-piecewise polynomial curves in the operator domain. Define

\[
C=(f,V(T))_\rho-\int_0^T(U(T-t),r_V(t))_\rho\,dt.
\]

Then

\[
|S(T)-C|\le
\left(\int_0^T\|A^{-1/2}r_U(t)\|_\rho^2dt\right)^{1/2}
\left(\int_0^T\|A^{-1/2}r_V(t)\|_\rho^2dt\right)^{1/2}.
\]

Indeed, \(e(t)=e^{-tA}\ell-V(t)\) has \(e(0)=0\) and \(e'+Ae=-r_V\). Differentiate \((U(T-t),e(t))_\rho\) and integrate to obtain

\[
S(T)-C=-\int_0^T(r_U(T-t),e(t))_\rho\,dt.
\]

The error equation with the energy test function \(e\) yields \((d/dt)\|e\|_\rho^2+\|A^{1/2}e\|_\rho^2\le\|A^{-1/2}r_V\|_\rho^2\). Integration and Cauchy–Schwarz give the claim. For the surrounding framework of goal-oriented parabolic residual estimation, see [Wu, van der Zee, Simsek and Van Brummelen](https://arxiv.org/abs/1706.04281).

For an exact spatial Galerkin solution with the same spaces for both curves, the residual is orthogonal to the trial space and the correction integral vanishes. To enclose a floating-point or quadrature-based approximate solution, use the correction integral and residuals evaluated for that actual approximation. This avoids an assumption of exact Galerkin orthogonality for the reported computed coefficients.

## A computable dual norm by weighted Fourier antidivergence

For any residual \(r\), write \(R_i=\rho_i r_i\), \(\bar R_i=\int R_i\), and \(Z_i=\nabla(-\Delta)^{-1}(R_i-\bar R_i)\) with mean-zero inverse Laplacian. Then \(R_i=\bar R_i-\nabla\cdot Z_i\). Integration by parts and Cauchy–Schwarz against the thermal and diffusion terms give

\[
\|A^{-1/2}r\|_\rho^2\le
\sum_i\left[\frac1D\int\frac{|Z_i|^2}{\rho_i}
+\frac{|\bar R_i|^2}{\lambda_0}\int\frac1{\rho_i}\right].
\]

The proof applies Cauchy–Schwarz once to the direct sum of the gradient and scalar terms. The source in exact spatial Galerkin has \(\bar R_i=0\) because the constant basis functions are included in every live state.

Finite Fourier validation can avoid placing an interval around every quadrature-grid point. Compute enclosures for \(\widehat R_i(k)\) through \(\|k\|_\infty\le M\), the corresponding weighted squared norm of \(Z_i^{(M)}\), and \(\int R_i^2\). Then

\[
\left\|\frac{Z_i-Z_i^{(M)}}{\sqrt{\rho_i}}\right\|_2^2
\le\frac{\int R_i^2-\sum_{\|k\|_\infty\le M}|\widehat R_i(k)|^2}
{\inf\rho_i\,[2\pi(M+1)]^2}.
\]

Use the upper end of the total norm and the lower ends of the resolved squared coefficients in this subtraction. Minkowski bounds the full weighted antidivergence norm by the sum of the resolved and unresolved norms. Apply the same inequality to the space-time norm if integrating both squared terms over time. The step supplies an explicit unresolved-frequency error, whereas merely measuring more Fourier modes does not bound the omitted ones.

All quantities in this proposal can be reduced to finite operations. Archive continuous, time-piecewise polynomial trial curves with exact dyadic coefficients, enforcing exact initial values and common values at piecewise junctions. For each time segment form the residual from the continuum coefficients and the exact derivative of the stored polynomial. Radial spatial coefficients and their products reduce to integrals of the form

\[
\widehat F(k)=4\pi\int_0^a r^2F(r)\operatorname{sinc}(2\pi |k|r)\,dr
\]

after subtracting the constant value outside contact. Split the radial integrals at \(a\) and at \(a\sqrt{1-\sqrt{\delta/E}}\), the switch in the acceptance formula. Polynomial factors in Cartesian position caused by \(\nabla U\) can be handled by differentiating the compact-support Fourier transform with respect to the frequency coordinates. For \(R_i^2\) these have degree at most two. Inside each radial branch the factors are polynomial and exponential functions and entire sinc functions. Interval Taylor integration on subdivided panels can bound each integral and its derivatives, treating the zero frequency with the sinc series. Integrate the time polynomials and their products exactly or outwardly, including the shifted correction integrand. Quantize the approximate trial curves before validation. Then their ODE propagation errors need no independent estimate because their exact residuals include those errors.

## Checking the whole fixed comparator interval

For the fixed family, \(A_\alpha=A_0+\alpha B\), with \(B\) nonnegative and bounded in the common Gibbs space. At fixed position each chemical paid edge gives an independent reversible two-state block. Its nonzero eigenvalue is the sum of the forward and reverse rates. Thus \(\|B\|_\rho\le k_1(N+1/8)=16.5\) for the stated rate ceilings. The same inequality holds for the complexification.

Writing \(F(x)=S((x+1)/2)\), complexify \(x\) to the Bernstein ellipse \(x=(z+z^{-1})/2,\ |z|=3/2\). The resulting \(\operatorname{Re}\alpha\ge-1/24\). The thermal bound yields

\[
\|e^{-TA_\alpha}\|_\rho\le \exp[-T(\lambda_0-16.5/24)]\le 1.
\]

This follows from the real-part energy estimate for the complex bounded perturbation. Analyticity in \(\alpha\) follows from the bounded-perturbation series on any compact complex set. Consequently \(|F(x)|\le\|f\|_\rho\|\ell\|_\rho<700\) on the ellipse. For this model use \(\|f\|_\rho^2\le(5/9)(1+e^4)\) and \(\|\ell\|_\rho^2=1/[(4/9)(1/9)^4]\).

For degree \(m\) Chebyshev–Lobatto interpolation, the exact interpolant error is bounded by \(4\cdot700(3/2)^{-m}/((3/2)-1)\) by the Chebyshev coefficient bound and aliasing of each omitted cosine into a bounded lower cosine. At \(m=64\) this is below \(3.1\times10^{-8}\). Enclose each of the 65 true node values using the residual method, convert the resulting interval interpolant to a polynomial with verified coefficient errors or use a validated Lebesgue factor, and bound its maximum on the whole interval by interval subdivision. A point comparison with \(\alpha=1\) alone cannot certify the supremum. This procedure would cover the full fixed family without a dense mesh based on a first-derivative bound. It has not been executed with validated continuum node enclosures.

## Exploratory calculation and what it identifies

The [replayable program](../scripts/explore_weighted_continuum_galerkin.py) builds a cosine Galerkin approximation on symmetric modes with coordinate frequencies up to \(K\), and uses midpoint spatial quadrature and Gauss time quadrature. Reflection and permutation invariance reduce the basis size. The program evaluates the dual norm using the sampled antidivergence. It does not bound the mid-point integration error at the contact and acceptance interfaces, the time quadrature error, or the Fourier frequencies absent from that sample. The sampled product is an indicator of the budget required by the theorem, not a validated error radius.

At \(K=10\) with 24 midpoint nodes per folded coordinate, estimated survival is \(0.2337171840\) for the spectator and \(0.2335705577\) for the fixed \(\alpha=1\) sample. The sampled dual-residual products are \(8.49\times10^{-6}\) and \(8.51\times10^{-6}\). With \(K=12\) and 32 midpoint nodes the two estimates are \(0.2337156043\) and \(0.2335688374\), a gap estimate of \(1.467668648\times10^{-4}\). The corresponding sampled products are \(2.4922\times10^{-6}\) and \(2.4934\times10^{-6}\). These magnitudes show the validation budget for this route. A rigorous continuum result requires the integral enclosures and fixed-family optimization described above. [Execution receipt](../results/continuum-residual-screen.json)
