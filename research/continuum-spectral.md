# Fourier checks and an error inequality for the continuum survival comparison

The independent Fourier calculations support a positive gap near 0.000147 in absolute survival probability at time two against the constant intensity \(\alpha=1\). The saved calculations use numerical quadrature and floating arithmetic. They do not give an interval proof of continuum superiority at time two or an optimizer certificate for all intensities.

## Numerical comparisons

| Discretization | Cutoff per spatial axis | Spectator survival | Constant intensity one | Difference |
|---|---:|---:|---:|---:|
| Gibbs conjugation, then project | 4 | 0.23296641 | 0.23274586 | 0.00022056 |
| Gibbs conjugation, then project | 6 | 0.23371085 | 0.23356318 | 0.00014767 |
| Gibbs conjugation, then project | 8 | 0.23370989 | 0.23356236 | 0.00014753 |
| Gibbs conjugation, then project | 10 | 0.23371131 | 0.23356400 | 0.00014731 |
| Gibbs conjugation, then project | 12 | 0.23371311 | 0.23356603 | 0.00014708 |
| Project in the original Gibbs space | 8 | 0.23371625 | 0.23356949 | 0.00014676 |
| Project in the original Gibbs space | 12 | 0.23371589 | 0.23356913 | 0.00014676 |

Commands from the repository root are `python -B scripts/run_continuum_fourier.py` and `python -B scripts/run_continuum_weighted.py`. Replay data include source hashes and quadrature tolerances or error estimates. The estimates from `quad_vec` are heuristic. [Conjugated receipt](../results/continuum-fourier.json) · [Weighted receipt](../results/continuum-weighted.json)

The flat operator comes from the similarity transformation by the square root of the Gibbs density. Its diffusion part is \(D\Delta+W_q\), where

\[
W_q=D\left(\frac{\Delta U_q}{2}-\frac{|\nabla U_q|^2}{4}\right).
\]

The transform of the chemical jump generator is symmetric pointwise. In particular the transformed off-diagonal coefficient for a switch with directional bare rates \(f,r\) is \(\sqrt{fr}\,g e^{-|U_1|/2}\). The diagonal loss terms retain their original rates. The distributional Laplacian has a bounded jump at the contact radius. There is no surface measure because \(U_1\) and its gradient are continuous. A radial multiplier compactly supported inside radius \(a\) has Fourier coefficient

\[
\widehat v(k)=4\pi\int_0^a r^2 v(r)\operatorname{sinc}(2\pi |k|r)\,dr,
\qquad \operatorname{sinc}(y)=\frac{\sin y}{y}.
\]

The continuum launch is uniform in one open chemical state, so its flat transformed vector contains only the zero mode. Retaining frequency components to \(K\) needs multiplier coefficients only through \(2K\). Multiplication of a \(2K\) field by a \(K\) field can reach \(3K\), so period \(4K+2\) of the padded transform cannot alias any resulting frequency to \([-K,K]\). This is a finite projection identity in exact arithmetic, independent of the quality of the input coefficient quadrature.

In the weighted projection, the diffusion stiffness matrix uses the Fourier coefficients of \(e^{-U_q}\) directly. A pair of complex basis functions of wavevectors \(i,j\) gives diffusion stiffness \(4\pi^2D(i\cdot j)\widehat{e^{-U_q}}(i-j)\). The mode-space mass matrix uses the same coefficients without \(4\pi^2D(i\cdot j)\). In the signed-permutation invariant sector write \(O_i\) for the orbit of the nonnegative sorted vector \(i\). A radial multiplication matrix element in normalized orbit coordinates is

\[
\sqrt{|O_i|/|O_j|}\sum_{k\in O_j}\widehat v(i-k).
\]

Cholesky conversion of each conformation mass block gives a symmetric matrix for the semigroup. Both the weighted launch density and the terminal payoff are spatially constant within each chemical state, so both belong exactly to the weighted trial space for all \(K\).

## A rigorous inequality for the weighted projection error

On the live states use the original Gibbs Hilbert space, or multiply each component by the square root of its chemical Gibbs weight and use spatial weights \(w_q=e^{-U_q}\). Let \(L\) be the killed generator and \(A=-L\). Reversible diffusion, chemical jumps, and killing make \(A\) positive. In fact the two-health-state damage block alone in each pair \((b,q)\) gives

\[
A\succeq \frac{7-\sqrt{17}}{4}\eta I.
\]

Indeed this block has the positive symmetric matrix \(\eta\begin{pmatrix}3/2&-1\\-1&2\end{pmatrix}\). Other terms are nonnegative quadratic forms. The closed form of \(A\) dominates \(D\sum_j\int w_{q_j}|\nabla v_j|^2\), where \(j\) runs through live chemical states. The compact spatial domain, bounded Lipschitz drift, bounded reaction rates, and positive finite weights allow the usual form realization and energy identity. [Background on duality-based estimates for parabolic problems](https://arxiv.org/abs/1706.04281)

Let \(P\) be the exact weighted orthogonal projection onto all cube Fourier modes \(|k_i|\le K\) in each component. One may equivalently restrict this to the signed-permutation invariant sector for the present data and operator. Let \(f\) be the initial density in the chosen Hilbert coordinates and \(g_0\) the terminal payoff. Both are in the trial space. Define \(a(t)=e^{tPLP}f\), \(b(t)=e^{tPLP}g_0\) in the trial space, and their spatial residuals \(r_a=La-a'\), \(r_b=Lb-b'\). These are orthogonal to the trial space at each time.

The error \(e(t)=e^{tL}f-a(t)\) has zero initial condition and satisfies \(e'+Ae=r_a\). The energy identity and Cauchy–Schwarz give

\[
\int_0^T\|A^{1/2}e\|^2\,dt\le\int_0^T\|A^{-1/2}r_a\|^2\,dt .
\]

To see the factor, integrate \(\tfrac12(\|e\|^2)' + \|A^{1/2}e\|^2=\langle A^{-1/2}r_a,A^{1/2}e\rangle\) and use Cauchy–Schwarz before dividing by the nonnegative energy norm. Differentiate \(\langle b(T-t),e(t)\rangle\). Self-adjointness cancels the generator terms, and \(\langle b(T-t),r_a(t)\rangle=0\). Thus the goal difference is \(\int_0^T\langle r_b(T-t),e(t)\rangle dt\). This proves

\[
\big|\langle g_0,e^{TL}f\rangle-\langle g_0,a(T)\rangle\big|
\le\|A^{-1/2} r_a\|_{L^2(0,T)}\;\|A^{-1/2} r_b\|_{L^2(0,T)}.
\]

For either residual put \(\rho_j=w_{q_j}r_j\). Galerkin orthogonality implies zero mean and vanishing retained Fourier coefficients. Define \(G_j=\nabla\Delta^{-1}\rho_j\), which satisfies \(\nabla\cdot G_j=\rho_j\). Integration by parts and the diffusion contribution to the form imply

\[
\|A^{-1/2}r\|^2\le \frac1D\sum_j\int_{\mathbb T^3}\frac{|G_j(x)|^2}{w_{q_j}(x)}\,dx .
\]

This bound uses the actual smoothness of the residual rather than imposing high-order differentiability on the Metropolis rate at its kink. The weighted prototype computes numerical diagnostics of the right-hand side at five sampled times. At \(K=12\) the largest sampled square roots are about \(0.0007322\) for the launch and \(0.003825\) for the payoff. These numbers are consistent with an eventual two-residual estimate below \(10^{-5}\) per curve. They are not bounds between the observation times, and the truncation of the spatial flux introduces additional uncertainty.

One can enclose the omitted flux explicitly after validating coefficient integrals. If \(P_M\) is the cube projection at cutoff \(M\), then

\[
\left\|\frac{(I-P_M)G_j}{\sqrt{w_{q_j}}}\right\|\le
\frac{1}{2\pi(M+1)\sqrt{\inf_x w_{q_j}(x)}}\;\|(I-P_M)\rho_j\| .
\]

Parseval bounds the last norm by the difference between an enclosed total square norm and lower bounds on resolved Fourier coefficient square sums. The radial coefficients of the required integrals split at the contact radius and the sphere where \(U_1=0\). Cartesian monomials in terms such as \(\nabla U_1\cdot\nabla a\) can be accounted for by differentiating the radial transform with respect to its frequency vector. The remaining numerical tasks are outward bounds for these integrals, reconstruction of the entire time paths with a validated matrix exponential or polynomial residual, and verification of roundoff in the projection and norm calculations. A direct norm bound for an inexact Galerkin path also needs to include the retained residual and its pairing with the approximate dual, which the exact Galerkin proof eliminates.

## Covering all constant intensities

For the fixed family \(L(\alpha)=L_0+\alpha B\), the full-space paid generator \(B\) consists of disjoint pairs within each \(x,h\). In the Gibbs norm its nonzero pointwise eigenvalue magnitudes are paid forward plus reverse rates, at most \(4(4+1/8)=16.5\) for the published parameters. Since all real \(\alpha\in[0,1]\) generate contractions, the iterated Duhamel formula bounds the derivatives of the survival objective:

\[
|S^{(n)}(\alpha)|\le\|f\|\,\|g_0\|(T\|B\|)^n .
\]

The launch can be normalized to \(\|f\|=1\) by choosing the chemical weight scale at its initial state. The chemical sum of live weights divided by the launch weight is \(8201.25\). Multiplying by \(1+e^4\) bounds the spatial-weight factor, giving \(\|f\|\|g_0\|<700\).

Interpolate at the \(n\) roots of the degree-\(n\) Chebyshev polynomial mapped to \([0,1]\). The monic nodal product is bounded by \(2^{1-2n}\), so at \(T=2\), the interpolation remainder is at most \(1400(8.25)^n/n!\). At \(n=38\) this is below \(1.8\times10^{-7}\). If every nodal enclosure has radius at most \(\epsilon\), the interpolation uncertainty is at most \(\sqrt{2n-1}\epsilon\). Indeed the discrete orthonormal cosine expansion at the roots gives coefficient squared sum at most \(\epsilon^2\), and the sum of squared basis values anywhere is at most \(2n-1\). Maximizing the interpolating polynomial with a verified interval subdivision then supplies an upper bound for the entire fixed family. This reduces the all-intensity comparison to finitely many validated evaluations with a quantitative accuracy requirement. The two numerical receipts above do not perform that enclosure.
