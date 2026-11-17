# Path-law error for the continuum simulation

The three Euler time steps in the [existing particle receipt](../results/organization-continuum-check.json) do not come with a bound on spatial discretization bias. The bound below applies to their exact-chemistry Euler interpolation and to every constant paid intensity in the published family. The bound is too wide to settle the sign at time two. It does give a continuum error bound whose inputs and error budget can be checked directly.

## A bound for a bounded step and every permitted intensity

Use \(D=1/40\), \(E=6\), \(a=8/25\), \(N=4\), and \(T=2\). The displaced position on the unit torus has drift zero while the interaction mode is open. In the closed mode its drift is \(F(x)=-D\nabla U_1(x)\). On the contact ball, write \(z=|x|^2/a^2\). Then

\[
F(x)=-\frac{4DE}{a^2}(1-z)x,\qquad
B:=\sup |F|=\frac{8DE}{3\sqrt3a},\qquad
L:=\operatorname{Lip}(F)\le\frac{8DE}{a^2}.
\]

The tangential and radial derivatives of the drift have eigenvalues \(-4DE(1-z)/a^2\) and \(-4DE(1-3z)/a^2\), respectively. The drift is zero outside the contact ball and continuous across its boundary. Integrating its derivative along a lifted periodic line segment gives the global Lipschitz bound. The maximum of \(r(1-r^2/a^2)\) occurs at \(r=a/\sqrt3\), giving \(B\).

Let \(P_h\) and \(P\) denote the Euler and exact laws on the tagged proposal paths. Each proposal waiting time is sampled from a rate bound determined by the current chemical state. Conditional on the state and the waiting time, that interval's position process has fixed interaction mode. It is propagated up to the proposal or observation boundary. The chemical mark then uses the same kernel under both laws, evaluated at the actual proposed position. This representation applies to the implementation's subdivisions of candidate intervals into steps of length at most \(h\). The waiting-time choice is independent of future Brownian increments conditional on the starting state.

On a closed interval with conditioned Euler step size \(\delta\), the interpolation has constant drift \(F(X_0)\). Under that law the conditional mean squared lifted displacement at time \(u\) is \(u^2|F(X_0)|^2+6Du\). The path relative entropy of this interpolation with respect to the true diffusion follows from Girsanov's formula. It is bounded by

\[
\frac{1}{4D}\int_0^\delta \mathbb E_h|F(X_0)-F(X_u)|^2\,du
\le L^2\left(\frac{B^2\delta^3}{12D}+\frac{3\delta^2}{4}\right).
\]

The same bound holds for a wrapped path by data processing. Both drifts are bounded, so the exponential martingales in this comparison are integrable. Conditional relative entropies add by the chain rule over proposal intervals and their Euler steps. Waiting times and chemical marks contribute zero entropy when conditioned on the same history. The result for the stopped path follows either by applying the chain rule to finite prefixes and passing to the limit or by forgetting the unused increments. Open-mode position updates are exact. With \(C_h\) the time spent alive in the closed mode under the Euler law, use \(\sum \delta^2\le h C_h\) and \(\sum\delta^3\le h^2 C_h\).

Pinsker's inequality bounds the difference in the probability of any event on these paths, including survival to the first-failure horizon. Hence for each policy \(\pi\) under consideration,

\[
|S_{h,\pi}(T)-S_\pi(T)|\le
L\sqrt{\mathbb E_h[C_h]\left(\frac{3h}{8}+\frac{B^2h^2}{24D}\right)}.
\]

This path-entropy method follows the same change-of-measure principle used for Euler approximations of Langevin diffusions. Here the chain rule accommodates independent proposals on intervals of random length and exact conditional chemical marks. [Dalalyan, 2017](https://arxiv.org/abs/1412.7392)

## Evaluating the occupation and error budgets

The initial open position is uniform. Before its first closure, the position is Brownian with zero drift. The total closure intensity is bounded by \(17g(x)\) for the spectator and all constants because \(k_0=1\), \(k_1N=16\), \(0\le m\le1\), and every paid factor is at most one. Use independent thinning on an un-killed free Brownian path to dominate the first closure by a Cox time with hazard \(17g\). The Euler scheme is exact up to this first closure. The uniform spatial average is

\[
\overline g=4\pi a^3\int_0^1 y^2(1-y^2)^4,dy
=\frac{512\pi a^3}{3465}.
\]

Let \(\gamma=17\overline g\). Jensen's inequality shows the dominating Cox time survives to \(t\) with probability at least \(e^{-\gamma t}\). The actual amount of closed, alive time is at most the time remaining after first closure. Therefore

\[
\mathbb E_h C_h\le c_T:=T-\frac{1-e^{-\gamma T}}{\gamma}.
\]

The bound is \(c_2<0.439\). For example, \(\gamma<0.2587\) follows from \(\pi<22/7\), and integrating the degree-five Taylor upper bound on \(1-e^{-\gamma t}\) at \(0.2587\) gives a rational value less than \(0.439\). Also \(B^2=25/48\) and \(L=375/32\). Substituting these rational upper bounds gives the following per-policy errors, rounded upward.

| Maximum step \(h\) | Upper bound on absolute continuum bias in survival |
|---:|---:|
| 0.01 | 0.481 |
| 0.005 | 0.339 |
| 0.0025 | 0.239 |

For any uniform per-policy bias bound \(e\), the absolute error in the spectator-minus-optimum gap is bounded by \(2e\). The inequality \(|\sup_\alpha f(\alpha)-\sup_\alpha g(\alpha)|\le\sup_\alpha|f(\alpha)-g(\alpha)|\) controls optimization over all constant intensities without a finite set of sampled intensities. A sampling error or numerical propagation error in the Euler objective would still require a separate budget.

For illustration, a per-policy bias budget of half the \(48^3\) certified lower gap, about \(7.30\times10^{-5}\), reaches the scale \(h\approx2.36\times10^{-10}\) in this sufficient bound. That spends the entire gap budget on bias and leaves none for sampling. Even without bias, a direct DKW band of that half-width at 99% confidence jointly for three independently simulated survival curves requires about \(6.0\times10^8\) paths per curve. This covers only those three curves. Verification of the optimizer over the full constant interval still requires further work. These figures quantify why the present path simulation does not settle the sign. They are not lower bounds on the complexity of every simulation or proof.

## An exact transition kernel for a different sampling route

The drift is a gradient and its potential is periodic with globally Lipschitz gradient. Let \(Y\) be Brownian motion with diffusion coefficient \(D\) on the torus. For a prescribed interval of closed conformation, Girsanov and Itô give the path likelihood relative to \(Y\), starting from the same \(x\), of

\[
\exp\left[\frac{U_1(x)-U_1(Y_t)}{2}+\int_0^t W(Y_s)\,ds\right],\qquad
W=\frac D2\Delta U_1-\frac D4|\nabla U_1|^2.
\]

Itô applies by smooth approximation. The bounded weak second derivatives agree with the classical derivatives away from the contact sphere. A Brownian path spends zero Lebesgue time on that sphere, so there is no ambiguity in the integral. On the contact ball,

\[
W=\frac{2DE}{a^2}(3-17z+24z^2-12z^3).
\]

The factor is strictly decreasing for \(0\le z\le1\), since its derivative is \(-17+48z-36z^2<0\). Off the ball \(W=0\), giving bounds \(W\in[-5.859375,8.7890625]\). The potential ranges from \(-4\) to \(2\).

This supplies a constructive exact sampler in the ideal real-arithmetic model. Propose a Brownian endpoint and skeleton at independently sampled Poisson times of rate \(W_{\max}-W_{\min}\). Accept the endpoint with factor \(\exp[(\min U_1-U_1(Y_t))/2]\). At every Poisson time \(s\), independently reject with probability \((W_{\max}-W(Y_s))/(W_{\max}-W_{\min})\). The probability of passing all Poisson tests conditional on the whole path is \(\exp[-W_{\max}t+\int_0^t W(Y_s)\,ds]\). Thus the net acceptance probability is the likelihood above multiplied by a path-independent constant. The Brownian skeleton is obtainable with successive independent Gaussian increments in chronological order, wrapped modulo one. For fixed durations, retrying has finite expectation. Long durations can be subdivided using the exact semigroup property. Swap this conditional kernel for Euler between chemical proposal times to remove the Euler bias. This construction adapts the rejection principle for exact gradient diffusions. It has not been run as part of this receipt and retains sampling uncertainty if used for the gap. [Beskos and Roberts, 2005](https://discovery.ucl.ac.uk/id/eprint/146969/1/Beskos_Beskos_Roberts_05.pdf)

## Independent short-horizon sign check

There is an exact continuum failure of superiority for short missions. This argument uses separate simple bounds and does not depend on any computed lattice result. Write \(L_\pi\) for the killed continuum generator acting on live-state functions and \(f=1\) on live states. The first generator action is \(-\eta\mathbf1_{h=1}\), and the second is \(\eta(3\eta/2+\kappa b g)\) for \(h=1\) and \(-\eta(2\eta+\kappa w g/2)\) for \(h=2\). Here \(g\) has at least two continuous derivatives across the contact boundary. Therefore \(f\in\operatorname{Dom}(L_\pi^3)\). A paid closure lowers \(b\) by one, immediately changing the second action in the launch state by \(-\eta\kappa g\). All other generator contributions are common at launch. The spectator has full paid intensity at \(b=N\). With \(\mu\) the uniform launch,

\[
S_{\rm spec}'''(0)-S_0'''(0)=-256\langle g^2m(U_1)\rangle_\mu<-13/10.
\]

For the last inequality, integrate only \(r\le a/2\), where \(U_1<0\) and \(m=1\). A rational lower bound is

\[
256\cdot12a^3\sum_{j=0}^{8}(-1)^j\binom 8j\frac{(1/2)^{3+2j}}{3+2j}>13/10,
\]

using \(\pi>3\) in the radial volume element. In the uniform open launch the diffusion contribution to \(\mu L_\pi\) integrates to zero on smooth periodic functions. The remaining signed measure has full total variation bounded by twice the mean outgoing rate, less than \(2(1.5+81\overline g)<6\). Dynkin's formula on the smooth generator core and closure imply \(\mu P_{\pi,t}-\mu=\int_0^t(\mu L_\pi)P_{\pi,s}\,ds\). This need not differentiate the acceptance factor \(m\) across its kink.

For \(L_\pi^3 f\), the diffusion term is at most 375 in magnitude using \(\|\Delta g\|_\infty\le24/a^2\) and the coefficient 64. The drift term is less than 281 using \(B<0.73\) and \(\|\nabla g\|_\infty<6\), whose square maximum occurs at \(z=1/7\). The jump term is at most \(82.5\cdot99.5\), using the outgoing rate cap and the range \([-34,65.5]\) of the second action, including the cemetery value. Thus \(\|L_\pi^3 f\|_\infty<9000\). Contraction gives \(|S_\pi'''(s)-S_\pi'''(0)|\le54000s\). With the first two derivatives common, integration gives

\[
S_{\rm spec}(t)-\sup_{\alpha\in[0,1]} S_\alpha(t)
\le S_{\rm spec}(t)-S_0(t)
\le -\frac{11}{300}t^3,\qquad 0<t\le \frac1{25000}.
\]

This sign follows in the Brownian continuum itself. It concerns an explicit short horizon. It leaves the published comparison at \(T=2\) to a separate bound.

## A strict long-time loss in the Brownian continuum

There is also an exact asymptotic obstruction for the same model, initial state and fixed comparator. In the passive policy \(\alpha=0\), repair preserves \(h+b\), and a passive tether switch preserves both \(h\) and \(b\). On live states \(h=1,2\), the killed thermal generator is

\[
K=\begin{pmatrix}-3/2&1/2\\2&-2\end{pmatrix}.
\]

Its smallest decay rate is \(\lambda_0=(7-\sqrt{17})/4\). Set \(r=3-2\lambda_0=(\sqrt{17}-1)/2>1\) and \(\phi(x,h,b,q)=r^{h+b}\). Then \((r,r^2)^\top\) is a positive eigenvector of \(K\) at \(-\lambda_0\). The passive generator \(L_0\) satisfies \(L_0\phi=-\lambda_0\phi\) because repair, tethering and spatial movement leave this quantity unchanged.

Use the common live restriction of the published Gibbs law \(\pi\) as reference measure. It has a positive bounded density in each state on the compact torus. Write \(\mathcal E_0(f)=\langle f,-L_0 f\rangle_\pi\). Integration by parts in each state and the pointwise paired detailed balance yield the ground-state identity

\[
\mathcal E_0(\phi u)-\lambda_0\|\phi u\|_\pi^2
=D\sum_i\int\pi_i\phi_i^2|\nabla u_i|^2\,dx
+\frac12\sum_{i,j}\int\pi_i a^0_{ij}\phi_i\phi_j|u_j-u_i|^2\,dx.
\]

The sum runs over live edges. The killing term cancels using \(-L_0\phi=\lambda_0\phi\). The identity can first be verified for smooth \(u\) by expansion of the Dirichlet form, and extends to its Sobolev form domain. Vanishing energy forces \(u_i\) constant in space in every state. Thermal recovery and damage connect the live values of \(h\) for each \((b,q)\). Passive tether edges connect the modes, and the repair edges link neighboring fuel counts. Each such connection has positive rate on a positive-volume region. Consequently equality requires \(u\) constant across all states, and the passive minimizers are exactly the multiples of \(\phi\).

For a paid policy \(\rho\), the form for \(-L_\rho\) is \(\mathcal E_0+\mathcal E_{\rm paid,\rho}\). The added form is a sum of nonnegative squared differences by the shared detailed balance. There is an active paid edge at full fuel both for the spectator and for any fixed \(\alpha>0\). Such an edge changes \(h+b\) by one. It follows that \(\mathcal E_{\rm paid,\rho}(\phi)>0\) in either case. The Rayleigh minimum \(\lambda_\rho\) is attained because the form controls the gradient norm and the finite vector Sobolev space embeds compactly into weighted \(L^2\). It must satisfy \(\lambda_\rho>\lambda_0\): equality would require simultaneously a passive minimizer proportional to \(\phi\) and zero added form.

The uniform launch has a bounded density \(\ell=d\mu/d\pi\), and both \(\ell\) and the unit live payoff belong to weighted \(L^2\). Self-adjoint spectral contraction gives \(S_\rho(t)\le \|\ell\|_\pi\|1\|_\pi e^{-\lambda_\rho t}\). Positivity and \(1\ge \phi/\|\phi\|_\infty\) give

\[
S_0(t)\ge \frac{\mu\phi}{\|\phi\|_\infty}e^{-\lambda_0 t}>0.
\]

Hence \(S_{\rm spec}(t)/S_0(t)\to0\) exponentially. Since the optimal constant arm can always choose \(\alpha=0\), the spectator loses to the full fixed family at all sufficiently long horizons.

### A conservative explicit bound on the exponent difference

The same argument gives an explicit long-time witness. Use the unnormalized density \(\pi_i=\binom2h2^{2-h}\binom4b8^{-b}e^{-U_q}\) for each live chemical state, with \(\chi_i=\pi_i\phi_i^2\). From \(12/5<r^2<5/2\) and the energy range, \(1/200<\chi_i<1100\). To verify the lower bound, use \(\binom2h2^{2-h}r^{2h}>144/25\), \(\binom4b(r^2/8)^b>(3/10)^4\), and \(e^{-U_q}>1/8\). To verify the upper bound, use the respective factors \(10\), \(5/4\) by inspection of the five fuel counts, and \(81\).

Let \(K\) be the ball \(|x|\le a/2\), of volume \(\nu>1/64\), and \(m_i\) the average of \(u_i\) on \(K\). Denote the right side of the ground-state identity by \(H\), and write \(P=\mathcal E_{\rm paid,spec}(\phi u)\). Then \(H\ge(1/8000)\sum_i\int |\nabla u_i|^2\). On \(K\), the undirected conductance in \(H\) is at least \(1/800\) for every thermal live edge, passive closure edge, and repair edge. Indeed use the direction \(h=1\to2\), \(q=0\to1\) or forward repair, respectively. The factor \(\phi_j/\phi_i\) is at least one in each direction, the rate is at least \((3/4)^4\), and the transformed starting weight is at least \(1/200\). Therefore, for any length-\(k\) simple path using those edges, Jensen and Cauchy–Schwarz give

\[
|m_i-m_j|^2\le k\frac{800}{\nu}H\le 51200kH.
\]

All 20 live states are connected, so \(k\le19\) suffices for any pair. In particular the full-fuel paid edge from \(i=(1,4,0)\) to \(j=(1,3,1)\) has a passive length-three path via \((2,3,0)\) and \((1,3,0)\). The paid conductance on \(K\) is at least \(1/256\), using the direct forward rate \(16g\). Thus \(|\phi_i m_i-\phi_jm_j|^2\le16384P\). Since \(\phi_i-\phi_j=r^4(r-1)>2\) and \(\phi_j^2=r^8<40\),

\[
m_i^2\le8192P+3{,}072{,}000H,\qquad
m_k^2\le16384P+8{,}089{,}600H\quad\text{for each state }k.
\]

For the torus mean \(\overline u_k\), Jensen gives \(|m_k-\overline u_k|^2\le \nu^{-1}\int|u_k-\overline u_k|^2\). Orthogonality to constants and the torus Poincare inequality give

\[
\sum_k\int |u_k-m_k|^2\le \frac{65\cdot8000}{36}H.
\]

The Poincare factor \(1/(4\pi^2)\le1/36\) follows from the first nonzero Fourier frequency on the unit torus. Decomposing \(u_k=(u_k-m_k)+m_k\) and using the upper weight bound gives

\[
\|\phi u\|_\pi^2\le1100\left[2\frac{65\cdot8000}{36}H+40(16384P+8{,}089{,}600H)\right]
\le 4\cdot10^{11}(H+P).
\]

The Rayleigh principle yields \(\lambda_{\rm spec}-\lambda_0\ge1/(4\cdot10^{11})\). The launch density in this unnormalized reference has \(\|\ell\|_\pi=32\), since the launch-state density of \(\pi\) is \(1/1024\). Also \(\|1\|_\pi^2<20\cdot324<81^2\), and \(\|\phi\|_\infty/(\mu\phi)=r<2\). The survival ratio therefore satisfies

\[
\frac{S_{\rm spec}(t)}{S_0(t)}\le 5184\exp[-t/(4\cdot10^{11})].
\]

In particular \(S_{\rm spec}(t)<S_0(t)\) for every \(t\ge4\cdot10^{12}\). That very large bound is a sufficient witness time, not a numerical crossover prediction. Neither spectral argument determines the sign at time two.
