# Stochastic calculus analysis: Brownian repair with finite fuel

This is the mathematical analysis behind the implemented spatial repair study. It exposes the stochastic differential equation, compensated reaction clocks, Itô generator, stopping-time martingale, backward survival equation and calibration problem. It uses established stochastic calculus; none of these tools is claimed as a new theorem. [Model, numerical results and verification](SPATIAL-REPAIR.html) · [Lessons and novelty assessment](SPATIAL-LESSONS.html)

## 1. The continuous and discrete state variables

Let \(X_t\in\mathbb T_L^3\) be the relative position of two independently diffusing complexes. Let \(C_t=(h_t,b_t)\), with \(h\in\{0,1,2\}\) healthy complexes and \(b\in\{0,\ldots,N\}\) fuel molecules. Damage count is \(2-h\) and waste count is \(N-b\). The Brownian equation is

\[
dX_t=\sqrt{2D}\,dW_t\pmod L,\qquad D=D_1+D_2.
\]

It is distinct from the older count-only models. Particle position is continuous, while chemistry changes through discrete jumps. Fuel is treated as fast-mixed, rather than individually diffusing.

Let \(g(x)=(1-|x|_{\rm per}^2/a^2)_+^4\). The feasible jumps and their rates are

\[
\begin{array}{lll}
\nu_1=(-1,0),&a_1(x,h,b)=\eta h,&h>0,\\
\nu_2=(1,0),&a_2(x,h,b)=\eta(2-h)/2,&h<2,\\
\nu_3=(1,-1),&a_3(x,1,b)=\kappa b g(x),&b>0,\\
\nu_4=(-1,1),&a_4(x,2,b)=\kappa(N-b)g(x)/2,&b<N.
\end{array}
\]

Rates vanish when the source conditions fail. The factor two in the aggregate reverse repair counts two labelled reverse channels. All rates are bounded because the chemical inventory is finite and \(0\le g\le1\).

## 2. Compensated reaction clocks

Using independent unit-rate Poisson processes \(Y_r\), independent of \(W\), the random time-change representation is

\[
C_t=C_0+\sum_{r=1}^4\nu_r
Y_r\!\left(\int_0^t a_r(X_s,C_{s-})\,ds\right).
\]

Write \(N_r(t)\) for each time-changed count and define

\[
\widetilde N_r(t)=N_r(t)-\int_0^t a_r(X_s,C_{s-})\,ds.
\]

These are compensated counting martingales. Together with Brownian motion they account for two distinct sources of noise: continuously fluctuating separations and discrete chemical events. The bounded-rate construction rules out infinitely many jumps in finite time. [Anderson and Kurtz, chemical-network formulation](https://people.math.wisc.edu/~dfanderson/papers/SURVEY_AndKurtz.pdf)

## 3. Itô's formula and the generator

For a function \(f(x,c)\) that is twice continuously differentiable and periodic in \(x\), let \(\Delta_r f=f(x,c+\nu_r)-f(x,c)\). Itô's formula with jumps gives

\[
df(X_t,C_t)=\mathcal Lf(X_t,C_{t-})\,dt
+\sqrt{2D}\,\nabla_xf(X_t,C_{t-})\cdot dW_t
+\sum_r\Delta_rf(X_t,C_{t-})\,d\widetilde N_r(t),
\]

where

\[
\mathcal Lf=D\Delta_x f+\sum_r a_r\Delta_r f.
\]

The stochastic terms have zero expectation under the stated boundedness and regularity conditions. The predictable quadratic variation of their sum is

\[
\langle M_f\rangle_t=\int_0^t\left[
2D|\nabla_xf|^2+\sum_r a_r(\Delta_rf)^2\right](X_s,C_{s-})\,ds.
\]

This identifies both fluctuation mechanisms explicitly. The generator and martingale also hold for the full uncensored process, where damaged complexes may later recover.

## 4. First collapse and the stopped martingale

Define \(\tau=\inf\{t:h_t=0\}\). For any admissible test function,

\[
M_f(t)=f(X_{t\wedge\tau},C_{t\wedge\tau})-f(X_0,C_0)
-\int_0^{t\wedge\tau}\mathcal Lf(X_s,C_{s-})\,ds
\]

is a martingale. We stop at first collapse, rather than count recoveries after collapse. The killed process has payoff zero at a cemetery state. Its live-state survival function satisfies

\[
\partial_tu=\mathcal L_{\rm kill}u,\qquad
u(0,x,h,b)=1\ (h>0),\qquad u(t,x,0,b)=0.
\]

Periodic spatial boundary conditions apply. For an initial position law \(\mu\),

\[
S_\mu(t)=\int u(t,x,h_0,b_0)\,\mu(dx).
\]

Only damage from \(h=1\) causes immediate collapse. Consequently,

\[
S_\mu'(t)=-\eta\Pr(\tau>t,h_t=1),\qquad -\eta\le S_\mu'(t)\le0.
\]

This identity supplies a rigorous temporal interpolation bound for the numerical survival curves. It does not supply a spatial mesh-error bound.

## 5. Exact short-time consequences

The killed generator applied to the constant live payoff is \(-\eta\) in HD and zero in HH. Applying the generator again gives

\[
\begin{aligned}
S_{\mu,HD,b}'(0)&=-\eta,\\
S_{\mu,HD,b}''(0)&=\eta\left(\frac32\eta+\kappa b\mathbb E_\mu g\right),\\
S_{\mu,HH,b}'(0)&=0,\\
S_{\mu,HH,b}''(0)&=-\eta\left(2\eta+\frac{\kappa(N-b)}2\mathbb E_\mu g\right).
\end{aligned}
\]

The diffusion term vanishes at this second application because the first generator image is position-independent. These are exact short-time formulas for the Brownian model. Position still affects higher-time behavior through the coupled dynamics. Different preparations with different expected contact activity cannot both be reproduced by one preparation-blind count model. Nearness can promote repair or its reverse, depending on the chemical launch.

## 6. Feynman–Kac calibration of a count rate

Freeze chemistry at HD,b and retain only the repair clock. Its survival is

\[
\Pr_\mu(\tau_b>t)=\mathbb E_\mu\exp\left[-\kappa b\int_0^tg(X_s)\,ds\right].
\]

The mean \(m_b(x)=\mathbb E_x\tau_b\) solves

\[
(-D\Delta+\kappa b g)m_b=1
\]

periodically. The implemented count calibration is \(c_b=1/[b\langle m_b\rangle]\), where the average is uniform. Each reverse edge receives the same multiplier as its forward partner, retaining the thermodynamic ratio. This matches one mean per fuel edge; it cannot match every encounter distribution.

For a uniform stationary Brownian launch, Jensen's inequality yields

\[
\Pr(\tau_b>t)\ge e^{-\kappa b\langle g\rangle t},
\qquad c_b\le\kappa\langle g\rangle.
\]

Thus matching the mean can already change the initial effective reaction intensity. This explains why the calibrated count model need not reproduce even the uniform short-time survival curvature.

## 7. What is proved, computed and still open

The stochastic equations, generator identities, event ledgers and short-time consequences are mathematical deductions within the model. The implemented solver discretizes space, exploits exact finite-grid symmetry and propagates a reversible killed generator with a controlled exact-arithmetic polynomial tail. Four spatial refinements provide numerical convergence evidence. Independent Brownian/Poisson thinning avoids a Brownian timestep approximation in this force-free spatial model. It would require modification for a future state-dependent binding force.

![Spatial curves and calibrated count comparison](../results/figures/spatial-repair-survival.png)

There are 51 passing tests, 200,000 formal continuum paths and 200,000 additional fresh-seed paths. Statistical confidence, arithmetic truncation and observed spatial refinement are separate. No rigorous continuum discretization enclosure, molecular calibration, or new general stochastic-calculus theorem is claimed. [Formal verification](../results/spatial-repair-verification.json) · [Independent review paths](../results/spatial-repair-review-check.json) · [Originality audit](../research/spatial-repair-originality-audit.md)
