# Model and derivations

All logarithms are natural. Energy is measured in units of k_B T_bath, and all rates share a time unit. T denotes task duration, not temperature. These are statements within a specified model family, not universal laws of agency.

## 1. Physical system and boundary

An externally imposed binary signal X flips symmetrically at rate r>0. A thermal two-state memory M has energy E(X,M)=A 1[X != M], where A>=0. For fixed X, memory transitions correct a mismatch at rate a and introduce a mismatch at rate b. Parameterize

\[
a=\frac{K}{1+e^{-A}},\quad b=\frac{K}{1+e^{A}},\quad a+b=K>0.
\]

Thus log(a/b)=A, satisfying local detailed balance for the energy difference and a single thermal bath. K is a kinetic capacity (the sum of directional rates), not the stationary number of transitions per unit time. Raising K can require hardware resources; their construction/maintenance costs are not included. A separately specified K ceiling is essential to avoid treating arbitrary speed as free.

Changing X changes the memory's energy landscape. The external work source pays +A for a match->mismatch signal transition and receives A back for mismatch->match. The imposed signal's own generation cost is outside the boundary. The model assumes ideal recovery of negative work; a nonregenerative driver costs more. The four-state chain is a driven thermal sensor, not a passive equilibrium joint X,M system and not an autonomous organism.

Take the effective action to be M. Its failure consequences are supplied by an explicit reliability layer below. There is no separate physical actuator, energy-depletion mechanism, or thermodynamic repair mechanism in this release. Average external work is NOT a finite-battery survival guarantee. These exclusions prevent attributing physical self-maintenance to a model that has only an effective failure criterion.

## 2. Stationary sensing and energetic balance

Let Z=1[X!=M]. This aggregation is lumpable because rates do not depend on which of the two matching (or mismatching) states is occupied. Define

\[
\begin{aligned}u&=r+b,\quad v=r+a,\quad s=u+v=2r+K,\\ q&=u/s,\quad d=1-2q=\frac{K\tanh(A/2)}{2r+K}.\end{aligned}
\]

The unfailed stationary joint distribution is [(1-q)/2,q/2,q/2,(1-q)/2] in states 00,01,10,11. Its mutual information is

\[
I=\ln 2+q\ln q+(1-q)\ln(1-q).
\]

Net memory heat to the bath is A[a q-b(1-q)]. Stationarity gives a q-b(1-q)=r(1-2q), so

\[
P=Ard.
\]

This equals net external work, because the mean internal energy is stationary. In physical units multiply P by k_B T_bath. The total steady-state entropy production of the driven four-state process is P in k_B/time units. The two-state Z chain alone is reversible: computing entropy production after aggregating its distinct signal and memory channels would incorrectly erase the dissipation.

The memory's learning rate, as a contribution to the mutual-information derivative, is

\[
\ell=rd\ln((1-q)/q)=2rd\operatorname{atanh}(d)\le P.
\]

The inequality follows because A > 2 atanh(d) for r>0, finite K, and d>0. ell is not transfer entropy, and it is not the stationary derivative of total mutual information, which is zero. These are elementary specializations of established bipartite information thermodynamics [L1,L2,L5], not claimed discoveries.

## 3. Exact kinetic-ceiling frontier

For required accuracy at least 1-q_star with 0<q_star<1/2, put d_star=1-2q_star. Restrict to 0<K<=Kmax and finite A. Feasibility requires

\[
d_\star<K_{\max}/(2r+K_{\max}).
\]

At equality, only an infinite energy gap would reach the target. Below the threshold, exact target accuracy requires

\[
A(K)=2\operatorname{atanh}[d_\star(1+2r/K)].
\]

**Proposition 1.** In this family, minimum stationary sensor work to reach the target under K<=Kmax is

\[
P_{\min}=2rd_\star\operatorname{atanh}[d_\star(1+2r/K_{\max})].
\]

**Proof.** For fixed d>0 the atanh argument, hence P, strictly decreases as K increases. For fixed K, d atanh[d(1+2r/K)] strictly increases with d on its feasible domain. Therefore select d=d_star and K=Kmax. Substitution proves the formula. The displayed rates construct a sensor achieving it. At q_star=1/2, A=0 gives zero cost. QED.

As Kmax->infinity, Pmin approaches 2 r d_star atanh(d_star)=ell; positive accuracy at finite environmental rate still dissipates. This is a bound-and-construction result in a two-parameter family, not a bound over arbitrary sensors, memory sizes, protocols, or devices. Hardware costs can move the optimum away from Kmax.

## 4. Failure law H: instantaneous mismatch hazard

Given a trajectory Z(t), death occurs with rate delta+h Z(t), with h>=0 and baseline delta>=0. This is an effective killing law, not a microscopically derived damage mechanism. The live-state row subgenerator is

\[
L=\begin{pmatrix}-u-\delta & u\\v & -v-\delta-h\end{pmatrix}.
\]

For any initial mismatch probability q0,

\[
S_H(T)=(1-q_0,q_0)e^{LT}(1,1)^\mathsf{T}.
\]

Unless stated otherwise q0=q (stationary **before** failure selection). Conditional distributions among survivors are generally not stationary.

Feynman-Kac gives S_H(T)=exp(-delta T) E exp[-h integral_0^T Z(t) dt]. With stationary initialization Jensen implies

\[
S_H(T)\ge e^{-(\delta+hq)T}.
\]

Thus replacing the fluctuating hazard by its stationary mean is a lower bound on survival in this model, not an equality. The difference expresses survivor selection and exposure heterogeneity; it does not establish that an initially endangered individual benefits from slow sensing.

## 5. Exact tradeoff at fixed accuracy and fixed environment

The asymptotic decay exponent is

\[
\lambda=\delta+\kappa,\qquad \kappa=\frac{s+h-\sqrt{(s+h)^2-4qsh}}{2}.
\]

For numerical stability, use kappa=2qsh/[s+h+sqrt((s+h)^2-4qsh)]. For 0<q<1, s,h>0 the slow eigenvalue is simple and 0<kappa<qh. It is an asymptotic exponent, not finite-time failure probability.

**Proposition 2.** At fixed r, 0<q<1/2, and h>0, increasing feasible K while adjusting A to maintain q strictly reduces stationary sensor power but strictly increases the asymptotic decay exponent under failure law H.

**Proof.** Proposition 1 gives dP/dK<0 at fixed q. The root satisfies kappa^2-(s+h)kappa+qsh=0. The polynomial is positive at zero and negative at qh because its value there is -q(1-q)h^2. Hence 0<kappa<qh, and implicit differentiation gives

\[
\frac{d\kappa}{ds}=\frac{qh-\kappa}{s+h-2\kappa}>0.
\]

Since ds/dK=1, the claimed ordering follows. QED.

The statistical meaning is important: slowly changing trajectories include a subset that remains matched for a long time. This subset disproportionately supplies long-lived survivors. Faster mixing makes all trajectories sample mismatch more evenly. This is NOT a theorem that energy dissipation universally improves survival. It is a conditional tradeoff derived from a particular observation/preparation protocol and damage law. Its algebra uses standard finite-state spectral theory; literature priority for this combined statement is unestablished.

An informative local expansion, as h->0 at fixed s,T, is

\[
\log S_H(T)=-\delta T-hqT+h^2q(1-q)\left[\frac{T}{s}-\frac{1-e^{-sT}}{s^2}\right]+O(h^3).
\]

This follows from Cov[Z(0),Z(t)]=q(1-q)exp(-s|t|) and the cumulant expansion of integrated exposure. Stationary information fixes q but not s: temporal correlations enter reliability already at second order. The expansion is explanatory, not used as an exact numerical solver.

## 6. Failure law D: multistage exposure with reset

Use n sequential damage stages, each advancing at z=n/tau only while mismatched. Reaching the end kills; returning to match resets damage instantly. The live states are C,D0,...,D(n-1). Rates are C->D0=u, Di->C=v, Di->D(i+1)=z, and D(n-1)->death=z. Add delta to all killing rates if required. In an uninterrupted mismatch, lethal exposure duration is Erlang(n,z), with mean tau. Reset is an idealized effective process; its work cost is not included in P.

Initialize [1-q0,q0,0,...,0]: the sensor state can be stationary while the damage clock is freshly prepared. This is not a stationary population of pre-damaged individuals. S_D(T) is obtained by the same matrix-exponential formula. When n=1 this is exactly law H with h=1/tau. This equality is a transfer check, not an approximate analogy.

## 7. Constructive finite-time reversal

Hold r=0.1, q=0.2, T=20, delta=0. Two designs:

| quantity | slow sensor | fast sensor |
|---|---:|---:|
| K | 0.4 | 5.8 |
| A | 2.9444389792 | 1.4522523289 |
| accuracy | 0.8 | 0.8 |
| mutual information (nats) | 0.1927447570 | 0.1927447570 |
| work rate (k_B T/time) | 0.1766663387 | 0.0871351397 |
| S_H(20), h=0.5 | 0.3026456742 | 0.1534855473 |
| S_D(20), n=4, tau=2 | 0.3825722145 | 0.8678821885 |

The environment is identical in both designs; only the sensor kinetics and gap differ. A fixed instantaneous information level cannot rank survival. Improved sensor thermodynamic efficiency cannot provide a task-independent survival ordering. Both failure laws have the same mean lethal uninterrupted mismatch duration of 2. Their exposure-duration distributions differ (exponential versus four-stage Erlang); equality of means does not make them the same hazard. The statement is about the within-law ordering of the same sensor pair.

This counterexample is calculated by exact small-matrix exponentials and checked independently by Gillespie trajectories. It is a candidate research synthesis, not an established discovery of new physics.

## 8. Scope of the barrier crossing

The link demonstrated is: explicit microscopic rates and energy levels -> information and response times -> path-dependent effective failure probabilities. The model exposes what stationary information leaves out and supplies an exact restricted physical cost frontier. It does not yet close the reverse loop in which energy reserves, repair, or actuator degradation change the sensing dynamics. Adding that loop, with local detailed balance throughout, is the next substantive research problem.
