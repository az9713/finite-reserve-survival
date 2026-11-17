## Main results and the logic connecting them

The results use two meanings of time: the integer round count \(N\) in the checkpoint model and continuous duration \(T\) in the autonomous model. All energetic statements specify what is inside the accounting boundary. No molecular mechanism is assumed merely because a Markov edge has a thermodynamically admissible rate ratio.

**Theorem A (finite-horizon checkpoint control).** Under the supplied-phase protocol below, the Bellman recurrence gives the supremum survival probability over all admissible history-dependent, randomized controllers. A deterministic maximizing policy and maximum-rate repair attain it. The proof is conditioning plus induction, with repair reduced to a two-state master equation.

**Theorem B (strict feedback value).** For two rounds, perfect readout, initial reserve two, unit sensing and engagement costs, capacity and harvest at least two, and IID opportunity probability \(q\in(0,1)\), adaptive sensing strictly exceeds every precommitted sensing word whenever \(0<\kappa<1\). At \(q=\kappa=1/2\) and capacity/harvest two the values are \(7/8\) and \(3/4\). With symmetric readout accuracy \(f\), this instance has a strict advantage precisely for \(f>2/3\). These constants belong to this architecture and parameter choice.

**Theorem C (all-word asymptotic comparison).** In the symmetric, perfect-readout, two-unit checkpoint model, the optimal adaptive and fixed-word survival root rates are

\[
r_{\rm ad}=\frac{1+\sqrt5}{4},\qquad
r_{\rm fixed}=\frac12\max\left\{\sqrt2,\frac{1+\sqrt{1+4\kappa}}2\right\}.
\]

The upper bounds apply even to nonperiodic words chosen separately for each horizon. Matching stationary or alternating constructions attain the rates. The proof reduces the fixed problem to products of two nonnegative matrices plus bounded additive terms. Without harvesting, finite reserves give the same root rate \(1-q(1-\rho)\) to both classes. Thus replenishment is essential to this example's exponential advantage.

**Theorem D (certified autonomous instance).** With the exact rational parameters and launch distribution specified below, one of 324 static binary circuits is uniquely optimal at \(T=4\). It senses at full reserve only. Writing \(D=288230376151711744\), its certified survival satisfies

\[
\frac{153847617132986168}{D}\le S_{\rm winner}\le
\frac{153847742532300097}{D}.
\]

Every other binary circuit has survival at most \(153763933422971466/D\). Allowing a constant sensing activity anywhere in \([0,1]\), with each actuator table reoptimized, gives the uniform comparator bound

\[
\sup_{\alpha,\,\mathrm{table}}S_{\rm uniform}(4)
\le\frac{145923831922267153}{D}
<\frac{153847617132986168}{D}.
\]

This is a certified gap of at least about 2.75 percentage points. Its proof combines nonnegative uniformization, outward integer rounding, and a rigorous interpolation remainder. The finite binary optimum is sharp within its declared class; the continuum upper bound need not be attained.

**Theorem E (the large-time objective selects differently).** Among the same binary autonomous circuits, the best reachable spectral decay exponent is attained with sensing off. The exponent is about \(0.24252027\), smaller than the best exponent with any sensing enabled. Exact positive-vector ratio bounds certify the ordering. This does not optimize all analogue circuits at infinite time.

**Theorem F (formation of the maintenance machinery).** In the closed six-species network of Study 05, the net assembly and repair counts obey \(J_A=h+d\) and \(J_A+J_R=N-b\) on every path. Reaching \(k\) complexes after \(r\) net fuel-consuming repairs therefore requires \(N\ge k+r\). For \(2\le k\le M\), an explicit assembly and damage–repair sequence attains the count. An exact first-assembly/persistence formula and a single-complex calculus theorem give additional sharp bounds and attaining static designs.

**Theorem G (a certified interior chemical design optimum).** For two pairs and four fuel molecules, formation by time one followed by four units of persistence has globally optimal probability in \([0.07486803450,0.07496800805]\). Every boundary design has probability at most \(0.07296461445\). A rational design attains a lower bound within \(10^{-4}\) of the continuum maximum. The complete parameter cover and its replay establish an interior optimum in the stated static kinetic family. Exact symbolic maximizing parameters and uniqueness remain open.

Information's survival value depends jointly on temporal error structure, resource state, fuel access, damage dynamics, and the horizon. In Study 05, the availability of assembled catalysts also determines whether cooperative maintenance can occur. The golden-ratio factor above comes from a two-state surviving reserve process. Detailed definitions and proofs follow.
