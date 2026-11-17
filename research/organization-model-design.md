# Finite-fuel tether formation and first collapse

## Benchmark and comparison class

Two preassembled repair complexes start uniformly dispersed in a unit periodic cube. One complex is healthy and one is damaged. Initially the tether is open and all four inventory molecules are fuel. The primary observable is the probability of avoiding simultaneous damage up to model time \(T=2\). These conditions and the coefficients below were fixed before computing the comparison. The problem is exploratory and dimensionless. The previous spatial model supplies the damage and repair kinetics. [Spatial maintenance model](../docs/spatial-repair.md)

A paid capture event converts one fuel to waste, closes the tether and leaves the relative position unchanged. A second, passive channel can close or open the same tether. In the fixed comparator, an arbitrary constant \(\alpha\in[0,1]\) scales both directions of the paid channel. In the candidate, a second fuel molecule must participate and comes back unchanged. Both designs have identical potential, passive binding, maximal paid capture rate, maximal reverse microscopic rate *ceilings*, exposure to fuel and waste, initial physical state and available repair chemistry. Their realized reverse rates can differ considerably. No comparison to all time-homogeneous networks or to every possible permanent spatial arrangement follows. The spectator network itself is time homogeneous.

| Quantity | Value |
|---|---:|
| Side \(L\), contact radius \(a\) | 1, 0.32 |
| Relative diffusivity \(D\), damage rate \(\eta\) | 0.025, 1 |
| Intrinsic repair coefficient \(\kappa\) | 16 |
| Inventory \(N\), initial healthy count \(h\), initial closed-tether indicator \(q\) | 4, 1, 0 |
| Passive and paid encounter coefficients \(k_0,k_1\) | 1, 4 |
| Tether offset \(\delta\), well depth parameter \(E\) in thermal units | 2, 6 |
| Primary time \(T\) | 2 |

The concentration-dependent activity models an effective concerted encounter involving a spectator. No sequence of elementary molecular recognition steps or synthesis cost for an actual engineered binding site has been established. Requiring a spectator also suppresses *reverse* paid events. The experiment therefore cannot assign a positive survival difference specifically to preservation of the final fuel molecule. Equal supplied hardware parameters are a benchmark convention. An actual spectator-binding motif may have different molecular costs.

## Motion, switching and inventory

Let \(x\in\mathbb T_L^3\) be relative position, \(r=|x|_{\rm per}\), \(h\in\{0,1,2\}\) the number of healthy complexes, \(b\in\{0,\ldots,N\}\) the fuel count, and \(q\in\{0,1\}\) the conformation. Write \(w=N-b,\ v=(1-r^2/a^2)_+^2,\ g=v^2,\ U_0=0,\ U_1=\delta-Ev\). The well is flat outside contact. Between jumps the diffusion is

\[
dx=-D\nabla U_q(x)\,dt+\sqrt{2D}\,dW.
\]

In contact, closure changes the force on relative position, even though the position at each chemical jump is preserved. Both conformations use exactly the same contact criterion, diffusivity, and repair rates. There is no geometric constraint pinning the pair to a site. In particular, “tether closed” and “inside contact” are different observables.

Thermal damage and restoration rates are \(\eta h\) and \(\eta(2-h)/2\). The contact repair edge \( (h=1,b)\to(2,b-1)\) has rate \(\kappa b g\). Its aggregated reverse at the destination has rate \(\kappa(w+1)g/2\), including two choices of partner to damage. Define \(m(u)=\min(1,e^{-u})\).

| Tether edge at fixed \(h,x\) | Forward rate | Reverse rate evaluated at the destination |
|---|---|---|
| Passive \((0,b)\leftrightarrow(1,b)\) | \(k_0g\,m(U_1)\) | \(k_0g\,m(-U_1)\) |
| Paid \((0,b)\leftrightarrow(1,b-1)\), \(b>0\) | \(k_1g\,b\,s(b)\,m(U_1)\) | \(k_1g\,(w+1)s(b)\,m(-U_1)/8\) |

For the comparator \(s(b)=\alpha\). For the candidate \(s(b)=(b-1)/(N-1)\), using the upper-fuel endpoint of the reversible edge. Its forward propensity is proportional to \(b(b-1)\), the count of ordered substrate/spectator choices. Its reverse sees exactly \(b-1\) available fuel spectators and \(w+1\) waste substrates. The corresponding stoichiometry is open \(+\,2F\leftrightarrow\) closed \(+\,F+W\). At \(b=1\) both rates on this paid edge vanish. Every permitted event has a reverse of positive propensity at its destination. Both designs can regenerate fuel.

## Thermodynamic and stochastic identities

The spatial counting measure supplies grid multiplicities on a lattice; in the continuum the reference is uniform volume. The uncensored equilibrium density is proportional to

\[
\binom2h\,2^{2-h}\binom Nb\,8^{-b}e^{-U_q(x)} .
\]

Diffusion is reversible with conditional spatial density proportional to \(e^{-U_q(x)}\). Across an active paid edge, the weight ratio is \(8b\,e^{-U_1}/(w+1)\). Dividing the displayed directional rates gives the same value. The passive edge ratio is \(e^{-U_1}\). The repair weight ratio from one healthy to two includes the lumping degeneracy and equals \(2b/(w+1)\). The shared equilibrium and shared launch imply identical initial nonequilibrium free energy \(k_BT D_{\rm KL}(p_0\|\pi)\), with the same additive equilibrium reference constant. Relaxation can redistribute the finite initial chemical free energy and dissipate it to the bath.

Let \(J_D,J_R,J_P,J_0\) be net forward damage, contact repair, paid tether and passive tether counts, respectively. Along any trajectory, including at a bounded stopping time,

\[
\Delta h=-J_D+J_R,\qquad \Delta b=-J_R-J_P,\qquad \Delta q=J_P+J_0.
\]

The state energy in thermal units is \(b\log8-(2-h)\log2+U_q(x)\). Its total change is \(-J_D\log2-J_R\log4-J_P\log8+\Delta U\); the potential change includes movement in the well and the instantaneous change in potential at each conformational event. Separating \(\Delta U\) into its continuous increments and its jump increments avoids treating installation of the potential as free external work. Factorials and health multiplicities in \(\pi\) account for entropy. Dynamics can bind passively, so tether occupancy alone cannot establish paid formation. Gross forward and reverse event counts furnish the additional evidence.

For smooth, periodic test functions \(\phi\) of the state and position, the generator is \(D(\Delta\phi-\nabla U_q\cdot\nabla\phi)+\sum_j k_j(\phi_{\mathrm{after}\ j}-\phi)\). These bounded jump propensities and the Lipschitz spatial drift on the compact torus define a nonexplosive switching diffusion. The generator martingale is \(\phi(Y_t)-\phi(Y_0)-\int_0^t\mathcal L\phi(Y_s)\,ds\). Stopping it at \(\tau=\inf\{t:h_t=0\}\) gives the killed backward equation with unit live initial payoff and zero failure boundary. This first-passage calculation allows the full chemistry to remain reversible before censoring.

## Finite grid, observables and optimization

The cubic periodic lattice has \(n^3\) positions before reflection and permutation lumping. For neighboring sites \(i,j\) the jump intensity is
\[
D (n/L)^2 \exp[-(U_q(j)-U_q(i))/2].
\]
The resulting spatial probability flux at equilibrium is symmetric on every undirected bond. The cubic symmetries preserve potentials and rates, so quotienting gives the exact spatial lattice chain for invariant observables and initial laws. This discretization changes the continuum generator. Lattice refinement can support convergence numerically but does not supply a quantitative continuum error enclosure.

For \(Q\) the killed *row* generator, \(\mu^\top e^{tQ}{\bf1}\) gives survival. Time-integrating a directional propensity gives its expected cumulative event count by compensation. In killed calculations endpoint fuel is weighted by survival and loses the residual fuel at a collapse. The only killing transition has rate \(\eta\) while \(h=1\), preserving \(b,q\). Integrating \(\eta b\mathbf1_{h=1}\) and \(\eta q\mathbf1_{h=1}\) against the live distribution records fuel and tether state at absorption. Adding these to the corresponding live moments gives \(\mathbb E[b_{t\wedge\tau}]\) and \(\mathbb E[q_{t\wedge\tau}]\) for the stopped ledgers. The implementation also reports uncensored cumulative counts and endpoint diagnostics when supplied the full generator. Tether, contact and gate observables of a killed process are unnormalized joint live moments. Division by live mass gives survivor-conditional values.

The fixed generator is affine, \(Q_\alpha=Q_0+\alpha B\). Both endpoints and all convex mixtures are valid killed generators. Iterated Duhamel expansions and positivity contraction in the sup norm bound the \(j\)th derivative of the survival objective by \(T^j\|B\|_\infty^j\). This supplies interpolation error bounds for a continuous-parameter comparator certificate. For this model \(\|B\|_\infty\le 2k_1N\). A sampled search alone supplies only a lower bound on the best fixed survival. The numerical certificate must name its particular lattice and floating-point treatment.

Two analytic controls clarify how to read this comparison.

1. From the primary launch, the first two survival derivatives are independent of paid tether conductance. For a live \(h=1\) state the second backward-generator action on the live payoff is \(\eta(3\eta/2+\kappa b g)\). A paid closure lowers \(b\) by one at fixed \(x\). Applying the generator once more gives \(S'''_\alpha(0)-S'''_0(0)=-\eta\kappa k_1N\alpha\langle g^2m(U_1)\rangle\). The spectator shares this derivative with \(\alpha=1\). Early paid closure depletes immediately available repair fuel. Its spatial drift can change contact exposure later.
2. In a separate ideal reservoir model with fixed dimensionless fuel and waste activities \(f,w>0\), the fixed paid forward/reverse coefficients are \(k_1\alpha f,k_1\alpha w/8\). Macroscopic mass action for the spectator route gives \(k_1 f^2/(N-1),k_1fw/[8(N-1)]\), before applying the same geometric gate and acceptance factors. Choosing \(f=N-1,w=1\) makes the two generators identical at \(\alpha=1\), hence every path observable has the same law. This diagnostic removes fuel fluctuations and feedback together and changes absolute rates from the finite-batch launch. Replacing the finite-batch counts mechanically by \(b=N,w=1\) does not give this identity: the reverse spectator sees \(b\) possible fuel particles, while the forward spectator sees \(b-1\) after choosing its substrate.

Chebyshev expansion projects the reversible semigroup using the full equilibrium restricted to represented states as its inner-product weight. The truncation estimate in the model code uses the symmetric Skellam tail and Cauchy–Schwarz. It is an exact-arithmetic series-truncation estimate. It does not account for floating-point arithmetic, for quadrature in accumulated event observables, or for spatial discretization. [Model implementation](../src/purposeful/organization.py)
