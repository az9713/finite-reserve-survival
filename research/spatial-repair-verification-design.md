# Independent verification of the two-complex spatial repair pilot

## Equilibrium and energy

Let \(B=b+w\) be the conserved resource count, and let \(h+d=2\) for the labelled complexes. In dimensionless energy units, choose \(E_H=E_W=0\), \(E_D=-\log 2\), and \(E_F=\log 8\). The equilibrium weight of a labelled chemical configuration is

\[
W(s)=\frac{2^d8^{-b}}{b!w!}.
\]

The normalizer is \(9(9/8)^B/B!\). The two labels are independent, with probability \(2/3\) of damage. Independently, \(b\) has a binomial law with parameters \(B,1/9\). The relative position has uniform torus measure. In an aggregated count distribution the label multiplicity \(\binom{2}{h}\) must be included. In a symmetry-reduced spatial grid the probability of each orbit includes its number of cells.

Along one damage pair the equilibrium weight ratio is 2 and so is the forward-to-reverse rate ratio. Along one repair pair both ratios are \(4b/(w+1)\). The same contact activity multiplies the two rates at each position, so the detailed-balance flux identity holds also at positions with zero activity. There are two distinct reverse repair channels from \(HH\), one for each choice of which complex becomes damaged. The unreduced diffusion edge rates are symmetric, and aggregation of such edges preserves equilibrium flux balance when orbit weights are used.

Let \(M_d\) and \(M_r\) be the signed counts of forward damage and forward repair, with a reverse event contributing minus one to its own count. Then

\[
\Delta h=-M_d+M_r,\qquad \Delta b=-M_r,\qquad
\Delta E=-M_d\log 2-M_r\log 4.
\]

In particular, \(\Delta(h+b)=-M_d\). Any path still functional with at least one healthy complex satisfies \(M_d\le h_0+b_0-1\). This is a bound on signed damage counts. Reversible cycles can produce arbitrarily large gross event counts. The factorial terms in \(W\) contribute counting entropy to the rate ratio, while \(E\) lists species energies.

## Failure and preparations

Failure means first entry into \(DD\). The killed generator retains the diagonal exit rate from \(HD\) and \(DH\) into \(DD\), with no onward transitions. The full chain includes returns from \(DD\) and has equilibrium \(\Pr(DD)=4/9\). Terminal mass on the live states of the full chain therefore measures occupancy rather than first-failure survival.

For a channel observed in the stationary process, the event-conditioned spatial density immediately before and after a repair or reverse-repair event is proportional to \(g(r)\), since both rates carry this factor and the event leaves the positions unchanged. The \(HD,B\) event-conditioned start can be prepared by observing a stationary reverse repair from \(HH,B-1\). An \(HH,B\) start under the same spatial density is a control preparation, since a forward repair would require \(B+1\) fuel molecules beforehand. Sampling an event after a specified transient history generally gives another spatial law. Each initialization comparison holds chemical state and initial fuel inventory fixed.

At any time \(S'(t)=-\eta\Pr(\tau>t,h_t=1)\) for a live start. The difference of two survival curves with the same \(\eta\) is therefore \(\eta\)-Lipschitz. On a mesh that includes both endpoints, the sampled maximum absolute difference lies at most \(\eta\) times half the maximum gap below the continuous-time maximum. This controls temporal sampling separately from spatial discretization and Monte Carlo uncertainty.

## Count calibration

Fix an \(HD,b\) chemical state and turn off other jumps for calibration. For the spatial diffusion generator \(A\), clamp the contact hazard to \(kb\,g\). Its mean reaction time \(m\) solves \((-A+kb\,\mathrm{diag}(g))m=\mathbf 1\). Average over the uniform spatial law. For \(b>0\), set the count forward rate to \(1/\overline m_b\), then obtain its reverse rate for each labelled orientation by multiplying by \((w+1)/(4b)\). Calibration is imposed separately for each positive fuel count or by an explicitly declared protocol. The calibrated reverse rate follows from the equilibrium ratio. Matching forward and reverse clamped waiting times independently generally gives two different multipliers. A reverse-reaction clock from \(HH\) has twice the per-orientation contact hazard when both channels are enabled.

For comparison, the fast-mixing rate of a single labelled forward channel is \(kb\langle g\rangle\). Averaging changes the clock shape as well as the mean. A finite-diffusion comparison uses the clamped mean calibration and can display the fast-mixing comparison alongside it. The mean-time estimate must remain separate from a fit to the pilot's actual survival curve.

## Independent numerical checks

- Enumerate labelled chemical jump channels independently for a small full periodic grid and compare a dense matrix exponential with the symmetry-reduced solver. Include positions at a wrapped contact, the reflecting fold edges, and the post-failure killing flux. Check finite-state per-edge flux balance under the exact orbit-weighted product law and exercise every type of energy and inventory change.
- Check separately that the finite-grid spatial solution approaches the exact stationary-average count generator when the diffusivity increases. Recompute a clamped-reaction mean by a small independent dense solve and check the intended effective forward and reverse rates.
- Use exact continuum candidate-event thinning as an independent cross-check for the periodic Brownian jump process. Between candidates, draw the relative displacement from its wrapped Gaussian law with variance \(2D_{\rm rel}\Delta t\) per coordinate. At the candidate, use an independent uniform mark against the state-specific aggregate bound and the actual channel rates. Because the bound dominates for all positions in the current chemical state, candidate thinning requires neither detection of contact entry nor discretization of elapsed time. Update position at rejected candidates as well. Stop at first \(DD\). At the prespecified comparison times, report pointwise Wilson intervals and compare with a 99% simultaneous Dvoretzky–Kiefer–Wolfowitz confidence band over time and four prespecified initial preparations. Fixed seeds enable replay. Report successive-grid differences separately from Monte Carlo sampling uncertainty.
- In the separate radial benchmark, take Doi hazard \(\kappa\) throughout the ball \(r\le a\) and reflection at \(R\). Let \(T_{\rm abs}\) be the Stage 1 absorbing-target mean for \(r>a\). Solving inside and matching mean and radial derivative at \(a\) gives

  \[
  T(r)=T_{\rm abs}(r)+\frac1\kappa+
  \frac{R^3-a^3}{3Da^2\,[\mu\coth(\mu a)-1/a]},\qquad \mu=\sqrt{\kappa/D}.
  \]

  A separate finite-volume integral balance on radial shells checks this mean at finite \(\kappa\). The exterior mean tends to the Stage 1 hard-capture mean as \(\kappa\to\infty\). This check concerns a reflecting sphere with a step hazard and is distinct from grid convergence in the periodic cube with a smooth gate.
