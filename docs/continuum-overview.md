## Brownian continuum proofs and the unresolved time-two comparison

The continuum analysis proves failure of the spectator strategy at short and long horizons in the original model. At the prespecified horizon \(T=2\), the sign against the best constant paid intensity remains unresolved by a rigorous continuum error bound. [Full continuum proofs, calculations and receipts](CONTINUUM-ORGANIZATION.html)

For \(0<t\le 1/1250=0.0008\),

\[
S_\star(t)\le S_0(t)-t^3/60,\qquad
S_\alpha(t)\le S_0(t)-\alpha t^3/60.
\]

Thus passive-only organization is the unique best constant intensity throughout this interval. The proof uses a bounded-perturbation Duhamel formula, a first-jump coupling and exactly replayable rational constants.

The passive generator has the positive ground state \(\phi=r^{h+b}\), with \(r=(\sqrt{17}-1)/2\), and exact decay rate \(\lambda_0=(7-\sqrt{17})/4\). Repair preserves \(h+b\). Paid mode changes add a positive reversible Dirichlet energy on this eigenfunction, strictly increasing the spectator decay rate. Consequently \(S_\star(t)/S_0(t)\to0\). A conservative bound guarantees loss for \(t\ge4\times10^{12}\). That is a sufficient witness, with no claim that it approximates the actual crossover.

![Continuum theorem, independent Fourier checks and remaining error budget](../results/figures/continuum-organization-bounds.png)

At \(T=2\), two independent Fourier formulations estimate a spectator advantage near 0.0147 percentage points over constant intensity one. The earlier finite-grid certificate covers the optimized constant family. The Fourier calculations still require validated spatial and temporal integration and a bound over that entire family. The rigorous Euler bias bound for the finest previously simulated step is 0.239 per policy. Transferring the lattice sign through equal per-policy bounds requires less than approximately 0.000073.

A weighted primal/dual residual identity and a parameter interpolation estimate provide a route to a tighter certificate. The retained residual measurements are numerical diagnostics. The receipt explicitly records `INCONCLUSIVE_CONTINUUM_SIGN` at \(T=2\).

[Short-time rational checks](../results/continuum-short-time.json) · [Path-law and long-time proofs](../research/continuum-coupling.md) · [Residual error theorem](../research/continuum-residual.md) · [Fourier derivation](../research/continuum-spectral.md) · [Verification receipt](../results/continuum-organization-verification.json) · [Full-size figure](../results/figures/continuum-organization-bounds.pdf)
