## Appendix — Paid repair contact and its sharp energy bound

Let a repair pair have damaged probability \(1-p(t)\), healthy probability \(p(t)\), energy drop \(\epsilon k_BT\), and sum of directional rates \(\gamma(t)\in[0,\Gamma]\). Fixed local detailed balance implies

\[
k_+(t)=\gamma(t)\sigma(\epsilon),\quad
k_-(t)=\gamma(t)[1-\sigma(\epsilon)],\quad
\sigma(x)=\frac1{1+e^{-x}}.
\]

With a damaged start, \(p(0)=0\) and the master equation is
\(p'=\gamma(t)[\sigma(\epsilon)-p]\). Integration gives

\[
p(\tau)=\sigma(\epsilon)\left[1-e^{-\int_0^\tau\gamma(t)\,dt}\right]
\le(1-e^{-\Gamma\tau})\sigma(\epsilon)=\kappa.
\]

Equality holds by taking \(\gamma\equiv\Gamma\). The protocol may choose the barrier schedule before contact but cannot latch the first repair using unpriced within-contact feedback. On a live outcome, the net battery expenditure is exactly one unit, irrespective of back-and-forth jumps. This proves the repair bound used by Theorem A.

For exogenous stresses, let \(M\) count the number of stressed checkpoints in \(N\) rounds, continued mathematically even beyond a device failure. With \(B\) identical reserve units, sequential conditioning bounds survival on every flag path by \(\kappa^M\mathbf1_{M\le B}\). Maximum-rate contact attains each factor, so

\[
S^*(B)=\sum_{m=0}^{\min(B,N)}\Pr(M=m)\kappa^m.
\]

This argument allows arbitrarily correlated stresses. With \(m\ge1\) certain stresses and variable positive battery gaps whose sum is at most \(F\), put \(a=1-e^{-\Gamma\tau}\). Since
\((\log\sigma)''=-\sigma(1-\sigma)<0\), Jensen's inequality yields

\[
\Pr(\mathrm{survive}\mid M=m)
\le a^m\prod_{i=1}^m\sigma(\epsilon_i)
\le [a\sigma(F/m)]^m.
\]

For \(F>0\), exactly \(m\) equal gaps and maximum-rate contact attain equality. This is a conditional allocation theorem for a known stress count. It does not optimize an unknown count with arbitrary gaps. The clock, load gate, and initial preparations remain supplied resources outside this reserve bound.
