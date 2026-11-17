# Study 02 — Sharp checkpoint survival with finite reserves and paid repair

Study 01 motivates a new checkpoint model with a finite battery and a
reversible repair reaction. This study proves an attainable survival
frontier for that class of protocols and shows that battery size alone can
reverse a survival ranking. All energies below are in units of k_B T.

## Physical model and protocol

A device has a binary stress flag Z, a damage bit D, and an integer battery
level b between 0 and B. The initial battery level is B and the device starts
healthy (D = 0). A round has an exposure interval of length Δ, a load operation,
and a repair contact of length τ. Failure occurs at the first end-of-contact
checkpoint with D = 1. A healthy device with zero remaining battery stays
operational until it needs a repair.

- During exposure, Z evolves independently of the repair device. For a physical
  realization, use a thermal two-state switch with transition rates 0 → 1 at
  q s and 1 → 0 at (1 − q)s, and energy difference E₁ − E₀ = ln((1 − q)/q).
  These rates satisfy local detailed balance. Start with P(Z = 1) = q. The
  sampled chain has correlation multiplier ρ = exp(−sΔ). It has exactly the same
  sampled mismatch law as the stationary Study 01 sensor at the same q and s.
  The thermal switch is one realization of the flag. Using the original driven
  sensor instead would require supplying its work separately.
- Freeze Z during loading and repair. The load operation applies the reversible
  permutation (Z,D,b) → (Z,D XOR Z,b) to a degenerate damage bit. On a surviving
  path this creates damage precisely at stressed checkpoints. The operation
  is an ideal prescribed, energy-neutral gate in this model. The flag is
  preserved and no damage bit is erased at a checkpoint.
- Only when Z = 1, the repair contact opens the pairs
  (D = 1, b) ↔ (D = 0, b − 1) for b ≥ 1. A forward repair lowers battery energy
  by ε and releases ε to the bath. Its reverse recreates damage, restores one
  unit to the battery, and absorbs ε from the bath. For jumps during contact,
  ΔE = −Q_bath, and log(k_forward/k_reverse) = ε. If the quanta vary with battery
  level, replace ε with the energy gap for the pair. The isolated damaged
  state at b = 0 has no repair transition.

An admissible protocol changes the sum γ(t) of the forward and reverse repair
rates in each pair, subject to 0 ≤ γ(t) ≤ Γ. The rate ratio stays fixed by the
energy gap. Both states of a pair share the same prescribed barrier schedule.
The schedule can depend on the pair and on information available at the start
of a round. Within-contact feedback that latches the first successful repair
would require additional states or a controller and lies outside the class.
Phase timing and barrier switching are ideal external operations. The energy
bound is for the initial battery above its empty state, with the device and
clock supplied. Failure is a first-passage event in the analysis and requires
no physical absorbing jump or erasure.

## Theorem 1 — Exact survival frontier

With a fixed quantum ε > 0, put

\[
a=1-e^{-\Gamma\tau},\qquad \kappa=a\sigma(\epsilon),\qquad \sigma(x)=\frac1{1+e^{-x}}.
\]

Let M count the stressed checkpoints among N rounds of the exogenous flag
trajectory, including its continuation after any device failure. Among all
admissible protocols, the maximum probability of surviving all N rounds with
B initial units is

\[
S^*(B)=\sum_{m=0}^{\min(B,N)}\Pr(M=m)\kappa^m,
\]

where the sum stops at N if B exceeds N. Opening every repair pair at the full
rate ceiling throughout each contact attains this bound simultaneously for
every B and N. The minimum initial energy for target survival η is ε times
the smallest B with S*(B) ≥ η. If S*(N) < η, additional units of the same
quantum cannot reach the target.

**Proof.** Starting damaged in a fueled pair, let p(t) be the probability of
being healthy at time t. Local detailed balance gives forward and reverse
rates γ(t)σ(ε) and γ(t)[1 − σ(ε)]. Thus p′ = γ(t)[σ(ε) − p] and p(0) = 0.
Integration gives

\[
p(\tau)=\sigma(\epsilon)\left[1-e^{-\int_0^\tau\gamma(t)\,dt}\right]\le\kappa.
\]

At a stressed checkpoint, ending healthy spends exactly one battery unit,
even if many forward and reverse jumps occurred. Unstressed checkpoints leave
the device healthy without spending a unit. Condition on a complete flag
path with M = m. Sequential conditioning on survival at earlier checkpoints
bounds its probability of survival by κ^m if m ≤ B, or zero if m > B.
Under the constant maximum schedule every factor equals κ. Averaging over
flag paths proves the formula and the matching construction. This argument
uses no independence assumption between successive stress flags.

## Theorem 2 — Sharp energy when the stress count is fixed

Allow arbitrary positive battery gaps and let F be the initial battery
energy above empty. Conditional on m ≥ 1 stressed checkpoints, every
admissible protocol satisfies

\[
\Pr(\mathrm{survive}\mid M=m)\le[a\sigma(F/m)]^m.
\]

When F > 0, a battery with exactly m equal gaps F/m and the maximum-rate
schedule attains the upper bound. Thus unequal allocation of the same
stored energy cannot improve survival for a fixed number of repairs.

**Proof.** When m repair pairs are available, let ε₁, …, εₘ be the gaps
of their successive forward transitions. Their sum is at most F. The
conditional product is bounded by a^m ∏ σ(εᵢ). The function log σ(x) is
increasing and strictly concave, with second derivative −σ(x)[1 − σ(x)].
Jensen's inequality gives the bound, with equality for m equal gaps
that exhaust F. If fewer than m pairs are available, survival is zero.

For a target η above (a/2)^m and below a^m, the exact minimum F is
m log[p/(1 − p)], where p = η^(1/m)/a. At or below (a/2)^m the infimum
is zero when strictly paid transitions are required. At or above a^m the
target is infeasible with any finite battery. For example, three certain
stresses, Γ = 2, and τ = 1 require 5.2542045 k_B T to reach 40% survival.
The construction uses three equal quanta of about 1.7514015 k_B T.

The energy theorem conditions on a fixed stress count. It does not solve
the allocation problem for arbitrary gaps when the count is unknown.
Theorem 1 gives an unconditional frontier for fixed-quantum batteries.

## A reserve-dependent ranking reversal

For N = 2, the count distribution at fixed stationary q has a simple
dependence on the checkpoint correlation ρ:

\[
\begin{aligned}p_0&=(1-q)^2+q(1-q)\rho,\\p_1&=2q(1-q)(1-\rho),\\p_2&=q^2+q(1-q)\rho.\end{aligned}
\]

Hence the derivatives of S*(B) with respect to ρ for B = 0, 1, and B ≥ 2
are q(1 − q) times 1, 1 − 2κ, and (1 − κ)², respectively. Whenever
κ > 1/2, faster mixing lowers optimal survival with zero battery, raises
it with one unit, and lowers it again with at least two units. With one
unit, two clustered stresses exhaust the reserve, and one isolated stress
can be repaired. This is an analytic ordering for 0 < q < 1 and κ > 1/2
under the two-round protocol.

Take q = 0.2, Δ = τ = 1, ε = ln 9, and Γ = ln 10, which give κ = 0.81.
Use the two mixing rates from Study 01, s = 0.6 and s = 6:

| Initial battery | Slow stress dynamics | Fast stress dynamics |
|---|---:|---:|
| 0 units | 72.78% | 64.04% |
| 1 unit | 84.48% | 89.90% |
| 2 units | 92.86% | 92.55% |

The fast design reaches 88% survival with one unit. The slow design needs
two. Both have the same marginal stress probability and identical repair
kinetics, energy per repair, and checkpoint protocol.

## Reproduce and scope

```bash
python scripts/run_reserve_study.py
```

The tests compare the count-based formula against a separate full-state matrix
exponential on (stress, damage, reserve), including both directions of every
repair edge and the projection at each checkpoint. They also check local
detailed balance, a throttled schedule, target inversion, and the ranking
reversal. The phase protocol is driven. The constructive result leaves open
the autonomous case where sensing and actuation draw from the same reserve,
action changes fuel access, and operation can fail between checkpoints.
Literature priority for these results has not been established.

Background on the thermal jump accounting appears in [Seifert's review](https://arxiv.org/abs/1205.4176).
