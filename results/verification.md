# Verification record

The analytic proofs are in `docs/model.md`. Numerical checks support their implementation; they do not replace the proofs or establish literature novelty.

## Independent event-driven simulations

Matrix exponentials and Gillespie paths use different numerical methods. These are survival probabilities for the stated effective models, not empirical device measurements. Intervals are 95% Wilson binomial intervals; a valid implementation need not land inside every such interval.

| Design / law | Exact | Monte Carlo | 95% interval | Samples | Seed |
|---|---:|---:|---|---:|---:|
| slow, hazard | 0.302645674 | 0.303366667 | [0.298190, 0.308594] | 30000 | 801 |
| slow, damage | 0.382572214 | 0.384133333 | [0.378645, 0.389652] | 30000 | 901 |
| fast, hazard | 0.153485547 | 0.155066667 | [0.151015, 0.159207] | 30000 | 802 |
| fast, damage | 0.867882189 | 0.868900000 | [0.865034, 0.872672] | 30000 | 902 |

## Energy and probability checks

- Maximum absolute disagreement among analytic power and independent four-state work, heat, and entropy-production calculations across 200 seeded parameter cases: 7.105e-15.
- Unit tests exercise local detailed balance, stationary distributions, full-to-reduced survival agreement, the information–dissipation inequality, the kinetic-ceiling frontier, exact damage/hazard equivalence at one stage, and numerical simulation.
- Test command: `PYTHONPATH=src python -m unittest discover -s tests -v`.
- Simulation seeds, package versions, and source hashes are recorded in `summary.json` and `run-manifest.json`.

## Transfer and counterexamples

- Ranking reversal found in 69 of 144 explicitly sampled cases. This is a grid count, not an estimated frequency in nature.
- The sweep varies accuracy, horizon, number of damage stages, and initial sensor preparation. Every case is retained in `transfer.csv`.
- Single-stage damage cannot reverse the ordering relative to the corresponding instantaneous hazard because the models coincide.
- Changing the initial condition can change conclusions based on stationary survivor selection. Fresh damage initialization is an assumption in every damage run.

## Remaining uncertainty

No experimental validation; no exhaustive novelty review; no full energy accounting for the effective repair/actuator layers; no energy-reserve feedback; and no guaranteed global optimization outside the explicit two-state sensor family. The interactive atlas is a finite sample, not an impossibility certificate.
