# Finite Reserve Survival

[**Open the live scientific dossier**](https://az9713.github.io/finite-reserve-survival/) · [Development journey](https://az9713.github.io/finite-reserve-survival/docs/DEVELOPMENT-JOURNEY.html) · [Paper draft](https://az9713.github.io/finite-reserve-survival/results/research-manuscript.pdf)

[![Scientific dossier preview](results/document-previews/scientific-dossier-1440.png)](https://az9713.github.io/finite-reserve-survival/)

Four stochastic models track how sensing, action, and repair affect the chance that a physical device keeps functioning. The strongest results are exact survival bounds with attaining policies and a certified example in which the value of sensing changes with the time horizon.

| Study | Finding | Significance |
|---|---|---|
| Thermal sensing | Two sensors with 80% stationary accuracy reverse their survival ranking when the damage mechanism changes. | Accuracy misses the timing of errors. Clustering can matter as much as frequency. |
| Paid repair | For a fixed energy quantum, the stress-count distribution gives a sharp survival frontier. With a known stress count and a fixed total reserve, equal repair quanta attain the optimal bound. | The model gives an exact energy requirement for a target survival probability within its protocol. |
| Shared reserve | Adaptive sensing achieves **87.5%** survival after two rounds versus **75%** for the best fixed sensing schedule. The advantage begins above **two-thirds** readout accuracy in this example. | A measurement can be worth paying for in a full battery state and wasteful in a depleted one. |
| Autonomous reactions | At time four, reserve-dependent sensing attains **53.38%** survival. Every constant sensing intensity in the specified comparator class is bounded above by **50.63%**. The best binary circuit for the asymptotic decay rate leaves sensing off. | Optimizing survival at a deadline can select a different architecture from optimizing the long-run decay rate. |

## Candidate contributions

- An exact long-run comparison against **every precommitted sensing schedule** in the two-unit checkpoint model, including an explicit switch in the best fixed strategy at repair success probability **2 − √2**.
- The exact **two-thirds** sensing accuracy threshold for a strict two-round feedback advantage in the symmetric instance.
- A rationally certified continuous-time benchmark that separates reserve-dependent sensing from a full continuum of constant sensing intensities, alongside a separate asymptotic comparison over binary circuits.

The [targeted literature review](https://az9713.github.io/finite-reserve-survival/research/priority-search-2026-11-17.md) found related work in costly sensing, energy harvesting, and autonomous information thermodynamics. It identified no exact match for the stated results. Literature priority and independent review remain open.

Tsimerman’s Golden Gate Institute talk motivated the project. The models connect microscopic energetics and temporal information to operational survival. The autonomous example uses a maintained fuel reservoir and designed reaction gates. A molecular realization and self-assembly remain open.

## Reproduce

Python 3.10+ with numpy, scipy, and matplotlib. PDF export additionally uses reportlab.

```bash
python -m pip install -e .
python -m pip install reportlab
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_study.py --samples 30000
python scripts/build_artifacts.py
```

`run_study.py` writes exact results, seeded event-driven simulations, 144 transfer cases, 1,875 sampled sensor designs at three horizons, figures, and environment/source provenance. `build_artifacts.py` builds the offline HTML and PDF. No external data or credentials are needed to reproduce calculations. Cached literature PDFs are excluded from git; source URLs and hashes remain available.

## What the theory does and does not establish

The physical layer satisfies local detailed balance and tracks external work and bath heat. Its power frontier is exact only within the stated sensor architecture and kinetic ceiling. Reliability is derived from explicit effective failure laws, whose damage/reset work is not included. No new law of physics, universal cost of agency, or biological validation is claimed.

The studies connect physical kinetics and energetic constraints to **temporal error structure** and task-dependent survival. Studies 02–04 add paid repair, a shared reserve, and autonomous reversible channels. A finite battery in the latter studies limits storage while an explicitly supplied external source replenishes it. Self-assembly, a chemical realization, and biological validation remain open. The dossier distinguishes each study's boundary and comparator class.

To rebuild the documentation, install `markdown` and run `python -B scripts/build_research_docs.py`. This reruns the 30 scientific tests and records source hashes. The HTML equations use the bundled, licensed KaTeX files and work offline. `scripts/build_manuscript.py` additionally uses `pypandoc_binary` and Tectonic to build the LaTeX/PDF and submission-source archive. See `results/dossier-ui-validation.json` for the browser and link checks.

