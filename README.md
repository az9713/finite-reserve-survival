# Finite Reserve Survival

[**Open the live scientific dossier**](https://az9713.github.io/finite-reserve-survival/) · [Development journey](https://az9713.github.io/finite-reserve-survival/docs/DEVELOPMENT-JOURNEY.html) · [Paper draft](https://az9713.github.io/finite-reserve-survival/results/research-manuscript.pdf)

[![Scientific dossier preview](results/document-previews/scientific-dossier-1440.png)](https://az9713.github.io/finite-reserve-survival/)

Five stochastic models track how sensing, action, and repair affect the chance that a physical device keeps functioning. The strongest results are exact survival bounds with attaining policies and a certified example in which the value of sensing changes with the time horizon.

| Study | Finding | Significance |
|---|---|---|
| Thermal sensing | Two sensors with 80% stationary accuracy reverse their survival ranking when the damage mechanism changes. | Accuracy misses the timing of errors. Clustering can matter as much as frequency. |
| Paid repair | For a fixed energy quantum, the stress-count distribution gives a sharp survival frontier. With a known stress count and a fixed total reserve, equal repair quanta attain the optimal bound. | The model gives an exact energy requirement for a target survival probability within its protocol. |
| Shared reserve | Adaptive sensing achieves **87.5%** survival after two rounds versus **75%** for the best fixed sensing schedule. The advantage begins above **two-thirds** readout accuracy in this example. | A measurement can be worth paying for in a full battery state and wasteful in a depleted one. |
| Autonomous reactions | At time four, reserve-dependent sensing attains **53.38%** survival. Every constant sensing intensity in the specified comparator class is bounded above by **50.63%**. The best binary circuit for the asymptotic decay rate leaves sensing off. | Optimizing survival at a deadline can select a different architecture from optimizing the long-run decay rate. |
| Self-assembled repair | A closed six-species network constructs active complexes from monomers and uses them to repair damaged complexes. Its globally optimized formation-and-maintenance probability lies between **7.486803% and 7.496801%**, with every boundary design below **7.296462%**. | An intermediate binding strength and a split between assembly and repair capacity are provably optimal within this static kinetic family. |

## Candidate contributions

- A sharp net-resource bound **N ≥ k + r** for forming k complexes and performing r net fuel-consuming repairs, with an attaining path for k ≥ 2.
- An analytic sharp allocation optimum for one complex, and a **globally certified interior optimum** for two complexes over a continuous two-parameter design square. A rational construction comes within **0.01 percentage points** of the optimum.

- An exact long-run comparison against **every precommitted sensing schedule** in the two-unit checkpoint model, including an explicit switch in the best fixed strategy at repair success probability **2 − √2**.
- The exact **two-thirds** sensing accuracy threshold for a strict two-round feedback advantage in the symmetric instance.
- A rationally certified continuous-time benchmark that separates reserve-dependent sensing from a full continuum of constant sensing intensities, alongside a separate asymptotic comparison over binary circuits.

The [targeted literature review](https://az9713.github.io/finite-reserve-survival/research/priority-search-2026-11-17.md) found related work in costly sensing, energy harvesting, and autonomous information thermodynamics. It identified no exact match for the stated results. Literature priority and independent review remain open.

Tsimerman’s Golden Gate Institute talk motivated the project. The models connect microscopic energetics and temporal information to operational survival. The autonomous example uses a maintained fuel reservoir and designed reaction gates. Study 05 adds an abstract chemical self-assembly and maintenance model. A molecular realization, replication, and heredity remain open.

## Self-assembly study and evidence

[Read the complete study](https://az9713.github.io/finite-reserve-survival/docs/SELF-ASSEMBLY.html) · [Global certificate](results/assembly.json) · [Certificate replay](results/assembly-certificate-check.json) · [Independent simulation](results/assembly-evidence.json) · [Literature comparison](research/assembly-literature.md)

The assembly-and-maintenance result uses two copies of each complementary monomer, four fuel molecules, formation by time one, and four time units of operation after formation. The bounds cover every static design in the specified square. The independent simulation uses 30,000 attempts. Every molecule remains inside the compartment, with a heat bath providing thermal fluctuations.

[![Assembly and maintenance tradeoff](results/figures/assembly-frontier.png)](https://az9713.github.io/finite-reserve-survival/docs/SELF-ASSEMBLY.html)

Intermediate binding optima and fuel-driven chemical persistence have earlier theoretical and experimental precedents. The candidate contribution is the particular finite-inventory reliability model, its constructive resource bound, and the replayable continuum certificate. Priority remains open.

## Reproduce

Python 3.10+ with numpy, scipy, and matplotlib. PDF export additionally uses reportlab.

```bash
python -m pip install -e .
python -m pip install reportlab
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_study.py --samples 30000
python scripts/build_artifacts.py
python -B scripts/run_assembly.py
python -B scripts/assembly_evidence.py
python -B scripts/check_assembly_certificate.py
```

`run_study.py` writes exact results, seeded event-driven simulations, 144 transfer cases, 1,875 sampled sensor designs at three horizons, figures, and environment/source provenance. `build_artifacts.py` builds the offline HTML and PDF. No external data or credentials are needed to reproduce calculations. Cached literature PDFs are excluded from git; source URLs and hashes remain available.

## What the theory does and does not establish

Study 01 tracks sensor work and heat under local detailed balance, with a power frontier restricted to its architecture and kinetic ceiling. Its effective damage/reset layer omits repair work. Later studies supply explicit repair and resource models. Study 05 has a closed finite-inventory species model with an exact equilibrium weight, a pathwise ledger, and a declared first-passage failure rule. Molecular realization and biological validation remain open.

The studies connect physical kinetics and energetic constraints to **temporal error structure** and task-dependent survival. Studies 02–04 add paid repair, a shared reserve, and autonomous reversible channels. A finite battery in the latter studies limits storage while an explicitly supplied external source replenishes it. Study 05 adds self-assembly with finite chemical inventories and exact detailed balance. Molecular realization and biological validation remain open. The dossier distinguishes each study's boundary and comparator class.

To rebuild the documentation, install `markdown` and run `python -B scripts/build_research_docs.py`. This reruns the 38 scientific tests and records source hashes. The HTML equations use the bundled, licensed KaTeX files and work offline. `scripts/build_manuscript.py` additionally uses `pypandoc_binary` and Tectonic to build the LaTeX/PDF and submission-source archive. See `results/dossier-ui-validation.json` for the browser and link checks.
