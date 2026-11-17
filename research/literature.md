# Literature audit and novelty boundary

Research date: 2026-11-16. This is a targeted initial audit, not an exhaustive priority search.

## Sources inspected

The arXiv abstracts and downloaded PDF text were inspected. Relevant sections are identified below; this does not assert independent verification of every result in each paper. Local PDF/text caches are ignored by git and excluded from the release archive. URLs and SHA-256 hashes are retained in `source-manifest.json`.

| ID | Reference | Relevant content and consequence |
|---|---|---|
| L1 | Hartich, Barato & Seifert (2014), [Stochastic thermodynamics of bipartite systems](https://arxiv.org/abs/1402.0419), DOI [10.1088/1742-5468/2014/02/P02016](https://doi.org/10.1088/1742-5468/2014/02/P02016) | Stationary bipartite systems, information flows, transfer-entropy inequalities, and four-state illustrations. Generic unification of information and thermodynamics is established. |
| L2 | Horowitz & Esposito (2014), [Thermodynamics with continuous information flow](https://arxiv.org/abs/1402.3276), DOI [10.1103/PhysRevX.4.031015](https://doi.org/10.1103/PhysRevX.4.031015) | Sections II–III: local detailed balance and subsystem entropy/information balances. Basis for keeping physical channels separate. |
| L3 | Still, Sivak, Bell & Crooks (2012), [The thermodynamics of prediction](https://arxiv.org/abs/1203.3271), DOI [10.1103/PhysRevLett.109.120604](https://doi.org/10.1103/PhysRevLett.109.120604) | Nonpredictive memory and dissipation, particularly the work-step relationship and total-dissipation bound. A claim that predictive information has energetic significance would not be new. |
| L4 | Kolchinsky & Wolpert (2018), [Semantic information, autonomous agency, and nonequilibrium statistical physics](https://arxiv.org/abs/1806.08053), DOI [10.1098/rsfs.2018.0041](https://doi.org/10.1098/rsfs.2018.0041) | Sections II–V: causal intervention, viability functions, syntactic/semantic information and information–viability curves. The information-to-existence bridge is already explicit here; our failure probabilities use a specified operational viability objective. |
| L5 | Hartich, Barato & Seifert (2016), [Sensory capacity](https://arxiv.org/abs/1509.02111), DOI [10.1103/PhysRevE.93.022116](https://doi.org/10.1103/PhysRevE.93.022116) | Sensory capacity, learning rates, memory and energetic efficiency. We do not identify learning rate, mutual information, transfer entropy, and sensory capacity with one another. |

## Search procedure and limitations

Crossref bibliographic queries were used to confirm titles and DOIs, then arXiv PDFs were retrieved and converted with pdftotext. Queries included the exact source titles and broader phrases `thermodynamic cost sensing survival first passage reliability` and `learning rate sensory capacity energy speed accuracy tradeoff`. The broad searches were noisy and returned largely irrelevant results. They do NOT establish an absence of prior work. One earlier unrelated Crossref request in this chat was rate-limited; the queries recorded here succeeded.

Additional adjacent titles surfaced: *Efficiency of cellular information processing* (10.1088/1367-2630/16/10/103024) and *Information-theoretic versus thermodynamic entropy production in autonomous sensory networks* (10.1103/PhysRevE.87.042104). These need full-text examination before publication. Further search should cover Markov-modulated killing, dynamic disorder, occupation-time distributions, intermittent exposure toxicology, thermodynamic speed limits, and risk-sensitive control. The matrix-exponential and Feynman–Kac methods themselves are standard.

## What this package can honestly claim

1. It reproduces established local-detailed-balance and information–dissipation relationships in a transparent sensor family.
2. It derives and verifies a kinetic-ceiling accuracy–power frontier within that family.
3. It proves a conditional asymptotic survival/power ordering and constructs a finite-time reversal when the error-to-damage mechanism changes.
4. It supplies reproducible numerical evidence and an interactive comparison.

Items 2–3 are derivations performed for this project. That is not evidence that their results are first in the literature. Their combination may be useful as a focused research note or benchmark; independent expert review and a broader novelty search remain necessary.

## Hypotheses explicitly rejected or deferred

- Rejected: average accuracy alone determines reliability.
- Rejected: higher learning efficiency necessarily means better survival under every task.
- Rejected: the coarse mismatch chain's zero entropy production means the sensor is thermodynamically free. Coarse-graining loses physical channels.
- Deferred: costs of arbitrary sensing or computation are bounded by this model's frontier.
- Deferred: the effective damage clock is a complete thermodynamic repair model.
- Deferred: this sensor constitutes an autonomous agent or derives biological purpose.

## Next novelty gate

Before claiming publication-level novelty: search adjacent reliability/occupation-time literature; examine L5 and the two adjacent sensory papers in depth; request review of the proposed joint ordering; and test whether the result survives an explicit energy reserve and reversible repair chemistry. A failure at this gate would still leave a useful reproduction and counterexample package.
