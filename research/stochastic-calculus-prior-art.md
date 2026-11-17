# Prior-work and parameter audit for the stochastic-calculus implementation plan

Date 17 November 2026. This receipt supports the implementation plan and its classical back-of-the-envelope calculation. It records inspected sources and unresolved comparisons. It makes no claim that diffusion-limited repair or its stochastic-calculus formulation is new.

## Claim-by-claim assessment

| Proposed element | Closest prior formulation | Assessment for this project |
|---|---|---|
| Poisson reaction-count representation, generator martingales | Anderson and Kurtz's chemical-network survey | Established machinery. Reuse to couple reaction jumps and positions. |
| Diffusing particle reaching a central target in a reflecting sphere | Agbanusi and Isaacson, 2013, Section 3 and Appendix A | The geometry is already explicit. Our mean-time and Laplace calculations are a classical benchmark reproduction. |
| Partially reactive target and diffusion-renormalized rate | Collins–Kimball theory, plus Berezhkovskii–Szabo, 2013 | Established. Effective forward and reverse rates require a consistent equilibrium ratio. |
| Validity of count dynamics and voxel approximations | Gillespie, Petzold and Seitaridou, 2014 | Established physical validity conditions. A large target-search time alone does not establish a failure of a calibrated count model. |
| Uniform-in-time small-target approximations | Isaacson, Mauro and Newby, 2016 | Already developed for Doi and partial absorption models. A generic request for a small-target approximation is insufficient as an originality claim. |
| Spatial detailed balance under association and dissociation | Zhang and Isaacson, 2022, Sections III–IV | Already formulated in detail. Binding/unbinding placement rules must obey a pointwise flux relation, including boundaries. |
| Efficient exact first-passage simulation | Oppelstrup et al., 2009 | Established method. Prefer reuse over inventing a new simulator without need. |
| Stop-time thermodynamic bounds | Li and Chou, 2026 preprint | Direct overlap with a proposed martingale route. Full theorem and assumptions comparison is required before any originality claim. |
| Finite-fuel depletion coupled to repair-generated spatial correlations | No complete theorem comparison performed in this pass | A question to audit. Neither novelty nor absence of prior work has been established. |

The deliberately narrow first experiment is a matched-rate comparison between count and spatial maintenance dynamics, using uniform versus localized initial conditions. Its output can be a validity finding or a reproduction. A substantive new claim must be identified before undertaking a broad theorem project.

## Sources and access depth

1. M. B. Elowitz, M. G. Surette, P.-E. Wolf, J. Stock, S. Leibler. *Protein Mobility in the Cytoplasm of Escherichia coli*. Journal of Bacteriology 181, 197–203 (1999). Primary seven-page PDF opened, abstract and measurement discussion inspected. GFP mobility is a tracer anchor. It does not calibrate the abstract catalyst, capture radius, or damage rate. [doi:10.1128/JB.181.1.197-203.1999](https://doi.org/10.1128/JB.181.1.197-203.1999)

2. D. F. Anderson, T. G. Kurtz. *Continuous Time Markov Chain Models for Chemical Reaction Networks*. Author-hosted survey formulation consulted, including the random time change and martingale problem. This is established mathematical infrastructure. [Author-hosted source](https://people.math.wisc.edu/~dfanderson/papers/SURVEY_AndKurtz.pdf)

3. I. C. Agbanusi, S. A. Isaacson. *A Comparison of Bimolecular Reaction Models for Stochastic Reaction Diffusion Systems*. Bulletin of Mathematical Biology (2013). Full HTML opened. Sections 2–3, the convergence statements in Section 5, and the mean-binding appendix were identified. The spherical target benchmark is already part of this paper. [Full text](https://arxiv.org/html/1301.0547v2), [doi:10.1007/s11538-013-9833-6](https://doi.org/10.1007/s11538-013-9833-6)

4. D. T. Gillespie, L. R. Petzold, E. Seitaridou. *Validity conditions for stochastic chemical kinetics in diffusion-limited systems*. Journal of Chemical Physics 140, 054111 (2014). Full primary article opened, abstract and propensity/voxel validity formulation inspected. [doi:10.1063/1.4863990](https://doi.org/10.1063/1.4863990)

5. A. M. Berezhkovskii, A. Szabo. *Effect of ligand diffusion on occupancy fluctuations of cell-surface receptors*. Journal of Chemical Physics 139, 121910 (2013). Primary article extraction of Section II and the diffusion-renormalized association and dissociation formulas inspected. Direct opening of PMC encountered a browser challenge. Bibliographic metadata checked against PubMed. [doi:10.1063/1.4816105](https://doi.org/10.1063/1.4816105)

6. Y. Zhang, S. A. Isaacson. *Detailed balance for particle models of reversible reactions in bounded domains*. Journal of Chemical Physics 156, 204105 (2022). Full HTML opened. Pointwise detailed balance and boundary-dependent unbinding kernels in Sections III–IV inspected. [Full text](https://arxiv.org/html/2201.04099v3), [doi:10.1063/5.0085296](https://doi.org/10.1063/5.0085296)

7. T. Oppelstrup et al. *First-passage kinetic Monte Carlo method*. Physical Review E 80, 066701 (2009). Primary publisher abstract and arXiv record inspected. Detailed algorithm implementation has yet to be reviewed for the proposed pilot. [doi:10.1103/PhysRevE.80.066701](https://doi.org/10.1103/PhysRevE.80.066701)

8. S. A. Isaacson, A. Mauro, J. Newby. *Uniform asymptotic approximation of diffusion to a small target: generalized reaction models*. Physical Review E 94, 042414 (2016). Abstract and bibliographic record inspected. This source directly anticipates uniform-in-time approximation as a general research direction. A full comparison remains necessary before proposing a new error-bound theorem. [doi:10.1103/PhysRevE.94.042414](https://doi.org/10.1103/PhysRevE.94.042414)

9. X. Li, T. Chou. *A generalized work theorem for stopped stochastic chemical reaction networks*. arXiv:2601.12854 (2026 preprint). Abstract inspected during the preceding discussion. It covers compensated Poisson martingales and bounded arbitrary stopping times. Full theorem comparison remains pending. [Preprint](https://arxiv.org/abs/2601.12854)

10. F. C. Collins, G. E. Kimball. *Diffusion-controlled reaction rates*. Journal of Colloid Science 4, 425–437 (1949). Publisher metadata confirmed. Formulas used here are also explicitly reproduced and discussed in the inspected later primary theoretical papers. [doi:10.1016/0095-8522(49)90023-9](https://doi.org/10.1016/0095-8522(49)90023-9)

## Search receipt

Queries included protein diffusion coefficients in cytoplasm, Smoluchowski capture and reflecting-sphere first-passage times, reversible diffusion-influenced reactions, detailed balance in bounded domains, first-passage kinetic Monte Carlo, and validity of well-mixed stochastic chemical kinetics. Exact source titles and DOIs were followed to primary records. The search selected practical model and verification dependencies. It is not an exhaustive novelty search for every possible coupled repair model.

## Parameter gaps that remain material

The capture radius and molecular orientation constraints are assumed. A relevant damage lifetime and intrinsic surface reactivity have not been measured for the proposed catalyst. A GFP diffusion coefficient cannot substitute for those data. The current chemical model uses dimensionless time and a reference volume. Mapping its rates into physical units requires volume-dependent propensities and a chosen molecular mechanism. These gaps determine whether the spatial extension is experimentally relevant.
