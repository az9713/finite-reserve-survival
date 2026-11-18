"""Update discovery links and the journey from an executed comparison receipt."""
from fractions import Fraction as F
import json
from pathlib import Path
import re
from build_validated_continuum_docs import decimal_bound

ROOT=Path(__file__).resolve().parents[1]


def main():
    data=json.loads((ROOT/'results/continuum-validated-comparison.json').read_text())
    cover=data['parameter_cover'];good=cover['strict_gap_arithmetically_verified']
    gap=decimal_bound(100*F(cover['spectator_minus_best_fixed_lower']),7)
    label='Validated Brownian continuum advantage at T=2' if good else 'Validated Brownian continuum bounds at T=2'
    verdict=(f'The spectator strategy beats **every constant paid intensity in [0,1]** at time two by at least **{gap} percentage points**. The proof includes spatial Fourier tails, exact polynomial time integration, roundoff and a complete parameter-interval cover.' if good else 'The Brownian observable and the whole constant-intensity interval now have rigorous residual enclosures. Their remaining overlap leaves the sign unresolved.')
    latest=f'**Latest continuum result:** [{label}](CONTINUUM-CERTIFICATE.html). {verdict} Short-time and eventual long-time loss remain proved; the comparison is limited to the declared shared-resource model and constant-intensity class.'
    builder=ROOT/'scripts/build_research_docs.py'
    text=builder.read_text(encoding='utf-8')
    text=re.sub(r'^\*\*Latest continuum result:\*\*.*$',lambda _:latest,text,flags=re.M)
    if '<a href="CONTINUUM-CERTIFICATE.html">Continuum certificate</a>' not in text:
        text=text.replace('<a href="PAID-ORGANIZATION.html">Paid organization</a>\')',
            '<a href="PAID-ORGANIZATION.html">Paid organization</a><a href="CONTINUUM-CERTIFICATE.html">Continuum certificate</a>\')')
    builder.write_text(text,encoding='utf-8')

    overview=f'''## {label}\n\n{verdict}\n\nThe comparison uses the original physical parameters, uniform dispersed start and four fuel molecules. A spectator-dependent reaction forms and maintains attraction without an external controller. [Complete proof and evidence](CONTINUUM-CERTIFICATE.html) · [Exact continuum comparison](../results/continuum-validated-comparison.json) · [Independent replay](../results/continuum-validated-independent-replay.json).\n\n![Validated continuum comparison and full-interval error budget](../results/figures/continuum-validated-certificate.png)\n\nThe earlier short-time theorem proves spectator loss for \\(0<t\\le0.0008\\). The passive ground state has exact decay rate \\((7-\\sqrt{{17}})/4\\); the spectator has a strictly larger decay rate and eventually loses. These statements coexist with a finite-horizon advantage when the certified time-two gap is positive. The result does not establish a molecular realization, a chemistry-to-biology transition, or confirmed publication priority.\n\n[Earlier continuum proofs and historical calculation](CONTINUUM-ORGANIZATION.html) · [Residual theorem](../research/continuum-residual.md) · [Validated spatial arithmetic](../research/validated-spatial.md) · [Exact time analysis](../research/validated-time.md) · [Final verification](../results/continuum-validated-verification.json)\n'''
    if (ROOT/'docs/continuum-assessment.md').exists():
        overview+='\n**What do we learn, and what could be novel next?** The result certifies a small continuum advantage within a specific family; it establishes no new general concept or chemistry-to-biology transition. The strongest proposed next step is a predictive mechanism theorem giving conditions for benefit and failure across a defined class of resource-matched reaction diffusions. Its novelty is unestablished. [Full assessment, literature comparison and mathematical next-step specification](CONTINUUM-CERTIFICATE.html#what-do-we-learn-from-the-certified-continuum-result).\n'
    (ROOT/'docs/continuum-overview.md').write_text(overview,encoding='utf-8')
    journey=ROOT/'docs/continuum-journey.md'
    text=journey.read_text(encoding='utf-8').split('\n## 15.')[0]
    text+=f'''\n## 15. Closing the time-two continuum gap\n\nThe user asked for validated residual bounds and complete coverage of constant intensities. Parallel agents developed Arb radial Fourier integrals, exact dyadic trial/time reconstruction, and an independent exact-rational parameter checker. All 40 retained trial archives were replayed: 104,104,000 integer coefficients passed initial-value, continuity, state and mode checks.\n\nThe first approach exposed two engineering bottlenecks: Fraction-per-entry time refinement and very small positive operands in roundoff-error matrix products. Exact integer refinement and monotone upward inflation in error products removed them. Both changes were independently checked; stopped attempts emitted no survival certificate and remain recorded as engineering diagnostics. The physical model and prespecified time were unchanged.\n\n{verdict}\n\nThe final calculation retains the nonzero corrected-observable integral, every unresolved Fourier mode, exact time integration, all node receipts and an independent full-interval Bernstein replay. The parameter checker proves a Chebyshev nodal amplification bound below 3.279 and an analytic interpolation remainder below 1.8e-7. An additional corollary bounds a nonempty advantage interval around time two whenever the certified gap is positive.\n\n![The completed continuum certificate](../results/figures/continuum-validated-certificate.png)\n\nThis closes a specific mathematical validation gap. The result remains a narrow designed example, with publication priority and an elementary chemical implementation unestablished. [Full new report](CONTINUUM-CERTIFICATE.html) · [Execution and verification receipts](../results/continuum-validated-verification.json) · [Engineering record](../results/continuum-validation-engineering.json)\n'''
    if (ROOT/'docs/continuum-assessment.md').exists():
        text+='\n## 16. Interpreting the certificate without inflating novelty\n\nThe next question asked what the continuum result teaches, whether it establishes a new concept, and which next step offers the strongest prospect of an original result. The assessment distinguishes a rigorous model-specific existence example from a general principle. A fresh focused literature comparison identifies precedents for interaction protocols, autonomous feedback, adjoint optimization, spatial detailed balance and residual estimation, recording access limits.\n\nThe proposed next target is a structural mechanism theorem for physically admitted fuel-dependent organization under a common reserve. A standard sensitivity/Taylor argument supplies an entry point; deriving predictive conditions, a useful quantitative bound, a null regime and fair stronger comparators remains unproved. An unspecified open neighborhood follows from continuity and is not treated as a new robustness theorem. [Interpretation and next-theorem section](CONTINUUM-CERTIFICATE.html#what-do-we-learn-from-the-certified-continuum-result) · [Literature check](../research/continuum-assessment-literature.md). No new model result or publication priority is claimed in this documentation step.\n'
    journey.write_text(text,encoding='utf-8')

    old=ROOT/'docs/continuum-organization.md'
    archive=ROOT/'docs/continuum-organization-initial.md'
    htmlarchive=ROOT/'docs/CONTINUUM-ORGANIZATION-INITIAL.html'
    if not archive.exists(): archive.write_bytes(old.read_bytes())
    if not htmlarchive.exists(): htmlarchive.write_bytes((ROOT/'docs/CONTINUUM-ORGANIZATION.html').read_bytes())
    initial=archive.read_text(encoding='utf-8')
    heading,rest=initial.split('\n',1)
    banner=f'\n\n> **Subsequent validation:** [{label}](CONTINUUM-CERTIFICATE.html). {verdict}\n> The derivations and numerical discussion below record the earlier unresolved stage. Its [original HTML snapshot](CONTINUUM-ORGANIZATION-INITIAL.html) and source are retained byte for byte so the historical receipt remains assessable.\n'
    old.write_text(heading+banner+rest,encoding='utf-8')

    org=ROOT/'docs/organization-overview.md'
    text=org.read_text(encoding='utf-8').replace('A rigorous continuum error enclosure remains open.',
        'The subsequent [validated Brownian comparison](CONTINUUM-CERTIFICATE.html) addresses the continuum error and the complete constant-intensity family.')
    org.write_text(text,encoding='utf-8')

    readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8')
    start=text.index('**[Brownian continuum proofs, error bounds and the unresolved T=2 question]') if '**[Brownian continuum proofs, error bounds and the unresolved T=2 question]' in text else text.index('<!-- continuum-certificate:start -->')
    end=text.index('[![Scientific dossier preview]',start)
    replacement=f'''<!-- continuum-certificate:start -->\n**[{label}](https://az9713.github.io/finite-reserve-survival/docs/CONTINUUM-CERTIFICATE.html)**\n\n{verdict}\n\nThe same strategy loses at sufficiently short and sufficiently long horizons. This is a rigorously enclosed finite-horizon example in one designed model; no new general concept, molecular realization or confirmed novelty is claimed. The validated calculation uses 38 fixed-intensity nodes, an exact full-interval cover and retained integer trial archives.\n\n[![Validated Brownian continuum comparison](results/figures/continuum-validated-certificate.png)](https://az9713.github.io/finite-reserve-survival/docs/CONTINUUM-CERTIFICATE.html)\n\n[Exact comparison receipt](results/continuum-validated-comparison.json) · [Independent replay](results/continuum-validated-independent-replay.json) · [Proof and implementation](research/validated-residual-implementation.md) · [Verification](results/continuum-validated-verification.json) · [Earlier short/long-time proofs](docs/CONTINUUM-ORGANIZATION.html)\n<!-- continuum-certificate:end -->\n\n'''
    text=text[:start]+replacement+text[end:]
    if (ROOT/'docs/continuum-assessment.md').exists():
        text=text.replace('<!-- continuum-certificate:end -->',
            '[What we learn, novelty limits, and the strongest next theorem](https://az9713.github.io/finite-reserve-survival/docs/CONTINUUM-CERTIFICATE.html#what-do-we-learn-from-the-certified-continuum-result)\n<!-- continuum-certificate:end -->')
    text=text.replace('The effect is small and has no rigorous Brownian-continuum error enclosure.',
        'The effect is small; the subsequent validated continuum certificate above now supplies a Brownian error enclosure and covers the full constant-intensity interval.')
    text=text.replace('Python 3.10+ with numpy, scipy, and matplotlib. PDF export additionally uses reportlab.',
        'Python 3.10+ with numpy, scipy, and matplotlib. The rigorous continuum calculation additionally requires python-flint >= 0.9 (the validation extra). PDF export additionally uses reportlab.')
    text=text.replace('python -m pip install -e .\n','python -m pip install -e ".[validation]"\n')
    if 'python -m pip install -e ".[validation]"' not in text:
        text+='\n## Reproduce the continuum certificate\n\n```powershell\npython -m pip install -e ".[validation]"\npython scripts/validate_continuum_nodes.py --policy spectator\n0..37 | ForEach-Object { python scripts/validate_continuum_nodes.py --node-index $_ }\npython scripts/complete_validated_continuum.py --cache results/continuum-residual-matrices-rebuild\npython scripts/check_validated_continuum.py --receipt results/continuum-validated-comparison.json --lebesgue --output results/continuum-validated-independent-replay.json\n```\n\nThe trial-generation solve proposes exact dyadic functions. Its errors are included in the residual certificate. The independent checker replays parameter arithmetic and the retained receipt chain; it does not independently recompute every spatial contraction. Matrix caches can be rebuilt from Arb integrals. The arithmetic contract and limits are stated in the proof notes.\n'
    readme.write_text(text,encoding='utf-8')
    print('Integrated the executed continuum result into README, dossier source and journey; preserved historical continuum snapshots.')


if __name__=='__main__': main()
