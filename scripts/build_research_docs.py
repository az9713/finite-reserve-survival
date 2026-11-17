"""Build the public journey, scientific dossier, and manuscript from audited sources."""
from pathlib import Path
import argparse
import hashlib
import html
import json
import re
import subprocess
import sys
from datetime import datetime, timezone

import markdown

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'

STYLE = '''
:root{color-scheme:dark;--bg:#0d141d;--panel:#152233;--ink:#e6edf7;--muted:#adc0d6;--line:#2b4159;--accent:#72dbc1}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--bg);color:var(--ink);font:17px/1.75 'Segoe UI',Arial,sans-serif}a{color:var(--accent);text-underline-offset:3px}nav{padding:18px 5vw;border-bottom:1px solid var(--line);display:flex;gap:24px;flex-wrap:wrap}main{max-width:1000px;margin:auto;padding:45px 30px 100px}h1{font:500 clamp(35px,5vw,65px)/1.12 Georgia,serif;letter-spacing:-.025em;margin:20px 0 30px}h2{font:500 32px/1.25 Georgia,serif;margin:65px 0 20px;padding-top:20px;border-top:1px solid var(--line)}h3{font-size:21px;margin:35px 0 12px}p{margin:16px 0}li{margin:9px 0}img{width:100%;height:auto;border:1px solid var(--line);border-radius:8px;margin:20px 0}pre{padding:20px;background:var(--panel);overflow:auto;border-radius:8px;font-size:14px}code{font-family:Consolas,monospace}table{width:100%;border-collapse:collapse;margin:25px 0;font-size:15px}td,th{padding:12px;text-align:left;vertical-align:top;border-bottom:1px solid var(--line)}th{background:var(--panel)}blockquote{margin:25px 0;padding:15px 25px;border-left:3px solid var(--accent);background:var(--panel)}.meta{color:var(--muted);font-size:14px;letter-spacing:.04em}.toc{background:var(--panel);padding:22px 30px;border-radius:9px}.toc ul{columns:2;column-gap:35px;padding-left:20px}.toc li{break-inside:avoid;font-size:14px}.toc ul ul{display:none}.math.display{overflow-x:auto;padding:12px 0}.katex{font-size:1.08em}.receipt{background:var(--panel);padding:16px 22px;border-left:3px solid var(--accent)}button{font:inherit;background:var(--panel);color:var(--ink);border:1px solid var(--line);border-radius:5px;padding:6px 14px;cursor:pointer}footer{border-top:1px solid var(--line);padding:25px 5vw;color:var(--muted)}@media(max-width:700px){main{padding:25px 18px 60px}body{font-size:16px}h2{font-size:27px}.toc ul{columns:1}table{display:block;overflow-x:auto}nav{gap:13px}.katex{font-size:.95em}}@media print{:root{--bg:white;--ink:#151515;--panel:#f1f3f5;--muted:#444;--line:#bbb;--accent:#165663;color-scheme:light}nav,button,.toc{display:none}main{max-width:none;padding:0}body{font-size:10pt;line-height:1.5}h1{font-size:28pt}h2{font-size:19pt;break-after:avoid}h3{break-after:avoid}img,table,pre{break-inside:avoid}a{color:inherit}footer{font-size:8pt}.katex{font-size:1em}}
'''

JOURNEY = r'''
# From a prediction about science to auditable survival models

The project began with the closing minutes of Jacob Tsimerman's Golden Gate Institute talk. Around 8:22–8:39, the supplied transcript records his prediction that the current landscape of mathematical problems would be resolved within five years, with several generations of progress compressed into that interval. Around 8:57 he asks researchers to reconsider their scientific ambition. At 9:05–9:28 he invokes the chain from mathematics to physics, chemistry, and biology, suggesting that computational difficulty might cease to enforce familiar disciplinary boundaries. At 9:36–9:56 he extends the challenge to philosophy, including attention, love, and value. These are a research invitation and predictions about the future. The experiments below test a small, explicitly defined part of that invitation. [Video from 8:57](https://www.youtube.com/watch?v=kbzIVh405r8&t=537s) · [Transcript provenance and closing excerpt](../research/talk-closing-excerpt.md)

> “we have to start revisiting how ambitious we are with our science” — supplied transcript, 8:57–9:05; audio verification pending.

The user reports that these comments inspired a search for projects that could cross those boundaries. The selected direction visible in the initial package is **Physical possibilities of purposeful activity**, operationalized as a thermal sensor whose errors affect survival. The subsequent research asks how sensing can help a device keep functioning when information, action, and repair all draw on physical resources.

## 1. The starting research package

The available record begins with a research package titled **Physical possibilities of purposeful activity**, together with the talk transcript. The first request was “Read the entire codebase in cwd.” This journey starts from those supplied artifacts without attributing their origin or reconstructing earlier project proposals. It then follows the assessment, the user's successive research questions, the resulting constructions, and the final scope of the evidence.

The local transcript contains rough automatic-transcription wording. This account paraphrases its closing argument and supplies timestamps. The video fetch was throttled during preparation of this dossier, so independent audio verification remains open. The transcript's earlier claims about contemporary mathematics are outside the scientific premises used here.

| Evidence | What it establishes |
|---|---|
| User's account of the motivation | The connection between the talk and the project exploration |
| Initial README, model, tests, figures, and paper | What the first study actually contains |
| This Codex conversation | Requests and implementation history of Studies 02–04 |
| Executable source and current test receipt | Current behavior within the implemented model classes |
| Source manifest and literature notes | Earlier claimed retrievals, with a fresh targeted comparison added here |

## 2. What the initial package delivered

Study 01 contains a four-state thermal sensor, a reduced mismatch chain, exact matrix-exponential survival, a Gillespie simulator, 15 original tests, an offline numerical atlas, and a research note. Two sensors have 80% stationary accuracy and equal mutual information. The faster one uses approximately half the stationary sensing work. Under instantaneous mismatch killing it survives less well, while under four-stage damage with reset it survives better.

| Quantity at horizon 20 | Slow | Fast |
|---|---:|---:|
| Sensor power in thermal energy per time | 0.176666 | 0.087135 |
| Instantaneous-hazard survival | 0.302646 | 0.153486 |
| Four-stage damage survival | 0.382572 | 0.867882 |

![Initial ranking reversal](../results/figures/ranking-reversal.png)

The numerical record includes 30,000 seeded trajectories for each sensor/law combination, 144 transfer cases with 69 ranking reversals, and 1,875 sampled designs at three horizons. These are computational experiments within a declared architecture. The proofs derive a sensing-power frontier and an asymptotic survival ordering. [Original proof](model.md) · [Saved results](../results/summary.json) · [Verification](../results/verification.md)

## 3. Codex's assessment and the change in ambition

The initial inspection reproduced the example and passed all 15 tests. It also found stale source hashes for the atlas template and artifact builder in the earlier numerical run manifest. The release manifest matched the files then present. Those two records answer different questions: file integrity and provenance of an earlier computation. Neither automatically dates a scientific discovery.

Codex's assessment was that the package is a serious tractable example, with established ingredients and unestablished priority. It handles temporal error structure and effective failure well. Free damage reset, externally supplied driving, and missing energy feedback limit the interpretation as a persistent physical agent. The proposed next step was a finite resource reserve with paid repair, followed by a single reserve for sensing, action, and maintenance.

This narrowed the original ambition to statements that could be refuted. What energy is necessary to survive? Which controller attains the bound? Does the advantage remain when the physical architecture changes? The user's subsequent requests drove each extension.

## 4. Paying for repair

The user asked for a minimal thermodynamically consistent model and a sharp survival bound with a matching construction. Study 02 uses an integer reserve. Repair consumes a unit, and the reverse reaction returns that unit while recreating damage. A rate ceiling and a fixed contact time bound the probability of successful repair. Conditioning on the number of stressed checkpoints gives an exact survival frontier. Concavity proves that equal energy allocation is optimal for a known number of repairs.

The resulting theorem bounds every contact schedule in the specified class and gives the attaining schedule. It also shows that reserve size can reverse the preferred error correlation. The construction remains externally phased. [Study 02 proof](finite-reserve-repair.md) · [Executable calculation](../scripts/run_reserve_study.py)

![Paid repair: reserve-dependent ordering and a sharp energy frontier](../results/figures/reserve-frontier.png)

## 5. Closing the operating-energy loop

The next request was to fund sensing, action, and repair from the same reserve, and compare adaptive sensing with the best fixed strategy. Study 03 adds paid measurements, paid engagement, and accessible fuel on correctly chosen engagement. Dynamic programming gives the exact finite-horizon frontier and an attaining controller. Fixed sensing schedules retain optimized actuation and a fallback when a scheduled measurement is unaffordable.

For the symmetric two-unit example, adaptive sensing achieves 7/8 survival in two rounds. The best fixed schedule achieves 3/4. The first probe creates two live reserve states that require opposite decisions in the final round. This makes the benefit of reserve feedback explicit.

![The two-round adaptive construction](../results/figures/closed-loop-witness.png)

Further derivations give an exact accuracy threshold of 2/3 for a strict two-round advantage, a switch in the optimal fixed long-run regime at repair success 2−√2, and distinct exponential survival factors. Removing harvesting eliminates the exponential advantage for a finite reserve. The long-run adaptive construction avoids repair on its surviving paths, which explains why its exponent does not depend on repair success. These distinctions matter when interpreting the theorem. [Study 03 full proofs](closed-loop-survival.md) · [Numerical record](../results/closed-loop.json)

## 6. Removing the clock changed the conclusion

Study 04 replaces rounds with autonomous competing reaction channels. The environment, memory, actuator, damage, and reserve evolve continuously. Forward and reverse rates obey local detailed balance. Survival ends on the first entry into terminal damage. Static gates use the available memory and battery level.

All 324 binary circuits are enumerated. Exact rational enclosures identify a unique time-four winner that senses only at full charge. A second-derivative bound extends the uniform-sensing comparison to every constant intensity between zero and one. The adaptive circuit reaches about 53.38%, while every such constant intensity is bounded above by 50.63%.

![Continuous-time comparisons](../results/figures/autonomous-survival.png)

The long-run binary optimum leaves sensing off. A positive-vector spectral certificate verifies this result. The checkpoint model and the autonomous network change several assumptions at once, including damage dynamics, reverse channels, and control architecture. Their different asymptotic behavior cannot be causally assigned to clock removal alone. [Study 04 proof](autonomous-survival.md) · [Exact certificate](../results/autonomous.json)

## 7. Tools, failures, and repairs

| Tool or method | Role in this project |
|---|---|
| Python, NumPy, SciPy | Markov generators, matrix exponentials, optimization, probability propagation |
| Fraction and integer arithmetic | Outward enclosures and spectral certificates |
| Matplotlib | Scientific figures generated from reproducible numerical calculations |
| Web research | Closest prior work and precise limits on novelty claims |
| PowerShell and file inspection | Workspace inventory, source checks, test execution |
| KaTeX and HTML | Offline equations and browsable evidence in this dossier |

The subsequent implementation was performed inline in Codex. No research subagents or independent expert reviewers were used. The initial package's internal production workflow is outside this record. Monetary costs and token usage were not recorded.

Several recoverable failures are visible in the conversation. A test attempted to multiply a Python list by a floating-point bound and raised `TypeError: can't multiply sequence by non-int of type 'float'`. It was corrected to a NumPy array. An exploratory rate search raised `ValueError: too many values to unpack (expected 7)` and was corrected before producing the retained benchmark. A documentation patch failed because expected lines did not match the current file; the actual file was read and the patch reapplied. A feedback arrow crossed text in a figure; visual inspection led to a smaller arc. The search also produced negative findings, including circuits whose sensing expense degraded survival and the loss of the long-run sensing advantage. Those outcomes shaped the final scientific claims.

An unrelated lecture-note request was explicitly withdrawn by the user as a wrong-project message. A prior request to commit and push was interrupted. The current workspace has no Git repository metadata, so this record does not claim a completed publication or push.

## 8. What this journey demonstrates

The project makes a mathematical chain assessable: transition rates determine error histories; error histories interact with resource-limited repair; observation and action change reserve distributions; those distributions determine the probability of reaching failure. It tests a narrow route through Tsimerman's proposed disciplinary boundaries.

Studies 01–04 begin with supplied architecture. Study 05 adds abstract chemical self-assembly and cooperative repair. Molecular synthesis, reproduction, evolving architecture, and biological measurements remain open. Terms such as purpose and survival here refer to designer-specified operation and stopping rules. The philosophy portion of the talk supplies motivation, while the mathematical claims concern explicitly defined stochastic systems.

The companion scientific document supplies full derivations, an updated literature comparison, and a source-linked receipt for every study. The accompanying manuscript draft concentrates on the mathematically strongest candidate contributions. Priority and independent proof review remain open. [Scientific dossier](SCIENTIFIC-DOSSIER.html) · [Manuscript PDF](../results/research-manuscript.pdf) · [Portable LaTeX sources and figures](../results/research-manuscript-source.zip)

## 9. Making the work assessable

The initial documentation pass reran its 30 tests and the full autonomous benchmark, and generated a new paid-repair figure directly from the model. The scientific dossier keeps the original mathematical assumptions visible and introduces the chemical vocabulary needed to interpret them. The literature comparison records which papers were inspected fully and which only at abstract level. The manuscript is a candidate paper draft with human authorship and independent review still pending.

The rendering checks exercise desktop and phone widths, every local link, every displayed equation, and every figure. A Windows newline-conversion problem initially stopped the LaTeX build; normalization fixed it. An absolutely positioned accessibility layer in the equation renderer initially widened the mobile page; anchoring that layer to its equation fixed the overflow while retaining accessible math. These are document-production corrections. They do not change the mathematical model. [Scientific execution receipt](../results/dossier-execution-receipt.json) · [Browser and source-integrity checks](../results/dossier-ui-validation.json) · [PDF build receipt](../results/manuscript-build.json)
'''

INTRO = r'''
# Information, energy, and the survival of designed systems

This dossier makes five small mathematical models and their spatial maintenance extension inspectable by researchers unfamiliar with the project. It assumes graduate mathematics and physics, and introduces the chemical and biological vocabulary explicitly. The motivating question comes from the closing minutes of Tsimerman's Golden Gate Institute talk: could stronger mathematical reasoning make explanations travel across physics, chemistry, biology, and eventually questions of value? The user selected the physical cost of persistent purposeful behavior as a tractable starting point. [Talk at 9:05](https://www.youtube.com/watch?v=kbzIVh405r8&t=545s) · [Development journey](DEVELOPMENT-JOURNEY.html) · [Spatial repair extension](SPATIAL-REPAIR.html)

## Abstract

We investigate how temporal error structure, finite reserves, and control architecture determine survival in finite-state stochastic models. A thermal binary sensor yields an exact architecture-specific accuracy–power frontier and a reversal of survival ranking across effective damage laws. A reversible paid-repair model yields an attainable reserve frontier. A shared-reserve checkpoint device admits exact optimal survival, an adaptive sensing advantage, an accuracy threshold, and an explicit long-run rate gap from all fixed sensing schedules. A continuous-time autonomous network admits a certified finite-time advantage over constant sensing intensity, while its best binary asymptotic circuit leaves sensing off. A closed six-species network then forms active repair catalysts from monomers and finite fuel. It yields a sharp resource bound, an analytic single-complex optimum, and a certified interior optimum for cooperative assembly and maintenance. The results quantify architecture-dependent conditions for persistence. Molecular realization and biological validation remain open.

## How to read and assess this dossier

Each study below contains its assumptions, mathematical derivation, construction or counterexample, figure, executable source, and numerical receipt. The proofs from the research notes are included in full. The final sections identify which ingredients are established in prior work, which exact results remain candidates for novelty, and what would be required before publication. Figures compare probabilities within the same study. Different studies use different dynamics and incompatible time units, so their percentages should not be compared as a common benchmark.

## Chemical and biological vocabulary

**Chemical potential** is the free-energy change per additional particle of a species. A reservoir maintaining a chemical-potential difference can drive a reaction cycle away from equilibrium. A **chemostat** maintains concentrations externally. Study 04 uses a maintained fuel reservoir, with an internal battery limiting storage. Study 05 includes every fuel and waste molecule inside a closed chemical compartment, coupled to a heat bath.

**Local detailed balance** specifies the log ratio of forward and reverse rates from the entropy transferred to the environment. It constrains ratios, leaving kinetic barriers and absolute rates to be specified. Multiplying both directions of a reaction by the same catalyst or gate factor preserves that ratio. Our reaction networks stipulate those gates. A particular molecule, catalytic pathway, or mass-action construction implementing them has yet to be supplied.

**Repair** means a transition from a designated damaged state to a designated operational state. **Survival** is the probability of avoiding a specified failure set up to a time or round count. **Homeostasis** in biology refers to regulation of internal conditions. **Self-assembly** would require components to form the functional organization through their own dynamics. Studies 01–04 supply the functional organization at initialization or fabrication. Study 05 starts with monomers and forms its active repair catalysts through the declared reactions. This distinction is essential when interpreting a link to biology.

## What crosses a boundary?

| Proposed crossing | Demonstrated achievement | Remaining gap |
|---|---|---|
| Mathematics to physics | Exact formulas and certified finite-state thermodynamic calculations | General many-body physical problems remain untouched |
| Physics to information | Rates and energy gaps determine accuracy, temporal errors, and sensing work | The physical architecture is restricted and largely drawn from established frameworks |
| Information to continued function | Proofs connect observations, reserves, actions, and first-passage survival | Failure and utility are externally specified |
| Physics to chemistry | Reversible reaction-like channels satisfy an explicit energy ledger | Study 05 specifies species and stoichiometry. Kinetic realization and experiments remain open |
| Chemistry to biology | A reaction network forms and maintains its own repair catalysts in Study 05 | Study 05 adds finite-inventory assembly and cooperative repair. Replication, membranes, and evolution remain open |
| Science to philosophy | Operational value of information can be quantified for a chosen task | No theory of love, subjective value, or consciousness follows |

The phrase “overcome a Tsimerman barrier” should be read as solving a declared toy crossing. The talk itself supplies no formal list of barriers, and these examples do not validate its broader forecast about the future of mathematics.
'''

LITERATURE = r'''
## Literature comparison and novelty assessment

The targeted search on 17 November 2026 used combinations of “survival”, “adaptive sensing”, “finite reserve”, “paid repair”, “golden ratio”, “spectral radius”, and autonomous information thermodynamics. Exact-expression searches for the two-thirds sensing threshold and the characteristic rate constants returned substantial unrelated material. A missing search hit supplies limited evidence. The search is not a proof that the formulas have never appeared in another notation or field.

| Prior source | Established overlap | Specific distinction requiring further priority review |
|---|---|---|
| Hartich, Barato and Seifert (2014), [bipartite stochastic thermodynamics](https://arxiv.org/abs/1402.0419) | Four-state examples, information flow, and thermodynamic inequalities | Our equal-accuracy task-ranking combination and restricted kinetic frontier |
| Kolchinsky and Wolpert (2018), [semantic information and agency](https://arxiv.org/abs/1806.08053) | Causal relation between information and viability | The reserve-constrained constructive bounds and specific feedback comparisons |
| Abad, Ercetin and Gündüz (2017), [costly channel sensing](https://arxiv.org/html/1703.10519v2) | Sensing spends energy, battery-dependent policies exploit correlated conditions | Our survival objective, paid repair, exact two-round separation, and explicit exponents |
| Vasconcelos et al. (2019), [energy-harvesting estimation](https://arxiv.org/abs/1908.06070) | Globally optimal scheduling with finite batteries | The present damage process and survival criterion |
| Guruacharya and Hossain (2018 version), [self-sustainability and outage](https://arxiv.org/abs/1702.01648) | Eventual outage, stochastic energy input, asymptotic bounds | Finite capacity with actions gating fuel and explicit repair |
| Tkachev et al. (2014), [controlled Markov reachability](https://arxiv.org/html/1407.5449v1) | Optimal event probabilities and dynamic programming | Our explicit closed-form example and model-specific constructions |
| Blondel and Nesterov (2009), [nonnegative matrix products](https://perso.uclouvain.be/vincent.blondel/publications/09BN.pdf) | Spectral growth bounds and special matrix families | The two matrices induced by this sensing task and its regime-switch formula |
| Ehrich and Sivak (2023), [autonomous energy and information flows](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2023.1108357/full) | Autonomous information thermodynamics and molecular-machine interpretation | Our shared-reserve first-passage comparison and its digital-circuit certificate |
| Fox and Glynn (1988), [Poisson probability computation](https://web.stanford.edu/~glynn/papers/1988/FG88.html) | Rigorous Poisson truncation and CTMC evaluation | Our small integer-grid implementation is a verification tool within established methodology |
| Chen, Chen and Ma (2026), [finite-time autonomous information machine](https://arxiv.org/abs/2604.15953) | Finite-time autonomous information/energy tradeoffs | The abstract studies erasure and refrigeration; our survival circuit differs. Full theorem comparison remains pending |
| Liu and Darabi (2005), [control with a dynamic sensing set](https://skoge.folk.ntnu.no/prost/proceedings/ifac2005/Fullpapers/03504.pdf) | Controller survival, changing observability, and sensor repair | Their finite-automaton formulation does not give this shared-reserve thermodynamic comparison |
| Hasidi and Kalech (2021), [anticipatory troubleshooting](https://www.mdpi.com/2076-3417/11/3/995) | Survival-aware diagnosis and maintenance decisions | Their reported experiments omit sensing costs; the objectives and physical constraints differ |

The Hartich, Kolchinsky, Vasconcelos, Guruacharya, and Chen entries were checked at abstract or overview level in this pass. The costly-sensing and reachability full HTML pages and the matrix-product paper were opened for their formulations. The Liu–Darabi introduction and the Hasidi–Kalech experimental discussion were also inspected. This is a targeted comparison, not an exhaustive full-text survey. The initial package contains a separate five-paper thermodynamic audit. [Earlier literature audit](../research/literature.md) · [Dated search receipt and access depths](../research/priority-search-2026-11-17.md)

**Assessment.** Generic information-to-viability connections, autonomous information engines, battery-dependent sensing, Bellman optimization, and uniformization are established. No exact prior match was identified for the present complete combination of assumptions and formulas. The most focused novelty candidates are the checkpoint rate formula and switch at 2−√2, the model-specific noise threshold, and the certified continuous-time separation with its asymptotic reversal. Their priority remains unestablished.

## Publication readiness and limitations

The accompanying manuscript is a candidate research draft with full mathematical definitions, proofs, references, and computational supplements. It carries no claim of established priority, peer review, or arXiv acceptance. Its author line identifies an AI-assisted workspace draft so that a human author can take responsibility before submission.

The sharp checkpoint bounds are within stated policy classes. The continuous-time sharp optimum is over binary circuit gates with restricted control tables; its uniform-intensity comparator permits a continuum only in the sensing multiplier. Arbitrary analogue reserve-dependent sensing, arbitrary reaction architectures, and every possible controller are outside that optimization. The finite-time autonomous benchmark was selected through exploratory search. The certificate verifies the exact chosen rates; it supplies no frequency estimate for how often such advantages occur in physical devices.

The direct-sensor-disconnection comparison changes its channel structure and may retain indirect information through the reserve. It demonstrates the performance difference for that intervention. It does not isolate a universal number of useful bits or a general thermodynamic price of meaning. Operating energy omits apparatus preparation and fabrication. The checkpoint model relies on supplied clocks and registers, while the autonomous model relies on a maintained fuel reservoir and designed barriers.

Before submission, a human reviewer should check the optimal-policy inequalities, the matrix-product proof for all fixed words, and the implementation of interval rounding. Priority review should extend to partially observed control, intermittent toxic exposure, occupation-time reliability, joint spectral radii, and autonomous molecular networks. An explicit mass-action realization would strengthen the chemical interpretation. Study 05 adds a self-assembly model with finite inventories. A molecular implementation remains open.
'''

REFERENCES = r'''
## References

1. U. Seifert. *Stochastic thermodynamics, fluctuation theorems and molecular machines*. Reports on Progress in Physics 75, 126001 (2012). [arXiv:1205.4176](https://arxiv.org/abs/1205.4176).
2. D. Hartich, A. C. Barato, U. Seifert. *Stochastic thermodynamics of bipartite systems: transfer entropy inequalities and a Maxwell's demon interpretation*. Journal of Statistical Mechanics, P02016 (2014). [arXiv:1402.0419](https://arxiv.org/abs/1402.0419).
3. J. M. Horowitz, M. Esposito. *Thermodynamics with continuous information flow*. Physical Review X 4, 031015 (2014). [arXiv:1402.3276](https://arxiv.org/abs/1402.3276).
4. S. Still, D. A. Sivak, A. J. Bell, G. E. Crooks. *Thermodynamics of prediction*. Physical Review Letters 109, 120604 (2012). [arXiv:1203.3271](https://arxiv.org/abs/1203.3271).
5. A. Kolchinsky, D. H. Wolpert. *Semantic information, autonomous agency, and nonequilibrium statistical physics*. Interface Focus 8, 20180041 (2018). [arXiv:1806.08053](https://arxiv.org/abs/1806.08053).
6. M. Salehi Heydar Abad, O. Ercetin, D. Gündüz. *Channel Sensing and Communication over a Time-Correlated Channel with an Energy Harvesting Transmitter* (2017). [arXiv:1703.10519](https://arxiv.org/abs/1703.10519).
7. M. M. Vasconcelos, M. Gagrani, A. Nayyar, U. Mitra. *Optimal scheduling strategy for networked estimation with energy harvesting* (2019). [arXiv:1908.06070](https://arxiv.org/abs/1908.06070).
8. S. Guruacharya, E. Hossain. *Self-Sustainability of Energy Harvesting Systems: Concept, Analysis, and Design* (2017, revised 2018). [arXiv:1702.01648](https://arxiv.org/abs/1702.01648).
9. I. Tkachev, A. Mereacre, J.-P. Katoen, A. Abate. *Quantitative model-checking of controlled discrete-time Markov processes* (2014). [arXiv:1407.5449](https://arxiv.org/abs/1407.5449).
10. V. D. Blondel, Y. Nesterov. *Polynomial-time computation of the joint spectral radius for some sets of nonnegative matrices*. SIAM Journal on Matrix Analysis and Applications 31, 865–876 (2009). [doi:10.1137/080723764](https://doi.org/10.1137/080723764).
11. J. Ehrich, D. A. Sivak. *Energy and information flows in autonomous systems*. Frontiers in Physics 11, 1108357 (2023). [doi:10.3389/fphy.2023.1108357](https://doi.org/10.3389/fphy.2023.1108357).
12. B. L. Fox, P. W. Glynn. *Computing Poisson Probabilities*. Communications of the ACM 31, 440–445 (1988). [Author page](https://web.stanford.edu/~glynn/papers/1988/FG88.html).
13. W. Chen, M. Chen, Y.-H. Ma. *Finite-Time Thermodynamics of an Autonomous Information Machine* (2026). [arXiv:2604.15953](https://arxiv.org/abs/2604.15953).
14. J. Liu, H. Darabi. *Reconfiguration of Discrete Event System Controllers with Dynamic Sensing Set*. IFAC World Congress (2005). [Proceedings PDF](https://skoge.folk.ntnu.no/prost/proceedings/ifac2005/Fullpapers/03504.pdf).
15. N. Hasidi, M. Kalech. *Anticipatory Troubleshooting*. Applied Sciences 11(3), 995 (2021). [doi:10.3390/app11030995](https://doi.org/10.3390/app11030995).
'''

def render(source):
    equations = []
    def hold(match):
        value = match.group(0)
        display = value.startswith('\\[')
        tex = value[2:-2]
        equations.append('<span class="math' + (' display' if display else '') + '" data-display="' + str(display).lower() + '">' + html.escape(tex) + '</span>')
        return f'MATHPLACEHOLDER{len(equations)-1}END'
    source = re.sub(r'\\\[.*?\\\]|\\\(.*?\\\)', hold, source, flags=re.S)
    engine = markdown.Markdown(extensions=['tables', 'fenced_code', 'toc'], extension_configs={'toc': {'toc_depth': '2-2'}})
    body = engine.convert(source)
    for i, equation in enumerate(equations):
        body = body.replace(f'MATHPLACEHOLDER{i}END', equation)
    return engine.toc, body

def page(filename, source, title):
    toc, body = render(source)
    first_section = body.find('<h2')
    if first_section >= 0:
        body = body[:first_section] + '<details class="contents"><summary>Contents and navigation</summary>' + toc + '</details>' + body[first_section:]
    js = '''document.querySelectorAll('.math').forEach(el=>{try{katex.render(el.textContent,el,{displayMode:el.dataset.display==='true',throwOnError:true,strict:'ignore'})}catch(e){el.dataset.error=e.message;console.error(e.message)}});document.getElementById('print').onclick=()=>window.print();'''
    extra_css = '.katex{position:relative}.math.display{display:block;max-width:100%}.math:not(.display){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}.math{scrollbar-width:none}.math::-webkit-scrollbar{display:none}.katex-display{margin:0;text-align:left}td{overflow-wrap:anywhere}.contents{margin:30px 0}.contents summary{cursor:pointer;color:var(--accent)}'
    content = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>{STYLE}{extra_css}</style><link rel="stylesheet" href="assets/katex/katex.min.css"></head><body><nav aria-label="Documents"><a href="DEVELOPMENT-JOURNEY.html">Development journey</a><a href="SCIENTIFIC-DOSSIER.html">Scientific dossier</a><a href="SELF-ASSEMBLY.html">Self-assembly study</a><a href="RESEARCH-MANUSCRIPT.html">Manuscript</a><button id="print">Print or save PDF</button></nav><main><p class="meta">Research dossier · 17 November 2026 · Candidate results, priority unestablished</p>{body}</main><footer>Equations and figures are local. Sources link to the research package and primary literature.</footer><script src="assets/katex/katex.min.js"></script><script>{js}</script></body></html>'''
    content = content.replace('<a href="SELF-ASSEMBLY.html">Self-assembly study</a>',
                              '<a href="SELF-ASSEMBLY.html">Self-assembly study</a><a href="SPATIAL-REPAIR.html">Spatial repair</a>')
    content = content.replace('</head>', '<link rel="icon" type="image/svg+xml" href="assets/favicon.svg"></head>')
    (DOCS / filename).write_text(content, encoding='utf-8')

def study(name, label):
    text = (DOCS / name).read_text(encoding='utf-8')
    text = re.sub(r'^# .+$', '## ' + label, text, count=1, flags=re.M)
    # Existing subsection headings become h3 to keep document navigation compact.
    lines = text.splitlines()
    lines = [line if i == 0 else ('#' + line if line.startswith('##') else line) for i, line in enumerate(lines)]
    lines.insert(1, '\n*Model boundary note: the derivation below retains the assumptions of this study. Statements about later extensions are historical and should be read alongside the subsequent studies.*\n')
    if name == 'model.md':
        lines.insert(2, '\nReference labels L1, L2, and L5 in this derivation refer to the [original literature audit](../research/literature.md).\n')
        lines.insert(3, '\n![Equal accuracy does not determine the survival ranking](../results/figures/ranking-reversal.png)\n')
    if name == 'finite-reserve-repair.md':
        lines.insert(2, '\n![Exact reserve frontier and equal-gap energy bound](../results/figures/reserve-frontier.png)\n\n[Figure data](../results/reserve.json) · [Reproduction script](../scripts/build_reserve_figure.py)\n')
    return '\n'.join(lines)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--site-only', action='store_true',
                        help='Refresh the dossier and journey while retaining the Studies 01–05 manuscript snapshot.')
    args = parser.parse_args()
    run = subprocess.run([sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests', '-v'], cwd=ROOT,
                         env={**__import__('os').environ, 'PYTHONPATH': str(ROOT / 'src')}, capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    receipt = ROOT / 'results' / 'dossier-verification.json'
    sources = sorted([*ROOT.glob('src/purposeful/*.py'), *ROOT.glob('tests/test_*.py'), *ROOT.glob('docs/*.md'),
                      *ROOT.glob('scripts/run_*.py'), *ROOT.glob('scripts/build_*.py'), *ROOT.glob('scripts/*assembly*.py'), *ROOT.glob('research/assembly-*.md'),
                      *[ROOT / 'results' / name for name in ['summary.json', 'closed-loop.json', 'autonomous.json', 'reserve.json', 'assembly.json', 'assembly-evidence.json', 'assembly-certificate-check.json', 'spatial-repair.json', 'spatial-repair-verification.json', 'spatial-repair-review-check.json']]])
    sources = sorted(set(p for p in sources if p.name != 'research-manuscript.md'))
    record = {'generated_utc': datetime.now(timezone.utc).isoformat(), 'command': 'python -B -m unittest discover -s tests -v',
              'test_output': run.stdout + run.stderr, 'return_code': run.returncode,
              'sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    receipt.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    receipt_md = '''\n## Evidence and reproduction receipts\n\nThe dossier build reran the existing test suite. The machine-readable receipt records the full output, UTC execution time, and SHA-256 hashes of source, tests, and mathematical notes. Passing tests supports implementation agreement within the tested cases. The proofs and certificates carry the mathematical claims. [Fresh verification receipt](../results/dossier-verification.json)\n\n| Study | Executable source | Reproduction | Numerical evidence |\n|---|---|---|---|\n| 01 | [model.py](../src/purposeful/model.py) | [run_study.py](../scripts/run_study.py) | [summary](../results/summary.json), [all 144 cases](../results/transfer.csv) |\n| 02 | [reserve.py](../src/purposeful/reserve.py) | [run_reserve_study.py](../scripts/run_reserve_study.py) | [independent full-state test](../tests/test_reserve.py) |\n| 03 | [closed_loop.py](../src/purposeful/closed_loop.py) | [run_closed_loop.py](../scripts/run_closed_loop.py) | [arrays and example](../results/closed-loop.json) |\n| 04 | [autonomous.py](../src/purposeful/autonomous.py) | [run_autonomous.py](../scripts/run_autonomous.py) | [integer and rational enclosures](../results/autonomous.json) |\n| 05 | [assembly.py](../src/purposeful/assembly.py) | [run_assembly.py](../scripts/run_assembly.py) | [certificate](../results/assembly.json), [full replay](../results/assembly-certificate-check.json), [simulation](../results/assembly-evidence.json) |\n\nThe earlier run-manifest contains historical hashes and package versions for Study 01. Two hashes were already stale on first inspection. The present receipt identifies the exact sources used for this dossier's fresh test run. The 38 tests are a local reproducibility check. Independent scientific review remains pending.\n'''
    guide = (DOCS / 'math-results-guide.md').read_text(encoding='utf-8')
    test_count = re.search(r'Ran (\d+) tests?', run.stdout + run.stderr).group(1)
    receipt_md = receipt_md.replace('The 38 tests', 'The ' + test_count + ' tests')
    spatial_overview = (DOCS / 'spatial-repair-overview.md').read_text(encoding='utf-8')
    spatial_journey = (DOCS / 'spatial-repair-journey.md').read_text(encoding='utf-8')
    repair = (DOCS / 'repair-lemma.md').read_text(encoding='utf-8')
    assembly = study('self-assembly.md', 'Study 05 — Self-assembly and cooperative maintenance')
    references = REFERENCES + '\n' + (ROOT / 'research/assembly-references.md').read_text(encoding='utf-8')
    assembly_results = (DOCS / 'assembly-results.md').read_text(encoding='utf-8')
    assembly += '\n\n' + assembly_results
    science = '\n\n'.join([INTRO, spatial_overview, guide, study('model.md', 'Study 01 — Thermal sensing, temporal errors, and effective survival'), study('finite-reserve-repair.md', 'Study 02 — Finite reserve and paid repair'), study('closed-loop-survival.md', 'Study 03 — A common reserve for sensing, actuation, and repair'), study('autonomous-survival.md', 'Study 04 — Autonomous reactions and continuous-time failure'), assembly, receipt_md, LITERATURE, references])
    manuscript = r'''# Finite-resource survival from sensing to self-assembled repair

AI-assisted research workspace draft · 17 November 2026

## Abstract

We study survival of devices whose sensing, control, and repair compete for finite stored energy. In a checkpoint model with endogenous harvesting, we construct a strict two-round sensing-feedback advantage and derive its accuracy threshold. For a symmetric two-unit device, we prove an exact asymptotic separation from every fixed sensing word and identify a switch in its optimal fixed regime. In an autonomous reversible reaction network, we certify a finite-time advantage over all constant sensing intensities within a specified control class, while the best binary asymptotic circuit leaves sensing off. A closed six-species reaction network extends the analysis to self-assembled repair catalysts. We derive a sharp net-resource bound with an attaining reaction path and an analytic single-complex optimum, then certify an interior two-complex design optimum over a full parameter square. The comparisons establish architecture-dependent survival tradeoffs. We provide explicit assumptions, attaining constructions, numerical certificates, and reproducible code. Priority and independent proof review remain open.

## Motivation and relation to previous work

Finite-inventory assembly and chemical persistence have extensive precedents [16–21]. Pathwise sensitivity theory [22] supplies an additional tool for our certified design comparison. Study 05 defines the complete network and proves the applicable resource, thermodynamic, and first-passage statements.

Information thermodynamics and information-dependent viability provide the conceptual setting [1–5,11]. Battery-dependent costly sensing and energy-harvesting control are established [6–8]. Controlled reachability and products of nonnegative matrices provide the mathematical tools [9,10]. Our candidate contributions are explicit survival and self-assembly examples with analytically or computationally certified optima. Each statement below quantifies over its declared policy or circuit class. The autonomous system uses a maintained fuel reservoir, and the checkpoint construction uses supplied phase timing and work registers.

The finite-time comparison and the asymptotic comparison are distinct objectives. Changing architectures changes several physical assumptions simultaneously. We therefore claim conditional results within each model, with no causal inference that removing a clock alone determines the outcome.
'''
    manuscript = '\n\n'.join([manuscript, guide, study('closed-loop-survival.md', 'Checkpoint architecture and proofs'), study('autonomous-survival.md', 'Autonomous architecture and certificates'), repair, assembly, receipt_md, LITERATURE, references])
    journey_extension = (DOCS / 'assembly-journey.md').read_text(encoding='utf-8')
    page('SELF-ASSEMBLY.html', (DOCS / 'self-assembly.md').read_text(encoding='utf-8') + '\n\n' + assembly_results + references, 'Self-assembly and maintenance with finite fuel')
    page('DEVELOPMENT-JOURNEY.html', JOURNEY + journey_extension + spatial_journey + receipt_md, 'Development journey — Tsimerman and physical survival')
    page('SCIENTIFIC-DOSSIER.html', science, 'Scientific dossier — information, energy, and survival')
    if not args.site_only:
        page('RESEARCH-MANUSCRIPT.html', manuscript, 'Research manuscript — finite-resource survival and assembly')
        (DOCS / 'research-manuscript.md').write_text(manuscript, encoding='utf-8')
    print('Built journey, scientific dossier, self-assembly page, and fresh verification receipt.'
          + (' Retained the Studies 01–05 manuscript snapshot.' if args.site_only
             else ' Also rebuilt manuscript HTML/source.'))

if __name__ == '__main__':
    main()
