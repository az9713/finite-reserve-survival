"""Render paid organization proofs, independent evidence and publication snippets."""
from datetime import datetime, timezone
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
from urllib.parse import unquote,urlsplit
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np
from playwright.sync_api import sync_playwright
from build_research_docs import page
from verify_research_docs import Links

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def bound(value,digits=8,upper=False,scale=1):
    """Print literal enclosure endpoints outward rather than rounding inward."""
    with localcontext() as context:
        context.prec=80
        scaled=Decimal.from_float(float(value))*Decimal(scale)
        return format(scaled.quantize(Decimal(1).scaleb(-digits),
            rounding=ROUND_CEILING if upper else ROUND_FLOOR),'f')


def current_sources(record):
    for name,expected in record['source_sha256'].items():
        assert sha(ROOT/name)==expected, f'Stale evidence input: {name}'


def main():
    data_path=ROOT/'results/organization.json'; proof_path=ROOT/'results/organization-verification.json'
    data=json.loads(data_path.read_bytes());proof=json.loads(proof_path.read_bytes())
    continuum=json.loads((ROOT/'results/organization-continuum-check.json').read_bytes())
    ablation=json.loads((ROOT/'results/organization-force-ablation.json').read_bytes())
    controls=json.loads((ROOT/'results/organization-certificate-controls.json').read_bytes())
    review=json.loads((ROOT/'results/organization-review.json').read_bytes())
    assert controls['status']=='PASS'
    for record in [data,proof,continuum,ablation,controls,review]:current_sources(record)
    assert data['verification_sha256']==sha(proof_path),'Experiment must be replayed after certificate changes.'
    assert all(c['passed_dkw_band'] for c in proof['paths']['checks'])
    fine=data['refinements'][-1];certificate=proof['certificates'][-1]
    assert certificate['spectator_minus_best_fixed_interval'][0]>0
    run=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s','tests','-v'],cwd=ROOT,
        env={**os.environ,'PYTHONPATH':str(ROOT/'src')},capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
    count=int(re.search(r'Ran (\d+) tests?',run.stdout+run.stderr)[1])
    source=(ROOT/'docs/paid-organization.md').read_text(encoding='utf-8')
    params='| Parameter | Value |\n|---|---:|\n'
    for key,value in data['model'].items():params+=f'| {key} | {value:g} |\n'
    certs='| Lattice n | Certified best constant survival interval | Spectator survival interval | Certified spectator minus best constant (percentage points) |\n|---|---:|---:|---:|\n'
    for row in proof['certificates']:
        lo,hi=row['spectator_minus_best_fixed_interval'];auto=row['spectator_interval'];fixed=row['fixed']
        certs+=f'| {row["n"]} | [{bound(fixed["lower"],11)}, {bound(fixed["upper"],11,True)}] | [{bound(auto[0],11)}, {bound(auto[1],11,True)}] | [{bound(lo,scale=100)}, {bound(hi,upper=True,scale=100)}] |\n'
    certs+='\nThese are absolute survival probabilities in the middle columns. The coarse n=4 grid reverses the result, warning against trusting a single coarse discretization. The certificate covers every real constant intensity in the specified interval, not just the plotted samples.\n'
    refinements='| n | Spectator survival (%) | Fixed attaining-sample survival (%) | Passive-only survival (%) | Spectator minus fixed sample (percentage points) |\n|---|---:|---:|---:|---:|\n'
    for row in data['refinements']:
        values={k:v[-1] for k,v in row['survival'].items()}
        refinements+=f'| {row["n"]} | {100*values["spectator"]:.8f} | {100*values["fixed"]:.8f} | {100*values["passive"]:.8f} | {100*(values["spectator"]-values["fixed"]):.8f} |\n'
    org=f'The following moments use n={data["observables_grid"]} and stop at first collapse or time two. They are unconditional per initial launch, except the explicitly conditional contact mean.\n\n'
    org+='| Policy | Gross paid closures | Gross paid openings | Stopped closed-mode probability | Conditional mean contact activity at time two | Expected stopped fuel |\n|---|---:|---:|---:|---:|---:|\n'
    for name,values in data['observables'].items():
        org+=f'| {name} | {values["paid_forward"][-1]:.6f} | {values["paid_reverse"][-1]:.6f} | {values["stopped_tether"][-1]:.6f} | {values["gate"][-1]/values["mass"][-1]:.6f} | {values["stopped_fuel"][-1]:.6f} |\n'
    maximum=max(max(v.values()) for v in data['stopped_ledger_checks'].values())
    org+=f'\nThe largest stopped expectation-ledger residual is {maximum:.3g}. Initial mode occupancy is zero. Uniform dispersion still permits natural chance contact; it does not require every launch to begin outside the contact ball.\n'
    paths=proof['paths'];n_paths=sum(c['paths'] for c in paths['checks'])
    verify=(f'The current full suite has **{count} passing tests**. Independent finite-lattice jump simulation uses **{n_paths:,} paths** on n={paths["n"]}, across three policies. Its 99% familywise DKW half-width is **{paths["dkw_probability_half_width"]:.6f}**. Channel counts, seeds and energy ledgers are retained in the [finite-grid receipt](../results/organization-verification.json).\n\n')
    verify+='| Lattice path policy | Maximum discrepancy from production solver | Maximum path energy-ledger residual |\n|---|---:|---:|\n'
    for row in paths['checks']:
        verify+=f'| {row["policy"]} | {row["maximum_grid_difference"]:.6f} | {row["max_energy_ledger_error"]:.3g} |\n'
    verify+='\n[Independent scientific review receipt](../results/organization-review.json) records the mathematical and reporting issues found, their resolutions, source fingerprints, and replay checks.\n'
    continuous=('A separately written continuum simulator uses free Brownian Gaussian increments and Euler–Maruyama integration for the attractive-mode force between Poisson reaction candidates. Three step sizes and three policies use **180,000 trajectories** in total. It checks survival, stopped mode/fuel and spatial contact. Its exact integer event ledgers pass at every recorded time. The [continuum receipt](../results/organization-continuum-check.json) retains all nine runs, seeds and source hashes. Sampling bands concern the numerical SDE laws; the Euler bias is unenclosed.\n\n')
    continuous+='```json\n'+json.dumps(continuum['confidence'],indent=2)+'\n```\n'
    al,ah=ablation['gap_interval']
    ablated=f'Without the attractive well, the n=16 certified spectator-minus-best-constant gap is **[{bound(al,scale=100)}, {bound(ah,upper=True,scale=100)}] percentage points**. The attaining fixed sample is alpha={ablation["fixed"]["attaining_sample_alpha"]:g}. In this ablation, the candidate loses. [Ablation receipt](../results/organization-force-ablation.json) · [Reproduction script](../scripts/run_organization_ablation.py)'
    for marker,value in {'@@PARAMETERS@@':params,'@@CERTIFICATES@@':certs,'@@REFINEMENT@@':refinements,
        '@@ORGANIZATION@@':org,'@@VERIFICATION@@':verify,'@@CONTINUUM@@':continuous,'@@ABLATION@@':ablated}.items():
        source=source.replace(marker,value)
    assert '@@' not in source
    page('PAID-ORGANIZATION.html',source,'Paid autonomous organization — shared fuel, best constant strategy and verification')
    gap=certificate['spectator_minus_best_fixed_interval']
    summary=f'''## Paid autonomous organization: the next verified benchmark

The system now starts uniformly dispersed and forms a finite-depth attractive interaction mode through reversible passive and fuel-paid reactions. Organization and repair share four fuel/waste molecules. A spectator-fuel reaction provides a local response to inventory without an external controller. [Full model, stochastic calculus, proofs and visuals](PAID-ORGANIZATION.html).

At the prespecified time two, the n={certificate['n']} finite-grid certificate bounds the gap to **the optimum over all constant paid-binding intensities alpha in [0,1]** between **{bound(gap[0],6,scale=100)} and {bound(gap[1],6,True,100)} percentage points**. The n={fine['n']} numerical values are {100*fine['survival']['spectator'][-1]:.6f}% survival for the candidate and {100*fine['survival']['fixed'][-1]:.6f}% for the fixed attaining sample. This is a small effect, not the earlier 23-point free-localization benefit.

![Paid organization comparison and refinements](../results/figures/organization-survival.png)

There are {count} passing tests, {n_paths:,} independent finite-lattice paths, 180,000 separate continuum SDE paths and stopped fuel/organization ledgers. The path precision cannot resolve the small policy gap; the interval certificate proves it for the stated finite grids. A rigorous continuum error enclosure remains open. Removing the attractive well reverses the n=16 comparison, but that changes both force and binding energetics. An open-bath limit gives an exact no-advantage control when the spectator generator equals a fixed-family member.

The comparison is against best constant paid-channel intensity with shared other parameters, not all fixed organizations or molecular architectures. Real molecular recognition/fabrication costs remain unverified. The audit records 29 queries and 26 sources; no new general concept or established publication priority is claimed. [Audit](../research/organization-originality-audit.md) · [Interval/path receipt](../results/organization-verification.json) · [Continuum step check](../results/organization-continuum-check.json).
'''
    (ROOT/'docs/organization-overview.md').write_text(summary,encoding='utf-8')
    journey='''## 13. Paying for organization instead of supplying localization

The next request required an originality audit before building paid autonomous organization. The audit found close precedents for autonomous batch-fuel assembly, reversible regeneration and improvement over constant interactions. The retained claim was narrowed to a reproducible comparison against all constant paid-channel prefactors in a common finite-resource model, not every possible fixed material.

The same dispersed launch, energy landscape, repair chemistry and rate ceilings were fixed before the experiment. A spectator fuel molecule changes paid binding kinetics locally. Both directions of every permitted transition retain their reverse, and motion in the activated conformation follows an attractive potential. No jump teleports particles into contact. Parallel work separated the model, primary-source audit and independent certificate/path verification.

'''+summary.replace('## Paid autonomous organization: the next verified benchmark','### Result and limitations')+'''

The certificate also showed a negative result on the coarsest grid, and the no-well ablation lost to passive-only organization. Short-time derivatives show that paying for formation initially hurts relative to the passive arm. These controls were retained rather than tuning the horizon or well depth to inflate a positive conclusion. The small final finite-grid advantage does not establish a new scientific concept, a continuum theorem or a molecular mechanism. The experiment advances the accounting question—when organization pays for itself—while making its unresolved assumptions explicit.
'''
    (ROOT/'docs/organization-journey.md').write_text(journey,encoding='utf-8')
    target=ROOT/'docs/PAID-ORGANIZATION.html';links=Links();links.feed(target.read_text(encoding='utf-8'))
    broken=[]
    for url in links.urls:
        part=urlsplit(url)
        if part.scheme or part.netloc:continue
        file=target.parent/unquote(part.path) if part.path else target
        if not file.exists():broken.append(url)
        elif file==target and part.fragment and unquote(part.fragment) not in links.ids:broken.append(url)
    assert not broken,broken
    checks=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        for width in [1440,390]:
            tab=browser.new_page(viewport={'width':width,'height':1000});errors=[]
            tab.on('pageerror',lambda e:errors.append(str(e)))
            tab.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
            tab.goto(target.as_uri(),wait_until='networkidle');tab.evaluate('document.fonts.ready')
            item=tab.evaluate('''() => ({equations:document.querySelectorAll('.math').length,
                rendered:document.querySelectorAll('.math .katex').length,
                mathErrors:document.querySelectorAll('.math[data-error]').length,
                brokenImages:[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length,
                width:innerWidth,scrollWidth:document.documentElement.scrollWidth})''')
            item['console_errors']=errors
            assert not errors and not item['mathErrors'] and not item['brokenImages']
            assert item['equations']==item['rendered'] and item['scrollWidth']<=width,item
            tab.screenshot(path=str(ROOT/f'results/document-previews/paid-organization-{width}.png'))
            checks.append(item);tab.close()
        browser.close()
    files=[*ROOT.glob('src/purposeful/*.py'),*ROOT.glob('tests/test_*.py'),
        *ROOT.glob('scripts/*organization*.py'),*ROOT.glob('research/organization*'),
        ROOT/'docs/paid-organization.md',target,ROOT/'docs/organization-overview.md',ROOT/'docs/organization-journey.md',
        data_path,proof_path,ROOT/'results/organization-continuum-check.json',ROOT/'results/organization-force-ablation.json',
        ROOT/'results/organization-certificate-controls.json',
        ROOT/'results/organization-review.json',
        *ROOT.glob('results/figures/organization*')]
    receipt={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS',
        'tests':{'count':count,'return_code':run.returncode,'output':run.stdout+run.stderr},
        'broken_local_links':broken,'browser_checks':checks,
        'sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in files if p.is_file()}}
    (ROOT/'results/organization-doc-verification.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','test_count':count,'browser_checks':checks},indent=2))


if __name__=='__main__':main()
