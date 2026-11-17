"""Render the spatial study and retain browser, links and provenance receipts."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
from urllib.parse import unquote, urlsplit

import numpy as np
from playwright.sync_api import sync_playwright
from build_research_docs import page
from verify_research_docs import Links

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    data_path=ROOT/'results/spatial-repair.json'
    verify_path=ROOT/'results/spatial-repair-verification.json'
    data=json.loads(data_path.read_text(encoding='utf-8'))
    verify=json.loads(verify_path.read_text(encoding='utf-8'))
    assert verify['passed'] is True, 'Independent verification did not pass.'
    assert verify['input_sha256']==sha(data_path), 'Verification used different experiment data.'
    for group in ['source_sha256','figure_sha256']:
        for relative,digest in verify[group].items():
            assert sha(ROOT/relative)==digest, f'Stale verification input or output: {relative}'
    for relative,digest in data['source_sha256'].items():
        assert sha(ROOT/relative)==digest, f'Stale computation source: {relative}'
    fine=data['refinements'][-1]
    source=(ROOT/'docs/spatial-repair.md').read_text(encoding='utf-8')
    params='| Parameter | Value | Meaning |\n|---|---:|---|\n'
    meanings={'fuel':'Conserved fuel plus waste inventory', 'eta':'Forward damage rate per healthy complex',
        'diffusivity':'Relative diffusivity D1 + D2', 'length':'Periodic box side',
        'radius':'Fixed physical gate radius in model length units',
        'activity':'Intrinsic bounded contact activity κ', 'local_width':'Localized preparation σ'}
    for key,value in data['parameters'].items():params+=f'| {key} | {value:g} | {meanings[key]} |\n'
    table='| Initial chemistry | Spatial start | Survival at t = 4 | Edge-calibrated count at t = 4 | Spatial minus count at sampled maximum absolute gap | Time of maximum |\n|---|---|---:|---:|---:|---:|\n'
    for initial,case in fine['initial_conditions'].items():
        for law,values in case['spatial'].items():
            summary=case['comparisons_to_edgewise'][law]
            table+=f'| {initial}, full fuel | {law} | {100*values[-1]:.5f}% | {100*case["count_edgewise"][-1]:.5f}% | {100*summary["signed_gap_at_sampled_max"]:+.5f} percentage points | {summary["time_of_sampled_max"]:.4f} |\n'
    refinement='| Grid comparison | HD uniform maximum change | HD localized maximum change | HD Palm maximum change |\n|---|---:|---:|---:|\n'
    for item in data['refinement_differences']:
        v=item['max_curve_differences']['HD']
        refinement+=f'| {item["grids"][0]} → {item["grids"][1]} | {v["uniform"]:.8f} | {v["localized"]:.8f} | {v["reaction"]:.8f} |\n'
    refinement+='\nThese entries are absolute probability differences, not percentages. Each n denotes the full periodic lattice points per axis before exact symmetry reduction. The finest grid has '+str(fine['spatial_orbits'])+' spatial orbits.\n'
    evidence='The independent verification receipt below retains test output, deterministic diagnostics, continuous-space path comparisons, random seeds, source hashes, and explicit acceptance checks.\n\n'
    evidence+='[Open machine-readable verification receipt](../results/spatial-repair-verification.json) · [Verification design and derivations](../research/spatial-repair-verification-design.md)\n\n'
    # Do not summarize unknown receipt schemas as a successful verification.
    evidence+='```json\n'+json.dumps(verify.get('summary',verify.get('status',{})),indent=2)+'\n```\n'
    evidence+='\n| Continuous-space launch | Maximum MC–lattice discrepancy | 99% simultaneous lower bound on maximum gap from fixed calibrated count curve |\n|---|---:|---:|\n'
    for case in verify['continuous_path_comparison']['cases']:
        bound=case['comparison_to_fixed_finest_grid_calibrated_count']['simultaneous_lower_bound_on_continuous_time_maximum_absolute_gap']
        evidence+=f'| {case["chemistry"]}, {case["preparation"]} | {case["maximum_dkw_check_discrepancy"]:.6f} | {100*bound:.3f} percentage points |\n'
    evidence+=('\nThese confidence lower bounds compare the true Brownian survival law with the **fixed count coefficients computed on the finest grid**. They do not enclose error in calibrating the exact continuum encounter mean. The localized HD lower bound exceeds 22 percentage points; the uniform HD lower bound is below the proposed one-point screening threshold. The deterministic uniform gap of 1.340 points should therefore not be presented as a statistically verified continuum gap above that threshold.\n')
    evidence+='\nA [supplemental fresh-seed review receipt](../results/spatial-repair-review-check.json) retains two further 100,000-path HD runs, their replay program, ledgers, sampled curves and source hashes. They are separate from the formal 200,000-path verification above.\n'
    hd=fine['initial_conditions']['HD']
    localized=hd['comparisons_to_edgewise']['localized']
    uniform=hd['comparisons_to_edgewise']['uniform']
    conclusion=(f'For the stated HD/full-fuel case, the largest sampled deviation from the edge-calibrated count chain is **{100*uniform["sampled_max_abs_gap"]:.3f} percentage points** for a uniform start and **{100*localized["sampled_max_abs_gap"]:.3f} percentage points** for a localized start. At time four the localized spatial survival is **{100*hd["spatial"]["localized"][-1]:.3f}%**, compared with **{100*hd["count_edgewise"][-1]:.3f}%** for the count chain. This is a substantial preparation effect within this exploratory parameter choice. Matching all uniform-start forward means and preserving equilibrium ratios does not retain the preparation information. The second reference calibration and HH launch are reported to make the scope assessable.\n\n'
        'The result is consistent with known localization and rebinding effects. Finite fuel, the chosen first-collapse event, and the replayable comparison define this benchmark; they do not by themselves establish originality. Rigorous continuum error enclosures and an experimentally calibrated molecular mechanism remain open.')
    for marker,value in {'@@PARAMETERS@@':params,'@@RESULT_TABLE@@':table,
        '@@REFINEMENT_TABLE@@':refinement,'@@VERIFICATION@@':evidence,'@@CONCLUSION@@':conclusion,
        '@@MC_BAND@@':f'The receipt uses {verify["summary"]["monte_carlo_paths_per_preparation"]:,} paths in each of four preparations and familywise confidence {100*verify["summary"]["confidence"]:g}%. Its simultaneous half-width is {verify["summary"]["dkw_half_width"]:.6f} in absolute probability.'}.items():
        source=source.replace(marker,value)
    assert '@@' not in source
    page('SPATIAL-REPAIR.html',source,'Two-complex spatial repair: model, proofs and verification')
    target=ROOT/'docs/SPATIAL-REPAIR.html'
    links=Links(); links.feed(target.read_text(encoding='utf-8'))
    broken=[]
    for url in links.urls:
        parts=urlsplit(url)
        if parts.scheme or parts.netloc:continue
        file=target.parent/unquote(parts.path) if parts.path else target
        if not file.exists():broken.append(url)
        elif file==target and parts.fragment and unquote(parts.fragment) not in links.ids:broken.append(url)
    assert not broken,broken
    views=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        for width in [1440,390]:
            tab=browser.new_page(viewport={'width':width,'height':1000})
            errors=[]
            tab.on('pageerror',lambda e:errors.append(str(e)))
            tab.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
            tab.goto(target.as_uri(),wait_until='networkidle')
            tab.evaluate('document.fonts.ready')
            stats=tab.evaluate('''() => ({equations:document.querySelectorAll('.math').length,
                rendered:document.querySelectorAll('.math .katex').length,
                mathErrors:document.querySelectorAll('.math[data-error]').length,
                brokenImages:[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length,
                scrollWidth:document.documentElement.scrollWidth,width:innerWidth})''')
            stats['console_errors']=errors
            assert stats['equations']==stats['rendered'] and stats['mathErrors']==0
            assert stats['brokenImages']==0 and not errors
            assert stats['scrollWidth']<=width
            tab.screenshot(path=str(ROOT/f'results/document-previews/spatial-repair-{width}.png'))
            views.append(stats);tab.close()
        browser.close()
    paths=[data_path,verify_path,target,ROOT/'docs/spatial-repair.md',
        ROOT/'scripts/build_spatial_repair_docs.py',
        *sorted((ROOT/'results/figures').glob('spatial-repair-*'))]
    record={'generated_utc':datetime.now(timezone.utc).isoformat(),
        'status':'PASS','broken_local_links':broken,'browser_checks':views,
        'source_and_artifact_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in paths}}
    (ROOT/'results/spatial-repair-ui-validation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
