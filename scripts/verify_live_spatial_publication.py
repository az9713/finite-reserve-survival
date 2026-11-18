"""Verify deployed scientific pages against local bytes and browser rendering.

The publication receipt is written beside the checkout, avoiding a commit/hash
cycle. Run after GitHub Pages reports the requested commit deployed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
BASE='https://az9713.github.io/finite-reserve-survival/'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--commit',required=True)
    parser.add_argument('--include-organization',action='store_true')
    parser.add_argument('--include-continuum',action='store_true')
    parser.add_argument('--include-validated-continuum',action='store_true')
    args=parser.parse_args()
    paths=['index.html','README.md','docs/SCIENTIFIC-DOSSIER.html',
           'docs/DEVELOPMENT-JOURNEY.html','docs/SPATIAL-REPAIR.html',
           'docs/STOCHASTIC-CALCULUS-PLAN.html','results/spatial-repair.json',
           'results/spatial-repair-verification.json','results/spatial-repair-review-check.json',
           'results/spatial-repair-ui-validation.json','results/spatial-repair-integrity.json',
           'results/spatial-repair-delivery-manifest.json','ARTIFACT_MANIFEST.json',
           'research/spatial-repair-originality-audit.md',
           'results/figures/spatial-repair-survival.png',
           'results/figures/spatial-repair-refinement.png',
           'results/figures/spatial-repair-independent.png']
    if args.include_organization:
        paths += ['docs/PAID-ORGANIZATION.html','docs/STOCHASTIC-CALCULUS.html','docs/SPATIAL-LESSONS.html',
            'src/purposeful/organization.py','research/organization-originality-audit.md',
            *['results/'+name+'.json' for name in ['organization','organization-verification',
                'organization-continuum-check','organization-force-ablation','organization-review',
                'organization-doc-verification','organization-integrity','organization-certificate-controls']],
            *['results/figures/organization-'+name+'.png' for name in ['survival','resources','fixed-certificate']]]
    if args.include_continuum:
        paths += ['docs/CONTINUUM-ORGANIZATION.html','docs/continuum-organization.md',
            'research/continuum-coupling.md','research/continuum-residual.md','research/continuum-spectral.md',
            *['results/continuum-'+name+'.json' for name in ['organization-verification','short-time',
                'fourier','weighted','residual-screen']],
            'results/figures/continuum-organization-bounds.png','results/figures/continuum-organization-bounds.pdf']
    if args.include_validated_continuum:
        paths += ['docs/CONTINUUM-CERTIFICATE.html','docs/continuum-certificate.md',
            'docs/CONTINUUM-ORGANIZATION-INITIAL.html','docs/continuum-organization-initial.md',
            'research/validated-residual-implementation.md','research/validated-spatial.md',
            'research/validated-time.md','research/validated-continuum-review.md',
            *['results/continuum-'+name+'.json' for name in ['validated-comparison',
                'validated-independent-replay','validated-verification','time-verification',
                'validated-spatial-tables','reproduction-check','validated-privacy']],
            'results/continuum-validated-spatial-tables.npz',
            'results/continuum-trials/K12-spectator.npz',
            'results/figures/continuum-validated-certificate.png',
            'results/figures/continuum-validated-certificate.pdf',
            'results/figures/continuum-certificate-pipeline.svg',
            'results/continuum-validated-nodes/spectator.json',
            *[f'results/continuum-validated-nodes/node{j:02d}.json' for j in range(38)]]
    paths=list(dict.fromkeys(paths))
    files=[]
    for relative in paths:
        request=Request(BASE+relative+'?publication='+args.commit,headers={'User-Agent':'Research-publication-verifier'})
        with urlopen(request,timeout=60) as response:
            payload=response.read();status=response.status
        actual=hashlib.sha256(payload).hexdigest()
        expected=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        assert status==200 and actual==expected, (relative,status,actual,expected)
        files.append({'path':relative,'status':status,'sha256':actual,'matches_local':True})
    pages=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        documents=['','docs/DEVELOPMENT-JOURNEY.html','docs/SPATIAL-REPAIR.html']
        if args.include_organization:documents+=['docs/PAID-ORGANIZATION.html','docs/STOCHASTIC-CALCULUS.html']
        if args.include_continuum:documents+=['docs/CONTINUUM-ORGANIZATION.html']
        if args.include_validated_continuum:documents+=['docs/CONTINUUM-CERTIFICATE.html']
        for relative in documents:
            for width in [1440,390]:
                page=browser.new_page(viewport={'width':width,'height':1000})
                errors=[]
                page.on('pageerror',lambda e:errors.append(str(e)))
                page.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
                page.goto(BASE+relative+'?publication='+args.commit,wait_until='networkidle')
                page.evaluate('document.fonts.ready')
                state=page.evaluate('''() => ({url:location.href,
                    equations:document.querySelectorAll('.math').length,
                    rendered:document.querySelectorAll('.math .katex').length,
                    mathErrors:document.querySelectorAll('.math[data-error]').length,
                    brokenImages:[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length,
                    hasSpatialResults:document.body.innerText.includes('23.208'),
                    width:innerWidth,scrollWidth:document.documentElement.scrollWidth})''')
                state['console_errors']=errors
                assert not errors and not state['mathErrors'] and not state['brokenImages'],state
                assert state['equations']==state['rendered'],state
                if relative not in ['docs/PAID-ORGANIZATION.html','docs/STOCHASTIC-CALCULUS.html','docs/CONTINUUM-ORGANIZATION.html','docs/CONTINUUM-CERTIFICATE.html']:
                    assert state['hasSpatialResults'],state
                if args.include_continuum and relative in ['', 'docs/DEVELOPMENT-JOURNEY.html','docs/CONTINUUM-ORGANIZATION.html']:
                    assert page.get_by_text('unresolved',exact=False).count()>0,state
                if args.include_validated_continuum and relative in ['', 'docs/DEVELOPMENT-JOURNEY.html','docs/CONTINUUM-CERTIFICATE.html']:
                    assert '0.0110051' in page.locator('body').inner_text(),state
                assert state['scrollWidth']<=width,state
                pages.append(state);page.close()
        browser.close()
    receipt={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS',
             'expected_commit':args.commit,'base_url':BASE,'files':files,'browser_checks':pages,
             'scope':'Live served bytes and rendering. Deployment commit is checked separately with GitHub Pages build metadata.'}
    output=ROOT.parent/'finite-reserve-survival-publication.json'
    output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','files_matched':len(files),'browser_checks':len(pages),'receipt':str(output)},indent=2))


if __name__=='__main__':main()
