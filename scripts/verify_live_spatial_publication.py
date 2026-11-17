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
        for relative in ['','docs/DEVELOPMENT-JOURNEY.html','docs/SPATIAL-REPAIR.html']:
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
                assert state['equations']==state['rendered'] and state['hasSpatialResults'],state
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
