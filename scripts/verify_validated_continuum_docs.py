"""Verify the final comparison, offline equations, links, figures and artifacts."""
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from urllib.parse import unquote,urlsplit
import json
import platform
import re
import sys
import numpy as np
import scipy
import flint
from playwright.sync_api import sync_playwright
from verify_research_docs import Links

ROOT=Path(__file__).resolve().parents[1]


def main():
    comp=json.loads((ROOT/'results/continuum-validated-comparison.json').read_text())
    replay=json.loads((ROOT/'results/continuum-validated-independent-replay.json').read_text())
    cover=comp['parameter_cover']
    assert replay['parameter_cover']['spectator_minus_best_fixed_lower']==cover['spectator_minus_best_fixed_lower']
    assert len(comp['node_records'])==38
    if comp['status']=='STRICT_CONTINUUM_ADVANTAGE_CERTIFIED':
        assert F(cover['spectator_minus_best_fixed_lower'])>0
        assert replay['parameter_cover']['strict_gap_arithmetically_verified']
    tests=(ROOT/'results/continuum-validated-full-tests.txt').read_text(encoding='utf-8-sig')
    count=re.search(r'Ran (\d+) tests',tests)
    assert count and re.search(r'^OK\s*$',tests,re.M)
    privacy_path=ROOT/'results/continuum-validated-privacy.json'
    if privacy_path.exists():
        privacy=json.loads(privacy_path.read_text())
        assert privacy['status']=='PASS' and not privacy['findings'],privacy
    half_tiny=float(np.finfo(float).tiny)/2
    probe_a=np.zeros((8,455));probe_a[:,0]=half_tiny
    probe_b=np.ones((455,64))
    arithmetic_probes={
        'binary64_round_to_nearest_ties_sample':float(np.float64(1)+np.float64(2**-53))==1.,
        'scalar_gradual_underflow_sample':float(np.float64(np.finfo(float).tiny)*np.float64(.5))==half_tiny and half_tiny>0,
        'blas_gradual_underflow_sample':bool(np.all(probe_a@probe_b==half_tiny))}
    assert all(arithmetic_probes.values()),arithmetic_probes
    historical=json.loads((ROOT/'results/continuum-organization-verification.json').read_text())
    relocations={'docs/CONTINUUM-ORGANIZATION.html':'docs/CONTINUUM-ORGANIZATION-INITIAL.html',
                 'docs/continuum-organization.md':'docs/continuum-organization-initial.md'}
    historical_replay=[]
    for name,expected in historical['sha256'].items():
        resolved=relocations.get(name,name)
        assert sha256((ROOT/resolved).read_bytes()).hexdigest()==expected,('historical artifact mismatch',name,resolved)
        historical_replay.append({'original_path':name,'retained_path':resolved,'sha256':expected})
    target=ROOT/'results/continuum-validated-verification.json'
    target.write_text('{}\n') # Make the self-referencing receipt link exist during checking.
    names=['CONTINUUM-CERTIFICATE.html','CONTINUUM-ORGANIZATION.html',
           'SCIENTIFIC-DOSSIER.html','DEVELOPMENT-JOURNEY.html']
    records=[];broken=[]
    shots=ROOT/'results/document-previews';shots.mkdir(exist_ok=True)
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        for name in names:
            file=ROOT/'docs'/name
            links=Links();links.feed(file.read_text(encoding='utf-8'))
            for url in links.urls:
                parts=urlsplit(url)
                if parts.scheme or parts.netloc: continue
                dst=file.parent/unquote(parts.path) if parts.path else file
                if not dst.exists(): broken.append([name,url])
                elif dst==file and parts.fragment and unquote(parts.fragment) not in links.ids: broken.append([name,url])
            for width in (1440,390):
                tab=browser.new_page(viewport={'width':width,'height':1000},device_scale_factor=1)
                errors=[]
                tab.on('pageerror',lambda error:errors.append(str(error)))
                tab.on('console',lambda msg:errors.append(msg.text) if msg.type=='error' else None)
                tab.goto(file.as_uri(),wait_until='networkidle');tab.evaluate('document.fonts.ready')
                item=tab.evaluate('''() => ({equations:document.querySelectorAll('.math').length,
                  rendered:document.querySelectorAll('.math .katex').length,
                  mathErrors:[...document.querySelectorAll('.math[data-error]')].map(e=>e.dataset.error),
                  brokenImages:[...document.images].filter(e=>!e.complete||!e.naturalWidth).map(e=>e.src),
                  width:document.documentElement.clientWidth,scrollWidth:document.documentElement.scrollWidth})''')
                item.update(file=name,viewport=width,console_errors=errors)
                records.append(item)
                assert item['equations']==item['rendered'] and not item['mathErrors'],item
                assert not errors and not item['brokenImages'],item
                assert item['scrollWidth']<=item['width']+1,item
                tab.screenshot(path=str(shots/f'{Path(name).stem.lower()}-{width}.png'))
                if name=='CONTINUUM-CERTIFICATE.html' and width==1440:
                    tab.locator('img').first.scroll_into_view_if_needed()
                    tab.screenshot(path=str(shots/'continuum-certificate-figure.png'))
                tab.close()
        browser.close()
    assert not broken,broken
    paths=set()
    for pattern in ('src/purposeful/*.py','tests/test_*.py',
                    'scripts/*validated*continuum*.py','scripts/validate_continuum*.py',
                    'scripts/build_validated_spatial_tables.py','research/validated*.md','research/continuum-*.md',
                    'results/continuum-*.json','results/continuum-*tests.txt',
                    'results/continuum-validated-nodes/*.json','results/continuum-trials/*.json',
                    'results/continuum-trials/*.npz','results/figures/continuum-validated-certificate.*',
                    'results/figures/continuum-certificate-pipeline.svg','results/continuum-time-tests.log',
                    'results/continuum-validated-*-progress.txt'):
        paths.update(ROOT.glob(pattern))
    paths.update(ROOT/p for p in ['README.md','pyproject.toml','docs/continuum-certificate.md',
        'docs/continuum-overview.md','docs/continuum-journey.md','results/continuum-validated-spatial-tables.npz',
        'docs/continuum-organization.md','docs/continuum-organization-initial.md',
        'docs/CONTINUUM-ORGANIZATION-INITIAL.html','results/research-manuscript.pdf',
        'scripts/build_research_docs.py','scripts/integrate_continuum_certificate.py','scripts/verify_research_docs.py'])
    paths.update(ROOT/'docs'/name for name in names)
    paths.discard(target);paths.discard(ROOT/'results/continuum-profile-local.json')
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'VERIFIED_LOCAL_CONTINUUM_DELIVERY',
        'comparison_status':comp['status'],'test_count':int(count.group(1)),
        'test_output':'results/continuum-validated-full-tests.txt','broken_local_links':broken,
        'browser_checks':records,'arithmetic_runtime_probes':arithmetic_probes,'environment':{'python':platform.python_version(),'numpy':np.__version__,
          'scipy':scipy.__version__,'python_flint':flint.__version__,'binary64_epsilon':float(np.finfo(float).eps)},
        'historical_continuum_receipt_replay':historical_replay,
        'source_and_artifact_sha256':{p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        'limits':'Computer-assisted enclosure under the stated IEEE754/Arb arithmetic contract. Independent parameter and receipt-chain replay; spatial proof-assistant formalization and external expert review not supplied. Local browser verification does not assert GitHub publication.'}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_and_artifact_sha256'},indent=2))


if __name__=='__main__':main()
