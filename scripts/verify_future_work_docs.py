"""Verify the planning document, without asserting its proposed research results."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import unquote,urlsplit
from playwright.sync_api import sync_playwright
from verify_research_docs import Links

ROOT=Path(__file__).resolve().parents[1]


def main():
    source=ROOT/'docs/future-work.md';html=ROOT/'docs/FUTURE-WORK.html'
    projects=re.findall(r'^### (FW\d{2}) — (.+)$',source.read_text(encoding='utf-8'),re.M)
    assert [p[0] for p in projects]==[f'FW{i:02d}' for i in range(1,25)]
    parser=Links();parser.feed(html.read_text(encoding='utf-8'));broken=[]
    for url in parser.urls:
        parts=urlsplit(url)
        if parts.scheme or parts.netloc:continue
        target=html.parent/unquote(parts.path) if parts.path else html
        if not target.exists():broken.append(url);continue
        if parts.fragment and target.suffix.lower()=='.html':
            other=Links();other.feed(target.read_text(encoding='utf-8'))
            if unquote(parts.fragment) not in other.ids:broken.append(url)
    assert not broken,broken
    checks=[];shots=ROOT/'results/document-previews';shots.mkdir(exist_ok=True)
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        for width in (1440,390):
            tab=browser.new_page(viewport={'width':width,'height':1000})
            errors=[];tab.on('pageerror',lambda e:errors.append(str(e)))
            tab.on('console',lambda msg:errors.append(msg.text) if msg.type=='error' else None)
            tab.goto(html.as_uri(),wait_until='networkidle');tab.evaluate('document.fonts.ready')
            data=tab.evaluate('''() => ({
              equations:document.querySelectorAll('.math').length,
              rendered:document.querySelectorAll('.math .katex').length,
              mathErrors:[...document.querySelectorAll('.math[data-error]')].map(e=>e.dataset.error),
              brokenImages:[...document.images].filter(e=>!e.complete||!e.naturalWidth).map(e=>e.src),
              width:document.documentElement.clientWidth,scrollWidth:document.documentElement.scrollWidth,
              projectHeadings:[...document.querySelectorAll('h3')].map(e=>e.textContent).filter(t=>/^FW[0-9]{2}/.test(t))})''')
            data.update(viewport=width,console_errors=errors)
            assert len(data['projectHeadings'])==24
            assert data['equations']==data['rendered'] and not data['mathErrors']
            assert not data['brokenImages'] and not errors
            assert data['scrollWidth']<=width
            assert 'proposed' in tab.locator('body').inner_text().lower()
            tab.screenshot(path=str(shots/f'future-work-{width}.png'))
            tab.get_by_role('heading',name=re.compile('^FW11')).scroll_into_view_if_needed()
            tab.screenshot(path=str(shots/f'future-work-self-construction-{width}.png'))
            if width==1440:
                tab.locator('img').first.screenshot(path=str(shots/'future-work-map.png'))
            checks.append(data);tab.close()
        browser.close()
    files=['docs/future-work.md','docs/FUTURE-WORK.html','research/future-work-sources.md',
           'results/figures/future-work-map.svg','scripts/build_future_work_docs.py',
           'scripts/verify_future_work_docs.py','scripts/build_research_docs.py',
           'scripts/integrate_continuum_certificate.py','README.md']
    record={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'FUTURE_WORK_DOCUMENTATION_VERIFIED',
        'scope':'Planning content, topic enumeration, local links, offline equations and layout. No proposed research outcome is asserted or experimentally tested.',
        'project_count':len(projects),'projects':[{'id':i,'title':t,'state':'PROPOSED_FOLLOWUP'} for i,t in projects],
        'broken_local_links':broken,'browser_checks':checks,
        'baseline_certificate_sha256':sha256((ROOT/'results/continuum-validated-comparison.json').read_bytes()).hexdigest(),
        'source_and_artifact_sha256':{f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files}}
    (ROOT/'results/future-work-verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'status':record['status'],'projects':len(projects),'browser_checks':len(checks),
                      'equations_per_view':checks[0]['equations'],'broken_links':broken},indent=2))


if __name__=='__main__':main()
