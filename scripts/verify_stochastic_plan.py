"""Check the plan's provenance, local links, equations, figures and layouts."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
from urllib.parse import urlsplit, unquote

from playwright.sync_api import sync_playwright
from verify_research_docs import Links

ROOT=Path(__file__).resolve().parents[1]


def main():
    document=ROOT/'docs/STOCHASTIC-CALCULUS-PLAN.html'
    data=json.loads((ROOT/'results/encounter-feasibility.json').read_text(encoding='utf-8'))
    receipt=json.loads((ROOT/'results/encounter-verification.json').read_text(encoding='utf-8'))
    test_match=re.search(r'Ran (\d+) tests?',receipt['test_output'])
    assert receipt['return_code']==0 and test_match and int(test_match[1])>=42
    test_count=int(test_match[1])
    assert receipt['data_sha256']==hashlib.sha256((ROOT/'results/encounter-feasibility.json').read_bytes()).hexdigest()
    mismatches=[name for name,digest in data['source_sha256'].items()
                if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest]
    assert not mismatches,mismatches
    links=Links()
    links.feed(document.read_text(encoding='utf-8'))
    broken=[]
    for url in links.urls:
        parts=urlsplit(url)
        if parts.scheme or parts.netloc:continue
        target=document.parent/unquote(parts.path) if parts.path else document
        if not target.exists():broken.append(url)
        elif target==document and parts.fragment and unquote(parts.fragment) not in links.ids:broken.append(url)
    assert not broken,broken
    records=[]
    previews=ROOT/'results/document-previews'
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        for width in [1440,390]:
            page=browser.new_page(viewport={'width':width,'height':1000})
            errors=[]
            page.on('pageerror',lambda error:errors.append(str(error)))
            page.on('console',lambda message:errors.append(message.text) if message.type=='error' else None)
            page.goto(document.as_uri(),wait_until='networkidle')
            page.evaluate('document.fonts.ready')
            record=page.evaluate('''() => ({equations:document.querySelectorAll('.math').length,
                rendered:document.querySelectorAll('.math .katex').length,
                mathErrors:document.querySelectorAll('.math[data-error]').length,
                brokenImages:[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length,
                scrollWidth:document.documentElement.scrollWidth,width:document.documentElement.clientWidth})''')
            record.update(viewport=width,console_errors=errors)
            assert record['equations']==record['rendered']>0
            assert not record['mathErrors'] and not record['brokenImages'] and not errors
            assert record['scrollWidth']<=width
            page.screenshot(path=str(previews/f'stochastic-plan-{width}.png'))
            records.append(record)
            page.close()
        browser.close()
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'source_hash_mismatches':mismatches,
            'broken_local_links':broken,'browser_checks':records,'test_count':test_count,
            'scope':'Classical encounter benchmark and historical implementation plan. The implemented maintenance model has its own report and verification.'}
    (ROOT/'results/stochastic-plan-ui-validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
