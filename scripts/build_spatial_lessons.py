"""Render and browser-check the spatial study assessment."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
from playwright.sync_api import sync_playwright
from build_research_docs import page

ROOT=Path(__file__).resolve().parents[1]


def main():
    source=ROOT/'docs/spatial-lessons.md'
    target=ROOT/'docs/SPATIAL-LESSONS.html'
    page(target.name,source.read_text(encoding='utf-8'),'Spatial model: lessons, novelty and next research claim')
    checks=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        for width in [1440,390]:
            tab=browser.new_page(viewport={'width':width,'height':1000})
            errors=[]
            tab.on('pageerror',lambda e:errors.append(str(e)))
            tab.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
            tab.goto(target.as_uri(),wait_until='networkidle')
            tab.evaluate('document.fonts.ready')
            record=tab.evaluate('''() => ({equations:document.querySelectorAll('.math').length,
                rendered:document.querySelectorAll('.math .katex').length,
                mathErrors:document.querySelectorAll('.math[data-error]').length,
                brokenImages:[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length,
                width:innerWidth,scrollWidth:document.documentElement.scrollWidth})''')
            record['console_errors']=errors
            assert record['equations']==record['rendered'] and not record['mathErrors']
            assert not record['brokenImages'] and not errors and record['scrollWidth']<=width
            checks.append(record);tab.close()
        browser.close()
    files=[source,target,Path(__file__)]
    receipt={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS','browser_checks':checks,
        'sha256':{f.relative_to(ROOT).as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}
    (ROOT/'results/spatial-lessons-validation.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
