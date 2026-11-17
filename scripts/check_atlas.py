"""Optional browser smoke check: requires playwright and Chromium."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright

root=Path(__file__).resolve().parents[1]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--no-sandbox'])
    page=browser.new_page(viewport={'width':1440,'height':1100},device_scale_factor=1)
    errors=[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    # Inline loading verifies the standalone artifact without a server or file-URL policy.
    page.set_content((root/'results/atlas.html').read_text())
    page.locator('#best').click()
    assert 'survival=' in page.locator('#inspection').inner_text()
    inspected={}
    for law in ['hazard','damage']:
        page.locator('#law').select_option(law)
        page.locator('#best').click()
        inspected[law]=page.locator('#inspection').inner_text()
    page.screenshot(path=str(root/'results/atlas-desktop.png'),full_page=True)
    page.set_viewport_size({'width':390,'height':844})
    page.wait_for_timeout(100)
    assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
    page.screenshot(path=str(root/'results/atlas-mobile.png'),full_page=True)
    assert not errors, errors
    report={'desktop':[1440,1100],'mobile':[390,844],'javascript_errors':errors,
            'mobile_horizontal_overflow':False,'inspected':inspected}
    (root/'results/browser-check.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    browser.close()
