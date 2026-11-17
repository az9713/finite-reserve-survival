"""Browser, link, PDF, and receipt checks for the three research documents."""
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urlsplit, unquote
from html.parser import HTMLParser
import hashlib
import json
import re

import fitz
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
RESULTS = ROOT / 'results'


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        for key in ['href', 'src']:
            if key in attrs:
                self.urls.append(attrs[key])
        if 'id' in attrs:
            self.ids.add(attrs['id'])


def main():
    shots = RESULTS / 'document-previews'
    shots.mkdir(exist_ok=True)
    names = ['DEVELOPMENT-JOURNEY.html', 'SCIENTIFIC-DOSSIER.html', 'RESEARCH-MANUSCRIPT.html']
    broken = []
    for name in names:
        document = DOCS / name
        parser = Links()
        parser.feed(document.read_text(encoding='utf-8'))
        for url in parser.urls:
            parts = urlsplit(url)
            if parts.scheme or parts.netloc:
                continue
            target = document.parent / unquote(parts.path) if parts.path else document
            if not target.exists():
                broken.append([name, url])
            elif target == document and parts.fragment and unquote(parts.fragment) not in parser.ids:
                broken.append([name, url])
    browser_records = []
    with sync_playwright() as p:
        browser = p.chromium.launch(channel='chrome', headless=True)
        for name in names:
            for width in [1440, 390]:
                page = browser.new_page(viewport={'width': width, 'height': 1000}, device_scale_factor=1)
                errors = []
                page.on('pageerror', lambda error: errors.append(str(error)))
                page.on('console', lambda message: errors.append(message.text) if message.type == 'error' else None)
                page.goto((DOCS / name).as_uri(), wait_until='networkidle')
                page.evaluate('document.fonts.ready')
                data = page.evaluate('''() => ({
                    equations: document.querySelectorAll('.math').length,
                    rendered: document.querySelectorAll('.math .katex').length,
                    mathErrors: [...document.querySelectorAll('.math[data-error]')].map(e=>e.dataset.error),
                    images: document.images.length,
                    brokenImages: [...document.images].filter(e=>!e.complete||e.naturalWidth===0).map(e=>e.src),
                    width: document.documentElement.clientWidth,
                    scrollWidth: document.documentElement.scrollWidth,
                    title: document.title
                })''')
                data.update(file=name, viewport=width, console_errors=errors)
                browser_records.append(data)
                page.screenshot(path=str(shots / f'{Path(name).stem.lower()}-{width}.png'))
                if name == 'SCIENTIFIC-DOSSIER.html' and width == 1440:
                    page.get_by_role('heading', name='Main results and the logic connecting them', exact=True).scroll_into_view_if_needed()
                    page.screenshot(path=str(shots / 'scientific-theorems.png'))
                page.close()
        browser.close()
    pdf_path = RESULTS / 'research-manuscript.pdf'
    pdf = fitz.open(pdf_path)
    pdf_pages = len(pdf)
    missing_text = [i + 1 for i, page in enumerate(pdf) if not page.get_text().strip()]
    for number in [0, min(3, pdf_pages - 1), pdf_pages - 1]:
        pdf[number].get_pixmap(matrix=fitz.Matrix(1.3, 1.3)).save(shots / f'manuscript-page-{number+1}.png')
    scientific_receipt = json.loads((RESULTS / 'dossier-verification.json').read_text(encoding='utf-8'))
    mismatches = [path for path, digest in scientific_receipt['sha256'].items()
                  if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != digest]
    record = {'generated_utc': datetime.now(timezone.utc).isoformat(), 'broken_local_links': broken,
              'browser_checks': browser_records, 'pdf_pages': pdf_pages, 'pdf_blank_text_pages': missing_text,
              'scientific_receipt_hash_mismatches': mismatches,
              'external_links': 'Not all external links are health-checked; primary sources were inspected during literature work.'}
    (RESULTS / 'dossier-ui-validation.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    assert not broken, broken
    assert not mismatches, mismatches
    for item in browser_records:
        assert not item['console_errors'], item
        assert not item['mathErrors'] and not item['brokenImages'], item
        assert item['equations'] == item['rendered'], item
        assert item['scrollWidth'] <= item['width'] + 1, item
    assert not missing_text
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
