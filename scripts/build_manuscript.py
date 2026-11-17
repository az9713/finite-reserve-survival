"""Compile the candidate manuscript and package portable LaTeX sources."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import zipfile

import pypandoc

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
RESULTS = ROOT / 'results'


def main():
    source = (DOCS / 'research-manuscript.md').read_text(encoding='utf-8')
    title = source.splitlines()[0].removeprefix('# ')
    source = source.split('## Abstract', 1)[1]
    abstract, body = source.split('## Motivation', 1)
    body = '## Motivation' + body
    body = re.sub(r'^(#{2,6}) \d+\. ', r'\1 ', body, flags=re.M)
    body = body.replace('2−√2', r'\(2-\sqrt2\)').replace('κ', r'\(\kappa\)')
    body = body.replace('| Environment | Sense | Control | Harvest | Wear | Terminal damage | Repair |',
                        '| Env. | Sense | Control | Harvest | Wear | Failure | Repair |')
    # Local links stay useful in the HTML. The PDF names them as supplement paths.
    body = re.sub(r'\[([^\]]+)\]\((?!https?://)([^)]+)\)',
                  lambda m: m.group(0) if m.start() and body[m.start()-1] == '!' else
                  f'{m.group(1)} (`{m.group(2).removeprefix("../")}`)', body)
    metadata = {'title': title,
                'author': 'AI-assisted working draft; human author attribution pending',
                'date': '17 November 2026',
                'abstract': abstract.strip()}
    with tempfile.TemporaryDirectory(prefix='survival-paper-') as staging:
        stage = Path(staging)
        (stage / 'metadata.json').write_text(json.dumps(metadata), encoding='utf-8')
        tex = pypandoc.convert_text(
            body, to='latex', format='markdown+tex_math_single_backslash',
            extra_args=['--standalone', '--number-sections', '--shift-heading-level-by=-1',
                        '--metadata-file=' + str(stage / 'metadata.json'),
                        '-V', 'documentclass=article', '-V', 'fontsize=11pt',
                        '-V', 'geometry:margin=1in',
                        '-V', 'colorlinks=true', '-V', 'linkcolor=blue', '-V', 'urlcolor=blue'])
        # Portable source archive places all figures beside the manuscript.
        tex = tex.replace('\r', '').replace('../results/figures/', 'figures/')
        tex = tex.replace(r'\real{0.1429}', r'\real{0.1428}')
        # Path strings in tables must allow line breaks at punctuation.
        tex = re.sub(r'\\texttt\{([^{}]*[/][^{}]*)\}',
                     lambda m: r'\path{' + m.group(1).replace(r'\_', '_') + '}', tex)
        (stage / 'manuscript.tex').write_text(tex, encoding='utf-8')
        (stage / 'figures').mkdir()
        used = sorted(set(re.findall(r'figures/([^}]+)', tex)))
        for name in used:
            shutil.copy2(RESULTS / 'figures' / name, stage / 'figures' / name)
        run = subprocess.run(['tectonic', '--keep-logs', 'manuscript.tex'], cwd=stage,
                             capture_output=True, text=True, encoding='utf-8', errors='replace')
        log = run.stdout + run.stderr
        (RESULTS / 'manuscript-build.log').write_text(log, encoding='utf-8')
        if run.returncode:
            raise RuntimeError(log)
        shutil.copy2(stage / 'manuscript.pdf', RESULTS / 'research-manuscript.pdf')
        # The readable workspace copy resolves figures relative to docs/.
        (DOCS / 'research-manuscript.tex').write_text(
            tex.replace('figures/', '../results/figures/'), encoding='utf-8')
        (stage / 'README.txt').write_text(
            'Candidate research draft; priority and independent review pending.\n'
            'Build: tectonic manuscript.tex (or XeLaTeX, twice).\n'
            'Fonts: standard Latin Modern supplied with TeX.\n'
            'Human authorship and responsibility must be settled before submission.\n'
            'Computational supplements are in the accompanying project package.\n', encoding='utf-8')
        with zipfile.ZipFile(RESULTS / 'research-manuscript-source.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
            for path in [stage / 'manuscript.tex', stage / 'README.txt', *sorted((stage / 'figures').glob('*'))]:
                archive.write(path, path.relative_to(stage))
        record = {'generated_utc': datetime.now(timezone.utc).isoformat(),
                  'compiler': 'tectonic', 'return_code': run.returncode,
                  'figures': used,
                  'warnings': [line for line in log.splitlines() if 'warning' in line.lower()],
                  'pdf_sha256': hashlib.sha256((RESULTS / 'research-manuscript.pdf').read_bytes()).hexdigest()}
        (RESULTS / 'manuscript-build.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
        print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
