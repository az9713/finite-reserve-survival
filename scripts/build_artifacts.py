"""Assemble offline HTML, verification note, source provenance, and PDF.

PDF export additionally requires reportlab; simulation does not.
"""
from pathlib import Path
import hashlib
import html
import json
import re
import sys
ROOT=Path(__file__).resolve().parents[1]


def make_pdf():
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Preformatted, Image, Table, TableStyle, PageBreak
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    import matplotlib
    fonts=Path(matplotlib.get_data_path())/'fonts/ttf'
    for name,filename in [('ResearchSans','DejaVuSans.ttf'),('ResearchSansBold','DejaVuSans-Bold.ttf'),('ResearchMono','DejaVuSansMono.ttf')]:
        pdfmetrics.registerFont(TTFont(name,str(fonts/filename)))
    pdfmetrics.registerFontFamily('ResearchSans',normal='ResearchSans',bold='ResearchSansBold',italic='ResearchSans',boldItalic='ResearchSansBold')
    styles=getSampleStyleSheet()
    styles['Title'].fontName='ResearchSansBold'
    styles['Heading2'].fontName='ResearchSansBold'
    styles.add(ParagraphStyle(name='BodyResearch',parent=styles['BodyText'],fontName='ResearchSans',fontSize=9,leading=13.5,spaceAfter=7))
    styles.add(ParagraphStyle(name='CodeResearch',fontName='ResearchMono',fontSize=7,leading=10,spaceAfter=8))
    styles.add(ParagraphStyle(name='TableResearch',fontName='ResearchSans',fontSize=7.3,leading=10))
    story=[]
    def inline(s):
        s=html.escape(s)
        s=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#087e82">\1</link>',s)
        s=re.sub(r'\[([^\]]+)\]\([^)]+\)',r'\1',s)
        s=re.sub(r'\*\*(.*?)\*\*',r'<b>\1</b>',s)
        s=re.sub(r'`([^`]+)`',r'<font name="ResearchMono">\1</font>',s)
        return s
    for docindex,name in enumerate(['docs/paper.md','docs/model.md','research/literature.md','results/verification.md']):
        if docindex:story.append(PageBreak())
        path=ROOT/name;lines=path.read_text().splitlines();i=0
        while i<len(lines):
            line=lines[i];i+=1
            if not line.strip():continue
            if line.startswith('!['):
                target=re.search(r'\]\(([^)]+)\)',line).group(1)
                imgpath=(path.parent/target).resolve()
                from reportlab.lib.utils import ImageReader
                iw,ih=ImageReader(str(imgpath)).getSize();width=465
                story.append(Image(str(imgpath),width=width,height=width*ih/iw));story.append(Spacer(1,10));continue
            if line.startswith('|'):
                rows=[line]
                while i<len(lines) and lines[i].startswith('|'):rows.append(lines[i]);i+=1
                cells=[[Paragraph(inline(c.strip()),styles['TableResearch']) for c in row.strip('|').split('|')]
                       for row in rows if not re.match(r'^\|[\s:|\-]+\|$',row)]
                if cells:
                    n=len(cells[0]); widths=[465/n]*n
                    if n==3:widths=[80,185,200]
                    t=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT')
                    t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8f2ef')),('LINEBELOW',(0,0),(-1,-1),.4,colors.HexColor('#d7e0df')),('BOTTOMPADDING',(0,0),(-1,-1),6)]));story.append(t);story.append(Spacer(1,10))
                continue
            if line.startswith('    '):
                code=[line[4:]]
                while i<len(lines) and lines[i].startswith('    '):code.append(lines[i][4:]);i+=1
                story.append(Preformatted('\n'.join(code),styles['CodeResearch'],maxLineLength=96));continue
            if line.startswith('#'):
                level=len(line)-len(line.lstrip('#'));style=styles['Title' if level==1 else 'Heading2'];story.append(Paragraph(inline(line.lstrip('#').strip()),style));continue
            story.append(Paragraph(inline(line),styles['BodyResearch']))
    def footer(canvas,doc):
        canvas.saveState();canvas.setFont('ResearchSans',7);canvas.setFillColor(colors.HexColor('#52626a'))
        canvas.drawString(48,28,'Purposeful activity | Research version 0.1 | Novelty unestablished')
        canvas.drawRightString(563,28,str(doc.page));canvas.restoreState()
    SimpleDocTemplate(str(ROOT/'results/research-note.pdf'),rightMargin=48,leftMargin=48,topMargin=44,bottomMargin=46,
        title='Physical sensing costs do not uniquely rank task survival',author='Mathematics End State research workspace').build(story,onFirstPage=footer,onLaterPages=footer)


def main():
    results=ROOT/'results';summary=json.loads((results/'summary.json').read_text())
    template=(ROOT/'docs/atlas-template.html').read_text()
    (results/'atlas.html').write_text(template.replace('__ATLAS_DATA__',(results/'atlas-data.json').read_text()))
    rows=[]
    for name,s in summary['sensors'].items():
        for law in ['hazard','damage']:
            mc=s[f'{law}_mc'];lo,hi=mc['ci95'];exact=s[f'{law}_survival']
            rows.append(f"| {name}, {law} | {exact:.9f} | {mc['estimate']:.9f} | [{lo:.6f}, {hi:.6f}] | {mc['samples']} | {mc['seed']} |")
    note='''# Verification record

The analytic proofs are in `docs/model.md`. Numerical checks support their implementation; they do not replace the proofs or establish literature novelty.

## Independent event-driven simulations

Matrix exponentials and Gillespie paths use different numerical methods. These are survival probabilities for the stated effective models, not empirical device measurements. Intervals are 95% Wilson binomial intervals; a valid implementation need not land inside every such interval.

| Design / law | Exact | Monte Carlo | 95% interval | Samples | Seed |
|---|---:|---:|---|---:|---:|
'''+ '\n'.join(rows)+f'''

## Energy and probability checks

- Maximum absolute disagreement among analytic power and independent four-state work, heat, and entropy-production calculations across 200 seeded parameter cases: {summary['max_flux_residual_200_cases']:.3e}.
- Unit tests exercise local detailed balance, stationary distributions, full-to-reduced survival agreement, the information–dissipation inequality, the kinetic-ceiling frontier, exact damage/hazard equivalence at one stage, and numerical simulation.
- Test command: `PYTHONPATH=src python -m unittest discover -s tests -v`.
- Simulation seeds, package versions, and source hashes are recorded in `summary.json` and `run-manifest.json`.

## Transfer and counterexamples

- Ranking reversal found in {summary['transfer']['reversals']} of {summary['transfer']['total']} explicitly sampled cases. This is a grid count, not an estimated frequency in nature.
- The sweep varies accuracy, horizon, number of damage stages, and initial sensor preparation. Every case is retained in `transfer.csv`.
- Single-stage damage cannot reverse the ordering relative to the corresponding instantaneous hazard because the models coincide.
- Changing the initial condition can change conclusions based on stationary survivor selection. Fresh damage initialization is an assumption in every damage run.

## Remaining uncertainty

No experimental validation; no exhaustive novelty review; no full energy accounting for the effective repair/actuator layers; no energy-reserve feedback; and no guaranteed global optimization outside the explicit two-state sensor family. The interactive atlas is a finite sample, not an impossibility certificate.
'''
    (results/'verification.md').write_text(note)
    manifest=[]
    for p in sorted((ROOT/'research/sources').glob('*.pdf')):
        manifest.append({'arxiv_id':p.stem,'url':'https://arxiv.org/pdf/'+p.stem,
            'retrieved_date':'2026-11-16','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
            'bytes':p.stat().st_size,'access':'downloaded PDF and locally extracted text; inspected relevant sections'})
    # Keep provenance when reproducing a release without the ignored PDF cache.
    if manifest:
        (ROOT/'research/source-manifest.json').write_text(json.dumps(manifest,indent=2))
    make_pdf()
    print('Built results/atlas.html, results/research-note.pdf, results/verification.md, source manifest.')


if __name__=='__main__':main()
