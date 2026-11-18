"""Render the proposed research agenda without rerunning numerical studies."""
from pathlib import Path
from build_research_docs import page

ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    source=(ROOT/'docs/future-work.md').read_text(encoding='utf-8')
    page('FUTURE-WORK.html',source,'Future work — self-maintaining organization')
    print('Built future-work HTML; all research directions remain proposals.')
