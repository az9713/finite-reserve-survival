"""Render the explicit stochastic-calculus chapter with existing offline assets."""
from pathlib import Path
from build_research_docs import page

ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    page('STOCHASTIC-CALCULUS.html',(ROOT/'docs/stochastic-calculus-analysis.md').read_text(encoding='utf-8'),
         'Stochastic calculus analysis — Brownian motion, reaction martingales and survival')
