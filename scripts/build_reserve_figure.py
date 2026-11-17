"""Generate Study 02 evidence and its two-panel explanatory figure."""
from pathlib import Path
import json
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]


def main():
    run = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/run_reserve_study.py')],
                         capture_output=True, text=True, check=True)
    data = json.loads(run.stdout)
    (ROOT / 'results/reserve.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    bg, ink, muted = '#0d141d', '#e6edf7', '#adc0d6'
    plt.rcParams.update({'figure.facecolor': bg, 'axes.facecolor': bg, 'text.color': ink,
                         'axes.labelcolor': ink, 'axes.edgecolor': muted,
                         'xtick.color': muted, 'ytick.color': muted, 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    for label, color, offset in [('slow', '#72dbc1', -.18), ('fast', '#edbd75', .18)]:
        values = np.array(data['designs'][label]['maximum_survival']) * 100
        bars = axes[0].bar(np.arange(3) + offset, values, width=.34, label=label.title() + ' stress', color=color)
        axes[0].bar_label(bars, fmt='%.2f', fontsize=9, color=ink, padding=3)
    axes[0].set(xticks=[0, 1, 2], xlabel='Initial reserve units', ylabel='Optimal two-round survival (%)',
                ylim=(0, 105), title='The preferred error correlation changes with reserve')
    axes[0].axhline(88, color=muted, linestyle=':', linewidth=1)
    axes[0].legend(loc='lower right', facecolor=bg, edgecolor=muted, labelcolor=ink)
    energy = np.linspace(0, 12, 401)
    survival = ((1-np.exp(-2))*expit(energy/3))**3
    minimum = data['three_certain_stresses']['minimum_initial_reserve_kBT']
    axes[1].plot(energy, survival, color='#72dbc1', linewidth=2.3)
    axes[1].axhline(.4, color=muted, linestyle=':', linewidth=1)
    axes[1].axvline(minimum, color='#edbd75', linestyle='--', linewidth=1)
    axes[1].scatter([minimum], [.4], color='#edbd75', zorder=3)
    axes[1].text(minimum+.3, .2, f'Minimum energy\n{minimum:.4f} kBT', color='#edbd75')
    axes[1].set(xlabel='Initial reserve energy F / (kBT)', ylabel='Conditional survival probability',
                ylim=(0, .7), title='Three certain repairs: the equal-gap bound is attained')
    fig.savefig(ROOT / 'results/figures/reserve-frontier.png', dpi=180)
    fig.savefig(ROOT / 'results/figures/reserve-frontier.pdf')
    plt.close(fig)
    print('Saved exact Study 02 data and reserve-frontier figure.')


if __name__ == '__main__':
    main()
