"""Build the numerical record and figures for the shared-reserve study."""
from dataclasses import replace
from itertools import product
import json
from math import log, sqrt
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from purposeful.closed_loop import (
    Device,
    STATIONARY_PRIOR,
    best_fixed_frontier,
    block_failure_bound,
    evaluate_fixed_schedule,
    minimum_reserve,
    optimal_policy,
)
from purposeful.reserve import RepairContact


BACKGROUND = "#101824"
PANEL = "#172437"
INK = "#edf3fc"
MUTED = "#b5c4d9"
BLUE = "#61c5f7"
AMBER = "#ffbe71"
GREEN = "#6de1ad"


def save(fig, name):
    fig.savefig(ROOT / "results" / "figures" / f"closed-loop-{name}.png", dpi=165)
    plt.close(fig)


def node(ax, x, y, label, color=BLUE, align="center"):
    ax.text(x, y, label, ha=align, va="center", color=INK, fontsize=11, linespacing=1.55,
            zorder=3, bbox={"boxstyle": "round,pad=0.55", "fc": PANEL, "ec": color})


def arrow(ax, start, end, color=MUTED, connection="arc3"):
    ax.annotate("", end, start, arrowprops={"arrowstyle": "->", "lw": 1.5, "color": color,
                                          "connectionstyle": connection})


def protocol_figure():
    fig, ax = plt.subplots(figsize=(11.7, 3.2), layout="constrained")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.set_title("One battery funds sensing, actuation, and repair", loc="left", pad=14, weight="bold")
    stages = [
        (0.13, "Sense or skip\nProbe costs c units"),
        (0.37, "Idle or engage\nEngagement costs d units"),
        (0.63, "Harvest or incur damage\nEngagement can refill"),
        (0.88, "Repair if damaged\nNet repair costs one unit"),
    ]
    for x, label in stages:
        node(ax, x, 0.61, label, AMBER if x == 0.88 else BLUE)
    for start, end in [(0.215, 0.275), (0.47, 0.525), (0.738, 0.785)]:
        arrow(ax, (start, 0.61), (end, 0.61))
    arrow(ax, (0.87, 0.41), (0.13, 0.41), GREEN, "arc3,rad=-0.13")
    ax.text(0.5, 0.085, "Survive only if healthy at the checkpoint. Carry the remaining reserve into the next round.",
            color=GREEN, ha="center", fontsize=11)
    save(fig, "protocol")


def witness_figure(fixed, adaptive):
    fig, (tree, bars) = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained",
                                    gridspec_kw={"width_ratios": [1.15, 1]})
    tree.set_title("Optimal feedback after round 1", loc="left", weight="bold")
    tree.set(xlim=(0, 1), ylim=(0, 1))
    tree.axis("off")
    node(tree, 0.5, 0.81, "Start with two units\nPay for the first probe")
    node(tree, 0.23, 0.24, "One unit left\nSkip probe · 75% survival", AMBER)
    node(tree, 0.77, 0.24, "Two units left\nProbe · 100% survival", GREEN)
    arrow(tree, (0.45, 0.7), (0.24, 0.36), AMBER)
    arrow(tree, (0.55, 0.7), (0.76, 0.36), GREEN)
    tree.text(0.27, 0.54, "Idle on X=0", ha="right", color=MUTED, fontsize=10)
    tree.text(0.72, 0.54, "Engage on X=1", ha="left", color=MUTED, fontsize=10)
    tree.text(0.5, 0.04, "Each branch occurs with probability 1/2", ha="center", color=MUTED)

    names = ["Skip both", "Skip, probe", "Probe, skip", "Probe both", "Adaptive"]
    values = [100 * probability for probability in (*fixed, adaptive)]
    x = np.arange(len(names))
    bars.bar(x, values, width=0.68, color=[PANEL, PANEL, PANEL, PANEL, GREEN], edgecolor=GREEN)
    bars.set(xticks=x, xticklabels=names, ylim=(0, 100), ylabel="Two-round survival (%)")
    bars.tick_params(axis="x", rotation=29)
    bars.set_title("Every fixed probe schedule", loc="left", weight="bold")
    for position, value in zip(x, values):
        bars.text(position, value + 2.2, f"{value:g}", ha="center", color=INK, fontsize=9)
    save(fig, "witness")


def robustness_figure(fidelity, persistence, gap, rounds, adaptive, fixed, no_harvest):
    fig, (heat, decay) = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    cmap = LinearSegmentedColormap.from_list("gap", [PANEL, "#2666a0", GREEN])
    image = heat.pcolormesh(fidelity, persistence, 100 * gap, cmap=cmap, vmin=0, vmax=12.5,
                           shading="nearest")
    heat.axvline(2 / 3, ymax=0.16, color=AMBER, ls="--", lw=1.2)
    heat.annotate("Exact IID threshold", (2 / 3, 0.025), (0.53, 0.23), color=AMBER,
                  arrowprops={"arrowstyle": "->", "color": AMBER}, fontsize=9)
    heat.set(xlabel="Sensing accuracy", ylabel="Environment persistence")
    heat.set_title("Two-round advantage with noise", loc="left", weight="bold")
    fig.colorbar(image, ax=heat, label="Adaptive gain (percentage points)")

    decay.semilogy(rounds, adaptive, color=GREEN, lw=2.4, label="Optimal adaptive")
    decay.semilogy(rounds, fixed, color=BLUE, lw=2.4, label="Best fixed schedule")
    decay.semilogy(rounds, no_harvest, color=MUTED, lw=1.5, label="Optimal with no harvest")
    decay.semilogy(rounds, ((1 + sqrt(5)) / 4) ** rounds, color=GREEN, ls="--", alpha=0.75,
                   label="Proven adaptive rate")
    decay.semilogy(rounds, (1 / sqrt(2)) ** rounds, color=BLUE, ls="--", alpha=0.75,
                   label="Proven fixed rate")
    decay.set(xlabel="Rounds", ylabel="Survival probability (log scale)")
    decay.set_title("Different exponential decay rates", loc="left", weight="bold")
    decay.legend(fontsize=8, frameon=False, loc="lower left")
    save(fig, "robustness")


def repair_rate_figure():
    kappa = np.linspace(0, 1, 201)
    adaptive = (1 + sqrt(5)) / 4
    alternating = 1 / sqrt(2)
    blind = (1 + np.sqrt(1 + 4 * kappa)) / 4
    threshold = 2 - sqrt(2)
    fixed = np.maximum(alternating, blind)
    fig, ax = plt.subplots(figsize=(8.5, 3.4), layout="constrained")
    ax.plot(kappa, fixed, color=BLUE, lw=2.5, label="Best fixed rate")
    ax.axhline(adaptive, color=GREEN, lw=2.5, label="Adaptive rate")
    ax.fill_between(kappa, fixed, adaptive, color=GREEN, alpha=0.12)
    ax.axvline(threshold, color=AMBER, ls="--", lw=1.2)
    ax.text(threshold - 0.025, 0.686, "Alternating", ha="right", color=MUTED)
    ax.text(threshold + 0.025, 0.686, "Always blind", ha="left", color=MUTED)
    ax.annotate("Switch at 2 − √2", xy=(threshold, alternating), xytext=(0.29, 0.753),
                color=AMBER, arrowprops={"arrowstyle": "->", "color": AMBER})
    ax.set(xlim=(0, 1), ylim=(0.665, 0.835), xlabel="Repair success probability per contact",
           ylabel="Asymptotic survival factor")
    ax.set_title("The best fixed schedule changes with repair reliability", loc="left", weight="bold")
    ax.legend(fontsize=9, frameon=False, loc="lower left", bbox_to_anchor=(0.01, 0.23))
    save(fig, "repair-rate")


def main():
    plt.rcParams.update({
        "figure.facecolor": BACKGROUND, "axes.facecolor": BACKGROUND, "text.color": INK,
        "axes.labelcolor": INK, "axes.titlecolor": INK, "axes.edgecolor": MUTED,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
        "axes.spines.right": False, "font.size": 10,
    })
    base = Device(2, 1, 1, 2, 0.5, RepairContact(log(3), log(3), 1.0))
    state = STATIONARY_PRIOR, base.capacity
    schedules = list(product([False, True], repeat=2))
    two_round = optimal_policy(base, 2)
    fixed_witness = [float(evaluate_fixed_schedule(base, word)[state]) for word in schedules]
    adaptive_witness = float(two_round.survival[2][state])
    protocol_figure()
    witness_figure(fixed_witness, adaptive_witness)
    fidelity = np.linspace(0.5, 1, 51)
    persistence = np.linspace(0, 1, 51)
    gap = np.zeros((len(persistence), len(fidelity)))
    for i, rho in enumerate(persistence):
        for j, accuracy in enumerate(fidelity):
            device = replace(base, fidelity=accuracy, persistence=rho)
            adaptive = optimal_policy(device, 2).survival[2][state]
            fixed = best_fixed_frontier(device, 2)[0][state]
            gap[i, j] = max(0.0, adaptive - fixed)
    horizons = np.arange(41)
    adaptive_series = optimal_policy(base, int(horizons[-1])).survival[:, state[0], state[1]]
    fixed_series = np.array([best_fixed_frontier(base, int(n))[0][state] for n in horizons])
    no_harvest = optimal_policy(replace(base, harvest=0), int(horizons[-1]))
    no_harvest_series = no_harvest.survival[:, state[0], state[1]]
    robustness_figure(fidelity, persistence, gap, horizons, adaptive_series, fixed_series, no_harvest_series)
    repair_rate_figure()
    results = {
        "witness": {
            "capacity": base.capacity,
            "sensing_cost_quanta": base.sensing_cost,
            "actuation_cost_quanta": base.actuation_cost,
            "harvest_quanta": base.harvest,
            "quantum_kBT": base.repair.token_energy,
            "opportunity_probability": base.opportunity,
            "repair_success": base.repair.success,
            "adaptive_survival": adaptive_witness,
            "best_fixed_survival": max(fixed_witness),
            "fixed_schedules": {"".join(str(int(bit)) for bit in word): value
                                for word, value in zip(schedules, fixed_witness)},
            "minimum_reserve_for_eighty_percent": minimum_reserve(two_round, 0.8),
        },
        "noise": {
            "iid_strict_gain_above_accuracy": 2 / 3,
            "iid_probe_first_above_accuracy": 1 / sqrt(2),
            "accuracy_grid": fidelity.tolist(),
            "persistence_grid": persistence.tolist(),
            "survival_gain": gap.tolist(),
        },
        "long_horizon": {
            "rounds": horizons.tolist(),
            "adaptive_survival": adaptive_series.tolist(),
            "fixed_survival": fixed_series.tolist(),
            "no_harvest_adaptive_survival": no_harvest_series.tolist(),
            "adaptive_exact_root_rate": (1 + sqrt(5)) / 4,
            "fixed_exact_root_rate": 1 / sqrt(2),
            "fixed_schedule_switch_repair_probability": 2 - sqrt(2),
            "no_harvest_exact_root_rate": 0.5,
            "uniform_block_failure_bound": block_failure_bound(base),
        },
    }
    (ROOT / "results" / "closed-loop.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(results["witness"], indent=2))


if __name__ == "__main__":
    main()
