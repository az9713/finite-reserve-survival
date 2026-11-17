"""Reproduce the certified comparisons and figures for the autonomous network."""
from fractions import Fraction
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from purposeful.autonomous import (
    Network, benchmark_kinetics, spectral_interval, survival_curves,
    survival_intervals, uniform_activity_enclosure,
)


BG, PANEL, INK, MUTED = "#101824", "#172437", "#edf3fc", "#b5c4d9"
BLUE, GREEN, AMBER = "#61c5f7", "#6de1ad", "#ffbe71"


def save(fig, name):
    fig.savefig(ROOT / "results" / "figures" / ("autonomous-" + name + ".png"), dpi=165)
    plt.close(fig)


def diagram():
    fig, ax = plt.subplots(figsize=(9.2, 4.4), layout="constrained")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.set_title("Autonomous reactions share a finite battery", loc="left", weight="bold", pad=15)

    def box(x, y, label, color):
        ax.text(x, y, label, color=INK, ha="center", va="center", fontsize=10, linespacing=1.55,
                zorder=3, bbox={"boxstyle": "round,pad=.65", "facecolor": PANEL, "edgecolor": color})

    def link(src, dst, color, kind):
        ax.annotate("", dst, src, arrowprops={"arrowstyle": kind, "color": color, "lw": 1.6})

    box(.5, .52, "Battery\n0, 1 or 2 units", GREEN)
    box(.15, .76, "Sensing\nUpdates memory Y", BLUE)
    box(.86, .76, "Control\nChanges actuator A", BLUE)
    box(.15, .26, "Harvest\nFuel enters when X=A=1", GREEN)
    box(.86, .26, "Repair\nRemoves one damage step", AMBER)
    link((.428, .555), (.264, .72), BLUE, "<->")
    link((.572, .555), (.752, .72), BLUE, "<->")
    link((.423, .47), (.277, .31), GREEN, "<->")
    link((.577, .47), (.734, .31), AMBER, "<->")
    ax.text(.5, .865, "Environment X flips continuously", ha="center", color=MUTED, fontsize=10)
    ax.text(.5, .09, "Mismatch advances wear. First entry into damage level 2 ends survival.",
            ha="center", color=AMBER, fontsize=10)
    save(fig, "reactions")


def plots(times, adaptive, fixed, uncoupled, mask_bounds, analogue, exponents):
    fig, (curve, bars) = plt.subplots(1, 2, figsize=(11.5, 4.1), layout="constrained")
    curve.plot(times, adaptive, lw=2.4, color=GREEN, label="Best reserve-dependent circuit")
    curve.plot(times, fixed, lw=2.0, color=BLUE, label="Best binary uniform circuit")
    curve.plot(times, uncoupled, lw=1.5, color=MUTED, ls="--", label="Sensor disconnected from X")
    curve.scatter([4], [analogue], marker="v", color=AMBER, s=42,
                  label="Upper bound for any constant intensity")
    curve.set(xlabel="Time", ylabel="Survival probability", xlim=(0, 12), ylim=(0, 1.03))
    curve.set_title("First-passage survival", loc="left", weight="bold")
    curve.legend(loc="upper right", frameon=False, fontsize=8)

    labels = ["Off", "Full charge\nonly", "One unit\nonly", "Always on"]
    colors = [BLUE, GREEN, AMBER, BLUE]
    bars.bar(labels, [bounds[1] * 100 for bounds in mask_bounds], color=colors, width=.65)
    bars.set(ylim=(0, 62), ylabel="Best survival at time 4 (%)")
    bars.set_title("All 324 binary circuits", loc="left", weight="bold")
    for i, bounds in enumerate(mask_bounds):
        bars.text(i, bounds[1] * 100 + 1.1, f"{bounds[1] * 100:.2f}", ha="center", fontsize=9, color=INK)
    save(fig, "survival")

    fig, ax = plt.subplots(figsize=(8.4, 3.6), layout="constrained")
    values = [-bounds[0] for bounds in exponents]
    ax.barh(labels, values, color=colors, height=.58)
    ax.invert_yaxis()
    ax.set(xlim=(0, .31), xlabel="Long-run decay exponent (lower is better)")
    for i, value in enumerate(values):
        ax.text(value + .004, i, f"{value:.5f}", va="center", fontsize=9, color=INK)
    ax.set_title("Turning sensing off gives the best long-run rate", loc="left", weight="bold")
    save(fig, "decay")


def main():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "text.color": INK,
        "axes.labelcolor": INK, "axes.titlecolor": INK, "axes.edgecolor": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "font.size": 10,
    })
    kinetics = benchmark_kinetics()
    network, disconnected = Network(kinetics), Network(kinetics, informative=False)
    circuits, matrices = network.digital_circuits()
    _, blank_matrices = disconnected.digital_circuits()
    initial = network.initial(2)
    time = 4
    lo, hi, scale = survival_intervals(matrices, network.denominator, time, initial)
    _, blank_hi, blank_scale = survival_intervals(
        blank_matrices, disconnected.denominator, time, initial
    )
    winner = int(lo.argmax())
    runner_up = int(np.delete(hi, winner).max())
    assert lo[winner] > runner_up
    reserve_frontier = []
    for battery in range(kinetics.capacity + 1):
        if battery == kinetics.capacity:
            low, high, divisor = lo, hi, scale
        else:
            low, high, divisor = survival_intervals(matrices, network.denominator, time, network.initial(battery))
        reserve_frontier.append([int(low.max()), int(high.max()), divisor])
    target = Fraction(13, 25)
    assert all(Fraction(high, divisor) < target for _, high, divisor in reserve_frontier[:-1])
    assert Fraction(reserve_frontier[-1][0], reserve_frontier[-1][2]) > target
    unique_masks = list(dict.fromkeys(circuit.sensing for circuit in circuits))
    groups = [np.array([c.sensing == mask for c in circuits]) for mask in unique_masks]
    finite = [(int(lo[m].max()) / scale, int(hi[m].max()) / scale) for m in groups]
    analogue = uniform_activity_enclosure(network, time, initial)
    assert Fraction(int(lo[winner]), scale) > analogue[1]
    assert Fraction(int(lo[winner]), scale) > Fraction(int(blank_hi.max()), blank_scale)
    spectral = [spectral_interval(matrix, network.denominator, initial) for matrix in matrices]
    by_mask = [
        (max(lo for (lo, _), keep in zip(spectral, group) if keep),
         max(hi for (_, hi), keep in zip(spectral, group) if keep))
        for group in groups
    ]
    assert by_mask[0][0] > max(pair[1] for pair in by_mask[1:])
    late_winner = max(range(len(circuits)), key=lambda i: spectral[i][0])
    print("Certified finite-time winner", circuits[winner], flush=True)
    times = np.linspace(0, 12, 121)
    curves = survival_curves(matrices, network.denominator, times, initial)
    digital_fixed = curves[groups[0] | groups[-1]].max(axis=0)
    blank_curves = survival_curves(blank_matrices, disconnected.denominator, times, initial)
    adaptive = curves.max(axis=0)
    blank = blank_curves.max(axis=0)
    diagram()
    plots(times, adaptive, digital_fixed, blank, finite, float(analogue[1]), by_mask)
    result = {
        "time": time,
        "binary_circuits": len(circuits),
        "design": {key: str(value) for key, value in vars(kinetics).items()},
        "winning_circuit": vars(circuits[winner]),
        "winner_survival_integer_enclosure": [int(lo[winner]), int(hi[winner]), scale],
        "runner_up_survival_upper_numerator": runner_up,
        "reserve_frontier_integer_enclosures": reserve_frontier,
        "minimum_initial_units_for_52_percent": kinetics.capacity,
        "sensor_disconnected_upper": float(blank_hi.max() / blank_scale),
        "analog_uniform_enclosure": [str(bound) for bound in analogue],
        "mask_survival_enclosures": {str(mask): interval for mask, interval in zip(unique_masks, finite)},
        "mask_spectral_enclosures": {str(mask): [str(v) for v in interval]
                                     for mask, interval in zip(unique_masks, by_mask)},
        "long_run_circuit": vars(circuits[late_winner]),
        "curves": {"time": times.tolist(), "adaptive": adaptive.tolist(),
                   "digital_uniform": digital_fixed.tolist(), "no_sensor_input": blank.tolist()},
    }
    (ROOT / "results" / "autonomous.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(f"Certified winner survival {lo[winner] / scale:.8f}–{hi[winner] / scale:.8f}")
    print(f"Best uniform intensity at most {float(analogue[1]):.8f}")
    print(f"Best sensor-disconnected circuit at most {blank_hi.max() / blank_scale:.8f}")
    print("Optimal asymptotic spectral abscissa", *map(float, by_mask[0]))


if __name__ == "__main__":
    main()
