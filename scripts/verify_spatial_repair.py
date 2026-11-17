"""Check the finite-grid experiment against exact continuous-space path samples."""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from math import comb, log, pi, sqrt
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.special import ndtr, ndtri

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from purposeful.spatial_repair import SpatialRepair
from test_spatial_repair import radial_verification


def initial_positions(model, preparation, size, rng):
    """Draw exactly from the stated continuous relative-position preparation."""
    half_length = model.length / 2
    if preparation == "uniform":
        return rng.uniform(-half_length, half_length, size=(size, 3))
    if preparation == "localized":
        minimum = ndtr(-half_length / model.local_width)
        maximum = ndtr(half_length / model.local_width)
        return model.local_width * ndtri(rng.uniform(minimum, maximum, size=(size, 3)))
    if preparation == "reaction":
        directions = rng.standard_normal((size, 3))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        # If u=(r/a)^2, radial volume contributes u^(1/2) and the gate (1-u)^4.
        radii = model.radius * np.sqrt(rng.beta(1.5, 5.0, size=size))
        return directions * radii[:, None]
    raise ValueError(f"Unknown preparation {preparation}")


def exact_failure_times(model, preparation, initial, horizon, size, seed):
    """Sample first collapse times using labelled clocks and periodic Brownian motion.

    A chemical state fixes the dominating marked-Poisson rate until the next
    accepted jump. A candidate interval can therefore be sampled independently
    of the Brownian increment. Only candidate endpoints need to be evaluated.
    Paths that reach the finite horizon without collapse are right-censored.
    """
    rng = np.random.default_rng(seed)
    positions = initial_positions(model, preparation, size, rng)
    h, b = initial
    healthy_first = np.full(size, h >= 1)
    healthy_second = np.full(size, h == 2)
    fuel = np.full(size, b, dtype=np.int64)
    net_damage = np.zeros(size, dtype=np.int64)
    net_repair = np.zeros(size, dtype=np.int64)
    elapsed = np.zeros(size)
    collapse = np.full(size, np.inf)
    active = np.arange(size)
    candidates = 0
    accepted = 0
    half_length = model.length / 2
    while active.size:
        first, second = healthy_first[active], healthy_second[active]
        fuel_left = fuel[active]
        first_flip = model.eta * np.where(first, 1, 0.5)
        second_flip = model.eta * np.where(second, 1, 0.5)
        damage_total = first_flip + second_flip
        forward = model.activity * np.where(first ^ second, fuel_left, 0)
        reverse_one = model.activity * np.where(first & second, model.fuel - fuel_left, 0) / 4
        bound = damage_total + forward + 2 * reverse_one
        duration = rng.exponential(1 / bound)
        scheduled = elapsed[active] + duration <= horizon
        active = active[scheduled]
        if not active.size:
            break
        candidates += active.size
        duration = duration[scheduled]
        first, second = first[scheduled], second[scheduled]
        elapsed[active] += duration
        brownian_step = rng.standard_normal((active.size, 3))
        brownian_step *= np.sqrt(2 * model.diffusivity * duration)[:, None]
        positions[active] = (positions[active] + brownian_step + half_length) % model.length - half_length
        squared_radius = np.sum(np.square(positions[active]), axis=1)
        gate = np.maximum(0, 1 - squared_radius / model.radius**2) ** 4
        mark = bound[scheduled] * rng.random(active.size)
        forward_edge = damage_total[scheduled] + forward[scheduled] * gate
        first_reverse_edge = forward_edge + reverse_one[scheduled] * gate
        second_reverse_edge = first_reverse_edge + reverse_one[scheduled] * gate
        first_damage = mark < first_flip[scheduled]
        second_damage = (mark >= first_flip[scheduled]) & (mark < damage_total[scheduled])
        repair = (mark >= damage_total[scheduled]) & (mark < forward_edge)
        first_reverse = (mark >= forward_edge) & (mark < first_reverse_edge)
        second_reverse = (mark >= first_reverse_edge) & (mark < second_reverse_edge)
        net_damage[active[first_damage]] += np.where(first[first_damage], 1, -1)
        net_damage[active[second_damage]] += np.where(second[second_damage], 1, -1)
        net_repair[active[repair]] += 1
        net_repair[active[first_reverse | second_reverse]] -= 1
        healthy_first[active[first_damage]] ^= True
        healthy_second[active[second_damage]] ^= True
        healthy_first[active[repair]] = True
        healthy_second[active[repair]] = True
        healthy_first[active[first_reverse]] = False
        healthy_second[active[second_reverse]] = False
        fuel[active[repair]] -= 1
        fuel[active[first_reverse | second_reverse]] += 1
        accepted += np.count_nonzero(mark < second_reverse_edge)
        has_collapsed = ~(healthy_first[active] | healthy_second[active])
        collapse[active[has_collapsed]] = elapsed[active[has_collapsed]]
        active = active[~has_collapsed]
    if np.any((fuel < 0) | (fuel > model.fuel)):
        raise RuntimeError("Continuous path checker violated the conserved inventory")
    final_h = healthy_first.astype(np.int64) + healthy_second
    if np.any(final_h - h != -net_damage + net_repair) or np.any(fuel - b != -net_repair):
        raise RuntimeError("Continuous path checker violated the signed jump ledger")
    if np.any(net_damage[np.isinf(collapse)] > h + b - 1):
        raise RuntimeError("A live path exceeded the signed damage bound")
    energy_change = (fuel - b) * log(8) + (final_h - h) * log(2)
    energy_ledger = -net_damage * log(2) - net_repair * log(4)
    energy_residual = float(np.max(np.abs(energy_change - energy_ledger)))
    if energy_residual > 1e-12:
        raise RuntimeError("Continuous path checker violated the energy ledger")
    ledger_examples = [{"path_index": i, "final_h": int(final_h[i]), "final_b": int(fuel[i]),
                        "net_damage": int(net_damage[i]), "net_repair": int(net_repair[i]),
                        "dimensionless_energy_change": float(energy_change[i]),
                        "first_collapse_time": float(collapse[i]) if np.isfinite(collapse[i]) else None}
                       for i in range(min(size, 5))]
    return collapse, {"candidates": int(candidates), "accepted_jumps": int(accepted),
                      "maximum_energy_ledger_residual": energy_residual,
                      "ledger_examples": ledger_examples}


def wilson_interval(survival, size):
    z = 1.959963984540054
    denominator = 1 + z**2 / size
    center = (survival + z**2 / (2 * size)) / denominator
    radius = z / denominator * np.sqrt(survival * (1 - survival) / size + z**2 / (4 * size**2))
    return np.clip(center - radius, 0, 1), np.clip(center + radius, 0, 1)


def compare_paths(data, paths):
    model = SpatialRepair(**data["parameters"])
    times = np.asarray(data["times"])
    previous, finest = data["refinements"][-2:]
    cases = [
        ("HD", 1, "uniform", 2026111701),
        ("HD", 1, "localized", 2026111702),
        ("HD", 1, "reaction", 2026111703),
        ("HH", 2, "localized", 2026111704),
    ]
    alpha = 0.01
    dkw = sqrt(log(2 * len(cases) / alpha) / (2 * paths))
    comparisons = []
    for name, h, preparation, seed in cases:
        collapse, events = exact_failure_times(model, preparation, (h, model.fuel), times[-1], paths, seed)
        survivors = paths - np.searchsorted(np.sort(collapse), times, side="right")
        estimate = survivors / paths
        lower, upper = wilson_interval(estimate, paths)
        fine_curve = np.asarray(finest["initial_conditions"][name]["spatial"][preparation])
        previous_curve = np.asarray(previous["initial_conditions"][name]["spatial"][preparation])
        count_curve = np.asarray(finest["initial_conditions"][name]["count_edgewise"])
        mismatch = np.abs(estimate - fine_curve)
        j = int(mismatch.argmax())
        count_gap = estimate - count_curve
        largest_gap = int(np.abs(count_gap).argmax())
        signed_gap = float(count_gap[largest_gap])
        spatial_step = float(np.max(np.abs(fine_curve - previous_curve)))
        within = float(mismatch[j]) <= dkw
        observed_failures = collapse[np.isfinite(collapse)]
        failure_quartiles = np.quantile(observed_failures, [0.25, 0.5, 0.75]).tolist() if observed_failures.size else []
        comparisons.append({
            "chemistry": name,
            "preparation": preparation,
            "rng_seed": seed,
            "events": events,
            "failures_by_end_of_horizon": int(observed_failures.size),
            "observed_first_collapse_time_quartiles": failure_quartiles,
            "surviving_paths": survivors.tolist(),
            "monte_carlo": estimate.tolist(),
            "wilson_95_lower": lower.tolist(),
            "wilson_95_upper": upper.tolist(),
            "finest_grid": finest["n"],
            "maximum_dkw_check_discrepancy": float(mismatch[j]),
            "time_at_maximum_discrepancy": float(times[j]),
            "within_familywise_dkw_band": within,
            "finest_previous_grid_curve_difference": spatial_step,
            "comparison_to_fixed_finest_grid_calibrated_count": {
                "time_of_largest_observed_absolute_gap": float(times[largest_gap]),
                "observed_signed_gap": signed_gap,
                "simultaneous_signed_gap_interval_at_that_time": [
                    max(-float(count_curve[largest_gap]), signed_gap - dkw),
                    min(1 - float(count_curve[largest_gap]), signed_gap + dkw),
                ],
                "simultaneous_lower_bound_on_continuous_time_maximum_absolute_gap": max(0, abs(signed_gap) - dkw),
            },
        })
        print(f"{name} {preparation}: MC difference {mismatch[j]:.6f}, DKW half-width "
              f"{dkw:.6f}, grid difference {spatial_step:.6g}", flush=True)
    return {"confidence": 1 - alpha, "paths_per_case": paths, "dkw_half_width": dkw,
            "time_mesh": times.tolist(), "cases": comparisons}


def verification_figure(data, comparison):
    background, ink, muted = "#0d141d", "#e6edf7", "#adc0d6"
    colors = ("#72dbc1", "#d899e8", "#edbd75", "#91b5ff")
    labels = ("HD uniform", "HD localized", "HD reverse-repair Palm", "HH localized")
    times = np.asarray(comparison["time_mesh"])
    step = max(1, (len(times) - 1) // 14)
    marker_sites = list(range(0, len(times), step))
    finest = data["refinements"][-1]
    with plt.rc_context({"figure.facecolor": background, "axes.facecolor": background,
                         "text.color": ink, "axes.labelcolor": ink, "axes.edgecolor": muted,
                         "xtick.color": muted, "ytick.color": muted,
                         "axes.spines.top": False, "axes.spines.right": False}):
        fig, (survival, residual) = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
        dkw = comparison["dkw_half_width"]
        residual_limit = 1.22 * max(dkw, *(case["maximum_dkw_check_discrepancy"] for case in comparison["cases"]))
        residual.axhspan(-dkw, dkw, color=muted, alpha=0.15, label="99% simultaneous sampling band")
        residual.axhline(0, color=muted, linewidth=0.8)
        for i, (case, label, color) in enumerate(zip(comparison["cases"], labels, colors)):
            fine = np.asarray(finest["initial_conditions"][case["chemistry"]]["spatial"][case["preparation"]])
            estimate = np.asarray(case["monte_carlo"])
            survival.plot(times, fine, color=color, label=label)
            indices = marker_sites[i::2] if i < 2 else marker_sites[(i - 2)::2]
            survival.plot(times[indices], estimate[indices], "o", color=color, markersize=4,
                          markeredgewidth=0.5, markeredgecolor=background)
            residual.plot(times, estimate - fine, color=color, linewidth=1.0)
        survival.set(title=f"Continuum samples and {finest['n']}³ lattice solution",
                     xlabel="Time (model units)", ylabel="Probability of no first collapse", ylim=(0, 1.02))
        survival.legend(title=f"Lattice curves · {comparison['paths_per_case']:,} paths per preparation",
                        facecolor=background, labelcolor=ink, fontsize=8, title_fontsize=8)
        residual.set(title="Difference from the lattice solution", xlabel="Time (model units)",
                     ylabel="Monte Carlo minus lattice survival", ylim=(-residual_limit, residual_limit))
        residual.legend(facecolor=background, labelcolor=ink, fontsize=8)
        paths = [ROOT / f"results/figures/spatial-repair-independent.{extension}" for extension in ("png", "pdf")]
        for path in paths:
            fig.savefig(path, dpi=180)
        plt.close(fig)
    return {path.relative_to(ROOT).as_posix(): file_hash(path) for path in paths}


def check_deterministic(data):
    probabilities = []
    for grid in data["refinements"]:
        for conditions in grid["initial_conditions"].values():
            probabilities.extend(conditions["spatial"].values())
            probabilities.extend((conditions["count_edgewise"], conditions["count_reference"]))
    curves = np.asarray(probabilities)
    if curves.ndim != 2 or curves.shape[1] != len(data["times"]) or not np.all(np.isfinite(curves)):
        raise ValueError("Every deterministic curve must have a finite value at every time")
    minimum = float(curves.min())
    maximum = float(curves.max())
    max_initial_error = float(np.max(np.abs(curves[:, 0] - 1)))
    max_increase = max(0., float(np.diff(curves, axis=1).max()))
    probability_tolerance = 1e-10
    previous, finest = data["refinements"][-2:]
    refinement_differences = {}
    for chemistry, conditions in finest["initial_conditions"].items():
        earlier_curves = previous["initial_conditions"][chemistry]["spatial"]
        refinement_differences[chemistry] = {
            preparation: float(np.max(np.abs(np.asarray(curve) - earlier_curves[preparation])))
            for preparation, curve in conditions["spatial"].items()
        }
    maximum_refinement_gap = max(gap for conditions in refinement_differences.values() for gap in conditions.values())
    maximum_calibration_residual = max(
        (float(value) for grid in data["refinements"] for value in grid["calibration_diagnostics"]["infinity_residual"]),
        default=0.,
    )
    maximum_polynomial_tail = max(
        grid["exponential_solver_diagnostics"]["exact_arithmetic_observable_tail_upper"]
        for grid in data["refinements"]
    )
    scale = data["parameters"]["radius"] / data["parameters"]["length"]
    continuum_gate = 4 * pi * scale**3 * sum((-1)**j * comb(4, j) / (2 * j + 3) for j in range(5))
    refinement_tolerance, residual_tolerance, tail_tolerance = 1e-4, 1e-8, 2e-14
    return {
        "probability_curves": {"checked": len(probabilities), "minimum": minimum,
                               "maximum": maximum, "maximum_initial_error": max_initial_error,
                               "maximum_increase": max_increase, "tolerance": probability_tolerance,
                               "passed": min(minimum, 1 - maximum) >= -probability_tolerance
                               and max(max_initial_error, max_increase) <= probability_tolerance},
        "last_grid_refinement": {"grids": [previous["n"], finest["n"]], "curves": refinement_differences,
                                 "maximum_difference": maximum_refinement_gap, "tolerance": refinement_tolerance,
                                 "passed": maximum_refinement_gap < refinement_tolerance},
        "forward_mean_linear_residual": {"maximum": maximum_calibration_residual, "tolerance": residual_tolerance,
                                         "passed": maximum_calibration_residual < residual_tolerance},
        "observable_polynomial_tail_exact_arithmetic": {"maximum": maximum_polynomial_tail,
                                                        "tolerance": tail_tolerance,
                                                        "passed": maximum_polynomial_tail <= tail_tolerance},
        "analytic_uniform_gate_mean": continuum_gate,
        "finest_grid_gate_mean_absolute_error": abs(finest["mean_gate"] - continuum_gate),
    }


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="results/spatial-repair.json")
    parser.add_argument("--output", default="results/spatial-repair-verification.json")
    parser.add_argument("--paths", type=int, default=50_000)
    args = parser.parse_args()
    if args.paths <= 0:
        parser.error("Supply a positive number of paths per preparation.")
    data_path = ROOT / args.input
    raw = data_path.read_bytes()
    data = json.loads(raw)
    input_hash = sha256(raw).hexdigest()
    mismatches = [name for name, digest in data["source_sha256"].items() if file_hash(ROOT / name) != digest]
    if mismatches:
        raise SystemExit(f"Rebuild spatial-repair data with the current sources: {mismatches}")
    sources = sorted([*ROOT.glob("src/purposeful/*.py"), *ROOT.glob("tests/test_*.py"),
                      ROOT / "scripts/verify_spatial_repair.py", ROOT / "scripts/run_spatial_repair.py",
                      ROOT / "research/spatial-repair-verification-design.md", ROOT / "pyproject.toml"])
    source_hashes = {source.relative_to(ROOT).as_posix(): file_hash(source) for source in sources}
    test_command = [sys.executable, "-B", "-m", "unittest", "discover", "-s", "tests", "-v"]
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    print("Running unit tests, including independently assembled dense and radial references", flush=True)
    tests = subprocess.run(test_command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, encoding="utf-8", errors="replace", check=False)
    test_count = re.search(r"Ran (\d+) tests?", tests.stdout)
    print(tests.stdout, end="", flush=True)
    deterministic = check_deterministic(data)
    deterministic_passed = all(item["passed"] for item in deterministic.values() if isinstance(item, dict))
    comparison = compare_paths(data, args.paths)
    paths_passed = all(case["within_familywise_dkw_band"] for case in comparison["cases"])
    passed = tests.returncode == 0 and deterministic_passed and paths_passed
    if file_hash(data_path) != input_hash:
        raise SystemExit("The experiment data changed during verification. Replay against a stable receipt.")
    if any(file_hash(ROOT / name) != digest for name, digest in source_hashes.items()):
        raise SystemExit("Sources changed during verification. Replay against stable sources.")
    case_results = comparison["cases"]
    previous, finest = data["refinements"][-2:]
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Finite-grid operator checks, independent spherical benchmark, and exact continuum path sampling.",
        "passed": passed,
        "summary": {"status": "passed" if passed else "failed",
                    "test_count": int(test_count[1]) if test_count else None,
                    "structural_checks": re.findall(r"^(test_\w+) \(test_spatial_repair\..*?\) \.\.\. ok$",
                                                     tests.stdout, flags=re.MULTILINE),
                    "monte_carlo_paths_per_preparation": args.paths,
                    "monte_carlo_paths_total": args.paths * len(case_results),
                    "confidence": comparison["confidence"],
                    "dkw_half_width": comparison["dkw_half_width"],
                    "maximum_mc_curve_discrepancy": max(x["maximum_dkw_check_discrepancy"] for x in case_results),
                    "grids_compared_for_refinement": [previous["n"], finest["n"]],
                    "maximum_last_grid_difference": deterministic["last_grid_refinement"]["maximum_difference"],
                    "deterministic_checks_passed": deterministic_passed},
        "input_sha256": input_hash,
        "source_sha256": source_hashes,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "tests": {"command": "python -B -m unittest discover -s tests -v",
                  "return_code": tests.returncode, "count": int(test_count[1]) if test_count else None,
                  "output": tests.stdout},
        "radial_benchmark": radial_verification(),
        "deterministic_checks": deterministic,
        "continuous_path_comparison": comparison,
        "figure_sha256": verification_figure(data, comparison),
        "uncertainty": "Wilson intervals are pointwise. The DKW envelope bounds sampling error simultaneously "
                       "over time and over the four preparations. Reported successive-grid differences are "
                       "observed refinement diagnostics. The count comparison treats the finest-grid calibrated curve "
                       "as fixed. It includes no enclosure for the limiting count calibration coefficients.",
    }
    (ROOT / args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Verification {'passed' if passed else 'failed'}. Receipt at {args.output}.", flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
