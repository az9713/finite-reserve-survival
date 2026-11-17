"""Non-validated Fourier Galerkin check of the Brownian organization model.

Projects the exact continuum differential operator in flat L2 after Gibbs
conjugation. Coefficient integration, floating arithmetic, and omission of high
Fourier modes are numerical approximations. This script does not bound them.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import scipy
from scipy.fft import dctn, idctn
from scipy.integrate import quad_vec
from scipy.special import ive

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from purposeful.organization import Organization
from purposeful.spatial_repair import _chebyshev_degree, chemical_generators, equilibrium

SPATIAL_AXES = (-3, -2, -1)


def radial_coefficients(model, maximum_mode, tolerance):
    """Return cube Fourier coefficients indexed by squared frequency length."""
    indices = np.arange(maximum_mode + 1)
    squared = indices[:, None, None] ** 2 + indices[None, :, None] ** 2 + indices[None, None, :] ** 2
    frequencies = 2 * np.pi * np.sqrt(np.arange(3 * maximum_mode**2 + 1))
    radius = model.radius
    offset_weight = np.exp(-model.offset / 2)

    def integrand(r):
        core = 1 - (r / radius) ** 2
        energy = model.offset - model.depth * core**2
        slope = 4 * model.depth * r * core / radius**2
        laplacian = 4 * model.depth * (3 - 5 * (r / radius) ** 2) / radius**2
        gate = core**4
        factors = np.array([
            gate, gate * np.exp(-max(0.0, energy)), gate * np.exp(min(0.0, energy)),
            gate * np.exp(-abs(energy) / 2),
            model.diffusivity * (laplacian / 2 - slope**2 / 4),
            np.exp(-energy / 2) - offset_weight,
        ])
        return (4 * np.pi * r**2) * factors[:, None] * np.sinc(frequencies * r / np.pi)

    kink = radius * np.sqrt(1 - np.sqrt(model.offset / model.depth))
    values, error_estimate = quad_vec(
        integrand, 0, radius, epsrel=tolerance, epsabs=tolerance / 20,
        points=[kink], limit=600,
    )
    values[-1, 0] += offset_weight
    return values[:, squared], float(error_estimate)


def solve(model, cutoff, policies, tolerance):
    start = perf_counter()
    states, damage, repair = chemical_generators(model)
    live = [i for i, (healthy, _) in enumerate(states) if healthy]
    states = [states[i] for i in live]
    locations = {s: i for i, s in enumerate(states)}
    launch = locations[1, model.fuel]
    weights = equilibrium(model)[live]
    weights = np.sqrt(weights / weights[launch])
    damage = weights[:, None] * damage[np.ix_(live, live)] / weights[None, :]
    repair = weights[:, None] * repair[np.ix_(live, live)] / weights[None, :]

    # A 2M-periodic grid with M=2K+1 aliases frequencies at 4K+2.
    # The multiplication V_{2K} u_K reaches 3K. Aliases cannot reach <=K.
    padded_intervals = 2 * cutoff + 1
    coefficients, integration_error = radial_coefficients(model, 2 * cutoff, tolerance)
    padded_shape = (padded_intervals + 1,) * 3
    padded_fields = np.zeros((5, *padded_shape))
    padded_fields[:, :-1, :-1, :-1] = coefficients[:-1]
    gate, on, off, geometric, well = idctn(
        padded_fields, type=1, norm="forward", axes=SPATIAL_AXES, workers=2,
    )

    k = np.arange(cutoff + 1)
    squared = k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2
    kinetic = 4 * np.pi**2 * model.diffusivity * squared
    kept = (slice(None), slice(None), *[slice(cutoff + 1)] * 3)
    padded = np.zeros((2, len(live), *padded_shape))
    output = np.zeros((2, len(live), *squared.shape))
    output[0, :, 0, 0, 0] = weights
    output[1] = weights[:, None, None, None] * coefficients[-1, :cutoff+1, :cutoff+1, :cutoff+1]
    initial = np.zeros_like(output)
    initial[0, launch, 0, 0, 0] = 1.0
    count = 1 + (k > 0)
    output *= count[:, None, None] * count[None, :, None] * count[None, None, :]

    # Minus twice the maximal exits of the killed reaction generator bounds its
    # spectrum in the Gibbs norm pointwise in x. Include the multiplication
    # potential's absolute analytic bound, using |dU/dr| <= 4*E/a.
    radius = model.radius
    well_norm_bound = model.diffusivity * (
        6 * model.depth / radius**2 + 4 * model.depth**2 / radius**2
    )
    maximum_chemical_rate = (
        2 * model.eta + model.activity * model.fuel + model.passive + model.paid * model.fuel
    )
    scale = (np.max(kinetic) + 2 * maximum_chemical_rate + well_norm_bound) / 2
    time = 2.0
    degree = _chebyshev_degree(scale * time, 1e-17)
    cheb_coefficients = ive(np.arange(degree + 1), scale * time)
    cheb_coefficients[1:] *= 2

    results = {}
    for policy in policies:
        pairs = []
        for i, (healthy, fuel) in enumerate(states):
            if fuel:
                j = locations[healthy, fuel - 1]
                multiplier = (fuel - 1) / (model.fuel - 1) if policy == "spectator" else float(policy)
                forward = model.paid * multiplier * fuel
                reverse = model.paid * multiplier * (model.fuel - fuel + 1) / 8
                pairs.append((i, j, forward, reverse, np.sqrt(forward * reverse)))

        def apply(vector):
            padded[kept] = vector
            physical = idctn(padded, type=1, norm="forward", axes=SPATIAL_AXES, workers=2)
            result = np.einsum("ij,qj...->qi...", repair, physical, optimize=True) * gate
            result[1] += well * physical[1]
            result[0] += model.passive * (geometric * physical[1] - on * physical[0])
            result[1] += model.passive * (geometric * physical[0] - off * physical[1])
            for i, j, forward, reverse, pair in pairs:
                result[0, i] += pair * geometric * physical[1, j] - forward * on * physical[0, i]
                result[1, j] += pair * geometric * physical[0, i] - reverse * off * physical[1, j]
            out = dctn(result, type=1, norm="forward", axes=SPATIAL_AXES, workers=2)[kept].copy()
            out += np.einsum("ij,qj...->qi...", damage, vector, optimize=True)
            return out - kinetic * vector

        previous = initial
        current = initial + apply(initial) / scale
        projections = np.empty(degree + 1)
        projections[0] = np.sum(previous * output)
        projections[1] = np.sum(current * output)
        for order in range(2, degree + 1):
            previous, current = current, 2 * (current + apply(current) / scale) - previous
            projections[order] = np.sum(current * output)
        results[str(policy)] = float(projections @ cheb_coefficients)

    return dict(cutoff=cutoff, values=results,
                paired_gap=results["spectator"] - results["1.0"],
                quadrature_error_estimate=integration_error, chebyshev_degree=degree,
                seconds=perf_counter() - start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoffs", nargs="+", type=int, default=[4, 6, 8, 10, 12])
    parser.add_argument("--fixed-intensities", nargs="+", type=float, default=[1.0])
    parser.add_argument("--tolerance", type=float, default=3e-12)
    args = parser.parse_args()
    if any(cutoff < 1 for cutoff in args.cutoffs):
        parser.error("cutoffs must be positive integers")
    if 1.0 not in args.fixed_intensities or any(not 0 <= x <= 1 for x in args.fixed_intensities):
        parser.error("fixed intensities must lie in [0,1] and include 1.0")
    model = Organization()
    rows = []
    for cutoff in args.cutoffs:
        row = solve(model, cutoff, ["spectator", *args.fixed_intensities], args.tolerance)
        print(json.dumps(row), flush=True)
        rows.append(row)
    files = ["scripts/run_continuum_fourier.py", "src/purposeful/organization.py",
             "src/purposeful/spatial_repair.py"]
    payload = dict(
        generated_utc=datetime.now(timezone.utc).isoformat(), model=asdict(model), horizon=2,
        method="Gibbs-conjugated Fourier Galerkin with de-aliased cosine transforms",
        scope="Numerical check; no rigorous continuum error enclosure or global parameter certificate.",
        quadrature_relative_tolerance=args.tolerance, numpy=np.__version__, scipy=scipy.__version__,
        rows=rows, source_sha256={name: sha256((ROOT / name).read_bytes()).hexdigest() for name in files},
    )
    (ROOT / "results" / "continuum-fourier.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8",
    )


if __name__ == "__main__":
    main()
