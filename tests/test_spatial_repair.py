"""Independent references for the periodic two-complex maintenance pilot."""

from itertools import product
from math import comb, log, pi, sqrt, tanh
import unittest

import numpy as np
from scipy.linalg import expm, solve_banded

from purposeful.encounter import mean_capture_time
from purposeful.spatial_repair import (
    SpatialRepair,
    build_grid,
    calibrate,
    chemical_generators,
    count_survival,
    equilibrium,
    spatial_generator,
    spatial_survival,
)


def labelled_channels(state, model):
    """One independent channel per labelled complex and reaction direction."""
    first, second, fuel = state
    labels = (first, second)
    for index, healthy in enumerate(labels):
        new_labels = list(labels)
        new_labels[index] = not healthy
        yield (*new_labels, fuel), model.eta * (1 if healthy else 0.5), "damage"
        if healthy and labels[1 - index] and fuel < model.fuel:
            yield (*new_labels, fuel + 1), model.activity * (model.fuel - fuel) / 4, "repair"
        if not healthy and labels[1 - index] and fuel:
            yield (*new_labels, fuel - 1), model.activity * fuel, "repair"


def periodic_reference(model, n):
    """Unfolded diffusion matrix and radial quantities at periodic grid nodes."""
    points = list(product(range(n), repeat=3))
    index = {position: i for i, position in enumerate(points)}
    diffusion = np.zeros((len(points), len(points)))
    hop = model.diffusivity * (n / model.length) ** 2
    for row, point in enumerate(points):
        for axis, direction in product(range(3), (-1, 1)):
            neighbor = list(point)
            neighbor[axis] = (neighbor[axis] + direction) % n
            diffusion[row, index[tuple(neighbor)]] += hop
            diffusion[row, row] -= hop
    coordinates = np.asarray(points) * (model.length / n)
    signed = np.where(coordinates >= model.length / 2, coordinates - model.length, coordinates)
    radius = np.linalg.norm(signed, axis=1)
    gate = np.maximum(0, 1 - np.square(radius / model.radius)) ** 4
    return diffusion, radius, gate


def labelled_reference(model, n, killed=True):
    """Reference generator in full spatial coordinates and full label space."""
    diffusion, radius, gate = periodic_reference(model, n)
    states = [(left, right, fuel) for left, right in product((False, True), repeat=2)
              for fuel in range(model.fuel + 1) if not killed or left or right]
    index = {state: row for row, state in enumerate(states)}
    sites = len(gate)
    generator = np.kron(np.eye(len(states)), diffusion)
    spatial_indices = np.arange(sites)
    for state, chem_index in index.items():
        rows = chem_index * sites + spatial_indices
        for destination, base_rate, reaction in labelled_channels(state, model):
            rate = base_rate * gate if reaction == "repair" else base_rate
            generator[rows, rows] -= rate
            if destination in index:
                columns = index[destination] * sites + spatial_indices
                generator[rows, columns] += rate
    return states, generator, radius, gate


def doi_shell_volume_mean(diffusion, inner, outer, activity, level):
    """Independent radial cell balance with the hazard jump on a cell face."""
    faces = np.r_[np.linspace(0, inner, level + 1), np.linspace(inner, outer, level + 1)[1:]]
    lower, upper = faces[:-1], faces[1:]
    volume = (upper**3 - lower**3) / 3
    centers = (upper**4 - lower**4) / (4 * volume)
    conductance = diffusion * faces[1:-1]**2 / np.diff(centers)
    band = np.zeros((3, len(centers)))
    band[0, 1:] = band[2, :-1] = -conductance
    band[1, :-1] += conductance
    band[1, 1:] += conductance
    band[1, :level] += activity * volume[:level]
    mean = solve_banded((1, 1), band, volume)
    return np.average(mean[level:], weights=volume[level:])


def radial_verification():
    diffusion, inner, outer, activity = 0.10, 0.25, 1.0, 2.2
    shell_volume = 4 * pi * (outer**3 - inner**3) / 3
    absorbing_mean = mean_capture_time(shell_volume, diffusion, inner)

    def finite_hazard_mean(rate):
        mu = sqrt(rate / diffusion)
        surface_mean = (outer**3 - inner**3) / (3 * diffusion * inner**2)
        return absorbing_mean + 1 / rate + surface_mean / (mu / tanh(mu * inner) - 1 / inner)

    analytic = finite_hazard_mean(activity)
    levels = (20, 40, 80, 160, 320, 640)
    numerical = [float(doi_shell_volume_mean(diffusion, inner, outer, activity, level)) for level in levels]
    return {
        "diffusivity": diffusion,
        "contact_radius": inner,
        "reflecting_radius": outer,
        "step_hazard": activity,
        "analytic_exterior_mean": analytic,
        "cells_per_region": list(levels),
        "numerical_exterior_mean": numerical,
        "absolute_errors": [abs(mean - analytic) for mean in numerical],
        "stage_1_absorbing_exterior_mean": absorbing_mean,
        "hard_contact_hazard": 1e10,
        "hard_contact_exterior_mean": finite_hazard_mean(1e10),
    }


class SpatialRepairTests(unittest.TestCase):
    def test_chemical_channels_and_energy_bookkeeping(self):
        model = SpatialRepair(fuel=3)
        states, damage, contact = chemical_generators(model)
        index = {state: i for i, state in enumerate(states)}
        for row, (healthy, fuel) in enumerate(states):
            label_state = (healthy > 0, healthy > 1, fuel)
            reference = {name: np.zeros(len(states)) for name in ("damage", "repair")}
            for destination, rate, reaction in labelled_channels(label_state, model):
                new_fuel = destination[2]
                new_healthy = sum(destination[:2])
                column = index[new_healthy, new_fuel]
                reference[reaction][column] += rate
                reference[reaction][row] -= rate
                dh, db = new_healthy - healthy, new_fuel - fuel
                old_counts = np.array([healthy, 2 - healthy, fuel, model.fuel - fuel])
                new_counts = np.array([new_healthy, 2 - new_healthy, new_fuel, model.fuel - new_fuel])
                forward = -dh if reaction == "damage" else dh
                expected = [-1, 1, 0, 0] if reaction == "damage" else [1, -1, -1, 1]
                np.testing.assert_array_equal(new_counts - old_counts, forward * np.array(expected))
                self.assertEqual(sum(new_counts[:2]), sum(old_counts[:2]))
                self.assertEqual(sum(new_counts[2:]), sum(old_counts[2:]))
                self.assertEqual(dh + db, -forward if reaction == "damage" else 0)
                old_energy = fuel * log(8) - (2 - healthy) * log(2)
                new_energy = new_fuel * log(8) - (2 - new_healthy) * log(2)
                change = -log(2) if reaction == "damage" else -log(4)
                self.assertAlmostEqual(new_energy - old_energy, forward * change, places=14)
            np.testing.assert_allclose(damage[row], reference["damage"], atol=1e-13)
            np.testing.assert_allclose(contact[row], reference["repair"], atol=1e-13)

    def test_uncensored_detailed_balance_and_killing(self):
        model = SpatialRepair(fuel=2)
        grid = build_grid(model, 6)
        states, full = spatial_generator(model, grid, killed=False)
        live_states, killed = spatial_generator(model, grid, killed=True)
        sites = len(grid.weight)
        exact_chem = np.array([
            comb(2, h) * (1 / 3)**h * (2 / 3)**(2 - h)
            * comb(model.fuel, b) * (1 / 9)**b * (8 / 9)**(model.fuel - b)
            for h, b in states
        ])
        np.testing.assert_allclose(equilibrium(model), exact_chem, atol=1e-14)
        self.assertAlmostEqual(sum(weight for (h, _), weight in zip(states, exact_chem) if h == 0), 4 / 9)
        rows = np.concatenate([i * sites + np.arange(sites) for i, (h, _) in enumerate(states) if h])
        np.testing.assert_allclose(killed.toarray(), full[rows, :][:, rows].toarray(), atol=1e-14)
        np.testing.assert_allclose(np.asarray(full.sum(axis=1)).ravel(), 0, atol=1e-12)
        np.testing.assert_allclose(np.asarray(killed.sum(axis=1)).ravel(),
                                   np.repeat([-model.eta if h == 1 else 0 for h, _ in live_states], sites),
                                   atol=1e-12)
        stationary = np.kron(exact_chem, grid.weight)
        np.testing.assert_allclose(stationary @ full, 0, atol=1e-13)
        flow = stationary[:, None] * full.toarray()
        np.testing.assert_allclose(flow, flow.T, atol=1e-13)

    def test_full_label_and_periodic_grid_reference(self):
        model = SpatialRepair(fuel=1, radius=0.40, activity=8.)
        n = 4
        grid = build_grid(model, n)
        states, matrix, radius, gate = labelled_reference(model, n)
        times = [0., .2, .7]
        reference_weights = {
            "uniform": np.ones_like(gate),
            "localized": np.exp(-0.5 * np.square(radius / model.local_width)),
            "reaction": gate,
        }
        site_count = n**3
        for initial, labelled in [((1, 1), (True, False, 1)), ((2, 1), (True, True, 1))]:
            native = spatial_survival(model, grid, times, initial)
            start = states.index(labelled) * site_count
            for j, time in enumerate(times):
                survival = expm(time * matrix).sum(axis=1)[start:start + site_count]
                for prep, weights in reference_weights.items():
                    self.assertAlmostEqual(native[prep][j], np.average(survival, weights=weights), places=11)

    def test_clamped_forward_calibration_and_count_reverses(self):
        model = SpatialRepair(fuel=2)
        n = 4
        grid = build_grid(model, n)
        diffusion, _, gate = periodic_reference(model, n)
        effective = calibrate(model, grid)
        self.assertTrue(np.isnan(effective[0]))
        for b in range(1, model.fuel + 1):
            mean = np.linalg.solve(np.diag(model.activity * b * gate) - diffusion, np.ones(n**3)).mean()
            self.assertAlmostEqual(effective[b], 1 / (b * mean), places=12)
            self.assertLess(effective[b], model.activity * np.mean(gate))
        states = [(h, b) for h in (1, 2) for b in range(model.fuel + 1)]
        index = {state: row for row, state in enumerate(states)}
        matrix = np.zeros((len(states), len(states)))
        for (h, b), row in index.items():
            for destination, rate, reaction in labelled_channels((h > 0, h > 1, b), model):
                new_b = destination[2]
                target = (sum(destination[:2]), new_b)
                if reaction == "repair":
                    rate *= effective[max(b, new_b)] / model.activity
                matrix[row, row] -= rate
                if target in index:
                    matrix[row, index[target]] += rate
        times = [0., .1, 1., 2.]
        for initial in [(1, 2), (2, 1)]:
            expected = [expm(time * matrix)[index[initial]].sum() for time in times]
            np.testing.assert_allclose(count_survival(model, times, effective, initial), expected, atol=1e-13)

    def test_fast_diffusion_approaches_stationary_average(self):
        times = [0., .1, .25, .5]
        errors = []
        for diffusion in (0.2, 2.0, 20.0):
            model = SpatialRepair(fuel=1, diffusivity=diffusion, activity=12.)
            grid = build_grid(model, 4)
            _, _, gate = periodic_reference(model, 4)
            rate = model.activity * gate.mean()
            target = count_survival(model, times, np.array([np.nan, rate]), (1, 1))
            actual = spatial_survival(model, grid, times, (1, 1))["uniform"]
            errors.append(np.max(np.abs(actual - target)))
        self.assertLess(errors[1], errors[0] / 3)
        self.assertLess(errors[2], errors[1] / 4)
        self.assertLess(errors[2], 2e-4)

    def test_no_fuel_and_disabled_repair_reduce_to_damage_only(self):
        spatial_model = SpatialRepair(fuel=0, eta=0.7)
        count_model = SpatialRepair(fuel=3, eta=spatial_model.eta)
        slow = spatial_model.eta * (-7 + sqrt(17)) / 4
        fast = spatial_model.eta * (-7 - sqrt(17)) / 4
        times = np.array([0, .01, .2, 1, 4])
        grid = build_grid(spatial_model, 4)
        np.testing.assert_array_equal(calibrate(spatial_model, grid), [np.nan])
        for h, initial_derivative in ((1, -spatial_model.eta), (2, 0)):
            slow_weight = (initial_derivative - fast) / (slow - fast)
            expected = slow_weight * np.exp(slow * times) + (1 - slow_weight) * np.exp(fast * times)
            actual = spatial_survival(spatial_model, grid, times, (h, 0))
            for curve in actual.values():
                np.testing.assert_allclose(curve, expected, atol=2e-14, rtol=0)
            for b in range(count_model.fuel + 1):
                curve = count_survival(count_model, times, np.zeros(count_model.fuel + 1), (h, b))
                np.testing.assert_allclose(curve, expected, atol=2e-14, rtol=0)

    def test_short_time_survival_curvature_resolves_contact(self):
        model = SpatialRepair(fuel=2, radius=0.40, activity=5.)
        n = 6
        grid = build_grid(model, n)
        _, radius, gate = periodic_reference(model, n)
        reference_weights = {
            "uniform": np.ones_like(gate),
            "localized": np.exp(-0.5 * np.square(radius / model.local_width)),
            "reaction": gate,
        }
        dt = 1e-4
        times = dt * np.arange(5)
        for h, b in [(1, 2), (2, 1)]:
            curves = spatial_survival(model, grid, times, (h, b))
            for law, curve in curves.items():
                mean_gate = np.average(gate, weights=reference_weights[law])
                first = np.dot([-25, 48, -36, 16, -3], curve) / (12 * dt)
                second = np.dot([35, -104, 114, -56, 11], curve) / (12 * dt**2)
                if h == 1:
                    expected = model.eta * (1.5 * model.eta + b * model.activity * mean_gate)
                    self.assertAlmostEqual(first, -model.eta, delta=1e-8)
                else:
                    expected = -model.eta * (2 * model.eta + (model.fuel - b) * model.activity * mean_gate / 2)
                    self.assertAlmostEqual(first, 0, delta=1e-8)
                self.assertAlmostEqual(second, expected, delta=2e-5)

    def test_invalid_parameters_mesh_and_failure_start(self):
        for parameters in ({"fuel": -1}, {"fuel": 2.5}, {"fuel": True}, {"eta": 0},
                           {"diffusivity": float("nan")}, {"radius": 0.6}, {"local_width": -1}):
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                SpatialRepair(**parameters)
        model = SpatialRepair()
        for n in (0, 1, 3, 4.0, True):
            with self.subTest(n=n), self.assertRaises(ValueError):
                build_grid(model, n)
        grid = build_grid(model, 4)
        with self.assertRaises(ValueError):
            spatial_survival(model, grid, [-1])
        with self.assertRaises(ValueError):
            spatial_survival(model, grid, [0], (0, model.fuel))
        with self.assertRaises(ValueError):
            count_survival(model, [0], np.full(model.fuel + 1, -1))

    def test_spherical_doi_reference_and_hard_contact_limit(self):
        record = radial_verification()
        errors = record["absolute_errors"]
        self.assertTrue(all(new < old / 3 for old, new in zip(errors, errors[1:])), errors)
        self.assertLess(errors[-1], 1e-4)
        self.assertGreater(record["analytic_exterior_mean"], record["hard_contact_exterior_mean"])
        self.assertLess(record["hard_contact_exterior_mean"] - record["stage_1_absorbing_exterior_mean"], 2e-4)


if __name__ == "__main__":
    unittest.main()
