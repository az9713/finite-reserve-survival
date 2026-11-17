"""Independent continuum particle check for the closed organization benchmark."""
import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Physics:
    length: float = 1.
    radius: float = .32
    diffusivity: float = .025
    depth: float = 6.
    offset: float = 2.
    fuel: int = 4
    eta: float = 1.
    repair: float = 16.
    passive: float = 1.
    paid: float = 4.


def advance(x, closed, duration, max_dt, p, rng):
    """Conditional propagation at fixed conformation, using minimum-image positions."""
    open_ids = np.flatnonzero(~closed)
    increments = rng.standard_normal((len(open_ids), 3))
    x[open_ids] += np.sqrt(2 * p.diffusivity * duration[open_ids])[:, None] * increments
    x[open_ids] = (x[open_ids] + p.length / 2) % p.length - p.length / 2
    closed_ids = np.flatnonzero(closed)
    if not len(closed_ids):
        return x
    count = np.maximum(1, np.ceil(duration[closed_ids] / max_dt).astype(int))
    dt = duration[closed_ids] / count
    for k in range(count.max()):
        mask = k < count
        ids = closed_ids[mask]
        pos = x[ids]
        step = dt[mask, None]
        core = np.maximum(0., 1. - np.sum(pos ** 2, axis=1) / p.radius ** 2)
        drift = (-4 * p.diffusivity * p.depth / p.radius ** 2) * core[:, None] * pos
        noise = np.sqrt(2 * p.diffusivity * step) * rng.standard_normal(pos.shape)
        x[ids] = (pos + step * drift + noise + p.length / 2) % p.length - p.length / 2
    return x


def run_paths(p, policy, max_dt, seed, count, times):
    rng = np.random.Generator(np.random.PCG64(seed))
    x = rng.uniform(-p.length / 2, p.length / 2, size=(count, 3))
    first = np.ones(count, dtype=bool)
    second = np.zeros(count, dtype=bool)
    closed = np.zeros(count, dtype=bool)
    fuel = np.full(count, p.fuel, dtype=np.int16)
    elapsed = np.zeros(count)
    net_damage = np.zeros(count, dtype=np.int32)
    channel_names = ('repair_forward', 'repair_reverse', 'passive_forward',
                     'passive_reverse', 'paid_forward', 'paid_reverse')
    gross = {name: np.zeros(count, dtype=np.int32) for name in channel_names}
    curve = {}
    errors = {}
    totals = {'poisson_candidates': 0, 'accepted_chemical_jumps': 0,
              'observed_path_ledger_checks': 0, 'max_path_ledger_error': 0}

    def persist(t):
        h = first.astype(int) + second
        alive = h > 0
        net_r = gross['repair_forward'] - gross['repair_reverse']
        net_p = gross['paid_forward'] - gross['paid_reverse']
        net_0 = gross['passive_forward'] - gross['passive_reverse']
        error = max(np.max(abs(fuel - p.fuel + net_r + net_p)),
                    np.max(abs(closed.astype(int) - net_p - net_0)),
                    np.max(abs(h - 1 + net_damage - net_r)))
        totals['observed_path_ledger_checks'] += count
        totals['max_path_ledger_error'] = max(totals['max_path_ledger_error'], int(error))
        if error or fuel.min() < 0 or fuel.max() > p.fuel:
            raise AssertionError('invalid inventory, health or tether trajectory ledger')
        r2 = np.sum(x ** 2, axis=1) / p.radius ** 2
        obs = {'survival': alive, 'live_contact': alive & (r2 < 1),
               'live_tether': alive & closed, 'live_gate': alive * np.maximum(0, 1-r2) ** 4,
               'stopped_fuel': fuel, 'stopped_tether': closed, **gross}
        for name, vector in obs.items():
            curve.setdefault(name, []).append(float(np.mean(vector)))
            errors.setdefault(name, []).append(float(np.std(vector, ddof=1) / np.sqrt(count)))

    persist(times[0])
    for horizon in times[1:]:
        ids = np.flatnonzero(first | second)
        while len(ids):
            h1, h2 = first[ids], second[ids]
            q, b = closed[ids], fuel[ids]
            w = p.fuel - b
            thermal1 = p.eta * np.where(h1, 1., .5)
            thermal2 = p.eta * np.where(h2, 1., .5)
            rep = p.repair * (h1 ^ h2) * b
            rev = p.repair * (h1 & h2) * w / 4  # One reverse per labelled target.
            if policy == 'spectator':
                # Choosing the forward substrate leaves b-1 spectator choices.
                # In reverse the substrate is waste, leaving all b fuel as spectators.
                on_scale = np.maximum(0, b - 1) / (p.fuel - 1)
                off_scale = b / (p.fuel - 1)
            else:
                on_scale = off_scale = policy
            paid_on = p.paid * (~q) * b * on_scale
            paid_off = p.paid * q * w * off_scale / 8
            bound = thermal1 + thermal2 + rep + 2 * rev + p.passive + paid_on + paid_off
            remaining = horizon - elapsed[ids]
            wait = rng.exponential(1 / bound)
            before = wait < remaining
            duration = np.minimum(wait, remaining)
            x[ids] = advance(x[ids], q, duration, max_dt, p, rng)
            elapsed[ids] += duration
            elapsed[ids[~before]] = horizon
            if not np.any(before):
                break
            ids = ids[before]
            proposals = len(ids)
            totals['poisson_candidates'] += proposals
            h1, h2, q = h1[before], h2[before], q[before]
            core = np.maximum(0., 1 - np.sum(x[ids] ** 2, axis=1) / p.radius ** 2)
            g = core ** 4
            potential = p.offset - p.depth * core ** 2
            on = g * np.exp(-np.maximum(0., potential))
            off = g * np.exp(np.minimum(0., potential))
            # Marks above all actual cumulative intensities are rejected candidates.
            rates = (thermal1[before], thermal2[before], rep[before] * g,
                     rev[before] * g, rev[before] * g,
                     p.passive * (~q) * on, p.passive * q * off,
                     paid_on[before] * on, paid_off[before] * off)
            bounds = np.cumsum(np.asarray(rates), axis=0)
            if np.any(bounds[-1] > bound[before] + 1e-12):
                raise AssertionError('chemical thinning bound violated')
            marks = rng.random(proposals) * bound[before]
            which = np.sum(marks >= bounds, axis=0)
            totals['accepted_chemical_jumps'] += int(np.sum(which < len(rates)))
            idx = ids[which == 0]
            net_damage[idx] += np.where(first[idx], 1, -1)
            first[idx] = ~first[idx]
            idx = ids[which == 1]
            net_damage[idx] += np.where(second[idx], 1, -1)
            second[idx] = ~second[idx]
            idx = ids[which == 2]
            first[idx] = second[idx] = True
            fuel[idx] -= 1
            gross['repair_forward'][idx] += 1
            r1, r2 = ids[which == 3], ids[which == 4]
            first[r1] = False
            second[r2] = False
            idx = ids[(which == 3) | (which == 4)]
            fuel[idx] += 1
            gross['repair_reverse'][idx] += 1
            for channel, code, new_q, fuel_delta in (
                    ('passive_forward', 5, True, 0), ('passive_reverse', 6, False, 0),
                    ('paid_forward', 7, True, -1), ('paid_reverse', 8, False, 1)):
                idx = ids[which == code]
                closed[idx] = new_q
                fuel[idx] += fuel_delta
                gross[channel][idx] += 1
            ids = ids[first[ids] | second[ids]]
        persist(horizon)
    return dict(policy=policy, max_euler_step=max_dt, seed=seed, trajectories=count,
                curve=curve, marginal_standard_error=errors, **totals)


def lattice_reference(n, times):
    """Read only the independent comparator. The SDE solver above has no lattice calls."""
    sys.path.insert(0, str(ROOT / 'src'))
    from purposeful.organization import Organization, build_system
    system = build_system(Organization(), n)
    mapping = {'survival': 'mass', 'live_contact': 'contact', 'live_tether': 'tether',
               'live_gate': 'gate', 'stopped_fuel': 'stopped_fuel',
               'stopped_tether': 'stopped_tether'}
    result = {}
    for name, policy in [('passive', 0.), ('fixed_one', 1.), ('spectator', 'spectator')]:
        values = system.observables(policy, times)
        result[name] = {key: values[v].tolist() for key, v in mapping.items()}
    return result


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paths', type=int, default=20000)
    parser.add_argument('--n', type=int, default=64)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/organization-continuum-check.json')
    args = parser.parse_args()
    if args.paths < 2:
        parser.error('--paths requires at least 2')
    p = Physics()
    times = np.linspace(0., 2., 9)
    dt_values = [.01, .005, .0025]
    ref = lattice_reference(args.n, times)
    alpha = .01
    family_count = len(dt_values) * 3
    band = float(np.sqrt(np.log(2 * family_count / alpha) / (2 * args.paths)))
    runs = []
    for d, dt in enumerate(dt_values):
        for k, (name, policy) in enumerate([('passive', 0.), ('fixed_one', 1.),
                                            ('spectator', 'spectator')]):
            seed = 2026111700 + 3 * d + k
            run = run_paths(p, policy, dt, seed, args.paths, times)
            run['name'] = name
            run['maximum_abs_survival_gap_to_lattice'] = float(np.max(np.abs(
                np.array(run['curve']['survival']) - ref[name]['survival'])))
            runs.append(run)
            print(name, dt, 'survival', run['curve']['survival'][-1],
                  'max lattice gap', run['maximum_abs_survival_gap_to_lattice'], flush=True)
    result = dict(physics=asdict(p), relative_position='3D unit torus, uniform initial displacement',
                  initial={'healthy_first': True, 'healthy_second': False, 'fuel': 4, 'closed': False},
                  random_generator='NumPy PCG64, independent seed for each of nine batches',
                  times=times.tolist(), lattice_points_per_axis=args.n, lattice_reference=ref, runs=runs,
                  confidence={'method': 'Two-sided DKW over survival time, union bound over 9 runs.',
                              'family_confidence': 1-alpha, 'simultaneous_half_width': band,
                              'scope': 'Sampling uncertainty for each numerical SDE law. Euler drift bias and lattice error unenclosed.'},
                  limitations=['Tiny survival differences between policies are unresolved at this sampling precision.',
                               'Free Brownian increments are exact. Tethered drift uses Euler-Maruyama between proposal/observation times.',
                               'Chemical thinning uses a valid dominating rate, evaluated on the numerically propagated positions.',
                               'Separate seeds permit sampling noise in differences between time steps.',
                               'Marginal standard errors for other observables are descriptive and are not simultaneous confidence bounds.'],
                  software={'numpy': np.__version__, 'scipy': scipy.__version__, 'python': sys.version.split()[0]},
                  source_sha256={'scripts/organization_continuum_check.py': digest(Path(__file__)),
                                 'src/purposeful/organization.py': digest(ROOT/'src/purposeful/organization.py'),
                                 'src/purposeful/spatial_repair.py': digest(ROOT/'src/purposeful/spatial_repair.py')})
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
