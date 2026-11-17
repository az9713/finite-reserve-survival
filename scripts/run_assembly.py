"""Reproduce and certify assembly/maintenance optimization before plotting."""
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import heapq
import json
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from purposeful.assembly import Assembly, success_interval, curvature_bounds


def certify(model, deadline=1, lifetime=4, tolerance=F(1, 10000)):
    """Cover the complete unit square with exact interpolation enclosures."""
    values, nodes = {}, []
    best = [F(0), (F(0), F(0))]
    def point(a, z):
        key = (a, z)
        if key not in values:
            values[key] = success_interval(model, deadline, lifetime, a, z)
            if values[key][0] > best[0]:
                best[:] = [values[key][0], key]
        return values[key]
    # This rational construction accelerates pruning. Its value is certified.
    point(F(33, 64), F(1, 4))
    def add(box):
        a0, a1, z0, z1 = box
        corner_upper = max(point(a, z)[1] for a in [a0, a1] for z in [z0, z1])
        ha, hz = curvature_bounds(model, deadline, lifetime, box)
        ea, ez = ha * (a1-a0)**2 / 8, hz * (z1-z0)**2 / 8
        upper = min(F(1), corner_upper + ea + ez)
        identifier = len(nodes)
        nodes.append({'box': list(map(str, box)), 'upper': str(upper), 'children': []})
        heapq.heappush(queue, (-upper, identifier, box, ea, ez))
        return identifier
    queue = []
    add((F(0), F(1), F(0), F(1)))
    iterations = 0
    while -queue[0][0] > best[0] + tolerance:
        _, parent, box, ea, ez = heapq.heappop(queue)
        a0, a1, z0, z1 = box
        if ea >= ez:
            middle = (a0+a1)/2
            children = [(a0, middle, z0, z1), (middle, a1, z0, z1)]
        else:
            middle = (z0+z1)/2
            children = [(a0, a1, z0, middle), (a0, a1, middle, z1)]
        nodes[parent]['children'] = [add(child) for child in children]
        iterations += 1
        if iterations % 1000 == 0:
            print(f'Certified {iterations} splits, {len(values)} points, gap {float(-queue[0][0]-best[0]):.7f}', flush=True)
        if iterations > 150000:
            raise RuntimeError('certificate exceeded its work limit')
    upper = -queue[0][0]
    # Leaves touching the square boundary provide a separate boundary bound.
    boundary_upper = F(0)
    for node in nodes:
        if node['children']:
            continue
        a0, a1, z0, z1 = map(F, node['box'])
        for fixed, value in [('a', F(0)), ('a', F(1)), ('z', F(0)), ('z', F(1))]:
            if fixed == 'a' and value in (a0, a1):
                endpoints = [point(value, z)[1] for z in [z0, z1]]
                _, curvature = curvature_bounds(model, deadline, lifetime, (value, value, z0, z1))
                bound = max(endpoints) + curvature*(z1-z0)**2/8
            elif fixed == 'z' and value in (z0, z1):
                endpoints = [point(a, value)[1] for a in [a0, a1]]
                curvature, _ = curvature_bounds(model, deadline, lifetime, (a0, a1, value, value))
                bound = max(endpoints) + curvature*(a1-a0)**2/8
            else:
                continue
            # The full leaf bound is also valid on its boundary segment.
            boundary_upper = max(boundary_upper, min(bound, F(node['upper'])))
    return {'lower': str(best[0]), 'upper': str(upper), 'tolerance': str(tolerance),
            'construction': list(map(str, best[1])), 'boundary_upper': str(boundary_upper),
            'nodes': nodes,
            'points': [{'alpha': str(a), 'z': str(z), 'lower': str(lo), 'upper': str(hi)}
                       for (a, z), (lo, hi) in sorted(values.items())]}


def main():
    start = time.monotonic()
    model = Assembly(fuel=4)
    certificate = certify(model)
    result = {'generated_utc': datetime.now(timezone.utc).isoformat(),
              'model': {key: str(value) for key, value in vars(model).items()},
              'deadline': 1, 'lifetime': 4, 'certificate': certificate,
              'elapsed_seconds': time.monotonic()-start,
              'source_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                for name in ['src/purposeful/assembly.py', 'scripts/run_assembly.py']}}
    (ROOT / 'results/assembly.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print({key: value for key, value in certificate.items() if key not in ['nodes', 'points']}, flush=True)
    print('seconds', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
