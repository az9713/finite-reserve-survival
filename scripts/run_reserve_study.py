"""Print reproducible results for the finite-reserve repair study."""
import json
from math import log
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from purposeful.model import sensor_for_accuracy
from purposeful.reserve import (
    RepairContact,
    minimum_energy_for_known_stress_count,
    minimum_tokens,
    stress_count_distribution,
    survival_frontier,
)


def main():
    rounds, interval, q = 2, 1.0, 0.2
    contact = RepairContact(log(9), log(10), 1.0)
    results = {
        "rounds": rounds,
        "exposure_interval": interval,
        "stationary_stress_probability": q,
        "token_energy_kBT": contact.token_energy,
        "repair_rate_ceiling": contact.rate_ceiling,
        "repair_duration": contact.duration,
        "maximum_success_per_repair": contact.success,
        "initial_tokens": list(range(rounds + 1)),
        "survival_target": 0.88,
        "designs": {},
    }
    for label, capacity in [("slow", 0.4), ("fast", 5.8)]:
        sensor = sensor_for_accuracy(0.1, capacity, q)
        mass = stress_count_distribution(rounds, q, sensor.s, interval)
        results["designs"][label] = {
            "mixing_rate": sensor.s,
            "challenge_count_distribution": mass.tolist(),
            "maximum_survival": survival_frontier(mass, contact).tolist(),
            "minimum_tokens_for_target": minimum_tokens(mass, contact, 0.88),
        }
    results["three_certain_stresses"] = {
        "survival_target": 0.4,
        "repair_rate_ceiling": 2.0,
        "repair_duration": 1.0,
        "minimum_initial_reserve_kBT": minimum_energy_for_known_stress_count(3, 0.4, 2.0, 1.0),
    }
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
