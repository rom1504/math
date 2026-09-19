"""Finite sanity checks for Theorem A in diagnostic_audit_report.md.

This is not a substitute for the one-line Markov proof.  It samples arbitrary
incidence-column masses and costs and verifies the explicit extraction
constants, including the two-cost variant.
"""

from __future__ import annotations

import random


def check_one_cost(rng: random.Random, trials: int = 100_000) -> None:
    for _ in range(trials):
        size = rng.randrange(1, 40)
        u = [rng.random() for _ in range(size)]
        cost = [10 ** rng.uniform(-3, 5) for _ in range(size)]
        raw = [rng.expovariate(1.0) for _ in range(size)]
        total = sum(raw)
        nu = [x / total for x in raw]
        z = sum(p * x for p, x in zip(nu, u))
        captured = sum(p * x * c for p, x, c in zip(nu, u, cost)) / z
        a = 1.0 + 10 ** rng.uniform(-2, 2)
        theta = 1.0 - 1.0 / a
        assert any(x >= theta * z and c <= a * captured
                   for x, c in zip(u, cost))


def check_two_costs(rng: random.Random, trials: int = 100_000) -> None:
    for _ in range(trials):
        size = rng.randrange(1, 40)
        u = [rng.random() for _ in range(size)]
        cost_1 = [10 ** rng.uniform(-3, 5) for _ in range(size)]
        cost_2 = [10 ** rng.uniform(-3, 5) for _ in range(size)]
        raw = [rng.expovariate(1.0) for _ in range(size)]
        total = sum(raw)
        nu = [x / total for x in raw]
        z = sum(p * x for p, x in zip(nu, u))
        mean_1 = sum(p * x * c for p, x, c in zip(nu, u, cost_1)) / z
        mean_2 = sum(p * x * c for p, x, c in zip(nu, u, cost_2)) / z
        inv_a = rng.uniform(0.01, 0.48)
        inv_b = rng.uniform(0.01, 0.49 - inv_a)
        a, b = 1.0 / inv_a, 1.0 / inv_b
        theta = 1.0 - inv_a - inv_b
        assert any(
            x >= theta * z and c1 <= a * mean_1 and c2 <= b * mean_2
            for x, c1, c2 in zip(u, cost_1, cost_2)
        )


if __name__ == "__main__":
    random_source = random.Random(0xA0D17)
    check_one_cost(random_source)
    check_two_costs(random_source)
    print("diagnostic extraction checks passed")
