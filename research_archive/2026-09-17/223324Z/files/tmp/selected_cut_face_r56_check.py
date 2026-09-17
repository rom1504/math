"""Finite checks for the selected-cut extraction identities in the Wave 56 memo."""

from __future__ import annotations

import math
import random


def kl(p: list[float], q: list[float]) -> float:
    return sum(x * math.log(x / y) for x, y in zip(p, q) if x)


def check_matching() -> None:
    # A nonuniform selector law on a perfect matching incidence.
    n = 37
    raw = [1.0 + ((11 * i) % 17) for i in range(n)]
    ps = [x / sum(raw) for x in raw]
    uniform = [1.0 / n] * n
    # D=S bijectively, hence I(S;D)=H(S).
    entropy = -sum(x * math.log(x) for x in ps)
    selector_kl = kl(ps, uniform)
    assert abs(entropy + selector_kl - math.log(n)) < 1e-12
    u = [1.0 / n] * n
    collision_rhs = sum(ps[i] * math.log(1.0 / u[i]) for i in range(n))
    assert abs(collision_rhs - math.log(n)) < 1e-12


def check_captured_extraction() -> None:
    rng = random.Random(20260730)
    for _ in range(1000):
        size = rng.randrange(2, 40)
        raw = [rng.random() for _ in range(size)]
        nu = [x / sum(raw) for x in raw]
        u = [rng.random() for _ in range(size)]
        row = [0.01 + 100.0 * rng.random() for _ in range(size)]
        deficit = [0.01 + 50.0 * rng.random() for _ in range(size)]
        z = sum(nu[i] * u[i] for i in range(size))
        mu = [nu[i] * u[i] / z for i in range(size)]
        rbar = sum(mu[i] * row[i] for i in range(size))
        dbar = sum(mu[i] * deficit[i] for i in range(size))
        a = b = 3.0
        theta = 1.0 - 1.0 / a - 1.0 / b
        good = [
            i
            for i in range(size)
            if row[i] <= a * rbar and deficit[i] <= b * dbar
        ]
        captured_good = sum(nu[i] * u[i] for i in good)
        assert captured_good + 1e-12 >= theta * z
        # Since nu(good)<=1, at least one good column has u>=theta*z.
        assert good and max(u[i] for i in good) + 1e-12 >= theta * z


if __name__ == "__main__":
    check_matching()
    check_captured_extraction()
    print("matching collision identity: PASS")
    print("selected-prior captured extraction (1000 trials): PASS")
