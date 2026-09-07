#!/usr/bin/env python3
"""Independent exact audit of the Wave 58 relative-cut formulas."""

from fractions import Fraction
from itertools import combinations, product
from math import comb
import random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def matvec(a, x):
    return [dot(row, x) for row in a]


def falling(a, k):
    ans = 1
    for j in range(k):
        ans *= a - j
    return ans


def trial(n, m, rng):
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = rng.choice((-1, 1))

    x = [rng.choice((-1, 1)) for _ in range(n)]
    z = [rng.choice((-1, 1)) for _ in range(n)]
    tau = rng.choice((-1, 1))
    upsilon = rng.choice((-1, 1))
    kappa = tau * upsilon

    s = [[0] * n for _ in range(n)]
    sp = [[0] * n for _ in range(n)]
    d_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            s[i][j] = s[j][i] = tau * a[i][j] * x[i] * x[j]
            sp[i][j] = sp[j][i] = upsilon * a[i][j] * z[i] * z[j]
            in_d = kappa * (x[i] * z[i]) * (x[j] * z[j]) == -1
            assert (sp[i][j] == -s[i][j]) == in_d
            if in_d:
                d_edges.append((i, j))

    r = [sum(row) for row in s]
    rp = [sum(row) for row in sp]
    b = [sum(s[i][j] for j in range(n) if i != j and sp[i][j] == -s[i][j])
         for i in range(n)]
    assert rp == [r[i] - 2 * b[i] for i in range(n)]
    assert sum(v * v for v in rp) - sum(v * v for v in r) == (
        4 * sum(v * v for v in b) - 4 * dot(r, b)
    )

    w = sum(s[i][j] for i, j in d_edges)
    energy = sum(x[i] * a[i][j] * x[j] for i in range(n) for j in range(n)) * tau
    energy_p = sum(z[i] * a[i][j] * z[j] for i in range(n) for j in range(n)) * upsilon
    # Delta=q-energy, so the unknown common q cancels.
    assert -(energy_p - energy) == 4 * w

    p2 = Fraction(falling(m, 2), falling(n, 2))
    values = []
    for ss in combinations(range(n), m):
        ss = set(ss)
        t = sum(s[i][j] for i, j in d_edges if i in ss and j in ss)
        local_energy = tau * sum(
            x[i] * a[i][j] * x[j] for i in ss for j in ss
        )
        local_energy_p = upsilon * sum(
            z[i] * a[i][j] * z[j] for i in ss for j in ss
        )
        # delta'=Q-local_energy', so local Q cancels.
        delta_change = -(local_energy_p - local_energy)
        assert delta_change == 4 * t
        # hat-ell change = delta change - p2 * Delta change.
        assert Fraction(delta_change) - p2 * 4 * w == 4 * (Fraction(t) - p2 * w)
        values.append(t)

    mean = Fraction(sum(values), len(values))
    assert mean == p2 * w
    variance = sum((Fraction(v) - mean) ** 2 for v in values) / len(values)
    p3 = Fraction(falling(m, 3), falling(n, 3)) if n >= 3 else Fraction(0)
    p4 = Fraction(falling(m, 4), falling(n, 4)) if n >= 4 else Fraction(0)
    predicted = (
        (p2 - 2 * p3 + p4) * len(d_edges)
        + (p3 - p4) * sum(v * v for v in b)
        + (p4 - p2 * p2) * w * w
    )
    assert variance == predicted, (variance, predicted)


def main():
    rng = random.Random(580058)
    for n in range(5, 10):
        for m in range(2, n):
            for _ in range(25):
                trial(n, m, rng)
    print("relative row, deficit, fibre-change, and exact variance checks passed")


if __name__ == "__main__":
    main()
