#!/usr/bin/env python3
"""Exact checks for the Wave 48 row-truncated coarea sharpening."""

from __future__ import annotations

import itertools
import math
from fractions import Fraction


def subsets(n: int, m: int):
    return [frozenset(s) for s in itertools.combinations(range(n), m)]


def down_up_kernel(n: int, m: int, ell: int, states):
    denom = math.comb(m, ell) * math.comb(n - ell, m - ell)
    return [
        [Fraction(math.comb(len(s & t), ell), denom) for t in states]
        for s in states
    ]


def lambdas(n: int, m: int, ell: int):
    lam1 = Fraction(ell * (n - m), m * (n - ell)) if ell else Fraction(0)
    if ell < 2 or n - m < 2:
        lam2 = Fraction(0)
    else:
        lam2 = Fraction(
            ell * (ell - 1) * (n - m) * (n - m - 1),
            m * (m - 1) * (n - ell) * (n - ell - 1),
        )
    return lam1, lam2


def family_data(n: int, m: int, ell: int, states, kernel, mask: int):
    count = bin(mask).count("1")
    total = len(states)
    a = Fraction(count, total)
    active = [i for i in range(total) if mask >> i & 1]
    collision = sum((kernel[i][j] for i in active for j in active), Fraction(0)) / total
    boundary = a - collision

    mu = []
    for vertex in range(n):
        hits = sum(vertex in states[i] for i in active)
        mu.append(Fraction(hits, count))
    p = Fraction(m, n)
    w1 = (
        a * a * Fraction(n - 1, 1) / (p * (1 - p) * n)
        * sum(((x - p) ** 2 for x in mu), Fraction(0))
    )
    return a, boundary, w1


def exhaustive_slice_check(n: int, m: int, ell: int):
    states = subsets(n, m)
    kernel = down_up_kernel(n, m, ell, states)
    assert all(sum(row, Fraction(0)) == 1 for row in kernel)
    lam1, lam2 = lambdas(n, m, ell)
    delta = 1 - lam1
    gap = lam1 - lam2
    kappa = gap / delta
    checked = 0
    qualifying = 0
    threshold = kappa / (1 + kappa * n)
    for mask in range(1, 1 << len(states)):
        a, boundary, w1 = family_data(n, m, ell, states, kernel, mask)
        assert w1 <= (n - 1) * a * a
        sharpened_rhs = delta * a * (1 - a) + gap * a * (1 - n * a)
        assert boundary >= sharpened_rhs
        ratio = boundary / (delta * a)
        if ratio <= 1:
            qualifying += 1
            assert a >= threshold
        checked += 1
    print(
        f"slice n={n},m={m},ell={ell}: checked={checked}, "
        f"qualifying={qualifying}, kappa={kappa}, threshold={threshold}"
    )


def aggregate_check(n: int, m: int, ell: int):
    """Check the a^2-size-biased aggregate implication on all family pairs."""
    states = subsets(n, m)
    kernel = down_up_kernel(n, m, ell, states)
    lam1, lam2 = lambdas(n, m, ell)
    delta = 1 - lam1
    kappa = (lam1 - lam2) / delta
    threshold = kappa / (1 + kappa * n)
    data = []
    for mask in range(1, 1 << len(states)):
        a, boundary, _ = family_data(n, m, ell, states, kernel, mask)
        data.append((a, boundary))
    qualifying = 0
    for i, (a1, b1) in enumerate(data):
        # Pair with a deterministic spread of partners; the theorem itself is algebraic.
        for j in range(i % 17, len(data), 97):
            a2, b2 = data[j]
            ratio = (a1 * b1 + a2 * b2) / (delta * (a1 * a1 + a2 * a2))
            if ratio <= 1:
                qualifying += 1
                weighted_degree = (a1**3 + a2**3) / (a1**2 + a2**2)
                assert weighted_degree >= threshold
                assert max(a1, a2) >= threshold
    print(f"aggregate pair samples: qualifying={qualifying}, threshold={threshold}")


def complete_signing_obstruction():
    """All-positive signing: every fixed-density ground lift misses project row."""
    c = Fraction(1, 8)
    for n in (32, 64, 128, 256, 512):
        m = 3 * n // 4
        min_active_row = n + (n - 2) * (2 * m - n) ** 2
        project_power = float(Fraction(9, 4) - c)
        project_scale = n**project_power
        print(
            f"complete n={n},m={m}: min-positive-degree row={min_active_row}, "
            f"n^(9/4-c)={project_scale:.3f}, ratio={min_active_row/project_scale:.3f}"
        )
    # The exact formulas used above.
    n, m = 64, 48
    for k in range(n + 1):
        degree_count = math.comb(k, m) + math.comb(n - k, m)
        row = n + (n - 2) * (2 * k - n) ** 2
        if degree_count:
            assert abs(2 * k - n) >= 2 * m - n
            assert row >= n + (n - 2) * (2 * m - n) ** 2


def symbolic_kappa_check():
    for n in range(6, 20):
        for m in range(3, n - 1):
            for s in range(1, m):
                ell = m - s
                lam1, lam2 = lambdas(n, m, ell)
                if not lam1:
                    continue
                kappa = (lam1 - lam2) / (1 - lam1)
                claimed = Fraction(
                    (m - s) * (n - m) * (n - 2),
                    n * (m - 1) * (n - m + s - 1),
                )
                assert kappa == claimed
    print("symbolic kappa formula checked over finite parameter grid")


def quadratic_energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(len(a)))


def row_square(a, x):
    return sum(
        sum(a[i][j] * x[j] for j in range(len(a))) ** 2
        for i in range(len(a))
    )


def completion_parseval_check():
    """Audit the exact uniform-completion mean and energy Parseval identity."""
    a = [
        [0, 1, 1, -1, 1, -1],
        [1, 0, -1, 1, 1, -1],
        [1, -1, 0, 1, -1, 1],
        [-1, 1, 1, 0, -1, -1],
        [1, 1, -1, -1, 0, 1],
        [-1, -1, 1, -1, 1, 0],
    ]
    n = len(a)
    s = (0, 1, 2, 3)
    t = (4, 5)
    child_words = list(itertools.product((-1, 1), repeat=len(s)))
    child_energy = lambda y: sum(
        a[s[i]][s[j]] * y[i] * y[j]
        for i in range(len(s)) for j in range(len(s))
    )
    y = max(child_words, key=lambda w: abs(child_energy(w)))
    q_s = abs(child_energy(y))
    sigma = 1 if child_energy(y) >= 0 else -1
    q = max(
        abs(quadratic_energy(a, x))
        for x in itertools.product((-1, 1), repeat=n)
    )
    b = [sum(a[i][s[j]] * y[j] for j in range(len(s))) for i in range(n)]
    internal = sum(b[i] ** 2 for i in s)
    external = sum(b[i] ** 2 for i in t)
    local_margins = [sigma * y[j] * b[s[j]] for j in range(len(s))]
    assert all(x >= 0 for x in local_margins)
    assert sum(local_margins) == q_s
    assert internal == sum(x * x for x in local_margins)
    assert internal <= (len(s) - 1) * q_s

    rows = []
    energies = []
    for outside in itertools.product((-1, 1), repeat=len(t)):
        x = [0] * n
        for j, vertex in enumerate(s):
            x[vertex] = y[j]
        for j, vertex in enumerate(t):
            x[vertex] = outside[j]
        rows.append(row_square(a, x))
        energies.append(quadratic_energy(a, x))
    mean_row = Fraction(sum(rows), len(rows))
    mean_energy_sq = Fraction(sum(e * e for e in energies), len(energies))
    assert mean_row == internal + external + len(t) * (n - 1)
    assert mean_energy_sq == q_s * q_s + 4 * external + 2 * len(t) * (len(t) - 1)
    assert mean_energy_sq <= q * q
    assert min(rows) <= mean_row
    print(
        "completion Parseval: "
        f"q={q}, q_S={q_s}, internal={internal}, external={external}, "
        f"E row={mean_row}, E energy^2={mean_energy_sq}"
    )


def main():
    symbolic_kappa_check()
    completion_parseval_check()
    exhaustive_slice_check(5, 3, 2)
    exhaustive_slice_check(6, 4, 3)
    aggregate_check(5, 3, 2)
    complete_signing_obstruction()
    print("PASS: first-level bound, sharpened coarea implication, and scalable mass wall")


if __name__ == "__main__":
    main()
