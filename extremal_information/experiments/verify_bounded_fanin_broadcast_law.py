#!/usr/bin/env python3
"""Finite wind tunnel for the bounded-fan-in broadcast law."""

from __future__ import annotations

import itertools
import math
import random


def q_value(edges, spin):
    return sum(value * spin[u] * spin[v] for u, v, value in edges)


def exact_q(n, edges):
    return max(
        abs(q_value(edges, spin))
        for spin in itertools.product((-1, 1), repeat=n)
    )


def check_parity_cells():
    n, g, t, d = 8, 2, 3, 8
    all_edges = list(itertools.combinations(range(n), 2))
    cells = [all_edges[j * d:(j + 1) * d] for j in range(g)]
    h = g * t
    for mask in range(1 << h):
        signs = []
        for j in range(g):
            parity = sum((mask >> (j * t + i)) & 1 for i in range(t)) & 1
            signs.append(-1 if parity else 1)
        for bit in range(h):
            j = bit // t
            changed = mask ^ (1 << bit)
            changed_parity = sum(
                (changed >> (j * t + i)) & 1 for i in range(t)
            ) & 1
            changed_sign = -1 if changed_parity else 1
            difference = [
                (u, v, changed_sign - signs[j]) for u, v in cells[j]
            ]
            assert exact_q(n, difference) == 2 * d
    assert h * 2 * d == 2 * t * g * d


def random_support_checks():
    rng = random.Random(20260817)
    for n in range(4, 10):
        universe = list(itertools.combinations(range(n), 2))
        for _ in range(50):
            sample = universe[:]
            rng.shuffle(sample)
            m = rng.randint(1, len(sample))
            edges = [(u, v, rng.choice((-2, 2))) for u, v in sample[:m]]
            assert exact_q(n, edges) + 1e-12 >= m / math.sqrt(2 * n)


def pairs(r):
    return [(i, j) for i in range(r) for j in range(i + 1, r)]


def wedge_mask(x, y, r):
    out = 0
    for bit, (i, j) in enumerate(pairs(r)):
        value = ((((x >> i) & 1) & ((y >> j) & 1))
                 ^ (((x >> j) & 1) & ((y >> i) & 1)))
        out |= value << bit
    return out


def fixed_sampler(r):
    rng = random.Random(10_000 * r + 400)
    return [rng.randrange(1 << r) for _ in range(4 * r * r)]


def gram_fanin_checks():
    reports = []
    for r in range(2, 6):
        labels = fixed_sampler(r)
        rows = [
            wedge_mask(labels[i], labels[j], r)
            for i in range(len(labels))
            for j in range(i + 1, len(labels))
        ]
        h = len(pairs(r))
        e = len(rows)
        degrees = [sum((row >> bit) & 1 for row in rows) for bit in range(h)]
        weights = [bin(row).count("1") for row in rows]
        assert min(degrees) >= e / 4
        assert sum(degrees) == sum(weights)
        assert max(weights) >= h / 4
        reports.append((r, h, e, min(degrees), sum(weights) / e, max(weights)))
    return reports


def local_switching_check():
    """Find a small switching code/cell and check its coherent witness."""
    rng = random.Random(17082026)
    n, t, d = 10, 2, 12
    words = []
    while len(words) < (1 << t):
        word = tuple(rng.choice((-1, 1)) for _ in range(n))
        if all(
            n // 4 <= sum(a != b for a, b in zip(word, old)) <= 3 * n // 4
            for old in words
        ):
            words.append(word)
    universe = list(itertools.combinations(range(n), 2))
    for _ in range(10000):
        cell = rng.sample(universe, d)
        if all(
            sum((a[u] * b[u]) != (a[v] * b[v]) for u, v in cell) >= 2
            for a, b in itertools.combinations(words, 2)
        ):
            break
    else:
        raise AssertionError("small switching cell not found")
    for a, b in itertools.combinations(words, 2):
        difference = []
        changed = 0
        for u, v in cell:
            old = a[u] * a[v]
            new = b[u] * b[v]
            if old != new:
                changed += 1
                difference.append((u, v, new - old))
        assert abs(q_value(difference, a)) == 2 * changed
        assert exact_q(n, difference) == 2 * changed


def main():
    check_parity_cells()
    random_support_checks()
    local_switching_check()
    reports = gram_fanin_checks()
    print("bounded-fan-in broadcast checks passed")
    for report in reports:
        print("Gram fan-in r,h,E,min-degree,average,max =", report)


if __name__ == "__main__":
    main()
