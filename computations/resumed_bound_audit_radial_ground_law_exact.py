#!/usr/bin/env python3
"""Exact finite falsifier of a radial common ground-state law; no LP."""

from fractions import Fraction
from itertools import combinations, product


def main():
    edges = list(combinations(range(4), 2))
    a = [(-1 if edge == (2, 3) else 1) for edge in edges]
    spins = [(1,) + rest for rest in product((-1, 1), repeat=3)]

    def features(x):
        return [x[i] * x[j] for i, j in edges]

    energies = [sum(ai * vi for ai, vi in zip(a, features(x))) for x in spins]
    assert max(map(abs, energies)) == 4
    assert sum(energy * energy for energy in energies) == 6 * len(spins)
    for signs in product((-1, 1), repeat=6):
        cap = max(abs(sum(ai * vi for ai, vi in zip(signs, features(x))))
                  for x in spins)
        assert cap >= 4

    for x in spins:
        separator = x[0]*x[1] + x[0]*x[3] + x[1]*x[2] - x[2]*x[3]
        assert abs(separator) == 2

    law = [(1, (1, 1, -1, 1)), (1, (1, 1, 1, -1)),
           (1, (1, 1, 1, 1)), (-1, (1, 1, -1, -1))]
    mean = [sum(Fraction(sigma * features(x)[e], 4) for sigma, x in law)
            for e in range(6)]
    assert mean == [Fraction(ai, 2) for ai in a]
    mean_energy = sum(ai * vi for ai, vi in zip(a, mean))
    assert mean_energy == 3
    assert Fraction(4, 6) > Fraction(1, 2)
    print("PASS: M_4=4, maximal radial coefficient=1/2, minimum mean gap=1")


if __name__ == "__main__":
    main()
