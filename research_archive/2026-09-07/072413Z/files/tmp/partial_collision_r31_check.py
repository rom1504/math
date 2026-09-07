#!/usr/bin/env python3
"""Exhaustive finite checks of the partial-collision sandwich (P1)."""

from fractions import Fraction
from itertools import combinations, product


def event_probability(sets, weights, r, s):
    total = Fraction(0)
    universe = range(len(weights))
    for sample in product(universe, repeat=r):
        if max(sum(x in hit for x in sample) for hit in sets) >= s:
            mass = Fraction(1)
            for x in sample:
                mass *= weights[x]
            total += mass
    return total


def main():
    omega = tuple(range(4))
    all_sets = [frozenset(c) for q in range(5) for c in combinations(omega, q)]
    laws = [
        (Fraction(1, 4),) * 4,
        (Fraction(1, 2), Fraction(1, 3), Fraction(1, 6), Fraction(0)),
        (Fraction(7, 10), Fraction(1, 10), Fraction(1, 10), Fraction(1, 10)),
    ]
    checked = 0
    for sets in combinations(all_sets, 3):
        for weights in laws:
            masses = [sum(weights[x] for x in hit) for hit in sets]
            h = max(masses)
            for r in range(1, 5):
                for s in range(1, r + 1):
                    prob = event_probability(sets, weights, r, s)
                    choose = Fraction(len(tuple(combinations(range(r), s))))
                    assert h**s <= prob
                    assert prob <= len(sets) * choose * h**s
                    checked += 1
    print("PARTIAL_COLLISION_CASES", checked)
    print("PASS partial_collision_r31")


if __name__ == "__main__":
    main()
