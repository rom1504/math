#!/usr/bin/env python3
"""Exact switching check for Theorem 17.1i's invariant rational grid."""

from fractions import Fraction as Q
from itertools import product


SHIFT = Q(1, 3)


def left_clamp(x):
    return max(Q(0), x - SHIFT)


def right_clamp(x):
    return min(Q(1), x + SHIFT)


MAPS = (left_clamp, right_clamp)


def follow(x, word):
    for letter in word:
        x = MAPS[letter](x)
    return x


def nearest_grid(x, denominator):
    scaled = x * denominator
    lower = scaled.numerator // scaled.denominator
    upper = min(lower + 1, denominator)
    candidates = (Q(lower, denominator), Q(upper, denominator))
    return min(candidates, key=lambda value: (abs(value - x), value))


def main():
    checks = 0
    maximum_error = Q(0)
    for refinement in range(1, 13):
        denominator = 3 * refinement
        grid = tuple(Q(i, denominator) for i in range(denominator + 1))

        # Exact invariance under both switching maps.
        for x in grid:
            for map_ in MAPS:
                assert map_(x) in grid
                checks += 1

        # Use half-grid test seeds so most starts are not themselves states.
        raw_seeds = tuple(Q(2 * i + 1, 2 * denominator)
                          for i in range(denominator))
        for raw in raw_seeds:
            center = nearest_grid(raw, denominator)
            initial_error = abs(raw - center)
            for depth in range(9):
                for word in product((0, 1), repeat=depth):
                    error = abs(follow(raw, word) - follow(center, word))
                    assert error <= initial_error
                    maximum_error = max(maximum_error, error)
                    checks += 1

    print(f"exact invariant-grid/switch-word checks: {checks}")
    print(f"largest shadow error: {maximum_error}")
    print("both generators have slope-one cells, so no strict contraction")


if __name__ == "__main__":
    main()
