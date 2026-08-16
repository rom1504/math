#!/usr/bin/env python3
"""Exact finite falsifiers for the normalized Max-Cut shell theorem.

The unrestricted projective lookup compiler is checked in
``verify_maxcut_projective_response.py``.  This script independently checks
the anisotropic Lipschitz condition, the outer max-plus shell identity, and
the coordinate-oscillation realization metric.
"""

from __future__ import annotations

from itertools import product
import json
import random


SpinWord = tuple[int, ...]


def canonical(x: SpinWord) -> SpinWord:
    negative = tuple(-a for a in x)
    return min(x, negative)


def projective_words(w: int) -> list[SpinWord]:
    return sorted({canonical(x) for x in product((-1, 1), repeat=w)})


def projective_distance(x: SpinWord, y: SpinWord, weights: tuple[int, ...]) -> int:
    direct = sum(c for a, b, c in zip(x, y, weights) if a != b)
    complement = sum(c for a, b, c in zip(x, y, weights) if a == b)
    return min(direct, complement)


def shell_response(
    f: dict[SpinWord, int], weights: tuple[int, ...], outer: SpinWord
) -> int:
    """Return the shell response with its boundary-independent constant removed."""
    return max(
        f[canonical(inner)]
        - sum(c for a, b, c in zip(outer, inner, weights) if a != b)
        for inner in product((-1, 1), repeat=len(weights))
    )


def coordinate_oscillations(f: dict[SpinWord, int], w: int) -> tuple[int, ...]:
    deltas = []
    for i in range(w):
        delta = 0
        for spin in product((-1, 1), repeat=w):
            flipped = list(spin)
            flipped[i] *= -1
            delta = max(
                delta,
                abs(f[canonical(spin)] - f[canonical(tuple(flipped))]),
            )
        deltas.append(delta)
    return tuple(deltas)


def check_trial(w: int, rng: random.Random) -> int:
    words = projective_words(w)
    weights = tuple(rng.randint(0, 4) for _ in range(w))

    # A maximum of weighted distance cones is automatically d_weights-Lipschitz.
    seeds = [(rng.choice(words), rng.randint(-7, 7)) for _ in range(8)]
    f = {
        x: max(bias - projective_distance(x, centre, weights) for centre, bias in seeds)
        for x in words
    }

    inequalities = 0
    for x, y in product(words, repeat=2):
        assert abs(f[x] - f[y]) <= projective_distance(x, y, weights)
        inequalities += 1

    for outer in product((-1, 1), repeat=w):
        assert shell_response(f, weights, outer) == f[canonical(outer)]

    deltas = coordinate_oscillations(f, w)
    for x, y in product(words, repeat=2):
        assert abs(f[x] - f[y]) <= projective_distance(x, y, deltas)

    return inequalities


def main() -> None:
    rng = random.Random(20260816)
    trials = 0
    inequalities = 0
    for w in range(2, 7):
        for _ in range(200):
            inequalities += check_trial(w, rng)
            trials += 1

    report = {
        "all_checks_passed": True,
        "seed": 20260816,
        "widths": [2, 3, 4, 5, 6],
        "trials": trials,
        "anisotropic_lipschitz_inequalities": inequalities,
        "shell_identity": True,
        "coordinate_oscillation_metric": True,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
