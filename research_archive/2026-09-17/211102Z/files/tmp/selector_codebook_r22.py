#!/usr/bin/env python3
"""Exact finite certificates for the Wave 22 selector-codebook memo."""

from fractions import Fraction
import importlib.util
import itertools
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "puncture_cycles_r20", ROOT / "tmp" / "puncture_cycles_r20.py"
)
PC = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PC)


def score(A, t, x, omitted=None):
    omitted = set(() if omitted is None else omitted)
    return t * sum(
        A[i][j] * x[i] * x[j]
        for i in range(len(A))
        for j in range(i + 1, len(A))
        if i not in omitted and j not in omitted
    )


def coverage_patterns(A):
    """All nonempty exact-child-ground incidence patterns of oriented cuts."""
    n = len(A)
    child_q = [PC.matrix_norm(PC.principal(A, (i,))[0])[0] for i in range(n)]
    patterns = {}
    for t in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            pattern = tuple(i for i in range(n) if score(A, t, x, (i,)) == child_q[i])
            if pattern:
                patterns.setdefault(pattern, []).append((t, x, score(A, t, x)))
    return patterns


CASES = {
    "A6": (
        PC.A6,
        {
            tuple(j for j in range(6) if j != i): Fraction(1, 6)
            for i in range(6)
        },
        (32, 26, 5),
        [Fraction(5, 6)] * 6,
    ),
    "A8": (
        PC.A8,
        {
            pattern: Fraction(1, 8)
            for pattern in (
                (3, 4, 5),
                (0, 3, 4),
                (0, 4, 5),
                (0, 3, 5),
                (1, 2, 6),
                (1, 2, 7),
                (2, 6, 7),
                (1, 6, 7),
            )
        },
        (24, 20, 3),
        [Fraction(3, 8)] * 8,
    ),
    "A9": (
        PC.A9,
        {
            (3, 6, 8): Fraction(1, 10),
            (3, 5, 7): Fraction(1, 20),
            (0, 4, 8): Fraction(1, 4),
            (1, 3, 8): Fraction(1, 20),
            (2, 3, 5): Fraction(1, 5),
            (1, 2, 7): Fraction(1, 5),
            (0, 6, 7): Fraction(3, 20),
        },
        (25, 22, 3),
        [
            Fraction(2, 5),
            Fraction(1, 4),
            Fraction(2, 5),
            Fraction(2, 5),
            Fraction(1, 4),
            Fraction(1, 4),
            Fraction(1, 4),
            Fraction(2, 5),
            Fraction(2, 5),
        ],
    ),
}


def verify_case(name, A, mass, expected_counts, expected_z):
    n = len(A)
    q, _ = PC.matrix_norm(A)
    patterns = coverage_patterns(A)
    state_count = sum(len(states) for states in patterns.values())
    counts = (state_count, len(patterns), max(map(len, patterns)))
    assert counts == expected_counts
    assert sum(mass.values()) == 1
    assert set(mass).issubset(patterns)

    # Every support word may be chosen to be an exact parent ground.
    for pattern in mass:
        assert max(parent_score for _, _, parent_score in patterns[pattern]) == q

    z = [sum(weight for pattern, weight in mass.items() if i in pattern) for i in range(n)]
    assert z == expected_z
    assert min(z) > 0

    # KKT certificate for minimizing -(1/n) sum_i log z_i over the simplex
    # of all oriented-cut incidence patterns.
    kkt = {
        pattern: sum(Fraction(1, n) / z[i] for i in pattern)
        for pattern in patterns
    }
    assert max(kkt.values()) == 1
    assert all(kkt[pattern] == 1 for pattern in mass)

    rate = -sum(math.log(float(value)) for value in z) / n
    print(
        f"{name}: states={state_count}, patterns={len(patterns)}, "
        f"max_cover={counts[2]}, z={[str(value) for value in z]}, "
        f"R_U(0)={rate:.15f} nats"
    )


def main():
    for name, data in CASES.items():
        verify_case(name, *data)
    print("all selector-codebook certificates passed")


if __name__ == "__main__":
    main()
