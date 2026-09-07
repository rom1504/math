#!/usr/bin/env python3
"""Exact finite checks for the Wave 24 joint row-square/coverage LP."""

from collections import Counter
from fractions import Fraction
import importlib.util
import itertools
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "selector_codebook_r22", ROOT / "tmp" / "selector_codebook_r22.py"
)
SC = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SC)


def ground_types(A):
    """Return exact (uniform exact-child coverage, row-square) types."""
    n = len(A)
    parent_q = SC.PC.matrix_norm(A)[0]
    child_q = [
        SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
        for i in range(n)
    ]
    types = []
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            if SC.score(A, orientation, x) != parent_q:
                continue
            rows = tuple(
                orientation
                * x[i]
                * sum(A[i][j] * x[j] for j in range(n) if j != i)
                for i in range(n)
            )
            assert min(rows) >= 0
            assert sum(rows) == 2 * parent_q
            exact_children = sum(
                SC.score(A, orientation, x, (i,)) == child_q[i]
                for i in range(n)
            )
            types.append((Fraction(exact_children, n), sum(r * r for r in rows)))
    return Counter(types)


EXPECTED = {
    "A6": Counter({(Fraction(5, 6), 30): 12}),
    "A8": Counter({(Fraction(3, 8), 64): 8}),
    "A9": Counter(
        {
            (Fraction(1, 9), 80): 2,
            (Fraction(2, 9), 96): 12,
            (Fraction(1, 3), 112): 4,
            (Fraction(2, 9), 112): 3,
            (Fraction(1, 3), 128): 4,
        }
    ),
}


def dual_value(types, budget):
    """Compute min_theta max_g u_g[1+theta(C-c_g)] exactly."""
    actions = list(types) + [(Fraction(0), Fraction(0))]
    candidates = {Fraction(0)}
    lines = [
        (u, u * (budget - cost))
        for u, cost in actions
    ]
    for a, (intercept_a, slope_a) in enumerate(lines):
        for intercept_b, slope_b in lines[a + 1 :]:
            if slope_a == slope_b:
                continue
            theta = (intercept_b - intercept_a) / (slope_a - slope_b)
            if theta >= 0:
                candidates.add(theta)
    return min(max(intercept + theta * slope for intercept, slope in lines) for theta in candidates)


def a9_expected(budget):
    if budget < 80:
        return Fraction(0)
    if budget <= 96:
        inverse = Fraction(9) - Fraction(9, 32) * (budget - 80)
        return 1 / inverse
    if budget <= 112:
        inverse = Fraction(9, 2) - Fraction(3, 32) * (budget - 96)
        return 1 / inverse
    return Fraction(1, 3)


def main():
    cases = {name: data[0] for name, data in SC.CASES.items()}
    found = {name: ground_types(A) for name, A in cases.items()}
    assert found == EXPECTED

    assert dual_value(found["A6"], Fraction(29)) == 0
    assert dual_value(found["A6"], Fraction(30)) == Fraction(5, 6)
    assert dual_value(found["A8"], Fraction(63)) == 0
    assert dual_value(found["A8"], Fraction(64)) == Fraction(3, 8)

    for budget in map(Fraction, (79, 80, 88, 96, 104, 112, 120)):
        value = dual_value(found["A9"], budget)
        assert value == a9_expected(budget)
        print(f"A9 C={budget}: Z*={value}")

    print("joint regular/captured-row LP certificates passed")


if __name__ == "__main__":
    main()
