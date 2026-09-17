#!/usr/bin/env python3
"""Exact finite sanity checks for the Wave 51 reverse-tail inequality."""

from fractions import Fraction
from itertools import combinations, product
from random import Random


def check_polynomial(k, linear, quadratic):
    vals = []
    for w in product((-1, 1), repeat=k):
        z = sum(linear[i] * w[i] for i in range(k))
        z += sum(quadratic[i, j] * w[i] * w[j]
                 for i, j in combinations(range(k), 2))
        vals.append(z)
    count = len(vals)
    assert sum(vals) == 0
    v = Fraction(sum(z * z for z in vals), count)
    if not v:
        return 0
    fourth = Fraction(sum(z ** 4 for z in vals), count)
    # Rational form of ||Z||_4 <= 3 ||Z||_2.
    assert fourth <= 81 * v * v
    # Check P{Z<=0} >= 1/324 exactly.
    nonpositive = sum(z <= 0 for z in vals)
    assert 324 * nonpositive >= count
    # Squared form of (A) for every attained nonnegative integer threshold r
    # below sqrt(V)/18.  Avoid irrational arithmetic: the claimed RHS is at
    # most 1/324, and direct floating evaluation is only a supplemental test.
    vf = float(v)
    for r in range(max(0, max(-z for z in vals)) + 1):
        rhs = max(0.0, 1 / 18 - r / (vf ** 0.5)) ** 2
        lhs = sum(z <= -r for z in vals) / count
        assert lhs + 1e-15 >= rhs
    return 1


def main():
    rng = Random(510051)
    tested = 0
    for k in range(1, 8):
        pairs = list(combinations(range(k), 2))
        for _ in range(250):
            linear = [rng.randrange(-4, 5) for _ in range(k)]
            quadratic = {ij: rng.randrange(-4, 5) for ij in pairs}
            tested += check_polynomial(k, linear, quadratic)
    print(f"checked {tested} nonzero centered degree-two Rademacher polynomials")
    print("all exact fourth-moment and finite reverse-tail checks passed")


if __name__ == "__main__":
    main()
