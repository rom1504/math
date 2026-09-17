"""Exact parity coupling and numerical corroboration of its sharp benchmark."""

import itertools
import json
from fractions import Fraction

import numpy as np
from scipy.special import logsumexp


def main():
    k = 5
    cube = list(itertools.product((-1, 1), repeat=k))
    code = [x for x in cube if np.prod(x) == 1]
    marginal = {x: Fraction(0) for x in cube}
    for x in code:
        conditional = {x: Fraction(1, 2)}
        for i in range(k):
            y = list(x)
            y[i] *= -1
            conditional[tuple(y)] = Fraction(1, 2 * k)
        for i in range(k):
            mean = sum(prob * y[i] for y, prob in conditional.items())
            assert mean == Fraction(k - 1, k) * x[i]
        for y, prob in conditional.items():
            marginal[y] += prob / len(code)
    assert all(p == Fraction(1, 2 ** k) for p in marginal.values())

    support = lambda y: max(sum(a * b for a, b in zip(x, y)) for x in code)
    left = sum(Fraction(support(x), len(code)) for x in code)
    right = sum(Fraction(support(y), len(cube)) for y in cube)
    assert left == k and right == k - 1
    assert left / right == Fraction(5, 4)

    # Exact convex-Taylor polynomial certificate used in the proof.
    for numerator in range(1001):
        u = Fraction(numerator, 1000)
        p = 5 * u ** 4 + 5 * u ** 2 - 12 * u + 5
        lower = (Fraction(17, 81) + Fraction(16, 27) * (u - Fraction(2, 3))
                 + 5 * (u - Fraction(2, 3)) ** 2)
        assert p >= lower >= Fraction(701, 3645)

    rng = np.random.default_rng(20260917)
    code_array = np.array(code, dtype=float)
    worst_gap = 0.0
    for _ in range(10000):
        theta = rng.normal(size=k) * 10 ** rng.uniform(-3, 1)
        log_mgf = logsumexp(code_array @ theta) - np.log(len(code))
        gap = float(log_mgf - np.dot(theta, theta) / 2)
        worst_gap = max(worst_gap, gap)
        assert gap <= 2e-12
    print(json.dumps({"status": "PASS", "dimension": k,
                      "code_size": len(code), "cube_size": len(cube),
                      "sharp_dilation": "5/4", "rational_polynomial_checks": 1001,
                      "mgf_directions_checked": 10000,
                      "largest_numerical_mgf_excess": worst_gap}, indent=2))


if __name__ == "__main__":
    main()
