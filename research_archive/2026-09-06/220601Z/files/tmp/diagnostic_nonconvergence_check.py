#!/usr/bin/env python3
"""Independent checks for diagnostic_nonconvergence_report.md.

This script checks only algebraic/numerical consequences that are used in the
report.  It does not claim to decide convergence or nonconvergence.
"""

from __future__ import annotations

import math
from functools import reduce
from math import gcd

import numpy as np


LOG2 = math.log(2.0)
C0 = math.sqrt(9.0 * LOG2 / 8.0)


def random_cap_bound(k: int) -> float:
    """The union-bound upper bound R_k for M_k."""
    if k <= 1:
        return 0.0
    return math.sqrt(k * (k - 1) * (k + 2) * LOG2)


def finite_relative_modulus(n: int, m: int) -> tuple[float, float, float]:
    """Bounds on a_m-a_n and a_n-a_m for 1 <= n <= m.

    Returns (upward_bound, downward_bound, max_bound), where a_k=M_k/k^(3/2).
    """
    assert 1 <= n <= m
    if n == m:
        return (0.0, 0.0, 0.0)
    h = m - n
    beta = h / m
    upward = (
        C0 * beta**1.5
        + math.sqrt(2.0 * LOG2 * beta * (1.0 - beta) * (1.0 + 2.0 / m))
    )
    downward = C0 * (1.0 - (1.0 - beta) ** 1.5)
    return upward, downward, max(upward, downward)


def delta_for_epsilon(eps: float) -> float:
    """A conservative relative-width parameter proved in the report."""
    assert eps > 0
    return min(
        0.5,
        (eps / (2.0 * C0)) ** (2.0 / 3.0),
        eps * eps / (16.0 * LOG2),
        2.0 * eps / (3.0 * C0),
    )


def augmented_cut_weight_gcd(n: int) -> int:
    N = n * (n - 1) // 2
    weights = set()
    for s in range(n + 1):
        w = s * (n - s)
        weights.add(w)
        weights.add(N - w)
    nonzero = sorted(w for w in weights if w)
    return reduce(gcd, nonzero)


def legendre(a: int, q: int) -> int:
    a %= q
    if a == 0:
        return 0
    z = pow(a, (q - 1) // 2, q)
    return 1 if z == 1 else -1


def paley_conference(q: int) -> np.ndarray:
    """Symmetric Paley conference matrix of order q+1, q prime and q=1 mod 4."""
    assert q % 4 == 1
    N = q + 1
    C = np.zeros((N, N), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for i in range(q):
        for j in range(q):
            if i != j:
                C[i + 1, j + 1] = legendre(i - j, q)
    return C


def check_random_bound_constant() -> None:
    ratios = []
    for k in range(2, 10000):
        ratios.append((random_cap_bound(k) / (k**1.5), k))
    value, argmax = max(ratios)
    assert argmax == 4
    assert abs(value - C0) < 1e-14
    print(f"random-bound constant C0={C0:.12f}, attained at k={argmax}")


def check_delta_rule() -> None:
    for eps in (0.01, 0.025, 0.05, 0.10):
        delta = delta_for_epsilon(eps)
        # beta <= delta is the hypothesis.  Check a fine grid and several m.
        worst = 0.0
        for m in (2, 3, 5, 10, 100, 10000):
            for j in range(1001):
                beta = delta * j / 1000.0
                up = C0 * beta**1.5 + math.sqrt(
                    2.0 * LOG2 * beta * (1.0 - beta) * (1.0 + 2.0 / m)
                )
                down = C0 * (1.0 - (1.0 - beta) ** 1.5)
                worst = max(worst, up, down)
        assert worst <= eps * (1.0 + 1e-12)
        print(
            f"epsilon={eps:.3f}: delta={delta:.6g}, "
            f"multiplicative radius=1/(1-delta)={1/(1-delta):.9f}, "
            f"checked worst={worst:.6g}"
        )


def check_code_gcd() -> None:
    for n in range(3, 101):
        expected = 2 if n % 4 == 1 else 1
        actual = augmented_cut_weight_gcd(n)
        assert actual == expected, (n, actual, expected)
    print("augmented cut-code weight gcd formula checked for 3 <= n <= 100")


def check_conference_minors() -> None:
    # All listed q are primes congruent to 1 mod 4.
    for q in (5, 13, 17, 29, 37, 41, 53, 61):
        C = paley_conference(q)
        N = q + 1
        assert np.array_equal(C @ C, q * np.eye(N, dtype=np.int64))
        # Delete a deterministic sublinear-sized tail (finite examples only).
        r = max(1, int(math.sqrt(N) // 2))
        n = N - r
        A = C[:n, :n]
        B = C[:n, n:]
        defect = A @ A - (n - 1) * np.eye(n, dtype=np.int64)
        assert np.array_equal(defect, r * np.eye(n, dtype=np.int64) - B @ B.T)
        defect_f = float(np.linalg.norm(defect, ord="fro"))
        certified = r * math.sqrt(n) + math.sqrt(n * r * (N - 1))
        op = float(np.linalg.norm(A, ord=2))
        assert defect_f <= certified * (1.0 + 1e-12)
        assert op <= math.sqrt(N - 1) * (1.0 + 1e-12)
        print(
            f"Paley N={N:2d}, minor n={n:2d}, r={r}: "
            f"op/sqrt(N-1)={op/math.sqrt(N-1):.9f}, "
            f"defect/n^(3/2)={defect_f/(n**1.5):.9f}, "
            f"certified={certified/(n**1.5):.9f}"
        )


def main() -> None:
    check_random_bound_constant()
    check_delta_rule()
    check_code_gcd()
    check_conference_minors()
    print("all diagnostic_nonconvergence checks passed")


if __name__ == "__main__":
    main()
