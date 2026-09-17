"""Exact small-law checks for the Hadamard-basis sparse-response construction."""

import itertools
import json
import math
from fractions import Fraction

import numpy as np
from scipy.special import logsumexp


def hadamard(n):
    h = np.ones((1, 1), dtype=int)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    assert len(h) == n
    return h


def absolute_walk_mean(p):
    return sum((Fraction(math.comb(p, j) * abs(p - 2 * j), 2 ** p)
                for j in range(p + 1)), Fraction(0))


def finite_law(k, p, rng):
    hk, hp = hadamard(k), hadamard(p)
    signs = np.array(list(itertools.product((-1, 1), repeat=p)), dtype=int)
    columns = np.array([np.kron(row, g) for row in hk for g in signs])
    n = k * p
    centers = np.kron(hk, hp)
    assert np.all(np.isin(columns, (-1, 1)))
    assert np.all(columns.sum(axis=0) == 0)
    assert np.array_equal(columns.T @ columns, len(columns) * np.eye(n, dtype=int))
    responses = centers @ columns.T
    target = absolute_walk_mean(p)
    for row in responses:
        assert Fraction(int(np.abs(row).sum()), len(columns)) == target
        assert sum(int(a) ** 2 for a in row) == n * len(columns)
    # Active probability is exactly 1/k before possible cancellation in S_p.
    for a in range(k):
        assert np.array_equal(hk @ hk[a], k * np.eye(k, dtype=int)[a])

    max_excess = 0.0
    for _ in range(1000):
        theta = rng.normal(size=n) * 10 ** rng.uniform(-2, 0.5)
        log_mgf = logsumexp(columns @ theta) - math.log(len(columns))
        excess = float(log_mgf - k * np.dot(theta, theta) / 2)
        max_excess = max(max_excess, excess)
        assert excess <= 1e-11
    return {"k": k, "p": p, "n": n, "law_atoms": len(columns),
            "center_absolute_mean": str(target),
            "normalized_center_mean": float(target) / math.sqrt(n),
            "covariance": "I exactly", "largest_mgf_excess": max_excess}


def all_order_law(k, p, leftover, rng):
    hk = hadamard(k)
    signs = list(itertools.product((-1, 1), repeat=p))
    extras = list(itertools.product((-1, 1), repeat=leftover))
    columns = np.array([np.concatenate((np.kron(row, g), extra))
                        for row in hk for g in signs for extra in extras])
    n = k * p + leftover
    assert leftover < k
    assert np.all(columns.sum(axis=0) == 0)
    assert np.array_equal(columns.T @ columns, len(columns) * np.eye(n, dtype=int))
    # Deliberately use an arbitrary, non-Hadamard short code.
    code = signs[:min(3, len(signs))]
    centers = np.array([np.concatenate((np.kron(row, g), extra))
                        for row in hk for g in code for extra in extras])
    bound = absolute_walk_mean(p) + leftover
    for values in centers @ columns.T:
        assert Fraction(int(np.abs(values).sum()), len(columns)) <= bound
    for _ in range(300):
        theta = rng.normal(size=n) * 10 ** rng.uniform(-2, 0.5)
        log_mgf = logsumexp(columns @ theta) - math.log(len(columns))
        assert log_mgf <= k * np.dot(theta, theta) / 2 + 1e-11
    return {"k": k, "p": p, "leftover": leftover, "n": n,
            "law_atoms": len(columns), "arbitrary_short_code_size": len(code),
            "center_absolute_mean_upper": str(bound), "covariance": "I exactly"}


def main():
    rng = np.random.default_rng(2026091702)
    cases = [finite_law(k, p, rng) for k, p in ((2, 2), (2, 4), (4, 4), (4, 8))]
    asymptotics = [{"p": p, "a_p_over_sqrt_p": float(absolute_walk_mean(p)) / math.sqrt(p)}
                   for p in (16, 64, 256, 1024)]
    # Verify the shell optimization and normalized limiting center constants.
    for k in (2, 4, 8):
        for kap in (0.01, 0.2, 1.0):
            eps, ratio = 0.003, 1.05
            eta_star = 2 * k * kap * ratio ** 2 * eps
            value = math.sqrt(8 * k * eps * kap * eta_star) - eta_star / ratio
            assert abs(value - 2 * k * kap * ratio * eps) < 1e-14
    c_lower = 0.4333221116640807
    all_orders = [all_order_law(k, p, ell, rng)
                  for k, p, ell in ((2, 3, 1), (4, 2, 3), (4, 3, 1))]
    print(json.dumps({"status": "PASS", "finite_laws": cases,
                      "all_order_arbitrary_code_checks": all_orders,
                      "walk_mean_asymptotics": asymptotics,
                      "gaussian_absolute_mean": math.sqrt(2 / math.pi),
                      "k2_center_slope": 1 / math.sqrt(math.pi),
                      "current_liminf_derivative_lower": 1.5 * c_lower,
                      "scope": "Cover theorem only; no minimizer-cover realization asserted."}, indent=2))


if __name__ == "__main__":
    main()
