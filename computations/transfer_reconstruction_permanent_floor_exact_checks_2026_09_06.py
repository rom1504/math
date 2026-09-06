"""Exact finite regressions for the universal full-row permanent floor.

At t=log(2) and integral coordinates, the folded Gaussian kernel is rational.
The permanent, deletion comparison, and Sinkhorn normalization are exact;
Gaussian/log comparisons use outward rational intervals.
"""

from fractions import Fraction as F
from math import factorial
import json
import random

from continued_feedback_conditional_variance_exact_certificate_2026_09_06 import gaussian_upper
from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import log_interval


def permanent(matrix):
    order = len(matrix)
    values = [F(0)] * (1 << order)
    values[0] = F(1)
    for mask in range(1 << order):
        row = bin(mask).count("1")
        if row == order:
            continue
        for col in range(order):
            if not (mask >> col) & 1:
                values[mask | (1 << col)] += values[mask] * matrix[row][col]
    return values[-1]


def folded_kernel(a, b):
    return (F(1, 2 ** ((a - b) ** 2)) + F(1, 2 ** ((a + b) ** 2))) / 2


def gaussian_floor_cases(rng):
    log2_low, _ = log_interval(2)
    checked = deletion_checked = 0
    least_gap = None
    for order in range(1, 9):
        vectors = [[0] * order, [1] * order]
        vectors.extend([[rng.randrange(-2, 3) for _ in range(order)] for _ in range(8)])
        for vector in vectors:
            kernel = [[folded_kernel(a, b) for b in vector] for a in vector]
            p = permanent(kernel) / factorial(order)
            variance = F(sum(value * value for value in vector), order)
            gaussian_up = gaussian_upper(variance, log2_low)
            permanent_log_low, _ = log_interval(p)
            gap = permanent_log_low - 2 * order * gaussian_up
            assert gap >= 0
            if variance:
                least_gap = gap if least_gap is None else min(least_gap, gap)
            for removed in range(order):
                minor = [[value for j, value in enumerate(row) if j != removed]
                         for i, row in enumerate(kernel) if i != removed]
                p_minor = permanent(minor) / factorial(order - 1)
                assert p <= kernel[removed][removed] * p_minor <= p_minor
                deletion_checked += 1
            checked += 1
    return {"vectors": checked, "coordinate_deletions": deletion_checked,
            "least_nonzero_log_floor_gap_lower": str(least_gap)}


def scaling_cases(rng):
    checked = 0
    for order in range(1, 9):
        for _ in range(5):
            # J/order mixed with symmetrized permutation matrices is positive
            # and exactly doubly stochastic, including repeated entries.
            d = [[F(1, 2 * order) for _ in range(order)] for _ in range(order)]
            for _ in range(3):
                perm = rng.sample(range(order), order)
                for i, j in enumerate(perm):
                    d[i][j] += F(1, 12)
                    d[j][i] += F(1, 12)
            assert all(sum(row) == 1 for row in d)
            assert all(sum(d[i][j] for i in range(order)) == 1 for j in range(order))
            factors = [F(rng.randrange(1, 7), rng.randrange(1, 7)) for _ in range(order)]
            k = [[order * d[i][j] / (factors[i] * factors[j]) for j in range(order)]
                 for i in range(order)]
            product_factor = F(1)
            for value in factors:
                product_factor *= value * value
            lhs = permanent(k) / factorial(order)
            scaled = permanent(d) * F(order**order, factorial(order)) / product_factor
            assert lhs == scaled
            assert lhs >= 1 / product_factor
            checked += 1
    return checked


def threshold_checks():
    log2_low, _ = log_interval(2)
    for index in range(1, 128):
        p = F(index, 128)
        lower = p * log2_low + log_interval(1 - F(15, 16) * p)[0] / 4
        assert lower > 0
    return {"interior_retention_points": 127,
            "threshold_squared": str(F(15, 64)),
            "endpoint_equalities": "p=0 and p=1 hold symbolically"}


def main():
    rng = random.Random(2026090603)
    print(json.dumps({
        "folded_gaussian_floor": gaussian_floor_cases(rng),
        "sinkhorn_normalization_cases": scaling_cases(rng),
        "cap_threshold_checks": threshold_checks(),
        "arithmetic": "exact rationals and outward rational intervals only",
        "scope": "finite regressions; universal result uses van der Waerden and Gaussian self-cost maximum",
    }, indent=2))


if __name__ == "__main__":
    main()
