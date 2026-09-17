#!/usr/bin/env python3
"""Exact binomial audit of the root's critical-size block counterexample.

Cell membership, normalization, and covariance are checked with rational
arithmetic. Constants involving sqrt(n) and the limiting normal law are
reported numerically, not as exact transcendental certificates.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from scipy.stats import norm


def exact_binomial(n: int) -> dict:
    counts = [0, 0, 0]
    first = [0, 0, 0]
    second = [0, 0, 0]
    total_first = 0
    binomial = 1
    for k in range(n + 1):
        s = abs(n - 2 * k)
        squared = 100 * s * s
        cell = None
        if squared <= n:
            cell = 0
        elif 81 * n <= squared <= 121 * n:
            cell = 1
        elif 361 * n <= squared <= 441 * n:
            cell = 2
        if cell is not None:
            counts[cell] += binomial
            first[cell] += binomial * s
            second[cell] += binomial * s * s
        total_first += binomial * s
        if k < n:
            binomial = binomial * (n - k) // (k + 1)
    if not all(counts):
        return {"n": n, "status": "at least one cell empty"}
    probabilities = [Fraction(c, 1 << n) for c in counts]
    a = [Fraction(v, n * c) for v, c in zip(second, counts)]
    b_times_root_n = [Fraction(v, c) for v, c in zip(first, counts)]
    delta = probabilities[1] / (2 * (a[2] - a[0]))
    f = [(a[2] - a[1]) / probabilities[0],
         -(a[2] - a[0]) / probabilities[1],
         (a[1] - a[0]) / probabilities[2]]
    weights = [1 + delta * v for v in f]
    normalization_error = sum(p * (w - 1) for p, w in zip(probabilities, weights))
    second_moment_error = sum(p * ai * (w - 1) for p, ai, w in zip(probabilities, a, weights))
    assert normalization_error == 0
    assert second_moment_error == 0
    assert weights[1] == Fraction(1, 2)
    assert min(weights) > 0
    first_change_times_root_n = sum(p * bi * (w - 1)
                                   for p, bi, w in zip(probabilities, b_times_root_n, weights))
    change = float(first_change_times_root_n) / math.sqrt(n)
    interval_bound = float(delta) * (-0.468)
    assert change <= interval_bound + 1e-12
    unbiased = float(Fraction(total_first, 1 << n)) / math.sqrt(n)
    weighted = unbiased + change
    return {"n": n, "status": "exact rational identities passed",
            "cell_probabilities": [float(v) for v in probabilities],
            "cell_conditional_second_moments": [float(v) for v in a],
            "cell_conditional_absolute_first_moments": [float(v) / math.sqrt(n) for v in b_times_root_n],
            "delta": float(delta), "density_weights": [float(v) for v in weights],
            "density_upper_bound": max(1.0, *(float(v) for v in weights)),
            "unbiased_normalized_absolute_mean": unbiased,
            "weighted_normalized_absolute_mean": weighted,
            "mean_change": change, "interval_certified_upper_bound_on_change": interval_bound,
            "gaussian_normalized_absolute_mean": math.sqrt(2 / math.pi),
            "weighted_minus_gaussian": weighted - math.sqrt(2 / math.pi),
            "exact_normalization_error": str(normalization_error),
            "exact_second_moment_error": str(second_moment_error)}


def normal_limit() -> dict:
    intervals = [(0.0, 0.1), (0.9, 1.1), (1.9, 2.1)]
    p, a, b = [], [], []
    for low, high in intervals:
        probability = 2 * (norm.cdf(high) - norm.cdf(low))
        first = 2 * (norm.pdf(low) - norm.pdf(high))
        second = probability + 2 * (low * norm.pdf(low) - high * norm.pdf(high))
        p.append(probability)
        a.append(second / probability)
        b.append(first / probability)
    delta = p[1] / (2 * (a[2] - a[0]))
    f = [(a[2] - a[1]) / p[0], -(a[2] - a[0]) / p[1], (a[1] - a[0]) / p[2]]
    weights = [1 + delta * v for v in f]
    change = delta * ((a[2] - a[1]) * b[0] - (a[2] - a[0]) * b[1] + (a[1] - a[0]) * b[2])
    return {"status": "numerical normal integrals, not interval certified",
            "cell_probabilities": p, "a": a, "b": b, "delta": delta,
            "density_weights": weights, "absolute_mean_change": change,
            "weighted_absolute_mean": math.sqrt(2 / math.pi) + change,
            "interval_proof_guaranteed_decrease": 0.468 * delta}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=[100, 400, 1000, 10000])
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/critical_block_audit.json"))
    args = parser.parse_args()
    result = {"exact_finite_audits": [exact_binomial(n) for n in args.orders],
              "normal_limit": normal_limit()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
