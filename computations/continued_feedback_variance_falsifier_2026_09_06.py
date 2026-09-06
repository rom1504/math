#!/usr/bin/env python3
"""Exact replay of the scalable cubic-feedback variance-floor falsifier.

The seed certificate and all finite lifted diagonal entries use integers
and Fraction. Decimal values are only for display. No solver is used.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

import numpy as np


SEED = np.array(
    [
        [-1, -1, 1, 1, 1],
        [-1, 1, 1, -1, 1],
        [1, 1, -1, -1, 1],
        [1, -1, -1, 1, 1],
        [1, 1, 1, 1, -1],
    ],
    dtype=np.int64,
)


def sylvester(order: int) -> np.ndarray:
    assert order >= 1 and order & (order - 1) == 0
    h = np.ones((1, 1), dtype=np.int64)
    while h.shape[0] < order:
        h = np.block([[h, h], [h, -h]])
    assert np.array_equal(h @ h, order * np.eye(order, dtype=np.int64))
    return h


def variance_numerators(a: np.ndarray) -> np.ndarray:
    q = a @ a
    return np.diag(a @ (q**3) @ a)


def exact_certificate(max_lift: int) -> dict:
    q = SEED @ SEED
    t = SEED @ (q**3) @ SEED
    assert np.array_equal(np.diag(t), [733, 733, 733, 733, 533])
    assert Fraction(int(t[4, 4]), 625) == Fraction(533, 625) < 1
    lifts = []
    order = 2
    while order <= max_lift:
        a = np.kron(SEED, sylvester(order))
        np.fill_diagonal(a, 0)
        n = a.shape[0]
        assert np.array_equal(a, a.T)
        assert np.all(np.diag(a) == 0)
        assert np.all(np.abs(a[np.triu_indices(n, 1)]) == 1)
        numerator = variance_numerators(a)
        denominator = (n - 1) ** 4
        fifth = [Fraction(int(v), denominator) for v in numerator[4 * order :]]
        lifts.append(
            {
                "hadamard_order": order,
                "signing_order": n,
                "fifth_block_min": str(min(fifth)),
                "fifth_block_max": str(max(fifth)),
                "fifth_block_below_one_count": sum(v < 1 for v in fifth),
                "fifth_block_mean": str(sum(fifth, Fraction()) / order),
                "fifth_block_mean_decimal": float(sum(fifth, Fraction()) / order),
            }
        )
        order *= 2
    return {
        "status": "exact integer and rational certificate",
        "seed": SEED.tolist(),
        "seed_square": q.tolist(),
        "cubic_feedback_numerator": t.tolist(),
        "cubic_feedback_denominator": 625,
        "fifth_block_limit": "533/625",
        "lifts": lifts,
    }


def apex_certificate(max_lift: int) -> list[dict]:
    """Apex plus opposite twins: no positive uniform pointwise floor."""
    result = []
    m = 2
    twin = np.array([[1, -1], [-1, 1]], dtype=np.int64)
    while m <= max_lift:
        h = sylvester(m)
        assert int(np.trace(h)) == 0
        a = np.ones((2 * m + 1, 2 * m + 1), dtype=np.int64)
        a[1:, 1:] = np.kron(twin, h)
        np.fill_diagonal(a, 0)
        q = a @ a
        numerator = int(np.sum(q[1:, 1:] ** 3))
        actual = Fraction(numerator, (2 * m) ** 4)
        claimed = Fraction(36 * m * m - 40 * m + 1, 8 * m**3)
        assert actual == claimed
        result.append(
            {
                "hadamard_order": m,
                "signing_order": 2 * m + 1,
                "apex_variance": str(actual),
                "apex_variance_decimal": float(actual),
                "verified_formula": "(36*m*m-40*m+1)/(8*m**3)",
            }
        )
        m *= 2
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-lift", type=int, default=64)
    args = parser.parse_args()
    # At these orders signed int64 arithmetic is comfortably exact.
    if args.max_lift > 256:
        raise ValueError("Use max-lift <=256 to retain the audited int64 range.")
    certificate = exact_certificate(args.max_lift)
    certificate["apex_vanishing_variance"] = apex_certificate(args.max_lift)
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()
