#!/usr/bin/env python3
"""Exact small checks for the fixed-rank upper-roof composition algebra."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import json


def response(states, t):
    """states are scalar-feature pairs (u,h), all exact Fractions."""
    return max(h + t * u for u, h in states)


def star(left, right):
    return [(u + v, h + k + u * v) for u, h in left for v, k in right]


def direct_response(left, right, t):
    return max(h + k + u * v + t * (u + v)
               for u, h in left for v, k in right)


def rank_one_landscape(n):
    scale = 2**n - 1
    states = []
    for x in product((-1, 1), repeat=n):
        raw = sum((2**i) * x[i] for i in range(n))
        u = Fraction(raw, scale)
        states.append((u, -(u * u)))
    return sorted(states)


def verify():
    landscapes = [
        [(Fraction(-1), Fraction(0)), (Fraction(1), Fraction(1))],
        [(Fraction(-2), Fraction(1)), (Fraction(0), Fraction(3)),
         (Fraction(2), Fraction(-1))],
        rank_one_landscape(3),
    ]
    fields = [Fraction(j, 3) for j in range(-9, 10)]
    binary_checks = 0
    associative_checks = 0
    for left in landscapes:
        for right in landscapes:
            composed = star(left, right)
            for t in fields:
                assert response(composed, t) == direct_response(left, right, t)
                binary_checks += 1
            for third in landscapes:
                lhs = star(star(left, right), third)
                rhs = star(left, star(right, third))
                for t in fields:
                    assert response(lhs, t) == response(rhs, t)
                    associative_checks += 1

    exposure_checks = 0
    for n in range(1, 9):
        states = rank_one_landscape(n)
        assert len({u for u, _ in states}) == 2**n
        for u, h in states:
            vals = [hh + 2 * u * uu for uu, hh in states]
            winner = max(range(len(vals)), key=vals.__getitem__)
            assert states[winner] == (u, h)
            assert vals.count(vals[winner]) == 1
            exposure_checks += 1

    return {
        "binary_response_identities": binary_checks,
        "three_way_associativity_responses": associative_checks,
        "rank_one_unique_exposures": exposure_checks,
        "max_rank_one_order": 8,
        "arithmetic": "exact rational",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
