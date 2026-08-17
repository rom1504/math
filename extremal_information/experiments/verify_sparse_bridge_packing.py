#!/usr/bin/env python3
"""Exact finite checks of the matching-bridge packing construction."""

from __future__ import annotations

from itertools import product
import json


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def hamming(x, y):
    return sum(a != b for a, b in zip(x, y))


def greedy_code(n, distance):
    remaining = list(product((-1, 1), repeat=n))
    code = []
    while remaining:
        c = remaining[0]
        code.append(c)
        remaining = [x for x in remaining if hamming(c, x) >= distance]
    return code


def response(code, labels, amplitude, query):
    vals = [amplitude * labels[x] + dot(x, query) for x in code]
    best = max(vals)
    return best, vals.count(best)


def verify():
    total_queries = 0
    minimum_margin = None
    rows = []
    for n in range(3, 9):
        distance = max(1, n // 3)
        code = greedy_code(n, distance)
        amplitude = distance
        # It suffices to test the two hostile label patterns separately for
        # each query: intended state 0 / every competitor 1, and vice versa.
        for c in code:
            hostile = {x: int(x != c) for x in code}
            value, multiplicity = response(code, hostile, amplitude, c)
            intended = n
            assert value == intended and multiplicity == 1
            competitor = max(
                (amplitude + dot(x, c) for x in code if x != c),
                default=-10**9,
            )
            margin = intended - competitor if competitor > -10**8 else None
            if margin is not None:
                minimum_margin = margin if minimum_margin is None else min(
                    minimum_margin, margin
                )
            friendly = {x: int(x == c) for x in code}
            value, multiplicity = response(code, friendly, amplitude, c)
            assert value == n + amplitude and multiplicity == 1
            total_queries += 2
        rows.append({"n": n, "distance": distance, "code_size": len(code)})
    return {
        "orders": rows,
        "hostile_and_friendly_queries_checked": total_queries,
        "minimum_observed_hostile_margin": minimum_margin,
        "arithmetic": "exact integer",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
