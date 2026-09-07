#!/usr/bin/env python3
"""Finite regressions for the fixed-seed mixed-orbit proof.

No asymptotic claim is inferred from these tests.  Exact integer/fraction
tests check the fourth moment and K5 cycle generation; finite Gaussian
kernel tests independently check reflection-twisted graph contraction.
"""
import itertools
import json
import math
from fractions import Fraction

import numpy as np


def hadamard(q):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < q:
        h = np.block([[h, h], [h, -h]])
    assert len(h) == q
    return h


def fourth_coefficient(row):
    coeff = [1, 0, 0, 0, 0]
    for x in row:
        for j in range(4, 0, -1):
            coeff[j] += int(x) * coeff[j - 1]
    return coeff[4]


def fourth_check(q, b):
    h = hadamard(q)
    signed_rows = list(h) + list(-h)
    total = 0
    count = 0
    for leaves in itertools.product(signed_rows, repeat=b):
        total += fourth_coefficient(np.concatenate(leaves))
        count += 1
    m = q * b
    observed = Fraction(total, count * math.comb(m, 4))
    target = Fraction(b * q * (q - 1) * (q - 2),
                      m * (m - 1) * (m - 2) * (m - 3))
    assert observed == target
    return {"q": q, "b": b, "row_words": count,
            "fourth_moment": str(observed), "k5_moment": str(observed**5),
            "status": "PASS"}


def unsigned_p(v, t):
    values = []
    for perm in itertools.permutations(range(len(v))):
        w = v[list(perm)]
        values.append(math.exp(-t * float(np.sum((v - w)**2))))
    return sum(values) / len(values)


def twisted_graph_check():
    rows = [np.array([-1.0, 2.0]), np.array([-3.0, 0.5]),
            np.array([1.0, -0.5])]
    edges = [(0, 1), (0, 2), (1, 2)]
    neighbors = [[1, 2], [0, 2], [0, 1]]
    results = []
    for t in [0.1, 0.5, 1.0, 3.0]:
        upper = math.prod(math.sqrt(unsigned_p(v, t)) for v in rows)
        ratios = []
        for signs in itertools.product([-1, 1], repeat=3):
            weights = []
            for flips in itertools.product([False, True], repeat=3):
                word = [v[::-1] if flip else v
                        for v, flip in zip(rows, flips)]
                energy = 0.0
                for (i, j), sign in zip(edges, signs):
                    a = word[i][neighbors[i].index(j)]
                    b = word[j][neighbors[j].index(i)]
                    energy += (a - sign * b)**2
                weights.append(math.exp(-t * energy))
            average = sum(weights) / len(weights)
            assert average <= upper + 1e-14
            ratios.append(average / upper)
        results.append({"t": t, "signings": 8,
                        "largest_ratio_to_product_norm": max(ratios),
                        "status": "PASS"})
    return results


def k5_edges(vertices):
    return {tuple(sorted(e)) for e in itertools.combinations(vertices, 2)}


def k5_span_check(m):
    count = 0
    for a, b, c, d in itertools.combinations(range(m), 4):
        rest = sorted(set(range(m)) - {a, b, c, d})[:3]
        assert len(rest) == 3
        result = set()
        for x, y in [(a, c), (b, c), (a, d), (b, d)]:
            result ^= k5_edges(rest + [x, y])
        expected = {tuple(sorted(e)) for e in
                    [(a, c), (b, c), (a, d), (b, d)]}
        assert result == expected
        count += 1
    return {"m": m, "four_cycle_identities": count, "status": "PASS"}


def main():
    output = {
        "status": "PASS",
        "scope": "finite regressions, not an asymptotic certificate",
        "fourth_moments": [fourth_check(4, b) for b in [1, 2, 4]]
                          + [fourth_check(8, b) for b in [1, 2]],
        "twisted_graph": twisted_graph_check(),
        "k5_span": [k5_span_check(m) for m in [7, 8, 9]],
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
