"""Rational moment checks for the recursive tensor-Hadamard ground law."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json

import numpy as np


PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def moment(a, indices):
    odd_indices = tuple(sorted(i for i, count in Counter(indices).items() if count % 2))
    return reduced_moment(a, odd_indices)


@lru_cache(None)
def reduced_moment(a, indices):
    if not indices:
        return Fraction(1)
    if len(indices) % 2:
        return Fraction(0)
    assert a > 0
    block_size = 4**(a-1)
    blocks = [[i % block_size for i in indices if i//block_size == b]
              for b in range(4)]
    answer = Fraction(0)
    for first, second in PAIRINGS:
        term = Fraction(1)
        for left, right in (first, second):
            term *= (-1)**len(blocks[right])
            term *= moment(a-1, tuple(blocks[left]+blocks[right]))
        answer += term/3
    return answer


def axis_line(a, indices):
    varying = 0
    for digit in range(a):
        values = {(i//4**digit) % 4 for i in indices}
        if len(values) > 1:
            if len(values) != 4:
                return False
            varying += 1
    return varying == 1


def sample(a, rng):
    if a == 0:
        return np.array([rng.choice([-1, 1])], dtype=np.int64)
    first, second = PAIRINGS[int(rng.integers(3))]
    blocks = [None]*4
    for left, right in (first, second):
        child = sample(a-1, rng)
        blocks[left], blocks[right] = child, -child
    return np.concatenate(blocks)


def main():
    bound = Fraction(11, 27)
    report = []
    for a in (1, 2):
        n = 4**a
        maximum = Fraction(0)
        exceptions = 0
        for indices in combinations(range(n), 4):
            value = abs(moment(a, indices))
            if axis_line(a, indices):
                assert value == 1
                exceptions += 1
            else:
                assert value <= bound
                maximum = max(maximum, value)
        assert exceptions == a*n//4
        for i, j in combinations(range(n), 2):
            differing_digits = sum((i//4**d) % 4 != (j//4**d) % 4 for d in range(a))
            assert moment(a, (i, j)) == Fraction(-1, 3)**differing_digits
        report.append({"n": n, "axis_line_exceptions": exceptions,
                       "largest_other_fourth_moment": str(maximum)})
    rng = np.random.default_rng(202609173)
    a, n = 3, 64
    h4 = np.ones((4, 4), dtype=np.int64)-2*np.eye(4, dtype=np.int64)
    h = np.kron(np.kron(h4, h4), h4)
    for _ in range(20):
        x = sample(a, rng)
        assert np.array_equal(h @ x, (-2)**a*x)
        permutation = rng.permutation(n)
        pairs = list(zip(permutation[::2], permutation[1::2]))
        signs = rng.choice([-1, 1], size=n//2)
        p = len(pairs)
        second_moment = Fraction(p)
        for i, j in combinations(range(p), 2):
            second_moment += 2*int(signs[i]*signs[j])*moment(a, pairs[i]+pairs[j])
        assert second_moment <= bound*p*p+2*(1-bound)*p
    print(json.dumps({"status": "all rational moment and integer eigenvector checks passed",
                      "exhaustive_fourth_moment_orders": report,
                      "random_signed_matchings_checked_at_n64": 20}, indent=2))


if __name__ == "__main__":
    main()
