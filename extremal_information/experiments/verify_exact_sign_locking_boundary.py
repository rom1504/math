#!/usr/bin/env python3
"""Finite checks for the sparse contrast identity and dense locking bound."""

from itertools import combinations, product
from math import comb


def spins(n):
    return list(product((-1, 1), repeat=n))


def check_sparse(k):
    xs = spins(k)
    edges = list(combinations(range(k), 2))
    ys = spins(len(edges))
    for tbits in product((-1, 1), repeat=len(edges)):
        # A deliberately nonquadratic test child makes the identity genuinely
        # pointwise rather than a cancellation special case.
        def f(x):
            return 3 * x[0] - 2 * x[-1] + int(sum(x) == k)

        lhs = max(
            f(x)
            + sum(y[e] * (x[i] - tbits[e] * x[j])
                  for e, (i, j) in enumerate(edges))
            for x in xs for y in ys
        )
        rhs = len(edges) + max(
            f(x) - sum(tbits[e] * x[i] * x[j]
                       for e, (i, j) in enumerate(edges))
            for x in xs
        )
        assert lhs == rhs


def check_locking(k):
    xs = spins(k)
    # Exhaust all bridges only through k=3; sample a deterministic rich slice
    # at k=4 so the script stays a wind tunnel rather than a proof substitute.
    matrices = product((-1, 1), repeat=k * k)
    cap = 2 ** (k * k) if k <= 3 else 512
    for idx, flat in enumerate(matrices):
        if idx >= cap:
            break
        R = [flat[i * k:(i + 1) * k] for i in range(k)]
        for d in spins(k):
            defects = []
            for x in xs:
                fields = [sum(R[i][j] * x[i] for i in range(k))
                          for j in range(k)]
                free = sum(abs(v) for v in fields)
                locked = sum(x[i] * R[i][j] * d[j] * x[j]
                             for i in range(k) for j in range(k))
                defects.append(free - locked)
            assert min(defects) >= 0
            # Exact average formula, stronger than merely checking the bound.
            mu_num = sum(abs(sum(x)) for x in xs)
            expected_sum = k * mu_num
            trace_term = len(xs) * sum(R[i][i] * d[i] for i in range(k))
            assert sum(defects) == expected_sum - trace_term
            assert max(defects) * len(xs) >= sum(defects)


def main():
    for k in (2, 3):
        check_sparse(k)
    for k in (2, 3, 4):
        check_locking(k)
    print("exact-sign locking-boundary checks passed")


if __name__ == "__main__":
    main()
