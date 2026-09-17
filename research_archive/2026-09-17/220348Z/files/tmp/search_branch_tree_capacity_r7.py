#!/usr/bin/env python3
"""Explore total one-orientation-per-edge capacity on endpoint-cut trees."""

from functools import lru_cache
from random import Random

from allocation_harvest_r7 import (
    endpoint_data,
    matrix_from_code,
    split_record,
    code_from_matrix,
)


@lru_cache(None)
def branch_capacity(n, code):
    if n <= 1:
        return 0, None
    P, N, ps, ns = endpoint_data(n, code)
    if P + N == 0:
        return 0, None
    A = matrix_from_code(n, code)
    best = (-1, None)
    for p in ps:
        for neg in ns:
            rec = split_record(A, p, neg)
            if rec is None:
                continue
            rows = []
            total = 0
            for shore in rec:
                local = max(0, *shore["mu"])
                child, _ = branch_capacity(len(shore["x"]), shore["code"])
                total += local + child
                rows.append((shore["x"], shore["mu"], local, child))
            if total > best[0]:
                best = (total, (p, neg, rows))
    return best


def random_code(n, rng):
    return rng.randrange(1 << (n * (n - 1) // 2))


def main():
    rng = Random(71060)
    for n in range(2, 12):
        trials = 1 << (n * (n - 1) // 2) if n <= 5 else 3000
        worst = (-1,)
        for _ in range(trials):
            code = _ if n <= 5 else random_code(n, rng)
            P, N, _, _ = endpoint_data(n, code)
            G, witness = branch_capacity(n, code)
            ratio = G / max(P, N)
            if ratio > worst[0]:
                worst = (ratio, G, P, N, code, witness)
        print("n", n, "max G/Q", worst[:5], flush=True)


if __name__ == "__main__":
    main()
