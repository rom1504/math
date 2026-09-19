#!/usr/bin/env python3
"""Exact support matching for the two-field stability kernel.

An elementary dominance-matching exchange algorithm, not a permanent
approximation or an asymptotic signing theorem. No novelty claim.
"""
import bisect
import itertools
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def dominance_matching(requests, columns):
    """Return a bijection with columns[j]>=requests[i] coordinatewise."""
    n = len(requests)
    assert len(columns) == n
    if not n:
        return []
    levels = sorted(set(c[1] for c in columns))
    buckets = [[] for _ in levels]
    tree = [0] * (len(levels) + 1)

    def add(i, value):
        i += 1
        while i < len(tree):
            tree[i] += value
            i += i & -i

    def prefix(count):
        value = 0
        while count:
            value += tree[count]
            count -= count & -count
        return value

    def kth(rank):
        index = 0
        bit = 1 << (len(levels).bit_length() - 1)
        while bit:
            candidate = index + bit
            if candidate < len(tree) and tree[candidate] < rank:
                rank -= tree[candidate]
                index = candidate
            bit >>= 1
        return index

    row_order = sorted(range(n), key=lambda i: requests[i][0], reverse=True)
    col_order = sorted(range(n), key=lambda j: columns[j][0], reverse=True)
    answer = [-1] * n
    pointer = 0
    for i in row_order:
        while pointer < n and columns[col_order[pointer]][0] >= requests[i][0]:
            j = col_order[pointer]
            level = bisect.bisect_left(levels, columns[j][1])
            buckets[level].append(j)
            add(level, 1)
            pointer += 1
        lower = bisect.bisect_left(levels, requests[i][1])
        before = prefix(lower)
        if before == prefix(len(levels)):
            return None
        level = kth(before + 1)
        answer[i] = buckets[level].pop()
        add(level, -1)
    return answer


def brute_feasible(requests, columns):
    return any(all(columns[p[i]][0] >= a and columns[p[i]][1] >= b
                   for i, (a, b) in enumerate(requests))
               for p in itertools.permutations(range(len(requests))))


def check(requests, columns):
    answer = dominance_matching(requests, columns)
    assert (answer is not None) == brute_feasible(requests, columns)
    if answer is not None:
        assert sorted(answer) == list(range(len(requests)))
        assert all(columns[j][0] >= requests[i][0] and
                   columns[j][1] >= requests[i][1] for i, j in enumerate(answer))


def main():
    rng = random.Random(2026091919)
    generic = 0
    for n in range(9):
        for _ in range(80):
            requests = [(rng.randrange(-3, 4), rng.randrange(-3, 4)) for _ in range(n)]
            columns = [(rng.randrange(-3, 4), rng.randrange(-3, 4)) for _ in range(n)]
            check(requests, columns)
            generic += 1
    actual = 0
    for n in range(2, 8):
        for _ in range(30):
            a = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i):
                    a[i][j] = a[j][i] = rng.choice([-1, 1])
            x, y, z = [[rng.choice([-1, 1]) for _ in range(n)] for _ in range(3)]
            p = list(range(n))
            rng.shuffle(p)
            w = [0] * n
            for i in range(n):
                w[p[i]] = x[i] * y[i] * z[p[i]]
            def apply(v):
                return [sum(a[i][j] * v[j] for j in range(n)) for i in range(n)]
            ax, ay, az, aw = map(apply, (x, y, z, w))
            for sector in (-1, 1):
                rows = [i for i in range(n) if x[i] * y[i] == sector]
                cols = [j for j in range(n) if z[j] * w[j] == sector]
                for threshold in (0, 2):
                    requests = [(threshold - x[i] * ax[i], threshold + y[i] * ay[i]) for i in rows]
                    columns = [(z[j] * aw[j], w[j] * az[j]) for j in cols]
                    check(requests, columns)
                    actual += 1
    result = dict(status="exact finite validation, not an asymptotic certificate",
                  seed=2026091919, generic_cases=generic,
                  actual_child_sector_threshold_cases=actual, all_pass=True,
                  complexity="O(n log n) comparisons/arithmetic via sorted thresholds and Fenwick counts")
    path = ROOT / "computations/results/twisted_chiral_dominance_matching_2026_09_19.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
