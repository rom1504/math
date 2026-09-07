"""Exact finite normalization checks; no numerical asymptotic claim."""
from collections import defaultdict
from itertools import combinations, product
from math import comb
from random import Random


def inner(a, b):
    return sum(v * b.get(k, 0) for k, v in a.items())


def check(c):
    m = len(c)
    n = 2 * m
    a = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            a[i][m + j] = a[m + j][i] = c[i][j]
    shore = [list(range(m)), list(range(m, n))]
    # V_i=(sqrt(2)/m)*v_i: own-spin, root degree two.
    # X_i=(sqrt(2)/m**1.5)*x_i: added external degree-one root.
    v = []
    x = []
    raw = []
    for i in range(n):
        opposite = shore[1 if i < m else 0]
        vi = {}
        for j, k in combinations(opposite, 2):
            vi[tuple(sorted((i, j, k)))] = a[i][j] * a[i][k]
        v.append(vi)
        xi = defaultdict(int)
        ri = defaultdict(int)
        own = shore[0 if i < m else 1]
        for center in opposite:
            for j, k in combinations(own, 2):
                key = tuple(sorted((center, j, k)))
                value = a[i][center] * a[center][j] * a[center][k]
                ri[key] += value
                if i not in (j, k):
                    xi[key] += value
        x.append(dict(xi))
        raw.append(dict(ri))
    for i in range(n):
        for j in range(n):
            assert inner(v[i], v[j]) == (comb(m, 2) if i == j else 0)
            aa = sum(a[i][z] * a[j][z] for z in range(n))
            if i == j:
                target = m * comb(m - 1, 2)
            elif (i < m) == (j < m):
                target = (comb(m - 2, 2) if m >= 4 else 0) * aa
            else:
                target = 0
            assert inner(x[i], x[j]) == target
            # The unexcluded outer product is literally A V.
            direct = defaultdict(int)
            for z in range(n):
                if a[i][z]:
                    for key, val in v[z].items():
                        direct[key] += a[i][z] * val
            assert dict(direct) == raw[i]
    beta = max(sum(abs(sum(c[i][j] * y[j] for j in range(m)))
                   for i in range(m))
               for y in product((-1, 1), repeat=m))
    high = -10 ** 9
    low = 10 ** 9
    for z in product((-1, 1), repeat=n):
        energy = sum(a[i][j] * z[i] * z[j]
                     for i in range(n) for j in range(i + 1, n))
        high = max(high, energy)
        low = min(low, energy)
    assert high == beta and low == -beta
    return n * n


def main():
    rng = Random(7309021)
    count = 0
    grams = 0
    for m in range(2, 7):
        samples = 16 if m == 2 else 12
        for sample in range(samples):
            if m == 2:
                bits = [(sample >> q) & 1 for q in range(m * m)]
                c = [[1 - 2 * bits[i * m + j] for j in range(m)]
                     for i in range(m)]
            else:
                c = [[rng.choice((-1, 1)) for _ in range(m)] for _ in range(m)]
            grams += check(c)
            count += 1
    print("PASS:", count, "bipartite matrices;", grams,
          "exact marked/output Gram entries;", count, "exact cap normalizations")


if __name__ == "__main__":
    main()
