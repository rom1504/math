#!/usr/bin/env python3
"""Exact finite checks for the Wave 14 cross-mosaic/entropy audit."""

from collections import Counter
from itertools import combinations, product
import math

import numpy as np


A9 = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int16)


def spins(n):
    return np.array([(1,) + z for z in product((-1, 1), repeat=n - 1)],
                    dtype=np.int16)


def energies(a):
    x = spins(len(a))
    return np.sum((x @ a) * x, axis=1)


def q(a):
    return int(np.max(np.abs(energies(a))))


def min_fixing_certificate(a):
    """Minimum number of fixed spins making z^T a z constant.

    The quadratic is multilinear on the cube.  After fixing F, constancy is
    equivalent to: no edge has both endpoints free, and every free vertex has
    zero signed field from F.
    """
    n = len(a)
    for k in range(n + 1):
        for fixed_tuple in combinations(range(n), k):
            fixed = set(fixed_tuple)
            free = [i for i in range(n) if i not in fixed]
            if any(a[i, j] != 0
                   for p, i in enumerate(free) for j in free[p + 1:]):
                continue
            for values in product((-1, 1), repeat=k):
                assignment = dict(zip(fixed_tuple, values))
                if all(sum(int(a[i, j]) * assignment[j] for j in fixed) == 0
                       for i in free):
                    return k, fixed_tuple, values
    raise AssertionError("fixing all variables must work")


def all_signings(n):
    edges = tuple(combinations(range(n), 2))
    for signs in product((-1, 1), repeat=len(edges)):
        a = np.zeros((n, n), dtype=np.int16)
        for (i, j), sign in zip(edges, signs):
            a[i, j] = a[j, i] = sign
        yield a


def minimizers(n):
    best = n * (n - 1)
    out = []
    for a in all_signings(n):
        value = q(a)
        if value < best:
            best, out = value, [a]
        elif value == best:
            out.append(a)
    return best, out


def random_cross(sizes, seed):
    rng = np.random.default_rng(seed)
    n = sum(sizes)
    a = np.zeros((n, n), dtype=np.int16)
    starts = np.cumsum((0,) + tuple(sizes))
    for b in range(len(sizes)):
        for c in range(b + 1, len(sizes)):
            for i in range(starts[b], starts[b + 1]):
                for j in range(starts[c], starts[c + 1]):
                    a[i, j] = a[j, i] = rng.choice((-1, 1))
    return a


def profile_row(a):
    es = energies(a)
    cap = int(np.max(np.abs(es)))
    return (cap, len(set(map(int, es))), min_fixing_certificate(a)[0],
            tuple(int(np.sum(np.abs(es) >= cap - gap))
                  for gap in (0, 4, 8)))


def check_a9():
    left, right = tuple(range(3)), tuple(range(3, 9))
    c = A9.copy()
    c[np.ix_(left, left)] = 0
    c[np.ix_(right, right)] = 0
    expected_histogram = Counter({
        -24: 1, -20: 4, -16: 8, -12: 16, -8: 31, -4: 44,
        0: 48,
        4: 44, 8: 31, 12: 16, 16: 8, 20: 4, 24: 1,
    })
    assert Counter(map(int, energies(c))) == expected_histogram
    assert q(c) == 24
    assert tuple(int(np.sum(np.abs(energies(c)) >= 24 - gap))
                 for gap in (0, 4, 8, 12, 16, 20)) == (2, 10, 26, 58, 120, 208)
    fixing = min_fixing_certificate(c)
    assert fixing[0] == 7
    assert fixing == (7, (0, 3, 4, 5, 6, 7, 8),
                      (-1, -1, -1, 1, -1, 1, -1))

    q3, gs3 = minimizers(3)
    q6, gs6 = minimizers(6)
    assert (q3, len(gs3), q6, len(gs6)) == (6, 8, 10, 384)
    b0 = q3 + q6
    distribution = Counter()
    for g3 in gs3:
        for g6 in gs6:
            h = c.copy()
            h[np.ix_(left, left)] = g3
            h[np.ix_(right, right)] = g6
            value = q(h)
            # The new two-sided comparison, checked on every local minimizer.
            assert abs(value - q(c)) <= b0
            distribution[value] += 1
    assert sum(distribution.values()) == 3072
    assert min(distribution) == 28

    # The Hoeffding sum in (10.549): t=4 fails badly, while t=B0 makes the
    # low-margin set empty.  These are diagnostics, not needed by the theorem.
    v0 = q3 * q3 + q6 * q6
    hoeffding_t4 = (2 * math.exp(-4**2 / (2 * v0))
                    + 8 * math.exp(-8**2 / (2 * v0))
                    + 16 * math.exp(-12**2 / (2 * v0)))
    assert hoeffding_t4 > 17
    return c, distribution, hoeffding_t4


def check_random_profiles():
    expected = {
        ((2, 2, 2, 2), 0): (24, 6, 6, (3, 3, 29)),
        ((2, 2, 2, 2), 1): (36, 8, 6, (1, 1, 2)),
        ((2, 2, 2, 2), 2): (24, 6, 6, (4, 4, 28)),
        ((3, 3, 3), 0): (30, 8, 6, (3, 5, 13)),
        ((3, 3, 3), 1): (26, 7, 7, (4, 10, 32)),
        ((3, 3, 3), 2): (34, 9, 7, (1, 2, 6)),
        ((3, 3, 4), 0): (42, 19, 8, (1, 1, 5)),
        ((3, 3, 4), 1): (38, 18, 7, (1, 7, 11)),
        ((3, 3, 4), 2): (30, 15, 8, (1, 14, 41)),
    }
    actual = {key: profile_row(random_cross(key[0], key[1])) for key in expected}
    assert actual == expected
    # Complete multipartite support alone proves fixing number >= n-s.
    assert all(row[2] >= sum(sizes) - max(sizes)
               for (sizes, _), row in actual.items())
    return actual


def check_random_sign_constant():
    # (m-1)(m+2) <= (9/8)m^2, with equality at m=4.
    for m in range(1, 1001):
        assert 8 * (m - 1) * (m + 2) <= 9 * m * m
    # The ledger's F(m) is the undoubled Hamiltonian optimum M_m, whereas
    # q_m=min Q=2M_m.  Thus the random-sign bound carries a factor two.
    kappa = 3 * math.sqrt(math.log(2)) / math.sqrt(2)
    assert 1.766 < kappa < 1.767
    # This normalization now safely dominates the exact q_3=6.
    q3_upper = 2 * math.sqrt(3 * 2 * 5 * math.log(2))
    assert q3_upper > 6
    return kappa


if __name__ == "__main__":
    cross, distribution, hoeffding_t4 = check_a9()
    random_rows = check_random_profiles()
    kappa = check_random_sign_constant()
    print("A9 cross Q / level count / exact robust dependence:",
          q(cross), len(set(map(int, energies(cross)))),
          min_fixing_certificate(cross)[0])
    print("A9 hybrid norm distribution:", sorted(distribution.items()))
    print("A9 Hoeffding sum at t=4:", hoeffding_t4)
    print("random rows ((sizes), seed): (Q, levels, fixing, cap counts gaps 0/4/8)")
    for key, row in random_rows.items():
        print(key, row)
    print("finite random-sign constant kappa:", kappa)
    print("PASS")
