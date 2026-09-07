#!/usr/bin/env python3
"""Exact finite checks for the Wave 55 joint-completion identities.

This checks only the algebra and the projective/full-cube distribution
identification.  The Hanson--Wright inequality and asymptotic comparison are
proved analytically in the ledger.
"""

from collections import Counter
from itertools import combinations, product

import numpy as np


def sign_matrix(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    upper = rng.choice((-1, 1), size=(n, n))
    upper = np.triu(upper, 1)
    return upper + upper.T


def quadratic(a: np.ndarray, x: tuple[int, ...]) -> int:
    v = np.asarray(x, dtype=int)
    return int(v @ a @ v)


checks = 0
for n in range(4, 9):
    a = sign_matrix(n, 1000 + n)
    aop = np.linalg.norm(a, ord=2)
    for m in range(2, n):
        # A few selectors of each size suffice because all assertions below
        # are exact identities, not a probabilistic test.
        for s_tuple in list(combinations(range(n), m))[: min(4, n)]:
            s = list(s_tuple)
            t = [i for i in range(n) if i not in s]
            ps = np.zeros((n, n), dtype=int)
            ps[np.ix_(s, s)] = a[np.ix_(s, s)]
            cmat = a - ps

            assert np.trace(cmat) == 0
            expected_frob_sq = 2 * m * (n - m) + (n - m) * (n - m - 1)
            assert int(np.sum(cmat * cmat)) == expected_frob_sq
            assert expected_frob_sq <= n * n
            assert np.linalg.norm(cmat, ord=2) <= 2 * aop + 1e-10

            # Gauge the projective local word at the first coordinate of S.
            anchor = s[0]
            free_s = s[1:]
            projective_values: list[int] = []
            for local_free in product((-1, 1), repeat=len(free_s)):
                y = {anchor: 1, **dict(zip(free_s, local_free))}
                for outside in product((-1, 1), repeat=len(t)):
                    x = tuple(y[i] if i in y else outside[t.index(i)] for i in range(n))
                    z_block = quadratic(cmat, x)
                    yv = np.asarray([y[i] for i in s], dtype=int)
                    wv = np.asarray(outside, dtype=int)
                    cross = 2 * int(yv @ a[np.ix_(s, t)] @ wv)
                    tail = int(wv @ a[np.ix_(t, t)] @ wv)
                    assert z_block == cross + tail
                    projective_values.append(z_block)

            full_values = [quadratic(cmat, x) for x in product((-1, 1), repeat=n)]
            # Every projective value occurs twice on the full cube because
            # x and -x have the same quadratic value.
            assert Counter(full_values) == Counter(projective_values * 2)
            checks += 1

print(f"joint completion identities verified on {checks} matrix/selector cases")
