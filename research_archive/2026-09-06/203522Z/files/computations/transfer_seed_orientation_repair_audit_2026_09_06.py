"""Exact finite replay of the director's orientation repair identities.

The probabilistic remainder bound and order-limit argument are proved in
the companion audit, not inferred from random finite trials.
"""

import itertools
from functools import lru_cache

import numpy as np


@lru_cache(None)
def spins(n):
    return np.array([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int64)


def extrema(a):
    x = spins(len(a))
    doubled = np.einsum("bi,ij,bj->b", x, a, x)
    assert not np.any(doubled % 2)
    return int(np.max(doubled)) // 2, int(-np.min(doubled)) // 2


rng = np.random.default_rng(2609062015)
cases = 0
for n in range(2, 8):
    for trial in range(12):
        upper = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
        a = upper + upper.T
        p, r = extrema(a)
        if r > p:
            a = -a
            p, r = r, p
        delta = p - r
        k = 0
        while k * (k - 1) // 2 < delta:
            k += 1
        assert delta == 0 and k == 0 or 0 <= k * (k - 1) // 2 - delta < k
        extra = trial % 2
        total = n + k + extra
        w = np.zeros((total, total), dtype=np.int64)
        w[:n, :n] = a
        w[n:n + k, n:n + k] = -1
        np.fill_diagonal(w, 0)
        wp, wr = extrema(w)
        assert wp == p + k // 2
        assert wr == r + k * (k - 1) // 2
        assert p <= wp <= p + k and p <= wr <= p + k
        remainder = np.zeros_like(w)
        edges = 0
        for i, j in itertools.combinations(range(total), 2):
            if w[i, j] == 0:
                remainder[i, j] = remainder[j, i] = rng.choice((-1, 1))
                edges += 1
        assert edges <= total * (k + extra)
        vp, vr = extrema(remainder)
        vcap = max(vp, vr)
        bp, br = extrema(w + remainder)
        assert p <= max(bp, br) <= p + k + vcap
        assert abs(bp - br) <= k + 2 * vcap
        assert abs(bp - wp) <= vcap and abs(br - wr) <= vcap
        cases += 1

print(f"PASS: {cases} exact negative-clique, padded-remainder and orientation checks.")
