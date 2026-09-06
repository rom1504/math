#!/usr/bin/env python3
"""Exact small-order search for boundary-compatible minimizing replacements."""

from functools import lru_cache
from itertools import combinations, product
import numpy as np


@lru_cache(None)
def spins(n):
    # Projective cube: first spin fixed to +1.
    return np.array([(1,) + z for z in product((-1, 1), repeat=n-1)], dtype=np.int16)


def gauge_signing(n, mask):
    a = np.zeros((n, n), dtype=np.int16)
    a[0, 1:] = a[1:, 0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i+1, n):
            a[i, j] = a[j, i] = 1 if (mask >> k) & 1 else -1
            k += 1
    return a


def energies(a):
    x = spins(len(a))
    return np.sum((x @ a) * x, axis=1)


def q(a):
    return int(np.max(np.abs(energies(a))))


@lru_cache(None)
def gauge_minimizers(n):
    ne = (n-1)*(n-2)//2
    best = n*(n-1)
    out = []
    for mask in range(1 << ne):
        a = gauge_signing(n, mask)
        qa = q(a)
        if qa < best:
            best, out = qa, [a]
        elif qa == best:
            out.append(a)
    return best, out


@lru_cache(None)
def all_minimizers(n):
    """Every labeled signing in every minimizing switching class."""
    qn, reps = gauge_minimizers(n)
    unique = {}
    for a in reps:
        for tail in product((-1, 1), repeat=n-1):
            s = np.array((1,) + tail, dtype=np.int16)
            b = a * np.outer(s, s)
            unique[b.tobytes()] = b
    return qn, list(unique.values())


def layer(c, d):
    y = spins(len(d))
    exposure = 2*np.sum(np.abs(y @ c.T), axis=1)
    payoff = exposure - q(d) + np.abs(energies(d))
    return int(max(0, np.max(payoff)))


def best_replacement_layer(c, m):
    q_m, mins = all_minimizers(m)
    y = spins(m)
    exposure = 2*np.sum(np.abs(y @ c.T), axis=1)
    best = None
    witness = None
    for g in mins:
        val = int(max(0, np.max(exposure - q_m + np.abs(energies(g)))))
        if best is None or val < best:
            best, witness = val, g
    return best, witness


def normalized_c_key(c):
    # Row signs and row order do not affect 2||Cy||_1.
    z = c * c[:, :1]
    rows = tuple(sorted(tuple(int(v) for v in row) for row in z))
    return rows


@lru_cache(None)
def best_replacement_layer_key(rows):
    c = np.array(rows, dtype=np.int16)
    return best_replacement_layer(c, c.shape[1])


def scan(n, min_tail=2, tail_sizes=None):
    qn, parents = gauge_minimizers(n)
    print(f'n={n}, q={qn}, gauge minimizers={len(parents)}')
    worst_gap = -10**9
    cert = None
    hist = {}
    for ai, a in enumerate(parents):
        verts = range(n)
        for m in (range(min_tail, n) if tail_sizes is None else tail_sizes):
            for tail_tuple in combinations(verts, m):
                tail = list(tail_tuple)
                head = [i for i in verts if i not in tail]
                d = a[np.ix_(tail, tail)]
                c = a[np.ix_(head, tail)]
                old = layer(c, d)
                new, g = best_replacement_layer_key(normalized_c_key(c))
                gap = new-old
                hist[gap] = hist.get(gap, 0)+1
                if gap > worst_gap:
                    worst_gap = gap
                    cert = (ai, tuple(head), tuple(tail), a.copy(), c.copy(), d.copy(), old, new, g.copy(), q(d))
    print('gap histogram', sorted(hist.items()))
    ai, head, tail, a, c, d, old, new, g, qd = cert
    print('worst gap',worst_gap,'parent',ai,'head',head,'tail',tail,'Q(D)',qd,'old/new',old,new)
    print('A=',a.tolist())
    print('C=',c.tolist())
    print('D=',d.tolist())
    print('G=',g.tolist())
    return cert


if __name__ == '__main__':
    for n in range(4,8):
        scan(n)
