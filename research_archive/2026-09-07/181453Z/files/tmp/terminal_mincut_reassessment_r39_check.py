#!/usr/bin/env python3
"""Exact checks for the Wave 39 subtree-flip/cut-code reassessment."""

from __future__ import annotations

import itertools

import numpy as np


def qnorm(a: np.ndarray) -> int:
    n = len(a)
    return max(
        abs(int(x @ a @ x))
        for tail in itertools.product((-1, 1), repeat=n - 1)
        for x in (np.array((1,) + tail, dtype=np.int64),)
    )


def edge_sets(x, parents):
    return {
        (parents[v], v): frozenset(np.flatnonzero(x[parents[v]] != x[v]).tolist())
        for v in range(1, len(parents))
    }


def descendants(parents, root):
    out = {root}
    changed = True
    while changed:
        changed = False
        for v in range(1, len(parents)):
            if parents[v] in out and v not in out:
                out.add(v)
                changed = True
    return out


def check_subtree_flip():
    rng = np.random.RandomState(3901)
    n = 6
    parents = (None, 0, 0, 1, 1, 3)
    raw = rng.choice((-1, 1), size=(n, n))
    a = np.triu(raw, 1)
    a = a + a.T
    x = [rng.choice((-1, 1), size=n).astype(np.int64) for _ in parents]
    selectors = [None] + [tuple(sorted(rng.choice(n, 4, replace=False))) for _ in parents[1:]]
    sectors = [None] + [int(rng.choice((-1, 1))) for _ in parents[1:]]
    before = edge_sets(x, parents)

    # Every edge is checked independently.  Flipping its disagreement
    # coordinates throughout the descendant subtree deletes exactly that
    # edge and leaves all other edge disagreement sets unchanged.
    for edge, u in before.items():
        _, child = edge
        shore = descendants(parents, child)
        xt = [w.copy() for w in x]
        for v in shore:
            if u:
                xt[v][list(u)] *= -1
        after = edge_sets(xt, parents)
        assert after[edge] == frozenset()
        for other in before:
            if other != edge:
                assert after[other] == before[other]
        assert sum(map(len, after.values())) == sum(map(len, before.values())) - len(u)

        # Check the exact signed-deficit variation at every descendant
        # terminal: delta(y^F)-delta(y)=4 C_y(F).
        for v in shore:
            if v == 0:
                continue
            s = selectors[v]
            loc = {j: k for k, j in enumerate(s)}
            f = tuple(j for j in u if j in loc)
            b = a[np.ix_(s, s)]
            y = x[v][list(s)]
            yt = y.copy()
            for j in f:
                yt[loc[j]] *= -1
            eps = sectors[v]
            q = qnorm(b)
            delta = q - eps * int(y @ b @ y)
            deltat = q - eps * int(yt @ b @ yt)
            fpos = {loc[j] for j in f}
            c = eps * sum(
                int(b[i, j] * y[i] * y[j])
                for i in fpos
                for j in range(len(s))
                if j not in fpos
            )
            assert deltat - delta == 4 * c

    # The edge-coordinate disagreement sets are an exact sparse cut code.
    xr = [x[0].copy()]
    for v in range(1, len(parents)):
        w = xr[parents[v]].copy()
        u = before[(parents[v], v)]
        if u:
            w[list(u)] *= -1
        xr.append(w)
    assert all(np.array_equal(u, v) for u, v in zip(x, xr))
    return sum(map(len, before.values()))


def check_complete_signing_wall():
    # Symmetric Hadamard R, H=[[R,-R],[-R,R]], and two blocks.  H 1=0,
    # ||K|| ||H|| / b = 1/2 <= 1-2/b, so the proof in the memo applies.
    r = np.array(
        [
            [1, 1, 1, 1],
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, -1, -1, 1],
        ],
        dtype=np.int64,
    )
    h = np.block([[r, -r], [-r, r]])
    b = len(h)
    internal = np.ones((b, b), dtype=np.int64) - np.eye(b, dtype=np.int64)
    a = np.block([[internal, h], [h, internal]])
    n = len(a)
    e0 = n * (b - 1)
    assert np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0)
    assert np.all(np.abs(a + np.eye(n, dtype=np.int64)) == 1)
    assert np.array_equal(h @ np.ones(b, dtype=np.int64), np.zeros(b, dtype=np.int64))
    assert abs(np.linalg.norm(h, 2) - 4.0) < 1e-10

    energies = []
    for tail in itertools.product((-1, 1), repeat=n - 1):
        z = np.array((1,) + tail, dtype=np.int64)
        energies.append(int(z @ a @ z))
    assert max(map(abs, energies)) == e0

    x = np.ones(n, dtype=np.int64)
    y = x.copy()
    y[:b] *= -1
    assert int(x @ a @ x) == int(y @ a @ y) == e0
    assert np.array_equal(a @ x, (b - 1) * x)
    assert np.array_equal(a @ y, (b - 1) * y)
    shore = sum(int(a[i, j] * x[i] * x[j]) for i in range(b) for j in range(b, n))
    assert shore == 0
    assert min(np.count_nonzero(x != y), np.count_nonzero(x != -y)) == n // 2
    return n, e0, int((a @ x) @ (a @ x))


if __name__ == "__main__":
    cut_size = check_subtree_flip()
    n, q, row = check_complete_signing_wall()
    print("PASS subtree flip and exact cut code; sample cut size =", cut_size)
    print("PASS complete-signing wall; n, Q(A), R2(ground) =", n, q, row)
