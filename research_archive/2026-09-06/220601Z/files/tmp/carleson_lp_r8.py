"""LP audit of the minimum normalized path-cover norm on endpoint trees.

Scratch only.  Energies use x^T A x.  A deterministic endpoint pair (or a
supplied tie selector later) generates the full binary tree.  The LP chooses
the conserved allocations and child obligations, rather than fixing them.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import random

import numpy as np
from scipy.optimize import linprog


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(len(a)))


def submatrix(a, ids):
    return tuple(tuple(a[i][j] for j in ids) for i in ids)


@lru_cache(None)
def extrema(a):
    if len(a) <= 1:
        return 0, 0, ((1,) if a else (),), ((1,) if a else (),)
    xs = tuple(product((-1, 1), repeat=len(a)))
    vals = tuple(energy(a, x) for x in xs)
    p = max(vals)
    mn = min(vals)
    return p, -mn, tuple(x for x, v in zip(xs, vals) if v == p), tuple(x for x, v in zip(xs, vals) if v == mn)


@dataclass
class Node:
    a: tuple
    p: int
    n: int
    cap: tuple[tuple[int, int], tuple[int, int]] | None = None
    kids: tuple["Node", "Node"] | None = None

    @property
    def r(self):
        return self.p + self.n

    @property
    def imb(self):
        return abs(self.p - self.n)


def build(a):
    pval, nval, ps, ns = extrema(a)
    node = Node(a, pval, nval)
    if len(a) <= 1:
        return node
    p, n = ps[0], ns[0]
    shores = (
        tuple(i for i in range(len(a)) if p[i] * n[i] == 1),
        tuple(i for i in range(len(a)) if p[i] * n[i] == -1),
    )
    assert all(shores), (a, p, n)
    caps = []
    kids = []
    for x in shores:
        y = tuple(i for i in range(len(a)) if i not in x)
        ax = submatrix(a, x)
        child = build(ax)
        px = tuple(p[i] for i in x)
        h = energy(ax, px)
        L = sum(abs(sum(a[j][i] * p[i] for i in x)) for j in y)
        q = max(child.p, child.n)
        caps.append((max(0, 2 * L - (q - h)), max(0, 2 * L - (q + h))))
        kids.append(child)
    node.cap = tuple(caps)
    node.kids = tuple(kids)
    return node


def solve(root, verbose=False):
    nodes = []
    def walk(v):
        idx = len(nodes)
        nodes.append(v)
        if v.kids:
            for ch in v.kids:
                walk(ch)
        return idx
    walk(root)
    ids = {id(v): i for i, v in enumerate(nodes)}

    names = []
    def var(name, ub=None):
        j = len(names)
        names.append(name)
        bounds.append((0, ub))
        return j
    bounds = []
    w = [var(f"w[{i}]", v.r if i == 0 else v.imb) for i, v in enumerate(nodes)]
    theta = [var(f"theta[{i}]") for i in range(len(nodes))]
    avar = {}
    zvar = {}
    for i, v in enumerate(nodes):
        if not v.kids:
            continue
        for side in range(2):
            zvar[i, side] = var(f"z[{i},{side}]")
            for sig in range(2):
                avar[i, side, sig] = var(f"a[{i},{side},{sig}]", v.cap[side][sig])

    eq_rows, eq_rhs = [], []
    ub_rows, ub_rhs = [], []
    def eq(co, rhs):
        row = np.zeros(len(names))
        for j, x in co.items(): row[j] = x
        eq_rows.append(row); eq_rhs.append(rhs)
    def le(co, rhs):
        row = np.zeros(len(names))
        for j, x in co.items(): row[j] = x
        ub_rows.append(row); ub_rhs.append(rhs)

    eq({w[0]: 1}, root.r)
    for i, v in enumerate(nodes):
        if not v.kids:
            eq({w[i]: 1}, 0)
            eq({theta[i]: 1}, 0)
            continue
        co = {w[i]: 1}
        for side, ch in enumerate(v.kids):
            ci = ids[id(ch)]
            co[w[ci]] = -1
            for sig in range(2): co[avar[i, side, sig]] = -1
        eq(co, 0)
        eq({theta[i]: 1, zvar[i, 0]: -1, zvar[i, 1]: -1}, 0)
        for side, ch in enumerate(v.kids):
            ci = ids[id(ch)]
            # z >= theta child
            le({theta[ci]: 1, zvar[i, side]: -1}, 0)
            # z >= sum a/c on this directed edge
            co = {zvar[i, side]: -1}
            for sig in range(2):
                c = v.cap[side][sig]
                if c:
                    co[avar[i, side, sig]] = 1 / c
            le(co, 0)

    obj = np.zeros(len(names)); obj[theta[0]] = 1
    res = linprog(obj, A_ub=np.array(ub_rows), b_ub=np.array(ub_rhs),
                  A_eq=np.array(eq_rows), b_eq=np.array(eq_rhs), bounds=bounds,
                  method="highs")
    assert res.success, res.message
    if verbose:
        print("K", res.fun, "R", root.r, "nodes", len(nodes))
        for name, x in zip(names, res.x):
            if x > 1e-7: print(name, x)
    return res.fun


def clique(n, sign=1):
    return tuple(tuple(0 if i == j else sign for j in range(n)) for i in range(n))


def random_signing(n, rng):
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = rng.choice((-1, 1))
    return tuple(tuple(r) for r in a)


if __name__ == "__main__":
    for n in (4, 8, 16):
        r = build(clique(n))
        print("clique", n, solve(r))
    rng = random.Random(8128)
    worst = (0, None)
    for n in range(3, 11):
        for trial in range(100):
            a = random_signing(n, rng)
            k = solve(build(a))
            if k > worst[0]:
                worst = (k, (n, trial, a))
                print("worst", worst[0], worst[1][:2])
    print("final", worst[0], worst[1][:2])
    solve(build(worst[1][2]), verbose=True)
