"""MILP search for the smallest parent embedding of the (10.409) gap.

The six-vertex core is fixed.  Its negative endpoint n (energy -10) is
required to be the restriction of a positive absolute ground of a larger
signed complete graph.  All edges touching k new vertices are free signs.
After switching by the target ground, the target is all-one and the core is
gauged by n.  Ground inequalities are linear in the free edge signs.
"""

from itertools import product

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


A = np.array(
    [
        [0, -1, 1, -1, -1, -1],
        [-1, 0, 1, -1, 1, 1],
        [1, 1, 0, -1, 1, -1],
        [-1, -1, -1, 0, 1, -1],
        [-1, 1, 1, 1, 0, -1],
        [-1, 1, -1, -1, -1, 0],
    ],
    dtype=int,
)
nspin = np.array([-1, -1, 1, -1, 1, 1], dtype=int)


def energy(a, x):
    return int(x @ a @ x)


def solve(k, C, label, verbose=True):
    total = k + 6
    # New vertices are 0..k-1, core vertices k..k+5.
    edges = [(i, j) for i in range(total) for j in range(i + 1, total) if i < k]
    m = len(edges)
    rows = []
    ub = []
    for prefix in product((-1, 1), repeat=total - 1):
        x = np.array(prefix + (1,), dtype=int)
        xc = x[k:]
        hc = energy(C, xc)
        q = np.array([x[i] * x[j] for i, j in edges], dtype=int)
        # s=2z-1.  H <= M, where M=-10+2 sum(s).
        # 2 sum s(q-1) <= -10-hc.
        coeff = 4 * (q - 1)
        rhs = -10 - hc + 2 * int(np.sum(q - 1))
        rows.append(coeff)
        ub.append(rhs)
        # -H <= M: -2 sum s(q+1) <= -10+hc.
        coeff = -4 * (q + 1)
        rhs = -10 + hc - 2 * int(np.sum(q + 1))
        rows.append(coeff)
        ub.append(rhs)
    # M >= 0: -10+2 sum(2z-1)>=0.
    rows.append(-4 * np.ones(m))
    ub.append(-10 - 2 * m)
    con = LinearConstraint(np.array(rows, dtype=float), -np.inf, np.array(ub, dtype=float))
    # Minimize target M, equivalently number of positive new edges.
    result = milp(
        c=np.ones(m),
        integrality=np.ones(m),
        bounds=Bounds(np.zeros(m), np.ones(m)),
        constraints=con,
        options={"time_limit": 120},
    )
    if verbose:
        print("k", k, label, "status", result.status, result.message)
    if result.x is None:
        return None
    z = np.rint(result.x).astype(int)
    s = 2 * z - 1
    parent = np.zeros((total, total), dtype=int)
    parent[k:, k:] = C
    for val, (i, j) in zip(s, edges):
        parent[i, j] = parent[j, i] = val
    vals = []
    grounds = []
    for prefix in product((-1, 1), repeat=total - 1):
        x = np.array(prefix + (1,), dtype=int)
        v = energy(parent, x)
        vals.append(v)
        if np.all(x == 1):
            target = v
    qabs = max(max(vals), -min(vals))
    if verbose:
        print("  target/Q/min/max", target, qabs, min(vals), max(vals), "new signs", s.tolist())
    assert target == qabs
    if verbose:
        print("  parent gauged")
        for row in parent:
            print(" ", row.tolist())
    return parent


def core_candidates():
    seen = set()
    # Positive parent orientation carrying a negative core endpoint.
    # Negative parent orientation carrying a positive endpoint is equivalent,
    # after negating/gauging the whole parent, to core -diag(p) A diag(p).
    rows = []
    for x in product((-1, 1), repeat=6):
        x = np.array(x, dtype=int)
        v = energy(A, x)
        if v == -10:
            rows.append((x[:, None] * A * x[None, :], ("negative", tuple(x))))
        if v == 10:
            rows.append((-x[:, None] * A * x[None, :], ("positive/opposite orientation", tuple(x))))
    for c, label in rows:
        assert energy(c, np.ones(6, dtype=int)) == -10
        key = tuple(c.flat)
        if key not in seen:
            seen.add(key)
            yield c, label


print("candidate gauged gap states", len(tuple(core_candidates())))
found = None
for k in range(1, 6):
    for idx, (core, label) in enumerate(core_candidates()):
        ans = solve(k, core, label, verbose=False)
        if ans is not None:
            print("FIRST FEASIBLE among all endpoint/orientation choices", "k", k, "candidate", idx, label)
            solve(k, core, label, verbose=True)
            found = (k, core, label)
            break
    if found:
        break
