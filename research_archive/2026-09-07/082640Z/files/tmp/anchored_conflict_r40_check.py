#!/usr/bin/env python3
"""Exact finite audit for Wave 40 anchored-conflict four-corner identities."""

from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A9


A5 = np.asarray([
    [0, -1, -1, -1, -1],
    [-1, 0, -1, -1, 1],
    [-1, -1, 0, 1, -1],
    [-1, -1, 1, 0, 1],
    [-1, 1, -1, 1, 0],
], dtype=np.int64)


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def q_value(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def complement(a: np.ndarray, s0: tuple[int, ...]) -> np.ndarray:
    b = -a.copy()
    s = np.asarray(s0, dtype=int)
    b[np.ix_(s, s)] = a[np.ix_(s, s)]
    return b


def cut(a: np.ndarray, x: np.ndarray, sigma: int, u: set[int]) -> int:
    return sum(
        sigma * int(a[i, j] * x[i] * x[j])
        for i in u for j in range(len(a)) if j not in u
    )


def cut_matrix(x: np.ndarray, sigma: int) -> np.ndarray:
    d = sigma * np.outer(x, x)
    np.fill_diagonal(d, 0)
    return d


def audit_pair(
    a: np.ndarray,
    s0: tuple[int, ...],
    t0: tuple[int, ...],
    x: np.ndarray,
    y: np.ndarray,
    sigma: int,
    anchor: int,
) -> dict[str, object]:
    n = len(a)
    q = q_value(a)
    assert x[anchor] == y[anchor] == 1
    bs, bt = complement(a, s0), complement(a, t0)
    ms, mt = q_value(bs), q_value(bt)
    ex = sigma * int(x @ a @ x)
    ey = sigma * int(y @ a @ y)
    lx = sigma * int(x @ bs @ x)
    ly = sigma * int(y @ bt @ y)
    dx, dy = q - ex, q - ey
    gx, gy = ms - lx, mt - ly
    rx = int((a @ x) @ (a @ x))
    ry = int((a @ y) @ (a @ y))

    global_disagreement = {i for i in range(n) if x[i] != y[i]}
    overlap = set(s0).intersection(t0)
    hset = (global_disagreement & overlap) - {anchor}
    rest = global_disagreement - hset

    cx = cut(a, x, sigma, hset)
    cy = cut(a, y, sigma, hset)
    j = sum(
        sigma * int(a[i, k] * x[i] * x[k])
        for i in hset for k in rest
    )
    assert cx + cy == 2 * j

    xh, yh = x.copy(), y.copy()
    for i in hset:
        xh[i] *= -1
        yh[i] *= -1
    assert sigma * int(xh @ a @ xh) + sigma * int(yh @ a @ yh) == ex + ey - 8 * j
    assert -(dx + dy) <= 8 * j <= 4 * q - dx - dy

    cbx = cut(bs, x, sigma, hset)
    cby = cut(bt, y, sigma, hset)
    xi = (cbx + cby) // 2
    assert cbx + cby == 2 * xi
    assert 8 * xi >= -(gx + gy)

    # Exact row crossover identity.
    uh = np.zeros(n, dtype=np.int64)
    ur = np.zeros(n, dtype=np.int64)
    uh[list(hset)] = x[list(hset)]
    ur[list(rest)] = x[list(rest)]
    row_cross = int((a @ uh) @ (a @ ur))
    rxh = int((a @ xh) @ (a @ xh))
    ryh = int((a @ yh) @ (a @ yh))
    assert rxh + ryh == rx + ry - 8 * row_cross

    d, e = cut_matrix(x, sigma), cut_matrix(y, sigma)
    cyclic = int(np.sum((bs - bt) * (d - e)))
    support_formula = 0
    for i in range(n):
        for k in range(i + 1, n):
            if (i in global_disagreement) == (k in global_disagreement):
                continue
            coeff = int(i in s0 and k in s0) - int(i in t0 and k in t0)
            support_formula += 8 * sigma * int(a[i, k] * x[i] * x[k]) * coeff
    assert cyclic == support_formula

    return {
        "q": q,
        "M": (ms, mt),
        "parent_energy": (ex, ey),
        "parent_deficit": (dx, dy),
        "complement_energy": (lx, ly),
        "complement_cap_slack": (gx, gy),
        "row": (rx, ry),
        "H": tuple(sorted(hset)),
        "global_disagreement": tuple(sorted(global_disagreement)),
        "parent_J": j,
        "complement_Xi": xi,
        "row_cross": row_cross,
        "cyclic": cyclic,
    }


def enumerate_q_min(n: int) -> int:
    edges = list(itertools.combinations(range(n), 2))
    answer = 10**9
    for signs in itertools.product((-1, 1), repeat=len(edges)):
        a = np.zeros((n, n), dtype=np.int64)
        for (i, j), value in zip(edges, signs):
            a[i, j] = a[j, i] = value
        answer = min(answer, q_value(a))
    return answer


def no_smaller_exact_pair(n: int) -> bool:
    """Exhaust proper high-ratio m=n-1 minimizers for the zero-gap pattern."""
    m = n - 1
    edges = list(itertools.combinations(range(n), 2))
    qn = enumerate_q_min(n)
    for signs in itertools.product((-1, 1), repeat=len(edges)):
        a = np.zeros((n, n), dtype=np.int64)
        for (i, j), value in zip(edges, signs):
            a[i, j] = a[j, i] = value
        if q_value(a) != qn:
            continue
        selectors = list(itertools.combinations(range(n), m))
        records: dict[tuple[int, ...], list[tuple[np.ndarray, int]]] = {}
        for s0 in selectors:
            b = complement(a, s0)
            mb = q_value(b)
            rows = []
            for x in projective_spins(n):
                for sigma in (-1, 1):
                    if sigma * int(x @ a @ x) != qn:
                        continue
                    if sigma * int(x @ b @ x) == mb and mb >= qn:
                        rows.append((x, sigma))
            records[s0] = rows
        for s0, t0 in itertools.combinations(selectors, 2):
            for anchor in set(s0).intersection(t0):
                for (x0, sigma), (y0, tau) in itertools.product(records[s0], records[t0]):
                    if sigma != tau:
                        continue
                    x, y = x0 * x0[anchor], y0 * y0[anchor]
                    h = (set(s0).intersection(t0) - {anchor})
                    if any(x[i] != y[i] for i in h):
                        return False
    return True


def main() -> None:
    assert enumerate_q_min(3) == 6
    assert enumerate_q_min(4) == 8
    assert enumerate_q_min(5) == 8
    assert q_value(A5) == 8
    assert no_smaller_exact_pair(3)
    assert no_smaller_exact_pair(4)

    s5 = (0, 1, 2, 3)
    t5 = (0, 1, 2, 4)
    # Parent cone is saturated (J=0), while overlap disagreement is maximal.
    a5_parent_flat = audit_pair(
        A5, s5, t5,
        np.asarray([1, 1, -1, -1, -1]),
        np.asarray([1, -1, 1, -1, -1]),
        1, 0,
    )
    # Complement cone is saturated (Xi=0), again with maximal disagreement.
    a5_complement_flat = audit_pair(
        A5, s5, t5,
        np.asarray([1, 1, -1, -1, 1]),
        np.asarray([1, -1, 1, 1, -1]),
        1, 0,
    )
    for result in (a5_parent_flat, a5_complement_flat):
        assert result["parent_deficit"] == (0, 0)
        assert result["complement_cap_slack"] == (0, 0)
        assert result["M"] == result["complement_energy"] == (8, 8)
        assert result["row"] == (24, 24)
        assert result["H"] == (1, 2)
    assert a5_parent_flat["global_disagreement"] == (1, 2)
    assert a5_parent_flat["parent_J"] == 0
    assert a5_complement_flat["global_disagreement"] == (1, 2, 3, 4)
    assert a5_complement_flat["complement_Xi"] == 0

    # Exact complement maximizers with zero cyclic-monotonicity gap but
    # nonzero overlap conflict on the known exact order-nine minimizer.
    a9_flat_normal = audit_pair(
        A9,
        (0, 1, 2, 3, 4, 5, 6),
        (0, 1, 2, 3, 4, 5, 7),
        np.asarray([1, -1, 1, 1, 1, 1, -1, 1, 1]),
        np.asarray([1, 1, 1, 1, 1, 1, -1, 1, 1]),
        1, 0,
    )
    assert a9_flat_normal["q"] == 24
    assert a9_flat_normal["M"] == a9_flat_normal["complement_energy"] == (28, 28)
    assert a9_flat_normal["complement_cap_slack"] == (0, 0)
    assert a9_flat_normal["H"] == (1,)
    assert a9_flat_normal["global_disagreement"] == (1,)
    assert a9_flat_normal["cyclic"] == 0

    print("A5 parent-flat", a5_parent_flat)
    print("A5 complement-flat", a5_complement_flat)
    print("A9 normal-fan-flat", a9_flat_normal)
    print("PASS anchored_conflict_r40_check")


if __name__ == "__main__":
    main()
