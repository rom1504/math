#!/usr/bin/env python3
"""Exact finite checks for drafts/mesoscopic_pair_query_visibility.md.

The proof is asymptotic, but its selector and completion identities are
finite.  This script checks those identities by exhaustive Boolean
optimization at orders small enough that the complete 2n-vertex parent can
also be enumerated.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from pathlib import Path

import numpy as np


def spins(n: int):
    return itertools.product((-1, 1), repeat=n)


def qform(a: list[list[int]], x: tuple[int, ...]) -> int:
    n = len(x)
    return sum(a[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))


def cross(b: list[list[int]], x: tuple[int, ...], y: tuple[int, ...]) -> int:
    return sum(b[i][j] * x[i] * y[j] for i in range(len(x)) for j in range(len(y)))


def random_hollow(n: int, rng: random.Random) -> list[list[int]]:
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = rng.choice((-1, 1))
    return a


def positive_ground(a: list[list[int]]):
    vals = [(qform(a, x), x) for x in spins(len(a))]
    value, x = max(vals, key=lambda item: abs(item[0]))
    sigma = 1 if value >= 0 else -1
    return abs(value), sigma, x


def bipartite_ground(c: list[list[int]]):
    best = None
    for x in spins(len(c)):
        fields = [sum(c[i][j] * x[i] for i in range(len(c))) for j in range(len(c[0]))]
        value = sum(abs(v) for v in fields)
        if best is None or value > best[0]:
            eta = tuple(1 if v >= 0 else -1 for v in fields)
            best = (value, x, eta, tuple(abs(v) for v in fields))
    assert best is not None
    return best


def run_case(n: int, seed: int):
    rng = random.Random(seed)
    a = random_hollow(n, rng)
    qa, sigma, u = positive_ground(a)
    k = min(2, n // 2)
    sset = set(range(k))
    v = tuple(-u[i] if i in sset else u[i] for i in range(n))
    outside = [i for i in range(n) if i not in sset]

    d = random_hollow(n, rng)
    qd, dsigma, y0 = positive_ground(d)
    if dsigma < 0:
        d = [[-entry for entry in row] for row in d]
    assert qform(d, y0) == qd

    raw_c = [[rng.choice((-1, 1)) for _ in range(n)] for _ in outside]
    lval, x0, eta, cabs = bipartite_ground(raw_c)
    row_gauge = [x0[ii] * u[i] for ii, i in enumerate(outside)]
    col_gauge = [eta[j] * y0[j] for j in range(n)]
    c = [
        [row_gauge[ii] * raw_c[ii][j] * col_gauge[j] for j in range(n)]
        for ii in range(len(outside))
    ]
    fields_u = [sum(c[ii][j] * u[i] for ii, i in enumerate(outside)) for j in range(n)]
    assert tuple(fields_u) == tuple(cabs[j] * y0[j] for j in range(n))

    b = [[0] * n for _ in range(n)]
    for i in range(n):
        if i in sset:
            b[i] = [u[i] * y0[j] for j in range(n)]
        else:
            ii = outside.index(i)
            b[i] = c[ii][:]

    fb = {}
    for x in spins(n):
        fb[x] = sum(abs(sum(b[i][j] * x[i] for i in range(n))) for j in range(n))
    fmax = max(fb.values())
    assert fb[u] == fmax == k * n + lval
    expected_v = sum(abs(cj - k) for cj in cabs)
    assert fb[v] == expected_v
    delta = fb[u] - fb[v]
    lower = 2.0 * k * lval * lval / (k * lval + sum(cj * cj for cj in cabs))
    assert delta + 1e-12 >= lower

    parent_cap = -1
    rooted_u = -1
    rooted_v = -1
    target = sigma * (qa + fb[u] + qd)
    for x in spins(n):
        rooted = -1
        for y in spins(n):
            value = qform(a, x) + sigma * cross(b, x, y) + sigma * qform(d, y)
            avalue = abs(value)
            parent_cap = max(parent_cap, avalue)
            rooted = max(rooted, avalue)
        if x == u:
            rooted_u = rooted
        if x == v:
            rooted_v = rooted
    assert abs(target) == qa + fb[u] + qd
    assert parent_cap == abs(target)
    assert rooted_u == parent_cap
    assert rooted_u - rooted_v >= delta

    op = float(np.linalg.svd(np.array(raw_c, dtype=float), compute_uv=False)[0])
    return {
        "n": n,
        "seed": seed,
        "k": k,
        "Q_A": qa,
        "Q_D": qd,
        "C_bipartite_cap": lval,
        "C_operator_norm": op,
        "B_cap": fmax,
        "pair_gap": delta,
        "moment_lower_bound": lower,
        "parent_cap": parent_cap,
        "rooted_rival_cap": rooted_v,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = [run_case(n, 1000 * n + seed) for n in range(4, 8) for seed in range(3)]
    out = {
        "description": "Exact checks of MQ.8 and MQ.14--MQ.25",
        "records": records,
        "all_checks_passed": True,
    }
    payload = json.dumps(out, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
