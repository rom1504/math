#!/usr/bin/env python3
"""Exact finite checks for the Wave 53 low-cross/box reduction.

The matrices are stored exact minimizers (the order-ten matrix is one
deterministic sample, not a classification).  Everything below is integer
enumeration.  No optimization solver or randomness is used.
"""

from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A as A10
from envelope_block_cover_r27 import A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def child_records(a: np.ndarray, m: int):
    n = len(a)
    records = []
    for s0 in itertools.combinations(range(n), m):
        s = np.asarray(s0, dtype=int)
        t = np.asarray(tuple(i for i in range(n) if i not in s0), dtype=int)
        b = a[np.ix_(s, s)]
        vals = [(y, int(y @ b @ y)) for y in projective_spins(m)]
        qs = max(abs(v) for _, v in vals)
        for y, raw in vals:
            if abs(raw) != qs:
                continue
            sigma = 1 if raw > 0 else -1
            r = sigma * y * (b @ y)
            assert np.all(r >= 0)
            assert int(np.sum(r)) == qs
            assert int(np.sum(r * r)) == int((b @ y) @ (b @ y))
            internal = int(np.sum(r * r))
            cross_vec = a[np.ix_(t, s)] @ y
            cross = int(cross_vec @ cross_vec)
            base = internal + cross + (n - m) * (n - 1)

            rows = []
            centered = []
            linear = []
            full_b = a[:, s] @ y
            full_c = a[:, t]
            h = a @ a - (n - 1) * np.eye(n, dtype=np.int64)
            d = (a @ a)[np.ix_(t, s)] @ y
            for w0 in itertools.product((-1, 1), repeat=n - m):
                w = np.asarray(w0, dtype=np.int64)
                z = np.empty(n, dtype=np.int64)
                z[s] = y
                z[t] = w
                row = int((a @ z) @ (a @ z))
                rows.append(row)
                polynomial = int(2 * d @ w + w @ h[np.ix_(t, t)] @ w)
                centered.append(polynomial)
                linear.append(int(d @ w))
                assert row == int((full_b + full_c @ w) @ (full_b + full_c @ w))
                assert row == base + polynomial
            assert sum(centered) == 0
            assert sum(rows) == len(rows) * base
            assert len(rows) * (base - min(rows)) >= 2 * sum(map(abs, linear))
            assert min(rows) <= base - np.sqrt(2.0) * np.linalg.norm(d) + 1e-12
            records.append({
                "s": tuple(map(int, s)),
                "y": tuple(map(int, y)),
                "qs": qs,
                "r": tuple(map(int, r)),
                "rmax": int(np.max(r)),
                "internal": internal,
                "cross": cross,
                "base": base,
                "vmin": min(rows),
            })
    return records


def audit(name: str, a: np.ndarray, m: int):
    n = len(a)
    q = qnorm(a)
    rec = child_records(a, m)
    qstar = max(x["qs"] for x in rec)
    minimal_cross = min(x["cross"] for x in rec)
    min_cross_records = [x for x in rec if x["cross"] == minimal_cross]
    spikiest_min_cross = max(min_cross_records, key=lambda x: x["internal"])
    best_box = min(rec, key=lambda x: x["vmin"])
    best_mean = min(rec, key=lambda x: x["base"])
    maximal = [x for x in rec if x["qs"] == qstar]
    print(
        f"{name},m={m}: q={q}, qstar={qstar}, incidences={len(rec)}, "
        f"minV={best_box['vmin']}, minMean={best_mean['base']}, "
        f"minX={minimal_cross}"
    )
    print(
        "  spikiest at minX:",
        (spikiest_min_cross["s"], spikiest_min_cross["y"],
         spikiest_min_cross["r"], spikiest_min_cross["internal"],
         spikiest_min_cross["cross"], spikiest_min_cross["vmin"],
         spikiest_min_cross["base"]),
    )
    print(
        "  best actual box:",
        (best_box["s"], best_box["y"], best_box["qs"], best_box["r"],
         best_box["internal"], best_box["cross"], best_box["vmin"],
         best_box["base"]),
    )
    print(
        "  maximal-selector ranges (rmax,I,X,V,mean):",
        (min(x["rmax"] for x in maximal), max(x["rmax"] for x in maximal)),
        (min(x["internal"] for x in maximal), max(x["internal"] for x in maximal)),
        (min(x["cross"] for x in maximal), max(x["cross"] for x in maximal)),
        (min(x["vmin"] for x in maximal), max(x["vmin"] for x in maximal)),
        (min(x["base"] for x in maximal), max(x["base"] for x in maximal)),
    )

    # The local-field estimate used in the analytic reduction is exact.
    for x in rec:
        assert x["internal"] <= x["rmax"] * x["qs"]
        assert x["rmax"] <= min(m - 1, x["qs"] // 2)
        assert x["vmin"] <= x["base"]
    return rec


def main() -> None:
    audit("A8", A8, 5)
    audit("A8", A8, 6)
    a9m6 = audit("A9", A9, 6)
    audit("A9", A9, 7)
    a10m6 = audit("A10", A10, 6)

    # Exact low-cross spikiness at maximal selectors: X=0 does not make the
    # internal fields uniform or small, even on stored exact minimizers.
    assert any(
        x["qs"] == 22 and x["cross"] == 0 and x["internal"] == 86
        and x["rmax"] == 5
        for x in a9m6
    )
    assert any(
        x["qs"] == 22 and x["cross"] == 0 and x["internal"] == 86
        and x["rmax"] == 5
        for x in a10m6
    )

    # The mean certificate is sufficient, not necessary: the best A10 box
    # completion has row 10 although its uniform-completion mean is 74.
    best10 = min(a10m6, key=lambda x: x["vmin"])
    assert (
        best10["qs"], best10["internal"], best10["cross"],
        best10["vmin"], best10["base"]
    ) == (10, 30, 8, 10, 74)
    print("PASS exact completion identities and local-field inequalities")


if __name__ == "__main__":
    main()
