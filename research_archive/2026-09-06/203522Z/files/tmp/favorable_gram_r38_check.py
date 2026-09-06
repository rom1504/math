#!/usr/bin/env python3
"""Exact finite checks for the Wave 38 completion/Gram bridge.

Only integer arithmetic is used.  For each exact child ground we orient the
whole signing so that the child energy is positive; this does not change any
Gram norm or the completion width.
"""

from __future__ import annotations

import itertools
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def all_spins(n: int):
    for bits in itertools.product((-1, 1), repeat=n):
        yield np.asarray(bits, dtype=np.int64)


def q_value(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def child_ground_records(a: np.ndarray, m: int):
    """Enumerate completion and Gram data for all projective child grounds."""
    n = len(a)
    q = q_value(a)
    records = []
    for s_tuple in itertools.combinations(range(n), m):
        s = np.asarray(s_tuple, dtype=int)
        t = np.asarray([i for i in range(n) if i not in s_tuple], dtype=int)
        child = a[np.ix_(s, s)]
        q_child = q_value(child)
        for y in projective_spins(m):
            raw_e = int(y @ child @ y)
            if abs(raw_e) != q_child:
                continue
            sigma = 1 if raw_e >= 0 else -1
            oriented = sigma * a
            e = abs(raw_e)
            h = oriented[np.ix_(t, s)] @ y
            ext = int(h @ h)
            internal_vector = oriented[np.ix_(s, s)] @ y
            internal = int(internal_vector @ internal_vector)

            vals = []
            for x in all_spins(len(t)):
                vals.append(
                    int(x @ oriented[np.ix_(t, t)] @ x + 2 * y @ oriented[np.ix_(s, t)] @ x)
                )
            zp = max(vals)
            zm = -min(vals)
            assert (zp + zm) % 2 == 0 and (zp - zm) % 2 == 0
            w = (zp + zm) // 2
            center = (zp - zm) // 2
            b = q - w

            # Exact Fourier/Parseval identity and interval bounds.
            assert sum(vals) == 0
            assert sum(v * v for v in vals) == len(vals) * (
                2 * len(t) * (len(t) - 1) + 4 * ext
            )
            assert 2 * len(t) * (len(t) - 1) + 4 * ext <= w * w - center * center
            assert w * w - center * center <= q * q - e * e
            assert abs(center + e) <= b

            records.append(
                {
                    "S": s_tuple,
                    "y": tuple(int(v) for v in y),
                    "q": q,
                    "q_child": q_child,
                    "e": e,
                    "w": w,
                    "a": center,
                    "b": b,
                    "ext": ext,
                    "internal": internal,
                }
            )
    return records


def audit_child_grounds(a: np.ndarray, m: int, name: str) -> None:
    records = child_ground_records(a, m)
    groups = defaultdict(lambda: {"external": set(), "internal": set(), "count": 0})
    for r in records:
        key = (r["w"], r["a"], r["b"])
        groups[key]["external"].add(r["ext"])
        groups[key]["internal"].add(r["internal"])
        groups[key]["count"] += 1

    print(name, "exact-child-ground records", len(records))
    for key in sorted(groups):
        g = groups[key]
        print(
            "  (w,a,b)=", key,
            "external=", sorted(g["external"]),
            "internal=", sorted(g["internal"]),
            "count=", g["count"],
        )

    if name == "A9,m=6":
        g = groups[(4, -2, 20)]
        assert g["external"] == {0}
        assert g["internal"] == {46, 86}
        witness = next(
            r for r in records
            if (r["w"], r["a"], r["b"], r["ext"], r["internal"])
            == (4, -2, 20, 0, 86)
        )
        print("  zero-external/high-internal witness", witness)


def audit_complement_certificates(a: np.ndarray, m: int, name: str) -> None:
    """Check the small-deficit certificate-to-external-Gram implication."""
    n = len(a)
    q = q_value(a)
    selectors = list(itertools.combinations(range(n), m))
    states = []
    for x in projective_spins(n):
        raw = int(x @ a @ x)
        for sigma in (-1, 1):
            energy = sigma * raw
            states.append((x, sigma, energy, q - energy))

    count = 0
    min_slack = None
    for s_tuple in selectors:
        s = np.asarray(s_tuple, dtype=int)
        t = np.asarray([i for i in range(n) if i not in s_tuple], dtype=int)
        child = a[np.ix_(s, s)]
        q_child = q_value(child)
        for x, sigma, energy, deficit in states:
            y = x[s]
            e = sigma * int(y @ child @ y)
            if 2 * e - energy < q:
                continue
            count += 1
            ext_vector = a[np.ix_(t, s)] @ y
            ext = int(ext_vector @ ext_vector)
            assert q_child - e <= deficit // 2 if deficit % 2 == 0 else q_child - e <= deficit / 2
            # 4 E_ext plus the outside quadratic Fourier floor is bounded by
            # q*Delta-Delta^2/4 (and hence by q*Delta).
            lhs = 4 * ext + 2 * len(t) * (len(t) - 1)
            rhs = q * deficit
            assert lhs <= rhs
            assert 4 * lhs <= 4 * q * deficit - deficit * deficit
            slack = rhs - lhs
            min_slack = slack if min_slack is None else min(min_slack, slack)

    print(
        name,
        "complement certificates", count,
        "minimum exact q*Delta-[4Eext+2t(t-1)] slack", min_slack,
    )


def main() -> None:
    audit_child_grounds(A6, 3, "A6,m=3")
    audit_child_grounds(A8, 5, "A8,m=5")
    audit_child_grounds(A9, 6, "A9,m=6")
    audit_complement_certificates(A6, 5, "A6,m=5")
    audit_complement_certificates(A8, 6, "A8,m=6")
    audit_complement_certificates(A9, 7, "A9,m=7")
    print("PASS favorable_gram_r38_check")


if __name__ == "__main__":
    main()
