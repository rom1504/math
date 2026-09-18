#!/usr/bin/env python3
"""Exact checks for the Wave 58 relative-fibre transfer and finite wall.

Integer energy, row, cut/uncut, and selector identities are checked exactly.
The bare threshold for m=6 contains sqrt(3), so its exhaustive column census
uses floating evaluation and asserts a positive distance from every boundary.
"""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A8, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def pk(m: int, n: int, k: int) -> Fraction:
    if m < k:
        return Fraction(0)
    ans = Fraction(1)
    for j in range(k):
        ans *= Fraction(m - j, n - j)
    return ans


def pair_data(a: np.ndarray, tau: int, x: np.ndarray,
              upsilon: int, z: np.ndarray):
    """Relative cut/uncut data, with each undirected edge counted once."""
    n = len(a)
    kappa = tau * upsilon
    y = x * z
    s = tau * a * np.outer(x, x)
    mask = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if kappa * y[i] * y[j] == -1:
                mask[i, j] = mask[j, i] = 1
    g = s * mask
    b = g @ np.ones(n, dtype=np.int64)
    w = int(np.triu(g, 1).sum())
    r = s @ np.ones(n, dtype=np.int64)
    sw = upsilon * a * np.outer(z, z)
    rw = sw @ np.ones(n, dtype=np.int64)
    return kappa, y, s, mask, g, b, w, r, rw


def global_pair_identity_check() -> None:
    """Check the orientation-unified version of the old shore-row identity."""
    a = A8
    q = qnorm(a)
    xs = list(projective_spins(len(a)))
    # A deterministic broad sample includes both relative sectors.
    indices = [(i, (37 * i + 11) % len(xs)) for i in range(len(xs))]
    for i, j in indices:
        x, z = xs[i], xs[j]
        for tau in (-1, 1):
            for upsilon in (-1, 1):
                kappa, y, s, mask, g, b, w, r, rw = pair_data(
                    a, tau, x, upsilon, z
                )
                # The mask is a cut for kappa=+1 and an uncut for kappa=-1.
                expected = np.zeros_like(mask)
                for u in range(len(a)):
                    for v in range(u + 1, len(a)):
                        cut = y[u] != y[v]
                        hit = cut if kappa == 1 else not cut
                        expected[u, v] = expected[v, u] = int(hit)
                assert np.array_equal(mask, expected)
                assert np.array_equal(rw, r - 2 * b)
                ed = tau * int(x @ a @ x)
                ew = upsilon * int(z @ a @ z)
                dd, dw = q - ed, q - ew
                assert dw - dd == 4 * w
                rd, roww = int(r @ r), int(rw @ rw)
                assert roww - rd == 4 * int(b @ b - r @ b)
                # G=(S-kappa YSY)/2, so masking cannot increase op norm.
                assert np.linalg.norm(g, 2) <= np.linalg.norm(a, 2) + 1e-10


def fibre_and_variance_check() -> None:
    """Check exact bare-fibre transfer and its fixed-slice variance formula."""
    a = A8
    n = len(a)
    xs = list(projective_spins(n))
    pairs = [
        (-1, xs[12], -1, xs[41]),
        (-1, xs[12], 1, xs[15]),
        (1, xs[3], -1, xs[79]),
    ]
    for m in (2, 3, 6):
        p2, p3, p4 = (pk(m, n, k) for k in (2, 3, 4))
        alpha = p2 - 2 * p3 + p4
        beta = p3 - p4
        gamma = p4 - p2 * p2
        assert gamma <= 0
        selectors = list(itertools.combinations(range(n), m))
        for tau, x, upsilon, z in pairs:
            _, _, s, mask, _, b, w, _, _ = pair_data(
                a, tau, x, upsilon, z
            )
            values = []
            for selector in selectors:
                ss = list(selector)
                ws = int(np.triu((s * mask)[np.ix_(ss, ss)], 1).sum())
                values.append(ws)
                # Q(A[S]) and B cancel from the difference, leaving this.
                local_d = tau * int(x[ss] @ a[np.ix_(ss, ss)] @ x[ss])
                local_w = upsilon * int(z[ss] @ a[np.ix_(ss, ss)] @ z[ss])
                delta_difference = local_d - local_w
                assert delta_difference == 4 * ws
                assert Fraction(delta_difference, 1) - 4 * p2 * w \
                    == 4 * (Fraction(ws, 1) - p2 * w)
            mean = sum(Fraction(v, 1) for v in values) / len(values)
            assert mean == p2 * w
            variance = sum((Fraction(v, 1) - mean) ** 2 for v in values) \
                / len(values)
            ecount = int(np.triu(mask, 1).sum())
            formula = alpha * ecount + beta * int(b @ b) + gamma * w * w
            assert variance == formula
            assert variance <= alpha * ecount + beta * int(b @ b)


def bare_columns(a: np.ndarray, m: int, t: float = 0.0):
    n = len(a)
    q = qnorm(a)
    selectors = list(itertools.combinations(range(n), m))
    qs = [qnorm(a[np.ix_(s, s)]) for s in selectors]
    p2 = m * (m - 1) / (n * (n - 1))
    allowance = ((m / n) ** 1.5 - p2) * q + t
    records = []
    margin = math.inf
    for x in projective_spins(n):
        raw = int(x @ a @ x)
        row = int((a @ x) @ (a @ x))
        for tau in (-1, 1):
            deficit = q - tau * raw
            good = []
            for selector, qs_value in zip(selectors, qs):
                ss = list(selector)
                local = tau * int(x[ss] @ a[np.ix_(ss, ss)] @ x[ss])
                residual = qs_value - local - p2 * deficit - allowance
                margin = min(margin, abs(residual))
                good.append(residual <= 0)
            count = sum(good)
            if count:
                records.append({
                    "x": x.copy(), "tau": tau, "deficit": deficit,
                    "row": row, "count": count, "good": np.asarray(good),
                    "h": -math.log(count / len(selectors)),
                })
    assert margin > 1e-6
    return records, selectors, margin


def find_record(records, tau, x):
    return next(r for r in records if r["tau"] == tau
                and np.array_equal(r["x"], x))


def finite_descent_wall_check() -> None:
    """A8 bare-column scalar optimum has migrated responses but no descent."""
    a = A8
    assert qnorm(a) == 20  # the stored exact order-eight minimizer
    records, selectors, margin = bare_columns(a, 6)
    lam = 1.0 / 1000.0
    objective = np.asarray([r["h"] + lam * r["row"] for r in records])
    best = float(objective.min())
    optimal = [r for r, f in zip(records, objective) if abs(f - best) < 1e-12]
    assert len(optimal) == 8
    assert {(r["count"], r["row"], r["deficit"]) for r in optimal} \
        == {(18, 72, 4)}
    second_gap = min(f - best for f in objective if f > best + 1e-12)
    assert abs(second_gap - math.log(9 / 8)) < 1e-12

    x = np.asarray([1, -1, -1, -1, 1, 1, -1, -1], dtype=np.int64)
    d = find_record(records, -1, x)
    fields = -x * (a @ x)
    assert fields.tolist() == [3, 5, 1, -1, -1, 3, 5, 1]
    assert (d["count"], d["row"], d["deficit"]) == (18, 72, 4)
    assert int(np.maximum(fields, 0).sum()) == 18
    # L=2 obeys the finite high-row premise R>=2(n-1)L.  Both edges are
    # positive at d and form the field-positive block used below.
    edge_block = [(0, 1), (0, 2)]
    assert d["row"] >= 2 * (len(a) - 1) * len(edge_block)
    assert all(-a[i, j] * x[i] * x[j] == 1 for i, j in edge_block)

    # An exact parent ground reverses one half of E and lowers row, but loses
    # seven of the eighteen favorable selectors, so the scalar cost rises.
    z = np.asarray([1, -1, 1, -1, 1, -1, -1, 1], dtype=np.int64)
    omega = find_record(records, -1, z)
    reversed_edges = sum(
        -a[i, j] * z[i] * z[j] == -1 for i, j in edge_block
    )
    assert reversed_edges == 1 and reversed_edges > len(edge_block) / 4
    assert (omega["count"], omega["row"], omega["deficit"]) == (11, 64, 0)
    assert int((d["good"] & omega["good"]).sum()) == 8
    objective_change = omega["h"] + lam * omega["row"] \
        - d["h"] - lam * d["row"]
    assert abs(objective_change - (math.log(18 / 11) - 8 / 1000)) < 1e-12
    assert objective_change > 0

    # A second equally near response reverses all E, yet does not lower row
    # and also has larger scalar cost.
    z2 = np.asarray([1, -1, -1, -1, 1, 1, 1, 1], dtype=np.int64)
    omega2 = find_record(records, 1, z2)
    reversed2 = sum(
        a[i, j] * z2[i] * z2[j] == -1 for i, j in edge_block
    )
    assert reversed2 == 2
    assert (omega2["count"], omega2["row"], omega2["deficit"]) == (16, 72, 4)
    assert omega2["h"] + lam * omega2["row"] > best

    print({
        "A8_m6_boundary_margin": margin,
        "scalar_optimum_type": (18, 72, 4),
        "second_objective_gap": second_gap,
        "ground_response": (11, 64, 0),
        "ground_response_objective_change": objective_change,
    })


def exponent_check() -> None:
    """Audit the powers in the project-scale relative-transfer estimate."""
    # T^2/R_* = n^(3-2c-(9/4-c)) = n^(3/4-c) = H.
    for c in (Fraction(1, 20), Fraction(1, 8), Fraction(1, 5)):
        t_exp = Fraction(3, 2) - c
        target_row_exp = Fraction(9, 4) - c
        h_exp = Fraction(3, 4) - c
        assert 2 * t_exp - target_row_exp == h_exp
        assert t_exp - Fraction(3, 4) == h_exp
        # With only the generic response row n^(5/2), the Gaussian exponent
        # drops strictly below H throughout 0<c<1/4.
        assert 2 * t_exp - Fraction(5, 2) == Fraction(1, 2) - 2 * c
        assert Fraction(1, 2) - 2 * c < h_exp


def main() -> None:
    global_pair_identity_check()
    fibre_and_variance_check()
    finite_descent_wall_check()
    exponent_check()
    print("PASS relative cut/uncut transfer, variance, and finite bare-row wall")


if __name__ == "__main__":
    main()
