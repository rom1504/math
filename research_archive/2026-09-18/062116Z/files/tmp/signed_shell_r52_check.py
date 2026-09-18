#!/usr/bin/env python3
"""Exact Wave 52 checks for level-two row shells and selector exchanges."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A, aggregate_retention, build_ground_relation
from envelope_block_cover_r27 import A6, A8, A9


def falling(x: int, j: int) -> int:
    ans = 1
    for i in range(j):
        ans *= x - i
    return ans


def energy(a: np.ndarray, x: np.ndarray) -> int:
    return int(x @ a @ x)


def row_square(a: np.ndarray, x: np.ndarray) -> int:
    return int(np.sum((a @ x) ** 2))


def qnorm(a: np.ndarray) -> int:
    n = len(a)
    return max(
        abs(energy(a, np.asarray((1,) + tail, dtype=np.int64)))
        for tail in itertools.product((-1, 1), repeat=n - 1)
    )


def adjacent(s: tuple[int, ...], t: tuple[int, ...]) -> bool:
    return len(set(s) & set(t)) == len(s) - 1


def projective_spins(n: int) -> list[np.ndarray]:
    return [
        np.asarray((1,) + tail, dtype=np.int64)
        for tail in itertools.product((-1, 1), repeat=n - 1)
    ]


def build_relation(a: np.ndarray, m: int):
    n = len(a)
    centers = projective_spins(n)
    selectors = list(itertools.combinations(range(n), m))
    child_caps = []
    favorable = np.zeros((len(centers), len(selectors)), dtype=bool)
    for js, s in enumerate(selectors):
        child = a[np.ix_(s, s)]
        cap = qnorm(child)
        child_caps.append(cap)
        for iz, z in enumerate(centers):
            favorable[iz, js] = abs(energy(child, z[list(s)])) == cap
    rows = np.asarray([row_square(a, z) for z in centers], dtype=np.int64)
    return centers, selectors, child_caps, favorable, rows


def retention_generic(
    favorable: np.ndarray,
    rows: np.ndarray,
    selectors: list[tuple[int, ...]],
    cap: int,
    ell: int,
) -> Fraction | None:
    n = max(max(s) for s in selectors) + 1
    m = len(selectors[0])
    d = math.comb(m, ell) * math.comb(n - ell, m - ell)
    h2 = Fraction(1, math.comb(n - 2, m - 2))
    numerator = 0
    denominator = 0
    for iz in np.flatnonzero(rows <= cap):
        fam_ids = list(map(int, np.flatnonzero(favorable[iz])))
        r = len(fam_ids)
        if not r:
            continue
        numerator += r * sum(
            math.comb(len(set(selectors[i]) & set(selectors[j])), ell)
            for i in fam_ids for j in fam_ids
        )
        denominator += r * r
    if not denominator:
        return None
    return Fraction(numerator, d * denominator)


def finite_shell_scan() -> None:
    """Search stored exact minimizers for shell-sign and prefix behavior."""
    reports = []
    for name, a, ms in (
        ("A6", A6, (3, 4, 5)),
        ("A8", A8, (4, 5, 6, 7)),
        ("A9", A9, (5, 6, 7, 8)),
        ("A10", A, (5, 6, 7, 8, 9)),
    ):
        n = len(a)
        for m in ms:
            _, selectors, _, favorable, rows = build_relation(a, m)
            active = favorable.any(axis=1)
            caps = sorted(set(map(int, rows[active])))
            ell = 2
            d = math.comb(m, ell) * math.comb(n - ell, m - ell)
            lam2 = Fraction(
                falling(ell, 2) * falling(n - m, 2),
                falling(m, 2) * falling(n - ell, 2),
            ) if n - m >= 2 else Fraction(0)
            shells = []
            prefixes = []
            for cap in caps:
                raw = Fraction(0)
                for iz in np.flatnonzero((rows == cap) & active):
                    fam_ids = list(map(int, np.flatnonzero(favorable[iz])))
                    r = len(fam_ids)
                    overlap = sum(
                        math.comb(len(set(selectors[i]) & set(selectors[j])), ell)
                        for i in fam_ids for j in fam_ids
                    )
                    raw += Fraction(r * overlap, d) - lam2 * r * r
                shells.append(raw)
                p = retention_generic(favorable, rows, selectors, cap, ell)
                assert p is not None
                prefixes.append(p - lam2)
            first_nonnegative = next((caps[i] for i, x in enumerate(prefixes) if x >= 0), None)
            h2 = Fraction(1, math.comb(n - 2, m - 2))
            # For every nonempty family p_2 is at least the kernel's diagonal
            # atom h_2.  Record whether this alone makes shell positivity
            # automatic; instances beyond that regime are only finite data.
            reports.append((name, m, caps[0], first_nonnegative,
                            sum(x < 0 for x in shells), min(prefixes),
                            h2 - lam2, h2 >= lam2))
    assert all(neg == 0 for _, _, _, _, neg, _, _, _ in reports)
    assert all(first is not None for _, _, _, first, _, _, _, _ in reports)
    print("stored exact-minimizer level-two shell scan")
    print("name m mincap first_nonnegative_prefix negative_shells min_prefix_gap h2-lambda2 automatic")
    for row in reports:
        print(" ".join(map(str, row)))

    # At fixed density the diagonal atom is exponentially small whereas
    # lambda_2 is polynomial.  The first rounded 3/5 example already occurs
    # just beyond the stored orders.
    n, m = 12, 7
    h2 = Fraction(1, math.comb(n - 2, m - 2))
    lam2 = Fraction(2 * (n - m) * (n - m - 1),
                    m * (m - 1) * (n - 2) * (n - 3))
    assert h2 == Fraction(1, 252) < lam2 == Fraction(2, 189)
    print(f"self-loop crossover witness n=12,m=7: h2={h2}<lambda2={lam2}")


def shell_atlas() -> None:
    n, m, ell = 10, 6, 2
    centers, selectors, child_caps, favorable, rows = build_ground_relation(m)
    d = math.comb(m, ell) * math.comb(n - ell, m - ell)
    h2 = Fraction(1, math.comb(n - 2, m - 2))
    lam2 = Fraction(
        falling(ell, 2) * falling(n - m, 2),
        falling(m, 2) * falling(n - ell, 2),
    )
    active_caps = sorted(set(map(int, rows[favorable.any(axis=1)])))
    shell_integer = {}
    prefix = {}
    for cap in active_caps:
        ids = np.flatnonzero((rows == cap) & favorable.any(axis=1))
        raw = Fraction(0)
        for iz in ids:
            fam_ids = list(map(int, np.flatnonzero(favorable[iz])))
            r = len(fam_ids)
            overlap = sum(
                math.comb(len(set(selectors[i]) & set(selectors[j])), ell)
                for i in fam_ids for j in fam_ids
            )
            off_diagonal = overlap - r * math.comb(m, 2)
            pz = Fraction(overlap, r * d)
            assert pz == h2 + Fraction(off_diagonal, r * d)
            assert pz <= r * h2
            # N^2 times the shell contribution a_z^2(p_2-lambda_2).
            raw += Fraction(r * overlap, d) - lam2 * r * r
        shell_integer[cap] = raw
        p = aggregate_retention(favorable, rows, selectors, cap, ell)
        prefix[cap] = p - lam2

    assert active_caps == [10, 42, 74, 106, 138]
    assert shell_integer[10] > 0
    assert all(v > 0 for v in shell_integer.values())
    assert all(prefix[c] >= 0 for c in active_caps)
    print("A10,m6 level-two shell atlas")
    print(f"lambda2={lam2}, active_caps={active_caps}")
    for cap in active_caps:
        print(
            f"cap={cap:3d} shell_scaled={shell_integer[cap]} "
            f"prefix_gap={prefix[cap]}"
        )

    low = list(map(int, np.flatnonzero(rows == 10)))
    assert len(low) == 2
    for iz in low:
        fam_ids = list(map(int, np.flatnonzero(favorable[iz])))
        fam = [selectors[j] for j in fam_ids]
        assert len(fam) == 5
        adj_pairs = sum(adjacent(s, t) for s in fam for t in fam if s < t)
        assert adj_pairs == 0
        hist = Counter(len(set(s) & set(t)) for s in fam for t in fam)
        print(
            f"row10 center={tuple(map(int, centers[iz]))} "
            f"family_size={len(fam)} adjacent_pairs={adj_pairs} "
            f"ordered_intersections={dict(sorted(hist.items()))}"
        )

    # Examine every one-coordinate outside flip of every row-ten incidence.
    center_index = {tuple(map(int, z)): i for i, z in enumerate(centers)}
    flip_stats = Counter()
    adjacency_stats = Counter()
    for iz in low:
        z = centers[iz]
        for js in np.flatnonzero(favorable[iz]):
            s = selectors[int(js)]
            outside = [j for j in range(n) if j not in s]
            for j in outside:
                zz = z.copy()
                zz[j] *= -1
                if zz[0] == -1:
                    zz *= -1
                kz = center_index[tuple(map(int, zz))]
                flip_stats[int(rows[kz])] += 1
                fam2 = [selectors[u] for u in np.flatnonzero(favorable[kz])]
                adjacency_stats[sum(adjacent(s, t) for t in fam2)] += 1
    assert max(flip_stats) <= 42
    print(f"row10-incidence one-flip row histogram={dict(sorted(flip_stats.items()))}")
    print(
        "number of adjacent favorable selectors to seed after one outside flip="
        f"{dict(sorted(adjacency_stats.items()))}"
    )


def maximal_exchange_identity() -> None:
    """Verify the exact maximal-selector port inequalities on A10,m6."""
    n, m = 10, 6
    centers, selectors, child_caps, favorable, rows = build_ground_relation(m)
    qstar = max(child_caps)
    maximal = [j for j, q in enumerate(child_caps) if q == qstar]
    checks = 0
    equality_pairs = 0
    strict_pairs = 0
    no_equality_incidence = 0
    row_hist = Counter()
    positive_star_min = None
    selector_index = {s: j for j, s in enumerate(selectors)}
    center_index = {tuple(map(int, z)): i for i, z in enumerate(centers)}
    for iz, z0 in enumerate(centers):
        for js in maximal:
            if not favorable[iz, js]:
                continue
            s = selectors[js]
            ss = set(s)
            tset = tuple(j for j in range(n) if j not in ss)
            raw = energy(A[np.ix_(s, s)], z0[list(s)])
            sigma = 1 if raw > 0 else -1
            z = z0.copy()
            # Oriented signed edge products; sigma is kept explicit.
            r = {
                i: sigma * int(z[i] * sum(A[i, u] * z[u] for u in s if u != i))
                for i in s
            }
            assert sum(r.values()) == qstar
            cross = sigma * sum(A[i, j] * z[i] * z[j] for i in s for j in tset)
            gap_sum = 0
            sq_gap_sum = 0
            eq_here = 0
            assigned = defaultdict(set)
            good_assignments = 0
            threshold = Fraction(4 * qstar, m)
            for i in s:
                for j in tset:
                    tij = sigma * int(
                        z[j] * sum(A[j, u] * z[u] for u in s if u != i)
                    )
                    assert abs(tij) <= r[i]
                    assert (r[i] - tij) % 2 == 0
                    gap_sum += r[i] - tij
                    sq_gap_sum += r[i] * r[i] - tij * tij
                    if abs(tij) == r[i]:
                        equality_pairs += 1
                        eq_here += 1
                    else:
                        strict_pairs += 1

                    # Choose the sign of the incoming port that maximizes the
                    # exchanged candidate energy.  It is realized at z or at
                    # the one-coordinate outside flip z^j.
                    use_flip = (tij < 0)
                    zz = z.copy()
                    if use_flip:
                        zz[j] *= -1
                    if zz[0] == -1:
                        zz *= -1
                    kz = center_index[tuple(map(int, zz))]
                    sprime = tuple(sorted((ss - {i}) | {j}))
                    jsp = selector_index[sprime]
                    candidate = abs(energy(A[np.ix_(sprime, sprime)], zz[list(sprime)]))
                    deficit = child_caps[jsp] - candidate
                    assert 0 <= deficit <= 2 * (r[i] - abs(tij))
                    if Fraction(deficit) <= threshold:
                        assigned[kz].add(sprime)
                        good_assignments += 1
                    checks += 1
            if eq_here == 0:
                no_equality_incidence += 1
            l = sigma * z * (A @ z)
            assert gap_sum == (n - 1) * qstar - (m - 1) * int(sum(l[list(s)]))
            gint = sum(v * v for v in r.values())
            cross_vec = A[np.ix_(tset, s)] @ z[list(s)]
            gext = int(cross_vec @ cross_vec)
            assert sq_gap_sum == (n - m) * gint - (m - 2) * gext - m * (n - m)
            assert sq_gap_sum >= 0
            assert good_assignments >= m * (n - m) // 2
            star = max(len(v) for v in assigned.values())
            assert star >= math.ceil(good_assignments / (n - m + 1))
            positive_star_min = star if positive_star_min is None else min(positive_star_min, star)
            for kz in assigned:
                assert rows[kz] <= (math.sqrt(rows[iz]) + 2 * math.sqrt(n)) ** 2 + 1e-9
            row_hist[int(rows[iz])] += 1

    assert checks > 0
    print(
        f"maximal-selector qstar={qstar}, port_checks={checks}, "
        f"equality_pairs={equality_pairs}, strict_pairs={strict_pairs}, "
        f"incidences_without_any_equality={no_equality_incidence}"
    )
    print(
        f"threshold=4qstar/m={float(Fraction(4*qstar,m)):.6f}, "
        f"minimum assigned near-ground star={positive_star_min}"
    )
    print(f"maximal-ground incidence row histogram={dict(sorted(row_hist.items()))}")


def main() -> None:
    assert qnorm(A) == 26
    shell_atlas()
    maximal_exchange_identity()
    finite_shell_scan()
    print("PASS: exact level-two shells and maximal-selector exchange identities")


if __name__ == "__main__":
    main()
