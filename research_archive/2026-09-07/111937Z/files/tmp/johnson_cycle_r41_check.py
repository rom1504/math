#!/usr/bin/env python3
"""Exact checks for Wave 41 Johnson-cycle averaging."""

from __future__ import annotations

import itertools
from fractions import Fraction

import numpy as np

from anchored_conflict_r40_check import (
    A5,
    complement,
    cut_matrix,
    enumerate_q_min,
    projective_spins,
    q_value,
)


def anchored_sets(n: int, m: int, v: int = 0):
    return [tuple((v,) + c) for c in itertools.combinations(
        [i for i in range(n) if i != v], m - 1)]


def lex_ground_map(a: np.ndarray, m: int, v: int = 0):
    states = sorted((sig, tuple(map(int, x)))
                    for sig in (-1, 1) for x in projective_spins(len(a)))
    answer = {}
    for s0 in anchored_sets(len(a), m, v):
        b = complement(a, s0)
        cap = q_value(b)
        sig, xt = next((sig, xt) for sig, xt in states
                       if sig * int(np.asarray(xt) @ b @ np.asarray(xt)) == cap)
        x = np.asarray(xt, dtype=np.int64)
        answer[s0] = (sig, x, cut_matrix(x, sig))
    return answer


def corrected_conflict(mapping, n: int, m: int, v: int = 0) -> Fraction:
    sets = list(mapping)
    p = Fraction(m - 1, n - 1)
    total = Fraction(0)
    for s0 in sets:
        for t0 in sets:
            x, y = mapping[s0][1], mapping[t0][1]
            total += sum(x[i] != y[i] for i in (set(s0) & set(t0)) - {v}) / p
    return total / len(sets) ** 2


def has_lex_conflict(a: np.ndarray, m: int, v: int = 0) -> bool:
    mapping = lex_ground_map(a, m, v)
    return corrected_conflict(mapping, len(a), m, v) > 0


def no_smaller_lex_example(n: int) -> bool:
    qn = enumerate_q_min(n)
    edges = list(itertools.combinations(range(n), 2))
    for signs in itertools.product((-1, 1), repeat=len(edges)):
        a = np.zeros((n, n), dtype=np.int64)
        for (i, j), sign in zip(edges, signs):
            a[i, j] = a[j, i] = sign
        if q_value(a) == qn and has_lex_conflict(a, n - 1):
            return False
    return True


def audit_a5() -> dict[str, object]:
    n, m, v = 5, 4, 0
    mapping = lex_ground_map(A5, m, v)
    gaps = []
    overlap_changes = []
    for s0, t0 in itertools.combinations(mapping, 2):
        ds, dt = mapping[s0][2], mapping[t0][2]
        gap = int(np.sum((complement(A5, s0) - complement(A5, t0)) * (ds - dt)))
        h = sum(mapping[s0][1][i] != mapping[t0][1][i]
                for i in (set(s0) & set(t0)) - {v})
        assert gap >= 0 and gap % 8 == 0
        if mapping[s0][0] != mapping[t0][0] or not np.array_equal(
                mapping[s0][1], mapping[t0][1]):
            assert gap >= 8
        gaps.append(gap)
        overlap_changes.append(h)
    assert q_value(A5) == enumerate_q_min(5) == 8
    assert corrected_conflict(mapping, n, m, v) == 1
    assert no_smaller_lex_example(3)
    assert no_smaller_lex_example(4)
    return {"conflict": 1, "gaps": tuple(gaps), "overlap": tuple(overlap_changes)}


def balanced_wall(h: int):
    """Canonical two-response wall; h=1 mod 4, n=2h, m=n-2."""
    assert h % 4 == 1 and h >= 5
    n, v = 2 * h, 0
    kset = set(range(1, h + 1))
    llist = [0] + list(range(h + 1, n))

    # A bipartite +/-1 matrix with every row and column sum +/-1.
    q = (h - 1) // 2
    degrees = [q + 1] * ((h + 1) // 2) + [q] * ((h - 1) // 2)
    remaining = degrees.copy()
    cross = {}
    for ii, i in enumerate(sorted(kset)):
        chosen = sorted(range(h), key=lambda j: (-remaining[j], j))[:degrees[ii]]
        chosen = set(chosen)
        for jj, j in enumerate(llist):
            cross[tuple(sorted((i, j)))] = 1 if jj in chosen else -1
        for jj in chosen:
            remaining[jj] -= 1
    assert remaining == [0] * h

    a = np.zeros((n, n), dtype=np.int64)
    for (i, j), sign in cross.items():
        a[i, j] = a[j, i] = sign

    # Balanced circulant signings inside both odd blocks.
    radius = (h - 1) // 4
    for block in (sorted(kset), llist):
        for ii, i in enumerate(block):
            for jj in range(ii + 1, h):
                j = block[jj]
                dist = min((jj - ii) % h, (ii - jj) % h)
                sign = 1 if dist <= radius else -1
                a[i, j] = a[j, i] = sign

    x = np.ones(n, dtype=np.int64)
    y = x.copy()
    y[list(kset)] = -1
    dx, dy = cut_matrix(x, 1), cut_matrix(y, 1)
    assert int((a @ x) @ (a @ x)) == n
    assert int((a @ y) @ (a @ y)) == n

    delta = a * (dx - dy)
    degree = delta.sum(axis=1)
    assert set(map(int, degree)) == {-2, 2}
    u = sum(int(delta[i, j]) for i in range(n) for j in range(i + 1, n))
    assert u == 2

    universe = list(range(1, n))
    records = {}
    for omitted in itertools.combinations(universe, 2):
        s0 = tuple(i for i in range(n) if i not in omitted)
        residual = -4 * sum(int(degree[o]) for o in omitted) + 4 * int(delta[omitted])
        # P(dx)-P(dy)=2u=4; ties use one fixed order preferring dy.
        chosen = dx if residual > 0 else dy
        label = x if residual > 0 else y
        records[omitted] = (s0, residual, chosen, label)
        assert residual in {-24, -16, -8, 0, 8, 16, 24}

    # At h=5, exhaust all oriented cuts and verify the common penalty rule.
    if h == 5:
        other_penalty = 4 * n * (n - 1) + 10
        for omitted, (s0, residual, chosen, _) in records.items():
            b = complement(a, s0)
            scored = []
            for sig in (-1, 1):
                for z in projective_spins(n):
                    dz = cut_matrix(z, sig)
                    if np.array_equal(dz, dy):
                        penalty, rank = 0, 0
                    elif np.array_equal(dz, dx):
                        penalty, rank = 2 * u, 1
                    else:
                        penalty, rank = other_penalty, 2
                    scored.append((sig * int(z @ b @ z) - penalty, -rank, dz))
            best = max(scored, key=lambda item: (item[0], item[1]))[2]
            assert np.array_equal(best, chosen)

    # Directed adjacent gaps and overlap changes.
    gap_sum = h_sum = edge_count = 0
    for omitted, (s0, residual, ds, xs) in records.items():
        oset = set(omitted)
        for keep in omitted:
            for new in universe:
                if new in oset:
                    continue
                target = tuple(sorted((keep, new)))
                t0, residual_t, dt, xt = records[target]
                gap = int(np.sum((complement(a, s0) - complement(a, t0)) * (ds - dt)))
                assert gap >= 0 and gap % 8 == 0 and gap <= 48
                if not np.array_equal(ds, dt):
                    assert gap == abs(residual - residual_t) >= 8
                gap_sum += gap
                h_sum += sum(xs[i] != xt[i] for i in (set(s0) & set(t0)) - {v})
                edge_count += 1

    # Corrected conflict by exact integer summation.
    p = Fraction(n - 3, n - 1)
    total = Fraction(0)
    for s0, _, _, xs in records.values():
        for t0, _, _, xt in records.values():
            total += sum(xs[i] != xt[i] for i in (set(s0) & set(t0)) - {v}) / p
    conflict = total / len(records) ** 2

    # Every elementary square obeys forward/reverse cyclicity, and its
    # symmetric part is exactly the sum of its four pair gaps.
    for a0, b0, c0, d0 in itertools.permutations(universe, 4):
        cycle_o = ((b0, d0), (a0, d0), (a0, c0), (b0, c0))
        cycle = [records[tuple(sorted(o))] for o in cycle_o]
        forward = reverse = pair_sum = 0
        for j, (s0, _, ds, _) in enumerate(cycle):
            sn, _, dn, _ = cycle[(j + 1) % 4]
            sp, _, dp, _ = cycle[(j - 1) % 4]
            forward += int(np.sum(complement(a, s0) * (ds - dn)))
            reverse += int(np.sum(complement(a, s0) * (ds - dp)))
            pair_sum += int(np.sum((complement(a, s0) - complement(a, sn)) * (ds - dn)))
        assert 0 <= forward <= 96 and 0 <= reverse <= 96
        assert forward + reverse == pair_sum

    return {
        "n": n,
        "conflict": conflict,
        "mean_adjacent_gap": Fraction(gap_sum, edge_count),
        "mean_adjacent_overlap_change": Fraction(h_sum, edge_count),
        "row": n,
    }


def main() -> None:
    print("A5 canonical", audit_a5())
    for h in (5,):
        print("wall", balanced_wall(h))
    print("PASS johnson_cycle_r41_check")


if __name__ == "__main__":
    main()
