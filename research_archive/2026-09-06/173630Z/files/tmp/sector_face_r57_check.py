#!/usr/bin/env python3
"""Exact checks for the Wave 57 sector-disagreement identity and A9 wall."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def cap_and_grounds(a: np.ndarray):
    records = [(x, int(x @ a @ x)) for x in projective_spins(len(a))]
    q = max(abs(v) for _, v in records)
    grounds = [
        (1 if v > 0 else -1, x)
        for x, v in records
        if abs(v) == q
    ]
    return q, grounds


def canonical_excess_block(a: np.ndarray, selector: np.ndarray):
    """Lex-first child ground and lex-first r_i positive edges at each i."""
    child = a[np.ix_(selector, selector)]
    q_s, grounds = cap_and_grounds(child)
    sigma, y = grounds[0]
    fields = sigma * y * (child @ y)
    assert np.all(fields >= 0)
    assert int(fields.sum()) == q_s

    edge_set: set[tuple[int, int]] = set()
    for ii, i in enumerate(selector):
        positive = []
        for jj, j in enumerate(selector):
            if ii == jj:
                continue
            if sigma * y[ii] * child[ii, jj] * y[jj] == 1:
                positive.append((int(i), int(j)))
        assert len(positive) == (len(selector) - 1 + int(fields[ii])) // 2
        for edge in positive[: int(fields[ii])]:
            edge_set.add(tuple(sorted(edge)))
    edges = sorted(edge_set)
    assert q_s / 2 <= len(edges) <= q_s
    return q_s, sigma, y, fields, edges


def all_legal_unions(child: np.ndarray, sigma: int, y: np.ndarray):
    """All unions obtained by choosing r_i positive incident edges per i."""
    m = len(child)
    edges = [(i, j) for i in range(m) for j in range(i + 1, m)]
    edge_index = {edge: k for k, edge in enumerate(edges)}
    c = {
        (i, j): int(sigma * y[i] * child[i, j] * y[j])
        for i, j in edges
    }
    fields = sigma * y * (child @ y)
    unions = {0}
    for i in range(m):
        positive_indices = []
        for j in range(m):
            if i == j:
                continue
            edge = tuple(sorted((i, j)))
            if c[edge] == 1:
                positive_indices.append(edge_index[edge])
        choices = []
        for choice in itertools.combinations(positive_indices, int(fields[i])):
            choices.append(sum(1 << k for k in choice))
        unions = {old | new for old in unions for new in choices}
    return edges, unions


def sector_identity_audit():
    # Stored exact minimizer and one exact negative-sector parent ground.
    a = A9
    n, m = 9, 7
    q, parent_grounds = cap_and_grounds(a)
    assert q == 24
    tau = -1
    x = np.asarray((1, -1, 1, 1, -1, -1, -1, -1, -1), dtype=np.int64)
    assert any(t == tau and np.array_equal(u, x) for t, u in parent_grounds)
    assert tau * int(x @ a @ x) == q
    assert int((a @ x) @ (a @ x)) == 96

    all_edge_sum = sum(
        tau * int(a[i, j]) * int(x[i]) * int(x[j])
        for i in range(n)
        for j in range(i + 1, n)
    )
    assert all_edge_sum == q // 2 == 12

    # The orientation-conditional correlation matrix is already rank one.
    correlation = np.outer(x, x)
    assert np.array_equal(np.diag(correlation), np.ones(n, dtype=np.int64))
    assert np.linalg.matrix_rank(correlation) == 1
    assert np.linalg.eigvalsh(correlation)[0] >= -1e-12

    escaped = []
    favorable = []
    orientation_count = Counter()
    deficit_count = Counter()
    sector_deficit_count = Counter()
    b_constant = 56 * math.sqrt(7) / 9 - 14
    assert 0 < b_constant < 4

    for selector0 in itertools.combinations(range(n), m):
        selector = np.asarray(selector0, dtype=int)
        q_s, sigma, y, _, edges = canonical_excess_block(a, selector)
        local = tau * int(x[selector] @ a[np.ix_(selector, selector)] @ x[selector])
        delta = q_s - local
        assert delta >= 0

        kappa = tau * sigma
        z = x[selector] * y
        position = {int(vertex): i for i, vertex in enumerate(selector)}
        disagreement = []
        c_sum = 0
        selected_sum = 0
        selected_disagreement = 0

        for ii in range(m):
            for jj in range(ii + 1, m):
                i, j = int(selector[ii]), int(selector[jj])
                c_e = int(sigma * y[ii] * a[i, j] * y[jj])
                s_e = int(tau * a[i, j] * x[i] * x[j])
                assert s_e == c_e * kappa * int(z[ii]) * int(z[jj])
                if kappa * int(z[ii]) * int(z[jj]) == -1:
                    disagreement.append((i, j))
                    c_sum += c_e

        for i, j in edges:
            ii, jj = position[i], position[j]
            c_e = int(sigma * y[ii] * a[i, j] * y[jj])
            assert c_e == 1
            s_e = int(tau * a[i, j] * x[i] * x[j])
            selected_sum += s_e
            if kappa * int(z[ii]) * int(z[jj]) == -1:
                selected_disagreement += 1

        # Unified identity, with all factors of two checked exactly.
        assert delta == 4 * c_sum
        assert selected_sum == len(edges) - 2 * selected_disagreement
        is_escape = 2 * selected_sum <= len(edges)
        assert is_escape == (4 * selected_disagreement >= len(edges))

        # At t=0 and Delta=0, B=(56 sqrt(7))/9-14 lies strictly in (0,4).
        widehat_ell = delta - b_constant
        is_favorable = widehat_ell <= 0
        assert is_favorable == (delta == 0)
        if is_escape:
            assert widehat_ell > 0
            escaped.append(selector0)
            orientation_count[kappa] += 1
            deficit_count[delta] += 1
            sector_deficit_count[(kappa, delta)] += 1
        if is_favorable:
            favorable.append(selector0)
        assert not (is_escape and is_favorable)

        # Flipping precisely the disagreement set raises this parent's energy.
        disagreement_sum = sum(
            tau * int(a[i, j]) * int(x[i]) * int(x[j])
            for i, j in disagreement
        )
        assert disagreement_sum == -delta // 4
        assert q - 4 * disagreement_sum == q + delta

    assert 56 * 56 * 7 > 14 * 14 * 9 * 9  # B>0, exact squared comparison.
    assert 56 * 56 * 7 < 18 * 18 * 9 * 9  # B<4, exact squared comparison.
    assert len(escaped) == 28
    assert len(favorable) == 5
    assert orientation_count == Counter({1: 17, -1: 11})
    assert deficit_count == Counter({8: 14, 4: 8, 12: 3, 16: 3})
    assert sector_deficit_count == Counter(
        {(-1, 4): 8, (-1, 12): 3, (1, 8): 14, (1, 16): 3}
    )

    print(
        "canonical A9,m=7 obstruction: parent energy=24 row=96; "
        "escape=28/36, favorable=5/36, intersection=0"
    )
    print("  escaped orientation sectors", dict(orientation_count))
    print("  escaped local deficits", dict(sorted(deficit_count.items())))
    print("  sector-by-deficit", dict(sorted(sector_deficit_count.items())))
    print(
        "  escaped hat-ell values",
        {delta: delta - b_constant for delta in sorted(deficit_count)},
    )


def choice_robust_audit():
    """The Fubini-size bad subfamily is not an artifact of lex choices."""
    a = A9
    n, m = 9, 7
    tau = -1
    x = np.asarray((1, -1, 1, 1, -1, -1, -1, -1, -1), dtype=np.int64)
    robust_bad = 0

    for selector0 in itertools.combinations(range(n), m):
        selector = np.asarray(selector0, dtype=int)
        child = a[np.ix_(selector, selector)]
        q_s, grounds = cap_and_grounds(child)
        local = tau * int(x[selector] @ child @ x[selector])
        delta = q_s - local
        every_choice_escapes = True

        for sigma, y in grounds:
            local_edges, unions = all_legal_unions(child, sigma, y)
            for mask in unions:
                count = bin(mask).count("1")
                selected_sum = 0
                for edge_index, (ii, jj) in enumerate(local_edges):
                    if (mask >> edge_index) & 1:
                        i, j = int(selector[ii]), int(selector[jj])
                        selected_sum += tau * int(a[i, j]) * int(x[i]) * int(x[j])
                if 2 * selected_sum > count:
                    every_choice_escapes = False
                    break
            if not every_choice_escapes:
                break

        if every_choice_escapes:
            assert delta >= 4
            robust_bad += 1

    assert robust_bad == 15
    assert robust_bad > 36 / 3
    print("choice-robust bad escapes: 15/36 (>1/3) for every legal child/block choice")


if __name__ == "__main__":
    sector_identity_audit()
    choice_robust_audit()
    print("PASS sector identity, PSD/full-energy audit, and exact finite wall")
