#!/usr/bin/env python3
"""Exact checks for the projective-cube Max-Cut profile construction.

This is a finite verifier, not evidence for the asymptotic theorem.  It checks
for small boundary widths that

* the pinning gadget has profile C-d_proj;
* the projective Hamming distance matrix is nonsingular; and
* its character eigenvalues agree with the closed formula in the proof.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random


def canonical(x: tuple[int, ...]) -> tuple[int, ...]:
    return x if x[0] == 1 else tuple(-a for a in x)


def projective_words(w: int) -> list[tuple[int, ...]]:
    return [(1,) + tail for tail in itertools.product((-1, 1), repeat=w - 1)]


def projective_distance(x: tuple[int, ...], y: tuple[int, ...]) -> int:
    d = sum(a != b for a, b in zip(x, y))
    return min(d, len(x) - d)


def gadget_score(sigma: tuple[int, ...], tau: tuple[int, ...]) -> int:
    """Maximize the unpadded pinning gadget over its anchor sign.

    tau_i=+1 is implemented by a two-edge equality path and tau_i=-1 by
    a direct inequality edge.  The path's internal spin has already been
    maximized out.
    """
    best = -1
    for anchor in (-1, 1):
        score = 0
        for s, t in zip(sigma, tau):
            if t == 1:
                score += 2 if s == anchor else 1
            else:
                score += 1 if s != anchor else 0
        best = max(best, score)
    return best


def formula_abs_eigenvalue(w: int, level: int) -> int:
    assert level % 2 == 0 and level > 0
    m = (w + 1) // 2 if w % 2 else w // 2
    j = level // 2
    numerator = math.factorial(2 * j - 2) * math.factorial(2 * m - 2 * j)
    denominator = (
        math.factorial(m - 1)
        * math.factorial(j - 1)
        * math.factorial(m - j)
    )
    value = numerator // denominator
    return value if w % 2 == 0 else value // 2


def character(word: tuple[int, ...], subset: tuple[int, ...]) -> int:
    out = 1
    for i in subset:
        out *= word[i]
    return out


def verify_width(w: int) -> None:
    words = projective_words(w)
    q = len(words)

    # The gadget constant depends on the representative before padding.
    for tau in words:
        c_tau = sum(2 if t == 1 else 1 for t in tau)
        for sigma in words:
            got = gadget_score(sigma, tau)
            want = c_tau - projective_distance(sigma, tau)
            assert got == want, (w, sigma, tau, got, want)

    distance_rows = [
        [projective_distance(x, y) for y in words]
        for x in words
    ]

    eigenvalues: list[int] = []
    for level in range(0, w + 1, 2):
        for subset in itertools.combinations(range(w), level):
            # Translation invariance makes it enough to evaluate D chi at
            # any row, then divide by the character at that row.  The first
            # representative in projective_words is not necessarily the
            # all-one identity.
            row_sum = sum(
                distance_rows[0][j] * character(words[j], subset)
                for j in range(q)
            )
            lam = row_sum // character(words[0], subset)
            eigenvalues.append(lam)
            if level > 0:
                assert abs(lam) == formula_abs_eigenvalue(w, level), (
                    w,
                    level,
                    lam,
                    formula_abs_eigenvalue(w, level),
                )

    assert len(eigenvalues) == q
    assert all(lam != 0 for lam in eigenvalues)
    print(
        f"w={w:2d} q={q:5d} "
        f"min_abs_eigenvalue={min(abs(x) for x in eigenvalues):5d} verified"
    )


def verify_lookup_universality(max_width: int, seed: int = 20260816) -> int:
    rng = random.Random(seed)
    checks = 0
    for w in range(1, min(max_width, 4) + 1):
        oriented = list(itertools.product((-1, 1), repeat=w))
        reps = projective_words(w)
        for _ in range(20):
            table = {rep: rng.randint(0, 7) for rep in reps}
            lifted = {a: table[canonical(a)] for a in oriented}
            for sigma in oriented:
                best = -math.inf
                for anchor in (-1, 1):
                    # Maximizing each t_a is the same as choosing y_a in {0,1}.
                    score = sum(
                        lifted[a]
                        * max(
                            0,
                            sum(ai * si * anchor for ai, si in zip(a, sigma))
                            - (w - 1),
                        )
                        for a in oriented
                    )
                    best = max(best, score)
                assert best == table[canonical(sigma)]
                checks += 1

            # Keeping every signed pair occurrence separate makes the
            # Max-Cut conversion offset depend only on sum(table).
            predicted_offset = (6 * w - 2) * sum(table.values())
            direct_offset = 0
            for a in oriented:
                lam = lifted[a]
                original_constant = -(w - 1) * lam / 2
                replacement_constant = 0.0
                for sign in a:
                    coupling = lam * sign / 2
                    replacement_constant += 2 * (
                        3 * coupling if coupling >= 0 else -coupling
                    )
                final_coupling = -(w - 1) * lam / 2
                replacement_constant += (
                    3 * final_coupling if final_coupling >= 0 else -final_coupling
                )
                direct_offset += replacement_constant - original_constant
            assert direct_offset == predicted_offset
            checks += 1
    return checks


def verify_explicit_maxcut_lookup(seed: int = 1618033) -> int:
    """Build the actual positive-edge graph and exhaust it through width two."""
    rng = random.Random(seed)
    checks = 0
    for w in (1, 2):
        reps = projective_words(w)
        oriented = list(itertools.product((-1, 1), repeat=w))
        for _ in range(6):
            table = {rep: rng.randint(0, 4) for rep in reps}
            lifted = {a: table[canonical(a)] for a in oriented}
            # Boundary vertices are 0,...,w-1; all remaining vertices are
            # private.  Pair terms are kept separate, as in the proof.
            next_vertex = w
            anchor = next_vertex
            next_vertex += 1
            selector = {}
            for a in oriented:
                selector[a] = next_vertex
                next_vertex += 1
            edges: list[tuple[int, int, int]] = []

            def implement_pair(u: int, v: int, twice_coupling: int) -> None:
                nonlocal next_vertex
                # twice_coupling=2J.  The direct edge or each path edge has
                # weight |2J|, so integer tables produce integer weights.
                if twice_coupling < 0:
                    edges.append((u, v, -twice_coupling))
                elif twice_coupling > 0:
                    middle = next_vertex
                    next_vertex += 1
                    edges.append((u, middle, twice_coupling))
                    edges.append((middle, v, twice_coupling))

            for a in oriented:
                lam = lifted[a]
                for i, sign in enumerate(a):
                    implement_pair(i, anchor, lam * sign)
                    implement_pair(selector[a], i, lam * sign)
                implement_pair(selector[a], anchor, -(w - 1) * lam)

            offset = (6 * w - 2) * sum(table.values())
            private_count = next_vertex - w
            for sigma in itertools.product((-1, 1), repeat=w):
                best = -1
                for private in itertools.product((-1, 1), repeat=private_count):
                    spins = sigma + private
                    score = sum(weight for u, v, weight in edges if spins[u] != spins[v])
                    best = max(best, score)
                assert best == table[canonical(sigma)] + offset
                checks += 1
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-width", type=int, default=10)
    args = parser.parse_args()
    for w in range(2, args.max_width + 1):
        verify_width(w)
    print(f"lookup_universality_checks={verify_lookup_universality(args.max_width)} verified")
    print(f"explicit_maxcut_lookup_checks={verify_explicit_maxcut_lookup()} verified")


if __name__ == "__main__":
    main()
