#!/usr/bin/env python3
"""Finite checks for drafts/approximate_residual_shell_law.md.

The script checks:
  * the tropical-hull estimate behind ARS.1 on random exact data;
  * the exact rank-one cyclic compatibility formula;
  * the sharp delta/4 one-state rate in ARS.3 through exhaustive words;
  * the uniform alpha row shell in ARS.4;
  * the mandatory empty-support-core/nonzero-semantic distinction.
"""

from __future__ import annotations

import itertools
import math
import random


NEG = -10**100


def mp_mul(a, b):
    n, k, m = len(a), len(b), len(b[0])
    return [
        [max(a[i][h] + b[h][j] for h in range(k)) for j in range(m)]
        for i in range(n)
    ]


def mp_word(mats, word):
    out = mats[word[0]]
    for e in word[1:]:
        out = mp_mul(out, mats[e])
    return out


def max_cycle_mean_2x2(a):
    return max(a[0][0], a[1][1], (a[0][1] + a[1][0]) / 2)


def rho_word(mats, word):
    return max_cycle_mean_2x2(mp_word(mats, word))


def rank_one(left, right):
    return [[x + y for y in right] for x in left]


def cyclic_compat(left, right, word):
    ans = 0.0
    for e, f in zip(word, word[1:] + word[:1]):
        ans += max(right[e][j] + left[f][j] for j in range(len(left[e])))
    return ans


def check_tropical_hull(rng):
    checks = 0
    for _ in range(200):
        n = 4
        p = [rng.uniform(-2, 2) for _ in range(n)]
        left = [rng.uniform(-2, 2) for _ in range(n)]
        eps = rng.uniform(0.001, 0.4)
        v = [
            [left[i] + p[j] + rng.uniform(-eps, eps) for j in range(n)]
            for i in range(n)
        ]
        prefix = [[rng.uniform(-3, 3) for _ in range(n)] for _ in range(n)]
        prod = mp_mul(prefix, v)
        for i in range(n):
            b = max(prefix[i][k] + left[k] for k in range(n))
            err = max(abs(prod[i][j] - b - p[j]) for j in range(n))
            assert err <= eps + 1e-10
            checks += 1
            for _ in range(5):
                z = [rng.uniform(-10, 10) for _ in range(n)]
                roof_row = max(prod[i][j] + z[j] for j in range(n)) - max(prod[i])
                roof_p = max(p[j] + z[j] for j in range(n)) - max(p)
                assert abs(roof_row - roof_p) <= 2 * eps + 1e-10
                checks += 1
    return checks


def check_rank_one_counterexample():
    delta = 1.0
    left = {"A": [0.0, 0.0], "B": [delta, 0.0]}
    right = {"A": [0.0, 0.0], "B": [0.0, delta]}
    mats = {e: rank_one(left[e], right[e]) for e in left}
    assert mats["A"] == [[0.0, 0.0], [0.0, 0.0]]
    assert mats["B"] == [[1.0, 2.0], [0.0, 1.0]]

    checks = 0
    worst = 0.0
    for t in range(1, 11):
        for bits in itertools.product("AB", repeat=t):
            word = "".join(bits)
            actual = rho_word(mats, word)
            formula = cyclic_compat(left, right, word)
            assert abs(actual - formula) < 1e-10
            naa = sum(word[i] == "A" and word[(i + 1) % t] == "A" for i in range(t))
            assert abs(actual - (t - naa)) < 1e-10
            prediction = 0.25 * word.count("A") + 1.25 * word.count("B")
            err_rate = abs(actual - prediction) / t
            assert err_rate <= 0.25 + 1e-10
            worst = max(worst, err_rate)
            checks += 3
    assert abs(worst - 0.25) < 1e-10

    # The three forced cycles give the matching LP lower bound.  A fine grid
    # independently locates its minimum at (1/4,5/4).
    grid_best = 100.0
    arg = None
    for ix in range(-40, 81):
        x = ix / 40
        for iy in range(0, 81):
            y = iy / 40
            val = max(abs(x), abs(y - 1), abs((x + y) / 2 - 1))
            if val < grid_best:
                grid_best, arg = val, (x, y)
    assert abs(grid_best - 0.25) < 1e-12
    assert arg == (0.25, 1.25)
    return checks + 2


def check_small_shell(rng):
    checks = 0
    alpha = 0.37
    mats = {
        e: [[alpha * rng.uniform(-1, 0) for _ in range(3)] for _ in range(3)]
        for e in "abc"
    }
    for t in range(1, 7):
        for word in itertools.product("abc", repeat=t):
            p = mp_word(mats, "".join(word))
            for row in p:
                assert max(row) - min(row) <= alpha + 1e-10
                midpoint = (max(row) + min(row)) / 2
                assert max(abs(x - midpoint) for x in row) <= alpha / 2 + 1e-10
            checks += 1
    return checks


def check_empty_core_guardrail():
    p = [0.0, -1.0]
    ta = [[0.0, -1.0], [-2.0, -3.0]]
    tb = [[-2.0, -3.0], [1.0, 0.0]]

    def row_times(row, mat):
        return [max(row[i] + mat[i][j] for i in range(2)) for j in range(2)]

    assert row_times(p, ta) == p
    assert row_times(p, tb) == p
    ra = {(i, j) for i in range(2) for j in range(2) if ta[i][j] >= 0}
    rb = {(i, j) for i in range(2) for j in range(2) if tb[i][j] >= 0}
    assert ra == {(0, 0)}
    assert rb == {(1, 0), (1, 1)}

    # Descending one-context core: K must be contained in both predecessor
    # images.  Starting at I, the first intersection is {0}; b then has no
    # source in {0}, so the next iterate is empty.
    k = {0, 1}
    for _ in range(3):
        image_a = {j for i, j in ra if i in k}
        image_b = {j for i, j in rb if i in k}
        k = k & image_a & image_b
    assert not k

    mats = {"a": ta, "b": tb}
    for t in range(1, 8):
        for word in itertools.product("ab", repeat=t):
            assert abs(rho_word(mats, "".join(word))) < 1e-10
    return 2 ** 8 - 2


def main():
    rng = random.Random(1729)
    checks = 0
    checks += check_tropical_hull(rng)
    checks += check_rank_one_counterexample()
    checks += check_small_shell(rng)
    checks += check_empty_core_guardrail()
    print(f"approximate residual-shell checks passed: {checks}")


if __name__ == "__main__":
    main()
