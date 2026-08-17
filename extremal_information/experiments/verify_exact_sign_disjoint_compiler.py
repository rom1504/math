#!/usr/bin/env python3
"""Finite checks for the character-preserving lock ceiling draft."""

from __future__ import annotations

import itertools
import math
import random


def spins(k: int):
    return itertools.product((-1, 1), repeat=k)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mat_t_vec(R, x):
    return [sum(R[i][j] * x[i] for i in range(len(x))) for j in range(len(R[0]))]


def quad(R, x, y):
    return sum(x[i] * R[i][j] * y[j] for i in range(len(x)) for j in range(len(y)))


def check_pair_set_star(k: int) -> int:
    """Every k-1 pairwise-intersecting distinct 2-set family is a star, k>=5."""
    edges = list(itertools.combinations(range(k), 2))
    checked = 0
    for fam in itertools.combinations(edges, k - 1):
        if all(set(a) & set(b) for a, b in itertools.combinations(fam, 2)):
            checked += 1
            common = set(fam[0]).intersection(*map(set, fam[1:]))
            assert common
    return checked


def check_square_lock(k: int, trials: int, rng: random.Random) -> int:
    xs = list(spins(k))
    mu = sum(abs(sum(x)) for x in xs) / len(xs)
    assert mu + 1e-12 >= math.sqrt(k / 2)
    checked = 0
    for _ in range(trials):
        R = [[rng.choice((-1, 1)) for _ in range(k)] for _ in range(k)]
        perm = list(range(k))
        rng.shuffle(perm)
        signs = [rng.choice((-1, 1)) for _ in range(k)]
        defects = []
        roofs = []
        intended = []
        for x in xs:
            y = [signs[j] * x[perm[j]] for j in range(k)]
            roof = sum(abs(v) for v in mat_t_vec(R, x))
            q = abs(quad(R, x, y))
            assert roof >= q
            defects.append(roof - q)
            roofs.append(roof)
            intended.append(q)
        assert abs(sum(roofs) / len(xs) - k * mu) < 1e-9
        assert sum(q * q for q in intended) / len(xs) <= 3 * k * k + 1e-9
        assert max(defects) + 1e-9 >= k * math.sqrt(k / 2) - math.sqrt(3) * k
        checked += len(xs)
    return checked


def check_replication(k: int, copies: int, trials: int, rng: random.Random) -> int:
    m = copies * k
    pi = [j for j in range(k) for _ in range(copies)]
    xs = list(spins(k))
    checked = 0
    for _ in range(trials):
        R = [[rng.choice((-1, 1)) for _ in range(m)] for _ in range(k)]
        ss = [rng.choice((-1, 1)) for _ in range(m)]
        defects = []
        q2 = 0.0
        for x in xs:
            y = [ss[a] * x[pi[a]] for a in range(m)]
            roof = sum(abs(v) for v in mat_t_vec(R, x))
            q = abs(quad(R, x, y))
            assert roof >= q
            defects.append(roof - q)
            q2 += q * q
        q2 /= len(xs)
        assert q2 <= m * m + 2 * k * copies * m + 1e-9
        rhs = m * math.sqrt(k / 2) - math.sqrt(m * m + 2 * k * copies * m)
        assert max(defects) + 1e-9 >= rhs
        checked += len(xs)
    return checked


def check_robust_pin_rank_one(k: int) -> int:
    """Exhaust the CD.5 flip inequalities for u=1 at small k."""
    assert k <= 3
    count = 0
    for bits in itertools.product((-1, 1), repeat=k * k):
        R = [list(bits[i * k : (i + 1) * k]) for i in range(k)]
        u = (1,) * k
        f_u = sum(abs(v) for v in mat_t_vec(R, u))
        robust = True
        for i in range(k):
            ui = list(u)
            ui[i] = -1
            f_i = sum(abs(v) for v in mat_t_vec(R, ui))
            robust &= f_u - f_i >= 2 * (k - 1)
        if robust:
            assert all(R[i] == R[0] for i in range(1, k))
        count += 1
    return count


def main() -> None:
    rng = random.Random(20260817)
    count = check_pair_set_star(5)
    for k in range(2, 9):
        count += check_square_lock(k, 16, rng)
    for k in range(2, 8):
        count += check_replication(k, 2, 12, rng)
    count += check_robust_pin_rank_one(2)
    count += check_robust_pin_rank_one(3)
    print(f"exact-sign disjoint compiler checks passed: {count}")


if __name__ == "__main__":
    main()
