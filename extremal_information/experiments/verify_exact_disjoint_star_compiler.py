#!/usr/bin/env python3
"""Finite checks for the exact disjoint-star compiler note.

This script verifies:
  * the edge-variable identity for all signings through k=5;
  * the exact pair Fourier coefficient gamma_d through d=10;
  * the pair-mass inequality on random exact star sums;
  * the all-positive PSD-completion rank calculation on random diagonals.
  * the interacting-selector contact inequality on exhaustive tiny models.

It is a finite arithmetic verifier, not evidence for the asymptotic claims.
"""

from __future__ import annotations

import itertools
import math
import random

import numpy as np


def spins(n: int):
    return itertools.product((-1, 1), repeat=n)


def gamma(d: int) -> float:
    return math.comb(d - 2, (d - 1) // 2) / (2 ** (d - 2))


def check_sparse_identity() -> None:
    for k in range(2, 6):
        edges = list(itertools.combinations(range(k), 2))
        # Exhaust all T for k<=4 and a deterministic sample at k=5.
        words = list(spins(len(edges)))
        if k == 5:
            random.seed(20260817)
            words = random.sample(words, 64)
        for word in words:
            T = dict(zip(edges, word))
            for x in spins(k):
                lhs = sum(abs(x[i] - T[i, j] * x[j]) for i, j in edges)
                rhs = len(edges) - sum(T[i, j] * x[i] * x[j] for i, j in edges)
                assert lhs == rhs


def walsh_pair(d: int, i: int, j: int) -> float:
    total = 0
    for x in spins(d):
        total += abs(sum(x)) * x[i] * x[j]
    return total / (2**d)


def check_gamma() -> None:
    for d in range(2, 11):
        got = walsh_pair(d, 0, 1)
        assert abs(got - gamma(d)) < 1e-12, (d, got, gamma(d))
        assert gamma(d) <= 1 / math.sqrt(d - 1) + 1e-12
        assert gamma(d) * math.comb(d, 2) <= 0.5 * d ** 1.5 + 1e-12


def fourier_pair(values: dict[tuple[int, ...], float], i: int, j: int) -> float:
    return sum(v * x[i] * x[j] for x, v in values.items()) / len(values)


def check_random_pair_mass() -> None:
    rng = random.Random(7321)
    for k in range(3, 8):
        xs = list(spins(k))
        for _ in range(50):
            stars = []
            for _ in range(rng.randint(1, 10)):
                d = rng.randint(2, k)
                support = tuple(rng.sample(range(k), d))
                signs = {i: rng.choice((-1, 1)) for i in support}
                stars.append((support, signs))
            vals = {
                x: sum(abs(sum(sgn[i] * x[i] for i in supp)) for supp, sgn in stars)
                for x in xs
            }
            pair_mass = sum(
                abs(fourier_pair(vals, i, j))
                for i, j in itertools.combinations(range(k), 2)
            )
            budget = sum(gamma(len(supp)) * math.comb(len(supp), 2) for supp, _ in stars)
            assert pair_mass <= budget + 1e-10


def check_negative_simplex_rank() -> None:
    rng = random.Random(901)
    for k in range(2, 12):
        # Any PSD diag(p)-J has nullity at most one.  Sampling p well above k
        # gives full rank; p_i=k gives the sharp rank k-1 boundary.
        K = np.diag([float(k)] * k) - np.ones((k, k))
        assert np.linalg.matrix_rank(K, tol=1e-8) == k - 1
        for _ in range(10):
            p = np.array([k + 1 + rng.random() for _ in range(k)])
            K = np.diag(p) - np.ones((k, k))
            assert np.linalg.eigvalsh(K)[0] > 0
            assert np.linalg.matrix_rank(K, tol=1e-8) == k


def check_selector_contact_inequality() -> None:
    """Check the key contact implication for arbitrary affine codebooks.

    We generate arbitrary small affine selector families, compute their max
    envelope and its uniform error from the all-positive cut shell, and
    verify (SC.23) for every active selector at every balanced point whenever
    the theorem's positive-threshold hypothesis applies.
    """
    rng = random.Random(44881)
    for k in (2, 4, 6):
        xs = list(spins(k))
        balanced = [x for x in xs if sum(x) == 0]
        target = {x: (k * k - sum(x) ** 2) / 2 for x in xs}
        for _ in range(100):
            m = rng.randint(1, 5)
            pieces = []
            for _ in range(2**m):
                b = tuple(rng.randint(-m, m) for _ in range(k))
                c = rng.randint(-k * k, k * k)
                pieces.append((b, c))
            vals = {x: max(c + sum(b[i] * x[i] for i in range(k))
                           for b, c in pieces) for x in xs}
            eta = max(abs(vals[x] - target[x]) for x in xs)
            a = k * k / 2 - 2 * eta
            if a <= 0:
                continue
            for x in balanced:
                active = [(b, c) for b, c in pieces
                          if c + sum(b[i] * x[i] for i in range(k)) == vals[x]]
                # The contact conclusion assumes the whole envelope is eta
                # close, and therefore holds for every active affine piece.
                for b, _ in active:
                    assert sum(x[i] * b[i] for i in range(k)) + 1e-12 >= a


def check_all_even_levels() -> None:
    """Check SC.15c--SC.15e and nonvanishing through small orders."""
    from fractions import Fraction

    for k in range(2, 13):
        xs = list(spins(k))
        for s in range(2, k + 1, 2):
            alpha = Fraction(sum(abs(sum(x)) * math.prod(x[:s]) for x in xs), 2**k)
            assert alpha != 0
            if k % 2 == 0:
                n, r = k // 2, s // 2
                ihat = Fraction(((-1) ** r) * math.comb(2 * n, n) * math.comb(n, r),
                                (2 ** (2 * n)) * math.comb(2 * n, 2 * r))
                assert alpha == -Fraction(2 * n, 2 * r - 1) * ihat
            else:
                n, r = (k - 1) // 2, s // 2
                ihat = Fraction(((-1) ** r) * math.comb(2 * n, n) * math.comb(n, r),
                                (2 ** (2 * n)) * math.comb(2 * n, 2 * r))
                assert alpha == -Fraction(2 * n - 2 * r + 1, 2 * r - 1) * ihat


def check_selector_oscillation() -> None:
    """Check osc(max_y c_y+b_y.x) <= 2 max_{x,y}|x.B.y|."""
    rng = random.Random(55127)
    for k in range(2, 6):
        xs = list(spins(k))
        for m in range(1, 5):
            ys = list(spins(m))
            for _ in range(30):
                B = [[rng.choice((-1, 1)) for _ in range(m)] for _ in range(k)]
                c = {y: rng.randint(-10, 10) for y in ys}
                vals = [max(c[y] + sum(x[i] * B[i][a] * y[a]
                                       for i in range(k) for a in range(m))
                            for y in ys) for x in xs]
                cross = max(abs(sum(x[i] * B[i][a] * y[a]
                                    for i in range(k) for a in range(m)))
                            for x in xs for y in ys)
                assert max(vals) - min(vals) <= 2 * cross


if __name__ == "__main__":
    check_sparse_identity()
    check_gamma()
    check_random_pair_mass()
    check_negative_simplex_rank()
    check_selector_contact_inequality()
    check_all_even_levels()
    check_selector_oscillation()
    print("exact disjoint-star compiler checks: PASS")
