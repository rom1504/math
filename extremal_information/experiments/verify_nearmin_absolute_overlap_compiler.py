#!/usr/bin/env python3
"""Finite identity checks for the absolute-overlap physical compiler.

This verifier checks only exact algebra and the deterministic spherical
inequality.  It does not simulate the asymptotic concentration theorem or
claim that exact minimizers satisfy the conditional overlap hypothesis.
"""

from __future__ import annotations

from itertools import product
import math

import numpy as np


SEED = 20260817


def edges(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def cuts(n: int) -> list[np.ndarray]:
    return [np.asarray(x, dtype=np.int8) for x in product((-1, 1), repeat=n)]


def cvec(x: np.ndarray, es: list[tuple[int, int]]) -> np.ndarray:
    return np.asarray([x[i] * x[j] for i, j in es], dtype=np.int8)


def hamiltonian(a: np.ndarray, x: np.ndarray) -> int:
    return int(sum(a[i, j] * x[i] * x[j]
                   for i in range(len(x)) for j in range(i + 1, len(x))))


def cap(a: np.ndarray) -> int:
    return max(abs(hamiltonian(a, x)) for x in cuts(len(a)))


def check_sparse_expectation(rng: np.random.Generator) -> int:
    checks = 0
    for n in range(3, 9):
        es = edges(n)
        for _ in range(20):
            a = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=len(es))
            u = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=n)
            sigma = int(rng.choice(np.asarray([-1, 1], dtype=np.int8)))
            zu = sigma * cvec(u, es)
            dmask = a * zu == -1
            p = float(rng.uniform(0.01, 0.7))
            expected_b = a.astype(float).copy()
            expected_b[dmask] *= 1.0 - 2.0 * p
            for y in cuts(n)[:: max(1, 2 ** (n - 5))]:
                for tau in (-1, 1):
                    z = tau * cvec(y, es)
                    lhs = float(expected_b @ z)
                    rhs = float((1.0 - p) * (a @ z) + p * (zu @ z))
                    assert abs(lhs - rhs) < 1e-10
                    checks += 1
    return checks


def check_edge_overlap_identity(rng: np.random.Generator) -> int:
    checks = 0
    for n in range(3, 25):
        es = edges(n)
        ecount = len(es)
        deterministic: list[tuple[np.ndarray, np.ndarray]] = []
        aligned = np.ones(n, dtype=np.int8)
        one_flip = aligned.copy()
        one_flip[0] = -1
        balanced = aligned.copy()
        balanced[: n // 2] = -1
        deterministic.extend([(aligned, aligned), (aligned, one_flip),
                              (aligned, balanced)])
        pairs = deterministic + [
            (
                rng.choice(np.asarray([-1, 1], dtype=np.int8), size=n),
                rng.choice(np.asarray([-1, 1], dtype=np.int8), size=n),
            )
            for _ in range(100)
        ]
        for u, v in pairs:
            su = int(rng.choice(np.asarray([-1, 1], dtype=np.int8)))
            sv = int(rng.choice(np.asarray([-1, 1], dtype=np.int8)))
            zu = (su * cvec(u, es)).astype(np.int64)
            zv = (sv * cvec(v, es)).astype(np.int64)
            lhs = int(zu @ zv)
            rhs = su * sv * (int(u @ v) ** 2 - n) // 2
            assert lhs == rhs
            absolute = abs(lhs) / ecount
            vertex = abs(int(u @ v)) / n
            formula = abs((n * vertex * vertex - 1) / (n - 1))
            assert abs(absolute - formula) < 1e-12
            checks += 1
    return checks


def check_spherical_inequality(rng: np.random.Generator) -> int:
    checks = 0
    for dimension in (2, 3, 5, 11):
        for _ in range(100_000):
            u = rng.normal(size=dimension)
            v = rng.normal(size=dimension)
            y = rng.normal(size=dimension)
            u /= np.linalg.norm(u)
            v /= np.linalg.norm(v)
            y /= np.linalg.norm(y)
            c = abs(float(u @ v))
            theta = 1.0 - c * c
            alpha = float(10 ** rng.uniform(-2, 1))
            lam = float(10 ** rng.uniform(-2, 1))
            lhs = 0.5 * alpha * float(v @ y) ** 2 + lam * abs(float(u @ y))
            delta = theta * min(alpha, lam) / 4.0
            rhs = 0.5 * alpha + lam - delta
            assert lhs <= rhs + 1e-11
            checks += 1
    return checks


def check_free_shore_identity(rng: np.random.Generator) -> int:
    checks = 0
    for n in range(3, 7):
        xs = cuts(n)
        for h in range(1, 5):
            etas = cuts(h)
            for _ in range(10):
                upper = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=(n, n))
                a = np.triu(upper, 1)
                a = a + a.T
                u = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=n)
                upper_c = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=(h, h))
                cmat = np.triu(upper_c, 1)
                cmat = cmat + cmat.T

                trust = max(
                    tau * hamiltonian(a, y) + h * abs(int(u @ y))
                    for y in xs for tau in (-1, 1)
                )
                parent = max(
                    abs(hamiltonian(a, y)
                        + int(u @ y) * int(np.sum(eta))
                        + hamiltonian(cmat, eta))
                    for y in xs for eta in etas
                )
                qc = cap(cmat) if h > 1 else 0
                assert abs(parent - trust) <= qc
                checks += 1
    return checks


def check_antipodal_countermodel() -> int:
    checks = 0
    for n in range(10, 500):
        k = math.ceil(math.sqrt(n))
        crossing = k * (n - k)
        ecount = n * (n - 1) // 2
        mean_l1 = crossing / ecount
        cross_r = 2.0 * crossing / ecount - 1.0
        expected_r = 0.5 * (1.0 + cross_r)
        assert abs(mean_l1 - expected_r) < 1e-14
        assert mean_l1 <= 4.0 / math.sqrt(n)
        assert abs(cross_r) >= 1.0 - 4.0 / math.sqrt(n)
        checks += 1
    return checks


def main() -> None:
    rng = np.random.default_rng(SEED)
    counts = {
        "sparse_expectation": check_sparse_expectation(rng),
        "edge_overlap": check_edge_overlap_identity(rng),
        "spherical": check_spherical_inequality(rng),
        "free_shore": check_free_shore_identity(rng),
        "antipodal_countermodel": check_antipodal_countermodel(),
    }
    print(f"PASS: {sum(counts.values())} checks; {counts}")


if __name__ == "__main__":
    main()
