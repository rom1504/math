#!/usr/bin/env python3
"""Exact small checks for the Wave-15 balanced-partition identities."""

from itertools import combinations, product
from fractions import Fraction


def equipartitions(vertices, s):
    vertices = tuple(vertices)
    if not vertices:
        yield ()
        return
    first, rest = vertices[0], vertices[1:]
    for tail in combinations(rest, s - 1):
        block = (first,) + tail
        remain = tuple(v for v in rest if v not in tail)
        for partition in equipartitions(remain, s):
            yield (block,) + partition


def energy(A, z, partition=None):
    n = len(A)
    if partition is None:
        allowed = {(i, j) for i in range(n) for j in range(i + 1, n)}
    else:
        allowed = {
            (min(i, j), max(i, j))
            for block in partition
            for i, j in combinations(block, 2)
        }
    return 2 * sum(A[i][j] * z[i] * z[j] for i, j in allowed)


def qnorm(A, partition=None, p=None):
    n = len(A)
    vals = []
    for tail in product((-1, 1), repeat=n - 1):
        z = (1,) + tail
        if p is None:
            vals.append(abs(energy(A, z, partition)))
        else:
            vals.append(abs(energy(A, z, partition) - p * energy(A, z)))
    return max(vals)


def matvec_sq(A, z):
    return sum(sum(A[i][j] * z[j] for j in range(len(A))) ** 2
               for i in range(len(A)))


def check(n=6, s=3):
    # A deterministic mixed signing.
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = A[j][i] = 1 if (i * i + 3 * j + i * j) % 5 < 3 else -1

    parts = list(equipartitions(range(n), s))
    p = Fraction(s - 1, n - 1)
    ca = -p * Fraction(n - s, (n - 1) * (n - 2))
    cd = 2 * p * Fraction(n - s, (n - 1) * (n - 2) * (n - 3))
    N = Fraction(n * (n - 1), 2)

    # Check indicator covariance directly on representative edge pairs.
    def same(P, i, j):
        return any(i in B and j in B for B in P)

    adj = sum(Fraction(same(P, 0, 1) * same(P, 0, 2), len(parts)) for P in parts) - p * p
    dis = sum(Fraction(same(P, 0, 1) * same(P, 2, 3), len(parts)) for P in parts) - p * p
    assert adj == ca, (adj, ca)
    assert dis == cd, (dis, cd)

    QA = qnorm(A)
    for tail in product((-1, 1), repeat=n - 1):
        z = (1,) + tail
        a = energy(A, z)
        samples = [Fraction(energy(A, z, P)) for P in parts]
        mean = sum(samples, Fraction()) / len(samples)
        assert mean == p * a
        var = sum((x - mean) ** 2 for x in samples) / len(samples)
        S = Fraction(a, 2)
        L = matvec_sq(A, z)
        predicted = 4 * (
            p * (1 - p) * N
            + (ca - cd) * (L - 2 * N)
            + cd * (S * S - N)
        )
        assert var == predicted, (z, var, predicted)

    # Check the centered-mask cross inequality for every partition.
    for P in parts:
        QC = qnorm(A, P)  # this is Q(D), not Q(C); form C explicitly below
        C = [row[:] for row in A]
        for B in P:
            for i, j in combinations(B, 2):
                C[i][j] = C[j][i] = 0
        left = qnorm(C)
        centered = qnorm(A, P, p)
        assert left <= (1 - p) * QA + centered
        assert QC >= 0

    # All-positive adaptive witness: half the blocks receive each sign.
    P = parts[0]
    assert len(P) % 2 == 0
    z = [0] * n
    for block_index, B in enumerate(P):
        for i in B:
            z[i] = 1 if block_index < len(P) // 2 else -1
    J = [[0 if i == j else 1 for j in range(n)] for i in range(n)]
    centered = abs(Fraction(energy(J, z, P)) - p * energy(J, z))
    assert centered == Fraction(n * n * (s - 1), n - 1)

    print(f"PASS n={n}, s={s}, partitions={len(parts)}, Q(A)={QA}")
    print(f"ca={ca}, cd={cd}; all-positive adaptive witness={centered}")


if __name__ == "__main__":
    check()
