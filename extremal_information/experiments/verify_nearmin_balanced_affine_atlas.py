#!/usr/bin/env python3
"""Finite verifier for AA.1--AA.2 and the AA.13 selector ingredient.

The run is exhaustive over every signing through n=5 and uses a fixed-seed
sample at n=6,7,8.  For every signing it checks every oriented spin atom and
every q, including the actual minimum-cost partition cell used in the proof.
It also enumerates the star-frame endpoints when the selected support permits
a nontrivial odd port frame.  The shell-law sampling statement AA.3 and the
response maximization in AA.13 are proved analytically, not instantiated by
this finite checker.
"""

from __future__ import annotations

import itertools
import math
import random


def energy(a, x):
    n = len(x)
    return sum(a[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))


def all_spins(n):
    return itertools.product((-1, 1), repeat=n)


def signing_from_bits(n, bits):
    a = [[0] * n for _ in range(n)]
    t = 0
    for i in range(n):
        for j in range(i + 1, n):
            v = 1 if (bits >> t) & 1 else -1
            a[i][j] = a[j][i] = v
            t += 1
    return a


def oriented_matrix(a, sigma, x):
    n = len(x)
    return [[0 if i == j else sigma * x[i] * a[i][j] * x[j] for j in range(n)] for i in range(n)]


def one_sided_negative_cap(d, vertices):
    if not vertices:
        return 0
    minimum = min(
        sum(
            d[vertices[i]][vertices[j]] * s[i] * s[j]
            for i in range(len(vertices))
            for j in range(i + 1, len(vertices))
        )
        for s in all_spins(len(vertices))
    )
    return -minimum


def balanced_partition(n, q):
    return [list(range(b, n, q)) for b in range(q)]


def check_signing(a):
    n = len(a)
    energies = [energy(a, x) for x in all_spins(n)]
    qcap = max(abs(v) for v in energies)
    atoms = 0
    frames = 0
    for sigma in (-1, 1):
        for x0 in all_spins(n):
            x = tuple(x0)
            p0 = sigma * energy(a, x)
            deficit = qcap - p0
            dmat = oriented_matrix(a, sigma, x)
            ell = [sum(row) for row in dmat]
            lminus = sum(max(-v, 0) for v in ell)
            assert lminus * lminus <= 4 * qcap * deficit
            for q in range(2, n + 1):
                cells = balanced_partition(n, q)
                costs = []
                for cell in cells:
                    cost = 2 * sum(max(ell[i], 0) for i in cell)
                    cost += 4 * one_sided_negative_cap(dmat, cell)
                    costs.append(cost)
                b = min(range(q), key=costs.__getitem__)
                k = n // q
                support = cells[b][:k]
                actual_cost = 2 * sum(max(ell[i], 0) for i in support)
                actual_cost += 4 * one_sided_negative_cap(dmat, support)
                analytic = (8 * qcap + 4 * math.sqrt(qcap * deficit)) / q
                assert actual_cost <= analytic + 1e-10
                for mask in range(1 << k):
                    y = list(x)
                    for t, i in enumerate(support):
                        if (mask >> t) & 1:
                            y[i] *= -1
                    assert sigma * energy(a, y) >= p0 - actual_cost

                # Check every endpoint of the largest even star subframe.
                kp = k - (k % 2)
                if kp:
                    chosen = support[:kp]
                    columns = [x]
                    for i in chosen:
                        y = list(x)
                        y[i] *= -1
                        columns.append(tuple(y))
                    for eps in all_spins(kp + 1):
                        field = [sum(eps[j] * columns[j][i] for j in range(kp + 1)) for i in range(n)]
                        assert all(v != 0 for v in field)
                        selector = tuple(1 if v > 0 else -1 for v in field)
                        assert sigma * energy(a, selector) >= p0 - actual_cost
                    frames += 1
            atoms += 1
    return atoms, frames


def main():
    rng = random.Random(20260817)
    signings = 0
    atoms = 0
    frames = 0
    for n in range(2, 6):
        for bits in range(1 << (n * (n - 1) // 2)):
            aa, ff = check_signing(signing_from_bits(n, bits))
            signings += 1
            atoms += aa
            frames += ff
    for n in (6, 7, 8):
        edge_count = n * (n - 1) // 2
        for _ in range(16):
            aa, ff = check_signing(signing_from_bits(n, rng.getrandbits(edge_count)))
            signings += 1
            atoms += aa
            frames += ff
    print({"status": "PASS", "signings": signings, "oriented_atoms": atoms, "star_frames": frames})


if __name__ == "__main__":
    main()
