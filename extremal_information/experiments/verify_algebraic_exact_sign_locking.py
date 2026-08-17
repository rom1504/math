#!/usr/bin/env python3
"""Exact wind tunnel for the algebraic exact-sign locking audit.

The checks are deliberately finite and integer-valued:

* enumerate every complete sign Hamiltonian on two paired blocks for k=2,3;
* classify those for which every coordinatewise duplicate is a ground state;
* verify the closed antibalanced energy formula and its mismatch ties;
* exercise the Boolean-eigenspace/star-row local-margin obstruction.
"""

from __future__ import annotations

import itertools
import math
import random

import numpy as np


def energy(A: np.ndarray, z: np.ndarray) -> int:
    return int(z @ A @ z // 2)


def matrices(n: int):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in itertools.product((-1, 1), repeat=len(edges)):
        A = np.zeros((n, n), dtype=int)
        for (i, j), a in zip(edges, bits):
            A[i, j] = A[j, i] = a
        yield A


def is_duplicate_ground_lock(A: np.ndarray, k: int) -> bool:
    cube = [np.array(z, dtype=int) for z in itertools.product((-1, 1), repeat=2 * k)]
    top = max(energy(A, z) for z in cube)
    return all(energy(A, np.r_[u, u]) == top
               for u in itertools.product((-1, 1), repeat=k))


def classified_sign(A: np.ndarray, k: int) -> np.ndarray | None:
    # Coordinates are x_0,...,x_(k-1),y_0,...,y_(k-1).
    if any(A[i, k + i] != 1 for i in range(k)):
        return None
    a = np.zeros((k, k), dtype=int)
    for i in range(k):
        for j in range(i + 1, k):
            v = A[i, j]
            block = np.array([[A[i, j], A[i, k + j]],
                              [A[k + i, j], A[k + i, k + j]]])
            if not np.array_equal(block, v * np.array([[1, -1], [-1, 1]])):
                return None
            a[i, j] = a[j, i] = v
    # Recover a_ij=-s_i s_j.  For k=1 the empty choice is harmless.
    if k == 1:
        return np.ones(1, dtype=int)
    s = np.ones(k, dtype=int)
    s[0] = -1
    for i in range(1, k):
        s[i] = a[0, i]
    if any(a[i, j] != -s[i] * s[j]
           for i in range(k) for j in range(i + 1, k)):
        return None
    return s


def check_exhaustive_classification() -> int:
    checks = 0
    for k, expected in ((2, 2), (3, 4)):
        found = 0
        for A in matrices(2 * k):
            ground = is_duplicate_ground_lock(A, k)
            s = classified_sign(A, k)
            assert ground == (s is not None)
            if ground:
                found += 1
            checks += 1
        assert found == expected
        checks += 1
    return checks


def lock_matrix(s: np.ndarray) -> np.ndarray:
    k = len(s)
    A = np.zeros((2 * k, 2 * k), dtype=int)
    for i in range(k):
        A[i, k + i] = A[k + i, i] = 1
    for i in range(k):
        for j in range(i + 1, k):
            a = -s[i] * s[j]
            block = a * np.array([[1, -1], [-1, 1]])
            inds_i, inds_j = (i, k + i), (j, k + j)
            for p in range(2):
                for q in range(2):
                    A[inds_i[p], inds_j[q]] = A[inds_j[q], inds_i[p]] = block[p, q]
    return A


def check_formula_and_ties() -> int:
    rng = random.Random(20260817)
    checks = 0
    for k in range(2, 9):
        s = np.array([rng.choice((-1, 1)) for _ in range(k)], dtype=int)
        A = lock_matrix(s)
        for zt in itertools.product((-1, 1), repeat=2 * k):
            z = np.array(zt, dtype=int)
            x, y = z[:k], z[k:]
            d = (x - y) // 2
            target = k - 2 * int(np.dot(s, d)) ** 2
            assert energy(A, z) == target
            checks += 1
        # Every pair has a nonduplicate two-mismatch ground state.
        x = np.ones(k, dtype=int)
        y = x.copy()
        y[0] = -x[0]
        # Choose x_1 so that s_0 d_0+s_1 d_1=0 after flipping y_1.
        x[1] = -s[0] * s[1] * x[0]
        y[1] = -x[1]
        assert not np.array_equal(x, y)
        assert energy(A, np.r_[x, y]) == k
        checks += 1
    return checks


def sylvester(n: int) -> np.ndarray:
    H = np.array([[1]], dtype=int)
    while len(H) < n:
        H = np.block([[H, H], [H, -H]])
    return H


def check_spectral_star_obstruction() -> int:
    rng = random.Random(271828)
    checks = 0
    for n in (4, 16):
        H = sylvester(n)
        lam = int(math.isqrt(n))
        code = [np.array(u, dtype=int)
                for u in itertools.product((-1, 1), repeat=n)
                if np.array_equal(H @ np.array(u), lam * np.array(u))]
        assert code
        u = code[0]
        i = 0
        A = np.zeros((n, n), dtype=int)
        for j in range(1, n):
            A[i, j] = A[j, i] = -u[i] * u[j]
        for a in range(1, n):
            for b in range(a + 1, n):
                A[a, b] = A[b, a] = rng.choice((-1, 1))
        assert u[i] * int((A @ u)[i]) == -(n - 1)
        # Repeated H-lock: the intended codeword fails the one-flip test
        # whenever the number of blocks is below the exact threshold.
        for blocks in range(1, max(2, (n - 1) // lam + 2)):
            margin = (blocks - 1) * lam + u[i] * int((A @ u)[i])
            assert (margin < 0) == ((blocks - 1) * lam < n - 1)
            checks += 1
        # The prescribed star costs only sqrt(n-1) in operator norm; the
        # random exact completion remains at the natural sqrt(n) scale here.
        assert np.linalg.norm(A, 2) < 6 * math.sqrt(n)
        checks += 1
    return checks


def main() -> None:
    count = 0
    count += check_exhaustive_classification()
    count += check_formula_and_ties()
    count += check_spectral_star_obstruction()
    print(f"algebraic exact-sign locking checks passed: {count}")


if __name__ == "__main__":
    main()
