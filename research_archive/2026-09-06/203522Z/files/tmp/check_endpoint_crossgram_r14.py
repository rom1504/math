#!/usr/bin/env python3
"""Exact checks for the Wave 14 balanced endpoint cross-Gram memo.

All decisions are integer/rational.  NumPy is used only as a convenient
integer-array container; no floating-point assertion is made.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product

import numpy as np


def mat_from_code(n: int, code: int) -> np.ndarray:
    """First-row-positive representative of a switching class."""
    A = np.zeros((n, n), dtype=np.int64)
    A[0, 1:] = 1
    A[1:, 0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            A[i, j] = A[j, i] = 1 if (code >> k) & 1 else -1
            k += 1
    return A


def projective_spins(n: int) -> np.ndarray:
    return np.asarray(
        [[1, *bits] for bits in product((-1, 1), repeat=n - 1)],
        dtype=np.int64,
    )


def endpoint_data(A: np.ndarray):
    X = projective_spins(len(A))
    energy = np.einsum("bi,ij,bj->b", X, A, X, optimize=True)
    P = int(energy.max())
    N = int(-energy.min())
    Gp = X[energy == P].T.copy()
    Gn = X[energy == -N].T.copy()
    return P, N, Gp, Gn


def rank_exact(A: np.ndarray) -> int:
    rows = [[Fraction(int(x)) for x in row] for row in A.tolist()]
    if not rows:
        return 0
    m, n = len(rows), len(rows[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        q = rows[r][c]
        rows[r] = [x / q for x in rows[r]]
        for i in range(m):
            if i != r and rows[i][c]:
                q = rows[i][c]
                rows[i] = [x - q * y for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == m:
            break
    return r


def det_bareiss(A: np.ndarray) -> int:
    M = [[int(x) for x in row] for row in A.tolist()]
    n = len(M)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if M[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            M[k], M[pivot] = M[pivot], M[k]
            sign *= -1
        p = M[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = M[i][j] * p - M[i][k] * M[k][j]
                assert numerator % prev == 0
                M[i][j] = numerator // prev
        prev = p
    return sign * M[-1][-1]


def charpoly_coeffs(A: np.ndarray) -> tuple[int, ...]:
    """Faddeev--LeVerrier coefficients of det(tI-A)."""
    n = len(A)
    Ao = A.astype(object)
    I = np.eye(n, dtype=object)
    B = I.copy()
    ans = [1]
    for k in range(1, n + 1):
        tr = int(np.trace(Ao @ B))
        assert tr % k == 0
        c = -tr // k
        ans.append(c)
        B = Ao @ B + c * I
    return tuple(ans)


def assert_zero(A: np.ndarray) -> None:
    assert np.count_nonzero(A) == 0, A


def endpoint_residual_identity(A: np.ndarray) -> int:
    """Check the general threshold formula on all endpoint pairs."""
    P, N, Gp, Gn = endpoint_data(A)
    I = abs(P - N)
    pairs = 0
    for p in Gp.T:
        Ap = A * np.outer(p, p)
        for n0 in Gn.T:
            n = p * n0  # negative endpoint in the p-gauge
            S = np.flatnonzero(n == 1)
            T = np.flatnonzero(n == -1)
            hS = int(Ap[np.ix_(S, S)].sum())
            hT = int(Ap[np.ix_(T, T)].sum())
            H = Fraction(P - N, 2)
            b = int(p @ A @ n0)
            assert hS + hT == H
            assert hS - hT == b
            R = max(Fraction(abs(hS)) - abs(H), 0) + max(
                Fraction(abs(hT)) - abs(H), 0
            )
            closed = Fraction(
                max(Fraction(abs(b)) - abs(H), 0)
                + max(Fraction(abs(b)) - 3 * abs(H), 0),
                2,
            )
            assert R == closed
            if P == N:
                assert R == abs(b)
                assert b % 4 == 0
            pairs += 1
    return pairs


def enumerate_small_orders():
    expected = {
        2: (2, 1, 1),
        3: (6, 2, 0),
        4: (8, 6, 6),
        5: (8, 12, 12),
        6: (10, 12, 12),
        7: (18, 3240, 0),
    }
    balanced_types = {}
    total_endpoint_pairs = 0
    for n in range(2, 8):
        best = n * n
        minimizers = []
        count = 1 << ((n - 1) * (n - 2) // 2)
        for code in range(count):
            A = mat_from_code(n, code)
            P, N, _, _ = endpoint_data(A)
            q = max(P, N)
            if q < best:
                best, minimizers = q, [(code, A, P, N)]
            elif q == best:
                minimizers.append((code, A, P, N))
        balanced = [row for row in minimizers if row[2] == row[3]]
        assert (best, len(minimizers), len(balanced)) == expected[n]
        types = Counter()
        for _, A, P, N in balanced:
            _, _, Gp, Gn = endpoint_data(A)
            C = Gp.T @ A @ Gn
            typ = (
                Gp.shape[1],
                Gn.shape[1],
                rank_exact(Gp),
                rank_exact(Gn),
                rank_exact(A),
                rank_exact(C),
                tuple(sorted(Counter(abs(int(x)) for x in C.flat).items())),
            )
            types[typ] += 1
            total_endpoint_pairs += endpoint_residual_identity(A)
        balanced_types[n] = types

    expected_types = {
        2: Counter({(1, 1, 1, 1, 2, 0, ((0, 1),)): 1}),
        3: Counter(),
        4: Counter({(1, 1, 1, 1, 4, 1, ((4, 1),)): 6}),
        5: Counter({(5, 5, 5, 5, 4, 4, ((0, 5), (4, 20))): 12}),
        6: Counter({(6, 6, 6, 6, 6, 6, ((0, 6), (4, 30))): 12}),
        7: Counter(),
    }
    assert balanced_types == expected_types
    return total_endpoint_pairs


EXAMPLES = {
    2: mat_from_code(2, 0),
    4: mat_from_code(4, 1),
    5: mat_from_code(5, 13),
    6: mat_from_code(6, 220),
    8: np.asarray(
        [
            [0, 1, 1, 1, 1, 1, 1, 1],
            [1, 0, 1, -1, 1, 1, -1, -1],
            [1, 1, 0, 1, -1, 1, -1, -1],
            [1, -1, 1, 0, -1, -1, -1, 1],
            [1, 1, -1, -1, 0, -1, 1, -1],
            [1, 1, 1, -1, -1, 0, 1, 1],
            [1, -1, -1, -1, 1, 1, 0, 1],
            [1, -1, -1, 1, -1, 1, 1, 0],
        ],
        dtype=np.int64,
    ),
}


def check_named_examples():
    expected = {
        # n: (P,N,g+,g-,d+,d-,rank A,rank C,det A)
        2: (2, 2, 1, 1, 1, 1, 2, 0, -1),
        4: (8, 8, 1, 1, 1, 1, 4, 1, 5),
        5: (8, 8, 5, 5, 5, 5, 4, 4, 0),
        6: (10, 10, 6, 6, 6, 6, 6, 6, -125),
        8: (20, 20, 4, 4, 4, 4, 8, 0, 729),
    }
    for n, A in EXAMPLES.items():
        P, N, Gp, Gn = endpoint_data(A)
        C = Gp.T @ A @ Gn
        got = (
            P,
            N,
            Gp.shape[1],
            Gn.shape[1],
            rank_exact(Gp),
            rank_exact(Gn),
            rank_exact(A),
            rank_exact(C),
            det_bareiss(A),
        )
        assert got == expected[n], (n, got)
        endpoint_residual_identity(A)

    I2 = np.eye(2, dtype=np.int64)
    I4 = np.eye(4, dtype=np.int64)
    I5 = np.eye(5, dtype=np.int64)
    I6 = np.eye(6, dtype=np.int64)
    I8 = np.eye(8, dtype=np.int64)
    A2, A4, A5, A6, A8 = (EXAMPLES[n] for n in (2, 4, 5, 6, 8))

    # These annihilating identities, ranks and traces certify all stated
    # singular-value multisets without numerical eigensolvers.
    assert_zero(A2 @ A2 - I2)
    assert_zero((A4 @ A4 - I4) @ (A4 @ A4 - 5 * I4))
    assert np.trace(A4 @ A4) == 12
    assert_zero(A5 @ (A5 @ A5 - 5 * I5))
    assert rank_exact(A5) == 4
    assert_zero(A6 @ A6 - 5 * I6)
    assert_zero((A8 @ A8 - I8) @ (A8 @ A8 - 9 * I8))
    assert np.trace(A8 @ A8) == 56

    expected_A_charpoly = {
        2: (1, 0, -1),
        4: (1, 0, -6, 0, 5),
        5: (1, 0, -10, 0, 25, 0),
        6: (1, 0, -15, 0, 75, 0, -125),
        8: (1, 0, -28, 0, 270, 0, -972, 0, 729),
    }
    assert {n: charpoly_coeffs(A) for n, A in EXAMPLES.items()} == expected_A_charpoly

    # Cross-Gram singular values and endpoint-frame spectra.
    for n in (2, 4, 5, 6, 8):
        A = EXAMPLES[n]
        _, _, Gp, Gn = endpoint_data(A)
        C = Gp.T @ A @ Gn
        if n in (2, 8):
            assert_zero(C)
        elif n == 4:
            assert C.shape == (1, 1) and abs(int(C[0, 0])) == 4
        elif n == 5:
            X = C @ C.T
            assert_zero(X @ (X - 80 * np.eye(5, dtype=np.int64)))
            assert rank_exact(C) == 4
        elif n == 6:
            assert_zero(C @ C.T - 80 * np.eye(6, dtype=np.int64))

        for G in (Gp, Gn):
            X = G @ G.T
            if n == 2:
                assert_zero(X @ (X - 2 * I2))
            elif n == 4:
                assert_zero(X @ (X - 4 * I4))
            elif n == 5:
                assert_zero(X @ (X - I5) @ ((X - 6 * I5) @ (X - 6 * I5) - 20 * I5))
            elif n == 6:
                assert_zero((X - 6 * I6) @ (X - 6 * I6) - 20 * I6)
            elif n == 8:
                assert_zero(X @ (X - 4 * I8) @ (X - 12 * I8))
                assert rank_exact(X) == 4
                assert np.trace(X) == 32
                assert np.trace(X @ X) == 320

    # A8 is the sharp complementary-sector wall: the endpoint spans are
    # Euclidean orthogonal as well as A-orthogonal and are A-invariant.
    _, _, Gp8, Gn8 = endpoint_data(A8)
    assert_zero(Gp8.T @ Gn8)
    assert_zero(Gp8.T @ A8 @ Gn8)
    assert rank_exact(np.concatenate((Gp8, Gn8), axis=1)) == 8
    # Restriction spectra are (-1,3,3,3) and (-3,-3,-3,1).
    for G, cp in (
        (Gp8, (1, -8, 18, 0, -27)),
        (Gn8, (1, 8, 18, 0, -27)),
    ):
        rows = None
        for ix in product(range(8), repeat=4):
            if len(set(ix)) == 4:
                B = G[list(ix), :]
                if det_bareiss(B):
                    rows = list(ix)
                    break
        assert rows is not None
        # Fractional inverse solve for the coordinate matrix M in AG=GM.
        B = [[Fraction(int(x)) for x in row] for row in G[rows, :].tolist()]
        Y = [[Fraction(int(x)) for x in row] for row in (A8 @ G)[rows, :].tolist()]
        aug = [b + y for b, y in zip(B, Y)]
        for c in range(4):
            p = next(i for i in range(c, 4) if aug[i][c])
            aug[c], aug[p] = aug[p], aug[c]
            q = aug[c][c]
            aug[c] = [x / q for x in aug[c]]
            for i in range(4):
                if i != c:
                    q = aug[i][c]
                    aug[i] = [x - q * y for x, y in zip(aug[i], aug[c])]
        Mq = [row[4:] for row in aug]
        assert all(x.denominator == 1 for row in Mq for x in row)
        M = np.asarray([[int(x) for x in row] for row in Mq], dtype=np.int64)
        assert_zero(A8 @ G - G @ M)
        assert charpoly_coeffs(M) == cp


def check_general_formula_sample():
    # Exhaust every switching class and every endpoint pair through order six,
    # including unbalanced matrices, to audit the threshold constants.
    pairs = 0
    for n in range(2, 7):
        count = 1 << ((n - 1) * (n - 2) // 2)
        for code in range(count):
            pairs += endpoint_residual_identity(mat_from_code(n, code))
    return pairs


if __name__ == "__main__":
    general_pairs = check_general_formula_sample()
    minimizer_pairs = enumerate_small_orders()
    check_named_examples()
    print("PASS balanced endpoint cross-Gram checks")
    print("general endpoint pairs through n=6:", general_pairs)
    print("balanced-minimizer endpoint pairs through n=7:", minimizer_pairs)
    print("A8: det=729, singular values 3^6,1^2, cross Gram zero")
