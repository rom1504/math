#!/usr/bin/env python3
"""Exact checks for the Wave 37 Gram-projection row transfer."""

from __future__ import annotations

import itertools
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def audit(a: np.ndarray, name: str, m: int) -> None:
    n = len(a)
    op = np.linalg.norm(a, 2)
    selectors = list(itertools.combinations(range(n), m))
    for z in spins(n):
        az = a @ z
        row = int(az @ az)
        kvals = []
        for s in selectors:
            ss = list(s)
            vec = a[:, ss] @ z[ss]
            kvals.append(int(vec @ vec))

            # Uniform outside completion identity.
            vals = []
            outside = [i for i in range(n) if i not in s]
            for tail in itertools.product((-1, 1), repeat=len(outside)):
                zz = z.copy()
                zz[outside] = tail
                vals.append(int((a @ zz) @ (a @ zz)))
            assert Fraction(sum(vals), len(vals)) == kvals[-1] + (n - m) * (n - 1)

            # Pointwise Hamming/Minkowski inequality, checked after squaring
            # only through the underlying vector triangle inequality.
            for ytail in itertools.product((-1, 1), repeat=m):
                y = np.asarray(ytail, dtype=np.int64)
                lhs = np.linalg.norm(a[:, ss] @ z[ss])
                rhs0 = np.linalg.norm(a[:, ss] @ y)
                e = int(np.sum(z[ss] != y))
                assert lhs <= rhs0 + 2 * op * np.sqrt(e) + 1e-9

        p = Fraction(m, n)
        p2 = Fraction(m * (m - 1), n * (n - 1))
        expected = Fraction(sum(kvals), len(kvals))
        assert expected == p2 * row + (p - p2) * n * (n - 1)

        # Anchored identity and scalar projection identities for every v.
        if row:
            u = az.astype(float) / np.sqrt(row)
            cols = [z[i] * a[:, i] for i in range(n)]
            coeff = np.asarray([u @ col for col in cols])
            assert abs(coeff.sum() - np.sqrt(row)) < 1e-9
            assert coeff @ coeff <= op ** 2 + 1e-8
        for v in range(n):
            anchored = [k for k, s in zip(kvals, selectors) if v in s]
            alpha1 = Fraction(m - 1, n - 1)
            alpha2 = Fraction((m - 1) * (m - 2), (n - 1) * (n - 2))
            tv = int(z[v] * (a @ az)[v] - (n - 1))
            rhs = m * (n - 1) + alpha2 * (row - n * (n - 1)) + 2 * (alpha1 - alpha2) * tv
            assert Fraction(sum(anchored), len(anchored)) == rhs

            if row:
                scalar_mean = np.mean([
                    sum(coeff[i] for i in s) for s in selectors if v in s
                ])
                rhs_scalar = float(alpha1) * np.sqrt(row) + (1 - float(alpha1)) * coeff[v]
                assert abs(scalar_mean - rhs_scalar) < 1e-9

    print(name, "n", n, "m", m, "selectors", len(selectors), "PASS")


def main() -> None:
    audit(A6, "A6", 3)
    audit(A8, "A8", 5)
    audit(A9, "A9", 6)
    print("gram_projection_row_r37_check: PASS")


if __name__ == "__main__":
    main()
