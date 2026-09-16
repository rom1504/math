#!/usr/bin/env python3
"""Exact Wave 50 table audit of the relaxed spectral-excess bound."""

from fractions import Fraction as F


rows = [
    (1, F(2, 27), F(1, 42), F(2, 25)),
    (2, F(1, 6), F(29, 1050), F(32, 175)),
    (3, F(2, 7), F(1, 25), F(8, 25)),
    (4, F(4, 9), F(17, 225), F(64, 125)),
    (5, F(2, 3), F(1, 5), F(4, 5)),
]


def main():
    n = 10
    actual_degree = F(1, 42)
    for ell, lam1, retention, kappa in rows:
        delta = 1 - lam1
        gap = kappa * delta
        lam2 = lam1 - gap
        bound = max(F(0), retention - lam2) / (delta + n * gap)
        assert bound <= actual_degree
        print(
            f"ell={ell}: lambda2={lam2}, P-lambda2={retention-lam2}, "
            f"degree_bound={bound}, actual={actual_degree}"
        )


if __name__ == "__main__":
    main()
