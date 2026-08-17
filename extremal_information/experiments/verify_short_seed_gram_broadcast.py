#!/usr/bin/env python3
"""Finite arithmetic checks for short-seed Gram broadcast.

This verifies the displayed constants, exhaustively checks the small-bias
radical estimate in tiny vector spaces, and constructs a toy polynomial-
trace t-wise-independent sample space.  It does not search the theorem-scale
good seed (whose guaranteed order is already at least 1024).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, exp, lgamma, log, log2


def dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def alternating_value(rows: tuple[int, ...], p: int, q: int) -> int:
    """B(p,q) for a symmetric zero-diagonal binary matrix given by rows."""
    out = 0
    r = len(rows)
    for i in range(r):
        if (p >> i) & 1:
            out ^= dot(rows[i], q)
    return out


def alternating_rows(r: int, mask: int) -> tuple[int, ...]:
    rows = [0] * r
    bit = 0
    for i in range(r):
        for j in range(i + 1, r):
            if (mask >> bit) & 1:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
            bit += 1
    return tuple(rows)


def bias_of_multiset(points: tuple[int, ...], r: int) -> Fraction:
    worst = Fraction(0)
    for ell in range(1, 1 << r):
        total = sum(1 if dot(ell, p) == 0 else -1 for p in points)
        worst = max(worst, abs(Fraction(total, len(points))))
    return worst


def check_radical_bound(points: tuple[int, ...], r: int) -> None:
    delta = bias_of_multiset(points, r)
    if delta > Fraction(1, 8):
        return
    h = comb(r, 2)
    for mask in range(1, 1 << h):
        rows = alternating_rows(r, mask)
        ones = sum(
            alternating_value(rows, p, q)
            for p in points
            for q in points
        )
        prob = Fraction(ones, len(points) ** 2)
        assert prob >= Fraction(1, 4), (r, mask, delta, prob)


# GF(2^3) with x^3+x+1; enough for a toy four-position, four-wise panel.
MOD = 0b1011


def gf_mul(a: int, b: int) -> int:
    out = 0
    x = a
    y = b
    while y:
        if y & 1:
            out ^= x
        y >>= 1
        x <<= 1
        if x & 0b1000:
            x ^= MOD
    return out & 0b111


def gf_pow(a: int, n: int) -> int:
    out = 1
    while n:
        if n & 1:
            out = gf_mul(out, a)
        a = gf_mul(a, a)
        n >>= 1
    return out


def gf_trace(a: int) -> int:
    return a ^ gf_pow(a, 2) ^ gf_pow(a, 4)


def poly_eval(coeffs: tuple[int, ...], x: int) -> int:
    out = 0
    for c in reversed(coeffs):
        out = gf_mul(out, x) ^ c
    return out


def check_twise_toy() -> None:
    positions = (0, 1, 2, 3)
    t = 4
    counts: dict[tuple[int, ...], int] = {}
    for coeffs in product(range(8), repeat=t):
        word = tuple(gf_trace(poly_eval(coeffs, x)) & 1 for x in positions)
        counts[word] = counts.get(word, 0) + 1
    assert len(counts) == 1 << t
    assert len(set(counts.values())) == 1


def check_constants() -> None:
    for r in range(2, 100):
        s = 256 * r
        k = s * r
        h = comb(r, 2)
        assert (r + 1) * log(2) - s / 128 < 0
        assert h >= k / 1024
        log_failure = k * log(9) + h * log(2) + 3 * k * log(Fraction(3, 8))
        assert log_failure < -0.743 * k
        assert 6 * k <= comb(k, 2)
        # Rao/orthogonal-character support bound is already Omega(k log k).
        E = comb(k, 2)
        log_choose = (lgamma(E + 1) - lgamma(3 * k + 1) - lgamma(E - 3 * k + 1)) / log(2)
        assert log_choose >= k * log2(k)
    assert exp(log(9) + log(2) / 512 + 3 * log(Fraction(3, 8))) < exp(-0.743)


def check_small_bias_examples() -> None:
    # Full spaces are zero-biased and permit exhaustive checks of SG.8--SG.10.
    for r in range(2, 5):
        points = tuple(range(1 << r))
        assert bias_of_multiset(points, r) == 0
        check_radical_bound(points, r)


def main() -> None:
    check_constants()
    check_small_bias_examples()
    check_twise_toy()
    print("short-seed Gram-broadcast checks passed")


if __name__ == "__main__":
    main()
