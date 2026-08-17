#!/usr/bin/env python3
"""Finite checks for drafts/flat_gram_exposed_entropy.md.

The proof in the draft is analytic.  These exact enumerations are only
regression checks for signs, normalizations, and the distinction between
signed full-cube exposure and absolute near-top exposure.
"""

from __future__ import annotations

from itertools import combinations, product
from math import sqrt
import random

import numpy as np


def cube(k: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=k)), dtype=np.int64)


def energies(coeffs: tuple[int, ...], k: int, xs: np.ndarray) -> np.ndarray:
    es = list(combinations(range(k), 2))
    out = np.zeros(len(xs), dtype=np.int64)
    for c, (i, j) in zip(coeffs, es):
        out += c * xs[:, i] * xs[:, j]
    return out


def exhaustive_sign_balance() -> None:
    # All {-2,0,2} quadratic coefficient arrays through k=5.
    for k in range(2, 6):
        xs = cube(k)
        e = k * (k - 1) // 2
        least_pos = 1.0
        least_neg = 1.0
        witness_pos = None
        witness_neg = None
        for coeffs in product((-2, 0, 2), repeat=e):
            if not any(coeffs):
                continue
            vals = energies(coeffs, k, xs)
            ppos = float(np.mean(vals > 0))
            pneg = float(np.mean(vals < 0))
            if ppos < least_pos:
                least_pos, witness_pos = ppos, coeffs
            if pneg < least_neg:
                least_neg, witness_neg = pneg, coeffs
        assert least_pos >= 1 / 324
        assert least_neg >= 1 / 324
        print(
            f"k={k}: min Pr(P>0)={least_pos:.6f} {witness_pos}; "
            f"min Pr(P<0)={least_neg:.6f} {witness_neg}"
        )


def alternating_value(mask: int, u: int, v: int, r: int) -> int:
    """Alternating form encoded by upper-triangular bit mask."""
    ans = 0
    bit = 0
    for i in range(r):
        for j in range(i + 1, r):
            ans ^= ((mask >> bit) & 1) & (
                (((u >> i) & 1) * ((v >> j) & 1))
                ^ (((u >> j) & 1) * ((v >> i) & 1))
            )
            bit += 1
    return ans


def small_gram_checks() -> None:
    rng = random.Random(20260817)
    for r in (2, 3, 4):
        k = 1 << r
        xs = cube(k)
        base = np.zeros((k, k), dtype=np.int64)
        for i in range(k):
            for j in range(i + 1, k):
                base[i, j] = base[j, i] = rng.choice((-1, 1))

        h = r * (r - 1) // 2
        min_nonnegative = 1.0
        min_abs_halfmax = 1.0
        for mask in range(1, 1 << h):
            d = np.zeros((k, k), dtype=np.int64)
            for i in range(k):
                for j in range(i + 1, k):
                    chi = -1 if alternating_value(mask, i, j, r) else 1
                    # A_mask - A_zero.
                    d[i, j] = d[j, i] = base[i, j] * (chi - 1)
            vals = np.einsum("bi,ij,bj->b", xs, d, xs, optimize=True) // 2
            q = int(np.max(np.abs(vals)))
            if q == 0:
                continue
            # Orient so the absolute extremum is the minimum.
            s = 1 if -int(np.min(vals)) >= int(np.max(vals)) else -1
            oriented = s * vals
            assert int(np.min(oriented)) == -q
            nonnegative = float(np.mean(oriented >= 0))
            abs_half = float(np.mean(np.abs(vals) >= q / 2))
            min_nonnegative = min(min_nonnegative, nonnegative)
            min_abs_halfmax = min(min_abs_halfmax, abs_half)
            assert nonnegative >= 1 / 324

            # Check the exact spectral Hamming Lipschitz inequality on all
            # pairs against one maximizer for these tiny examples.
            zidx = int(np.argmax(np.abs(vals)))
            z = xs[zidx]
            op = float(np.linalg.norm(d.astype(float), ord=2))
            for row, value in zip(xs, vals):
                dist = int(np.sum(row != z))
                assert abs(int(value) - int(vals[zidx])) <= (
                    2 * op * sqrt(dist * k) + 1e-8
                )

        print(
            f"r={r}, k={k}: min oriented nonnegative mass="
            f"{min_nonnegative:.6f}; min absolute half-max mass="
            f"{min_abs_halfmax:.6f}"
        )


def block_absolute_counterexample() -> None:
    for s in (2, 4):
        k = s * s
        xs = cube(k)
        d = np.zeros((k, k), dtype=np.int64)
        for b in range(s):
            block = range(b * s, (b + 1) * s)
            for i, j in combinations(block, 2):
                d[i, j] = d[j, i] = 2
        vals = np.einsum("bi,ij,bj->b", xs, d, xs, optimize=True) // 2
        assert int(np.max(vals)) == s**3 - s**2
        assert int(np.min(vals)) == -(s**2)
        assert abs(float(np.linalg.norm(d.astype(float), ord=2)) - 2 * (s - 1)) < 1e-8
        print(
            f"block-flat s={s}, k={k}: q={int(np.max(np.abs(vals)))}, "
            f"half-max mass={float(np.mean(np.abs(vals) >= np.max(vals)/2)):.6f}"
        )


def main() -> None:
    exhaustive_sign_balance()
    small_gram_checks()
    block_absolute_counterexample()
    print("flat Gram exposed-entropy checks passed")


if __name__ == "__main__":
    main()
