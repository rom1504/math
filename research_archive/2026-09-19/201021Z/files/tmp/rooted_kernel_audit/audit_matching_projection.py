#!/usr/bin/env python3
"""Audit the partial-matching moving-projection kernel.

This is a scratch research checker.  It verifies the exact matching-layer
transition path and evaluates the rooted cut-code inequality on saved small
signings.  It writes no project files and uses no temporary directory.
"""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path("/home/math/quadra")


def matching_count(n: int, ell: int) -> int:
    return math.comb(n, 2 * ell) * math.prod(range(1, 2 * ell, 2))


def matching_path(n: int, L: int) -> tuple[float, np.ndarray, np.ndarray]:
    """Return Perron value, Perron vector, and overlap layer weights alpha."""
    N = math.comb(n, 2)
    J = np.zeros((L + 1, L + 1), dtype=float)
    for ell in range(L):
        dplus = math.comb(n - 2 * ell, 2)
        c = math.sqrt((ell + 1) * dplus) / N
        J[ell, ell + 1] = J[ell + 1, ell] = c
    vals, vecs = np.linalg.eigh(J)
    v = vecs[:, -1]
    if v.sum() < 0:
        v = -v
    assert np.all(v > -1e-12)
    dims = np.asarray([matching_count(n, ell) for ell in range(L + 1)], float)
    w = np.sqrt(dims) * v
    alpha = w / w.sum()
    return float(vals[-1]), v, alpha


def standard_matching_path(
    n: int, L: int
) -> tuple[float, np.ndarray, np.ndarray]:
    """The common standard S^(n-1,1) module on matching layers 1,...,L."""
    N = math.comb(n, 2)
    if not (1 <= L <= (n - 1) // 2):
        raise ValueError("standard path requires 1 <= L <= floor((n-1)/2)")
    J = np.zeros((L, L), dtype=float)
    for ell in range(1, L):
        dplus = math.comb(n - 2 * ell, 2)
        up_sq = ell * (n - 2 * ell - 1) * (n - 2 * ell - 2) / 2
        c = up_sq / (N * math.sqrt((ell + 1) * dplus))
        J[ell - 1, ell] = J[ell, ell - 1] = c
    vals, vecs = np.linalg.eigh(J)
    v = vecs[:, -1]
    if v.sum() < 0:
        v = -v
    dims = np.asarray([matching_count(n, ell) for ell in range(1, L + 1)], float)
    w = np.sqrt(dims) * v
    alpha = w / w.sum()
    return float(vals[-1]), v, alpha


def elementary_sign_sums(n: int, negatives: int) -> list[int]:
    """Coefficients of (1+t)^(n-r)(1-t)^r."""
    coeff = [1]
    for sign, count in ((1, n - negatives), (-1, negatives)):
        for _ in range(count):
            coeff.append(0)
            for j in range(len(coeff) - 1, 0, -1):
                coeff[j] += sign * coeff[j - 1]
    return coeff


def code_budget(n: int, alpha: np.ndarray, lam: float) -> float:
    """Compute E_c[(tau(c)-lambda) kappa(c)] exactly up to FP summation."""
    total = 0.0
    for r in range(n + 1):
        e = elementary_sign_sums(n, r)
        even = 0.0
        odd = 0.0
        for ell, aell in enumerate(alpha):
            term = aell * e[2 * ell] / math.comb(n, 2 * ell)
            if ell & 1:
                odd += term
            else:
                even += term
        # Average over the global sign sigma.
        kappa = even * even + odd * odd
        tau_kappa = 2.0 * even * odd * (
            ((n - 2 * r) ** 2 - n) / (n * (n - 1))
        )
        total += math.comb(n, r) * (tau_kappa - lam * kappa)
    return total / (2**n)


def standard_layer_operator(n: int, ell: int, x: np.ndarray) -> np.ndarray:
    """Return phi_l^* D_cut phi_l on the standard vertex module."""
    size = 2 * ell
    norm_b_sq = size * (n - size) / n
    scale = (n - 1) / (math.comb(n, size) * norm_b_sq)
    out = np.zeros((n, n), dtype=float)
    one = np.ones(n, dtype=float)
    for subset in combinations(range(n), size):
        b = -size / n * one
        b = b.copy()
        b[list(subset)] += 1.0
        xprod = int(np.prod(x[list(subset)]))
        out += xprod * np.outer(b, b)
    return scale * out


def standard_code_budget(n: int, alpha: np.ndarray, lam: float) -> float:
    """Compute the standard-module complete-Gram code budget."""
    total = 0.0
    for r in range(n + 1):
        x = np.ones(n, dtype=np.int8)
        x[:r] = -1
        even = np.zeros((n, n), dtype=float)
        odd = np.zeros((n, n), dtype=float)
        for offset, aell in enumerate(alpha):
            ell = offset + 1
            term = aell * standard_layer_operator(n, ell, x)
            if ell & 1:
                odd += term
            else:
                even += term
        kappa = float(np.sum(even * even) + np.sum(odd * odd))
        tau_kappa = 2.0 * float(np.sum(even * odd)) * (
            ((n - 2 * r) ** 2 - n) / (n * (n - 1))
        )
        total += math.comb(n, r) * (tau_kappa - lam * kappa)
    return total / (2**n)


def all_hafnians(matrix: np.ndarray) -> list[int]:
    """Return unsigned hafnian polynomial values for every vertex subset."""
    n = len(matrix)

    @lru_cache(maxsize=None)
    def haf(mask: int) -> int:
        if mask == 0:
            return 1
        if bin(mask).count("1") & 1:
            return 0
        i_bit = mask & -mask
        i = i_bit.bit_length() - 1
        rest = mask ^ i_bit
        ans = 0
        probe = rest
        while probe:
            j_bit = probe & -probe
            j = j_bit.bit_length() - 1
            ans += int(matrix[i, j]) * haf(rest ^ j_bit)
            probe ^= j_bit
        return ans

    return [haf(mask) for mask in range(1 << n)]


def root_mass(matrix: np.ndarray, alpha: np.ndarray) -> tuple[float, float]:
    """Return exact cut average of kappa and its parity-only lower floor."""
    n = len(matrix)
    haf = all_hafnians(matrix)
    by_layer = [0] * len(alpha)
    for mask, value in enumerate(haf):
        size = bin(mask).count("1")
        if size % 2 == 0 and size // 2 < len(alpha):
            by_layer[size // 2] += value * value
    actual = 0.0
    floor = 0.0
    for ell, aell in enumerate(alpha):
        m = matching_count(n, ell)
        actual += aell * aell * by_layer[ell] / (m * m)
        floor += aell * aell * math.comb(n, 2 * ell) / (m * m)
    return actual, floor


def standard_root_mass(matrix: np.ndarray, alpha: np.ndarray) -> tuple[float, float]:
    """Standard-module root mass; exact multiple of the scalar hafnian sum."""
    scalar_actual, scalar_floor = root_mass(
        matrix, np.concatenate((np.asarray([0.0]), alpha))
    )
    factor = (len(matrix) - 1) ** 2
    return factor * scalar_actual, factor * scalar_floor


def cap(matrix: np.ndarray) -> int:
    n = len(matrix)
    best = 0
    for mask in range(1 << (n - 1)):
        x = np.ones(n, dtype=np.int64)
        for i in range(1, n):
            if mask >> (i - 1) & 1:
                x[i] = -1
        best = max(best, abs(int(x @ matrix @ x // 2)))
    return best


def saved_matrices() -> list[tuple[str, np.ndarray]]:
    rows: list[tuple[str, np.ndarray]] = []
    for n in range(3, 11):
        path = ROOT / f"computations/results/m{n}_minimizer_orbits.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text())
        for entry in payload["classes"]:
            rows.append(
                (
                    f"m{n}_class{entry['class']}",
                    np.asarray(entry["representative_matrix"], dtype=np.int8),
                )
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=10)
    args = parser.parse_args()
    output = []
    for name, matrix in saved_matrices():
        n = len(matrix)
        if n > args.max_n:
            continue
        for L in range(1, n // 2 + 1):
            lam, _, alpha = matching_path(n, L)
            budget = code_budget(n, alpha, lam)
            root, floor = root_mass(matrix, alpha)
            output.append(
                {
                    "name": name,
                    "n": n,
                    "L": L,
                    "cap": cap(matrix),
                    "lambda": lam,
                    "root_mass": root,
                    "parity_floor": floor,
                    "code_budget": budget,
                    "bound_normalized": lam - budget / root,
                    "floor_bound_normalized": lam - budget / floor,
                }
            )
        for L in range(2, (n - 1) // 2 + 1):
            lam, _, alpha = standard_matching_path(n, L)
            budget = standard_code_budget(n, alpha, lam)
            root, floor = standard_root_mass(matrix, alpha)
            output.append(
                {
                    "name": name,
                    "family": "standard_module",
                    "n": n,
                    "L": L,
                    "cap": cap(matrix),
                    "lambda": lam,
                    "root_mass": root,
                    "parity_floor": floor,
                    "code_budget": budget,
                    "bound_normalized": lam - budget / root,
                    "floor_bound_normalized": lam - budget / floor,
                }
            )
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
