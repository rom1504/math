#!/usr/bin/env python3
"""Numerical/exact-energy audit of the soft selected-child identity.

All energies are enumerated exactly as integers.  Only log-sum-exp evaluation
is floating point.  The factorization identity is also checked directly from
the cavity fields to near machine precision.
"""

from __future__ import annotations

from itertools import combinations, product
import math

import numpy as np
from scipy.special import logsumexp


QMIN = {1: 0, 2: 2, 3: 6, 4: 8, 5: 8, 6: 10, 7: 18, 8: 20, 9: 24}


def mat_from_code(n: int, code: int) -> np.ndarray:
    a = np.zeros((n, n), dtype=np.int64)
    a[0, 1:] = a[1:, 0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            a[i, j] = a[j, i] = 1 if (code >> k) & 1 else -1
            k += 1
    return a


A9 = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int64)

A8 = np.array([
    [0, 1, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, -1, 1, 1, -1, -1],
    [1, 1, 0, 1, -1, 1, -1, -1],
    [1, -1, 1, 0, -1, -1, -1, 1],
    [1, 1, -1, -1, 0, -1, 1, -1],
    [1, 1, 1, -1, -1, 0, 1, 1],
    [1, -1, -1, -1, 1, 1, 0, 1],
    [1, -1, -1, 1, -1, 1, 1, 0],
], dtype=np.int64)

A7 = np.array([
    [0, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, 1, -1, -1, 1],
    [1, 1, 0, 1, -1, 1, -1],
    [1, 1, 1, 0, 1, -1, -1],
    [1, -1, -1, 1, 0, -1, -1],
    [1, -1, 1, -1, -1, 0, -1],
    [1, 1, -1, -1, -1, -1, 0],
], dtype=np.int64)

A5 = np.array([
    [0, -1, 1, -1, 1],
    [-1, 0, -1, 1, 1],
    [1, -1, 0, 1, 1],
    [-1, 1, 1, 0, 1],
    [1, 1, 1, 1, 0],
], dtype=np.int64)

NAMED = {
    "A2": mat_from_code(2, 0),
    "A4": mat_from_code(4, 1),
    "A5": A5,
    "A6": mat_from_code(6, 220),
    "A7": A7,
    "A8": A8,
    "A9": A9,
}


def spins(n: int) -> np.ndarray:
    if n == 0:
        return np.zeros((1, 0), dtype=np.int64)
    return np.asarray([(1,) + t for t in product((-1, 1), repeat=n - 1)], dtype=np.int64)


def projective_energies(a: np.ndarray) -> np.ndarray:
    x = spins(len(a))
    return np.einsum("bi,ij,bj->b", x, a, x, optimize=True)


def oriented_energies(a: np.ndarray) -> np.ndarray:
    e = projective_energies(a)
    return np.concatenate((e, -e))


def p0(a: np.ndarray, beta: float) -> float:
    e = oriented_energies(a).astype(float)
    return float((logsumexp(beta * e) - math.log(len(e))) / beta)


def child(a: np.ndarray, i: int) -> np.ndarray:
    ids = [j for j in range(len(a)) if j != i]
    return a[np.ix_(ids, ids)]


def cavity_kappa(a: np.ndarray, i: int, beta: float) -> float:
    ids = [j for j in range(len(a)) if j != i]
    c = a[np.ix_(ids, ids)]
    y = spins(len(c))
    ec = np.einsum("bi,ij,bj->b", y, c, y, optimize=True).astype(float)
    f = (y @ a[i, ids]).astype(float)
    # Sum over sigma; cosh does not depend on sigma.  Stable log numerator.
    oriented_ec = np.concatenate((ec, -ec))
    oriented_f = np.concatenate((f, f))
    logw = beta * oriented_ec
    log_num = logsumexp(logw + np.logaddexp(2 * beta * oriented_f, -2 * beta * oriented_f) - math.log(2))
    return float((log_num - logsumexp(logw)) / beta)


def audit_matrix(name: str, a: np.ndarray, beta: float) -> tuple[float, float]:
    r = len(a)
    if r <= 1:
        return math.inf, math.inf
    delta_q = QMIN[r] - QMIN[r - 1]
    pb = p0(a, beta) - QMIN[r]
    threshold = delta_q + pb * (1 - ((r - 1) / r) ** 1.5)
    kappas = []
    norm_drops = []
    for i in range(r):
        c = child(a, i)
        k_pressure = p0(a, beta) - p0(c, beta)
        k_cavity = cavity_kappa(a, i, beta)
        assert abs(k_pressure - k_cavity) < 2e-11, (name, i, beta, k_pressure, k_cavity)
        kappas.append(k_pressure)
        pc = p0(c, beta) - QMIN[r - 1]
        lhs = pc / ((r - 1) ** 1.5)
        rhs = pb / (r ** 1.5)
        norm_drops.append(rhs - lhs)
    return max(kappas) - delta_q, max(norm_drops)


def main():
    betas = (0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0)
    print("margin=max_i kappa_i-(q_r-q_{r-1}); softdrop=max_i(Phi_parent-Phi_child)")
    for beta in betas:
        print("beta", beta)
        required = []
        for name, a in NAMED.items():
            margin, softdrop = audit_matrix(name, a, beta)
            r = len(a)
            dvr = 1 / math.sqrt(r - 1) - 1 / math.sqrt(r)
            c_req = math.log(2) + beta * max(0.0, -softdrop) / dvr
            required.append(c_req)
            print(f"  {name}: margin={margin:+.9f} softdrop={softdrop:+.9f}")
        print(f"  named fixed-beta coefficient c required: {max(required):.9f}")

    # Every restriction of A9, grouped by order: worst best-child margins.
    for beta in (0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
        print("A9 restriction census beta", beta)
        all_required = []
        for r in range(3, 10):
            vals = []
            soft = []
            for ids in combinations(range(9), r):
                a = A9[np.ix_(ids, ids)]
                m, d = audit_matrix(str(ids), a, beta)
                vals.append(m)
                soft.append(d)
            dvr = 1 / math.sqrt(r - 1) - 1 / math.sqrt(r)
            c_req = math.log(2) + beta * max(0.0, -min(soft)) / dvr
            all_required.append(c_req)
            print(f"  r={r}: worst cavity margin={min(vals):+.9f}; worst selected softdrop={min(soft):+.9f}")
            print(f"       fixed-beta coefficient c required={c_req:.9f}")
        print(f"  all A9 restrictions coefficient c required: {max(all_required):.9f}")

    # Locate the A9 threshold(s) for the simpler sufficient cavity condition.
    grid = np.geomspace(1e-4, 20.0, 2000)
    margins = np.asarray([audit_matrix("A9", A9, float(b))[0] for b in grid])
    good = grid[margins >= -1e-10]
    print("A9 cavity-good beta interval on grid:", (float(good.min()), float(good.max())) if len(good) else None)
    j = int(np.argmax(margins))
    print("A9 best cavity margin/grid beta:", float(margins[j]), float(grid[j]))


if __name__ == "__main__":
    main()
