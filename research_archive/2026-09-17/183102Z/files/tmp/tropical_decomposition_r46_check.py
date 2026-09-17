#!/usr/bin/env python3
"""Exact checks for Wave 46C tropical principal decomposition.

All calculations are integer/exact except the optional finite-beta normalization
check, whose tolerance is explicitly asserted.
"""

from __future__ import annotations

import itertools
from collections import Counter

import numpy as np
from scipy.special import logsumexp


def spins(n: int) -> np.ndarray:
    """Projective spin representatives, first spin fixed to +1."""
    if n == 0:
        return np.ones((1, 0), dtype=np.int8)
    return np.asarray([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int8)


def edges(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def matrix_from_normalized_code(n: int, code: int) -> np.ndarray:
    """Switching-normalized signing: every edge from vertex zero is +1."""
    A = np.zeros((n, n), dtype=np.int8)
    bit = 0
    for i in range(1, n):
        A[0, i] = A[i, 0] = 1
    for i in range(1, n):
        for j in range(i + 1, n):
            A[i, j] = A[j, i] = 1 if ((code >> bit) & 1) else -1
            bit += 1
    return A


def raw_energies(A: np.ndarray) -> np.ndarray:
    X = spins(len(A)).astype(np.int64)
    return np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)


def qnorm(A: np.ndarray) -> int:
    return int(np.max(np.abs(raw_energies(A))))


def enumerate_normalized_minimizers(n: int) -> tuple[int, list[np.ndarray]]:
    """Enumerate one representative per switching orbit."""
    var_edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    nv = len(var_edges)
    count = 1 << nv
    X = spins(n).astype(np.int16)
    # Raw energy / 2 is sum over upper-triangle edges.
    fixed = np.sum(X[:, [0]] * X[:, 1:], axis=1, dtype=np.int16)
    products = np.asarray([X[:, i] * X[:, j] for i, j in var_edges], dtype=np.int16)
    best = 10**9
    best_codes: list[int] = []
    block = 4096
    shifts = np.arange(nv, dtype=np.uint64)
    for start in range(0, count, block):
        cc = np.arange(start, min(start + block, count), dtype=np.uint64)
        signs = (((cc[:, None] >> shifts[None, :]) & 1).astype(np.int16) * 2 - 1)
        half = fixed[None, :] + signs @ products
        qs = 2 * np.max(np.abs(half), axis=1)
        local = int(np.min(qs))
        if local < best:
            best = local
            best_codes = []
        if local == best:
            best_codes.extend(int(cc[j]) for j in np.flatnonzero(qs == best))
    return best, [matrix_from_normalized_code(n, c) for c in best_codes]


def selector_profile(A: np.ndarray, m: int) -> dict[str, object]:
    """Test the common-ground/tight-decomposition property for one slice."""
    n = len(A)
    X = spins(n).astype(np.int64)
    full_raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
    qn = int(np.max(np.abs(full_raw)))
    parent = [(x, 1 if e > 0 else -1) for x, e in zip(X, full_raw) if abs(int(e)) == qn]
    selectors = list(itertools.combinations(range(n), m))
    data: list[tuple[tuple[int, ...], int, list[tuple[int, ...]]]] = []
    for S in selectors:
        AS = A[np.ix_(S, S)].astype(np.int64)
        Y = spins(m).astype(np.int64)
        vals = np.einsum("bi,ij,bj->b", Y, AS, Y, optimize=True)
        qS = int(np.max(np.abs(vals)))
        # Oriented ground keys (sigma, all pair products on S).
        grounds: list[tuple[int, ...]] = []
        for y, val in zip(Y, vals):
            if abs(int(val)) == qS:
                orientations = (-1, 1) if int(val) == 0 else ((1,) if val > 0 else (-1,))
                for sigma in orientations:
                    key = (sigma,) + tuple(int(sigma * y[i] * y[j]) for i in range(m) for j in range(i + 1, m))
                    grounds.append(key)
        data.append((S, qS, grounds))
    qstar = max(qS for _, qS, _ in data)
    max_data = [(S, grounds) for S, qS, grounds in data if qS == qstar]

    lifted: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    common_parents = 0
    for x, sigma in parent:
        hit = False
        for S, grounds in max_data:
            key = (sigma,) + tuple(int(sigma * x[S[i]] * x[S[j]]) for i in range(m) for j in range(i + 1, m))
            if key in grounds:
                lifted.add((S, key))
                hit = True
        common_parents += int(hit)
    total_pairs = sum(len(grounds) for _, grounds in max_data)
    return {
        "q_n": qn,
        "m": m,
        "q_star": qstar,
        "max_selectors": len(max_data),
        "max_child_ground_pairs": total_pairs,
        "lifted_pairs": len(lifted),
        "common_parent_grounds": common_parents,
        "common_exists": common_parents > 0,
        "every_max_child_ground_lifts": len(lifted) == total_pairs,
    }


def actual_tropical_data(A: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Return E_d, a_d=min_S max external energy, and q_* for oriented cuts."""
    n = len(A)
    X = spins(n).astype(np.int64)
    states = [(sigma, x) for sigma in (-1, 1) for x in X]
    E = np.asarray([sigma * int(x @ A @ x) for sigma, x in states], dtype=np.int64)
    a = np.full(len(states), 10**9, dtype=np.int64)
    qstar = -1
    for S in itertools.combinations(range(n), m):
        S = tuple(S)
        AS = A[np.ix_(S, S)].astype(np.int64)
        qS = qnorm(AS)
        qstar = max(qstar, qS)
        groups: dict[tuple[int, ...], list[int]] = {}
        for j, (sigma, x) in enumerate(states):
            key = (sigma,) + tuple(int(sigma * x[S[i]] * x[S[k]]) for i in range(m) for k in range(i + 1, m))
            groups.setdefault(key, []).append(j)
        for inds in groups.values():
            first = inds[0]
            sigma, x = states[first]
            child = sigma * int(x[list(S)] @ AS @ x[list(S)])
            extmax = max(int(E[j] - child) for j in inds)
            a[inds] = np.minimum(a[inds], extmax)
    return E, a, qstar


def finite_beta_cancellation_check(A: np.ndarray, m: int, beta: float = 0.73) -> float:
    """Check q(S)h_S=e^{-F_S}/sum z and the resulting f exactly numerically."""
    n = len(A)
    X = spins(n).astype(np.int64)
    states = [(sigma, x) for sigma in (-1, 1) for x in X]
    E = np.asarray([sigma * int(x @ A @ x) for sigma, x in states], dtype=float)
    logZA = float(logsumexp(beta * E))
    eF_rows = []
    z_rows = []
    for S in itertools.combinations(range(n), m):
        S = tuple(S)
        groups: dict[tuple[int, ...], list[int]] = {}
        for j, (sigma, x) in enumerate(states):
            key = (sigma,) + tuple(int(sigma * x[S[i]] * x[S[k]]) for i in range(m) for k in range(i + 1, m))
            groups.setdefault(key, []).append(j)
        row = np.empty(len(states))
        for inds in groups.values():
            sigma, x = states[inds[0]]
            AS = A[np.ix_(S, S)].astype(np.int64)
            child = sigma * int(x[list(S)] @ AS @ x[list(S)])
            ext = E[inds] - child
            logKnorm = float(logsumexp(beta * ext) - (n - m) * np.log(2.0))
            row[inds] = np.exp(-logKnorm)
        eF_rows.append(row)
        z_rows.append(float(np.sum(np.exp(beta * E - logZA) * row)))
    eF = np.asarray(eF_rows)
    z = np.asarray(z_rows)
    q = z / np.sum(z)
    h = eF / z[:, None]
    f1 = q @ h
    f2 = np.sum(eF, axis=0) / np.sum(z)
    return float(np.max(np.abs(f1 - f2)))


def batched_common_exhaustion(n: int, qn: int, minimizers: list[np.ndarray]) -> dict[int, int]:
    """Fast all-class test of common tight decomposition at every slice."""
    X = spins(n).astype(np.int16)
    ee = edges(n)
    products = np.asarray([X[:, i] * X[:, j] for i, j in ee], dtype=np.int16).T
    subsets = [S for m in range(1, n) for S in itertools.combinations(range(n), m)]
    sizes = np.asarray([len(S) for S in subsets])
    masks = np.asarray(
        [[int(i in S and j in S) for i, j in ee] for S in subsets], dtype=np.int16
    )
    avec = np.asarray([[A[i, j] for i, j in ee] for A in minimizers], dtype=np.int16)
    failures = {m: 0 for m in range(1, n)}
    for start in range(0, len(minimizers), 32):
        av = avec[start : start + 32]
        raw = 2 * (av @ products.T)
        parent = np.abs(raw) == qn
        sigma = np.where(raw >= 0, 1, -1).astype(np.int16)
        internal = 2 * np.einsum("xe,be,se->bxs", products, av, masks, optimize=True)
        qs = np.max(np.abs(internal), axis=1)
        for m in range(1, n):
            ids = np.flatnonzero(sizes == m)
            qstar = np.max(qs[:, ids], axis=1)
            eligible = qs[:, ids] == qstar[:, None]
            oriented = internal[:, :, ids] * sigma[:, :, None]
            hit = np.any(
                parent[:, :, None]
                & eligible[:, None, :]
                & (oriented == qstar[:, None, None]),
                axis=(1, 2),
            )
            failures[m] += int(np.sum(~hit))
    return failures


def main() -> None:
    summary = []
    stronger_failures = []
    # Detailed enumeration through six, then a faster all-class exhaustion
    # at seven and eight (the largest switching-normalized space practical
    # here without a separate classification).
    for n in range(2, 7):
        qn, minimizers = enumerate_normalized_minimizers(n)
        failures = []
        for ai, A in enumerate(minimizers):
            for m in range(1, n):
                z = selector_profile(A, m)
                if not z["common_exists"]:
                    failures.append((ai, z))
                if not z["every_max_child_ground_lifts"]:
                    stronger_failures.append((n, ai, z))
                # The direct line computation is redundant and substantially
                # more expensive; audit it on one representative per order.
                if ai == 0:
                    E, a, qstar = actual_tropical_data(A, m)
                    assert int(np.max(E - a)) == qstar
                    common_line = np.any((E == qn) & ((E - a) == qstar))
                    assert common_line == z["common_exists"]
        assert not failures
        summary.append({"n": n, "q_n": qn, "switching_classes": len(minimizers), "slice_failures": len(failures)})

    for n in (7, 8):
        qn, minimizers = enumerate_normalized_minimizers(n)
        failures = batched_common_exhaustion(n, qn, minimizers)
        assert not any(failures.values())
        # Audit the full tropical line identity on a representative.
        for m in range(1, n):
            E, a, qstar = actual_tropical_data(minimizers[0], m)
            assert int(np.max(E - a)) == qstar
            assert np.any((E == qn) & ((E - a) == qstar))
        summary.append(
            {
                "n": n,
                "q_n": qn,
                "switching_classes": len(minimizers),
                "failures_by_m": failures,
            }
        )

    # Also cover the stored order-eight/order-nine exact minimizers.
    import sys
    sys.path.insert(0, "/home/math/quadra/tmp")
    from endpoint_transport_r25 import A8, A9
    large = []
    for name, A in (("A8", A8), ("A9", A9)):
        rows = [selector_profile(A, m) for m in range(1, len(A))]
        assert all(z["common_exists"] for z in rows)
        large.append((name, rows))

    err = finite_beta_cancellation_check(A8, 4)
    assert err < 2e-13
    print("normalized exact-minimizer exhaustion", summary)
    print("stronger all-pairs failures", len(stronger_failures), stronger_failures[:8])
    print("stored large minimizers", large)
    print("finite-beta cancellation max error", err)
    print("PASS tropical_decomposition_r46_check")


if __name__ == "__main__":
    main()
