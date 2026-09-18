#!/usr/bin/env python3
"""Exact small-order audit of the symplectic switching-cell construction."""

from __future__ import annotations

import itertools

import numpy as np


def dot2(x: int, y: int) -> int:
    return bin(x & y).count("1") & 1


def build(k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = k * k
    K = np.empty((n, n), dtype=np.int64)
    for u, v, x, y in itertools.product(range(k), repeat=4):
        K[u * k + v, x * k + y] = (-1) ** (dot2(v, x) ^ dot2(u, y))
    A = K - np.eye(n, dtype=np.int64)
    D = np.zeros_like(A)
    for u in range(k):
        sl = slice(u * k, (u + 1) * k)
        D[sl, sl] = A[sl, sl]
    return A, A - D, D


def spins(n: int) -> np.ndarray:
    masks = np.arange(1 << (n - 1), dtype=np.uint64)[:, None]
    bits = (masks >> np.arange(n - 1, dtype=np.uint64)) & 1
    return np.concatenate((np.ones((len(masks), 1), dtype=np.int64), 1 - 2 * bits.astype(np.int64)), axis=1)


def audit(k: int) -> None:
    A, C, D = build(k)
    n = k * k
    Z = spins(n)
    ec = np.einsum("bi,ij,bj->b", Z, C, Z, optimize=True)
    ed = np.einsum("bi,ij,bj->b", Z, D, Z, optimize=True)
    ea = ec + ed
    qc = int(np.max(np.abs(ec)))
    qa = int(np.max(np.abs(ea)))
    assert qa == n * (k + 1)
    assert qc == n * k
    assert int(np.max(np.abs(ed))) == n * (k - 1)
    print(f"k={k} n={n} Q(A)={qa} Q(C)={qc} Q(D)={int(np.max(np.abs(ed)))}")
    print("C extrema", int(ec.max()), int(ec.min()), "counts", int(np.sum(ec == ec.max())), int(np.sum(ec == ec.min())))

    # G consists of arbitrary signs constant on the k character cells.  For
    # every oriented C-ground, compute its best opposite state in its orbit.
    cell_switches = []
    for mask in range(1 << k):
        signs = np.array([(-1 if (mask >> u) & 1 else 1) for u in range(k)], dtype=np.int64)
        cell_switches.append(np.repeat(signs, k))
    gaps = []
    variances = []
    witnesses = []
    for idx in np.flatnonzero(np.abs(ec) == qc):
        z = Z[idx]
        sigma = 1 if ec[idx] > 0 else -1
        orbit = np.array([sigma * int((z * g) @ C @ (z * g)) for g in cell_switches])
        gap = qc + int(orbit.min())
        gaps.append(gap)
        variances.append(int(orbit @ orbit) // len(orbit))
        if gap == min(gaps):
            witnesses.append((idx, sigma, tuple(map(int, orbit))))
    print("ground orbit gaps", sorted(set(gaps)), "best", min(gaps), "worst", max(gaps))
    print("ground orbit variances", sorted(set(variances)))
    print("best witness orbit", min(witnesses, key=lambda row: qc + min(row[2])))

    pos = np.flatnonzero(ec == qc)
    neg = np.flatnonzero(ec == -qc)
    assert np.all(ed[pos] == -n) and np.all(ed[neg] == -n)
    # Every opposite exact cross-ground pair is orthogonal, hence its full-A
    # mutual cross-Gram and paired raw-shore functional are zero.
    grams = Z[pos] @ A @ Z[neg].T
    assert np.all(grams == 0)
    full_pos_deficit = qa - (qc - n)
    full_neg_deficit = qa - (qc + n)
    assert (full_pos_deficit, full_neg_deficit) == (2 * n, 0)
    paired = np.maximum(np.abs(grams) - full_pos_deficit - full_neg_deficit, 0)
    assert np.all(paired == 0)
    # At t=0 the negative C-ground, under negative orientation, has D
    # response +n, so the active defect is exactly zero.
    response = max(int(ed[pos].max()), int((-ed[neg]).max()))
    assert response == n
    print("active response", response, "Gamma", max(-response, 0), "paired-functional max", int(paired.max()))


def main() -> None:
    audit(2)
    audit(4)


if __name__ == "__main__":
    main()
