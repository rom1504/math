#!/usr/bin/env python3
"""Finite audit of the OpenAI Boolean moving-projection kernel on cut codes.

Scratch only.  The ambient Hamming length is N=binom(v,2).  For a Boolean
moving-projection path (N,k,L), construct the overlap kernel K(t), evaluate
the internal augmented-cut-code Gram remainder, and test the rooted annular
criterion derived from equation (27) of Chapter 2 of ten-proofs-oai.pdf.
"""

from __future__ import annotations

import argparse
import itertools
import math

import numpy as np
import scipy.linalg


def subsets(n: int, r: int) -> list[tuple[int, ...]]:
    return list(itertools.combinations(range(n), r))


def incidence_up(n: int, r: int) -> np.ndarray:
    """U: V_r -> V_{r+1}."""
    lo = subsets(n, r)
    hi = subsets(n, r + 1)
    lo_index = {s: j for j, s in enumerate(lo)}
    out = np.zeros((len(hi), len(lo)), dtype=float)
    for i, t in enumerate(hi):
        for q in range(r + 1):
            s = t[:q] + t[q + 1 :]
            out[i, lo_index[s]] = 1.0
    return out


def primitive_basis(n: int, k: int) -> np.ndarray:
    if k == 0:
        return np.ones((1, 1))
    dmat = incidence_up(n, k - 1).T
    return scipy.linalg.null_space(dmat, rcond=1e-11)


def embeddings(n: int, k: int, L: int) -> dict[int, np.ndarray]:
    q = primitive_basis(n, k)
    out = {k: q}
    cur = q
    for level in range(k, L):
        cur = incidence_up(n, level) @ cur
        # Formula (15): ||U^(i-k)h||^2.
        i = level + 1
        norm2 = (
            math.factorial(i - k)
            * math.factorial(n - 2 * k)
            / math.factorial(n - i - k)
        )
        out[i] = cur / math.sqrt(norm2)
    for phi in out.values():
        err = np.linalg.norm(phi.T @ phi - np.eye(phi.shape[1]), ord=2)
        if err > 2e-8:
            raise AssertionError(f"embedding isometry error {err}")
    return out


def jacobi(n: int, k: int, L: int) -> tuple[float, np.ndarray]:
    size = L - k + 1
    jmat = np.zeros((size, size), dtype=float)
    for level in range(k, L):
        c = ((level - k + 1) * (n - level - k)) / (
            n * math.sqrt((level + 1) * (n - level))
        )
        q = level - k
        jmat[q, q + 1] = jmat[q + 1, q] = c
    vals, vecs = np.linalg.eigh(jmat)
    lam = float(vals[-1])
    v = vecs[:, -1]
    if v.sum() < 0:
        v = -v
    return lam, v


def kernel_by_distance(n: int, k: int, L: int) -> tuple[float, int, np.ndarray]:
    phis = embeddings(n, k, L)
    lam, v = jacobi(n, k, L)
    levels = list(range(k, L + 1))
    dims = np.array([math.comb(n, i) for i in levels], dtype=float)
    w = np.sqrt(dims) * v
    if np.any(w <= 0):
        raise AssertionError("Perron/dimension weights are not positive")
    amplitudes2 = w / w.sum()
    rank = phis[k].shape[1]
    kvals = np.empty(n + 1, dtype=float)
    for dist in range(n + 1):
        flipped = set(range(dist))
        overlap = np.zeros((rank, rank), dtype=float)
        for amp2, level in zip(amplitudes2, levels):
            ss = subsets(n, level)
            signs = np.array(
                [(-1.0) ** len(flipped.intersection(s)) for s in ss]
            )
            phi = phis[level]
            overlap += amp2 * (phi.T @ (signs[:, None] * phi))
        kvals[dist] = float(np.sum(overlap * overlap))
    if abs(kvals[0] - rank) > 2e-7:
        raise AssertionError((kvals[0], rank))
    return lam, rank, kvals


def augmented_cut_distances(v: int) -> list[int]:
    edges = list(itertools.combinations(range(v), 2))
    out: list[int] = []
    # Fix the first vertex sign to +1 to quotient the global spin flip.
    for tail in itertools.product((-1, 1), repeat=v - 1):
        x = (1,) + tail
        base = [x[i] * x[j] for i, j in edges]
        for sigma in (-1, 1):
            out.append(sum(1 for z in base if sigma * z == -1))
    assert len(out) == 2**v
    return out


def audit(v: int, k: int, L: int) -> None:
    n = math.comb(v, 2)
    lam, rank, kvals = kernel_by_distance(n, k, L)
    dists = augmented_cut_distances(v)
    correlations = np.array([1 - 2 * d / n for d in range(n + 1)])
    remainder = sum((correlations[d] - lam) * kvals[d] for d in dists)
    print(f"vertices={v} ambient={n} k={k} L={L}")
    print(f"lambda={lam:.12g} rank={rank} remainder={remainder:.12g}")
    print("d K(d)")
    for d, val in enumerate(kvals):
        print(d, f"{val:.12g}")
    best = None
    # If max correlation <=s for the antipodal code, all correlations lie in [-s,s].
    for q in range(n + 1):
        s = abs(correlations[q])
        if s >= lam:
            continue
        central = [kvals[d] for d in range(n + 1) if abs(correlations[d]) <= s + 1e-12]
        if not central:
            continue
        kmin = min(central)
        lhs = (lam - s) ** 2 * len(dists) * kmin**2
        rhs = (1 - lam) * rank * remainder
        if lhs > rhs + 1e-8:
            best = max(best or 0.0, s)
        print(
            f"s={s:.8g} kmin={kmin:.6g} criterion={lhs-rhs:.6g} "
            f"{'FORCES' if lhs > rhs + 1e-8 else ''}"
        )
    print("best forced s", best)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("v", type=int)
    ap.add_argument("k", type=int)
    ap.add_argument("L", type=int)
    args = ap.parse_args()
    audit(args.v, args.k, args.L)
