#!/usr/bin/env python3
"""Tiny exact wind-tunnel for the double-partial rank-one channel."""

from __future__ import annotations

import itertools
import math
import numpy as np


def signs(d: int) -> np.ndarray:
    return np.asarray(list(itertools.product((-1, 1), repeat=d)), dtype=np.int8)


def softmax(v: np.ndarray) -> np.ndarray:
    z = v - np.max(v)
    w = np.exp(z)
    return w / w.sum()


def kl(a: np.ndarray, b: np.ndarray) -> float:
    keep = a > 0
    return float(np.sum(a[keep] * (np.log(a[keep]) - np.log(b[keep]))))


def product_tensor(rows: list[np.ndarray]) -> np.ndarray:
    out = rows[0]
    for p in rows[1:]:
        out = np.multiply.outer(out, p)
    return out


def mf_reverse_kl(logq: np.ndarray, starts: list[list[np.ndarray]], rounds=1000):
    m = logq.ndim
    best = (math.inf, None)
    for init in starts:
        rows = [p.copy() for p in init]
        old = math.inf
        for _ in range(rounds):
            for i in range(m):
                # Contract every other row against its current factor.
                vals = np.empty(logq.shape[i])
                for a in range(logq.shape[i]):
                    sl = np.take(logq, a, axis=i)
                    others = rows[:i] + rows[i + 1 :]
                    vals[a] = float(np.sum(sl * product_tensor(others)))
                rows[i] = softmax(vals)
            pt = product_tensor(rows)
            val = float(np.sum(pt * (np.log(np.maximum(pt, 1e-300)) - logq)))
            if abs(old - val) < 1e-13:
                break
            old = val
        if val < best[0]:
            best = (val, rows)
    return best


def run(m: int, n: int, k: int, r: int, beta: float, lam: float):
    ell, s = m - k, n - r
    nscale = m + n
    u = beta / math.sqrt(nscale)
    row_words = signs(n)
    row_count = len(row_words)

    xs = np.asarray(
        [(sig,) * k + xi for sig in (-1, 1) for xi in itertools.product((-1, 1), repeat=ell)],
        dtype=np.int8,
    )
    ys = np.asarray(
        [(tau,) * r + eta for tau in (-1, 1) for eta in itertools.product((-1, 1), repeat=s)],
        dtype=np.int8,
    )
    qwords = np.unique(np.asarray([np.outer(x, y).ravel() for x in xs for y in ys]), axis=0)

    bridge_rows = np.asarray(list(itertools.product(range(row_count), repeat=m)), dtype=np.int16)
    bridges = row_words[bridge_rows].reshape((-1, m * n)).astype(np.float64)
    p = np.mean(np.exp(u * bridges @ qwords.T), axis=1) / math.cosh(u) ** (m * n)
    q = p ** (-lam)
    q /= q.sum()
    qt = q.reshape((row_count,) * m)

    sd = row_words[:, :r].sum(axis=1)
    nu = np.cosh(u * sd) ** (-lam)
    nu /= nu.sum()
    un = np.full(row_count, 1 / row_count)
    can_rows = [nu.copy() for _ in range(m)]
    comp_rows = [un.copy() for _ in range(k)] + [nu.copy() for _ in range(ell)]
    can = product_tensor(can_rows)
    comp = product_tensor(comp_rows)
    uniform = np.full_like(qt, 1 / qt.size)
    j = kl(can.ravel(), q)
    compkl = kl(comp.ravel(), q)
    fairkl = kl(uniform.ravel(), q)

    rng = np.random.default_rng(20260818)
    starts = [can_rows, comp_rows, [un.copy() for _ in range(m)]]
    for _ in range(20):
        starts.append([softmax(rng.normal(size=row_count)) for _ in range(m)])
    ikl, _ = mf_reverse_kl(np.log(qt), starts)
    print(
        dict(
            m=m,
            n=n,
            k=k,
            r=r,
            beta=beta,
            lam=lam,
            support=len(qwords),
            J=j,
            I_mf=ikl,
            J_minus_I=j - ikl,
            comparison=compkl,
            J_minus_comparison=j - compkl,
            fair=fairkl,
        )
    )


if __name__ == "__main__":
    for beta in (0.5, 1.0, 2.0, 4.0):
        run(4, 4, 2, 2, beta, 1.0)
