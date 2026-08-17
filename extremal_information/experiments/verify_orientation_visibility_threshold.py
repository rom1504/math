#!/usr/bin/env python3
"""Finite checks for drafts/orientation_visibility_threshold.md.

The asymptotic construction is probabilistic.  This script fixes seeds,
checks the exact Boolean caps at n=4,16, and verifies the scalar envelopes
used in the proof.  It is diagnostic; the theorem uses the norm-existence
lemma rather than these finite samples.
"""

from __future__ import annotations

from itertools import product
from math import floor, sqrt

import numpy as np

from verify_bcx_two_port_holonomy import regular_hadamard


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int8)


def exact_caps(m_walsh: int, seed: int) -> tuple[int, int, float]:
    q, h, _ = regular_hadamard(m_walsh)
    n = q * q
    shore = floor(n ** 0.75)
    a = shore / n
    rng = np.random.default_rng(seed)
    b = np.where(rng.random((n, shore)) < (1 + a) / 2, 1, -1).astype(np.int8)
    residual_norm = float(
        np.linalg.norm(b.astype(float) - a * np.ones((n, shore)), ord=2)
    )

    xs, ys = cube(n), cube(shore)
    old = np.einsum("bi,ij,bj->b", xs, h, xs, optimize=True) // 2
    sums = np.sum(ys, axis=1, dtype=np.int64)
    clique = (sums * sums - shore) // 2

    qp = qm = 0
    for start in range(0, len(xs), 1024):
        xchunk = xs[start : start + 1024].astype(np.int64)
        cross = xchunk @ b.astype(np.int64) @ ys.astype(np.int64).T
        hp = old[start : start + len(xchunk), None]
        qp = max(qp, int(np.max(np.abs(hp + cross + clique[None, :]))))
        qm = max(qm, int(np.max(np.abs(-hp + cross + clique[None, :]))))

    assert abs(qp - qm) <= shore * (shore - 1)
    print(
        f"n={n}, m={shore}, seed={seed}: Q+={qp}, Q-={qm}, "
        f"gap={qp-qm}, ||B-aJ||={residual_norm:.6f}"
    )
    return qp, qm, residual_norm


def scalar_envelopes() -> None:
    grid = np.linspace(-1, 1, 1601)
    for lam in (0.1, 0.5, 0.9, 1.0):
        p, s = np.meshgrid(grid, grid, indexing="ij")
        plus = 0.5 - p * p + lam * p * s + 0.5 * lam * s * s
        minus = 0.5 + lam * (-p * s - 0.5 * s * s)
        assert float(np.max(plus)) <= 1.25 + 1e-9
        assert float(np.max(minus)) <= 1.0 + 1e-9
    print("scalar orientation envelopes: PASS")


def main() -> None:
    scalar_envelopes()
    exact_caps(1, 20260817)
    exact_caps(2, 20260817)
    print("orientation visibility threshold checks: PASS")


if __name__ == "__main__":
    main()
