"""Exact all-diagonal audit of A9 x the order-two conference factor."""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
import numpy as np

from mixed_product_r17 import A9


def spins(n: int, projective: bool = False) -> np.ndarray:
    width = n - 1 if projective else n
    z = np.arange(1 << width, dtype=np.uint32)[:, None]
    b = ((z >> np.arange(width, dtype=np.uint32)) & 1).astype(np.int16)
    x = 2 * b - 1
    if projective:
        x = np.concatenate((np.ones((len(x), 1), dtype=np.int16), x), axis=1)
    return x


def main() -> None:
    # The small certificate used by the growing-Hadamard proof.  There are ten
    # positive and fifteen negative projective A9 grounds.  For every diagonal
    # completion, one opposite-sign pair has completed cross-Gram at least 15.
    all9p = spins(9, projective=True)
    a9e = np.sum((all9p @ A9) * all9p, axis=1)
    POS = all9p[a9e == 24]
    NEG = all9p[a9e == -24]
    assert (len(POS), len(NEG)) == (10, 15)
    cross_max_hist = Counter()
    for d in spins(9):
        P = A9 + np.diag(d)
        cross_max = max(abs(int(p @ P @ n)) for p in POS for n in NEG)
        cross_max_hist[cross_max] += 1
    print("max opposite-ground completed cross-Gram histogram", sorted(cross_max_hist.items()))
    assert min(cross_max_hist) == 15

    X = spins(9, projective=True)  # global antipodal quotient in first fibre
    Y = spins(9)
    # All 2^17 global projective states, indexed by (x,y).
    XX = np.repeat(X, len(Y), axis=0)
    YY = np.tile(Y, (len(X), 1))
    ax = np.sum((XX @ A9) * XX, axis=1).astype(np.int16)
    ay = np.sum((YY @ A9) * YY, axis=1).astype(np.int16)
    cross0 = (2 * np.sum((XX @ A9) * YY, axis=1)).astype(np.int16)
    coeff = (2 * XX * YY).astype(np.int16)

    DS = spins(9).astype(np.int16)
    hist = Counter()
    by_e = defaultdict(Counter)
    minima = {}
    # Process diagonal completions in chunks to keep memory modest.
    for e1, e2 in itertools.product((-1, 1), repeat=2):
        vals = np.full(len(DS), -1, dtype=np.int32)
        base = cross0.astype(np.int32) + e1 * ax.astype(np.int32) + e2 * ay.astype(np.int32)
        for start in range(0, len(DS), 64):
            dchunk = DS[start:start + 64]
            energies = base[:, None] + coeff.astype(np.int32) @ dchunk.T.astype(np.int32)
            vals[start:start + len(dchunk)] = np.max(np.abs(energies), axis=0)
        for d, q in zip(DS, vals):
            delta = int(d.sum())
            hist[int(q)] += 1
            by_e[(e1, e2)][(delta, int(q))] += 1
        minima[(e1, e2)] = (int(vals.min()), int(vals.max()))

    print("all 2048 completion Q histogram", sorted(hist.items()))
    print("per E min/max", minima)
    for e, h in by_e.items():
        print("E", e, "delta,Q counts", sorted(h.items()))
    assert sum(hist.values()) == 2048
    assert min(hist) > 18 * np.sqrt(17)  # order-18 conference spectral upper
    assert min(hist) == 78
    print("all mixed order-two products are rigorously nonminimal: PASS")

    H2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    H4 = np.kron(H2, H2)
    w = np.array([1, 1, 1, -1], dtype=np.int64)
    assert np.array_equal(H4 @ w, 2 * w)
    ratio = 78 / (18 * np.sqrt(18))
    assert ratio > 1
    print("growing Sylvester witness ratio", ratio, ": PASS")


if __name__ == "__main__":
    main()
