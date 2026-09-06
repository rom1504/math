#!/usr/bin/env python3
"""Explore high-replica absolute mixed pressure on A6/A8/A9."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a):
    xs = np.asarray(list(projective_spins(len(a))), dtype=np.int64)
    return int(np.max(np.abs(np.einsum("bi,ij,bj->b", xs, a, xs))))


def geometry(a, m):
    n = len(a)
    z = np.asarray(list(projective_spins(n)), dtype=np.int64)
    sels = list(itertools.combinations(range(n), m))
    qn = qnorm(a)
    parent_e = np.einsum("bi,ij,bj->b", z, a, z)
    row = np.einsum("bi,ij,bj->b", z, a @ a, z)
    deficits = np.empty((len(z), len(sels)), dtype=np.int16)
    for j, ss in enumerate(sels):
        s = list(ss)
        child = a[np.ix_(s, s)]
        qs = qnorm(child)
        cs = np.einsum("bi,ij,bj->b", z[:, s], child, z[:, s])
        deficits[:, j] = qs - np.abs(cs)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    h = (p**1.5 - p2) * qn
    return z, sels, qn, parent_e, row, deficits, h


def audit(name, a, m, uvals=(0.25, 0.5, 1, 2, 4, 8)):
    z, sels, qn, pe, row, d, h = geometry(a, m)
    n = len(a)
    c2 = row <= 2 * n * (n - 1)
    print(f"\n{name} m={m} ncenters={len(z)} C2={int(c2.sum())} H={h:.9f}")
    print("deficit histogram all pairs", sorted(Counter(map(int, d.ravel())).items()))
    for u in uvals:
        theta = u / h
        q = np.mean(np.exp(-theta * d), axis=1)
        elig = np.flatnonzero(c2)
        qm = float(q[elig].max())
        maxim = elig[np.isclose(q[elig], qm, rtol=1e-12, atol=1e-14)]
        kmax = u + math.log(qm)
        attrs = Counter((int(abs(pe[i])), int(row[i]), int(np.mean(d[i] <= h) * len(sels))) for i in maxim)
        prs = []
        for r in (1, 2, 4, 8, 16, 32, 64):
            vals = r * (u + np.log(q[elig]))
            vmax = float(vals.max())
            pr = (vmax + math.log(float(np.mean(np.exp(vals - vmax))))) / r
            prs.append(pr)
        print(
            f"u={u:g} maxK={kmax:.8f} mult={len(maxim)} "
            f"attrs(|E|,R,hardcount)={dict(attrs)} P_r="
            + ",".join(f"{x:.5f}" for x in prs)
        )


def main():
    for name, a, m in (("A6", A6, 5), ("A8", A8, 6), ("A9", A9, 7)):
        audit(name, a, m)


if __name__ == "__main__":
    main()
