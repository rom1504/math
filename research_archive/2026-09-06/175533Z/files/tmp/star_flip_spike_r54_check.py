#!/usr/bin/env python3
"""Exact star-flip coverage diagnostics for Wave 54.

Checks the simultaneous-flip certificate on positive internal stars of
spiky child grounds in the stored exact minimizers.  All arithmetic is
integer and exhaustive at these orders.
"""

from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A as A10
from envelope_block_cover_r27 import A8, A9


def oriented_states(a: np.ndarray):
    n = len(a)
    states = []
    q = 0
    for tail in itertools.product((-1, 1), repeat=n - 1):
        x = np.asarray((1,) + tail, dtype=np.int64)
        e = int(x @ a @ x)
        q = max(q, abs(e))
        states.append((x, e))
    out = []
    for x, e in states:
        for sigma in (-1, 1):
            out.append((sigma, x, q - sigma * e))
    return q, out


def spikiest_child(a: np.ndarray, m: int):
    n = len(a)
    best = None
    for s0 in itertools.combinations(range(n), m):
        s = np.asarray(s0, dtype=int)
        b = a[np.ix_(s, s)]
        vals = []
        q_s = 0
        for tail in itertools.product((-1, 1), repeat=m - 1):
            y = np.asarray((1,) + tail, dtype=np.int64)
            e = int(y @ b @ y)
            q_s = max(q_s, abs(e))
            vals.append((y, e))
        for y, e in vals:
            if abs(e) != q_s:
                continue
            sigma = 1 if e > 0 else -1
            r = sigma * y * (b @ y)
            rec = (int(np.max(r)), int(np.sum(r * r)), tuple(s0), y, sigma, r)
            if best is None or rec[:2] > best[:2]:
                best = rec
    assert best is not None
    return best


def analyze(name: str, a: np.ndarray, m: int):
    q, states = oriented_states(a)
    rmax, internal, s0, y, sigma_y, r = spikiest_child(a, m)
    i_local = int(np.argmax(r))
    i = s0[i_local]
    positive = []
    for j_local, j in enumerate(s0):
        if j == i:
            continue
        child_sign = int(a[i, j] * sigma_y * y[i_local] * y[j_local])
        if child_sign == 1:
            positive.append((min(i, j), max(i, j)))
    h = len(positive)
    assert 2 * h - (m - 1) == rmax

    # For every subset F of P, verify some oriented state obeys
    # Delta + 4 sum_{e in F} s_e <= 0.  Record the best deficit and the
    # number of negative positive-star edges for the full-star flip.
    full_best = None
    full_witnesses = []
    by_k = []
    for k in range(h + 1):
        worst_best_margin = None
        union_witnesses = set()
        for idxs in itertools.combinations(range(h), k):
            f = [positive[t] for t in idxs]
            best_margin = None
            witnesses = []
            for state_idx, (sig, x, delta) in enumerate(states):
                sf = sum(int(a[u, v] * sig * x[u] * x[v]) for u, v in f)
                margin = delta + 4 * sf
                if best_margin is None or margin < best_margin:
                    best_margin = margin
                    witnesses = [(state_idx, delta, sf)]
                elif margin == best_margin:
                    witnesses.append((state_idx, delta, sf))
            assert best_margin is not None and best_margin <= 0
            worst_best_margin = (
                best_margin if worst_best_margin is None
                else max(worst_best_margin, best_margin)
            )
            union_witnesses.update(t[0] for t in witnesses if t[1] + 4 * t[2] <= 0)
            if k == h:
                full_best = best_margin
                full_witnesses = witnesses
        by_k.append((k, worst_best_margin, len(union_witnesses)))

    # Near-ground states and their negative density on P.
    near = []
    for state_idx, (sig, x, delta) in enumerate(states):
        neg = sum(
            int(a[u, v] * sig * x[u] * x[v] == -1)
            for u, v in positive
        )
        if delta <= 4 * h:
            near.append((delta, neg, state_idx))

    print(f"{name},m={m}: q={q}, S={s0}, child sigma={sigma_y}, r={tuple(map(int,r))}")
    print(f"  I={internal}, spike vertex={i}, rmax={rmax}, |P|={h}")
    print(f"  by k (k,worst optimum,union of optimum witnesses)={by_k}")
    print(f"  full P best margin={full_best}, best witnesses={full_witnesses[:8]}")
    print(f"  near count={len(near)}, max negatives on P={max(t[1] for t in near)}, profile={sorted(set((d,z) for d,z,_ in near))}")


def main():
    analyze("A8", A8, 5)
    analyze("A9", A9, 6)
    analyze("A10", A10, 6)
    print("PASS star-flip certificates")


if __name__ == "__main__":
    main()
