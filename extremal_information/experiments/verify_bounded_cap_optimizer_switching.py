#!/usr/bin/env python3
"""Exact finite checks for bounded_cap_optimizer_switching.md."""

from __future__ import annotations

from itertools import combinations, product
from math import comb
import random

import numpy as np


def cube(k: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=k)), dtype=np.int64)


def all_children(k: int):
    edges = list(combinations(range(k), 2))
    for signs in product((-1, 1), repeat=len(edges)):
        a = np.zeros((k, k), dtype=np.int64)
        for sign, (i, j) in zip(signs, edges):
            a[i, j] = a[j, i] = sign
        yield a


def hvals(a: np.ndarray, xs: np.ndarray) -> np.ndarray:
    return np.einsum("bi,ij,bj->b", xs, a, xs, optimize=True) // 2


def huge_common_fibre_check() -> None:
    for k in range(2, 6):
        xs = cube(k)
        u = np.ones(k, dtype=np.int64)
        uid = int(np.flatnonzero(np.all(xs == u, axis=1))[0])
        count = 0
        worst_parent = 0
        for a in all_children(k):
            hv = hvals(a, xs)
            q = int(np.max(np.abs(hv)))
            if q > 2 * k ** 1.5 + 1e-9:
                continue
            if int(hv[uid]) != int(np.max(hv)):
                continue
            g = np.abs(xs @ u)
            assert int(hv[uid] + g[uid]) == int(np.max(hv + g))
            # Max over the one new spin of absolute parent energy.
            parent_q = int(np.max(np.abs(hv[:, None] + (xs @ u)[:, None] * np.asarray([-1, 1]))))
            assert parent_q <= q + k
            worst_parent = max(worst_parent, parent_q)
            count += 1
        lower = 2 ** (comb(k, 2) - k)
        assert count + 1e-9 >= lower
        print(
            f"k={k}: common-u bounded children={count}, "
            f"theorem lower={lower:g}, worst parent cap={worst_parent}"
        )


def exact_future(B: np.ndarray, C: np.ndarray, xs: np.ndarray) -> np.ndarray:
    ys = cube(B.shape[1])
    hc = hvals(C, ys)
    values = xs @ B @ ys.T + hc[None, :]
    return np.max(values, axis=1)


def witness_cover_inequality_check() -> None:
    rng = random.Random(266101)
    k, m = 4, 3
    xs = cube(k)
    children = list(all_children(k))
    h = np.stack([hvals(a, xs) for a in children])
    futures = []
    for _ in range(7):
        b = np.asarray(
            [[rng.choice((-1, 1)) for _ in range(m)] for _ in range(k)],
            dtype=np.int64,
        )
        c = np.zeros((m, m), dtype=np.int64)
        for i, j in combinations(range(m), 2):
            c[i, j] = c[j, i] = rng.choice((-1, 1))
        futures.append(exact_future(b, c, xs))

    # A deliberately incomplete reusable dictionary.
    uids = [0, 3, 5, 10]
    full = np.stack([np.max(h + g[None, :], axis=1) for g in futures], axis=1)
    trunc = np.stack(
        [np.max(h[:, uids] + g[None, uids], axis=1) for g in futures], axis=1
    )
    error = full - trunc
    assert np.all(error >= 0)
    tau_raw = int(np.max(error))

    evals = h[:, uids]
    for i in range(len(children)):
        for j in range(i + 1, len(children)):
            response_distance = int(np.max(np.abs(full[i] - full[j])))
            eval_distance = int(np.max(np.abs(evals[i] - evals[j])))
            assert response_distance <= eval_distance + 2 * tau_raw

            response_proj = (int(np.max(full[i] - full[j])) - int(np.min(full[i] - full[j]))) / 2
            eval_proj = (int(np.max(evals[i] - evals[j])) - int(np.min(evals[i] - evals[j]))) / 2
            assert response_proj <= eval_proj + tau_raw + 1e-9

    print(
        f"witness-cover core inequalities passed: children={len(children)}, "
        f"futures={len(futures)}, |U|={len(uids)}, raw tau={tau_raw}"
    )


def main() -> None:
    huge_common_fibre_check()
    witness_cover_inequality_check()
    print("bounded-cap optimizer-switching checks passed")


if __name__ == "__main__":
    main()
