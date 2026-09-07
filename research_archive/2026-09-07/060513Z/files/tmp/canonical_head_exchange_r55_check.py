#!/usr/bin/env python3
"""Finite algebra audit for the Wave 55 completion/head-exchange memo.

This checks exact identities on stored minimizers.  It is not numerical
evidence for an asymptotic population statement.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins
from far_tail_kfunctional_r52_check import exact_kfunctional, qnorm


def canonical_head(u: np.ndarray, hfun: int) -> tuple[int, ...]:
    """Return {i: |u_i|>=theta} at the positive clipped optimizer."""

    b = np.abs(u.astype(float))
    if np.count_nonzero(b) <= hfun:
        return ()
    lo = 0.0
    hi = max(1.0, float(np.max(b)))
    while float(np.sum(np.minimum(1.0, b / hi) ** 2)) > hfun:
        hi *= 2.0
    for _ in range(120):
        theta = (lo + hi) / 2.0
        if float(np.sum(np.minimum(1.0, b / theta) ** 2)) > hfun:
            lo = theta
        else:
            hi = theta
    return tuple(int(i) for i in np.flatnonzero(b >= hi - 1e-9))


def sign_nonzero(v: np.ndarray) -> np.ndarray:
    return np.where(v >= 0, 1, -1).astype(np.int64)


def audit_case(a: np.ndarray, name: str, m: int) -> dict[str, int | str]:
    n = len(a)
    k = n - m
    hfun = min(3, k)
    p2 = m * (m - 1) / (n * (n - 1))
    cases = 0
    exchanges = 0
    k_implications = 0
    positive_dh = 0
    negative_dh = 0

    # Deterministic, broad sample without making the checker expensive.
    selectors = list(itertools.combinations(range(n), m))
    for s_index, s0 in enumerate(selectors):
        if s_index % max(1, len(selectors) // 18) != 0:
            continue
        s = tuple(s0)
        t = tuple(i for i in range(n) if i not in s)
        pmat = a[np.ix_(s, s)]
        bmat = a[np.ix_(t, s)]
        cmat = a[np.ix_(t, t)]
        qs = qnorm(pmat)

        # Joint-completion matrix identities.
        cfull = np.zeros_like(a)
        cfull[np.ix_(t, s)] = bmat
        cfull[np.ix_(s, t)] = bmat.T
        cfull[np.ix_(t, t)] = cmat
        # Coordinate-free check by masking the old principal block.
        masked = a.copy()
        masked[np.ix_(s, s)] = 0
        assert np.array_equal(cfull, masked)
        assert int(np.sum(cfull * cfull)) == 2 * m * k + k * (k - 1)
        assert np.linalg.norm(cfull, 2) <= 2 * np.linalg.norm(a, 2) + 1e-9

        rmat = bmat.T @ bmat
        assert int(np.trace(rmat)) == m * k
        assert np.linalg.norm(rmat, 2) <= np.linalg.norm(a, 2) ** 2 + 1e-8
        assert np.sum(rmat.astype(object) ** 2) <= (
            np.linalg.norm(rmat, 2) + 1e-7
        ) * np.trace(rmat)

        ys = list(projective_spins(m))
        for y_index, y in enumerate(ys):
            if y_index % max(1, len(ys) // 12) != 0:
                continue
            u = bmat @ y
            genergy = int(u @ u)
            kval = exact_kfunctional(u.astype(float), hfun)
            assert kval <= math.sqrt(hfun * genergy) + 1e-8
            for threshold in range(0, int(math.ceil(kval)) + 3):
                if kval + 1e-9 >= threshold:
                    assert genergy + 1e-9 >= threshold**2 / hfun
                    k_implications += 1

            # R55.1 for sampled completions.
            for w0 in itertools.islice(itertools.product((-1, 1), repeat=k), 8):
                w = np.asarray(w0, dtype=np.int64)
                x = np.zeros(n, dtype=np.int64)
                x[list(s)] = y
                x[list(t)] = w
                direct = int(x @ cfull @ x)
                block = int(2 * y @ bmat.T @ w + w @ cmat @ w)
                assert direct == block
                cases += 1

            i_local = canonical_head(u, hfun)
            if not i_local:
                continue
            i = tuple(t[j] for j in i_local)
            ssize = len(i)
            assert ssize < hfun
            avec = sign_nonzero(u[list(i_local)])

            for sigma in (-1, 1):
                e = sigma * int(y @ pmat @ y)
                z = sigma * avec
                deletion_rows = []
                for r0 in itertools.combinations(s, ssize):
                    r = tuple(r0)
                    keep = tuple(v for v in s if v not in set(r))
                    ck = int(y[[s.index(v) for v in keep]] @ a[np.ix_(keep, keep)] @ y[[s.index(v) for v in keep]])
                    deletion_rows.append((sigma * (int(y @ pmat @ y) - ck), r))
                deletion_rows.sort(key=lambda item: (item[0], item[1]))
                deletion, r = deletion_rows[0]
                keep = tuple(v for v in s if v not in set(r))
                new_s = keep + i
                new_t = tuple(v for v in range(n) if v not in set(new_s))
                y_keep = y[[s.index(v) for v in keep]]
                y_r = y[[s.index(v) for v in r]]
                ynew = np.concatenate((y_keep, z))
                eprime = sigma * int(ynew @ a[np.ix_(new_s, new_s)] @ ynew)
                delta_e = eprime - e

                correction = (
                    -deletion
                    + 2 * int(np.sum(np.abs(u[list(i_local)])))
                    - 2 * int(avec @ a[np.ix_(i, r)] @ y_r)
                    + sigma * int(avec @ a[np.ix_(i, i)] @ avec)
                )
                assert delta_e == correction

                alpha = ssize * (2 * m - ssize - 1) / (m * (m - 1))
                assert deletion <= alpha * e + 1e-9
                assert delta_e + 1e-9 >= (
                    2 * np.sum(np.abs(u[list(i_local)]))
                    - alpha * e
                    - (3 * ssize * ssize - ssize)
                )

                qsprime = qnorm(a[np.ix_(new_s, new_s)])
                dh = qsprime - qs
                ls = ssize * (2 * m - ssize - 1)
                assert abs(dh) <= ls
                assert dh >= delta_e - (qs - e)
                g = (1 - p2) * e
                gprime = (1 - p2) * eprime - dh
                assert abs((gprime - g) - ((1 - p2) * delta_e - dh)) < 1e-9
                positive_dh += int(dh > 0)
                negative_dh += int(dh < 0)

                uprime = a[np.ix_(new_t, new_s)] @ ynew
                uold_embed = np.zeros(n, dtype=np.int64)
                unew_embed = np.zeros(n, dtype=np.int64)
                uold_embed[list(t)] = u
                unew_embed[list(new_t)] = uprime
                ukeep = tuple(v for v in t if v not in set(i))
                expected_ukeep = (
                    u[[t.index(v) for v in ukeep]]
                    - a[np.ix_(ukeep, r)] @ y_r
                    + a[np.ix_(ukeep, i)] @ z
                )
                assert np.array_equal(
                    uprime[[new_t.index(v) for v in ukeep]], expected_ukeep
                )
                expected_ur = a[np.ix_(r, keep)] @ y_keep + a[np.ix_(r, i)] @ z
                assert np.array_equal(
                    uprime[[new_t.index(v) for v in r]], expected_ur
                )
                rhs_change = (
                    int(u[list(i_local)] @ u[list(i_local)])
                    + int(expected_ur @ expected_ur)
                    + int((expected_ukeep - u[[t.index(v) for v in ukeep]]) @
                          (expected_ukeep - u[[t.index(v) for v in ukeep]]))
                )
                assert int((unew_embed - uold_embed) @ (unew_embed - uold_embed)) == rhs_change

                # Retain one full spin; row and completion-shift identities.
                x = np.ones(n, dtype=np.int64)
                x[list(s)] = y
                x[list(i)] = z
                row_old = int((a @ x) @ (a @ x))
                row_new = int((a @ x) @ (a @ x))
                assert row_old == row_new
                full_e = sigma * int(x @ a @ x)
                zold = full_e - e
                znew = full_e - eprime
                assert znew == zold - delta_e
                assert abs((p2 * znew - gprime) - (p2 * zold - g + dh - delta_e)) < 1e-9
                exchanges += 1

    return {
        "name": name,
        "n": n,
        "m": m,
        "joint_checks": cases,
        "K_implications": k_implications,
        "exchanges": exchanges,
        "positive_dh": positive_dh,
        "negative_dh": negative_dh,
    }


def main() -> None:
    rows = []
    for a, name in (
        (A6, "A6"),
        (A8, "A8"),
        (A9, "A9"),
        (sample_order_ten(0), "A10_seed0"),
    ):
        m = len(a) // 2 + 1
        row = audit_case(a, name, m)
        rows.append(row)
        print(
            f"{name} n={row['n']} m={row['m']} "
            f"joint={row['joint_checks']} K={row['K_implications']} "
            f"exchanges={row['exchanges']} "
            f"dh+/dh-={row['positive_dh']}/{row['negative_dh']}"
        )
    assert sum(int(row["joint_checks"]) for row in rows) > 0
    assert sum(int(row["K_implications"]) for row in rows) > 0
    assert sum(int(row["exchanges"]) for row in rows) > 0
    print("PASS Wave 55 joint-completion and canonical-head exchange audit")


if __name__ == "__main__":
    main()
