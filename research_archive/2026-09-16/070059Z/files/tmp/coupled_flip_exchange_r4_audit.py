#!/usr/bin/env python3
"""Exact finite checks for wave-4 coupled-flip exchange results."""

from __future__ import annotations

import itertools

import numpy as np


def spins(n: int) -> np.ndarray:
    return np.array(
        [[1, *tail] for tail in itertools.product((1, -1), repeat=n - 1)],
        dtype=int,
    )


def norm_single(a: np.ndarray) -> int:
    x = spins(len(a))
    return int(np.abs(np.einsum("bi,ij,bj->b", x, a, x)).max() // 2)


def norm_doubled(a: np.ndarray) -> int:
    return 2 * norm_single(a)


def flip_edges(a: np.ndarray, edges: set[tuple[int, int]]) -> np.ndarray:
    out = a.copy()
    for i, j in edges:
        out[i, j] *= -1
        out[j, i] *= -1
    return out


def order_five_check() -> None:
    a = np.array(
        [
            [0, -1, 1, -1, 1],
            [-1, 0, -1, 1, 1],
            [1, -1, 0, 1, 1],
            [-1, 1, 1, 0, 1],
            [1, 1, 1, 1, 0],
        ],
        dtype=int,
    )
    s = np.array([True, True, False, False, False])
    t = ~s
    r = a.sum(axis=1)
    beta = a[np.ix_(s, t)].sum(axis=1)
    c = int(beta.sum())
    d_block = a[np.ix_(s, s)]
    e_block = a[np.ix_(t, t)]
    q = norm_doubled(a)
    g = 4 * c
    h_d = int(d_block.sum())
    h_e = int(e_block.sum())
    gamma = norm_doubled(e_block) - h_e
    d = q - norm_doubled(e_block)
    ell = 2 * int(np.abs(beta).sum()) - gamma
    covariance = int(r[s] @ beta)
    assert (q, tuple(r), tuple(beta), c, g) == (8, (0, 0, 2, 2, 4), (1, 1), 2, 8)
    assert (h_d, h_e, gamma, d, ell, covariance) == (-2, 6, 0, 2, 4, 0)
    assert norm_single(d_block) + norm_single(e_block) == norm_single(a) == 4
    print("n=5", q, r.tolist(), beta.tolist(), c, g, h_d, h_e, gamma, d, ell, covariance)


def order_nine_matrix() -> np.ndarray:
    b = np.ones((4, 4), dtype=int) - np.eye(4, dtype=int)
    d = np.array(
        [
            [0, 1, -1, 1, 1],
            [1, 0, 1, 1, -1],
            [-1, 1, 0, 1, 1],
            [1, 1, 1, 0, 1],
            [1, -1, 1, 1, 0],
        ],
        dtype=int,
    )
    c = np.array(
        [
            [1, 1, -1, 1, -1],
            [-1, 1, 1, -1, -1],
            [1, -1, 1, -1, 1],
            [-1, -1, 1, 1, -1],
        ],
        dtype=int,
    )
    return np.block([[b, c], [c.T, d]])


def coupled_switching_check(a: np.ndarray) -> None:
    shore_s = set(range(4))
    shore_t = set(range(4, 9))
    histogram: dict[int, int] = {}
    for rmask in range(1 << 4):
        chosen_r = {i for i in shore_s if (rmask >> i) & 1}
        internal_cut = {
            (min(i, j), max(i, j))
            for i in chosen_r
            for j in shore_s - chosen_r
        }
        for umask in range(1 << 5):
            chosen_u = {4 + j for j in range(5) if (umask >> j) & 1}
            f = internal_cut | {(i, j) for i in chosen_r for j in chosen_u}
            h = {(i, j) for i in chosen_r for j in shore_t - chosen_u}
            switched = np.array([-1 if i in chosen_r else 1 for i in range(9)])
            af = flip_edges(a, f)
            ah = flip_edges(a, h)
            assert np.array_equal(af * switched[:, None] * switched[None, :], ah)
            assert norm_single(af) == norm_single(ah)
            value = norm_single(af)
            histogram[value] = histogram.get(value, 0) + 1
    assert sum(histogram.values()) == 512
    assert min(histogram) == 12
    print("n=9 coupled histogram", dict(sorted(histogram.items())))


def qnorm_submatrix(e: np.ndarray) -> int:
    if len(e) <= 1:
        return 0
    return norm_doubled(e)


def order_nine_all_endpoint_check(a0: np.ndarray) -> None:
    n = len(a0)
    x_table = spins(n)
    scores = np.einsum("bi,ij,bj->b", x_table, a0, x_table) // 2
    records: list[tuple[int, int, int, int]] = []
    for endpoint_mask, x in enumerate(x_table):
        for orientation in (1, -1):
            if orientation * scores[endpoint_mask] != 12:
                continue
            a = orientation * a0 * x[:, None] * x[None, :]
            r = a.sum(axis=1)
            assert np.all(r >= 0)
            for mask in range(1, 1 << (n - 1)):
                s = np.array([(mask >> i) & 1 for i in range(n)], dtype=bool)
                t = ~s
                cross = a[np.ix_(s, t)]
                beta = cross.sum(axis=1)
                c = int(beta.sum())
                if c <= 0:
                    continue
                g = 4 * c
                covariance = int(r[s] @ beta)
                e = a[np.ix_(t, t)]
                gamma = qnorm_submatrix(e) - int(e.sum())
                ell = 2 * int(np.abs(beta).sum()) - gamma
                records.append((covariance, c, g, ell))

    nonpositive = [z for z in records if z[0] <= 0]
    bad_eighth = [z for z in records if z[3] < z[2] / 8]
    negative_ell = [z for z in records if z[3] < 0]
    assert len(nonpositive) == 18
    assert all(cov == 0 and ell / g == 0.5 for cov, _, g, ell in nonpositive)
    assert len(bad_eighth) == 636
    assert min(cov / (8 * c) for cov, c, _, _ in bad_eighth) == 0.25
    assert len(negative_ell) == 4
    assert min(cov / (8 * c) for cov, c, _, _ in negative_ell) == 0.625
    print(
        "n=9 endpoint audit",
        "nonpositive=", len(nonpositive),
        "ell<g/8=", len(bad_eighth),
        "ell<0=", len(negative_ell),
    )


def margin_identity_check(a: np.ndarray) -> None:
    n = len(a)
    q = norm_doubled(a)
    r = a.sum(axis=1)
    for mask in range(1, (1 << n) - 1):
        s = np.array([(mask >> i) & 1 for i in range(n)], dtype=bool)
        t = ~s
        d = a[np.ix_(s, s)]
        beta = a[np.ix_(s, t)].sum(axis=1)
        c = int(beta.sum())
        covariance = int(r[s] @ beta)
        rs = r[s]
        scale = n - 1
        delta = q / 2 - c
        rhs = (
            float(rs @ rs)
            - float(rs @ d @ rs) / scale
            - scale * delta
        ) / 2
        assert covariance >= rhs - 1e-9
    print("n=9 margin inequalities checked on all nontrivial shores")


def main() -> None:
    order_five_check()
    a9 = order_nine_matrix()
    assert norm_single(a9) == 12
    coupled_switching_check(a9)
    order_nine_all_endpoint_check(a9)
    margin_identity_check(a9)


if __name__ == "__main__":
    main()
