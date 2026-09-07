#!/usr/bin/env python3
"""Exact finite checks for the Wave 16 common-mosaic response dual.

The script uses integer arithmetic throughout.  Fractional mixed matrices are
stored after multiplication by two, so their Boolean norms are also checked
exactly.  Complete minimizer lists are obtained by enumerating every
first-row-positive switching representative and expanding all switchings.
"""

from functools import lru_cache
from itertools import product

import numpy as np


def spins(n: int) -> np.ndarray:
    """One representative of every projective Boolean state."""
    return np.array(
        [(1,) + tail for tail in product((-1, 1), repeat=n - 1)],
        dtype=np.int64,
    )


def energies(matrix: np.ndarray) -> np.ndarray:
    z = spins(len(matrix))
    return np.einsum("bi,ij,bj->b", z, matrix, z)


def qnorm(matrix: np.ndarray) -> int:
    return int(np.max(np.abs(energies(matrix))))


def gauge_signing(n: int, mask: int) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=np.int64)
    matrix[0, 1:] = matrix[1:, 0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            matrix[i, j] = matrix[j, i] = 1 if (mask >> bit) & 1 else -1
            bit += 1
    return matrix


@lru_cache(None)
def all_minimizers(n: int):
    best = n * (n - 1)
    representatives = []
    free_edges = (n - 1) * (n - 2) // 2
    for mask in range(1 << free_edges):
        matrix = gauge_signing(n, mask)
        value = qnorm(matrix)
        if value < best:
            best, representatives = value, [matrix]
        elif value == best:
            representatives.append(matrix)

    labelled = {}
    for matrix in representatives:
        for tail in product((-1, 1), repeat=n - 1):
            switch = np.array((1,) + tail, dtype=np.int64)
            switched = matrix * np.outer(switch, switch)
            labelled[switched.tobytes()] = switched
    return best, tuple(labelled.values())


A9 = np.array(
    [
        [0, 1, -1, 1, 1, 1, 1, -1, -1],
        [1, 0, 1, -1, -1, -1, 1, -1, -1],
        [-1, 1, 0, 1, 1, 1, 1, 1, -1],
        [1, -1, 1, 0, 1, 1, 1, -1, 1],
        [1, -1, 1, 1, 0, 1, -1, 1, -1],
        [1, -1, 1, 1, 1, 0, -1, -1, -1],
        [1, 1, 1, 1, -1, -1, 0, 1, 1],
        [-1, -1, 1, -1, 1, -1, 1, 0, -1],
        [-1, -1, -1, 1, -1, -1, 1, -1, 0],
    ],
    dtype=np.int64,
)


A8 = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, -1, 1, 1, -1, -1],
        [1, 1, 0, 1, -1, 1, -1, -1],
        [1, -1, 1, 0, -1, -1, -1, 1],
        [1, 1, -1, -1, 0, -1, 1, -1],
        [1, 1, 1, -1, -1, 0, 1, 1],
        [1, -1, -1, -1, 1, 1, 0, 1],
        [1, -1, -1, 1, -1, 1, 1, 0],
    ],
    dtype=np.int64,
)


A5 = np.array(
    [
        [0, -1, 1, -1, 1],
        [-1, 0, -1, 1, 1],
        [1, -1, 0, 1, 1],
        [-1, 1, 1, 0, 1],
        [1, 1, 1, 1, 0],
    ],
    dtype=np.int64,
)


A7 = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, 1, -1, -1, 1],
        [1, 1, 0, 1, -1, 1, -1],
        [1, 1, 1, 0, 1, -1, -1],
        [1, -1, -1, 1, 0, -1, -1],
        [1, -1, 1, -1, -1, 0, -1],
        [1, 1, -1, -1, -1, -1, 0],
    ],
    dtype=np.int64,
)


def embedded(size: int, block, local: np.ndarray) -> np.ndarray:
    matrix = np.zeros((size, size), dtype=np.int64)
    matrix[np.ix_(block, block)] = local
    return matrix


def check_a9():
    block = (0, 1, 2, 4, 5, 6)
    singletons = (3, 7, 8)
    d6 = A9[np.ix_(block, block)]
    cross = A9.copy()
    cross[np.ix_(block, block)] = 0
    q6, minimizers6 = all_minimizers(6)

    assert qnorm(A9) == 24
    assert qnorm(d6) == 22
    assert qnorm(cross) == 22
    assert q6 == 10 and len(minimizers6) == 384
    x = qnorm(d6) - q6
    assert x == 12

    # Pure game: exactly 40 completions have norm 24, so its response gap is 12.
    best = 10**9
    attaining = []
    for local in minimizers6:
        hybrid = cross + embedded(9, block, local)
        value = qnorm(hybrid)
        if value < best:
            best, attaining = value, [(local, hybrid)]
        elif value == best:
            attaining.append((local, hybrid))
    assert best == 24 and len(attaining) == 40
    pure_gap = x + best - 24
    assert pure_gap == 12

    # Every optimal pure completion has a persistent common ground.  At such
    # a state the old and new internal energies agree, so the response change
    # is entirely the falling local deficit baseline.
    e_a = energies(A9)
    e_c = energies(cross)
    e_d = energies(embedded(9, block, d6))
    persistent_counts = []
    for local, hybrid in attaining:
        e_h = energies(hybrid)
        e_g = energies(embedded(9, block, local))
        profiles = []
        for sigma in (-1, 1):
            mask = (sigma * e_a == 24) & (sigma * e_h == 24)
            for index in np.flatnonzero(mask):
                profiles.append(
                    (
                        int(sigma * e_c[index]),
                        int(sigma * e_d[index]),
                        int(sigma * e_g[index]),
                    )
                )
        assert profiles
        assert set(profiles) <= {(14, 10, 10), (18, 6, 6)}
        persistent_counts.append(len(profiles))
    assert min(persistent_counts) >= 4 and max(persistent_counts) == 13

    # Mixed primal upper certificate.  The average of these two exact local
    # minimizers gives a fractional hybrid of norm 20.  We multiply by two.
    ga = np.array(
        [
            [0, -1, 1, 1, -1, -1],
            [-1, 0, -1, 1, -1, 1],
            [1, -1, 0, 1, 1, 1],
            [1, 1, 1, 0, -1, 1],
            [-1, -1, 1, -1, 0, 1],
            [-1, 1, 1, 1, 1, 0],
        ],
        dtype=np.int64,
    )
    gb = np.array(
        [
            [0, -1, 1, -1, 1, -1],
            [-1, 0, -1, 1, 1, -1],
            [1, -1, 0, 1, -1, -1],
            [-1, 1, 1, 0, -1, -1],
            [1, 1, -1, -1, 0, -1],
            [-1, -1, -1, -1, -1, 0],
        ],
        dtype=np.int64,
    )
    minimizer_keys = {matrix.tobytes() for matrix in minimizers6}
    assert ga.tobytes() in minimizer_keys and gb.tobytes() in minimizer_keys
    doubled_fractional = 2 * cross + embedded(9, block, ga + gb)
    assert qnorm(doubled_fractional) == 40

    # Mixed dual lower certificate.  The local restrictions coincide while
    # the orientations are opposite, so every possible local matrix cancels.
    z_plus = np.array([1, 1, -1, -1, -1, 1, -1, -1, -1], dtype=np.int64)
    z_minus = np.array([1, 1, -1, 1, -1, 1, -1, 1, 1], dtype=np.int64)
    assert np.array_equal(z_plus[list(block)], z_minus[list(block)])
    assert z_plus @ cross @ z_plus == 18
    assert -(z_minus @ cross @ z_minus) == 22
    assert z_plus @ A9 @ z_plus == 16
    assert -(z_minus @ A9 @ z_minus) == 24
    for local in minimizers6:
        local_plus = z_plus[list(block)] @ local @ z_plus[list(block)]
        local_minus = -(z_minus[list(block)] @ local @ z_minus[list(block)])
        assert local_plus + local_minus == 0
        # The expected block deficit difference is exactly X=12.
        delta_d_sum = (22 - z_plus[list(block)] @ d6 @ z_plus[list(block)])
        delta_d_sum += (22 + z_minus[list(block)] @ d6 @ z_minus[list(block)])
        delta_g_sum = (10 - local_plus) + (10 - local_minus)
        assert delta_d_sum - delta_g_sum == 24

    mixed_norm = qnorm(doubled_fractional) // 2
    average_global_slack = ((24 - 16) + (24 - 24)) // 2
    block_margin = 12
    mixed_gap = x + mixed_norm - 24
    assert mixed_norm == 20
    assert average_global_slack == 4
    assert mixed_gap == block_margin - average_global_slack == 8

    print(
        "A9:",
        {
            "pure_norm": best,
            "pure_gap": pure_gap,
            "optimal_pure_completions": len(attaining),
            "mixed_norm": mixed_norm,
            "mixed_gap": mixed_gap,
            "dual_block_margin": block_margin,
            "dual_global_slack": average_global_slack,
        },
    )


def check_a8():
    blocks = ((0, 3, 4, 5), (1, 2, 6, 7))
    d_blocks = [A8[np.ix_(block, block)] for block in blocks]
    cross = A8.copy()
    for block in blocks:
        cross[np.ix_(block, block)] = 0
    q4, minimizers4 = all_minimizers(4)

    assert qnorm(A8) == 20
    assert q4 == 8 and len(minimizers4) == 48
    assert [qnorm(local) for local in d_blocks] == [12, 12]
    assert qnorm(cross) == 16
    x_parts = [qnorm(local) - q4 for local in d_blocks]
    assert x_parts == [4, 4]
    x = sum(x_parts)

    # Pure response problem.
    best = 10**9
    count = 0
    for g1 in minimizers4:
        for g2 in minimizers4:
            hybrid = cross + embedded(8, blocks[0], g1) + embedded(8, blocks[1], g2)
            value = qnorm(hybrid)
            if value < best:
                best, count = value, 1
            elif value == best:
                count += 1
    assert best == 20 and count == 16
    pure_gap = x + best - 20
    assert pure_gap == 8

    # Mixed upper: 0 is in each local-minimizer convex hull (average G and -G),
    # so the cross matrix itself is feasible and has norm 16.
    minimizer_keys = {matrix.tobytes() for matrix in minimizers4}
    assert all((-matrix).tobytes() in minimizer_keys for matrix in minimizers4)

    # Mixed dual: two opposite orientations with block restrictions equal up
    # to sign cancel every local quadratic form and expose cross energy 16.
    z_plus = np.ones(8, dtype=np.int64)
    z_minus = np.array([1, -1, -1, 1, 1, 1, -1, -1], dtype=np.int64)
    assert z_plus @ cross @ z_plus == 16
    assert -(z_minus @ cross @ z_minus) == 16
    assert z_plus @ A8 @ z_plus == 12
    assert -(z_minus @ A8 @ z_minus) == 20
    for block, old in zip(blocks, d_blocks):
        zp = z_plus[list(block)]
        zm = z_minus[list(block)]
        assert np.array_equal(zm, zp) or np.array_equal(zm, -zp)
        assert (zp @ old @ zp) - (zm @ old @ zm) == 0
        for local in minimizers4:
            assert (zp @ local @ zp) - (zm @ local @ zm) == 0

    mixed_norm = qnorm(cross)
    average_global_slack = ((20 - 12) + (20 - 20)) // 2
    mixed_gap = x + mixed_norm - 20
    assert mixed_norm == 16
    assert average_global_slack == 4
    assert mixed_gap == sum(x_parts) - average_global_slack == 4

    # Endpoint neutrality: enumerate all 4x4 positive/negative endpoint pairs.
    z = spins(8)
    e = energies(A8)
    positive = z[e == 20]
    negative = z[e == -20]
    assert len(positive) == len(negative) == 4
    directed_shores = 0
    for p in positive:
        for n in negative:
            sign = p * n
            shores = (tuple(np.flatnonzero(sign == 1)), tuple(np.flatnonzero(sign == -1)))
            assert tuple(map(len, shores)) == (4, 4)
            for child in shores:
                sibling = tuple(i for i in range(8) if i not in child)
                child_matrix = A8[np.ix_(child, child)]
                y = p[list(child)]
                h = int(y @ child_matrix @ y)
                ell = int(np.sum(np.abs(A8[np.ix_(sibling, child)] @ y)))
                child_norm = qnorm(child_matrix)
                decrement = qnorm(A8) - child_norm
                assert child_norm == 8 and h == 0 and ell == 10
                for sigma in (-1, 1):
                    capacity = max(0, 2 * ell - (child_norm - sigma * h))
                    residual = max(0, capacity - decrement)
                    assert capacity == decrement == 12 and residual == 0
                directed_shores += 1
    assert directed_shores == 32

    print(
        "A8:",
        {
            "pure_norm": best,
            "pure_gap": pure_gap,
            "optimal_pure_completions": count,
            "mixed_norm": mixed_norm,
            "mixed_gap": mixed_gap,
            "dual_block_margins": x_parts,
            "dual_global_slack": average_global_slack,
            "endpoint_pairs": len(positive) * len(negative),
            "all_root_endpoint_residuals": 0,
        },
    )


def check_timing_walls():
    # A5: the endpoint pair in (10.518) has zero residual on both shores.
    assert qnorm(A5) == 8
    p = np.ones(5, dtype=np.int64)
    n = np.array([1, 1, 1, 1, -1], dtype=np.int64)
    assert p @ A5 @ p == 8 and n @ A5 @ n == -8
    for child in ((0, 1, 2, 3), (4,)):
        sibling = tuple(i for i in range(5) if i not in child)
        child_matrix = A5[np.ix_(child, child)]
        y = p[list(child)]
        h = int(y @ child_matrix @ y)
        ell = int(np.sum(np.abs(A5[np.ix_(sibling, child)] @ y)))
        decrement = 8 - qnorm(child_matrix)
        for sigma in (-1, 1):
            capacity = max(0, 2 * ell - (qnorm(child_matrix) - sigma * h))
            assert max(0, capacity - decrement) == 0

    # Exact signs of the temporal quantities used in the A5/A7 walls.
    assert 8 - 64 * np.sqrt(5) / 25 > 0
    assert 6 - 24 * np.sqrt(15) / 25 > 0
    x1 = 10 - 108 * np.sqrt(42) / 49
    x2 = (108 * np.sqrt(42) - 90 * np.sqrt(35)) / 49 - 2
    assert x1 < 0 < x2 and x1 + x2 < 0

    # A7 deterministic suffix: 7 -> 6 -> 5 are all exact minimizers.
    child6 = A7[1:, 1:]
    assert qnorm(A7) == 18 and qnorm(child6) == 10
    assert all(
        qnorm(np.delete(np.delete(child6, i, axis=0), i, axis=1)) == 8
        for i in range(6)
    )
    print(
        "timing walls:",
        {
            "A5_first_demand_positive": True,
            "A5_completion_demand_positive": True,
            "A7_first_atom_negative": True,
            "A7_suffix_atom_positive": True,
            "A7_root_sum_negative": True,
        },
    )


def main():
    check_a9()
    check_a8()
    check_timing_walls()
    print("all exact response-dual checks: PASS")


if __name__ == "__main__":
    main()
