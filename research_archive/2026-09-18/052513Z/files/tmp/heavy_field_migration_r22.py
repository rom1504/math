#!/usr/bin/env python3
"""Exact audits for the Wave 22 edge-cube/star-cover identities.

This is deliberately a finite verifier.  The all-edge-set theorem in the
memo is algebraic and does not require enumerating 2^N perturbations.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np

from check_response_dual_r16 import A8, A9, spins


def oriented_states(A: np.ndarray):
    """One copy of each augmented-cut word (sigma, projective x)."""
    n = len(A)
    X0 = spins(n)
    X = np.tile(X0, (2, 1))
    sigma = np.repeat(np.array((-1, 1), dtype=np.int64), len(X0))
    edge_list = [(i, j) for i in range(n) for j in range(i + 1, n)]
    edge_sign = np.column_stack(
        [sigma * A[i, j] * X[:, i] * X[:, j] for i, j in edge_list]
    )
    energy = 2 * edge_sign.sum(axis=1)
    q = int(energy.max())
    deficit = q - energy
    negative = edge_sign == -1
    return X, sigma, edge_list, edge_sign, negative, energy, q, deficit


def star_audit(A: np.ndarray, name: str, beta: float) -> None:
    n = len(A)
    m = n - 1
    X, sigma, edges, signs, negative, energy, q, deficit = oriented_states(A)
    edge_index = {edge: j for j, edge in enumerate(edges)}
    N = len(edges)
    assert (2 * N - q) % 4 == 0
    R = (2 * N - q) // 4

    # Complete-energy/coset identity, state by state.
    neg_count = negative.sum(axis=1)
    assert np.array_equal(deficit, 4 * (neg_count - R))
    assert int(neg_count.min()) == R

    weights = np.exp(-beta * deficit)
    D = float(weights.sum())
    all_subset_min_margins = []
    all_subset_min_witness_deficits = []
    star_cover_density = []
    star_bcs = []
    average_competitor_D = []

    for i in range(n):
        H = np.array(
            [edge_index[tuple(sorted((i, j)))] for j in range(n) if j != i],
            dtype=int,
        )
        d_i = negative[:, H].sum(axis=1)
        h_i = signs[:, H].sum(axis=1)
        rho = d_i - deficit // 4
        assert np.array_equal(deficit + 2 * h_i, 2 * m - 4 * rho)

        # The restricted balls cover every row replacement.  Record both the
        # best distance margin and the least deficit of any certifying state.
        competitor_D = []
        local_margins = []
        local_witness_deficits = []
        cover_multiplicities = []
        for bits in itertools.product((0, 1), repeat=m):
            target = np.asarray(bits, dtype=bool)
            distances = np.logical_xor(negative[:, H], target).sum(axis=1)
            margins = rho - distances
            certifying = margins >= 0
            assert np.any(certifying)
            local_margins.append(int(margins.max()))
            local_witness_deficits.append(int(deficit[certifying].min()))
            cover_multiplicities.append(int(certifying.sum()))

            selected_sign_sum = signs[:, H][:, target].sum(axis=1)
            shifted_deficit = deficit + 4 * selected_sign_sum
            ds = float(np.exp(-beta * shifted_deficit).sum())
            assert ds >= 1.0 - 1e-12
            competitor_D.append(ds)

        local_margins = np.asarray(local_margins)
        all_subset_min_margins.extend(local_margins.tolist())
        all_subset_min_witness_deficits.extend(local_witness_deficits)

        # Direct Hellinger/field formula and the uniform row-replacement
        # transform.  Complementary star perturbations are switch-equivalent.
        bc = float((weights * np.exp(-2 * beta * h_i)).sum() / D)
        avg_ds = float(np.mean(competitor_D))
        transformed = D * math.cosh(2 * beta) ** m * bc
        assert math.isclose(avg_ds, transformed, rel_tol=1e-11, abs_tol=1e-11)
        for mask in range(1 << m):
            assert math.isclose(
                competitor_D[mask], competitor_D[((1 << m) - 1) ^ mask],
                rel_tol=1e-11,
                abs_tol=1e-11,
            )
        assert bc + 1e-12 >= 1.0 / (D * math.cosh(2 * beta) ** m)
        star_bcs.append(bc)
        average_competitor_D.append(avg_ds)

        # A volume union bound is another exact consequence of the same
        # restricted cover; negative radii certify no target.
        volume_sum = sum(
            sum(math.comb(m, k) for k in range(min(int(rad), m) + 1))
            for rad in rho
            if rad >= 0
        )
        assert volume_sum >= 2**m
        star_cover_density.append(volume_sum / 2**m)
        assert min(cover_multiplicities) >= 1

    # Ground-field migration census: which coordinates can carry the largest
    # positive field among exact ground states?
    ground = deficit == 0
    ground_argmax = Counter()
    ground_profiles = Counter()
    for row in signs[ground]:
        hs = np.zeros(n, dtype=int)
        for eidx, (i, j) in enumerate(edges):
            hs[i] += row[eidx]
            hs[j] += row[eidx]
        maximum = int(hs.max())
        ground_profiles[tuple(sorted(map(int, hs), reverse=True))] += 1
        for i in np.flatnonzero(hs == maximum):
            ground_argmax[int(i)] += 1

    print(
        name,
        f"beta={beta}",
        {
            "N": N,
            "Q": q,
            "R": R,
            "oriented_states": len(deficit),
            "ground_states": int(ground.sum()),
            "star_best_margin_hist": sorted(Counter(all_subset_min_margins).items()),
            "max_min_witness_deficit": max(all_subset_min_witness_deficits),
            "star_cover_density_range": (
                min(star_cover_density),
                max(star_cover_density),
            ),
            "BC_range": (min(star_bcs), max(star_bcs)),
            "avg_competitor_D_range": (
                min(average_competitor_D),
                max(average_competitor_D),
            ),
            "ground_heavy_coordinate_coverage": sorted(ground_argmax),
            "ground_field_profiles": sorted(ground_profiles.items()),
        },
    )


def block_transform_audit(A: np.ndarray, name: str, beta: float) -> None:
    """Check the arbitrary-block transform on several non-star blocks."""
    n = len(A)
    _, _, edges, signs, _, _, q, deficit = oriented_states(A)
    weights = np.exp(-beta * deficit)
    D = float(weights.sum())
    edge_index = {edge: j for j, edge in enumerate(edges)}

    blocks = []
    # A triangle, a four-clique, and a short matching.  These test internal
    # blocks as well as disjoint edges without expensive 2^N enumeration.
    for vertices in ((0, 1, 2), (0, 1, 2, 3)):
        blocks.append(
            tuple(
                edge_index[(i, j)]
                for i, j in itertools.combinations(vertices, 2)
            )
        )
    blocks.append(
        tuple(
            edge_index[(2 * k, 2 * k + 1)]
            for k in range(min(n // 2, 4))
        )
    )

    reports = []
    for H in blocks:
        m = len(H)
        z = signs[:, H].sum(axis=1)
        moment = float((weights * np.exp(-2 * beta * z)).sum() / D)
        competitor_D = []
        for bits in itertools.product((0, 1), repeat=m):
            target = np.asarray(bits, dtype=bool)
            selected = signs[:, H][:, target].sum(axis=1)
            ds = float(np.exp(-beta * (deficit + 4 * selected)).sum())
            assert ds >= 1.0 - 1e-12
            competitor_D.append(ds)
        rhs = D * math.cosh(2 * beta) ** m * moment
        assert math.isclose(np.mean(competitor_D), rhs, rel_tol=1e-11)
        assert moment + 1e-12 >= 1.0 / (D * math.cosh(2 * beta) ** m)
        reports.append((m, min(competitor_D), max(competitor_D), moment))
    print(name, f"beta={beta}", "block reports (|H|, minD, maxD, moment)=", reports)


def main() -> None:
    for beta in (0.5, 1.0):
        for A, name in ((A8, "A8"), (A9, "A9")):
            star_audit(A, name, beta)
            block_transform_audit(A, name, beta)
    print("PASS: edge-cube cover, every star replacement, and block transforms")


if __name__ == "__main__":
    main()
