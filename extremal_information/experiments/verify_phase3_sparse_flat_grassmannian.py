#!/usr/bin/env python3
"""Exhaustive finite audit of the sparse-flat Grassmannian ball identity."""

from __future__ import annotations

import json
from math import prod
from pathlib import Path


def all_binary_subspaces(d: int) -> tuple[frozenset[int], ...]:
    spaces = {frozenset({0})}
    frontier = [frozenset({0})]
    while frontier:
        space = frontier.pop()
        for vector in range(1 << d):
            if vector in space:
                continue
            enlarged = frozenset(space | {x ^ vector for x in space})
            if enlarged not in spaces:
                spaces.add(enlarged)
                frontier.append(enlarged)
    return tuple(spaces)


def dim(space: frozenset[object]) -> int:
    return len(space).bit_length() - 1


def gaussian_binary(n: int, r: int) -> int:
    if r < 0 or r > n:
        return 0
    if r == 0:
        return 1
    return prod((2 ** (n - i) - 1) for i in range(r)) // prod(
        (2 ** (r - i) - 1) for i in range(r)
    )


def distance_to_space(x: int, space: frozenset[int]) -> int:
    return min(bin(x ^ y).count("1") for y in space)


def directed_distance(space: frozenset[int], center: frozenset[int]) -> int:
    return max(distance_to_space(x, center) for x in space)


def hausdorff(left: frozenset[int], right: frozenset[int]) -> int:
    return max(directed_distance(left, right), directed_distance(right, left))


def quotient_coset(x: int, center: frozenset[int]) -> frozenset[int]:
    return frozenset(x ^ c for c in center)


def sparse_flat_spectrum(
    d: int,
    center: frozenset[int],
    delta: int,
    all_spaces: tuple[frozenset[int], ...],
) -> dict[int, int]:
    cosets = {quotient_coset(x, center) for x in range(1 << d)}
    leader_ball = {
        coset for coset in cosets if min(bin(x).count("1") for x in coset) <= delta
    }
    quotient_subspaces = {
        frozenset(quotient_coset(x, center) for x in space)
        for space in all_spaces
    }
    spectrum: dict[int, int] = {}
    for quotient_space in quotient_subspaces:
        if quotient_space <= leader_ball:
            degree = dim(quotient_space)
            spectrum[degree] = spectrum.get(degree, 0) + 1
    return spectrum


def audit(d: int, k: int) -> dict[str, object]:
    all_spaces = all_binary_subspaces(d)
    grassmannian = tuple(space for space in all_spaces if dim(space) == k)
    checked = 0
    minimum_radius_one_ball = len(grassmannian)
    for center in grassmannian:
        minimum_radius_one_ball = min(
            minimum_radius_one_ball,
            sum(hausdorff(space, center) <= 1 for space in grassmannian),
        )
        for delta in range(d + 1):
            spectrum = sparse_flat_spectrum(d, center, delta, all_spaces)
            predicted = sum(
                gaussian_binary(k, ell) * 2 ** (ell * ell) * count
                for ell, count in spectrum.items()
                if ell <= min(k, d - k)
            )
            observed = sum(
                directed_distance(space, center) <= delta
                for space in grassmannian
            )
            assert predicted == observed, (d, k, center, delta, predicted, observed)
            checked += 1
    if d - k >= 1:
        assert minimum_radius_one_ball >= 2**k
    return {
        "ambient_dimension": d,
        "subspaces_total": len(all_spaces),
        "grassmannian_size": len(grassmannian),
        "center_threshold_pairs_checked": checked,
        "minimum_symmetric_radius_one_ball": minimum_radius_one_ball,
        "systematic_deformation_lower_bound": 2**k if d - k >= 1 else 1,
    }


def rooted_lift_counterexample() -> dict[str, object]:
    d, k, delta = 4, 1, 1
    all_spaces = all_binary_subspaces(d)
    grassmannian = tuple(space for space in all_spaces if dim(space) == k)
    weight_two = frozenset({0, 0b0011})
    weight_one = frozenset({0, 0b0010})

    spectra = [
        sparse_flat_spectrum(d, center, delta, all_spaces)
        for center in (weight_two, weight_one)
    ]
    quotient_weight_lists = []
    symmetric_ball_sizes = []
    for center in (weight_two, weight_one):
        cosets = {quotient_coset(x, center) for x in range(1 << d)}
        quotient_weight_lists.append(
            sorted(min(bin(x).count("1") for x in coset) for coset in cosets)
        )
        symmetric_ball_sizes.append(
            sum(hausdorff(space, center) <= delta for space in grassmannian)
        )
    assert spectra[0] == spectra[1] == {0: 1, 1: 3}
    assert quotient_weight_lists[0] == quotient_weight_lists[1]
    assert symmetric_ball_sizes == [5, 7]

    # Direct sums preserve the quotient Hamming metric, while the maximum
    # weight inside the two kernels is respectively 2r and r.
    direct_sum_diameters = []
    for blocks in range(1, 9):
        direct_sum_diameters.append(
            {"blocks": blocks, "weight_two": 2 * blocks, "weight_one": blocks}
        )
    return {
        "quotient_leader_weights": quotient_weight_lists[0],
        "common_sparse_flat_spectrum_at_one": spectra[0],
        "symmetric_radius_one_ball_sizes": symmetric_ball_sizes,
        "direct_sum_kernel_diameters": direct_sum_diameters,
    }


def main() -> None:
    audits = [audit(4, 2), audit(5, 2)]
    result = {
        "claim": "SF.4 exact on every audited center and integer threshold",
        "audits": audits,
        "rooted_lift_counterexample": rooted_lift_counterexample(),
    }
    output = Path(__file__).with_name("phase3_sparse_flat_grassmannian.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
