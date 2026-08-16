#!/usr/bin/env python3
"""Exact finite falsifiers for the phase-2 response-geometry results.

The analytic proofs live in ``../theorems.md`` after promotion.  This script
checks only finite identities:

* the inverse-Hamming modulus of binary endpoint response kernels;
* the equal outer/different inner spectra of the four-bit code pair;
* the strict loss incurred by retaining only covering radius; and
* the rare matching-fibre conditional variance and uniform response gap in
  the deterministic-synchronization counterexample.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path


def hamming(a: int, b: int) -> int:
    return bin(a ^ b).count("1")


def code_profiles(code: tuple[int, ...], dimension: int) -> dict:
    roots = range(1 << dimension)
    outer = Counter(min(hamming(x, c) for c in code) for x in roots)
    inner = Counter(hamming(a, b) for a in code for b in code)
    return {
        "outer": dict(sorted(outer.items())),
        "inner": dict(sorted(inner.items())),
        "radius": max(outer),
    }


def boundary_kernel_modulus(q: int) -> dict:
    """Exhaustively compute kappa for all binary q-by-q kernels.

    The endpoint-pinning response is the flattened kernel itself, with the
    uniform measure on the q^2 endpoint pairs.
    """

    size = q * q
    kernels = list(product((0, 1), repeat=size))
    best = None
    witnesses = 0
    for left, right in combinations(kernels, 2):
        changed = sum(a != b for a, b in zip(left, right))
        squared_l2 = Fraction(changed, size)
        ratio = squared_l2 / changed
        if best is None or ratio < best:
            best = ratio
            witnesses = 1
        elif ratio == best:
            witnesses += 1
    assert best == Fraction(1, size)
    return {
        "q": q,
        "latent_bits": size,
        "kernels_checked": len(kernels),
        "pairwise_modulus": str(best),
        "minimizing_pairs": witnesses,
    }


def rooted_code_modulus(m: int) -> dict:
    """Exhaustively compute kappa for anchored arbitrary codes in Q_m."""

    root_count = 1 << m
    tables = []
    for mask in range(1 << (root_count - 1)):
        code = (0,) + tuple(
            u for u in range(1, root_count) if (mask >> (u - 1)) & 1
        )
        distances = tuple(
            min(hamming(u, c) for c in code) for u in range(root_count)
        )
        tables.append((mask, distances))
    best = None
    for (left_mask, left), (right_mask, right) in combinations(tables, 2):
        changed = hamming(left_mask, right_mask)
        squared_l2 = Fraction(
            sum((a - b) ** 2 for a, b in zip(left, right)), root_count
        )
        ratio = squared_l2 / changed
        best = ratio if best is None else min(best, ratio)
    assert best == Fraction(1, root_count)
    return {
        "dimension": m,
        "anchored_codes_checked": len(tables),
        "inverse_hamming_modulus": str(best),
    }


def quadratic_response_modulus(n: int) -> dict:
    """Check the exact shifted-Ising inverse-Hamming modulus for a=1."""

    edges = list(combinations(range(n), 2))
    spins = list(product((-1, 1), repeat=n))
    responses = []
    for coefficients in product((-1, 1), repeat=len(edges)):
        energy = tuple(
            sum(
                coefficients[k] * x[i] * x[j]
                for k, (i, j) in enumerate(edges)
            )
            for x in spins
        )
        cap = max(energy)
        responses.append((coefficients, tuple(value - cap for value in energy)))
    best = None
    for (left_coeff, left), (right_coeff, right) in combinations(responses, 2):
        changed = sum(a != b for a, b in zip(left_coeff, right_coeff))
        squared_l2 = Fraction(
            sum((a - b) ** 2 for a, b in zip(left, right)), len(spins)
        )
        ratio = squared_l2 / changed
        best = ratio if best is None else min(best, ratio)
    assert best == 4
    return {
        "order": n,
        "signings_checked": len(responses),
        "inverse_hamming_modulus_for_a_1": str(best),
    }


def maxcut_response_modulus(n: int) -> dict:
    """Check the uncentered counterfactual Max-Cut response modulus."""

    edges = list(combinations(range(n), 2))
    spins = list(product((-1, 1), repeat=n))
    responses = []
    for graph in product((0, 1), repeat=len(edges)):
        cut_values = tuple(
            sum(
                graph[k] * (1 - x[i] * x[j]) // 2
                for k, (i, j) in enumerate(edges)
            )
            for x in spins
        )
        responses.append((graph, cut_values))
    best = None
    for (left_graph, left), (right_graph, right) in combinations(responses, 2):
        changed = sum(a != b for a, b in zip(left_graph, right_graph))
        squared_l2 = Fraction(
            sum((a - b) ** 2 for a, b in zip(left, right)), len(spins)
        )
        ratio = squared_l2 / changed
        best = ratio if best is None else min(best, ratio)
    expected = Fraction(1, 2) if n == 2 else Fraction(1, 4)
    assert best == expected
    return {
        "order": n,
        "graphs_checked": len(responses),
        "inverse_hamming_modulus": str(best),
    }


def matching_fibre(m: int) -> dict:
    """Check the rho=1 instance of the rare matching-fibre example."""

    # Conditional on a uniformly chosen matching edge, R_1 is one on one
    # edge and zero on the other m-1 edges.
    mean = Fraction(1, m)
    variance = mean * (1 - mean)
    values = [1] + [0] * (m - 1)
    best_constant_error = min(
        max(abs(Fraction(v) - c) for v in values)
        for c in (Fraction(0), Fraction(1, 2), Fraction(1))
    )
    assert variance == Fraction(m - 1, m * m)
    assert best_constant_error == Fraction(1, 2)
    return {
        "states": 2 * m,
        "matching_edges": m,
        "conditional_variance": str(variance),
        "best_uniform_constant_error": str(best_constant_error),
    }


def tropical_resonance(r: int) -> dict:
    """Check the exact error of replacing 0-diagonal/1-offdiag by all ones."""

    squared_errors = [1 if i == j else 0 for i in range(r) for j in range(r)]
    normalized_mse = Fraction(sum(squared_errors), r * r)
    assert normalized_mse == Fraction(1, r)
    return {
        "order": r,
        "distinguished_diagonal_cells": r,
        "normalized_rank_one_mse": str(normalized_mse),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "extremal_information/experiments/phase2_response_geometry_results.json"
        ),
    )
    args = parser.parse_args()

    code_a = (0b0000, 0b0001, 0b0010, 0b0011)
    code_b = (0b0000, 0b0001, 0b0010, 0b0101)
    profile_a = code_profiles(code_a, 4)
    profile_b = code_profiles(code_b, 4)
    assert profile_a["outer"] == profile_b["outer"]
    assert profile_a["inner"] != profile_b["inner"]

    radius_code_a = code_profiles((0b00, 0b01), 2)
    radius_code_b = code_profiles((0b00, 0b01, 0b10), 2)
    assert radius_code_a["radius"] == radius_code_b["radius"] == 1
    assert radius_code_a["outer"] != radius_code_b["outer"]

    result = {
        "schema": "extremal-information-phase2-response-geometry-v1",
        "boundary_kernel_moduli": [
            boundary_kernel_modulus(q) for q in (1, 2, 3)
        ],
        "rooted_code_moduli": [rooted_code_modulus(m) for m in (1, 2, 3)],
        "shifted_quadratic_moduli": [
            quadratic_response_modulus(n) for n in (2, 3, 4)
        ],
        "maxcut_moduli": [maxcut_response_modulus(n) for n in (2, 3, 4)],
        "outer_spectrum_collision": {
            "code_a": profile_a,
            "code_b": profile_b,
        },
        "radius_is_strictly_coarser": {
            "code_a": radius_code_a,
            "code_b": radius_code_b,
        },
        "rare_matching_fibres": [matching_fibre(m) for m in range(2, 9)],
        "tropical_rank_resonances": [
            tropical_resonance(r) for r in (2, 4, 8, 16, 32)
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
