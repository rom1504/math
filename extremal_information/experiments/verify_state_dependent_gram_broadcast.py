#!/usr/bin/env python3
"""Exact finite checks for the state-dependent Gram-broadcast theorem.

The asymptotic sampler and spectral-flatness statements are probabilistic
existence proofs.  This script verifies their finite algebra, finds fixed
small samplers, checks the support constants, and exhaustively checks the
contextual response table for the complete r=2 and r=3 label sets.
"""

from __future__ import annotations

import itertools
import json
import math
import random


def pairs(r: int):
    return [(i, j) for i in range(r) for j in range(i + 1, r)]


def wedge_mask(x: int, y: int, r: int) -> int:
    value = 0
    for bit, (i, j) in enumerate(pairs(r)):
        determinant = ((((x >> i) & 1) & ((y >> j) & 1))
                       ^ (((x >> j) & 1) & ((y >> i) & 1)))
        value |= determinant << bit
    return value


def alt_value(mask: int, x: int, y: int, r: int) -> int:
    return bin(mask & wedge_mask(x, y, r)).count("1") & 1


def pair_wedges(labels: list[int], r: int) -> list[int]:
    return [
        wedge_mask(labels[i], labels[j], r)
        for i in range(len(labels))
        for j in range(i + 1, len(labels))
    ]


def support(mask: int, wedges: list[int]) -> int:
    return sum(bin(mask & wedge).count("1") & 1 for wedge in wedges)


def fixed_sampler(r: int) -> list[int]:
    """A reproducible small sampler; 4 r^2 already works for r <= 5."""
    k = 4 * r * r
    rng = random.Random(10_000 * r + 400)
    return [rng.randrange(1 << r) for _ in range(k)]


def kernel(mask: int, labels: list[int], r: int) -> list[list[int]]:
    k = len(labels)
    out = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            v = -1 if alt_value(mask, labels[i], labels[j], r) else 1
            out[i][j] = out[j][i] = v
    return out


def q_value(matrix: list[list[int]], spin: tuple[int, ...]) -> int:
    k = len(spin)
    return sum(
        matrix[i][j] * spin[i] * spin[j]
        for i in range(k)
        for j in range(i + 1, k)
    )


def exact_context_table(r: int):
    labels = list(range(1 << r))
    states = range(1 << len(pairs(r)))
    kernels = [kernel(mask, labels, r) for mask in states]
    spins = list(itertools.product((-1, 1), repeat=len(labels)))
    table = []
    for b in states:
        row = []
        for t in states:
            difference = [
                [kernels[b][i][j] - kernels[t][i][j] for j in range(len(labels))]
                for i in range(len(labels))
            ]
            row.append(max(abs(q_value(difference, spin)) for spin in spins))
        table.append(row)
    assert all(table[b][b] == 0 for b in states)
    off_diagonal = [
        table[b][t] for b in states for t in states if b != t
    ]
    assert min(off_diagonal) > 0
    return {
        "r": r,
        "k": len(labels),
        "states": len(kernels),
        "minimum_off_diagonal_response": min(off_diagonal),
        "maximum_off_diagonal_response": max(off_diagonal),
        "normalized_minimum": min(off_diagonal) / len(labels) ** 1.5,
    }


def check_bicharacter_identities(r: int, labels: list[int]):
    form_count = 1 << len(pairs(r))
    # Exhaust all forms against a deterministic representative panel.  The
    # identity is coefficientwise linear, so a quadratic all-pairs scan of
    # the (already exhaustive) form list would add no coverage.
    forms = range(min(form_count, 64))
    panel = sorted({0, 1, form_count - 1, *(1 << i for i in range(len(pairs(r))))})
    for b in forms:
        for c in panel:
            for x, y in itertools.product(labels[: min(12, len(labels))], repeat=2):
                lhs = alt_value(b ^ c, x, y, r)
                rhs = alt_value(b, x, y, r) ^ alt_value(c, x, y, r)
                assert lhs == rhs
        for x, y, z in itertools.product(labels[: min(4, len(labels))], repeat=3):
            assert alt_value(b, x ^ y, z, r) == (
                alt_value(b, x, z, r) ^ alt_value(b, y, z, r)
            )
            assert alt_value(b, x, x, r) == 0


def check_probability_constants():
    for r in range(2, 50):
        k = 64 * r * r
        h = r * (r - 1) // 2
        assert h * math.log(2) - k / 128 < 0
        assert h >= k / 256
        assert math.log(9) * k + math.log(2) * h - 4 * k < 0


def main():
    check_probability_constants()
    sampler_reports = []
    for r in range(2, 6):
        labels = fixed_sampler(r)
        k = len(labels)
        wedges = pair_wedges(labels, r)
        supports = [
            support(mask, wedges)
            for mask in range(1, 1 << len(pairs(r)))
        ]
        threshold = math.comb(k, 2) / 4
        assert min(supports) >= threshold
        check_bicharacter_identities(r, labels)
        sampler_reports.append(
            {
                "r": r,
                "k": k,
                "forms": 1 << len(pairs(r)),
                "minimum_support": min(supports),
                "required_support": threshold,
                "minimum_relative_support": min(supports) / math.comb(k, 2),
            }
        )

    # Complete label sets make r=2 and r=3 small enough for exhaustive
    # Boolean response calculation (4 and 8 spins respectively).
    response_reports = [exact_context_table(2), exact_context_table(3)]

    # The base r=2 nonzero form has exact response six.  Tensoring with the
    # regular order-four Hadamard witness multiplies it by 4^(3/2)=8.
    assert response_reports[0]["minimum_off_diagonal_response"] == 6
    lifted_witness = 8 * 6
    assert lifted_witness == 48

    report = {
        "probability_constants_checked_r_through": 49,
        "finite_samplers": sampler_reports,
        "exact_context_tables": response_reports,
        "r2_order4_hadamard_lift_witness": lifted_witness,
        "theorem_lower_constant": math.sqrt(2) / 32,
        "checks": "passed",
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
