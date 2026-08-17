#!/usr/bin/env python3
"""Exhaustive GF(2) checks for the parity-trellis response benchmark.

A parity-check column is stored as an integer bit mask in the syndrome
space.  Every past and future word is enumerated; no randomness or external
package is used.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class CodeCut:
    """A split parity-check matrix H=[H_P H_F], stored by columns."""

    name: str
    syndrome_bits: int
    past_columns: Tuple[int, ...]
    future_columns: Tuple[int, ...]


def syndrome(columns: Sequence[int], word: int) -> int:
    """Return the XOR of columns selected by the binary word."""

    value = 0
    for index, column in enumerate(columns):
        if (word >> index) & 1:
            value ^= column
    return value


def image(columns: Sequence[int]) -> FrozenSet[int]:
    """Enumerate the image of the matrix having the given columns."""

    return frozenset(syndrome(columns, word) for word in range(1 << len(columns)))


def gf2_rank(vectors: Iterable[int]) -> int:
    """Compute the rank of bit-mask vectors by exact XOR elimination."""

    basis: Dict[int, int] = {}
    for vector in vectors:
        reduced = vector
        while reduced:
            pivot = reduced.bit_length() - 1
            if pivot in basis:
                reduced ^= basis[pivot]
            else:
                basis[pivot] = reduced
                break
    return len(basis)


def subspace_dimension(size: int) -> int:
    """Return log2(size), asserting that size is a positive power of two."""

    assert size > 0 and size & (size - 1) == 0
    return size.bit_length() - 1


def verify_cut(cut: CodeCut) -> Dict[str, int]:
    """Exhaust every word and verify all quotient and dimension identities."""

    limit = 1 << cut.syndrome_bits
    assert cut.past_columns and cut.future_columns
    assert all(0 <= column < limit for column in cut.past_columns)
    assert all(0 <= column < limit for column in cut.future_columns)

    past_words = range(1 << len(cut.past_columns))
    future_words = range(1 << len(cut.future_columns))
    past_syndrome = {
        word: syndrome(cut.past_columns, word) for word in past_words
    }
    future_syndrome = {
        word: syndrome(cut.future_columns, word) for word in future_words
    }

    codewords = frozenset(
        (past, future)
        for past in past_words
        for future in future_words
        if past_syndrome[past] == future_syndrome[future]
    )
    compatible_futures = {
        past: frozenset(
            future
            for future in future_words
            if (past, future) in codewords
        )
        for past in past_words
    }
    reachable_pasts = tuple(
        past for past in past_words if compatible_futures[past]
    )

    past_supported = frozenset(
        past for past in past_words if past_syndrome[past] == 0
    )
    future_supported = frozenset(
        future for future in future_words if future_syndrome[future] == 0
    )
    past_coset = {
        past: frozenset(past ^ supported for supported in past_supported)
        for past in reachable_pasts
    }

    pair_checks = 0
    for first in reachable_pasts:
        for second in reachable_pasts:
            same_future_coset = (
                compatible_futures[first] == compatible_futures[second]
            )
            same_partial_syndrome = (
                past_syndrome[first] == past_syndrome[second]
            )
            same_past_quotient_coset = past_coset[first] == past_coset[second]
            assert same_future_coset == ((first ^ second) in past_supported)
            assert same_future_coset == same_partial_syndrome
            assert same_future_coset == same_past_quotient_coset
            pair_checks += 1

    past_image = image(cut.past_columns)
    future_image = image(cut.future_columns)
    intersection = past_image & future_image
    reachable_syndromes = frozenset(
        past_syndrome[past] for past in reachable_pasts
    )
    assert reachable_syndromes == intersection

    # A hard-fixed future word exposes exactly one partial-syndrome state.
    probe_checks = 0
    for state in intersection:
        witness = next(
            future for future in future_words if future_syndrome[future] == state
        )
        accepted_pasts = frozenset(
            past for past in reachable_pasts if (past, witness) in codewords
        )
        expected_pasts = frozenset(
            past for past in reachable_pasts if past_syndrome[past] == state
        )
        assert accepted_pasts == expected_pasts
        assert accepted_pasts
        probe_checks += 1

    rank_past = gf2_rank(cut.past_columns)
    rank_future = gf2_rank(cut.future_columns)
    rank_total = gf2_rank(cut.past_columns + cut.future_columns)
    intersection_dimension = subspace_dimension(len(intersection))
    assert intersection_dimension == rank_past + rank_future - rank_total

    code_dimension = subspace_dimension(len(codewords))
    past_supported_dimension = subspace_dimension(len(past_supported))
    future_supported_dimension = subspace_dimension(len(future_supported))
    quotient_dimension = (
        code_dimension - past_supported_dimension - future_supported_dimension
    )
    assert quotient_dimension == intersection_dimension

    distinct_future_cosets = frozenset(
        compatible_futures[past] for past in reachable_pasts
    )
    distinct_past_cosets = frozenset(
        past_coset[past] for past in reachable_pasts
    )
    state_count = len(intersection)
    assert state_count == 1 << quotient_dimension
    assert state_count == len(distinct_future_cosets)
    assert state_count == len(distinct_past_cosets)
    assert len(reachable_pasts) // len(past_supported) == state_count
    assert (
        len(codewords) // (len(past_supported) * len(future_supported))
        == state_count
    )

    return {
        "ambient_syndromes": limit,
        "code_dimension": code_dimension,
        "codewords_checked": len(codewords),
        "future_supported_dimension": future_supported_dimension,
        "hard_future_probe_checks": probe_checks,
        "intersection_dimension": intersection_dimension,
        "past_future_assignments_checked": len(past_words) * len(future_words),
        "past_supported_dimension": past_supported_dimension,
        "quotient_dimension": quotient_dimension,
        "reachable_past_pair_checks": pair_checks,
        "reachable_pasts": len(reachable_pasts),
        "state_count": state_count,
    }


def exhaustive_small_matrix_census() -> Dict[str, int]:
    """Check every split matrix of the declared one- and two-bit sizes."""

    totals = {
        "codewords_checked": 0,
        "hard_future_probe_checks": 0,
        "matrix_cases": 0,
        "past_future_assignments_checked": 0,
        "reachable_past_pair_checks": 0,
    }
    for syndrome_bits in (1, 2):
        columns = range(1 << syndrome_bits)
        for past_width in (1, 2):
            for future_width in (1, 2):
                for past_columns in itertools.product(columns, repeat=past_width):
                    for future_columns in itertools.product(
                        columns, repeat=future_width
                    ):
                        stats = verify_cut(
                            CodeCut(
                                name="exhaustive-census-case",
                                syndrome_bits=syndrome_bits,
                                past_columns=past_columns,
                                future_columns=future_columns,
                            )
                        )
                        totals["matrix_cases"] += 1
                        for key in (
                            "codewords_checked",
                            "hard_future_probe_checks",
                            "past_future_assignments_checked",
                            "reachable_past_pair_checks",
                        ):
                            totals[key] += stats[key]

    # 36 one-check-bit matrices and 400 two-check-bit matrices.
    assert totals["matrix_cases"] == 436
    return totals


def named_examples() -> List[Dict[str, object]]:
    """Audit examples exhibiting full, compressed, and trivial interfaces."""

    cuts = (
        CodeCut(
            name="full_two_bit_identity",
            syndrome_bits=2,
            past_columns=(0b01, 0b10),
            future_columns=(0b01, 0b10),
        ),
        CodeCut(
            name="four_check_bits_compress_to_two",
            syndrome_bits=4,
            past_columns=(0b0001, 0b0010, 0b0100),
            future_columns=(0b0001, 0b0010, 0b1000),
        ),
        CodeCut(
            name="disjoint_images_leave_only_zero_state",
            syndrome_bits=4,
            past_columns=(0b0001, 0b0010),
            future_columns=(0b0100, 0b1000),
        ),
        CodeCut(
            name="nontrivial_past_and_future_supported_subcodes",
            syndrome_bits=3,
            past_columns=(0b001, 0b010, 0b011),
            future_columns=(0b010, 0b100, 0b110),
        ),
    )

    reports: List[Dict[str, object]] = []
    for cut in cuts:
        report: Dict[str, object] = {
            "name": cut.name,
            "past_columns": list(cut.past_columns),
            "future_columns": list(cut.future_columns),
            "syndrome_bits": cut.syndrome_bits,
        }
        report.update(verify_cut(cut))
        reports.append(report)

    assert reports[1]["ambient_syndromes"] == 16
    assert reports[1]["state_count"] == 4
    assert reports[2]["state_count"] == 1
    assert reports[3]["past_supported_dimension"] == 1
    assert reports[3]["future_supported_dimension"] == 1
    assert reports[3]["state_count"] == 2
    return reports


def main() -> None:
    print(
        json.dumps(
            {
                "exhaustive_small_matrix_census": exhaustive_small_matrix_census(),
                "named_examples": named_examples(),
                "status": "all parity-trellis quotient checks passed",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
