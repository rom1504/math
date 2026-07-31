#!/usr/bin/env python3
"""Exact local-field-pattern decision for a cyclic order-38 conference.

For target energy 113 the conference moment and local-field penalty identities
leave only six possible oriented row-sum histograms.  This script decides each
histogram separately as a Seidel-switching degree-sequence problem and removes
the exact order-19 cyclic symmetry by lexicographic orbit normalization.

A status of INFEASIBLE for every pattern is a certificate that the selected
signing has no energy of absolute value 113.  UNKNOWN patterns remain open.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

from two_fiber_cyclic_conference import circulant, verify


def matrix_from_sequences(a: list[int], c: list[int]) -> np.ndarray:
    verify(3, a, c)
    matrix_a = circulant(a)
    matrix_c = circulant(c)
    return np.block([[matrix_a, matrix_c], [matrix_c.T, -matrix_a]])


def local_field_patterns(target: int) -> list[dict[int, int]]:
    penalty_target = 24 * (114 - target)
    values = [
        value
        for value in range(-37, 38, 2)
        if (value - 5) * (value - 7) <= penalty_target
    ]
    patterns = []

    def recurse(
        position: int,
        remaining_count: int,
        remaining_sum: int,
        remaining_square: int,
        remaining_penalty: int,
        counts: list[int],
    ) -> None:
        if position == len(values):
            if (
                remaining_count == 0
                and remaining_sum == 0
                and remaining_square == 0
                and remaining_penalty == 0
            ):
                patterns.append(
                    {
                        value: count
                        for value, count in zip(values, counts)
                        if count
                    }
                )
            return
        value = values[position]
        penalty = (value - 5) * (value - 7)
        for count in range(remaining_count + 1):
            if count * penalty > remaining_penalty:
                break
            recurse(
                position + 1,
                remaining_count - count,
                remaining_sum - count * value,
                remaining_square - count * value * value,
                remaining_penalty - count * penalty,
                counts + [count],
            )

    recurse(0, 38, 2 * target, 38 * 37, penalty_target, [])
    return patterns


def equality_indicator(
    model: cp_model.CpModel,
    left: cp_model.IntVar,
    right: cp_model.IntVar,
    name: str,
) -> cp_model.IntVar:
    same = model.new_bool_var(name)
    model.add(left == right).only_enforce_if(same)
    model.add(left != right).only_enforce_if(same.negated())
    return same


def add_cyclic_lex_maximum(
    model: cp_model.CpModel, bits: list[cp_model.IntVar]
) -> None:
    """Require the two-fiber word to be lex-maximal under 19 shifts."""

    if len(bits) != 38:
        raise AssertionError(len(bits))
    for shift in range(1, 19):
        rotated = [
            bits[(index + shift) % 19] if index < 19
            else bits[19 + ((index - 19 + shift) % 19)]
            for index in range(38)
        ]
        prefix = model.new_bool_var(f"shift_{shift}_prefix_0")
        model.add(prefix == 1)
        for index, (left, right) in enumerate(zip(bits, rotated)):
            model.add(left >= right).only_enforce_if(prefix)
            same = equality_indicator(
                model, left, right, f"shift_{shift}_same_{index}"
            )
            next_prefix = model.new_bool_var(
                f"shift_{shift}_prefix_{index + 1}"
            )
            model.add_implication(next_prefix, prefix)
            model.add_implication(next_prefix, same)
            model.add_bool_or(
                [next_prefix, prefix.negated(), same.negated()]
            )
            prefix = next_prefix


def decide_pattern(
    matrix: np.ndarray,
    pattern: dict[int, int],
    time_limit: float,
    workers: int,
) -> dict[str, object]:
    model = cp_model.CpModel()
    n = len(matrix)
    spin_bits = [model.new_bool_var(f"spin_{index}") for index in range(n)]
    # Global negation is free; cyclic lex maximum then removes all 19 shifts.
    model.add(spin_bits[0] == 1)
    add_cyclic_lex_maximum(model, spin_bits)

    disagreements: dict[tuple[int, int], cp_model.IntVar] = {}
    for left in range(n):
        for right in range(left + 1, n):
            disagree = model.new_bool_var(f"xor_{left}_{right}")
            model.add(disagree >= spin_bits[left] - spin_bits[right])
            model.add(disagree >= spin_bits[right] - spin_bits[left])
            model.add(disagree <= spin_bits[left] + spin_bits[right])
            model.add(disagree <= 2 - spin_bits[left] - spin_bits[right])
            disagreements[left, right] = disagree

    values = sorted(pattern)
    field_indicators: dict[int, list[cp_model.IntVar]] = {
        value: [] for value in values
    }
    for vertex in range(n):
        field = model.new_int_var(min(values), max(values), f"field_{vertex}")
        row_constant = int(np.sum(matrix[vertex]))
        terms = []
        for other in range(n):
            if other == vertex:
                continue
            edge = (min(vertex, other), max(vertex, other))
            terms.append(-2 * int(matrix[vertex, other]) * disagreements[edge])
        model.add(field == row_constant + sum(terms))
        indicators = []
        for value in values:
            indicator = model.new_bool_var(f"field_{vertex}_is_{value}")
            model.add(field == value).only_enforce_if(indicator)
            model.add(field != value).only_enforce_if(indicator.negated())
            indicators.append(indicator)
            field_indicators[value].append(indicator)
        model.add(sum(indicators) == 1)
    for value, count in pattern.items():
        model.add(sum(field_indicators[value]) == count)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 20260731
    status = solver.solve(model)
    record: dict[str, object] = {
        "status": solver.status_name(status),
        "wall_time_seconds": solver.wall_time,
        "branches": solver.num_branches,
        "conflicts": solver.num_conflicts,
        "pattern": {str(value): count for value, count in sorted(pattern.items())},
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        spin = np.asarray(
            [1 if solver.value(bit) else -1 for bit in spin_bits], dtype=np.int64
        )
        energy = int(spin @ matrix @ spin // 2)
        fields = spin * (matrix @ spin)
        if energy != sum(value * count for value, count in pattern.items()) // 2:
            raise AssertionError("energy check failed")
        if Counter(map(int, fields)) != Counter(pattern):
            raise AssertionError("field histogram check failed")
        record["energy"] = energy
        record["spin_positive_bits_little_endian"] = format(
            sum(1 << index for index, value in enumerate(spin) if value > 0),
            "x",
        )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orbit-source",
        type=Path,
        default=Path("computations/results/two_fiber_cyclic_conference.json"),
    )
    parser.add_argument("--orbit-index", type=int, default=1)
    parser.add_argument("--target", type=int, default=113)
    parser.add_argument("--time-per-pattern", type=float, default=120.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.orbit_source.read_text())
    orbit = source["k3_exhaustive_natural_orbits"]["orbits"][args.orbit_index]
    matrix = matrix_from_sequences(orbit["a"], orbit["c"])
    reversal = [(-index) % 19 for index in range(19)]
    permutation = [19 + reversal[index] for index in range(19)] + reversal
    switching = np.asarray([1] * 19 + [-1] * 19, dtype=np.int64)
    switched = switching[:, None] * matrix * switching[None, :]
    if not np.array_equal(switched[np.ix_(permutation, permutation)], -matrix):
        raise AssertionError("saved -S equivalence formula failed")
    patterns = local_field_patterns(args.target)
    print(f"target={args.target} exact patterns={len(patterns)}", flush=True)
    records = []
    for index, pattern in enumerate(patterns):
        record = decide_pattern(
            matrix, pattern, args.time_per_pattern, args.workers
        )
        records.append(record)
        print(f"pattern {index + 1}/{len(patterns)} {record}", flush=True)
        if record["status"] in ("OPTIMAL", "FEASIBLE"):
            break

    all_infeasible = len(records) == len(patterns) and all(
        record["status"] == "INFEASIBLE" for record in records
    )
    payload = {
        "schema": "quadratic-signing-order38-local-field-pattern-decision-v1",
        "classification": (
            "certified target exclusion" if all_infeasible else
            "partial exact pattern search; UNKNOWN is not exclusion"
        ),
        "orbit_index": args.orbit_index,
        "orbit_clique_counts": orbit["clique_counts"],
        "target_positive_energy": args.target,
        "absolute_value_reduction": (
            "-S is switching/permutation equivalent to S by second-fiber "
            "switching, fiber swap, and cyclic reversal"
        ),
        "exact_local_field_pattern_count": len(patterns),
        "cyclic_symmetry_reduction": (
            "global negation fixed at spin_0=+1; lex maximum removes all 19 "
            "simultaneous cyclic shifts"
        ),
        "records": records,
        "all_patterns_infeasible": all_infeasible,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
