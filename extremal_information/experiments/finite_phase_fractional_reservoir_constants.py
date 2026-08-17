#!/usr/bin/env python3
"""Exact basic-support constants for the finite-phase reservoir lemma.

For each 1 <= s <= max_k, enumerate invertible s-by-s sign matrices up to
row permutation.  A basic vertex has mass

    q = B^{-1} 1 > 0,   U(B) = 1^T q.

The largest U(B) is a universal upper bound for every capped instance.
When the maximizing B also has B^{-T}1 >= 0, the capacity q is an exact
lower example by LP duality, so the bound is sharp.

The implementation uses integer permutation expansions for every
determinant; no floating-point decision enters the certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction


def permutation_terms(n: int) -> list[tuple[tuple[int, ...], int]]:
    terms: list[tuple[tuple[int, ...], int]] = []
    for perm in itertools.permutations(range(n)):
        inversions = sum(
            perm[i] > perm[j] for i in range(n) for j in range(i + 1, n)
        )
        terms.append((perm, -1 if inversions % 2 else 1))
    return terms


def determinant(
    rows: tuple[tuple[int, ...], ...],
    terms: list[tuple[tuple[int, ...], int]],
) -> int:
    total = 0
    for perm, sign in terms:
        product = sign
        for i, j in enumerate(perm):
            product *= rows[i][j]
        total += product
    return total


def replace_column(
    rows: tuple[tuple[int, ...], ...], column: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(1 if j == column else row[j] for j in range(len(row)))
        for row in rows
    )


def transpose(rows: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(rows[i][j] for i in range(len(rows))) for j in range(len(rows)))


def solve_ones(
    rows: tuple[tuple[int, ...], ...],
    terms: list[tuple[tuple[int, ...], int]],
) -> tuple[Fraction, ...] | None:
    denominator = determinant(rows, terms)
    if denominator == 0:
        return None
    return tuple(
        Fraction(determinant(replace_column(rows, j), terms), denominator)
        for j in range(len(rows))
    )


def enumerate_order(s: int) -> dict[str, object]:
    sign_rows = tuple(itertools.product((-1, 1), repeat=s))
    terms = permutation_terms(s)
    best_mass = Fraction(-1)
    best_matrix: tuple[tuple[int, ...], ...] | None = None
    best_primal: tuple[Fraction, ...] | None = None
    best_dual: tuple[Fraction, ...] | None = None
    invertible = 0
    positive = 0

    # An invertible matrix has distinct rows, and row permutation leaves
    # B^{-1}1 unchanged.  Hence row subsets are exhaustive.
    for row_indices in itertools.combinations(range(2**s), s):
        matrix = tuple(sign_rows[i] for i in row_indices)
        primal = solve_ones(matrix, terms)
        if primal is None:
            continue
        invertible += 1
        if not all(value > 0 for value in primal):
            continue
        positive += 1
        mass = sum(primal, Fraction())
        if mass > best_mass:
            best_mass = mass
            best_matrix = matrix
            best_primal = primal
            dual = solve_ones(transpose(matrix), terms)
            best_dual = dual

    return {
        "order": s,
        "invertible_row_subsets": invertible,
        "positive_basic_vertices": positive,
        "maximum_vertex_mass": (
            None
            if best_matrix is None
            else [best_mass.numerator, best_mass.denominator]
        ),
        "matrix": best_matrix,
        "primal_B_inverse_one": (
            None
            if best_primal is None
            else [[x.numerator, x.denominator] for x in best_primal]
        ),
        "dual_B_inverse_transpose_one": (
            None
            if best_dual is None
            else [[x.numerator, x.denominator] for x in best_dual]
        ),
        "dual_is_nonnegative": (
            best_dual is not None and all(value >= 0 for value in best_dual)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-k", type=int, default=5)
    parser.add_argument("--output", type=str)
    args = parser.parse_args()
    if not 1 <= args.max_k <= 5:
        raise SystemExit("exact exhaustive mode is intended for 1 <= max-k <= 5")

    records = [enumerate_order(s) for s in range(1, args.max_k + 1)]
    cumulative: list[list[int]] = []
    current = Fraction(0)
    for record in records:
        raw = record["maximum_vertex_mass"]
        if raw is not None:
            current = max(current, Fraction(raw[0], raw[1]))
        cumulative.append([current.numerator, current.denominator])

    result = {
        "max_k": args.max_k,
        "orders": records,
        "cumulative_constants": cumulative,
        "interpretation": (
            "A nonnegative dual at each cumulative maximizer makes the "
            "corresponding capacity vector a matching exact lower example."
        ),
    }
    rendered = json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
