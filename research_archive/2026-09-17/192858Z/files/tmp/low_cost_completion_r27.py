#!/usr/bin/env python3
"""Exact finite audit for the Wave 27 low-cost completion memo.

The enumerations use one representative of every projective child spin and
the full (not projectively quotiented) outside cube.  All row-square,
variance, determinant, and least-squares calculations are exact integers or
Fractions.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import sys

import numpy as np

sys.path.insert(0, "tmp")
from check_response_dual_r16 import A8, A9, qnorm, spins  # noqa: E402
from verify_compatible_replacement_r12 import A6  # noqa: E402


def bareiss_det(a: np.ndarray) -> int:
    """Exact determinant by fraction-free elimination."""
    b = [[int(value) for value in row] for row in a]
    n = len(b)
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        if b[pivot_index][pivot_index] == 0:
            swap = next(
                (j for j in range(pivot_index + 1, n)
                 if b[j][pivot_index] != 0),
                None,
            )
            if swap is None:
                return 0
            b[pivot_index], b[swap] = b[swap], b[pivot_index]
            sign *= -1
        pivot = b[pivot_index][pivot_index]
        for i in range(pivot_index + 1, n):
            for j in range(pivot_index + 1, n):
                numerator = b[i][j] * pivot - b[i][pivot_index] * b[pivot_index][j]
                assert numerator % previous == 0
                b[i][j] = numerator // previous
        previous = pivot
    return sign * b[-1][-1]


def one_deletion_histogram(a: np.ndarray) -> tuple[Counter, Counter]:
    """Return exact (mean, variance, min, max) child-ground fiber types."""
    a = a.astype(np.int64)
    n = len(a)
    m = n - 1
    parent_q = qnorm(a)
    global_row_ceiling = 2 * (n - 1) * parent_q
    histogram: Counter = Counter()
    distance_histogram: Counter = Counter()

    for selected_tuple in combinations(range(n), m):
        selected = np.array(selected_tuple, dtype=int)
        outside = np.array(
            [i for i in range(n) if i not in selected_tuple], dtype=int
        )
        child = a[np.ix_(selected, selected)]
        child_spins = spins(m)
        child_energies = [int(y @ child @ y) for y in child_spins]
        child_q = max(abs(value) for value in child_energies)

        for y, child_energy in zip(child_spins, child_energies):
            if abs(child_energy) != child_q:
                continue
            orientation = 1 if child_energy > 0 else -1
            b = a[:, selected] @ y
            vmat = a[:, outside]
            h = vmat.T @ b
            gram = vmat.T @ vmat

            rows = []
            parents = []
            for outside_tuple in product((-1, 1), repeat=n - m):
                x = np.empty(n, dtype=np.int64)
                x[selected] = y
                x[outside] = np.array(outside_tuple, dtype=np.int64)
                rows.append(int((a @ x) @ (a @ x)))
                parents.append(orientation * int(x @ a @ x))

            mean_row = Fraction(sum(rows), len(rows))
            variance_row = sum(
                (Fraction(value) - mean_row) ** 2 for value in rows
            ) / len(rows)
            predicted_mean = int(b @ b) + (n - m) * (n - 1)
            predicted_variance = 4 * int(h @ h) + 4 * sum(
                int(gram[i, j]) ** 2
                for i in range(n - m)
                for j in range(i + 1, n - m)
            )
            assert mean_row == predicted_mean
            assert variance_row == predicted_variance

            # Parent-energy mean/variance and row/parent covariance.
            cross = a[np.ix_(outside, selected)] @ y
            mean_parent = Fraction(sum(parents), len(parents))
            variance_parent = sum(
                (Fraction(value) - mean_parent) ** 2 for value in parents
            ) / len(parents)
            covariance = sum(
                (Fraction(row) - mean_row) * (Fraction(parent) - mean_parent)
                for row, parent in zip(rows, parents)
            ) / len(rows)
            predicted_parent_variance = (
                4 * int(cross @ cross) + 2 * (n - m) * (n - m - 1)
            )
            predicted_covariance = 4 * orientation * int(h @ cross) + 4 * orientation * sum(
                int(gram[i, j]) * int(a[outside[i], outside[j]])
                for i in range(n - m)
                for j in range(i + 1, n - m)
            )
            assert mean_parent == child_q
            assert variance_parent == predicted_parent_variance
            assert covariance == predicted_covariance

            # Every exact-child completion has favorable retained loss because
            # its parent oriented payoff is at most Q(A).
            assert max(parents) <= parent_q

            # Conditional expectation, Bhatia--Davis, and the one-sided
            # hypercontractive extraction with the safe constant 1/10.
            minimum = min(rows)
            maximum = max(rows)
            assert minimum <= mean_row
            if mean_row < global_row_ceiling:
                assert (
                    (mean_row - minimum) * (global_row_ceiling - mean_row)
                    >= variance_row
                )
            assert float(minimum) <= float(mean_row) - float(variance_row) ** 0.5 / 10 + 1e-12

            # With one outside bit, greedy conditional expectation is exact.
            assert minimum == mean_row - 2 * abs(int(h[0]))

            # Exact biased least-squares rounding and zonotope distance.
            # Here V has one column of squared norm n-1.
            bias = Fraction(-int(h[0]), n - 1)
            assert -1 <= bias <= 1
            distance_sq = Fraction(int(b @ b)) - Fraction(int(h[0]) ** 2, n - 1)
            biased_mean = distance_sq + (n - 1) * (1 - bias * bias)
            enumerated_biased_mean = (
                Fraction(1 + bias, 2) * rows[1]
                + Fraction(1 - bias, 2) * rows[0]
            )
            assert biased_mean == enumerated_biased_mean

            histogram[(int(mean_row), int(variance_row), minimum, maximum)] += 1
            distance_histogram[distance_sq] += 1

    return histogram, distance_histogram


def main() -> None:
    matrices = tuple(np.asarray(a, dtype=np.int64) for a in (A6, A8, A9))
    determinants = tuple(bareiss_det(a) for a in matrices)
    assert determinants == (-125, 729, 808)

    expected = {
        6: Counter({(30, 0, 30, 30): 60}),
        8: Counter({(68, 16, 64, 72): 24}),
        9: Counter({
            (88, 64, 80, 96): 2,
            (96, 0, 96, 96): 6,
            (104, 64, 96, 112): 4,
            (112, 0, 112, 112): 6,
            (112, 256, 96, 128): 6,
            (120, 64, 112, 128): 2,
            (128, 0, 128, 128): 2,
        }),
    }
    expected_distances = {
        6: Counter({Fraction(25): 60}),
        8: Counter({Fraction(423, 7): 24}),
        9: Counter({
            Fraction(78): 2,
            Fraction(88): 6,
            Fraction(94): 4,
            Fraction(96): 6,
            Fraction(104): 6,
            Fraction(110): 2,
            Fraction(120): 2,
        }),
    }
    for matrix in matrices:
        observed, distances = one_deletion_histogram(matrix)
        assert observed == expected[len(matrix)]
        assert distances == expected_distances[len(matrix)]
        print(f"A{len(matrix)}: {dict(sorted(observed.items()))}")
        print(f"  squared zonotope distances: {dict(sorted(distances.items()))}")

    print({"determinants": determinants})
    print("PASS: exact low-cost completion identities and finite obstructions")


if __name__ == "__main__":
    main()
