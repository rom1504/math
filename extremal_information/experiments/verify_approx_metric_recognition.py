#!/usr/bin/env python3
"""Exact finite checks for approximate metric recognition and power drift.

The core checks use only the Python standard library.  If NumPy and SciPy
are importable, a second check solves the nearest-pseudometric linear
program and verifies the sharp constants for q <= 10.
"""

from __future__ import annotations

from fractions import Fraction
import json
from typing import Sequence, Union


Number = Union[int, Fraction]
Matrix = list[list[Number]]


def minplus(left: Sequence[Sequence[Number]], right: Sequence[Sequence[Number]]) -> Matrix:
    size = len(left)
    return [
        [min(left[i][k] + right[k][j] for k in range(size)) for j in range(size)]
        for i in range(size)
    ]


def chain_kernel(size: int, slope: int = 5, defect: int = 1) -> Matrix:
    assert size >= 2 and slope > defect > 0
    return [
        [0 if i == j else slope * abs(i - j) - defect for j in range(size)]
        for i in range(size)
    ]


def predicted_chain_power(
    size: int, power: int, slope: int = 5, defect: int = 1
) -> Matrix:
    return [
        [
            0
            if i == j
            else slope * abs(i - j) - defect * min(power, abs(i - j))
            for j in range(size)
        ]
        for i in range(size)
    ]


def triangle_defect(kernel: Sequence[Sequence[Number]]) -> Number:
    size = len(kernel)
    return max(
        kernel[i][j] - kernel[i][k] - kernel[k][j]
        for i in range(size)
        for j in range(size)
        for k in range(size)
    )


def sup_distance(left: Sequence[Sequence[Number]], right: Sequence[Sequence[Number]]) -> Number:
    return max(
        abs(left[i][j] - right[i][j])
        for i in range(len(left))
        for j in range(len(left))
    )


def shape_distance(first: Sequence[Number], second: Sequence[Number]) -> Fraction:
    differences = [a - b for a, b in zip(first, second)]
    return Fraction(max(differences) - min(differences), 2)


def shortest_path_repair(kernel: Matrix, shift: Fraction) -> Matrix:
    size = len(kernel)
    distance: Matrix = [
        [0 if i == j else kernel[i][j] + shift for j in range(size)]
        for i in range(size)
    ]
    for middle in range(size):
        for i in range(size):
            for j in range(size):
                distance[i][j] = min(
                    distance[i][j], distance[i][middle] + distance[middle][j]
                )
    return distance


def assert_pseudometric(distance: Sequence[Sequence[Number]]) -> None:
    size = len(distance)
    for i in range(size):
        assert distance[i][i] == 0
        for j in range(size):
            assert distance[i][j] >= 0
            assert distance[i][j] == distance[j][i]
            for k in range(size):
                assert distance[i][j] <= distance[i][k] + distance[k][j]


def verify_chain_family() -> dict[str, int]:
    defect_checks = 0
    power_checks = 0
    drift_checks = 0
    repair_checks = 0

    for size in range(3, 11):
        kernel = chain_kernel(size)
        assert triangle_defect(kernel) == 1
        assert sup_distance(kernel, minplus(kernel, kernel)) == 1
        defect_checks += size**3 + size**2

        current = kernel
        for power in range(1, size):
            if power > 1:
                current = minplus(current, kernel)
            predicted = predicted_chain_power(size, power)
            assert current == predicted
            assert sup_distance(current, kernel) == power - 1
            assert shape_distance(current[0], kernel[0]) == Fraction(power - 1, 2)
            power_checks += size**2
            drift_checks += size**2 + size

        # For q >= 4 this family attains the sharp universal repair factor
        # c_q=(q-2)/q.  The theorem's constructive shifted-shortest-path
        # repair attains the same value on this example.
        if size >= 4:
            sharp = Fraction(size - 2, size)
            repaired = shortest_path_repair(kernel, sharp)
            assert_pseudometric(repaired)
            assert sup_distance(kernel, repaired) == sharp
            repair_checks += size**3 + size**2

    return {
        "defect_scalar_checks": defect_checks,
        "power_entry_checks": power_checks,
        "projective_drift_checks": drift_checks,
        "constructive_repair_checks": repair_checks,
    }


def negative_edge_kernel(size: int) -> list[list[float]]:
    kernel = [[0.0 for _ in range(size)] for _ in range(size)]
    kernel[0][1] = kernel[1][0] = -0.5
    return kernel


def lp_distance_to_pseudometrics(kernel: Sequence[Sequence[Number]]) -> float:
    """Return the exact LP objective numerically; called only with SciPy."""

    import numpy as np
    from scipy.optimize import linprog

    size = len(kernel)
    pairs = [(i, j) for i in range(size) for j in range(i + 1, size)]
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    variable_count = len(pairs) + 1
    error_index = variable_count - 1

    objective = np.zeros(variable_count)
    objective[error_index] = 1.0
    inequalities: list[object] = []
    bounds: list[float] = []

    for index, (i, j) in enumerate(pairs):
        # d_ij - epsilon <= K_ij and -d_ij - epsilon <= -K_ij.
        row = np.zeros(variable_count)
        row[index] = 1.0
        row[error_index] = -1.0
        inequalities.append(row)
        bounds.append(float(kernel[i][j]))

        row = np.zeros(variable_count)
        row[index] = -1.0
        row[error_index] = -1.0
        inequalities.append(row)
        bounds.append(-float(kernel[i][j]))

    for i, j in pairs:
        for middle in range(size):
            if middle in (i, j):
                continue
            row = np.zeros(variable_count)
            row[pair_index[(i, j)]] = 1.0
            row[pair_index[tuple(sorted((i, middle)))]] -= 1.0
            row[pair_index[tuple(sorted((j, middle)))]] -= 1.0
            inequalities.append(row)
            bounds.append(0.0)

    result = linprog(
        objective,
        A_ub=np.asarray(inequalities),
        b_ub=np.asarray(bounds),
        bounds=[(0.0, None)] * variable_count,
        method="highs",
    )
    assert result.success, result.message
    return float(result.fun)


def verify_lp_sharpness() -> dict[str, object]:
    try:
        import numpy  # noqa: F401
        import scipy  # noqa: F401
    except ImportError:
        return {"status": "skipped (NumPy/SciPy unavailable)", "cases": 0}

    observed: dict[int, float] = {}
    for size in range(2, 11):
        if size <= 3:
            kernel: Sequence[Sequence[Number]] = negative_edge_kernel(size)
            expected = 0.5
        else:
            kernel = chain_kernel(size)
            expected = (size - 2) / size
        optimum = lp_distance_to_pseudometrics(kernel)
        assert abs(optimum - expected) <= 1e-8, (size, optimum, expected)
        observed[size] = optimum

    return {"status": "passed", "cases": len(observed), "optima": observed}


def main() -> None:
    print(
        json.dumps(
            {
                "core": verify_chain_family(),
                "optional_lp": verify_lp_sharpness(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
