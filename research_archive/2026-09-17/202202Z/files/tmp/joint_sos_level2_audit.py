#!/usr/bin/env python3
"""Degree-four Boolean moment-SOS audit for selected saved signing seeds."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import cvxpy as cp
import numpy as np


ROOT = Path("/home/math/quadra")


def load(path: str, key: str) -> np.ndarray:
    return np.asarray(json.loads((ROOT / path).read_text())[key], dtype=float)


def masks_upto(n: int, degree: int) -> list[int]:
    result = []
    for size in range(degree + 1):
        for subset in itertools.combinations(range(n), size):
            mask = 0
            for index in subset:
                mask |= 1 << index
            result.append(mask)
    return result


def endpoint(matrix: np.ndarray, maximize: bool) -> tuple[float, str]:
    n = len(matrix)
    basis = masks_upto(n, 2)
    moments = masks_upto(n, 4)
    position = {mask: index for index, mask in enumerate(moments)}
    y = cp.Variable(len(moments))
    moment_matrix = cp.bmat(
        [[y[position[left ^ right]] for right in basis] for left in basis]
    )
    objective = sum(
        matrix[i, j] * y[position[(1 << i) | (1 << j)]]
        for i in range(n)
        for j in range(i + 1, n)
    )
    problem = cp.Problem(
        cp.Maximize(objective) if maximize else cp.Minimize(objective),
        [y[position[0]] == 1, moment_matrix >> 0],
    )
    value = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-8,
        tol_feas=1e-8,
        tol_gap_rel=1e-8,
        max_iter=1000,
    )
    return float(value), str(problem.status)


def exact_cap(matrix: np.ndarray) -> int:
    n = len(matrix)
    best = 0
    integer = matrix.astype(np.int64)
    for code in range(1 << (n - 1)):
        spin = np.ones(n, dtype=np.int64)
        for i in range(1, n):
            if code & (1 << (i - 1)):
                spin[i] = -1
        best = max(best, abs(int(spin @ integer @ spin) // 2))
    return best


def main() -> None:
    m8 = json.loads((ROOT / "computations/results/m8_minimizer_orbits.json").read_text())
    cases = [
        ("M6_conference", load("computations/results/exact_m6.json", "matrix")),
        ("M8_class0", np.asarray(m8["classes"][0]["representative_matrix"], dtype=float)),
        ("M8_class1", np.asarray(m8["classes"][1]["representative_matrix"], dtype=float)),
        ("M10_nonconference", load("computations/results/exact_m10.json", "matrix")),
        ("GF9_conference_n10", load("computations/results/conference_order10_gf9.json", "conference_matrix")),
        ("M12_witness", load("computations/results/heuristic_m12_seed20260731.json", "matrix")),
        ("M14_conference", load("computations/results/conference_completion_m13.json", "conference_matrix")),
    ]
    for name, matrix in cases:
        upper, upper_status = endpoint(matrix, True)
        lower, lower_status = endpoint(matrix, False)
        width = (upper - lower) / 2
        print(json.dumps({
            "name": name,
            "n": len(matrix),
            "cap": exact_cap(matrix),
            "U2": upper,
            "L2": lower,
            "W2": width,
            "R2": max(upper, -lower),
            "width_gap": width - exact_cap(matrix),
            "statuses": [upper_status, lower_status],
        }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
