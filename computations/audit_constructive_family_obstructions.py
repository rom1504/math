#!/usr/bin/env python3
"""Reproduce finite checks used in the phase-2 constructive-family audit."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from conference_prime_square import PrimeSquare


def energy(matrix: np.ndarray, spins: np.ndarray) -> int:
    return int(spins.astype(np.int64) @ matrix.astype(np.int64) @ spins // 2)


def exact_cap(matrix: np.ndarray) -> tuple[int, int, int]:
    n = len(matrix)
    spins = np.ones((1 << (n - 1), n), dtype=np.int8)
    codes = np.arange(1 << (n - 1), dtype=np.uint64)
    for j in range(1, n):
        spins[:, j] = 1 - 2 * ((codes >> (j - 1)) & 1).astype(np.int8)
    values = np.einsum(
        "bi,ij,bj->b", spins, matrix.astype(np.int64), spins, optimize=True
    ) // 2
    return int(np.abs(values).max()), int(values.min()), int(values.max())


def characteristic_coefficients_from_traces(traces: list[int]) -> list[int]:
    elementary = [Fraction(1)]
    for k in range(1, len(traces) + 1):
        value = sum(
            (-1) ** (i - 1) * elementary[k - i] * traces[i - 1]
            for i in range(1, k + 1)
        ) / k
        elementary.append(value)
    if any(value.denominator != 1 for value in elementary):
        raise AssertionError(elementary)
    return [int((-1) ** k * elementary[k]) for k in range(len(elementary))]


def matrix_from_root_mask(n: int, mask: int) -> np.ndarray:
    matrix = np.ones((n, n), dtype=np.int8)
    np.fill_diagonal(matrix, 0)
    edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    for index, (i, j) in enumerate(edges):
        if mask & (1 << index):
            matrix[i, j] = matrix[j, i] = -1
    return matrix


def audit_cospectral_pair() -> dict[str, object]:
    records = []
    for mask in (6875, 6887):
        matrix = matrix_from_root_mask(8, mask)
        power = np.eye(8, dtype=np.int64)
        traces = []
        for _ in range(8):
            power = power @ matrix.astype(np.int64)
            traces.append(int(np.trace(power)))
        cap, minimum, maximum = exact_cap(matrix)
        records.append(
            {
                "root_internal_edge_mask": mask,
                "cap": cap,
                "min_energy": minimum,
                "max_energy": maximum,
                "power_traces_1_through_8": traces,
                "characteristic_coefficients_descending": (
                    characteristic_coefficients_from_traces(traces)
                ),
                "matrix_sha256": hashlib.sha256(matrix.tobytes()).hexdigest(),
                "matrix": matrix.astype(int).tolist(),
            }
        )
    if records[0]["power_traces_1_through_8"] != records[1][
        "power_traces_1_through_8"
    ]:
        raise AssertionError("pair is not cospectral")
    if records[0]["cap"] == records[1]["cap"]:
        raise AssertionError("pair does not separate cap")
    return {"order": 8, "records": records}


def line_union_checks() -> list[dict[str, object]]:
    records = []
    for r in (3, 5, 7):
        matrix = PrimeSquare(r).conference()
        for labels in itertools.chain.from_iterable(
            itertools.combinations(range(r), s) for s in range(1, r + 1)
        ):
            keep = [1 + a * r + b for a in range(r) for b in labels]
            principal = matrix[np.ix_(keep, keep)]
            s = len(labels)
            best = -10**18
            for signs in itertools.product((-1, 1), repeat=s):
                spins = np.array(
                    [signs[index] for _a in range(r) for index in range(s)],
                    dtype=np.int8,
                )
                best = max(best, energy(principal, spins))
            predicted = r * (s * r - (s % 2)) // 2
            if best != predicted:
                raise AssertionError((r, labels, best, predicted))
        records.append(
            {
                "r": r,
                "subsets_checked": 2**r - 1,
                "status": "all line-union constant-spin energies match formula",
            }
        )
    return records


def deletion_identity_checks() -> dict[str, object]:
    r = 5
    matrix = PrimeSquare(r).conference().astype(np.int64)
    line_signs = np.ones(r, dtype=np.int64)
    line_signs[(r + 1) // 2 :] = -1
    spins = np.ones(r * r + 1, dtype=np.int64)
    for a in range(r):
        for b in range(r):
            spins[1 + a * r + b] = line_signs[b]
    if not np.array_equal(matrix @ spins, r * spins):
        raise AssertionError("explicit Boolean eigenvector failed")
    subsets = [tuple(range(d)) for d in range(0, 11)]
    for deleted in subsets:
        kept = tuple(i for i in range(len(matrix)) if i not in deleted)
        full = energy(matrix, spins)
        deleted_energy = energy(
            matrix[np.ix_(deleted, deleted)], spins[list(deleted)]
        )
        kept_energy = energy(matrix[np.ix_(kept, kept)], spins[list(kept)])
        if kept_energy != full - r * len(deleted) + deleted_energy:
            raise AssertionError((deleted, kept_energy, full, deleted_energy))
    return {
        "r": r,
        "conference_order": r * r + 1,
        "deletion_prefix_sizes_checked": list(range(0, 11)),
        "status": "H_K=H_full-r|D|+H_D exactly",
    }


def polynomial_fractional_double_checks() -> list[dict[str, object]]:
    repo = Path(__file__).resolve().parents[1]
    structure = json.loads(
        (repo / "computations/results/algebraic_m12_structure.json").read_text()
    )
    source = json.loads((repo / structure["source"]).read_text())
    parent = np.asarray(source.get("matrix", source.get("parent_matrix")), dtype=np.int8)
    left = structure["partition"][0]
    cases = [
        ("order_6_conference_minimizer", parent[np.ix_(left, left)]),
        ("order_8_cospectral_cap_14", matrix_from_root_mask(8, 6875)),
        ("order_8_cospectral_cap_12", matrix_from_root_mask(8, 6887)),
    ]
    records = []
    for name, child in cases:
        n = len(child)
        zero_diagonal_bridge = child.copy()
        double = np.block(
            [[child, zero_diagonal_bridge], [zero_diagonal_bridge.T, -child]]
        )
        child_cap, _, _ = exact_cap(child)
        double_cap, _, _ = exact_cap(double)
        spins = np.asarray(
            [
                [1 if not (code >> i) & 1 else -1 for i in range(n)]
                for code in range(1 << n)
            ],
            dtype=np.int64,
        )
        child_energies = np.einsum(
            "bi,ij,bj->b", spins, child.astype(np.int64), spins, optimize=True
        ) // 2
        diagonal_invariant_witness = None
        for i, x in enumerate(spins):
            for j, y in enumerate(spins):
                if int(x @ y) != 0:
                    continue
                value = int(child_energies[i] - child_energies[j] + x @ child @ y)
                if abs(value) == double_cap:
                    diagonal_invariant_witness = {
                        "x_code": i,
                        "y_code": j,
                        "energy_for_bridge_S_plus_tI_all_t": value,
                        "x_dot_y": 0,
                    }
                    break
            if diagonal_invariant_witness is not None:
                break
        product = double.astype(np.int64) @ double.astype(np.int64)
        expected = np.block(
            [
                [2 * child.astype(np.int64) @ child.astype(np.int64), np.zeros((n, n), dtype=np.int64)],
                [np.zeros((n, n), dtype=np.int64), 2 * child.astype(np.int64) @ child.astype(np.int64)],
            ]
        )
        if not np.array_equal(product, expected):
            raise AssertionError(name)
        records.append(
            {
                "name": name,
                "order": n,
                "child_cap": child_cap,
                "ideal_equal_child_energy_target": float(2 * np.sqrt(2) * child_cap),
                "zero_diagonal_fractional_double_cap": double_cap,
                "double_square_identity_verified": True,
                "diagonal_invariant_extremizer": diagonal_invariant_witness,
            }
        )
    return records


def sylvester(order: int) -> np.ndarray:
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    if len(matrix) != order or not np.array_equal(matrix @ matrix, order * np.eye(order, dtype=np.int64)):
        raise AssertionError(order)
    return matrix


def hadamard_tensor_checks() -> list[dict[str, object]]:
    child = PrimeSquare(3).conference().astype(np.int64)
    n = len(child)
    records = []
    for t in (2, 4):
        hadamard = sylvester(t)
        fractional = np.kron(hadamard, child)
        expected_square = t * (n - 1) * np.eye(t * n, dtype=np.int64)
        if not np.array_equal(fractional @ fractional, expected_square):
            raise AssertionError((t, "tensor square"))
        integral = fractional.copy()
        missing = []
        for a in range(t):
            for b in range(a + 1, t):
                for i in range(n):
                    u, v = a * n + i, b * n + i
                    if integral[u, v] != 0:
                        raise AssertionError((t, a, b, i))
                    integral[u, v] = integral[v, u] = 1
                    missing.append((u, v))
        off_diagonal = integral[~np.eye(t * n, dtype=bool)]
        if not np.all(np.abs(off_diagonal) == 1):
            raise AssertionError((t, "not an integral signing"))
        perturbation = integral - fractional
        perturbation_norm = float(np.linalg.norm(perturbation.astype(float), ord=2))
        if perturbation_norm > t - 1 + 1e-9:
            raise AssertionError((t, perturbation_norm))
        records.append(
            {
                "conference_child_order": n,
                "hadamard_order": t,
                "parent_order": t * n,
                "tensor_square_identity": f"P0^2={t * (n - 1)}I",
                "missing_cross_edges": len(missing),
                "predicted_missing_cross_edges": t * (t - 1) * n // 2,
                "filled_perturbation_operator_norm": perturbation_norm,
                "proved_operator_norm_bound": t - 1,
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {
        "schema": "quadratic-signing-constructive-family-audit-v1",
        "classification": "exact finite arithmetic checks",
        "deletion_identity": deletion_identity_checks(),
        "line_unions": line_union_checks(),
        "cospectral_cap_collision": audit_cospectral_pair(),
        "polynomial_fractional_double": polynomial_fractional_double_checks(),
        "hadamard_tensor_composition": hadamard_tensor_checks(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
