#!/usr/bin/env python3
"""Exact small-order optimized-bridge responses of minimizer orbit classes.

For a hollow signing ``A`` of order ``n`` and a hollow signing ``C`` of
order ``k``, compute

    F_C(A) = min_B Q([[A, B], [B^T, C]])

over all ``n*k``-entry sign bridges.  Orders ``k <= 2`` are independently
exhausted in vectorized chunks.  Order ``k = 3`` is solved by deterministic
CP-SAT feasibility, starting from the certified global lower bound
``M_{n+3}`` and advancing by the exact parity step.

The source class files are exhaustive signed-permutation/global-sign orbit
classifications.  The output compares response vectors only after quotienting
the residual global-sign action correctly: at order three it swaps the two
triangle-product query classes.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "extremal_information"
    / "experiments"
    / "results"
    / "exact_minimizer_optimized_bridge_response.json"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins with coordinate zero fixed to +1."""

    return np.asarray(
        [(1,) + tail for tail in itertools.product((1, -1), repeat=n - 1)],
        dtype=np.int16,
    )


def all_spins(n: int) -> np.ndarray:
    return np.asarray(list(itertools.product((1, -1), repeat=n)), dtype=np.int16)


def energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    n = len(spin)
    return int(
        sum(
            int(matrix[i, j]) * int(spin[i]) * int(spin[j])
            for i in range(n)
            for j in range(i + 1, n)
        )
    )


def cap(matrix: np.ndarray) -> int:
    return max(abs(energy(matrix, spin)) for spin in projective_spins(len(matrix)))


def parent_matrix(A: np.ndarray, C: np.ndarray, B: np.ndarray) -> np.ndarray:
    n, k = B.shape
    parent = np.zeros((n + k, n + k), dtype=np.int8)
    parent[:n, :n] = A
    parent[n:, n:] = C
    parent[:n, n:] = B
    parent[n:, :n] = B.T
    return parent


def contexts(A: np.ndarray, C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return constants and bridge coefficients for every projective context."""

    n, k = len(A), len(C)
    rows: list[list[int]] = []
    constants: list[int] = []
    for x in projective_spins(n):
        h_a = energy(A, x)
        for y in all_spins(k):
            constants.append(h_a + energy(C, y))
            rows.append([int(x[i]) * int(y[j]) for i in range(n) for j in range(k)])
    return np.asarray(constants, dtype=np.int16), np.asarray(rows, dtype=np.int16)


def exhaustive_bridge_minimum(
    A: np.ndarray, C: np.ndarray, chunk_size: int
) -> dict[str, object]:
    """Exhaust every bridge; intended for k <= 2."""

    constants, coefficients = contexts(A, C)
    variables = coefficients.shape[1]
    total = 1 << variables
    bit_positions = np.arange(variables, dtype=np.uint64)
    best = (len(A) + len(C)) * (len(A) + len(C) - 1) // 2
    best_count = 0
    witness: np.ndarray | None = None
    started = time.monotonic()
    for start in range(0, total, chunk_size):
        masks = np.arange(start, min(start + chunk_size, total), dtype=np.uint64)
        bridges = 1 - 2 * ((masks[:, None] >> bit_positions) & 1).astype(np.int16)
        values = bridges @ coefficients.T + constants[None, :]
        caps = np.max(np.abs(values), axis=1)
        chunk_best = int(caps.min())
        if chunk_best < best:
            best = chunk_best
            best_count = int(np.count_nonzero(caps == best))
            witness = bridges[int(np.flatnonzero(caps == best)[0])].copy()
        elif chunk_best == best:
            best_count += int(np.count_nonzero(caps == best))
    assert witness is not None
    B = witness.reshape(len(A), len(C)).astype(np.int8)
    assert cap(parent_matrix(A, C, B)) == best
    return {
        "method": "complete vectorized enumeration",
        "classification": "proved exhaustive finite computation",
        "bridge_count": total,
        "optimal_bridge_count": best_count,
        "objective": best,
        "witness_bridge": B.tolist(),
        "elapsed_seconds": time.monotonic() - started,
    }


def signed_permutation_automorphisms(C: np.ndarray) -> list[np.ndarray]:
    """Signed permutation matrices U satisfying U^T C U = C."""

    k = len(C)
    answer: list[np.ndarray] = []
    for permutation in itertools.permutations(range(k)):
        P = np.eye(k, dtype=np.int8)[:, permutation]
        for signs in itertools.product((1, -1), repeat=k):
            U = P @ np.diag(np.asarray(signs, dtype=np.int8))
            if np.array_equal(U.T @ C @ U, C):
                answer.append(U)
    return answer


def first_row_orbit_representatives(C: np.ndarray) -> list[tuple[int, ...]]:
    """Canonical first rows modulo query automorphisms and B -> -B."""

    k = len(C)
    automorphisms = signed_permutation_automorphisms(C)
    unseen = set(itertools.product((1, -1), repeat=k))
    representatives: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        orbit: set[tuple[int, ...]] = set()
        row = np.asarray(seed, dtype=np.int8)
        for U in automorphisms:
            image = tuple(int(value) for value in row @ U)
            orbit.add(image)
            orbit.add(tuple(-value for value in image))
        representatives.append(max(orbit))
        unseen -= orbit
    # CP variables encode b=1-2z.
    return sorted(
        tuple((1 - value) // 2 for value in representative)
        for representative in representatives
    )


def cp_sat_decision(
    A: np.ndarray, C: np.ndarray, target: int, time_limit: float
) -> tuple[str, dict[str, object], np.ndarray | None]:
    """Decide whether some bridge has cap at most target."""

    n, k = len(A), len(C)
    constants, coefficients = contexts(A, C)
    model = cp_model.CpModel()
    z = [model.NewBoolVar(f"z_{i}_{j}") for i in range(n) for j in range(k)]

    allowed_first_rows = first_row_orbit_representatives(C)
    model.AddAllowedAssignments(z[:k], allowed_first_rows)
    for constant, coefficient in zip(constants.tolist(), coefficients.tolist()):
        # b=1-2z, hence constant + coefficient.b.
        expression = int(constant + sum(coefficient)) - 2 * sum(
            int(value) * variable for value, variable in zip(coefficient, z)
        )
        model.Add(expression <= target)
        model.Add(expression >= -target)

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.max_time_in_seconds = time_limit
    started = time.monotonic()
    status_code = solver.Solve(model)
    elapsed = time.monotonic() - started
    status = solver.StatusName(status_code)
    witness = None
    if status_code in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        witness = np.asarray(
            [1 - 2 * int(solver.Value(variable)) for variable in z], dtype=np.int8
        ).reshape(n, k)
        assert cap(parent_matrix(A, C, witness)) <= target
    record = {
        "target": target,
        "status": status,
        "elapsed_seconds": elapsed,
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "best_objective_bound": solver.BestObjectiveBound(),
        "wall_time_seconds": solver.WallTime(),
        "first_row_orbit_representatives_bits": [list(row) for row in allowed_first_rows],
        "query_signed_permutation_automorphism_count": len(
            signed_permutation_automorphisms(C)
        ),
        "response_stats": solver.ResponseStats(),
    }
    return status, record, witness


def cpsat_bridge_minimum(
    A: np.ndarray,
    C: np.ndarray,
    global_lower_bound: int,
    time_limit: float,
) -> dict[str, object]:
    """Find the exact minimum by parity-stepped feasibility decisions."""

    total_edges = (len(A) + len(C)) * (len(A) + len(C) - 1) // 2
    target = global_lower_bound
    if (target - total_edges) % 2:
        target += 1
    decisions: list[dict[str, object]] = []
    witness: np.ndarray | None = None
    while target <= total_edges:
        status, record, candidate = cp_sat_decision(A, C, target, time_limit)
        decisions.append(record)
        if status in ("OPTIMAL", "FEASIBLE"):
            witness = candidate
            break
        if status != "INFEASIBLE":
            raise RuntimeError(
                f"CP-SAT did not decide target {target}: {status}; increase --time-limit"
            )
        target += 2
    if witness is None:
        raise AssertionError("no bridge found")
    exact_cap = cap(parent_matrix(A, C, witness))
    if exact_cap != target:
        raise AssertionError((exact_cap, target))
    return {
        "method": "parity-stepped deterministic CP-SAT feasibility",
        "classification": "solver-certified exact finite computation; no standalone proof object",
        "global_parent_lower_bound": global_lower_bound,
        "objective": target,
        "witness_bridge": witness.tolist(),
        "decisions": decisions,
    }


def query_representatives() -> list[tuple[str, np.ndarray]]:
    return [
        ("k1_unique", np.zeros((1, 1), dtype=np.int8)),
        ("k2_unique", np.asarray([[0, 1], [1, 0]], dtype=np.int8)),
        (
            "k3_triangle_product_plus",
            np.asarray([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=np.int8),
        ),
        (
            "k3_triangle_product_minus",
            np.asarray([[0, 1, 1], [1, 0, -1], [1, -1, 0]], dtype=np.int8),
        ),
    ]


def verify_query_orbits() -> dict[str, object]:
    """Exhaustively verify the switching/permutation query classes through k=3."""

    counts: dict[str, int] = {}
    products: dict[str, list[int]] = {}
    for k in (1, 2, 3):
        edges = list(itertools.combinations(range(k), 2))
        matrices = []
        for edge_signs in itertools.product((1, -1), repeat=len(edges)):
            C = np.zeros((k, k), dtype=np.int8)
            for value, (i, j) in zip(edge_signs, edges):
                C[i, j] = C[j, i] = value
            matrices.append(C)
        unseen = {matrix.tobytes(): matrix for matrix in matrices}
        orbits: list[list[np.ndarray]] = []
        while unseen:
            _, seed = next(iter(unseen.items()))
            orbit: dict[bytes, np.ndarray] = {}
            for permutation in itertools.permutations(range(k)):
                P = np.eye(k, dtype=np.int8)[:, permutation]
                permuted = P.T @ seed @ P
                for signs in itertools.product((1, -1), repeat=k):
                    S = np.diag(np.asarray(signs, dtype=np.int8))
                    image = S @ permuted @ S
                    orbit[image.tobytes()] = image
            for key in orbit:
                unseen.pop(key, None)
            orbits.append(list(orbit.values()))
        counts[str(k)] = len(orbits)
        if k == 3:
            products[str(k)] = sorted(
                {
                    int(orbit[0][0, 1] * orbit[0][0, 2] * orbit[0][1, 2])
                    for orbit in orbits
                }
            )
    if counts != {"1": 1, "2": 1, "3": 2} or products["3"] != [-1, 1]:
        raise AssertionError((counts, products))
    return {
        "switching_and_permutation_orbit_counts": counts,
        "k3_orbits_classified_by_triangle_product": products["3"],
    }


def verify_switch_transport(
    A: np.ndarray, C: np.ndarray, B: np.ndarray
) -> dict[str, object]:
    """Verify witness transport over every projective child switch."""

    parent = parent_matrix(A, C, B)
    reference = cap(parent)
    n, k = len(A), len(C)
    checks = 0
    for s in projective_spins(n):
        S = np.diag(s.astype(np.int8))
        switched_A = S @ A @ S
        for t in projective_spins(k):
            T = np.diag(t.astype(np.int8))
            switched_C = T @ C @ T
            switched_B = S @ B @ T
            G = np.block(
                [
                    [S, np.zeros((n, k), dtype=np.int8)],
                    [np.zeros((k, n), dtype=np.int8), T],
                ]
            )
            transported = parent_matrix(switched_A, switched_C, switched_B)
            if not np.array_equal(transported, G @ parent @ G):
                raise AssertionError("switch transport failed")
            checks += 1

    reverse_n = np.eye(n, dtype=np.int8)[:, ::-1]
    reverse_k = np.eye(k, dtype=np.int8)[:, ::-1]
    permuted = parent_matrix(
        reverse_n.T @ A @ reverse_n,
        reverse_k.T @ C @ reverse_k,
        reverse_n.T @ B @ reverse_k,
    )
    reverse_parent = np.block(
        [
            [reverse_n, np.zeros((n, k), dtype=np.int8)],
            [np.zeros((k, n), dtype=np.int8), reverse_k],
        ]
    )
    if not np.array_equal(permuted, reverse_parent.T @ parent @ reverse_parent):
        raise AssertionError("permutation matrix identity failed")
    if cap(permuted) != reference:
        raise AssertionError("permutation transport failed")
    negated = parent_matrix(-A, -C, -B)
    if not np.array_equal(negated, -parent):
        raise AssertionError("global-sign matrix identity failed")
    if cap(negated) != reference:
        raise AssertionError("global-sign transport failed")
    return {
        "all_projective_child_switch_pairs_checked": checks,
        "reverse_permutation_checked": True,
        "simultaneous_global_sign_transport_checked": True,
    }


def load_global_exact_values() -> tuple[dict[int, int], dict[str, object]]:
    values: dict[int, int] = {}
    sources: dict[str, object] = {}
    for n in range(3, 11):
        path = ROOT / "computations" / "results" / f"exact_m{n}.json"
        payload = json.loads(path.read_text())
        values[n] = int(payload["profile"]["M"])
        sources[str(n)] = {"path": str(path.relative_to(ROOT)), "sha256": sha256_file(path)}
    certified_path = ROOT / "computations" / "results" / "certified_m11_m12.json"
    certified = json.loads(certified_path.read_text())
    values[11] = int(certified["values"]["11"])
    sources["11"] = {
        "path": str(certified_path.relative_to(ROOT)),
        "sha256": sha256_file(certified_path),
    }
    expected = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: 12, 10: 13, 11: 17}
    if values != expected:
        raise AssertionError(values)
    return values, sources


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--time-limit", type=float, default=300.0)
    args = parser.parse_args()

    exact_values, exact_sources = load_global_exact_values()
    queries = query_representatives()
    query_orbit_check = verify_query_orbits()
    records: list[dict[str, object]] = []
    source_classes: dict[str, object] = {}

    for n in range(3, 9):
        class_path = ROOT / "computations" / "results" / f"m{n}_minimizer_orbits.json"
        payload = json.loads(class_path.read_text())
        source_classes[str(n)] = {
            "path": str(class_path.relative_to(ROOT)),
            "sha256": sha256_file(class_path),
            "classification": payload["classification"],
            "class_count": payload["signed_permutation_and_global_sign_class_count"],
        }
        for class_row in payload["classes"]:
            A = np.asarray(class_row["representative_matrix"], dtype=np.int8)
            response: dict[str, int] = {}
            query_records: list[dict[str, object]] = []
            for query_name, C in queries:
                k = len(C)
                if k <= 2:
                    computation = exhaustive_bridge_minimum(A, C, args.chunk_size)
                else:
                    computation = cpsat_bridge_minimum(
                        A, C, exact_values[n + k], args.time_limit
                    )
                objective = int(computation["objective"])
                B = np.asarray(computation["witness_bridge"], dtype=np.int8)
                if objective != cap(parent_matrix(A, C, B)):
                    raise AssertionError("saved witness cap mismatch")
                response[query_name] = objective
                query_records.append(
                    {
                        "query": query_name,
                        "query_order": k,
                        "query_matrix": C.tolist(),
                        "objective": objective,
                        "computation": computation,
                        "transport_verification": verify_switch_transport(A, C, B),
                    }
                )

            canonical_signature = [
                response["k1_unique"],
                response["k2_unique"],
                *sorted(
                    [
                        response["k3_triangle_product_plus"],
                        response["k3_triangle_product_minus"],
                    ]
                ),
            ]
            records.append(
                {
                    "order": n,
                    "class": class_row["class"],
                    "canonical_orbit_sha256": class_row["canonical_orbit_sha256"],
                    "representative_matrix_sha256": class_row[
                        "representative_matrix_sha256"
                    ],
                    "source_self_complementary": class_row["self_complementary"],
                    "oriented_response": response,
                    "canonical_response_signature_mod_global_sign": canonical_signature,
                    "queries": query_records,
                }
            )

    partitions: dict[str, object] = {}
    for n in range(3, 9):
        groups: defaultdict[tuple[int, ...], list[int]] = defaultdict(list)
        for record in records:
            if record["order"] == n:
                groups[tuple(record["canonical_response_signature_mod_global_sign"])].append(
                    int(record["class"])
                )
        partitions[str(n)] = [
            {"signature": list(signature), "classes": classes}
            for signature, classes in sorted(groups.items())
        ]

    output = {
        "schema": "exact-minimizer-optimized-bridge-response-v1",
        "classification": (
            "proved exhaustive finite computation for k<=2; solver-certified exact "
            "finite computation without standalone proof object for k=3"
        ),
        "normalization": "Q(A)=max_x |sum_{i<j} a_ij x_i x_j|",
        "definition": "F_C(A)=min_{B in {+-1}^{n by k}} Q([[A,B],[B^T,C]])",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "ortools": __import__("ortools").__version__,
            "cp_sat_num_search_workers": 1,
            "cp_sat_random_seed": 0,
        },
        "source_minimizer_classes": source_classes,
        "source_global_exact_values": exact_sources,
        "global_exact_values": {str(key): value for key, value in exact_values.items()},
        "query_orbit_check": query_orbit_check,
        "global_sign_quotient_rule": (
            "F_C(-A)=F_{-C}(A); k=3 triangle-product coordinates swap, so the "
            "two k=3 values are sorted in the canonical class signature"
        ),
        "records": records,
        "response_partitions_by_order": partitions,
        "frozen_conclusion": {
            "order_7": (
                "The optimized k=1 response separates class 1 from classes 0 and 2; "
                "the tested k<=3 language does not separate classes 0 and 2."
            ),
            "order_8": (
                "The k<=2 responses coincide, while optimized k=3 responses separate "
                "the two exhaustive exact-minimizer classes: class 0 has sorted pair "
                "[17,17] and class 1 has [19,19]."
            ),
            "asymptotic_scope": (
                "These are fixed-order separations only; the experiment provides no "
                "scaling law or asymptotic response gap."
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(partitions, indent=2, sort_keys=True))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
