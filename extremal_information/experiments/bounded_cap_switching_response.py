#!/usr/bin/env python3
"""Exact finite probe of bounded-cap switching-orbit bridge responses.

For a hollow sign matrix ``A`` and a dense sign bridge ``B``, this script
studies the projective response functions

    R_s(y) = max_x { H_{D_s A D_s}(x) + x^T B y }.

All maxima and response distances reported here use exact integer arithmetic.
For orders at most ten the seeded-random-bridge experiment uses the complete
projective switching orbit and every projective Boolean query.  At larger
orders it uses a declared query-linked subset; the resulting distances and
packing sizes are rigorous lower bounds for the full response metric.

This is experimental mathematics, not asymptotic evidence by itself.  The
falsifiable scaling conjecture tested by the program is recorded in the JSON
output and in the accompanying draft.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

CASES = (
    (4, "exact_m4.json", "matrix", "exact minimizer"),
    (6, "exact_m6.json", "matrix", "exact minimizer"),
    (8, "exact_m8.json", "matrix", "exact minimizer"),
    (10, "conference_order10_gf9.json", "conference_matrix", "conference"),
    (12, "heuristic_m12.json", "matrix", "certified-cap heuristic"),
    (14, "conference_double_p13.json", "conference_matrix", "conference"),
)

EPSILONS = (0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25)


def projective_spins(order: int) -> np.ndarray:
    """Return one representative of each pair ``{x,-x}``, with x_0=1."""

    masks = np.arange(1 << (order - 1), dtype=np.uint32)[:, None]
    bits = (masks >> np.arange(order - 1, dtype=np.uint32)) & 1
    tail = (1 - 2 * bits).astype(np.int16)
    return np.concatenate(
        [np.ones((len(masks), 1), dtype=np.int16), tail], axis=1
    )


def load_case(order: int) -> tuple[np.ndarray, dict]:
    matches = [case for case in CASES if case[0] == order]
    if not matches:
        raise ValueError(f"no configured matrix at order {order}")
    _, filename, key, classification = matches[0]
    path = ROOT / "computations" / "results" / filename
    payload = json.loads(path.read_text(encoding="utf-8"))
    matrix = np.asarray(payload[key], dtype=np.int16)
    if matrix.shape != (order, order):
        raise ValueError((path, matrix.shape))
    if not np.array_equal(matrix, matrix.T):
        raise ValueError(f"{path} is not symmetric")
    if np.any(np.diag(matrix)) or not np.all(
        np.isin(matrix[~np.eye(order, dtype=bool)], (-1, 1))
    ):
        raise ValueError(f"{path} is not a hollow sign matrix")
    return matrix, {
        "source": str(path.relative_to(ROOT)),
        "source_classification": classification,
    }


def quadratic_energies(matrix: np.ndarray, states: np.ndarray) -> np.ndarray:
    wide_states = states.astype(np.int32)
    wide_matrix = matrix.astype(np.int32)
    return (
        np.einsum(
            "bi,ij,bj->b", wide_states, wide_matrix, wide_states, optimize=True
        )
        // 2
    ).astype(np.int32)


def canonical(vector: np.ndarray) -> tuple[int, ...]:
    normalized = vector if vector[0] == 1 else -vector
    return tuple(int(value) for value in normalized)


def sha256_matrix(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(matrix, dtype=np.int8).tobytes()).hexdigest()


def sylvester_hadamard(order: int) -> np.ndarray:
    if order <= 0 or order & (order - 1):
        raise ValueError("Sylvester order must be a power of two")
    matrix = np.ones((1, 1), dtype=np.int16)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    assert np.array_equal(matrix @ matrix.T, order * np.eye(order, dtype=int))
    return matrix


def dense_bridges(matrix: np.ndarray, seed: int) -> list[tuple[str, np.ndarray]]:
    order = len(matrix)
    rng = np.random.default_rng(seed + 1009 * order)
    random_bridge = rng.choice(
        np.asarray((-1, 1), dtype=np.int16), size=(order, order)
    )
    self_bridge = matrix + np.eye(order, dtype=np.int16)
    bridges = [("seeded Rademacher", random_bridge), ("A plus identity", self_bridge)]
    if order & (order - 1) == 0:
        bridges.append(("Sylvester Hadamard", sylvester_hadamard(order)))
    return bridges


def linked_family(
    states: np.ndarray,
    fields: np.ndarray,
    top_state: np.ndarray,
    target_count: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Select queries and distinct switches linked by ``s=u* sign(By)``."""

    rng = np.random.default_rng(seed)
    order = rng.permutation(len(states))
    query_rows: list[np.ndarray] = []
    field_rows: list[np.ndarray] = []
    switch_rows: list[np.ndarray] = []
    seen_switches: set[tuple[int, ...]] = set()
    for index in order:
        field = fields[index]
        field_sign = np.where(field >= 0, 1, -1).astype(np.int16)
        switch = top_state * field_sign
        key = canonical(switch)
        if key in seen_switches:
            continue
        seen_switches.add(key)
        query_rows.append(states[index])
        field_rows.append(field)
        switch_rows.append(np.asarray(key, dtype=np.int16))
        if len(query_rows) == target_count:
            break
    return (
        np.asarray(query_rows, dtype=np.int16),
        np.asarray(field_rows, dtype=np.int16),
        np.asarray(switch_rows, dtype=np.int16),
    )


def response_table(
    states: np.ndarray,
    energies: np.ndarray,
    switches: np.ndarray,
    query_fields: np.ndarray,
) -> np.ndarray:
    """Evaluate all displayed responses exactly, exploiting evenness of H."""

    table = np.empty((len(switches), len(query_fields)), dtype=np.int32)
    for row, switch in enumerate(switches):
        switched_fields = query_fields * switch
        linear_scores = states.astype(np.int32) @ switched_fields.astype(np.int32).T
        table[row] = np.max(
            energies[:, None] + np.abs(linear_scores), axis=0
        )
    return table


def projective_distance_matrix(table: np.ndarray) -> np.ndarray:
    count = len(table)
    distances = np.zeros((count, count), dtype=np.float32)
    for first in range(count):
        differences = table[first + 1 :] - table[first]
        if len(differences):
            row = (
                np.max(differences, axis=1) - np.min(differences, axis=1)
            ) / 2.0
            distances[first, first + 1 :] = row
            distances[first + 1 :, first] = row
    return distances


def upper_triangle(distances: np.ndarray) -> np.ndarray:
    return distances[np.triu_indices(len(distances), k=1)]


def linked_deficits(
    table: np.ndarray, energies: np.ndarray, fields: np.ndarray
) -> np.ndarray:
    """Return the directed diagonal-minus-off-diagonal deficit matrix."""

    top_energy = int(np.max(energies))
    predicted_diagonal = top_energy + np.sum(np.abs(fields), axis=1)
    if not np.array_equal(np.diag(table), predicted_diagonal):
        raise AssertionError("query-linked diagonal exposure identity failed")
    deficits = predicted_diagonal[None, :] - table
    if np.min(deficits) < 0:
        raise AssertionError("an off-diagonal response exceeded the roof")
    if np.any(np.diag(deficits)):
        raise AssertionError("linked diagonal deficits must vanish")
    return deficits.astype(np.int32)


def quantiles(values: np.ndarray) -> dict[str, float | None]:
    if not len(values):
        return {name: None for name in ("minimum", "q01", "q10", "median", "q90", "maximum")}
    levels = np.quantile(values.astype(float), (0, 0.01, 0.10, 0.50, 0.90, 1))
    return {
        name: float(round(float(value), 12))
        for name, value in zip(
            ("minimum", "q01", "q10", "median", "q90", "maximum"), levels
        )
    }


def greedy_pack(
    distances: np.ndarray, threshold: float, trials: int, seed: int
) -> int:
    """Return a certified packing lower bound, optimized over seeded orders."""

    rng = np.random.default_rng(seed)
    best = 0
    orders = [np.arange(len(distances), dtype=int)]
    orders.extend(
        rng.permutation(len(distances)) for _ in range(max(0, trials - 1))
    )
    for order in orders:
        chosen: list[int] = []
        for candidate in order:
            if not chosen or np.all(distances[candidate, chosen] >= threshold):
                chosen.append(int(candidate))
        best = max(best, len(chosen))
    return best


def analyze_linked(
    matrix: np.ndarray,
    bridge: np.ndarray,
    states: np.ndarray,
    energies: np.ndarray,
    seed: int,
) -> dict:
    order = len(matrix)
    all_fields = states.astype(np.int32) @ bridge.astype(np.int32).T
    target_count = min(len(states), 1 << (order // 2))
    queries, fields, switches = linked_family(
        states, all_fields, states[int(np.argmax(energies))], target_count, seed
    )
    table = response_table(states, energies, switches, fields)
    distance_matrix = projective_distance_matrix(table)
    distances = upper_triangle(distance_matrix)
    deficit_matrix = linked_deficits(table, energies, fields)
    deficits = deficit_matrix[~np.eye(len(table), dtype=bool)]
    two_way_bounds = np.asarray(
        [
            (int(deficit_matrix[second, first]) + int(deficit_matrix[first, second]))
            / 2
            for first in range(len(table))
            for second in range(first + 1, len(table))
        ],
        dtype=np.float64,
    )
    two_way_matrix = (deficit_matrix + deficit_matrix.T).astype(np.float32) / 2
    if np.any(two_way_bounds > distances + 1e-9):
        raise AssertionError("two-way diagonal exposure exceeded response metric")
    scale = order**1.5
    packing = []
    two_way_packing = []
    for index, epsilon in enumerate(EPSILONS):
        size = greedy_pack(
            distance_matrix,
            epsilon * scale,
            trials=12,
            seed=seed + 7919 * index,
        )
        packing.append(
            {
                "epsilon": epsilon,
                "threshold": epsilon * scale,
                "certified_greedy_packing_lower_bound": size,
                "packing_bits": math.log2(size),
                "packing_bits_per_vertex": math.log2(size) / order,
            }
        )
        diagonal_size = greedy_pack(
            two_way_matrix,
            epsilon * scale,
            trials=12,
            seed=seed + 15485863 * index,
        )
        two_way_packing.append(
            {
                "epsilon": epsilon,
                "certified_packing_from_two_matched_queries": diagonal_size,
                "packing_bits": math.log2(diagonal_size),
                "packing_bits_per_vertex": math.log2(diagonal_size) / order,
            }
        )
    return {
        "family": "query-linked switching orbit",
        "selection_rule": "s_y = u_star coordinatewise-times sign(B y)",
        "target_family_size": target_count,
        "distinct_family_size": len(switches),
        "query_count": len(queries),
        "metric_scope": "displayed linked queries; hence a lower bound for all futures",
        "diagonal_exposure_identity_verified": True,
        "directed_exposure_deficit": quantiles(deficits),
        "two_way_diagonal_exposure_lower_bound": quantiles(two_way_bounds),
        "normalized_two_way_diagonal_exposure_lower_bound": quantiles(
            two_way_bounds / scale
        ),
        "projective_distance": quantiles(distances),
        "normalized_projective_distance": quantiles(distances / scale),
        "packing_lower_bounds": packing,
        "two_matched_query_packing_lower_bounds": two_way_packing,
    }


def analyze_full_orbit(
    bridge: np.ndarray,
    states: np.ndarray,
    energies: np.ndarray,
    seed: int,
) -> dict:
    order = bridge.shape[0]
    fields = states.astype(np.int32) @ bridge.astype(np.int32).T
    table = response_table(states, energies, states, fields)
    distance_matrix = projective_distance_matrix(table)
    distances = upper_triangle(distance_matrix)
    unique_rows = len(np.unique(table, axis=0))
    scale = order**1.5
    packing = []
    for index, epsilon in enumerate(EPSILONS):
        size = greedy_pack(
            distance_matrix,
            epsilon * scale,
            trials=12,
            seed=seed + 104729 * index,
        )
        packing.append(
            {
                "epsilon": epsilon,
                "threshold": epsilon * scale,
                "certified_greedy_packing_lower_bound": size,
                "packing_bits": math.log2(size),
                "packing_bits_per_vertex": math.log2(size) / order,
            }
        )
    return {
        "family": "complete projective switching orbit",
        "family_size": len(states),
        "query_count": len(states),
        "metric_scope": "every projective Boolean query; exact response metric",
        "distinct_response_rows": unique_rows,
        "projective_distance": quantiles(distances),
        "normalized_projective_distance": quantiles(distances / scale),
        "packing_lower_bounds": packing,
    }


def analyze_order(order: int, seed: int) -> dict:
    matrix, provenance = load_case(order)
    states = projective_spins(order)
    energies = quadratic_energies(matrix, states)
    positive_cap = int(np.max(energies))
    negative_cap = int(-np.min(energies))
    result = {
        "order": order,
        **provenance,
        "matrix_sha256_int8": sha256_matrix(matrix),
        "positive_cap": positive_cap,
        "negative_cap": negative_cap,
        "absolute_cap": max(positive_cap, negative_cap),
        "normalized_absolute_cap": max(positive_cap, negative_cap) / order**1.5,
        "projective_spin_count": len(states),
        "bridges": [],
    }
    for bridge_index, (name, bridge) in enumerate(dense_bridges(matrix, seed)):
        bridge_seed = seed + order * 65537 + bridge_index * 8191
        row = {
            "bridge": name,
            "bridge_sha256_int8": sha256_matrix(bridge),
            "operator_norm_over_sqrt_n": float(
                round(
                    float(
                        np.linalg.svd(bridge.astype(float), compute_uv=False)[0]
                        / math.sqrt(order)
                    ),
                    12,
                )
            ),
            "linked_probe": analyze_linked(
                matrix, bridge, states, energies, bridge_seed
            ),
        }
        if order <= 10 and name == "seeded Rademacher":
            row["full_orbit_probe"] = analyze_full_orbit(
                bridge, states, energies, bridge_seed
            )
        result["bridges"].append(row)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orders", type=int, nargs="+", default=[case[0] for case in CASES]
    )
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "extremal_information"
        / "experiments"
        / "results"
        / "bounded_cap_switching_response.json",
    )
    args = parser.parse_args()

    output = {
        "schema": "bounded-cap-switching-response-v1",
        "classification": (
            "exact finite integer maxima and distances; seeded finite experiment; "
            "greedy packing sizes are certified lower bounds, not maxima"
        ),
        "falsifiable_scaling_conjecture": (
            "For some fixed epsilon,c,C>0, a bounded-cap sequence "
            "Q(A_n)<=C n^(3/2) and dense sign bridges B_n admit at least "
            "exp(c n) switching-orbit children with pairwise projective "
            "response distance at least epsilon n^(3/2)."
        ),
        "warning": (
            "Finite packing ratios do not establish the scaling conjecture; "
            "failure of these particular bridges would not disprove it."
        ),
        "seed": args.seed,
        "orders": [analyze_order(order, args.seed) for order in args.orders],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    for row in output["orders"]:
        random_row = row["bridges"][0]
        linked = random_row["linked_probe"]
        packing = next(
            item for item in linked["packing_lower_bounds"] if item["epsilon"] == 0.1
        )
        print(
            f"n={row['order']:2d} cap/n^1.5={row['normalized_absolute_cap']:.4f} "
            f"linked={linked['distinct_family_size']:4d} "
            f"median-gap/n^1.5="
            f"{linked['normalized_projective_distance']['median']:.4f} "
            f"pack(0.1)={packing['certified_greedy_packing_lower_bound']:4d}"
        )


if __name__ == "__main__":
    main()
