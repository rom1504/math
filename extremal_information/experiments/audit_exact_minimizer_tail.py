#!/usr/bin/env python3
"""Finite, counterexample-first audit of the ``L_tail`` proposal.

The matrix-bearing low-cap corpus, random seed, and control sizes come only
from the already frozen ``nearmin_blind_structural_results.json``.  No
minimizer search is performed.  See ``exact_minimizer_tail_finite_protocol.md``
for the preregistered orientation and threshold rules.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "extremal_information"
    / "experiments"
    / "nearmin_blind_structural_results.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "extremal_information"
    / "experiments"
    / "exact_minimizer_tail_finite_results.json"
)
D0_GRID = (Fraction(1, 64), Fraction(1, 32), Fraction(1, 16), Fraction(1, 8))
EXPECTED_INPUT_SHA256 = (
    "2c086cf7523ead804942948e800c6231eac33d954e049b5aa113c9fb0cca47a5"
)
EXPECTED_SMALL_EXACT_COUNTS = {3: 2, 4: 6, 5: 12, 6: 12, 7: 3240}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(matrix: np.ndarray) -> str:
    payload = json.dumps(matrix.astype(int).tolist(), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


@lru_cache(maxsize=None)
def landscape(n: int) -> tuple[np.ndarray, tuple[tuple[int, int], ...], np.ndarray]:
    rows = np.arange(1 << (n - 1), dtype=np.uint32)[:, None]
    bits = ((rows >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int16)
    spins = np.concatenate(
        [np.ones((len(rows), 1), dtype=np.int16), 1 - 2 * bits], axis=1
    )
    edges = tuple(itertools.combinations(range(n), 2))
    products = np.column_stack(
        [spins[:, i] * spins[:, j] for i, j in edges]
    ).astype(np.int16)
    return spins, edges, products


def matrix_from_signs(n: int, signs: np.ndarray) -> np.ndarray:
    _, edges, _ = landscape(n)
    matrix = np.zeros((n, n), dtype=np.int16)
    for (i, j), value in zip(edges, signs):
        matrix[i, j] = matrix[j, i] = int(value)
    return matrix


def signs_from_matrix(matrix: np.ndarray) -> np.ndarray:
    _, edges, _ = landscape(len(matrix))
    return np.asarray([matrix[i, j] for i, j in edges], dtype=np.int16)


def energies(matrix: np.ndarray) -> np.ndarray:
    _, _, products = landscape(len(matrix))
    return products @ signs_from_matrix(matrix)


def cyclic_matrix(n: int, pattern: int) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=np.int16)
    for i in range(n):
        for j in range(i + 1, n):
            distance = min((j - i) % n, (i - j) % n)
            sign = -1 if (pattern >> (distance - 1)) & 1 else 1
            matrix[i, j] = matrix[j, i] = sign
    return matrix


def root_gauged_population(n: int) -> Iterable[np.ndarray]:
    variable_edges = tuple(itertools.combinations(range(1, n), 2))
    for code in range(1 << len(variable_edges)):
        matrix = np.ones((n, n), dtype=np.int16)
        np.fill_diagonal(matrix, 0)
        for bit, (i, j) in enumerate(variable_edges):
            if (code >> bit) & 1:
                matrix[i, j] = matrix[j, i] = -1
        yield matrix


def exact_shell_count(deficits: np.ndarray, n: int, d0: Fraction) -> int:
    """Count ``deficit < d0*n^(3/2)`` using integer arithmetic."""
    lhs = (d0.denominator * deficits.astype(np.int64)) ** 2
    rhs = d0.numerator * d0.numerator * n**3
    return int(np.sum(lhs < rhs))


def d0_key(d0: Fraction) -> str:
    return f"{d0.numerator}/{d0.denominator}"


def matrix_record(matrix: np.ndarray, supplied_m: dict[int, int]) -> dict[str, Any]:
    matrix = np.asarray(matrix, dtype=np.int16)
    n = len(matrix)
    values = energies(matrix).astype(np.int64)
    p_plus = int(np.max(values))
    p_minus = int(-np.min(values))
    cap = max(p_plus, p_minus)
    spectral = float(np.max(np.abs(np.linalg.eigvalsh(matrix.astype(float)))) / math.sqrt(n))
    tails: dict[str, Any] = {}
    total_projective = len(values)
    for d0 in D0_GRID:
        plus_count = (
            exact_shell_count(cap - values, n, d0) if p_plus == cap else None
        )
        minus_count = (
            exact_shell_count(cap + values, n, d0) if p_minus == cap else None
        )
        eligible = []
        if plus_count is not None:
            eligible.append((plus_count, 1))
        if minus_count is not None:
            eligible.append((minus_count, -1))
        # Counterexample-first tie rule: maximize density, then use + on an
        # exact tie only for deterministic output.
        selected_count, selected_orientation = max(eligible, key=lambda x: (x[0], x[1]))
        density = selected_count / total_projective
        rate_bits = -math.log2(density) / n
        rate_nats = rate_bits * math.log(2)
        tails[d0_key(d0)] = {
            "positive_orientation_count": plus_count,
            "negative_orientation_count": minus_count,
            "selected_count_projective": selected_count,
            "selected_density": density,
            "selected_orientation": selected_orientation,
            "rate_bits_per_vertex": rate_bits,
            "rate_nats_per_vertex": rate_nats,
            "full_count_log_per_vertex_nats": math.log(2) - rate_nats,
        }
    return {
        "n": n,
        "matrix_sha256": stable_hash(matrix),
        "cap": cap,
        "cap_delta": cap - supplied_m[n],
        "positive_max": p_plus,
        "negative_max": p_minus,
        "absolute_cap_orientation_tie": p_plus == p_minus,
        "operator_norm_over_sqrt_n": spectral,
        "tails": tails,
    }


def quantiles(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    array = np.asarray(values, dtype=float)
    q = np.quantile(array, [0, 0.1, 0.5, 0.9, 1])
    return {
        "min": float(q[0]),
        "q10": float(q[1]),
        "median": float(q[2]),
        "q90": float(q[3]),
        "max": float(q[4]),
        "mean": float(np.mean(array)),
    }


def summarize_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    if not records:
        return {"count": 0}
    spectral_worst = max(records, key=lambda r: r["operator_norm_over_sqrt_n"])
    out: dict[str, Any] = {
        "count": len(records),
        "cap_histogram": {
            str(k): int(v) for k, v in sorted(Counter(r["cap"] for r in records).items())
        },
        "operator_norm_over_sqrt_n": quantiles(
            [r["operator_norm_over_sqrt_n"] for r in records]
        ),
        "operator_norm_worst_sha256": spectral_worst["matrix_sha256"],
        "orientation_tie_fraction": float(
            np.mean([r["absolute_cap_orientation_tie"] for r in records])
        ),
        "tails": {},
    }
    for d0 in D0_GRID:
        key = d0_key(d0)
        worst = min(records, key=lambda r: r["tails"][key]["rate_nats_per_vertex"])
        out["tails"][key] = {
            "rate_nats_per_vertex": quantiles(
                [r["tails"][key]["rate_nats_per_vertex"] for r in records]
            ),
            "rate_bits_per_vertex": quantiles(
                [r["tails"][key]["rate_bits_per_vertex"] for r in records]
            ),
            "density": quantiles(
                [r["tails"][key]["selected_density"] for r in records]
            ),
            "worst_sha256": worst["matrix_sha256"],
            "worst_selected_count_projective": worst["tails"][key][
                "selected_count_projective"
            ],
            "worst_selected_density": worst["tails"][key]["selected_density"],
            "worst_rate_nats_per_vertex": worst["tails"][key][
                "rate_nats_per_vertex"
            ],
        }
    return out


def deduplicate_labelled(rows: list[dict[str, Any]]) -> list[np.ndarray]:
    unique: dict[str, np.ndarray] = {}
    for row in rows:
        matrix = np.asarray(row["matrix"], dtype=np.int16)
        unique.setdefault(stable_hash(matrix), matrix)
    return list(unique.values())


def deduplicate_matrices(matrices: list[np.ndarray]) -> list[np.ndarray]:
    unique: dict[str, np.ndarray] = {}
    for matrix in matrices:
        unique.setdefault(stable_hash(matrix), matrix)
    return list(unique.values())


def compact_records(matrices: list[np.ndarray], supplied_m: dict[int, int]) -> list[dict[str, Any]]:
    return [matrix_record(matrix, supplied_m) for matrix in matrices]


def summarize_by_order(records: list[dict[str, Any]]) -> dict[str, Any]:
    groups: defaultdict[int, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        groups[record["n"]].append(record)
    return {str(n): summarize_records(rows) for n, rows in sorted(groups.items())}


def reconstruct_uniform_controls(
    frozen: dict[str, Any], supplied_m: dict[int, int]
) -> tuple[list[np.ndarray], dict[str, Any]]:
    seed = int(frozen["parameters"]["seed"])
    sample_count = int(frozen["parameters"]["random_samples_per_order"])
    max_order = max(supplied_m)
    rng = np.random.default_rng(seed)
    matrices: list[np.ndarray] = []
    checks: dict[str, Any] = {}
    for n in range(3, max_order + 1):
        order_rows = []
        edge_count = n * (n - 1) // 2
        for _ in range(sample_count):
            signs = rng.choice(np.asarray([-1, 1], dtype=np.int16), size=edge_count)
            matrix = matrix_from_signs(n, signs)
            order_rows.append(matrix)
            matrices.append(matrix)
        actual = Counter(int(np.max(np.abs(energies(matrix)))) for matrix in order_rows)
        expected = {
            int(k): int(v)
            for k, v in frozen["population_summaries"]["uniform_random"][str(n)][
                "cap_histogram"
            ].items()
        }
        checks[str(n)] = {
            "actual_cap_histogram": {str(k): v for k, v in sorted(actual.items())},
            "expected_cap_histogram": {str(k): v for k, v in sorted(expected.items())},
            "passed": dict(actual) == expected,
        }
        if dict(actual) != expected:
            raise AssertionError(f"uniform-control reconstruction failed at n={n}")
    return matrices, checks


def reconstruct_cyclic_controls(
    frozen: dict[str, Any], supplied_m: dict[int, int]
) -> tuple[list[np.ndarray], dict[str, Any]]:
    matrices: list[np.ndarray] = []
    checks: dict[str, Any] = {}
    for n in range(3, max(supplied_m) + 1):
        order_rows = [cyclic_matrix(n, pattern) for pattern in range(1 << (n // 2))]
        matrices.extend(order_rows)
        actual = Counter(int(np.max(np.abs(energies(matrix)))) for matrix in order_rows)
        expected = {
            int(k): int(v)
            for k, v in frozen["population_summaries"]["cyclic_distance"][str(n)][
                "cap_histogram"
            ].items()
        }
        checks[str(n)] = {
            "actual_cap_histogram": {str(k): v for k, v in sorted(actual.items())},
            "expected_cap_histogram": {str(k): v for k, v in sorted(expected.items())},
            "passed": dict(actual) == expected,
        }
        if dict(actual) != expected:
            raise AssertionError(f"cyclic-control reconstruction failed at n={n}")
    return matrices, checks


def witness_payload(record: dict[str, Any], matrix: np.ndarray) -> dict[str, Any]:
    return {"record": record, "matrix": matrix.astype(int).tolist()}


def global_witnesses(
    matrices: list[np.ndarray], supplied_m: dict[int, int]
) -> dict[str, Any]:
    pairs = [(matrix_record(matrix, supplied_m), matrix) for matrix in matrices]
    if not pairs:
        return {}
    output: dict[str, Any] = {}
    spectral_record, spectral_matrix = max(
        pairs, key=lambda pair: pair[0]["operator_norm_over_sqrt_n"]
    )
    output["maximum_operator_norm_over_sqrt_n"] = witness_payload(
        spectral_record, spectral_matrix
    )
    output["minimum_tail_rates"] = {}
    for d0 in D0_GRID:
        key = d0_key(d0)
        record, matrix = min(
            pairs, key=lambda pair: pair[0]["tails"][key]["rate_nats_per_vertex"]
        )
        output["minimum_tail_rates"][key] = witness_payload(record, matrix)
    return output


def exhaustive_small_orders(
    supplied_m: dict[int, int]
) -> tuple[dict[str, Any], dict[str, Any], list[np.ndarray], list[np.ndarray]]:
    summaries: dict[str, Any] = {}
    checks: dict[str, Any] = {}
    exact_matrices: list[np.ndarray] = []
    one_step_matrices: list[np.ndarray] = []
    for n in range(3, 8):
        strata: defaultdict[int, list[dict[str, Any]]] = defaultdict(list)
        matrices_by_hash: dict[str, np.ndarray] = {}
        for matrix in root_gauged_population(n):
            record = matrix_record(matrix, supplied_m)
            strata[record["cap"]].append(record)
            matrices_by_hash[record["matrix_sha256"]] = matrix.copy()
            if record["cap"] == supplied_m[n]:
                exact_matrices.append(matrix.copy())
            elif record["cap"] == supplied_m[n] + 2:
                one_step_matrices.append(matrix.copy())
        # Make both predeclared low-cap strata explicit even when lattice
        # parity leaves one empty (this occurs for n=3, cap M_n+2).
        strata.setdefault(supplied_m[n], [])
        strata.setdefault(supplied_m[n] + 2, [])
        exact_count = len(strata[supplied_m[n]])
        checks[str(n)] = {
            "root_gauged_population_size": sum(len(rows) for rows in strata.values()),
            "exact_count": exact_count,
            "expected_exact_count": EXPECTED_SMALL_EXACT_COUNTS[n],
            "passed": exact_count == EXPECTED_SMALL_EXACT_COUNTS[n],
        }
        if exact_count != EXPECTED_SMALL_EXACT_COUNTS[n]:
            raise AssertionError(f"exhaustive minimizer count failed at n={n}")
        order_summary: dict[str, Any] = {}
        for cap, records in sorted(strata.items()):
            summary = summarize_records(records)
            summary["classification"] = (
                "exact" if cap == supplied_m[n] else
                "one_step_near" if cap == supplied_m[n] + 2 else
                "other_cap"
            )
            # Retain the exact finite worst matrices for the two principal
            # low-cap strata only.
            if records and cap in {supplied_m[n], supplied_m[n] + 2}:
                hashes = {summary["operator_norm_worst_sha256"]}
                hashes.update(item["worst_sha256"] for item in summary["tails"].values())
                summary["worst_witness_matrices"] = {
                    digest: matrices_by_hash[digest].astype(int).tolist() for digest in sorted(hashes)
                }
            order_summary[str(cap)] = summary
        summaries[str(n)] = order_summary
    return summaries, checks, exact_matrices, one_step_matrices


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    input_sha = file_sha256(args.input)
    if input_sha != EXPECTED_INPUT_SHA256:
        raise AssertionError(
            f"frozen input hash mismatch: expected {EXPECTED_INPUT_SHA256}, got {input_sha}"
        )
    frozen = json.loads(args.input.read_text())
    supplied_m = {int(n): int(value) for n, value in frozen["supplied_exact_M"].items()}

    uniform, uniform_checks = reconstruct_uniform_controls(frozen, supplied_m)
    cyclic, cyclic_checks = reconstruct_cyclic_controls(frozen, supplied_m)

    stratum_matrices: dict[str, list[np.ndarray]] = {
        "repository_exact": deduplicate_labelled(frozen["repository_exact_representatives"]),
        "repository_one_step_near": deduplicate_labelled(
            frozen["repository_one_step_near_representatives"]
        ),
        "adversarial_low_cap": deduplicate_labelled(
            frozen["cap_constrained_adversarial_samples"]
        ),
        "greedy_low_cap": deduplicate_labelled(
            frozen["independently_generated_greedy_low_cap"]
        ),
        "uniform_random": deduplicate_matrices(uniform),
        "uniform_random_low_cap": deduplicate_labelled(
            frozen["random_draws_that_are_low_cap"]
        ),
        "cyclic_structured": deduplicate_matrices(cyclic),
    }

    compact: dict[str, list[dict[str, Any]]] = {}
    summaries: dict[str, Any] = {}
    witnesses: dict[str, Any] = {}
    for name, matrices in stratum_matrices.items():
        records = compact_records(matrices, supplied_m)
        compact[name] = records
        summaries[name] = {
            "pooled": summarize_records(records),
            "by_order": summarize_by_order(records),
        }
        if name in {"repository_exact", "repository_one_step_near", "adversarial_low_cap"}:
            witnesses[name] = global_witnesses(matrices, supplied_m)

    exhaustive, exhaustive_checks, exhaustive_exact, exhaustive_one_step = (
        exhaustive_small_orders(supplied_m)
    )
    witnesses["combined_exact_including_exhaustive_n_le_7"] = global_witnesses(
        deduplicate_matrices(stratum_matrices["repository_exact"] + exhaustive_exact),
        supplied_m,
    )
    witnesses["combined_one_step_including_exhaustive_n_le_7"] = global_witnesses(
        deduplicate_matrices(
            stratum_matrices["repository_one_step_near"] + exhaustive_one_step
        ),
        supplied_m,
    )

    result = {
        "status": "finite experimental audit only; no asymptotic theorem",
        "protocol": "exact_minimizer_tail_finite_protocol.md",
        "frozen_input": {
            "path": str(args.input.relative_to(ROOT)),
            "sha256": input_sha,
        },
        "d0_grid": [d0_key(d0) for d0 in D0_GRID],
        "orientation_policy": (
            "orient so max H=Q; on a two-sided absolute-cap tie use the orientation "
            "with larger upper-tail density separately at each d0"
        ),
        "rate_definition": (
            "I_nats=-(1/n) log p, where p is projective/full-cube upper-tail density; "
            "full count=exp((log 2-I_nats)n)"
        ),
        "stratum_summaries": summaries,
        "compact_records": compact,
        "global_worst_witnesses": witnesses,
        "exhaustive_root_gauged_orders_3_to_7": exhaustive,
        "checks": {
            "uniform_control_reconstruction": uniform_checks,
            "cyclic_control_reconstruction": cyclic_checks,
            "exhaustive_exact_counts": exhaustive_checks,
            "all_passed": True,
        },
        "limitations": [
            "orders are at most fourteen",
            "repository exact representatives are not exhaustive above order eight",
            "byte deduplication is not a switching/permutation orbit quotient",
            "finite positive rates do not imply a uniform asymptotic kappa",
            "finite spectral ratios do not imply a uniform operator-norm theorem",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
