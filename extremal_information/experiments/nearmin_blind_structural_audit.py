#!/usr/bin/env python3
"""Blind, reproducible structural audit of low-cap quadratic signings.

The observable set was frozen in ``nearmin_blind_observable_freeze.md`` before
the near-minimizer campaign prompt or theory drafts were read.  This script
uses only matrix-bearing computation results and independently generated
controls.  It never treats a heuristic search as an optimality certificate;
all caps reported here are recomputed exhaustively over projective spins.

The default run is intentionally modest: exhaustive root-gauged populations
through order 7, deterministic random and cyclic controls through order 14,
single-edge perturbation neighborhoods, local descent, and cap-constrained
adversarial walks.  No finite pattern emitted by this program is a theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "computations" / "results"
DEFAULT_OUTPUT = (
    ROOT / "extremal_information" / "experiments" / "nearmin_blind_structural_results.json"
)
M_EXACT = dict(zip(range(3, 15), (3, 4, 4, 5, 9, 10, 12, 13, 17, 18, 20, 21)))
SUMMARY_KEYS = (
    "cap_delta",
    "cap_over_sqrt_edges",
    "cap_over_n32",
    "d4_conference_defect",
    "dop_singular_deviation",
    "spectral_effective_rank_fraction",
    "pair_correlation_mean_abs",
    "pair_correlation_rms",
    "pair_correlation_max",
    "pair_correlation_zero_fraction",
    "triangle_bias_abs",
    "energy_m4_over_edges2",
    "energy_m6_over_edges3",
    "energy_entropy_bits",
    "boundary_mass_cap",
    "boundary_mass_cap_minus_2",
    "boundary_mass_cap_minus_4",
    "active_density",
    "active_frame_defect",
    "active_frame_operator",
    "local_field_effective_support_fraction",
    "local_field_max_share",
    "deletion_exact_fraction",
    "deletion_max_gap",
    "edge_improving_fraction",
    "edge_neutral_fraction",
)


def stable_hash(matrix: np.ndarray) -> str:
    payload = json.dumps(matrix.astype(int).tolist(), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


@lru_cache(maxsize=None)
def landscape(n: int) -> tuple[np.ndarray, tuple[tuple[int, int], ...], np.ndarray]:
    """Projective spins, edges, and edge characters for order ``n``."""
    rows = np.arange(1 << (n - 1), dtype=np.uint32)[:, None]
    bits = ((rows >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int16)
    spins = np.concatenate(
        [np.ones((len(rows), 1), dtype=np.int16), 1 - 2 * bits], axis=1
    )
    edges = tuple(combinations(range(n), 2))
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


def is_signing(value: Any) -> bool:
    if not isinstance(value, list) or not 3 <= len(value) <= 14:
        return False
    n = len(value)
    if not all(isinstance(row, list) and len(row) == n for row in value):
        return False
    try:
        matrix = np.asarray(value, dtype=np.int16)
    except (TypeError, ValueError):
        return False
    upper = matrix[np.triu_indices(n, 1)]
    return bool(
        np.array_equal(matrix, matrix.T)
        and np.all(np.diag(matrix) == 0)
        and np.all(np.abs(upper) == 1)
    )


def energies(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _, _, products = landscape(len(matrix))
    signs = signs_from_matrix(matrix)
    return products @ signs, signs


def exact_cap(matrix: np.ndarray) -> int:
    values, _ = energies(matrix)
    return int(np.max(np.abs(values)))


def quantile_summary(values: Iterable[float]) -> dict[str, float] | None:
    array = np.asarray(list(values), dtype=float)
    if not len(array):
        return None
    qs = np.quantile(array, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1])
    return {
        "min": float(qs[0]),
        "q10": float(qs[1]),
        "q25": float(qs[2]),
        "median": float(qs[3]),
        "q75": float(qs[4]),
        "q90": float(qs[5]),
        "max": float(qs[6]),
        "mean": float(np.mean(array)),
    }


def json_safe(value: Any) -> Any:
    """Replace non-finite floats so the output is strict JSON, not JSON+NaN."""
    if isinstance(value, dict):
        return {key: json_safe(child) for key, child in value.items()}
    if isinstance(value, list):
        return [json_safe(child) for child in value]
    if isinstance(value, tuple):
        return [json_safe(child) for child in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    if isinstance(value, np.integer):
        return int(value)
    return value


def metric_vector(matrix: np.ndarray, include_edge_response: bool = True) -> dict[str, Any]:
    """Compute all scalar frozen observables needed for population summaries."""
    matrix = np.asarray(matrix, dtype=np.int16)
    n = len(matrix)
    spins, edges, products = landscape(n)
    values, signs = energies(matrix)
    cap = int(np.max(np.abs(values)))
    edge_count = len(edges)
    abs_values = np.abs(values)
    active_mask = abs_values == cap
    active = spins[active_mask].astype(np.float64)
    active_count = int(np.sum(active_mask))

    matrix64 = matrix.astype(np.float64)
    square = matrix64 @ matrix64
    defect = square - (n - 1) * np.eye(n)
    d4 = float(np.sum(defect * defect) / (n * (n - 1) ** 2))
    eigenvalues = np.linalg.eigvalsh(matrix64)
    abs_eigenvalues = np.abs(eigenvalues)
    dop = float(np.max(np.abs(abs_eigenvalues / math.sqrt(n - 1) - 1)))
    spectral_weights = eigenvalues * eigenvalues / (n * (n - 1))
    positive_weights = spectral_weights[spectral_weights > 0]
    effective_rank = math.exp(float(-np.sum(positive_weights * np.log(positive_weights))))
    effective_rank_fraction = effective_rank / n

    upper = np.triu_indices(n, 1)
    pair_correlations = np.abs(square[upper]) / (n - 2)
    triangle_bias = float(
        abs(np.trace(matrix64 @ matrix64 @ matrix64)) / (n * (n - 1) * (n - 2))
    )

    values64 = values.astype(np.float64)
    histogram = Counter(int(v) for v in values)
    probabilities = np.asarray(list(histogram.values()), dtype=float) / len(values)
    entropy_bits = float(-np.sum(probabilities * np.log2(probabilities)))

    # The projective active frame is invariant up to signed permutation.
    active_covariance = active.T @ active / active_count
    active_frame_defect = float(
        np.sum((active_covariance - np.eye(n)) ** 2) / n
    )
    active_frame_operator = float(np.linalg.eigvalsh(active_covariance)[-1])

    fields = spins * (spins @ matrix)
    active_fields = fields[active_mask].astype(np.float64)
    active_signs = np.sign(values[active_mask]).astype(np.float64)
    z = active_signs[:, None] * active_fields / cap
    sum_z2 = np.sum(z * z, axis=1)
    effective_support_fraction = 4.0 / (n * sum_z2)

    deletion_caps = np.max(np.abs(values[:, None] - fields), axis=0).astype(int)
    child_target = M_EXACT.get(n - 1)
    if child_target is None:
        deletion_gaps: list[int] = []
        deletion_exact_fraction = float("nan")
        deletion_max_gap = float("nan")
    else:
        deletion_gaps = [int(v - child_target) for v in deletion_caps]
        deletion_exact_fraction = float(np.mean(deletion_caps == child_target))
        deletion_max_gap = float(np.max(deletion_caps - child_target))

    if include_edge_response:
        neighbor_values = values[:, None] - 2 * products * signs[None, :]
        neighbor_caps = np.max(np.abs(neighbor_values), axis=0).astype(int)
        neighbor_deltas = neighbor_caps - cap
        edge_improving_fraction = float(np.mean(neighbor_deltas < 0))
        edge_neutral_fraction = float(np.mean(neighbor_deltas == 0))
        edge_response_histogram = {
            str(int(k)): int(v) for k, v in sorted(Counter(neighbor_deltas.tolist()).items())
        }
    else:
        neighbor_caps = np.asarray([], dtype=int)
        edge_improving_fraction = float("nan")
        edge_neutral_fraction = float("nan")
        edge_response_histogram = {}

    determinant_scale = (n - 1) ** (n / 2)
    return {
        "n": n,
        "matrix_sha256": stable_hash(matrix),
        "cap": cap,
        "cap_delta": cap - M_EXACT[n],
        "cap_over_sqrt_edges": cap / math.sqrt(edge_count),
        "cap_over_n32": cap / (n ** 1.5),
        "d4_conference_defect": d4,
        "dop_singular_deviation": dop,
        "spectral_effective_rank_fraction": effective_rank_fraction,
        "pair_correlation_mean_abs": float(np.mean(pair_correlations)),
        "pair_correlation_rms": float(math.sqrt(np.mean(pair_correlations**2))),
        "pair_correlation_max": float(np.max(pair_correlations)),
        "pair_correlation_zero_fraction": float(np.mean(pair_correlations == 0)),
        "triangle_bias_abs": triangle_bias,
        "energy_m4_over_edges2": float(np.mean(values64**4) / edge_count**2),
        "energy_m6_over_edges3": float(np.mean(values64**6) / edge_count**3),
        "energy_entropy_bits": entropy_bits,
        "boundary_mass_cap": float(np.mean(abs_values == cap)),
        "boundary_mass_cap_minus_2": float(np.mean(abs_values >= cap - 2)),
        "boundary_mass_cap_minus_4": float(np.mean(abs_values >= cap - 4)),
        "active_count": active_count,
        "active_density": active_count / len(values),
        "active_frame_defect": active_frame_defect,
        "active_frame_operator": active_frame_operator,
        "local_field_effective_support_fraction": float(np.mean(effective_support_fraction)),
        "local_field_effective_support_fraction_min": float(
            np.min(effective_support_fraction)
        ),
        "local_field_max_share": float(np.max(z) / 2),
        "local_field_z_min": float(np.min(z)),
        "local_field_z_max": float(np.max(z)),
        "local_field_sum_max_error": float(np.max(np.abs(np.sum(z, axis=1) - 2))),
        "deletion_caps": [int(v) for v in deletion_caps],
        "deletion_gaps": deletion_gaps,
        "deletion_exact_fraction": deletion_exact_fraction,
        "deletion_max_gap": deletion_max_gap,
        "edge_improving_fraction": edge_improving_fraction,
        "edge_neutral_fraction": edge_neutral_fraction,
        "edge_response_histogram": edge_response_histogram,
        "absolute_determinant_normalized": float(
            abs(np.linalg.det(matrix64)) / determinant_scale
        ),
    }


def full_observables(matrix: np.ndarray) -> dict[str, Any]:
    """Scalar observables plus exact histograms and active-code geometry."""
    result = metric_vector(matrix, include_edge_response=True)
    n = len(matrix)
    spins, _, _ = landscape(n)
    values, _ = energies(matrix)
    cap = result["cap"]
    active = spins[np.abs(values) == cap].astype(np.int16)
    histogram = Counter(int(v) for v in values)
    result["energy_histogram"] = {
        str(k): int(v) for k, v in sorted(histogram.items())
    }
    result["absolute_energy_quantiles_over_cap"] = {
        str(q): float(v)
        for q, v in zip(
            (0.5, 0.75, 0.9, 0.95, 0.99),
            np.quantile(np.abs(values) / cap, (0.5, 0.75, 0.9, 0.95, 0.99)),
        )
    }
    result["eigenvalues"] = [float(v) for v in np.linalg.eigvalsh(matrix.astype(float))]

    overlap_histogram: Counter[int] = Counter()
    overlap_sum = 0.0
    overlap_pairs = 0
    overlap_max = 0.0
    # Chunking avoids a large active-by-active temporary for unusual controls.
    for start in range(0, len(active), 512):
        gram = np.abs(active[start : start + 512] @ active.T)
        for local_i, row in enumerate(gram):
            global_i = start + local_i
            if global_i + 1 >= len(active):
                continue
            tail = row[global_i + 1 :]
            overlap_histogram.update(int(v) for v in tail)
            overlap_sum += float(np.sum(tail)) / n
            overlap_pairs += len(tail)
            if len(tail):
                overlap_max = max(overlap_max, float(np.max(tail)) / n)
    result["active_overlap_abs_histogram"] = {
        str(k): int(v) for k, v in sorted(overlap_histogram.items())
    }
    result["active_overlap_abs_mean"] = (
        overlap_sum / overlap_pairs if overlap_pairs else None
    )
    result["active_overlap_abs_max_off_diagonal"] = (
        overlap_max if overlap_pairs else None
    )
    positive = int(np.sum(values == cap))
    negative = int(np.sum(values == -cap))
    result["cap_sign_counts_unordered"] = sorted((positive, negative))
    return result


@dataclass
class Extracted:
    matrix: np.ndarray
    sources: list[dict[str, str]]


def recursive_matrices(value: Any, keypath: str = "") -> Iterable[tuple[str, np.ndarray]]:
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{keypath}.{key}" if keypath else key
            if is_signing(child):
                yield path, np.asarray(child, dtype=np.int16)
            else:
                yield from recursive_matrices(child, path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from recursive_matrices(child, f"{keypath}[{index}]")


def extract_repository_matrices() -> dict[str, Extracted]:
    extracted: dict[str, Extracted] = {}
    for path in sorted(RESULTS.glob("*.json")):
        try:
            payload = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        classification = str(payload.get("classification", "unlabelled computation result"))
        for keypath, matrix in recursive_matrices(payload):
            digest = stable_hash(matrix)
            source = {
                "file": str(path.relative_to(ROOT)),
                "keypath": keypath,
                "classification": classification,
            }
            if digest in extracted:
                extracted[digest].sources.append(source)
            else:
                extracted[digest] = Extracted(matrix=matrix, sources=[source])
    return extracted


def invariant_signature(observation: dict[str, Any]) -> str:
    """A separating signature, explicitly not a complete orbit invariant."""
    histogram = observation["energy_histogram"]
    reflected = {str(-int(k)): v for k, v in histogram.items()}
    histogram_pair = min(
        json.dumps(histogram, sort_keys=True), json.dumps(reflected, sort_keys=True)
    )
    payload = {
        "abs_eigenvalues": sorted(round(abs(v), 8) for v in observation["eigenvalues"]),
        "energy_histogram_up_to_reflection": histogram_pair,
        "deletion_caps": sorted(observation["deletion_caps"]),
        "d4": round(observation["d4_conference_defect"], 12),
        "triangle": round(observation["triangle_bias_abs"], 12),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def summarize_group(records: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {"count": len(records)}
    for key in SUMMARY_KEYS:
        values = [r[key] for r in records if key in r and math.isfinite(float(r[key]))]
        summary[key] = quantile_summary(values)
    summary["cap_histogram"] = {
        str(k): int(v) for k, v in sorted(Counter(r["cap"] for r in records).items())
    }
    return summary


def root_gauged_population(n: int) -> Iterable[np.ndarray]:
    variable_edges = tuple(combinations(range(1, n), 2))
    for code in range(1 << len(variable_edges)):
        matrix = np.ones((n, n), dtype=np.int16)
        np.fill_diagonal(matrix, 0)
        for bit, (i, j) in enumerate(variable_edges):
            if (code >> bit) & 1:
                matrix[i, j] = matrix[j, i] = -1
        yield matrix


def exhaustive_small_orders(max_n: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output: dict[str, Any] = {}
    low_cap_records: list[dict[str, Any]] = []
    for n in range(3, min(max_n, 7) + 1):
        strata: defaultdict[int, list[dict[str, Any]]] = defaultdict(list)
        extreme: dict[str, tuple[float, np.ndarray, dict[str, Any]]] = {}
        for matrix in root_gauged_population(n):
            record = metric_vector(matrix)
            strata[record["cap"]].append(record)
            if record["cap"] <= M_EXACT[n] + 2:
                low_cap_records.append(record)
                candidates = {
                    "max_d4": record["d4_conference_defect"],
                    "min_boundary": -record["boundary_mass_cap"],
                    "max_active_frame_defect": record["active_frame_defect"],
                    "min_deletion_exact_fraction": -record["deletion_exact_fraction"]
                    if math.isfinite(record["deletion_exact_fraction"])
                    else -2,
                    "max_edge_improving_fraction": record["edge_improving_fraction"],
                }
                for name, score in candidates.items():
                    if name not in extreme or score > extreme[name][0]:
                        extreme[name] = (score, matrix.copy(), record)
        output[str(n)] = {
            "root_gauged_population_size": sum(len(v) for v in strata.values()),
            "cap_strata": {str(cap): summarize_group(rows) for cap, rows in sorted(strata.items())},
            "low_cap_extremes": {
                name: {
                    "matrix": matrix.astype(int).tolist(),
                    "observables": full_observables(matrix),
                }
                for name, (_, matrix, _) in extreme.items()
            },
        }
    return output, low_cap_records


def cyclic_matrix(n: int, pattern: int) -> np.ndarray:
    max_distance = n // 2
    matrix = np.zeros((n, n), dtype=np.int16)
    for i in range(n):
        for j in range(i + 1, n):
            distance = min((j - i) % n, (i - j) % n)
            sign = -1 if (pattern >> (distance - 1)) & 1 else 1
            matrix[i, j] = matrix[j, i] = sign
    assert max_distance == max(min((j - i) % n, (i - j) % n) for i in range(n) for j in range(i + 1, n))
    return matrix


def steepest_descent(matrix: np.ndarray) -> np.ndarray:
    n = len(matrix)
    _, edges, products = landscape(n)
    signs = signs_from_matrix(matrix).copy()
    values = products @ signs
    while True:
        neighbor_values = values[:, None] - 2 * products * signs[None, :]
        neighbor_caps = np.max(np.abs(neighbor_values), axis=0)
        current_cap = int(np.max(np.abs(values)))
        edge = int(np.argmin(neighbor_caps))
        if int(neighbor_caps[edge]) >= current_cap:
            break
        values = neighbor_values[:, edge].copy()
        signs[edge] *= -1
    return matrix_from_signs(n, signs)


def cap_constrained_walk(
    seed_matrix: np.ndarray,
    threshold: int,
    objective: str,
    steps: int,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Search for atypical observables while never leaving the low-cap set."""
    n = len(seed_matrix)
    spins, edges, products = landscape(n)
    signs = signs_from_matrix(seed_matrix).copy()
    values = products @ signs
    matrix = seed_matrix.copy()

    def score(candidate_matrix: np.ndarray, candidate_values: np.ndarray) -> float:
        cap = int(np.max(np.abs(candidate_values)))
        active_mask = np.abs(candidate_values) == cap
        if objective == "d4_max":
            square = candidate_matrix @ candidate_matrix
            defect = square - (n - 1) * np.eye(n)
            return float(np.sum(defect * defect) / (n * (n - 1) ** 2))
        if objective == "boundary_min":
            return -float(np.mean(active_mask))
        if objective == "active_frame_max":
            active = spins[active_mask].astype(float)
            covariance = active.T @ active / len(active)
            return float(np.sum((covariance - np.eye(n)) ** 2) / n)
        if objective == "field_support_min":
            fields = spins * (spins @ candidate_matrix)
            z = (
                np.sign(candidate_values[active_mask])[:, None]
                * fields[active_mask]
                / cap
            )
            support = 4 / (n * np.sum(z * z, axis=1))
            return -float(np.mean(support))
        if objective == "random":
            return 0.0
        raise ValueError(objective)

    current_score = score(matrix, values)
    best_score = current_score
    best = matrix.copy()
    accepted = 0
    samples = [matrix.copy()]
    temperature = 0.025
    for _ in range(steps):
        edge = int(rng.integers(len(edges)))
        proposed_values = values - 2 * products[:, edge] * signs[edge]
        if int(np.max(np.abs(proposed_values))) > threshold:
            continue
        i, j = edges[edge]
        matrix[i, j] *= -1
        matrix[j, i] *= -1
        proposed_score = score(matrix, proposed_values)
        delta = proposed_score - current_score
        accept = objective == "random" or delta >= 0 or rng.random() < math.exp(
            delta / temperature
        )
        if accept:
            signs[edge] *= -1
            values = proposed_values
            current_score = proposed_score
            accepted += 1
            if proposed_score > best_score:
                best_score = proposed_score
                best = matrix.copy()
            if accepted % 75 == 0 and len(samples) < 25:
                samples.append(matrix.copy())
        else:
            matrix[i, j] *= -1
            matrix[j, i] *= -1
    samples.append(best)
    unique: dict[str, np.ndarray] = {}
    for candidate in samples:
        unique[stable_hash(candidate)] = candidate
    return list(unique.values())


def transformed(matrix: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = len(matrix)
    permutation = rng.permutation(n)
    switches = rng.choice(np.asarray([-1, 1], dtype=np.int16), size=n)
    sign = int(rng.choice(np.asarray([-1, 1], dtype=np.int16)))
    permuted = matrix[np.ix_(permutation, permutation)]
    return sign * switches[:, None] * permuted * switches[None, :]


def invariance_checks(matrices: list[np.ndarray], rng: np.random.Generator) -> dict[str, Any]:
    numeric_keys = [
        key
        for key in SUMMARY_KEYS
        if key not in {"cap_delta", "deletion_exact_fraction", "deletion_max_gap"}
    ]
    maximum_error = 0.0
    checked = 0
    exact_failures: list[str] = []
    for matrix in matrices[:12]:
        first = metric_vector(matrix)
        second = metric_vector(transformed(matrix, rng))
        checked += 1
        if first["cap"] != second["cap"]:
            exact_failures.append(f"cap:{first['matrix_sha256']}")
        if sorted(first["deletion_caps"]) != sorted(second["deletion_caps"]):
            exact_failures.append(f"deletion:{first['matrix_sha256']}")
        if first["edge_response_histogram"] != second["edge_response_histogram"]:
            exact_failures.append(f"edge:{first['matrix_sha256']}")
        for key in numeric_keys:
            a, b = float(first[key]), float(second[key])
            if math.isfinite(a) and math.isfinite(b):
                maximum_error = max(maximum_error, abs(a - b))
    return {
        "matrices_checked": checked,
        "maximum_numeric_error": maximum_error,
        "exact_failures": exact_failures,
        "passed": not exact_failures and maximum_error < 1e-10,
    }


def select_extreme(
    labelled: list[dict[str, Any]], key: str, maximize: bool
) -> dict[str, Any] | None:
    candidates = [x for x in labelled if math.isfinite(float(x["observables"][key]))]
    if not candidates:
        return None
    chosen = max(candidates, key=lambda x: x["observables"][key]) if maximize else min(
        candidates, key=lambda x: x["observables"][key]
    )
    return chosen


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--random-samples", type=int, default=48)
    parser.add_argument("--greedy-restarts", type=int, default=10)
    parser.add_argument("--walk-steps", type=int, default=3000)
    parser.add_argument("--max-order", type=int, default=14)
    args = parser.parse_args()
    started = time.time()
    rng = np.random.default_rng(args.seed)

    extracted = extract_repository_matrices()
    repository_exact: list[dict[str, Any]] = []
    repository_near: list[dict[str, Any]] = []
    exact_seed_by_n: dict[int, np.ndarray] = {}
    for digest, item in extracted.items():
        n = len(item.matrix)
        if n > args.max_order:
            continue
        cap = exact_cap(item.matrix)
        if cap not in {M_EXACT[n], M_EXACT[n] + 2}:
            continue
        observation = full_observables(item.matrix)
        labelled = {
            "label": "repository exact minimizer" if cap == M_EXACT[n] else "repository one-step near-minimizer",
            "matrix": item.matrix.astype(int).tolist(),
            "observables": observation,
            "sources": item.sources,
        }
        if cap == M_EXACT[n]:
            repository_exact.append(labelled)
            exact_seed_by_n.setdefault(n, item.matrix.copy())
        else:
            repository_near.append(labelled)

    # Authoritative exhaustive orbit inventory through order 8.
    orbit_inventory: dict[str, Any] = {}
    orbit_representatives: list[dict[str, Any]] = []
    for n in range(3, min(8, args.max_order) + 1):
        path = RESULTS / f"m{n}_minimizer_orbits.json"
        payload = json.loads(path.read_text())
        classes = []
        for row in payload["classes"]:
            matrix = np.asarray(row["representative_matrix"], dtype=np.int16)
            observation = full_observables(matrix)
            labelled = {
                "label": "authoritative exhaustive orbit representative",
                "matrix": matrix.astype(int).tolist(),
                "observables": observation,
                "sources": [
                    {
                        "file": str(path.relative_to(ROOT)),
                        "keypath": f"classes[{row['class']}].representative_matrix",
                        "classification": payload["classification"],
                    }
                ],
                "orbit_class": int(row["class"]),
                "orbit_sha256": row["canonical_orbit_sha256"],
                "root_gauged_labeled_count": int(row["root_gauged_labeled_count"]),
                "self_complementary": bool(row["self_complementary"]),
            }
            orbit_representatives.append(labelled)
            classes.append(labelled)
        orbit_inventory[str(n)] = {
            "classification": payload["classification"],
            "target_cap": int(payload["target_cap"]),
            "root_gauged_population_size": int(payload["root_gauged_signing_count"]),
            "minimizing_signing_count": int(payload["minimizing_signing_count"]),
            "signed_permutation_and_global_sign_class_count": int(
                payload["signed_permutation_and_global_sign_class_count"]
            ),
            "classes": classes,
        }

    exact_small, exhaustive_low_cap = exhaustive_small_orders(args.max_order)
    exhaustive_extreme_labelled: list[dict[str, Any]] = []
    for n, row in exact_small.items():
        for objective, item in row["low_cap_extremes"].items():
            exhaustive_extreme_labelled.append(
                {
                    "label": "exhaustive root-gauged low-cap extreme",
                    "objective": objective,
                    "matrix": item["matrix"],
                    "observables": item["observables"],
                    "population_order": int(n),
                }
            )

    # Random and cyclic controls.
    control_records: defaultdict[str, defaultdict[int, list[dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    control_extremes: list[dict[str, Any]] = []
    random_low_cap_labelled: list[dict[str, Any]] = []
    cyclic_low_cap_labelled: list[dict[str, Any]] = []
    for n in range(3, args.max_order + 1):
        edge_count = n * (n - 1) // 2
        for sample_index in range(args.random_samples):
            signs = rng.choice(np.asarray([-1, 1], dtype=np.int16), size=edge_count)
            matrix = matrix_from_signs(n, signs)
            record = metric_vector(matrix)
            control_records["uniform_random"][n].append(record)
            if record["cap"] <= M_EXACT[n] + 2:
                random_low_cap_labelled.append(
                    {
                        "label": "uniform-random draw that happened to be low-cap",
                        "sample_index": sample_index,
                        "matrix": matrix.astype(int).tolist(),
                        "observables": full_observables(matrix),
                    }
                )
        cyclic_rows: list[tuple[np.ndarray, dict[str, Any]]] = []
        for pattern in range(1 << (n // 2)):
            matrix = cyclic_matrix(n, pattern)
            record = metric_vector(matrix)
            control_records["cyclic_distance"][n].append(record)
            cyclic_rows.append((matrix, record))
            if record["cap"] <= M_EXACT[n] + 2:
                cyclic_low_cap_labelled.append(
                    {
                        "label": "cyclic-distance low-cap control",
                        "pattern": pattern,
                        "matrix": matrix.astype(int).tolist(),
                        "observables": full_observables(matrix),
                    }
                )
        best_matrix, best_record = min(cyclic_rows, key=lambda item: item[1]["cap"])
        control_extremes.append(
            {
                "label": "best cyclic-distance control",
                "matrix": best_matrix.astype(int).tolist(),
                "observables": full_observables(best_matrix),
            }
        )
        all_one = np.ones((n, n), dtype=np.int16) - np.eye(n, dtype=np.int16)
        control_records["all_one_switching_class"][n].append(metric_vector(all_one))

    # Independent random-restart steepest local minima.
    greedy_labelled: list[dict[str, Any]] = []
    for n in range(8, args.max_order + 1):
        edge_count = n * (n - 1) // 2
        for restart in range(args.greedy_restarts):
            signs = rng.choice(np.asarray([-1, 1], dtype=np.int16), size=edge_count)
            matrix = steepest_descent(matrix_from_signs(n, signs))
            observation = metric_vector(matrix)
            control_records["independent_greedy_local_minimum"][n].append(observation)
            if observation["cap"] <= M_EXACT[n] + 2:
                greedy_labelled.append(
                    {
                        "label": "independently generated greedy low-cap witness",
                        "matrix": matrix.astype(int).tolist(),
                        "observables": full_observables(matrix),
                        "seed": args.seed,
                        "restart": restart,
                    }
                )

    # Every single edge flip of an exact minimizer lies within cap M_n+2.
    one_edge_records: defaultdict[int, list[dict[str, Any]]] = defaultdict(list)
    one_edge_labelled: list[dict[str, Any]] = []
    for n, seed_matrix in sorted(exact_seed_by_n.items()):
        if n < 7:
            continue
        _, edges, _ = landscape(n)
        for edge_index, (i, j) in enumerate(edges):
            matrix = seed_matrix.copy()
            matrix[i, j] *= -1
            matrix[j, i] *= -1
            observation = metric_vector(matrix)
            if observation["cap"] > M_EXACT[n] + 2:
                raise AssertionError("single-edge Lipschitz check failed")
            one_edge_records[n].append(observation)
            one_edge_labelled.append(
                {
                    "label": "one-edge perturbation of an exact minimizer",
                    "matrix": matrix.astype(int).tolist(),
                    "observables": observation,
                    "parent_sha256": stable_hash(seed_matrix),
                    "flipped_edge": [i, j],
                    "flipped_edge_index": edge_index,
                }
            )

    # Adversarial walks stay entirely inside cap <= M_n+2.
    walk_labelled: list[dict[str, Any]] = []
    objectives = (
        "d4_max",
        "boundary_min",
        "active_frame_max",
        "field_support_min",
        "random",
    )
    for n, seed_matrix in sorted(exact_seed_by_n.items()):
        if n < 8:
            continue
        for objective in objectives:
            matrices = cap_constrained_walk(
                seed_matrix,
                M_EXACT[n] + 2,
                objective,
                args.walk_steps,
                rng,
            )
            for index, matrix in enumerate(matrices):
                observation = full_observables(matrix)
                if observation["cap"] > M_EXACT[n] + 2:
                    raise AssertionError("walk escaped cap threshold")
                walk_labelled.append(
                    {
                        "label": "cap-constrained adversarial walk sample",
                        "objective": objective,
                        "sample_index": index,
                        "matrix": matrix.astype(int).tolist(),
                        "observables": observation,
                        "parent_sha256": stable_hash(seed_matrix),
                    }
                )

    # Distinguishable signatures among available exact representatives.  This
    # is deliberately weaker than orbit classification beyond order 8.
    available_exact_inventory: dict[str, Any] = {}
    for n in range(3, args.max_order + 1):
        rows = [x for x in repository_exact if x["observables"]["n"] == n]
        signatures = defaultdict(list)
        for row in rows:
            signature = invariant_signature(row["observables"])
            signatures[signature].append(row["observables"]["matrix_sha256"])
        available_exact_inventory[str(n)] = {
            "byte_distinct_representatives": len(rows),
            "distinguishable_invariant_signature_count": len(signatures),
            "warning": "signature collisions are not asserted to be the same switching orbit",
            "signature_classes": dict(signatures),
        }

    population_summaries: dict[str, Any] = {
        "repository_exact": {},
        "repository_one_step_near": {},
        "authoritative_orbit_representatives": {},
        "one_edge_neighborhood": {},
        "cap_constrained_walk": {},
        "independent_greedy_low_cap": {},
        "random_draw_low_cap": {},
        "cyclic_distance_low_cap": {},
        "exhaustive_low_cap_extremes": {},
    }
    labelled_groups = {
        "repository_exact": repository_exact,
        "repository_one_step_near": repository_near,
        "authoritative_orbit_representatives": orbit_representatives,
        "one_edge_neighborhood": one_edge_labelled,
        "cap_constrained_walk": walk_labelled,
        "independent_greedy_low_cap": greedy_labelled,
        "random_draw_low_cap": random_low_cap_labelled,
        "cyclic_distance_low_cap": cyclic_low_cap_labelled,
        "exhaustive_low_cap_extremes": exhaustive_extreme_labelled,
    }
    for group, rows in labelled_groups.items():
        for n in range(3, args.max_order + 1):
            records = [x["observables"] for x in rows if x["observables"]["n"] == n]
            if records:
                population_summaries[group][str(n)] = summarize_group(records)
    for group, by_n in control_records.items():
        population_summaries[group] = {
            str(n): summarize_group(records) for n, records in sorted(by_n.items())
        }

    all_low_cap_labelled = (
        repository_exact
        + repository_near
        + orbit_representatives
        + one_edge_labelled
        + walk_labelled
        + greedy_labelled
        + random_low_cap_labelled
        + cyclic_low_cap_labelled
        + exhaustive_extreme_labelled
    )
    exact_only_labelled = [
        x for x in all_low_cap_labelled if x["observables"]["cap_delta"] == 0
    ]
    near_only_labelled = [
        x for x in all_low_cap_labelled if x["observables"]["cap_delta"] == 2
    ]
    falsifiers = {
        "H1_max_d4_among_exact": select_extreme(
            exact_only_labelled, "d4_conference_defect", True
        ),
        "H1_max_d4_among_one_step_near": select_extreme(
            near_only_labelled, "d4_conference_defect", True
        ),
        "H2_min_boundary_mass_among_exact": select_extreme(
            exact_only_labelled, "boundary_mass_cap", False
        ),
        "H2_min_boundary_mass_among_one_step_near": select_extreme(
            near_only_labelled, "boundary_mass_cap", False
        ),
        "H2_max_active_frame_defect_among_one_step_near": select_extreme(
            near_only_labelled, "active_frame_defect", True
        ),
        "H3_min_deletion_exact_fraction_among_exact": select_extreme(
            exact_only_labelled, "deletion_exact_fraction", False
        ),
        "H3_min_deletion_exact_fraction_among_one_step_near": select_extreme(
            near_only_labelled, "deletion_exact_fraction", False
        ),
        "H4_max_improving_edge_fraction_among_one_step_near": select_extreme(
            near_only_labelled, "edge_improving_fraction", True
        ),
        "local_field_min_effective_support_among_one_step_near": select_extreme(
            near_only_labelled, "local_field_effective_support_fraction", False
        ),
    }

    exact_matrices_for_checks = [
        np.asarray(x["matrix"], dtype=np.int16) for x in repository_exact
    ]
    checks = {
        "switching_permutation_global_sign_invariance": invariance_checks(
            exact_matrices_for_checks, rng
        ),
        "all_exact_caps_match_supplied_M": all(
            x["observables"]["cap"] == M_EXACT[x["observables"]["n"]]
            for x in repository_exact
        ),
        "all_one_edge_perturbations_within_M_plus_2": all(
            x["observables"]["cap"] <= M_EXACT[x["observables"]["n"]] + 2
            for x in one_edge_labelled
        ),
        "all_local_field_constraints_hold": all(
            x["observables"]["local_field_z_min"] >= -1e-12
            and x["observables"]["local_field_z_max"] <= 1 + 1e-12
            and x["observables"]["local_field_sum_max_error"] <= 1e-12
            for x in all_low_cap_labelled
        ),
    }
    orbit_enumeration_count_checks = {}
    for n, row in exact_small.items():
        expected = orbit_inventory[n]["minimizing_signing_count"]
        actual = row["cap_strata"][str(M_EXACT[int(n)])]["count"]
        orbit_enumeration_count_checks[n] = {
            "expected": expected,
            "actual": actual,
            "passed": expected == actual,
        }
    checks["independent_exhaustive_minimizer_counts"] = orbit_enumeration_count_checks

    payload = {
        "schema": "nearmin-blind-structural-audit-v1",
        "status": "FINITE EXPERIMENTAL AUDIT; NOT A THEOREM OR ASYMPTOTIC CLAIM",
        "blind_freeze": "extremal_information/experiments/nearmin_blind_observable_freeze.md",
        "objective": "Q(A)=max_x |sum_{i<j} a_ij x_i x_j|",
        "supplied_exact_M": {str(k): v for k, v in M_EXACT.items()},
        "parameters": {
            "seed": args.seed,
            "random_samples_per_order": args.random_samples,
            "greedy_restarts_per_order_8_plus": args.greedy_restarts,
            "walk_steps_per_objective_and_order": args.walk_steps,
            "max_order": args.max_order,
        },
        "classification_policy": {
            "cap": "exhaustively recomputed over all 2^(n-1) projective spins",
            "exact": "cap equals the supplied exact M_n",
            "one_step_near": "cap equals M_n+2, the next parity-compatible cap",
            "orbit_warning": "only computations/results/m3..m8_minimizer_orbits.json claim exhaustive orbit classification",
            "search_warning": "greedy and cap-constrained-walk discovery is heuristic even though each saved cap is exact",
        },
        "authoritative_orbit_inventory": orbit_inventory,
        "available_exact_representative_inventory": available_exact_inventory,
        "repository_exact_representatives": repository_exact,
        "repository_one_step_near_representatives": repository_near,
        "independently_generated_greedy_low_cap": greedy_labelled,
        "random_draws_that_are_low_cap": random_low_cap_labelled,
        "cyclic_distance_low_cap_controls": cyclic_low_cap_labelled,
        "exhaustive_low_cap_extreme_representatives": exhaustive_extreme_labelled,
        "one_edge_low_cap_neighborhood": one_edge_labelled,
        "cap_constrained_adversarial_samples": walk_labelled,
        "control_extremes": control_extremes,
        "population_summaries": population_summaries,
        "exhaustive_root_gauged_orders_3_to_7": exact_small,
        "falsifier_shortlist": falsifiers,
        "checks": checks,
        "elapsed_seconds": time.time() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(json_safe(payload), indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(f"repository exact representatives: {len(repository_exact)}")
    print(f"repository one-step near representatives: {len(repository_near)}")
    print(f"one-edge low-cap samples: {len(one_edge_labelled)}")
    print(f"walk samples: {len(walk_labelled)}")
    print(f"greedy low-cap samples: {len(greedy_labelled)}")
    print(f"elapsed_seconds: {time.time() - started:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
