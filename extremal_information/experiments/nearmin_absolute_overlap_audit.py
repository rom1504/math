#!/usr/bin/env python3
"""Finite absolute-overlap audit for positive near-top shells.

The protocol was frozen in ``nearmin_absolute_overlap_observable_freeze.md``
before this program was run.  Every cap and shell is enumerated exactly over
projective Boolean spins.  Maximum packings are certified only when CP-SAT
returns OPTIMAL; otherwise the output retains explicit lower/upper bounds.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "extremal_information/experiments/nearmin_blind_structural_results.json"
DEFAULT_OUTPUT = ROOT / "extremal_information/experiments/nearmin_absolute_overlap_results.json"
DEFAULT_SUMMARY_OUTPUT = (
    ROOT / "extremal_information/experiments/nearmin_absolute_overlap_summary.json"
)
SEED = 20260817
RANDOM_SIGNINGS_PER_ORDER = 24
MATCHED_SUBSETS_PER_SHELL = 32
SUPPLEMENT_PER_STRATUM_ORDER = 8
EXACT_PACKING_SIZE_LIMIT = 120
EXACT_PACKING_SECONDS = 1.0


def edge_list(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(n), 2))


def projective_spins(n: int) -> np.ndarray:
    masks = np.arange(1 << (n - 1), dtype=np.uint32)
    bits = ((masks[:, None] >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int8)
    spins = np.ones((len(masks), n), dtype=np.int8)
    spins[:, 1:] = 1 - 2 * bits
    return spins


def upper_signs(matrix: np.ndarray) -> np.ndarray:
    return np.asarray([matrix[i, j] for i, j in edge_list(len(matrix))], dtype=np.int8)


def matrix_from_signs(n: int, signs: np.ndarray) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=np.int8)
    for (i, j), value in zip(edge_list(n), signs):
        matrix[i, j] = matrix[j, i] = int(value)
    return matrix


def matrix_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(upper_signs(matrix).tobytes()).hexdigest()


def energy_vector(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    signs = upper_signs(matrix).astype(np.int16)
    products = np.stack(
        [spins[:, i] * spins[:, j] for i, j in edge_list(len(matrix))], axis=1
    ).astype(np.int16)
    return products @ signs


def quantiles(values: Iterable[float]) -> dict[str, float] | None:
    array = np.asarray(list(values), dtype=float)
    if not len(array):
        return None
    result = np.quantile(array, [0, 0.25, 0.5, 0.75, 1])
    return {
        "min": float(result[0]),
        "q25": float(result[1]),
        "median": float(result[2]),
        "q75": float(result[3]),
        "max": float(result[4]),
        "mean": float(np.mean(array)),
    }


def shell_points(energies: np.ndarray, cap: int, deficit: int) -> tuple[np.ndarray, np.ndarray]:
    threshold = cap - deficit
    if threshold <= 0:
        raise ValueError("frozen shells require a positive threshold")
    plus = np.flatnonzero(energies >= threshold)
    minus = np.flatnonzero(energies <= -threshold)
    indices = np.concatenate([plus, minus]).astype(np.int32)
    orientations = np.concatenate(
        [np.ones(len(plus), dtype=np.int8), -np.ones(len(minus), dtype=np.int8)]
    )
    if len(np.unique(indices)) != len(indices):
        raise AssertionError("positive threshold should give one orientation per projective spin")
    return indices, orientations


def pair_geometry(spins: np.ndarray) -> dict[str, Any]:
    """Exact pair data and compatibility matrix for one projective subset."""
    count, n = spins.shape
    if count <= 1:
        return {
            "pair_count": 0,
            "edge_abs_overlap": None,
            "vertex_projective_overlap": None,
            "edge_fraction_gt_half": None,
            "edge_fraction_ge_three_quarters": None,
            "edge_fraction_ge_nine_tenths": None,
            "compatible": np.ones((count, count), dtype=bool),
        }
    gram = spins.astype(np.int16) @ spins.astype(np.int16).T
    iu = np.triu_indices(count, 1)
    dot = gram[iu].astype(np.int64)
    edge_num = np.abs(dot * dot - n)
    edge_den = n * (n - 1)
    edge_overlap = edge_num.astype(float) / edge_den
    vertex_overlap = np.abs(dot).astype(float) / n
    reconstructed = np.abs((n * vertex_overlap * vertex_overlap - 1) / (n - 1))
    if not np.allclose(edge_overlap, reconstructed, atol=1e-14, rtol=0):
        raise AssertionError("edge/vertex overlap identity failed")
    # The comparison is performed with integer arithmetic: 2*num <= den.
    all_num = np.abs(gram.astype(np.int64) ** 2 - n)
    compatible = 2 * all_num <= edge_den
    np.fill_diagonal(compatible, True)
    return {
        "pair_count": int(len(dot)),
        "edge_abs_overlap": quantiles(edge_overlap),
        "vertex_projective_overlap": quantiles(vertex_overlap),
        "edge_fraction_gt_half": float(np.mean(2 * edge_num > edge_den)),
        "edge_fraction_ge_three_quarters": float(np.mean(4 * edge_num >= 3 * edge_den)),
        "edge_fraction_ge_nine_tenths": float(np.mean(10 * edge_num >= 9 * edge_den)),
        "compatible": compatible,
    }


def greedy_packing(compatible: np.ndarray) -> list[int]:
    """Deterministic maximal packing, repeatedly choosing least conflict."""
    count = len(compatible)
    remaining = np.ones(count, dtype=bool)
    chosen: list[int] = []
    while np.any(remaining):
        live = np.flatnonzero(remaining)
        # Among live vertices, minimizing conflicts is the same as maximizing
        # compatibility degree.  ``argmax`` gives the frozen least-index tie
        # break because ``live`` is sorted.
        degrees = compatible[np.ix_(live, live)].sum(axis=1)
        vertex = int(live[int(np.argmax(degrees))])
        chosen.append(vertex)
        remaining &= compatible[vertex]
        remaining[vertex] = False
    return chosen


def certified_packing(compatible: np.ndarray, greedy: list[int]) -> dict[str, Any]:
    count = len(compatible)
    if count <= 1:
        return {
            "greedy_lower_bound": count,
            "certified_lower_bound": count,
            "certified_upper_bound": count,
            "exact": True,
            "optimum": count,
        }
    if len(greedy) == count:
        # Every pair selected by a full-size greedy packing is compatible, so
        # no solver is needed for an exact certificate.
        return {
            "greedy_lower_bound": count,
            "certified_lower_bound": count,
            "certified_upper_bound": count,
            "exact": True,
            "optimum": count,
            "status": "ALL_PAIRS_COMPATIBLE",
        }
    if count > EXACT_PACKING_SIZE_LIMIT:
        return {
            "greedy_lower_bound": len(greedy),
            "certified_lower_bound": len(greedy),
            "certified_upper_bound": count,
            "exact": False,
            "optimum": None,
            "status": "SKIPPED_SIZE_LIMIT",
        }
    model = cp_model.CpModel()
    variables = [model.NewBoolVar(f"x{i}") for i in range(count)]
    for i in range(count):
        for j in range(i + 1, count):
            if not compatible[i, j]:
                model.Add(variables[i] + variables[j] <= 1)
    model.Maximize(sum(variables))
    greedy_set = set(greedy)
    for i, variable in enumerate(variables):
        model.AddHint(variable, int(i in greedy_set))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = EXACT_PACKING_SECONDS
    solver.parameters.num_search_workers = 8
    solver.parameters.random_seed = SEED
    status = solver.Solve(model)
    names = {
        cp_model.OPTIMAL: "OPTIMAL",
        cp_model.FEASIBLE: "FEASIBLE",
        cp_model.INFEASIBLE: "INFEASIBLE",
        cp_model.MODEL_INVALID: "MODEL_INVALID",
        cp_model.UNKNOWN: "UNKNOWN",
    }
    lower = len(greedy)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        lower = max(lower, int(round(solver.ObjectiveValue())))
    upper = count
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE, cp_model.UNKNOWN):
        upper = min(upper, int(math.floor(solver.BestObjectiveBound() + 1e-7)))
    upper = max(lower, upper)
    exact = status == cp_model.OPTIMAL
    return {
        "greedy_lower_bound": len(greedy),
        "certified_lower_bound": lower,
        "certified_upper_bound": upper,
        "exact": exact,
        "optimum": lower if exact else None,
        "status": names.get(status, str(status)),
    }


def stable_rng(*parts: object) -> np.random.Generator:
    digest = hashlib.sha256(":".join(map(str, parts)).encode()).digest()
    seed = int.from_bytes(digest[:8], "little")
    return np.random.default_rng(seed)


@lru_cache(maxsize=None)
def matched_controls(n: int, size: int) -> dict[str, Any]:
    universe = 1 << (n - 1)
    if size > universe:
        raise AssertionError("shell exceeds projective cube")
    spins = projective_spins(n)
    minimum_edge: list[float | None] = []
    packing: list[int] = []
    for repetition in range(MATCHED_SUBSETS_PER_SHELL):
        # Use a common standardized null panel for equal ``(n,size)``.  Every
        # shell still receives 32 fixed-seed uniform cardinality-matched
        # subsets, while caching prevents identical null experiments from
        # dominating run time.
        rng = stable_rng(SEED, "matched", n, size, repetition)
        indices = rng.choice(universe, size=size, replace=False)
        geometry = pair_geometry(spins[indices])
        edge = geometry["edge_abs_overlap"]
        minimum_edge.append(None if edge is None else float(edge["min"]))
        packing.append(len(greedy_packing(geometry.pop("compatible"))))
    finite_minimum = [x for x in minimum_edge if x is not None]
    return {
        "replicates": MATCHED_SUBSETS_PER_SHELL,
        "minimum_edge_abs_overlap": quantiles(finite_minimum),
        "greedy_packing_lower_bound": quantiles(packing),
    }


def audit_shell(matrix: np.ndarray, cap: int, deficit: int, include_matched: bool,
                attempt_exact_packing: bool = True) -> dict[str, Any]:
    n = len(matrix)
    spins = projective_spins(n)
    energies = energy_vector(matrix, spins)
    indices, orientations = shell_points(energies, cap, deficit)
    selected = spins[indices]
    selected_cuts = np.stack(
        [selected[:, i] * selected[:, j] for i, j in edge_list(n)], axis=1
    ).astype(np.int8)
    signed_words = orientations[:, None] * selected_cuts
    signed_barycenter = np.mean(signed_words.astype(float), axis=0)
    geometry = pair_geometry(selected)
    compatible = geometry.pop("compatible")
    greedy = greedy_packing(compatible)
    result = {
        "deficit": deficit,
        "threshold": cap - deficit,
        "shell_size": int(len(indices)),
        "orientation_plus": int(np.sum(orientations == 1)),
        "orientation_minus": int(np.sum(orientations == -1)),
        **geometry,
        "packing_edge_abs_at_most_half": (
            certified_packing(compatible, greedy)
            if attempt_exact_packing
            else {
                "greedy_lower_bound": len(greedy),
                "certified_lower_bound": len(greedy),
                "certified_upper_bound": len(compatible),
                "exact": len(compatible) <= 1,
                "optimum": len(compatible) if len(compatible) <= 1 else None,
                "status": "CONTROL_GREEDY_ONLY",
            }
        ),
        # This auxiliary was added transparently after the frozen primary run
        # to compare the observed geometry with the already-proved signed
        # barycentre dichotomy in AO.21.  It is not a preregistered outcome.
        "posthoc_auxiliary_signed_balance": {
            "edge_l1_mean": float(np.mean(np.abs(signed_barycenter))),
            "independent_signed_overlap_mean": float(
                np.mean(signed_barycenter * signed_barycenter)
            ),
        },
    }
    if include_matched:
        result["cardinality_matched_projective_null"] = matched_controls(n, len(indices))
    return result


def row_matrix(row: dict[str, Any]) -> np.ndarray:
    matrix = np.asarray(row["matrix"], dtype=np.int8)
    if not (
        np.array_equal(matrix, matrix.T)
        and np.all(np.diag(matrix) == 0)
        and np.all(np.abs(matrix[np.triu_indices(len(matrix), 1)]) == 1)
    ):
        raise AssertionError("invalid matrix in source audit")
    return matrix


def deduplicate_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    unique: dict[str, dict[str, Any]] = {}
    for row in rows:
        unique.setdefault(matrix_hash(row_matrix(row)), row)
    return [unique[key] for key in sorted(unique)]


def source_note(row: dict[str, Any]) -> list[str]:
    notes: list[str] = []
    for source in row.get("sources", []):
        if isinstance(source, dict) and "file" in source:
            notes.append(str(source["file"]))
    return sorted(set(notes))


def selected_physical_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    # Orders 7--8 have authoritative complete orbit inventories.
    for n in (7, 8):
        for row in data["authoritative_orbit_inventory"][str(n)]["classes"]:
            selected.append({
                **row,
                "audit_stratum": "exact_authoritative_orbit_representative",
                "provenance_strength": "exhaustive orbit classification and exact cap",
            })
    # At larger orders retain all available byte-distinct repository witnesses.
    for row in deduplicate_rows(data["repository_exact_representatives"]):
        n = int(row["observables"]["n"])
        if 9 <= n <= 14:
            selected.append({
                **row,
                "audit_stratum": "exact_available_repository_witness",
                "provenance_strength": "exact cap equals certified M_n; not orbit-uniform",
            })
    for row in deduplicate_rows(data["repository_one_step_near_representatives"]):
        n = int(row["observables"]["n"])
        if 7 <= n <= 14:
            selected.append({
                **row,
                "audit_stratum": "one_step_near_repository_witness",
                "provenance_strength": "exact enumerated cap M_n+2; witness source varies",
            })
    supplements = (
        ("independently_generated_greedy_low_cap", "one_step_near_greedy_search_witness"),
        ("cap_constrained_adversarial_samples", "one_step_near_cap_walk_witness"),
    )
    for source_key, label in supplements:
        by_n: defaultdict[int, list[dict[str, Any]]] = defaultdict(list)
        for row in deduplicate_rows(data[source_key]):
            n = int(row["observables"]["n"])
            if 7 <= n <= 14 and int(row["observables"]["cap_delta"]) == 2:
                by_n[n].append(row)
        for n, rows in sorted(by_n.items()):
            for row in sorted(rows, key=lambda item: matrix_hash(row_matrix(item)))[:SUPPLEMENT_PER_STRATUM_ORDER]:
                selected.append({
                    **row,
                    "audit_stratum": label,
                    "provenance_strength": "heuristic discovery; exact enumerated cap M_n+2",
                })
    # A matrix may occur in several semantically distinct source strata.  Keep
    # the strongest/earliest stratum only to avoid duplicate outcome weight.
    final: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in selected:
        digest = matrix_hash(row_matrix(row))
        if digest not in seen:
            seen.add(digest)
            final.append(row)
    return final


def audit_physical(data: dict[str, Any]) -> list[dict[str, Any]]:
    exact_m = {int(k): int(v) for k, v in data["supplied_exact_M"].items()}
    records: list[dict[str, Any]] = []
    for row in selected_physical_rows(data):
        matrix = row_matrix(row)
        n = len(matrix)
        spins = projective_spins(n)
        energies = energy_vector(matrix, spins)
        cap = int(np.max(np.abs(energies)))
        expected_cap = int(row["observables"]["cap"])
        if cap != expected_cap:
            raise AssertionError("source cap did not reproduce")
        if cap not in (exact_m[n], exact_m[n] + 2):
            raise AssertionError("selected matrix is not exact or one-step near")
        deficits = sorted(set((0, 2, 4, 2 * math.floor(math.sqrt(n)))))
        records.append({
            "n": n,
            "matrix_sha256_upper": matrix_hash(matrix),
            "audit_stratum": row["audit_stratum"],
            "provenance_strength": row["provenance_strength"],
            "source_files": source_note(row),
            "cap": cap,
            "M_n": exact_m[n],
            "cap_delta": cap - exact_m[n],
            "shells": [audit_shell(matrix, cap, d, include_matched=True) for d in deficits],
        })
    return records


def summarize_numeric(values: list[float]) -> dict[str, float] | None:
    return quantiles(values)


def audit_random_controls(data: dict[str, Any]) -> dict[str, Any]:
    exact_m = {int(k): int(v) for k, v in data["supplied_exact_M"].items()}
    rng = np.random.default_rng(SEED)
    result: dict[str, Any] = {}
    for n in range(7, 15):
        per_deficit: defaultdict[int, list[dict[str, Any]]] = defaultdict(list)
        caps: list[int] = []
        hashes: list[str] = []
        for _ in range(RANDOM_SIGNINGS_PER_ORDER):
            signs = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=n * (n - 1) // 2)
            matrix = matrix_from_signs(n, signs)
            energies = energy_vector(matrix, projective_spins(n))
            cap = int(np.max(np.abs(energies)))
            caps.append(cap)
            hashes.append(matrix_hash(matrix))
            for deficit in sorted(set((0, 2, 4, 2 * math.floor(math.sqrt(n))))):
                per_deficit[deficit].append(
                    audit_shell(
                        matrix, cap, deficit, include_matched=False,
                        attempt_exact_packing=False,
                    )
                )
        summaries: dict[str, Any] = {}
        for deficit, rows in sorted(per_deficit.items()):
            def field(path: tuple[str, ...]) -> list[float]:
                values: list[float] = []
                for row in rows:
                    value: Any = row
                    for key in path:
                        value = value.get(key) if isinstance(value, dict) else None
                    if value is not None:
                        values.append(float(value))
                return values
            summaries[str(deficit)] = {
                "shell_size": summarize_numeric(field(("shell_size",))),
                "minimum_edge_abs_overlap": summarize_numeric(field(("edge_abs_overlap", "min"))),
                "minimum_vertex_projective_overlap": summarize_numeric(field(("vertex_projective_overlap", "min"))),
                "greedy_packing_lower_bound": summarize_numeric(
                    field(("packing_edge_abs_at_most_half", "greedy_lower_bound"))
                ),
                "exact_packing_certified_count": sum(
                    bool(row["packing_edge_abs_at_most_half"]["exact"]) for row in rows
                ),
                "diffuse_certificate_packing_at_least_3_count": sum(
                    int(row["packing_edge_abs_at_most_half"]["certified_lower_bound"]) >= 3
                    for row in rows
                ),
            }
        result[str(n)] = {
            "sample_count": RANDOM_SIGNINGS_PER_ORDER,
            "M_n_reference_only": exact_m[n],
            "caps": quantiles(caps),
            "matrix_hashes": hashes,
            "shell_summaries": summaries,
        }
    return result


def invariance_check(records: list[dict[str, Any]], data: dict[str, Any]) -> dict[str, Any]:
    source_by_hash = {
        matrix_hash(row_matrix(row)): row_matrix(row) for row in selected_physical_rows(data)
    }
    failures: list[str] = []
    for record in records[:16]:
        matrix = source_by_hash[record["matrix_sha256_upper"]]
        n = len(matrix)
        rng = stable_rng(SEED, "invariance", record["matrix_sha256_upper"])
        switch = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=n)
        permutation = rng.permutation(n)
        transformed = (switch[:, None] * matrix * switch[None, :])[np.ix_(permutation, permutation)]
        if rng.integers(2):
            transformed = -transformed
        cap = int(np.max(np.abs(energy_vector(matrix, projective_spins(n)))))
        transformed_cap = int(np.max(np.abs(energy_vector(transformed, projective_spins(n)))))
        if cap != transformed_cap:
            failures.append(f"cap:{record['matrix_sha256_upper']}")
            continue
        for shell in record["shells"]:
            first = audit_shell(
                matrix, cap, int(shell["deficit"]), include_matched=False,
                attempt_exact_packing=False,
            )
            second = audit_shell(
                transformed, cap, int(shell["deficit"]), include_matched=False,
                attempt_exact_packing=False,
            )
            keys = (
                "shell_size",
                "orientation_plus",
                "orientation_minus",
                "pair_count",
                "edge_abs_overlap",
                "vertex_projective_overlap",
            )
            # Global matrix sign swaps the two orientation counts, so compare sorted.
            if first["shell_size"] != second["shell_size"] or first["pair_count"] != second["pair_count"]:
                failures.append(f"size:{record['matrix_sha256_upper']}:{shell['deficit']}")
            if sorted((first["orientation_plus"], first["orientation_minus"])) != sorted(
                (second["orientation_plus"], second["orientation_minus"])
            ):
                failures.append(f"orientation:{record['matrix_sha256_upper']}:{shell['deficit']}")
            for key in keys[4:]:
                first_values = first[key]
                second_values = second[key]
                if first_values is None or second_values is None:
                    equal = first_values is second_values
                else:
                    equal = first_values.keys() == second_values.keys() and all(
                        math.isclose(
                            float(first_values[name]), float(second_values[name]),
                            rel_tol=0, abs_tol=1e-14,
                        )
                        for name in first_values
                    )
                if not equal:
                    failures.append(f"{key}:{record['matrix_sha256_upper']}:{shell['deficit']}")
    return {"passed": not failures, "failure_count": len(failures), "failures": failures}


def compact_summary(payload: dict[str, Any]) -> dict[str, Any]:
    records = payload["physical_nearmin_records"]
    by_order: dict[str, Any] = {}
    for n in range(7, 15):
        order: dict[str, Any] = {}
        for cap_delta, label in ((0, "exact"), (2, "one_step_near")):
            chosen = [row for row in records if row["n"] == n and row["cap_delta"] == cap_delta]
            if not chosen:
                continue
            shell_rows: dict[str, Any] = {}
            deficits = sorted({shell["deficit"] for row in chosen for shell in row["shells"]})
            for deficit in deficits:
                shells = [
                    next(shell for shell in row["shells"] if shell["deficit"] == deficit)
                    for row in chosen
                ]
                packing = [
                    shell["packing_edge_abs_at_most_half"]["certified_lower_bound"]
                    for shell in shells
                ]
                minimum_overlap = [
                    shell["edge_abs_overlap"]["min"]
                    for shell in shells if shell["edge_abs_overlap"] is not None
                ]
                matched_above = 0
                matched_below = 0
                for shell, value in zip(shells, packing):
                    null = shell["cardinality_matched_projective_null"][
                        "greedy_packing_lower_bound"
                    ]
                    matched_above += value > null["max"]
                    matched_below += value < null["min"]
                shell_rows[str(deficit)] = {
                    "record_count": len(shells),
                    "shell_size_range": [
                        min(shell["shell_size"] for shell in shells),
                        max(shell["shell_size"] for shell in shells),
                    ],
                    "minimum_edge_abs_overlap_range": (
                        [min(minimum_overlap), max(minimum_overlap)]
                        if minimum_overlap else None
                    ),
                    "packing_certified_lower_bound_range": [min(packing), max(packing)],
                    "exact_packing_certificate_count": sum(
                        shell["packing_edge_abs_at_most_half"]["exact"] for shell in shells
                    ),
                    "packing_at_least_three_count": sum(value >= 3 for value in packing),
                    "above_all_32_matched_controls_count": matched_above,
                    "below_all_32_matched_controls_count": matched_below,
                    "posthoc_signed_barycenter_l1_range": [
                        min(
                            shell["posthoc_auxiliary_signed_balance"]["edge_l1_mean"]
                            for shell in shells
                        ),
                        max(
                            shell["posthoc_auxiliary_signed_balance"]["edge_l1_mean"]
                            for shell in shells
                        ),
                    ],
                }
            order[label] = {
                "record_count": len(chosen),
                "stratum_counts": {
                    stratum: sum(row["audit_stratum"] == stratum for row in chosen)
                    for stratum in sorted({row["audit_stratum"] for row in chosen})
                },
                "shells": shell_rows,
            }
        by_order[str(n)] = order
    random_controls = {
        n: {
            "sample_count": row["sample_count"],
            "caps": row["caps"],
            "shell_summaries": row["shell_summaries"],
        }
        for n, row in payload["uniform_random_signing_controls"].items()
    }
    return {
        "schema": "nearmin-absolute-overlap-audit-compact-summary-v1",
        "status": payload["status"],
        "frozen_protocol": payload["frozen_protocol"],
        "detailed_output": "extremal_information/experiments/nearmin_absolute_overlap_results.json",
        "parameters": payload["parameters"],
        "classification": payload["classification"],
        "checks": payload["checks"],
        "physical_by_order": by_order,
        "uniform_random_signing_controls": random_controls,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    physical = audit_physical(data)
    payload = {
        "schema": "nearmin-absolute-overlap-audit-v1",
        "status": "FINITE AUDIT ONLY; NO ASYMPTOTIC CLAIM",
        "frozen_protocol": "extremal_information/experiments/nearmin_absolute_overlap_observable_freeze.md",
        "input": str(args.input.relative_to(ROOT)),
        "parameters": {
            "seed": SEED,
            "orders": [7, 8, 9, 10, 11, 12, 13, 14],
            "random_signings_per_order": RANDOM_SIGNINGS_PER_ORDER,
            "matched_subsets_per_physical_shell": MATCHED_SUBSETS_PER_SHELL,
            "supplement_per_stratum_order": SUPPLEMENT_PER_STRATUM_ORDER,
            "exact_packing_size_limit": EXACT_PACKING_SIZE_LIMIT,
            "exact_packing_seconds": EXACT_PACKING_SECONDS,
            "packing_rule": "pairwise absolute edge overlap <= 1/2",
        },
        "classification": {
            "exact": "exhaustively recomputed cap equals supplied certified M_n",
            "one_step_near": "exhaustively recomputed cap equals M_n+2",
            "search_warning": "greedy/cap-walk source is heuristic even when the saved cap is exact",
            "orbit_warning": "orders 9--14 are available witnesses, not an orbit-uniform sample",
        },
        "physical_nearmin_records": physical,
        "uniform_random_signing_controls": audit_random_controls(data),
        "checks": {
            "all_caps_classified": all(row["cap_delta"] in (0, 2) for row in physical),
            "overlap_identity": True,
            "switching_permutation_global_sign_invariance": invariance_check(physical, data),
        },
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    args.summary_output.write_text(
        json.dumps(compact_summary(payload), indent=2, sort_keys=True) + "\n"
    )
    print(f"wrote {args.output}")
    print(f"wrote {args.summary_output}")
    print(f"physical records: {len(physical)}")
    print(
        "shells:",
        sum(len(row["shells"]) for row in physical),
        "exact packing certificates:",
        sum(
            shell["packing_edge_abs_at_most_half"]["exact"]
            for row in physical for shell in row["shells"]
        ),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
