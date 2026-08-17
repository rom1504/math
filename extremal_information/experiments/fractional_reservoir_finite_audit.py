#!/usr/bin/env python3
"""Finite LP audit of simultaneous fractional reservoirs.

The selection and validation protocol was frozen in
``fractional_reservoir_finite_protocol.md`` before this program was run.
Every claim emitted here is finite numerical evidence only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import scipy
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT / "extremal_information" / "experiments"
    / "nearmin_blind_structural_results.json"
)
DEFAULT_OUTPUT = (
    ROOT / "extremal_information" / "experiments"
    / "fractional_reservoir_finite_results.json"
)
PROTOCOL = (
    ROOT / "extremal_information" / "experiments"
    / "fractional_reservoir_finite_protocol.md"
)
SUPPORT_TOL = 1e-9
VALIDATION_TOL = 1e-7


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_matrix(matrix: np.ndarray) -> str:
    payload = json.dumps(matrix.astype(int).tolist(), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def signing_data(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(matrix)
    edges = tuple(combinations(range(n), 2))
    rows = np.arange(1 << (n - 1), dtype=np.uint32)[:, None]
    bits = ((rows >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int8)
    spins = np.concatenate(
        [np.ones((len(rows), 1), dtype=np.int8), 1 - 2 * bits], axis=1
    )
    products = np.column_stack(
        [spins[:, i] * spins[:, j] for i, j in edges]
    ).astype(np.int8)
    signs = np.asarray([matrix[i, j] for i, j in edges], dtype=np.int8)
    return products, signs, products @ signs


def normalize_entry(entry: dict[str, Any], stratum: str) -> dict[str, Any]:
    matrix = np.asarray(entry["matrix"], dtype=np.int8)
    observables = entry["observables"]
    matrix_hash = sha256_matrix(matrix)
    expected_hash = observables.get("matrix_sha256")
    if expected_hash is not None and matrix_hash != expected_hash:
        raise ValueError(f"matrix hash mismatch for {entry.get('label')}")
    return {
        "stratum": stratum,
        "label": entry.get("label", stratum),
        "matrix": matrix,
        "matrix_sha256": matrix_hash,
        "n": int(observables["n"]),
        "cap": int(observables["cap"]),
        "cap_delta": int(observables["cap_delta"]),
        "sources": entry.get("sources", []),
    }


def first_by_stratum(
    entries: Iterable[dict[str, Any]], stratum: str
) -> list[dict[str, Any]]:
    normalized = [normalize_entry(entry, stratum) for entry in entries]
    groups: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for entry in normalized:
        groups[(entry["n"], entry["cap_delta"])].append(entry)
    chosen = []
    for key in sorted(groups):
        unique = {entry["matrix_sha256"]: entry for entry in groups[key]}
        chosen.append(unique[sorted(unique)[0]])
    return chosen


def select_corpus(data: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    candidates: list[dict[str, Any]] = []

    # Exact: all certified orbit classes through order 8.
    for n_text in sorted(data["authoritative_orbit_inventory"], key=int):
        n = int(n_text)
        if n > 8:
            continue
        entries = [
            normalize_entry(entry, "exact")
            for entry in data["authoritative_orbit_inventory"][n_text]["classes"]
        ]
        candidates.extend(sorted(entries, key=lambda entry: entry["matrix_sha256"]))

    # Exact: at most two SHA-first repository representatives per order 9--14.
    exact_by_n: dict[int, dict[str, dict[str, Any]]] = defaultdict(dict)
    for raw in data["repository_exact_representatives"]:
        entry = normalize_entry(raw, "exact")
        if entry["n"] >= 9:
            exact_by_n[entry["n"]][entry["matrix_sha256"]] = entry
    for n in sorted(exact_by_n):
        for matrix_hash in sorted(exact_by_n[n])[:2]:
            candidates.append(exact_by_n[n][matrix_hash])

    # Every distinct repository one-step-near representative.
    near = {
        entry["matrix_sha256"]: entry
        for entry in (
            normalize_entry(raw, "one_step_near")
            for raw in data["repository_one_step_near_representatives"]
        )
    }
    candidates.extend(near[matrix_hash] for matrix_hash in sorted(near))

    # One SHA-first item per (order, cap delta) for the remaining sources.
    candidates.extend(
        first_by_stratum(
            data["independently_generated_greedy_low_cap"], "heuristic_low_cap"
        )
    )
    candidates.extend(
        first_by_stratum(data["random_draws_that_are_low_cap"], "random_low_cap")
    )
    candidates.extend(
        first_by_stratum(
            data["cyclic_distance_low_cap_controls"], "structured_control"
        )
    )
    controls = {
        entry["matrix_sha256"]: entry
        for entry in (
            normalize_entry(raw, "structured_control")
            for raw in data["control_extremes"]
        )
    }
    candidates.extend(controls[matrix_hash] for matrix_hash in sorted(controls))

    # Global priority-order deduplication.
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    duplicate_count = 0
    for entry in candidates:
        if entry["matrix_sha256"] in seen:
            duplicate_count += 1
            continue
        seen.add(entry["matrix_sha256"])
        selected.append(entry)

    inventory: dict[str, Any] = {
        "candidate_count_before_global_deduplication": len(candidates),
        "duplicate_count": duplicate_count,
        "selected_count": len(selected),
        "by_stratum": {},
        "by_order": {},
        "by_cap_delta": {},
    }
    for field in ("stratum", "n", "cap_delta"):
        counts: dict[str, int] = defaultdict(int)
        for entry in selected:
            counts[str(entry[field])] += 1
        inventory[f"by_{field}" if field != "n" else "by_order"] = dict(
            sorted(counts.items(), key=lambda item: item[0])
        )
    return selected, inventory


def quantiles(values: Iterable[float]) -> dict[str, float] | None:
    array = np.asarray(list(values), dtype=float)
    if not len(array):
        return None
    q = np.quantile(array, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1])
    return {
        "min": float(q[0]),
        "q10": float(q[1]),
        "q25": float(q[2]),
        "median": float(q[3]),
        "q75": float(q[4]),
        "q90": float(q[5]),
        "max": float(q[6]),
        "mean": float(np.mean(array)),
    }


def solve_shell(entry: dict[str, Any], shell_type: str) -> dict[str, Any]:
    matrix = entry["matrix"]
    products, signs, energies = signing_data(matrix)
    cap = int(np.max(np.abs(energies)))
    if cap != entry["cap"]:
        raise ValueError(f"recomputed cap mismatch for {entry['matrix_sha256']}")
    if shell_type == "active":
        threshold = cap
        mask = np.abs(energies) == cap
    elif shell_type == "deficit_2":
        threshold = cap - 2
        mask = np.abs(energies) >= threshold
    else:
        raise ValueError(shell_type)
    if threshold <= 0:
        raise ValueError(f"nonpositive threshold {threshold}")

    shell_products = products[mask]
    shell_energies = energies[mask]
    orientation = np.sign(shell_energies).astype(np.int8)
    if np.any(orientation == 0):
        raise ValueError("zero energy entered positive shell")
    oriented_words = shell_products * orientation[:, None]
    response = oriented_words * signs[None, :]
    # Defensive deduplication; cut words should already be distinct here.
    response = np.unique(response, axis=0)
    shell_size = len(response)
    edge_count = len(signs)

    c = np.ones(edge_count, dtype=float)
    a_ub = -response.astype(float)
    b_ub = -threshold * np.ones(shell_size, dtype=float)
    bounds = [(0.0, 1.0)] * edge_count
    options = {
        "primal_feasibility_tolerance": 1e-9,
        "dual_feasibility_tolerance": 1e-9,
    }
    ds = linprog(c, A_ub=a_ub, b_ub=b_ub, bounds=bounds, method="highs-ds", options=options)
    ipm = linprog(c, A_ub=a_ub, b_ub=b_ub, bounds=bounds, method="highs-ipm", options=options)
    if not ds.success or not ipm.success:
        raise RuntimeError(
            f"LP failure {entry['matrix_sha256']} {shell_type}: "
            f"DS={ds.message}; IPM={ipm.message}"
        )

    w = np.asarray(ds.x, dtype=float)
    weighted_responses = response @ w
    objective = float(np.sum(w))
    objective_cross_error = abs(float(ds.fun) - float(ipm.fun))
    constraint_violation = max(0.0, float(threshold - np.min(weighted_responses)))
    lower_violation = max(0.0, float(-np.min(w)))
    upper_violation = max(0.0, float(np.max(w) - 1.0))

    # For min c*x subject to A*x<=b and 0<=x<=1, HiGHS marginals
    # reconstruct b*y + upper*marginal (lower bounds contribute zero).
    dual_objective = float(
        np.dot(b_ub, ds.ineqlin.marginals) + np.sum(ds.upper.marginals)
    )
    duality_gap = abs(objective - dual_objective)
    scale = max(1.0, objective)
    accepted = bool(
        objective_cross_error <= VALIDATION_TOL * scale
        and constraint_violation <= VALIDATION_TOL
        and lower_violation <= VALIDATION_TOL
        and upper_violation <= VALIDATION_TOL
        and duality_gap <= VALIDATION_TOL * scale
    )
    if not accepted:
        raise RuntimeError(
            f"validation failure {entry['matrix_sha256']} {shell_type}: "
            f"cross={objective_cross_error}, primal={constraint_violation}, "
            f"bounds={lower_violation, upper_violation}, gap={duality_gap}"
        )

    positive = w > SUPPORT_TOL
    unit = w >= 1.0 - SUPPORT_TOL
    fractional = positive & ~unit
    square_mass = float(np.dot(w, w))
    probabilities = w[positive] / objective
    entropy_effective_support = float(
        math.exp(-np.dot(probabilities, np.log(probabilities)))
    )
    common_correct = int(np.sum(np.all(response == 1, axis=0)))
    pattern_count = int(len(np.unique(response.T, axis=0)))
    column_sums = response.astype(np.int64).sum(axis=0)
    # Uniformly averaging all shell constraints is an exact finite dual
    # certificate for W>=E when every averaged edge response is <=m/E.
    uniform_scaled_excess = edge_count * column_sums - shell_size * threshold
    uniform_full_mass_certificate = bool(np.all(uniform_scaled_excess <= 0))
    full_edge_mass = bool(abs(objective - edge_count) <= VALIDATION_TOL)
    if uniform_full_mass_certificate and not full_edge_mass:
        raise RuntimeError("exact uniform dual certificate contradicts LP objective")

    return {
        "stratum": entry["stratum"],
        "label": entry["label"],
        "matrix_sha256": entry["matrix_sha256"],
        "n": entry["n"],
        "edge_count": edge_count,
        "cap": cap,
        "cap_delta": entry["cap_delta"],
        "shell_type": shell_type,
        "threshold": threshold,
        "shell_size": shell_size,
        "shell_density": shell_size / (1 << (entry["n"] - 1)),
        "edge_phase_pattern_count": pattern_count,
        "full_response_min": int(np.min(response.sum(axis=1))),
        "full_response_max": int(np.max(response.sum(axis=1))),
        "common_correct_edge_count": common_correct,
        "common_correct_mass_over_threshold": common_correct / threshold,
        "uniform_constraint_average_full_mass_certificate": uniform_full_mass_certificate,
        "uniform_average_scaled_excess_min": int(np.min(uniform_scaled_excess)),
        "uniform_average_scaled_excess_max": int(np.max(uniform_scaled_excess)),
        "objective_weight": objective,
        "C_inst": objective / threshold,
        "weight_over_edges": objective / edge_count,
        "full_edge_mass_solution": full_edge_mass,
        "support_count": int(np.sum(positive)),
        "support_fraction": float(np.mean(positive)),
        "unit_weight_count": int(np.sum(unit)),
        "fractional_weight_count": int(np.sum(fractional)),
        "max_weight": float(np.max(w)),
        "max_normalized_atom": float(np.max(w) / objective),
        "inverse_herfindahl_support": objective * objective / square_mass,
        "inverse_herfindahl_support_fraction": objective * objective / (square_mass * edge_count),
        "entropy_effective_support": entropy_effective_support,
        "entropy_effective_support_fraction": entropy_effective_support / edge_count,
        "weighted_response_min": float(np.min(weighted_responses)),
        "weighted_response_max": float(np.max(weighted_responses)),
        "active_constraint_count": int(
            np.sum(np.abs(weighted_responses - threshold) <= 1e-7)
        ),
        "validation": {
            "accepted": accepted,
            "dual_simplex_status": int(ds.status),
            "interior_point_status": int(ipm.status),
            "objective_cross_solver_error": objective_cross_error,
            "constraint_violation": constraint_violation,
            "lower_bound_violation": lower_violation,
            "upper_bound_violation": upper_violation,
            "dual_objective": dual_objective,
            "duality_gap": duality_gap,
        },
    }


SUMMARY_FIELDS = (
    "C_inst",
    "shell_size",
    "shell_density",
    "edge_phase_pattern_count",
    "common_correct_mass_over_threshold",
    "full_edge_mass_solution",
    "uniform_constraint_average_full_mass_certificate",
    "weight_over_edges",
    "support_fraction",
    "fractional_weight_count",
    "max_normalized_atom",
    "inverse_herfindahl_support_fraction",
    "entropy_effective_support_fraction",
)


def summarize(records: list[dict[str, Any]], keys: tuple[str, ...]) -> dict[str, Any]:
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        groups[tuple(str(record[key]) for key in keys)].append(record)
    output: dict[str, Any] = {}
    for group, items in sorted(groups.items()):
        label = " | ".join(f"{key}={value}" for key, value in zip(keys, group))
        output[label] = {
            "count": len(items),
            **{
                field: quantiles(item[field] for item in items)
                for field in SUMMARY_FIELDS
            },
        }
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = json.loads(args.input.read_text())
    selected, inventory = select_corpus(data)
    records = [
        solve_shell(entry, shell_type)
        for entry in selected
        for shell_type in ("active", "deficit_2")
    ]
    result = {
        "schema": "fractional-reservoir-finite-audit-v1",
        "status": "FINITE NUMERICAL LP EVIDENCE; NOT A THEOREM OR ASYMPTOTIC CLAIM",
        "objective": (
            "minimize sum_e w_e subject to all positively oriented shell "
            "responses sum_e w_e a_e z_e >= threshold and 0<=w_e<=1"
        ),
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "protocol_sha256": sha256_file(PROTOCOL),
        "input": str(args.input.relative_to(ROOT)),
        "input_sha256": sha256_file(args.input),
        "script": str(Path(__file__).resolve().relative_to(ROOT)),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "solver_methods": ["highs-ds", "highs-ipm"],
            "support_reporting_tolerance": SUPPORT_TOL,
            "validation_tolerance": VALIDATION_TOL,
        },
        "inventory": inventory,
        "all_records_validated": all(
            record["validation"]["accepted"] for record in records
        ),
        "record_count": len(records),
        "records": records,
        "summaries": {
            "by_shell": summarize(records, ("shell_type",)),
            "by_stratum_and_shell": summarize(records, ("stratum", "shell_type")),
            "by_order_and_shell": summarize(records, ("n", "shell_type")),
            "by_cap_delta_and_shell": summarize(records, ("cap_delta", "shell_type")),
            "by_stratum_cap_delta_and_shell": summarize(
                records, ("stratum", "cap_delta", "shell_type")
            ),
        },
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "selected_matrices": len(selected),
                "records": len(records),
                "validated": result["all_records_validated"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
