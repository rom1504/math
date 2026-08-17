#!/usr/bin/env python3
"""Exact finite audit of local edge-edit response roofs.

The observable definitions were frozen in
``nearmin_radial_response_observable_freeze.md`` before this script was run.
All landscape and edit-context enumerations are exact.  Set-cover entries are
labelled exact only when CP-SAT proves optimality; otherwise certified bounds
are retained.  Greedy metric packing/covering entries are only bounds.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "extremal_information/experiments/nearmin_blind_structural_results.json"
DEFAULT_OUTPUT = ROOT / "extremal_information/experiments/nearmin_radial_response_results.json"
SEED = 20260817
M_EXACT = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: 12, 10: 13, 11: 17,
           12: 18, 13: 20, 14: 23}


def matrix_hash(a: np.ndarray) -> str:
    return hashlib.sha256(a.astype(np.int8).tobytes()).hexdigest()


def edge_list(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(n), 2))


def projective_spins(n: int) -> np.ndarray:
    masks = np.arange(1 << (n - 1), dtype=np.uint32)
    bits = ((masks[:, None] >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int8)
    spins = np.ones((len(masks), n), dtype=np.int8)
    spins[:, 1:] = 1 - 2 * bits
    return spins


def edge_products(a: np.ndarray, spins: np.ndarray, edges: list[tuple[int, int]]) -> np.ndarray:
    """q[x,e]=a_e x_i x_j."""
    return np.stack(
        [a[i, j] * spins[:, i] * spins[:, j] for i, j in edges], axis=1
    ).astype(np.int8)


def contexts_through_radius(edge_count: int, maximum_radius: int) -> tuple[list[tuple[int, ...]], dict[int, int]]:
    contexts: list[tuple[int, ...]] = [()]
    ends = {0: 1}
    for r in range(1, maximum_radius + 1):
        contexts.extend(combinations(range(edge_count), r))
        ends[r] = len(contexts)
    return contexts, ends


def edit_sums(q: np.ndarray, contexts: list[tuple[int, ...]]) -> np.ndarray:
    chunks = [np.zeros((q.shape[0], 1), dtype=np.int16)]
    start = 1
    while start < len(contexts):
        size = len(contexts[start])
        end = start
        while end < len(contexts) and len(contexts[end]) == size:
            end += 1
        indices = np.asarray(contexts[start:end], dtype=np.int16)
        # Advanced indexing gives (spin, context, radius); radius <= 3.
        chunks.append(q[:, indices].sum(axis=2, dtype=np.int16))
        start = end
    return np.concatenate(chunks, axis=1)


def entropy_bits(weights: np.ndarray) -> float:
    positive = weights[weights > 0]
    if not len(positive):
        return 0.0
    return float(-np.sum(positive * np.log2(positive)))


def reduce_requirements(requirements: Iterable[tuple[int, ...]]) -> list[tuple[int, ...]]:
    """Drop duplicate and superset constraints from a monotone set cover."""
    unique = sorted(set(requirements), key=lambda row: (len(row), row))
    kept: list[tuple[int, ...]] = []
    kept_sets: list[frozenset[int]] = []
    for row in unique:
        s = frozenset(row)
        if any(old.issubset(s) for old in kept_sets):
            continue
        kept.append(row)
        kept_sets.append(s)
    return kept


def greedy_set_cover(requirements: list[tuple[int, ...]], candidate_count: int) -> list[int]:
    remaining = set(range(len(requirements)))
    incidence: list[set[int]] = [set() for _ in range(candidate_count)]
    for c, row in enumerate(requirements):
        for z in row:
            incidence[z].add(c)
    chosen: list[int] = []
    while remaining:
        z = max(range(candidate_count), key=lambda i: (len(incidence[i] & remaining), -i))
        hit = incidence[z] & remaining
        if not hit:
            raise AssertionError("uncovered response context")
        chosen.append(z)
        remaining.difference_update(hit)
    return chosen


def solve_cover(requirements: list[tuple[int, ...]], candidate_count: int,
                time_limit: float) -> dict[str, object]:
    requirements = reduce_requirements(requirements)
    greedy = greedy_set_cover(requirements, candidate_count)
    model = cp_model.CpModel()
    variables = [model.NewBoolVar(f"z{i}") for i in range(candidate_count)]
    for row in requirements:
        model.AddBoolOr([variables[i] for i in row])
    objective = sum(variables)
    model.Minimize(objective)
    # A concrete feasible cover makes the upper bound independent of solver status.
    for i in range(candidate_count):
        model.AddHint(variables[i], int(i in set(greedy)))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
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
    feasible_value = len(greedy)
    chosen: list[int] | None = None
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        chosen = [i for i, var in enumerate(variables) if solver.Value(var)]
        feasible_value = min(feasible_value, len(chosen))
    lower = int(math.ceil(solver.BestObjectiveBound() - 1e-7)) if status != cp_model.MODEL_INVALID else 1
    lower = max(1, min(lower, feasible_value))
    return {
        "status": names.get(status, str(status)),
        "lower_bound": lower,
        "upper_bound": feasible_value,
        "exact": bool(status == cp_model.OPTIMAL),
        "optimum": feasible_value if status == cp_model.OPTIMAL else None,
        "reduced_constraint_count": len(requirements),
        "greedy_upper_bound": len(greedy),
        "selected_candidate_indices": chosen if status == cp_model.OPTIMAL else None,
    }


def contextual_metric(t: np.ndarray, radius: int) -> np.ndarray:
    """Exact d_r between affine edit-response functions.

    t[z,e]=a_e z_e is signed.  For w=t_z-t_z', editing an edge subtracts
    2w_e.  Since w_e is in {-2,0,2}, the two extrema over |F|<=r depend only
    on the oriented mismatch counts.
    """
    positive = (t == 1).astype(np.int16)
    negative = (t == -1).astype(np.int16)
    pos_count = positive @ negative.T
    neg_count = pos_count.T
    base = t.sum(axis=1, dtype=np.int16)
    base_diff = base[:, None] - base[None, :]
    maximum = base_diff + 4 * np.minimum(neg_count, radius)
    minimum = base_diff - 4 * np.minimum(pos_count, radius)
    return np.maximum(np.abs(maximum), np.abs(minimum)).astype(np.int16)


def greedy_metric_bounds(distance: np.ndarray, tolerance: int) -> dict[str, int]:
    count = len(distance)
    # Deterministic maximal separated set: inspect high-eccentricity points first.
    order = sorted(range(count), key=lambda i: (-int(np.sum(distance[i] > tolerance)), i))
    packing: list[int] = []
    for i in order:
        if all(distance[i, j] > tolerance for j in packing):
            packing.append(i)

    uncovered = set(range(count))
    cover: list[int] = []
    balls = [set(np.flatnonzero(distance[i] <= tolerance).tolist()) for i in range(count)]
    while uncovered:
        i = max(range(count), key=lambda j: (len(balls[j] & uncovered), -j))
        cover.append(i)
        uncovered.difference_update(balls[i])
    return {"greedy_packing_lower_bound": len(packing), "greedy_cover_upper_bound": len(cover)}


def radius_audit(q: np.ndarray, adjusted: np.ndarray, cap: int, eta: int,
                 context_count: int, radius: int, solver_seconds: float) -> dict[str, object]:
    adj = adjusted[:, :context_count]
    response_abs = np.max(np.abs(adj), axis=0)
    spin_ties = np.abs(adj) == response_abs[None, :]

    exposed_ids: set[int] = set()
    tie_id_rows: list[tuple[int, ...]] = []
    probability_mass: Counter[int] = Counter()
    tie_logs = 0.0
    tie_counts: list[int] = []
    for c in range(context_count):
        xs = np.flatnonzero(spin_ties[:, c])
        ids: list[int] = []
        for x in xs:
            value = int(adj[x, c])
            if value >= 0:
                ids.append(2 * int(x) + 1)  # plus orientation
            if value <= 0:
                ids.append(2 * int(x))      # minus orientation
        row = tuple(sorted(ids))
        tie_id_rows.append(row)
        exposed_ids.update(row)
        tie_counts.append(len(row))
        tie_logs += math.log2(len(row))
        mass = 1.0 / (context_count * len(row))
        for z in row:
            probability_mass[z] += mass

    exposed = sorted(exposed_ids)
    local_index = {z: i for i, z in enumerate(exposed)}
    xs = np.asarray([z // 2 for z in exposed], dtype=np.int32)
    signs = np.asarray([1 if z % 2 else -1 for z in exposed], dtype=np.int16)
    candidate_values = signs[:, None] * adj[xs, :]

    covers: dict[str, object] = {}
    for delta in (0, 2, 4):
        good = candidate_values >= response_abs[None, :] - delta
        requirements = [tuple(np.flatnonzero(good[:, c]).tolist()) for c in range(context_count)]
        covers[str(delta)] = solve_cover(requirements, len(exposed), solver_seconds)

    weights = np.asarray([probability_mass[z] for z in exposed], dtype=float)
    weights /= weights.sum()
    hz = entropy_bits(weights)
    conditional = tie_logs / context_count

    # The common RS shell is counted in the full augmented family, including
    # the rarely relevant opposite orientation when the threshold is broad.
    base = q.sum(axis=1, dtype=np.int16)
    base_scores = np.concatenate((-base, base))
    shell_threshold = eta + 2 * radius
    shell_size = int(np.count_nonzero(cap - base_scores <= shell_threshold))

    t = signs[:, None] * q[xs, :]
    metric = contextual_metric(t, radius)
    metric_summary: dict[str, object] = {"diameter": int(np.max(metric))}
    for tolerance in (2, 4):
        metric_summary[str(tolerance)] = greedy_metric_bounds(metric, tolerance)

    histogram = Counter(int(v - cap) for v in response_abs)
    response_probs = np.asarray(list(histogram.values()), dtype=float) / context_count
    response_entropy = entropy_bits(response_probs)
    return {
        "radius": radius,
        "context_count": context_count,
        "rs_shell_width": shell_threshold,
        "rs_shell_size": shell_size,
        "exposed_witness_count": len(exposed),
        "exposed_fraction_of_rs_shell": len(exposed) / shell_size,
        "response_cover": covers,
        "optimizer_information": {
            "H_Z_bits": hz,
            "H_Z_given_F_bits": conditional,
            "I_F_Z_bits": hz - conditional,
            "effective_support": 2.0 ** hz,
            "mean_tie_count": float(np.mean(tie_counts)),
            "max_tie_count": int(max(tie_counts)),
        },
        "response_distribution": {
            "response_minus_base_cap_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "entropy_bits": response_entropy,
            "normalized_entropy_by_log_contexts": (
                response_entropy / math.log2(context_count) if context_count > 1 else 0.0
            ),
        },
        "contextual_affine_metric": metric_summary,
        "normalized": {
            # Per-vertex normalizations are injected by the caller, which has
            # the vertex count; this helper sees only the edge table.
            "log2_rs_shell": math.log2(shell_size),
            "log2_exposed": math.log2(len(exposed)),
            "optimizer_mutual_information_bits": hz - conditional,
        },
    }


def cap_of_q(q: np.ndarray) -> int:
    return int(np.max(np.abs(q.sum(axis=1, dtype=np.int16))))


def audit_matrix(a: np.ndarray, label: str, source: str, solver_seconds: float) -> dict[str, object]:
    n = len(a)
    edges = edge_list(n)
    spins = projective_spins(n)
    q = edge_products(a, spins, edges)
    cap = cap_of_q(q)
    contexts, ends = contexts_through_radius(len(edges), 3)
    sums = edit_sums(q, contexts)
    adjusted = q.sum(axis=1, dtype=np.int16)[:, None] - 2 * sums
    radii = [radius_audit(q, adjusted, cap, cap - M_EXACT[n], ends[r], r, solver_seconds)
             for r in (1, 2, 3)]
    for row in radii:
        row["normalized"]["log2_rs_shell_over_n"] = row["normalized"]["log2_rs_shell"] / n
        row["normalized"]["log2_exposed_over_n"] = row["normalized"]["log2_exposed"] / n
        row["normalized"]["optimizer_mutual_information_over_n"] = (
            row["normalized"]["optimizer_mutual_information_bits"] / n
        )
        for delta, cover in row["response_cover"].items():
            cover["log2_upper_bound_over_n"] = math.log2(cover["upper_bound"]) / n
            cover["log2_lower_bound_over_n"] = math.log2(cover["lower_bound"]) / n
    return {
        "n": n,
        "edge_count": len(edges),
        "label": label,
        "source": source,
        "matrix_sha256": matrix_hash(a),
        "matrix": a.astype(int).tolist(),
        "cap": cap,
        "M_n": M_EXACT[n],
        "cap_excess": cap - M_EXACT[n],
        "radii": radii,
    }


def unique_sorted(records: list[dict[str, object]], n: int, cap: int | None = None) -> list[dict[str, object]]:
    by_hash: dict[str, dict[str, object]] = {}
    for row in records:
        a = np.asarray(row["matrix"], dtype=np.int8)
        if len(a) != n:
            continue
        observed_cap = int(row.get("observables", {}).get("cap", -999))
        if cap is not None and observed_cap != cap:
            continue
        by_hash.setdefault(matrix_hash(a), row)
    return [by_hash[key] for key in sorted(by_hash)]


def select_population(payload: dict[str, object], min_n: int, max_n: int,
                      per_class: int) -> list[tuple[np.ndarray, str, str]]:
    selected: list[tuple[np.ndarray, str, str]] = []
    rng = np.random.default_rng(SEED)
    for n in range(min_n, max_n + 1):
        exact = unique_sorted(payload["repository_exact_representatives"], n, M_EXACT[n])[:per_class]
        near = unique_sorted(payload["one_edge_low_cap_neighborhood"], n, M_EXACT[n] + 2)[:per_class]
        for row in exact:
            selected.append((np.asarray(row["matrix"], dtype=np.int8), "exact minimizer",
                             json.dumps(row.get("sources", "repository"), sort_keys=True)))
        for row in near:
            selected.append((np.asarray(row["matrix"], dtype=np.int8), "one-step near-minimizer",
                             f"blind audit one-edge neighborhood; edge={row.get('flipped_edge')}"))
        for draw in range(per_class):
            upper = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=(n, n))
            a = np.triu(upper, 1)
            a = a + a.T
            selected.append((a, "uniform random", f"seed={SEED}; order-local draw={draw}"))
    return selected


def quantiles(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"count": 0, "min": None, "median": None, "max": None}
    arr = np.asarray(values, dtype=float)
    return {"count": len(values), "min": float(np.min(arr)),
            "median": float(np.median(arr)), "max": float(np.max(arr))}


def summarize(records: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in sorted({int(x["n"]) for x in records}):
        for label in ("exact minimizer", "one-step near-minimizer", "uniform random"):
            group = [x for x in records if x["n"] == n and x["label"] == label]
            for radius in (1, 2, 3):
                rr = [x["radii"][radius - 1] for x in group]
                row: dict[str, object] = {"n": n, "label": label, "radius": radius,
                                           "matrix_count": len(group)}
                for name, getter in {
                    "cap_excess": lambda x, y: x["cap_excess"],
                    "log2_shell_over_n": lambda x, y: y["normalized"]["log2_rs_shell_over_n"],
                    "log2_exposed_over_n": lambda x, y: y["normalized"]["log2_exposed_over_n"],
                    "mutual_information_over_n": lambda x, y: y["normalized"]["optimizer_mutual_information_over_n"],
                    "exact_cover0_log2_over_n": lambda x, y: y["response_cover"]["0"]["log2_upper_bound_over_n"],
                    "cover2_log2_over_n": lambda x, y: y["response_cover"]["2"]["log2_upper_bound_over_n"],
                    "cover4_log2_over_n": lambda x, y: y["response_cover"]["4"]["log2_upper_bound_over_n"],
                    "exposed_fraction_shell": lambda x, y: y["exposed_fraction_of_rs_shell"],
                    "response_entropy_normalized": lambda x, y: y["response_distribution"]["normalized_entropy_by_log_contexts"],
                }.items():
                    row[name] = quantiles([float(getter(x, y)) for x, y in zip(group, rr)])
                rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--min-order", type=int, default=7)
    parser.add_argument("--max-order", type=int, default=11)
    parser.add_argument("--per-class", type=int, default=3)
    parser.add_argument("--solver-seconds", type=float, default=4.0)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    population = select_population(payload, args.min_order, args.max_order, args.per_class)
    records = []
    for index, (a, label, source) in enumerate(population, 1):
        print(f"[{index}/{len(population)}] n={len(a)} {label} {matrix_hash(a)[:10]}", flush=True)
        records.append(audit_matrix(a, label, source, args.solver_seconds))

    output = {
        "schema": "nearmin-radial-response-roof-audit-v1",
        "status": "EXACT FINITE ENUMERATION; COVER OPTIMALITY LABELLED PER ENTRY",
        "observable_freeze": "extremal_information/experiments/nearmin_radial_response_observable_freeze.md",
        "seed": SEED,
        "parameters": vars(args) | {"input": str(args.input), "output": str(args.output)},
        "records": records,
        "summary": summarize(records),
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(records)}, indent=2))


if __name__ == "__main__":
    main()
