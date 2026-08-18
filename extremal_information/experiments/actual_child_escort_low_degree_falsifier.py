#!/usr/bin/env python3
"""Complete-cube low-degree cavity audit for actual thermal minimizers.

For the certified minimizing children at beta=4 and balanced/comparable
splits through total order ten, this program computes the exact finite bridge
cube (with floating transcendental evaluation), the negative escort q with
lambda=1, and the *best* degree-at-most-K polynomial approximation to each
edge-cavity response in weighted L2(q).

Child/signing and bridge enumeration are complete.  Thermal comparisons use
high-precision evaluations of exact energy histograms.  The final least-
squares residuals are numerical, so the output is finite evidence rather
than an asymptotic or interval-certified theorem.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_orbit_posterior_quotient as orbit  # noqa: E402


BETA = 4.0
LAMBDA = 1.0
DEGREES = (1, 3)
PLANS = ((6, 3, 3), (8, 4, 4), (10, 3, 7))


def fwht(values: np.ndarray) -> np.ndarray:
    """Unnormalized floating Walsh transform."""

    result = np.asarray(values, dtype=np.float64).copy()
    width = 1
    while width < len(result):
        view = result.reshape(-1, 2 * width)
        left = view[:, :width].copy()
        right = view[:, width:].copy()
        view[:, :width] = left + right
        view[:, width:] = left - right
        width *= 2
    return result


def vertex_orbits(matrix: np.ndarray) -> list[list[int]]:
    """Vertex orbits of the projective signed-automorphism group."""

    group = orbit.signed_automorphisms(matrix.astype(np.int16))
    unseen = set(range(len(matrix)))
    answer: list[list[int]] = []
    while unseen:
        representative = min(unseen)
        cell = sorted({permutation[representative] for permutation, _ in group})
        unseen.difference_update(cell)
        answer.append(cell)
    return answer


def child_record(
    order: int, total_order: int, cache: dict[tuple[int, int], tuple[np.ndarray, dict]]
) -> tuple[np.ndarray, dict]:
    key = (order, total_order)
    if key in cache:
        return cache[key]
    space = exact.build_signing_space(order)
    classes, certificate = exact.thermal_minimizer_classes(
        space, format(BETA, ".12g"), total_order
    )
    if len(classes) != 1:
        raise AssertionError(
            f"expected one thermal-minimizer class at order {order}, got {len(classes)}"
        )
    row = classes[0]
    matrix = np.asarray(row["representative_matrix"], dtype=np.int8)
    winner = space.absolute_histograms[int(row["representative_mask"])]
    mp.iv.dps = 80
    interval_t = mp.iv.mpf(format(BETA, ".12g")) / mp.iv.sqrt(
        mp.iv.mpf(total_order)
    )
    interval_cosh = [
        (mp.iv.exp(interval_t * value) + mp.iv.exp(-interval_t * value)) / 2
        for value in range(space.unique_histograms.shape[1])
    ]
    strict_differences = []
    for histogram in space.unique_histograms:
        if np.array_equal(histogram, winner):
            continue
        difference = sum(
            int(other - chosen) * value
            for other, chosen, value in zip(histogram, winner, interval_cosh)
        ) / int(len(space.spins))
        if not difference.a > 0:
            raise AssertionError(
                f"interval child-pressure comparison did not separate order {order}"
            )
        strict_differences.append(difference)
    minimum_interval = (
        min(strict_differences, key=lambda value: value.a)
        if strict_differences
        else None
    )
    result = {
        "order": order,
        "class_id": row["class_id"],
        "representative_sha256": row["representative_sha256"],
        "root_gauged_member_count": row["root_gauged_member_count"],
        "thermal_selection": certificate,
        "interval_pressure_certificate": {
            "interval_decimal_digits": mp.iv.dps,
            "competitor_histogram_count": len(strict_differences),
            "all_competitor_partition_sums_strictly_larger": True,
            "minimum_partition_sum_gap_interval": (
                None if minimum_interval is None else str(minimum_interval)
            ),
        },
        "signed_automorphism_group_order_projective": len(
            orbit.signed_automorphisms(matrix.astype(np.int16))
        ),
        "vertex_orbits": vertex_orbits(matrix),
    }
    cache[key] = (matrix, result)
    return cache[key]


def cavity_table(
    pressure: np.ndarray, edge: int, raw_t: float
) -> np.ndarray:
    """Exact flip-identity cavity table, numerically evaluated."""

    size = len(pressure)
    indices = np.arange(size, dtype=np.uint64)
    bit = 1.0 - 2.0 * (
        ((indices >> np.uint64(edge)) & np.uint64(1)).astype(np.float64)
    )
    flipped = indices ^ np.uint64(1 << edge)
    response = (
        bit
        * np.tanh(0.5 * (pressure - pressure[flipped]))
        / math.tanh(raw_t)
    )
    if np.max(np.abs(response - response[flipped])) > 3e-12:
        raise AssertionError("deleted-edge response depends on its deleted bit")
    if np.max(np.abs(response)) > 1.0 + 3e-10:
        raise AssertionError("cavity response left [-1,1]")
    return response


def best_weighted_degree(
    moment: np.ndarray,
    response_moment: np.ndarray,
    response_energy: float,
    popcount: np.ndarray,
    edge: int,
    degree: int,
) -> dict:
    """Solve the exact finite weighted normal equations in double precision."""

    indices = np.arange(len(moment), dtype=np.uint64)
    eligible = (popcount <= degree) & (
        (indices & np.uint64(1 << edge)) == 0
    )
    basis = np.flatnonzero(eligible).astype(np.int64)
    gram = moment[np.bitwise_xor(basis[:, None], basis[None, :])]
    target = response_moment[basis]
    coefficients = np.linalg.solve(gram, target)
    residual = response_energy - float(target @ coefficients)
    if residual < -2e-10:
        raise AssertionError(f"negative least-squares residual {residual}")
    residual = max(0.0, residual)
    normal_error = float(np.max(np.abs(gram @ coefficients - target)))
    condition = float(np.linalg.cond(gram))
    if not math.isfinite(condition) or condition > 200.0:
        raise AssertionError(f"ill-conditioned weighted Gram matrix {condition}")
    if normal_error > 1e-12:
        raise AssertionError(f"weighted normal equations lost accuracy {normal_error}")
    return {
        "degree": degree,
        "basis_dimension": len(basis),
        "weighted_l2_residual": residual,
        "fraction_of_cavity_energy_unexplained": residual / response_energy,
        "gram_condition_number_2": condition,
        "maximum_normal_equation_residual": normal_error,
    }


def run_orientation(
    total_order: int,
    left: np.ndarray,
    right: np.ndarray,
    orientation: int,
) -> dict:
    rows, columns = len(left), len(right)
    edge_count = rows * columns
    raw_t = BETA / math.sqrt(total_order)
    pressure, pressure_audit = exact.bridge_pressures(
        left, right, BETA, total_order, orientation
    )
    shifted = -LAMBDA * pressure
    shifted -= float(np.max(shifted))
    escort = np.exp(shifted)
    escort /= float(np.sum(escort))

    indices = np.arange(len(pressure), dtype=np.uint64)
    complement = indices ^ np.uint64(len(pressure) - 1)
    evenness_error = float(np.max(np.abs(escort - escort[complement])))
    moment = fwht(escort)
    popcount = np.bitwise_count(indices)

    left_cells = vertex_orbits(left)
    right_cells = vertex_orbits(right)
    edge_orbits = []
    for left_cell in left_cells:
        for right_cell in right_cells:
            edge = left_cell[0] * columns + right_cell[0]
            response = cavity_table(pressure, edge, raw_t)
            oddness_error = float(
                np.max(np.abs(response + response[complement]))
            )
            response_moment = fwht(escort * response)
            even_mask = (popcount % 2) == 0
            parity_correlation_error = float(
                np.max(np.abs(response_moment[even_mask]))
            )
            response_energy = float(np.dot(escort, response * response))
            fits = [
                best_weighted_degree(
                    moment,
                    response_moment,
                    response_energy,
                    popcount,
                    edge,
                    degree,
                )
                for degree in DEGREES
            ]
            if (
                fits[1]["weighted_l2_residual"]
                > fits[0]["weighted_l2_residual"] + 1e-12
            ):
                raise AssertionError("larger polynomial space increased residual")
            edge_orbits.append(
                {
                    "representative_edge": [left_cell[0], right_cell[0]],
                    "orbit_size": len(left_cell) * len(right_cell),
                    "left_vertex_orbit": left_cell,
                    "right_vertex_orbit": right_cell,
                    "weighted_cavity_energy": response_energy,
                    "global_oddness_max_error": oddness_error,
                    "maximum_even_degree_weighted_correlation": (
                        parity_correlation_error
                    ),
                    "optimal_polynomial_fits": fits,
                }
            )

    if sum(row["orbit_size"] for row in edge_orbits) != edge_count:
        raise AssertionError("edge-orbit sizes do not cover the bridge")
    return {
        "N": total_order,
        "split": [rows, columns],
        "beta": BETA,
        "raw_t": raw_t,
        "lambda": LAMBDA,
        "orientation": orientation,
        "bridge_escort_global_evenness_max_error": evenness_error,
        "bridge_pressure_audit": pressure_audit,
        "edge_orbits": edge_orbits,
    }


def main() -> None:
    mp.mp.dps = 80
    cache: dict[tuple[int, int], tuple[np.ndarray, dict]] = {}
    children: dict[str, dict] = {}
    records = []
    for total_order, left_order, right_order in PLANS:
        left, left_record = child_record(left_order, total_order, cache)
        right, right_record = child_record(right_order, total_order, cache)
        children[f"N{total_order}_left_{left_order}"] = left_record
        children[f"N{total_order}_right_{right_order}"] = right_record
        for orientation in (-1, 1):
            record = run_orientation(total_order, left, right, orientation)
            records.append(record)
            minimum = min(
                fit["weighted_l2_residual"]
                for edge in record["edge_orbits"]
                for fit in edge["optimal_polynomial_fits"]
                if fit["degree"] == 3
            )
            maximum = max(
                fit["weighted_l2_residual"]
                for edge in record["edge_orbits"]
                for fit in edge["optimal_polynomial_fits"]
                if fit["degree"] == 3
            )
            if total_order == 6 and not (0.015 < minimum <= maximum < 0.021):
                raise AssertionError("order-six residual regression failed")
            if total_order >= 8 and not (0.08 < minimum <= maximum < 0.11):
                raise AssertionError("order-eight/ten residual regression failed")
            print(
                f"N={total_order} split={left_order}+{right_order} "
                f"eps={orientation:+d} degree3=[{minimum:.12g},{maximum:.12g}]",
                flush=True,
            )

    result = {
        "schema": "actual-child-escort-low-degree-falsifier-v1",
        "classification": (
            "complete child/signing and bridge enumeration; high-precision "
            "thermal comparison; numerical transcendental and linear-algebra "
            "evaluation"
        ),
        "parameters": {
            "beta": BETA,
            "lambda": LAMBDA,
            "degrees": list(DEGREES),
            "plans": [list(plan) for plan in PLANS],
            "mp_dps_for_child_selection": mp.mp.dps,
        },
        "weighted_projection_definition": (
            "For each deleted edge e, minimize E_q[(r_e-f(B_-e))^2] "
            "over Walsh polynomials f of total degree at most K. The Gram "
            "matrix is G_ST=E_q chi_(S xor T), and b_S=E_q r_e chi_S."
        ),
        "parity_identity": (
            "Antipodal latent symmetry gives q(B)=q(-B) and "
            "r_e(-B_-e)=-r_e(B_-e). Hence every even-degree weighted "
            "correlation vanishes exactly, so degree 2k gives no improvement "
            "over degree 2k-1. Reported nonzero values are roundoff audits."
        ),
        "children": children,
        "records": records,
        "scope": {
            "proved": (
                "the parity orthogonality identity for every antipodal actual "
                "child channel"
            ),
            "certified_finite_input": (
                "thermal minimizers are selected by exhaustive signing and "
                "exact energy-histogram enumeration"
            ),
            "numerical_finite_output": (
                "complete bridge cubes and weighted least-squares residuals"
            ),
            "not_claimed": [
                "an interval certificate for the displayed residuals",
                "an asymptotic lower bound",
                "failure of degree o(N)",
                "a recurrence or Level-6 implication",
            ],
        },
    }
    output = ROOT / "computations/results/actual_child_escort_low_degree_falsifier.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
