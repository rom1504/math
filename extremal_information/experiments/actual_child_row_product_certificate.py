#!/usr/bin/env python3
"""Global row-product certificates for three actual optimized-child laws.

The script does four things, all on the actual contracted-temperature child
minimizers.

1. It groups every bridge by an exact integer pressure signature.
2. It evaluates selected four-bit row-feature images of the negative escort
   with outward ``mpmath.iv`` intervals.
3. It proves the N=8 coarse reverse-product KL lower bound by interval
   branch-and-bound.
4. It checks the rectangle-Hessian criterion on the coarse and full row
   alphabets and combines the certified lower bounds with the existing
   feasible row-product candidates.

No Monte Carlo sampling is used.  The only floating computation not used as
a proof certificate is the displayed full-alphabet rectangle spectral radius;
its failure margins are large and it is explicitly classified as a numerical
falsification test.  The coarse probabilities and KL bounds are outward
interval certificates.
"""

from __future__ import annotations

import argparse
import functools
import heapq
import itertools
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_row_product_shadow as shadow  # noqa: E402


@dataclass(frozen=True)
class Case:
    name: str
    total_n: int
    beta_text: str
    lambda_text: str
    feature_masks: tuple[int, ...]


CASES = (
    Case("N8_beta4_target", 8, "4", "5.382104195764755", (2, 1, 4, 8)),
    Case("N9_beta2_lambda1", 9, "2", "1", (1, 1, 2, 2)),
    Case("N9_beta4_lambda1", 9, "4", "1", (1, 1, 4, 4)),
)


def iv_cosh(value: mp.iv.mpf) -> mp.iv.mpf:
    return (mp.iv.exp(value) + mp.iv.exp(-value)) / 2


def outward_float_bounds(value: mp.iv.mpf) -> list[float]:
    lower = mp.libmp.to_float(value._mpi_[0], rnd="f")
    upper = mp.libmp.to_float(value._mpi_[1], rnd="c")
    return [
        math.nextafter(lower, -math.inf),
        math.nextafter(upper, math.inf),
    ]


def interval_hist_pressure(
    histogram: np.ndarray, beta_text: str, normalization_order: int
) -> mp.iv.mpf:
    t = mp.iv.mpf(beta_text) / mp.iv.sqrt(
        mp.iv.mpf(str(normalization_order))
    )
    partition = sum(
        (
            int(count) * iv_cosh(t * value)
            for value, count in enumerate(histogram)
            if count
        ),
        mp.iv.mpf("0"),
    ) / int(np.sum(histogram))
    return mp.iv.log(partition)


def certify_histogram_minimum(
    space: exact.SigningSpace,
    winner_matrix: np.ndarray,
    beta_text: str,
    normalization_order: int,
) -> dict:
    spins = exact.projective_spins(len(winner_matrix))
    energies = exact.energies_for_matrix(winner_matrix, spins)
    edge_count = len(winner_matrix) * (len(winner_matrix) - 1) // 2
    winner_histogram = np.bincount(
        np.abs(energies), minlength=edge_count + 1
    )
    winner = interval_hist_pressure(
        winner_histogram, beta_text, normalization_order
    )
    winner_bounds = outward_float_bounds(winner)
    competitor_lower_bounds = []
    for histogram in space.unique_histograms:
        if np.array_equal(histogram, winner_histogram):
            continue
        competitor_lower_bounds.append(
            outward_float_bounds(
                interval_hist_pressure(
                    histogram, beta_text, normalization_order
                )
            )[0]
        )
    next_lower = min(competitor_lower_bounds, default=math.inf)
    if not winner_bounds[1] < next_lower:
        raise AssertionError(
            (winner_bounds, next_lower, beta_text, normalization_order)
        )
    return {
        "order": len(winner_matrix),
        "temperature_normalization_order": normalization_order,
        "minimum_pressure_outward_interval": winner_bounds,
        "next_distinct_histogram_pressure_lower_bound": next_lower,
        "certified_gap_lower_bound": math.nextafter(
            next_lower - winner_bounds[1], -math.inf
        ),
    }


def actual_children(
    case: Case,
) -> tuple[np.ndarray, np.ndarray, list[str], list[dict]]:
    m = case.total_n // 2
    n = case.total_n - m
    matrices = []
    hashes = []
    interval_certificates = []
    for order in (m, n):
        space = exact.build_signing_space(order)
        classes, certificate = exact.thermal_minimizer_classes(
            space, case.beta_text, case.total_n
        )
        if certificate["signed_permutation_global_sign_class_count"] != 1:
            raise AssertionError((case, order, certificate))
        matrices.append(
            np.asarray(classes[0]["representative_matrix"], dtype=np.int16)
        )
        hashes.append(classes[0]["representative_sha256"])
        interval_certificates.append(
            certify_histogram_minimum(
                space,
                matrices[-1],
                case.beta_text,
                case.total_n,
            )
        )
    return matrices[0], matrices[1], hashes, interval_certificates


def same_temperature_target_certificate(case: Case) -> dict:
    component_certificates = []
    for order in (case.total_n // 2, case.total_n - case.total_n // 2):
        space = exact.build_signing_space(order)
        classes, selector = exact.thermal_minimizer_classes(
            space, case.beta_text, order
        )
        if selector["minimizing_histogram_count"] != 1:
            raise AssertionError((case, order, selector))
        matrix = np.asarray(
            classes[0]["representative_matrix"], dtype=np.int16
        )
        component_certificates.append(
            certify_histogram_minimum(
                space, matrix, case.beta_text, order
            )
        )
    lower = math.nextafter(
        sum(
            item["minimum_pressure_outward_interval"][0]
            for item in component_certificates
        ),
        -math.inf,
    )
    upper = math.nextafter(
        sum(
            item["minimum_pressure_outward_interval"][1]
            for item in component_certificates
        ),
        math.inf,
    )
    return {
        "components": component_certificates,
        "same_temperature_minimum_child_target_outward_interval": [
            lower,
            upper,
        ],
    }


def enumerate_pressure_signatures(
    case: Case,
    left: np.ndarray,
    right: np.ndarray,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return exact signature vectors and signature-by-feature counts."""

    total_n = case.total_n
    m = total_n // 2
    n = total_n - m
    d = m * n
    x = exact.projective_spins(m).astype(np.int16)
    y = exact.projective_spins(n).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
    # All three audited cases use relative orientation -1.
    internal = np.abs(ex[:, None] - ey[None, :]).reshape(-1)
    max_internal = int(np.max(internal))
    width = d + 1
    signature_width = (max_internal + 1) * width
    rank_one = np.asarray(
        [
            (xx[:, None] * yy[None, :]).reshape(-1)
            for xx in x
            for yy in y
        ],
        dtype=np.int16,
    )

    counts: dict[bytes, np.ndarray] = {}
    cube_size = 1 << d
    started = time.time()
    for lower in range(0, cube_size, batch_size):
        masks = np.arange(
            lower, min(lower + batch_size, cube_size), dtype=np.uint64
        )
        bits = (
            (masks[:, None] >> np.arange(d, dtype=np.uint64)) & 1
        ).astype(np.int16)
        bridges = 1 - 2 * bits
        cross = np.abs(bridges @ rank_one.T).astype(np.int16)
        bins = internal[None, :] * width + cross
        for local_index, mask_value in enumerate(masks):
            signature = np.bincount(
                bins[local_index], minlength=signature_width
            ).astype(np.uint8).tobytes()
            row_values = [
                (int(mask_value) >> (row * n)) & ((1 << n) - 1)
                for row in range(m)
            ]
            code = 0
            for row, feature_mask in enumerate(case.feature_masks):
                parity = bin(row_values[row] & feature_mask).count("1") & 1
                code |= parity << row
            if signature not in counts:
                counts[signature] = np.zeros(1 << m, dtype=np.int64)
            counts[signature][code] += 1

    signatures = np.frombuffer(b"".join(counts), dtype=np.uint8).reshape(
        len(counts), signature_width
    )
    count_matrix = np.stack(list(counts.values()))
    if int(np.sum(count_matrix)) != cube_size:
        raise AssertionError((np.sum(count_matrix), cube_size))
    metadata = {
        "bridge_cube_size": cube_size,
        "rank_one_character_count": len(rank_one),
        "exact_pressure_signature_count": len(signatures),
        "signature_width": signature_width,
        "internal_bin_width": width,
        "wall_time_seconds": time.time() - started,
    }
    return signatures, count_matrix, metadata


def interval_coarse_law(
    case: Case,
    signatures: np.ndarray,
    counts: np.ndarray,
) -> tuple[list[list[float]], list[list[float]], dict]:
    total_n = case.total_n
    m = total_n // 2
    n = total_n - m
    d = m * n
    width = d + 1
    t = mp.iv.mpf(case.beta_text) / mp.iv.sqrt(mp.iv.mpf(str(total_n)))
    lam = mp.iv.mpf(case.lambda_text)
    max_internal = signatures.shape[1] // width - 1
    cross_cosh = [iv_cosh(t * value) for value in range(d + 1)]
    internal_cosh = [
        iv_cosh(t * value) for value in range(max_internal + 1)
    ]
    rank_one_count = (1 << (m - 1)) * (1 << (n - 1))
    weights = []
    for signature in signatures:
        partition = mp.iv.mpf("0")
        for index, multiplicity in enumerate(signature):
            if multiplicity:
                partition += (
                    int(multiplicity)
                    * internal_cosh[index // width]
                    * cross_cosh[index % width]
                )
        partition /= rank_one_count
        weights.append(partition ** (-lam))

    numerators = [
        sum(
            (
                int(counts[index, code]) * weights[index]
                for index in range(len(weights))
            ),
            mp.iv.mpf("0"),
        )
        for code in range(1 << m)
    ]
    denominator = sum(numerators, mp.iv.mpf("0"))
    probabilities = [value / denominator for value in numerators]
    log_probabilities = [mp.iv.log(value) for value in probabilities]
    coefficients = []
    for mask in range(1 << m):
        coefficient = sum(
            (
                (-1) ** (bin(code & mask).count("1"))
                * log_probabilities[code]
                for code in range(1 << m)
            ),
            mp.iv.mpf("0"),
        ) / (1 << m)
        coefficients.append(coefficient)
    soft_pressure = -mp.iv.log(denominator / (1 << d)) / lam
    return (
        [outward_float_bounds(value) for value in probabilities],
        [outward_float_bounds(value) for value in coefficients],
        {
            "negative_moment_soft_pressure_outward_interval": (
                outward_float_bounds(soft_pressure)
            ),
        },
    )


def interval_product(left: tuple[float, float], right: tuple[float, float]):
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return (min(products), max(products))


def interval_square(interval: tuple[float, float]):
    lower, upper = interval
    return (
        0.0 if lower <= 0 <= upper else min(lower * lower, upper * upper),
        max(lower * lower, upper * upper),
    )


@functools.lru_cache(maxsize=None)
def negative_binary_entropy_outward(mean: float) -> tuple[float, float]:
    """Validated enclosure of ``-h((1+mean)/2)`` for a binary64 mean.

    Every branch-and-bound endpoint is dyadic.  Reconstructing that rational
    from ``as_integer_ratio`` prevents the platform ``libm`` implementation
    from entering the proof certificate.
    """

    if mean <= -1.0 or mean >= 1.0:
        return (0.0, 0.0)
    numerator, denominator = mean.as_integer_ratio()
    value = mp.iv.mpf(numerator) / mp.iv.mpf(denominator)
    probability = (1 + value) / 2
    entropy = probability * mp.iv.log(probability) + (
        1 - probability
    ) * mp.iv.log(1 - probability)
    bounds = outward_float_bounds(entropy)
    return (bounds[0], bounds[1])


def n8_two_dimensional_bounds(
    coefficient_bounds: list[list[float]],
    target: float,
) -> dict:
    # Exact signature-count equalities prove the coefficient identifications;
    # take interval hulls to retain outward rounding.
    c0 = (
        -coefficient_bounds[0][1],
        -coefficient_bounds[0][0],
    )
    j = (
        min(coefficient_bounds[3][0], coefficient_bounds[12][0]),
        max(coefficient_bounds[3][1], coefficient_bounds[12][1]),
    )
    k = (
        min(coefficient_bounds[index][0] for index in (5, 6, 9, 10)),
        max(coefficient_bounds[index][1] for index in (5, 6, 9, 10)),
    )
    four = (
        -coefficient_bounds[15][1],
        -coefficient_bounds[15][0],
    )
    if j[0] - four[1] <= 0 or k[0] <= 0:
        raise AssertionError((j, k, four))

    # The entropy endpoint itself is enclosed by arbitrary-precision interval
    # arithmetic.  On this fixed certificate all remaining basic binary64
    # intermediates have magnitude below 16 and there are fewer than 64
    # operations per box.  The 1e-12 subtraction is more than four times the
    # standard worst-case 64*16*2^-52 accumulation bound.
    arithmetic_safety = 1e-12
    if max(
        abs(endpoint)
        for interval in (c0, j, k, four)
        for endpoint in interval
    ) >= 8:
        raise AssertionError("binary64 safety proof range exceeded")

    def lower_bound(box: tuple[tuple[float, float], ...]) -> float:
        u, v = box
        u2 = interval_square(u)
        v2 = interval_square(v)
        uv = interval_product(u, v)
        u2v2 = interval_product(u2, v2)
        value = c0[0]
        value -= j[1] * (u2[1] + v2[1])
        value -= 4 * k[1] * uv[1]
        value += four[0] * u2v2[0]
        for interval in box:
            closest = (
                0.0
                if interval[0] <= 0 <= interval[1]
                else min(interval, key=abs)
            )
            value += 2 * negative_binary_entropy_outward(closest)[0]
        return math.nextafter(value - arithmetic_safety, -math.inf)

    initial = ((0.0, 1.0), (0.0, 1.0))
    queue: list[tuple[float, int, tuple[tuple[float, float], ...]]] = [
        (lower_bound(initial), 0, initial)
    ]
    serial = 1
    processed = 0
    minimum_pruned_bound = math.inf
    while queue:
        lower, _, box = heapq.heappop(queue)
        processed += 1
        if lower >= target:
            minimum_pruned_bound = min(minimum_pruned_bound, lower)
            continue
        widths = [interval[1] - interval[0] for interval in box]
        axis = 0 if widths[0] >= widths[1] else 1
        endpoint_lower, endpoint_upper = box[axis]
        midpoint = (endpoint_lower + endpoint_upper) / 2
        for segment in (
            (endpoint_lower, midpoint),
            (midpoint, endpoint_upper),
        ):
            child = box[:axis] + (segment,) + box[axis + 1 :]
            child_lower = lower_bound(child)
            if child_lower < target:
                heapq.heappush(queue, (child_lower, serial, child))
                serial += 1
            else:
                minimum_pruned_bound = min(
                    minimum_pruned_bound, child_lower
                )

    candidate_mean = 0.74237662
    candidate_negative_entropy = negative_binary_entropy_outward(
        candidate_mean
    )
    candidate_lower = (
        c0[0]
        - 2 * j[1] * candidate_mean**2
        - 4 * k[1] * candidate_mean**2
        + four[0] * candidate_mean**4
        + 4 * candidate_negative_entropy[0]
    )
    candidate_upper = (
        c0[1]
        - 2 * j[0] * candidate_mean**2
        - 4 * k[0] * candidate_mean**2
        + four[1] * candidate_mean**4
        + 4 * candidate_negative_entropy[1]
    )
    return {
        "certified_lower_bound": target,
        "feasible_coarse_product_mean": candidate_mean,
        "feasible_coarse_product_value_interval": [
            math.nextafter(candidate_lower - arithmetic_safety, -math.inf),
            math.nextafter(candidate_upper + arithmetic_safety, math.inf),
        ],
        "branch_and_bound_boxes_processed": processed,
        "minimum_pruned_box_lower_bound": minimum_pruned_bound,
        "per_box_binary64_safety_subtraction": arithmetic_safety,
        "dimension_reduction": (
            "exact signature classes imply twin pairs; first-order "
            "monotonicity forces s1=s2 and s3=s4, then complement/sign "
            "symmetry reduces to u,v in [0,1]"
        ),
    }


def coarse_rectangle_upper(
    coefficient_bounds: list[list[float]],
) -> tuple[np.ndarray, float]:
    """Upper-bound the binary potential's rectangle matrix."""

    m = int(round(math.log2(len(coefficient_bounds))))
    # Potential coefficients are the negatives of log Q coefficients.
    potential = [(-bounds[1], -bounds[0]) for bounds in coefficient_bounds]
    matrix = np.zeros((m, m), dtype=np.float64)
    for i in range(m):
        for j in range(i):
            rest = [index for index in range(m) if index not in (i, j)]
            maximum = 0.0
            for signs in itertools.product((-1, 1), repeat=len(rest)):
                total = (0.0, 0.0)
                for mask, coefficient in enumerate(potential):
                    if not ((mask >> i) & 1 and (mask >> j) & 1):
                        continue
                    sign = math.prod(
                        signs[rest.index(index)]
                        for index in rest
                        if (mask >> index) & 1
                    )
                    term = (
                        coefficient
                        if sign > 0
                        else (-coefficient[1], -coefficient[0])
                    )
                    total = (total[0] + term[0], total[1] + term[1])
                maximum = max(maximum, 4 * abs(total[0]), 4 * abs(total[1]))
            matrix[i, j] = matrix[j, i] = math.nextafter(
                maximum + 1e-12, math.inf
            )
    return matrix, float(np.max(np.sum(matrix, axis=1)))


def full_rectangle_matrix(
    pressure: np.ndarray, rows: int, columns: int
) -> np.ndarray:
    row_size = 1 << columns
    tensor = shadow.pressure_tensor(pressure, rows, row_size)
    matrix = np.zeros((rows, rows), dtype=np.float64)
    pairs = [(a, b) for a in range(row_size) for b in range(a)]
    for i in range(rows):
        for j in range(i):
            axes = [i, j] + [
                index for index in range(rows) if index not in (i, j)
            ]
            view = np.transpose(tensor, axes).reshape(
                row_size, row_size, -1
            )
            maximum = 0.0
            for a, b in pairs:
                difference = view[a] - view[b]
                ranges = np.max(difference, axis=0) - np.min(
                    difference, axis=0
                )
                maximum = max(maximum, float(np.max(ranges)))
            matrix[i, j] = matrix[j, i] = maximum
    return matrix


def find_shadow_record(case: Case) -> dict:
    if case.name == "N8_beta4_target":
        path = (
            ROOT
            / "computations/results/actual_child_row_product_shadow_target_threshold_n8.json"
        )
    else:
        path = ROOT / "computations/results/actual_child_row_product_shadow.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    m = case.total_n // 2
    n = case.total_n - m
    beta = float(case.beta_text)
    lam = float(case.lambda_text)
    for record in payload["records"]:
        if (
            record["N"] == case.total_n
            and record["split"] == [m, n]
            and record["beta"] == beta
            and record["relative_child_orientation"] == -1
        ):
            for law in record["laws"]:
                if abs(law["lambda"] - lam) < 1e-12:
                    return law
    raise LookupError(case)


def run_case(case: Case, batch_size: int) -> dict:
    left, right, child_hashes, child_interval_certificates = actual_children(case)
    signatures, counts, signature_metadata = enumerate_pressure_signatures(
        case, left, right, batch_size
    )
    probability_bounds, coefficient_bounds, interval_soft = interval_coarse_law(
        case, signatures, counts
    )
    target_certificate = same_temperature_target_certificate(case)
    m = case.total_n // 2
    n = case.total_n - m
    complement_equal = all(
        np.array_equal(counts[:, code], counts[:, code ^ ((1 << m) - 1)])
        for code in range(1 << m)
    )
    if not complement_equal:
        raise AssertionError((case, "coarse complement symmetry failed"))

    pressure, pressure_audit = exact.bridge_pressures(
        left, right, float(case.beta_text), case.total_n, -1
    )
    full_rectangle = full_rectangle_matrix(pressure, m, n)
    full_spectral_radius = float(np.linalg.eigvalsh(full_rectangle)[-1])
    lam = float(case.lambda_text)
    shadow_record = find_shadow_record(case)
    # The imported coordinate-Gibbs files record binary64 feasible values but
    # not outward intervals.  Inflate every imported upper endpoint by 1e-6;
    # this is far above their reported residuals and roundoff audits.  The
    # soft pressure is recomputed with outward interval arithmetic below.
    feasible_upper_safety = 1e-6

    result = {
        "case": case.name,
        "N": case.total_n,
        "split": [m, n],
        "beta": float(case.beta_text),
        "lambda": lam,
        "relative_child_orientation": -1,
        "actual_child_representative_sha256": child_hashes,
        "actual_child_minimality_interval_certificates": (
            child_interval_certificates
        ),
        "same_temperature_target_interval_certificate": target_certificate,
        "row_Walsh_feature_masks": list(case.feature_masks),
        "signature_enumeration": signature_metadata,
        "coarse_probability_outward_intervals": probability_bounds,
        "coarse_log_probability_Walsh_coefficient_outward_intervals": (
            coefficient_bounds
        ),
        "negative_moment_soft_pressure_interval_certificate": interval_soft,
        "exact_signature_count_complement_symmetry": complement_equal,
        "full_row_rectangle_matrix_numerical": full_rectangle.tolist(),
        "full_row_rectangle_spectral_radius_numerical": full_spectral_radius,
        "full_row_strong_convexity_left_side": lam * full_spectral_radius,
        "full_row_strong_convexity_required_upper_bound": 4.0,
        "full_row_condition_passes": lam * full_spectral_radius < 4.0,
        "pressure_cube_audit": pressure_audit,
        **interval_soft,
        "exact_negative_moment_soft_pressure": shadow_record[
            "exact_negative_moment_soft_pressure"
        ],
        "feasible_full_row_product_objective": shadow_record[
            "best_evaluated_product_objective"
        ],
        "imported_feasible_upper_endpoint_safety": feasible_upper_safety,
        "feasible_full_reverse_projection_upper_bound": shadow_record[
            "candidate_reverse_projection_upper_bound"
        ],
    }

    if case.name == "N8_beta4_target":
        groups = (
            (0, 15),
            (3, 12),
            (1, 2, 4, 7, 8, 11, 13, 14),
            (5, 6, 9, 10),
        )
        equalities = [
            all(
                np.array_equal(counts[:, group[0]], counts[:, code])
                for code in group[1:]
            )
            for group in groups
        ]
        if not all(equalities):
            raise AssertionError((case, groups, equalities))
        proof = n8_two_dimensional_bounds(coefficient_bounds, 1.075)
        lower_i = proof["certified_lower_bound"]
        upper_i = (
            shadow_record["candidate_reverse_projection_upper_bound"]
            + feasible_upper_safety
        )
        soft_interval = interval_soft[
            "negative_moment_soft_pressure_outward_interval"
        ]
        target_interval = target_certificate[
            "same_temperature_minimum_child_target_outward_interval"
        ]
        vrow_lower = outward_float_bounds(
            mp.iv.mpf([soft_interval[0], soft_interval[1]])
            + mp.iv.mpf(str(lower_i)) / mp.iv.mpf(case.lambda_text)
        )[0]
        result.update(
            {
                "exact_signature_count_equality_groups": [
                    list(group) for group in groups
                ],
                "coarse_reverse_product_projection_certificate": proof,
                "full_reverse_projection_certified_interval": [
                    lower_i,
                    upper_i,
                ],
                "V_row_certified_interval": [
                    vrow_lower,
                    shadow_record["best_evaluated_product_objective"]
                    + feasible_upper_safety,
                ],
                "target_excess_certified_lower_bound": math.nextafter(
                    vrow_lower - target_interval[1], -math.inf
                ),
            }
        )
    else:
        coarse_matrix, row_sum_upper = coarse_rectangle_upper(
            coefficient_bounds
        )
        if row_sum_upper >= 4:
            raise AssertionError((case, coarse_matrix, row_sum_upper))
        constant_potential = (
            -coefficient_bounds[0][1],
            -coefficient_bounds[0][0],
        )
        log_two = outward_float_bounds(mp.iv.log(mp.iv.mpf(2)))
        uniform_gap = [
            math.nextafter(constant_potential[0] - m * log_two[1], -math.inf),
            math.nextafter(constant_potential[1] - m * log_two[0], math.inf),
        ]
        soft_interval = interval_soft[
            "negative_moment_soft_pressure_outward_interval"
        ]
        result.update(
            {
                "coarse_rectangle_matrix_outward_upper": coarse_matrix.tolist(),
                "coarse_rectangle_max_row_sum_upper": row_sum_upper,
                "coarse_uniform_is_exact_unique_product_minimizer": True,
                "coarse_reverse_product_projection_outward_interval": uniform_gap,
                "full_reverse_projection_certified_interval": [
                    uniform_gap[0],
                    shadow_record["candidate_reverse_projection_upper_bound"]
                    + feasible_upper_safety,
                ],
                "V_row_certified_interval": [
                    math.nextafter(
                        soft_interval[0] + uniform_gap[0] / lam,
                        -math.inf,
                    ),
                    shadow_record["best_evaluated_product_objective"]
                    + feasible_upper_safety,
                ],
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--iv-dps", type=int, default=70)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "computations/results/actual_child_row_product_certificate.json"
        ),
    )
    args = parser.parse_args()
    mp.mp.dps = args.mp_dps
    mp.iv.dps = args.iv_dps
    negative_binary_entropy_outward.cache_clear()
    started = time.time()
    results = [run_case(case, args.batch_size) for case in CASES]
    payload = {
        "schema": "actual-child-row-product-global-certificate-v1",
        "classification": (
            "complete actual-child bridge enumeration; exact integer "
            "pressure signatures; outward interval coarse laws; rigorous "
            "coarse reverse-product lower certificates; feasible full "
            "row-product upper bounds; full-alphabet rectangle radii are "
            "floating falsification diagnostics only"
        ),
        "theorem_note": (
            "extremal_information/drafts/"
            "actual_child_row_product_global_certificate.md"
        ),
        "parameters": {
            "mp_dps": args.mp_dps,
            "iv_dps": args.iv_dps,
            "batch_size": args.batch_size,
        },
        "cases": results,
        "wall_time_seconds": time.time() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
