#!/usr/bin/env python3
"""Audit the concrete asymmetric law X uniform, Y = sign(A X).

All quadratic energies use doubled normalization ``z.T @ A @ z``.  At a
zero field this checker uses the switching-equivariant convention ``Y_i=X_i``.
The response ``X.T @ A @ Y = ||A X||_1`` is unaffected by that convention.
The augmented greedy witness uses the original free-shore spin to break a
zero cross-field tie, making that certificate switching covariant as well.

For orders at most 14 every projective X is enumerated.  The order-18
conference response and quadratic cap are enumerated exactly, while its shore
defects are sampled.  Larger conference cases use reproducible Monte Carlo and
the polynomial spectral shore certificates only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "computations" / "results"
ARCSIN_KAPPA = math.pi / 2 - 1
CURRENT_DOUBLED_CONSTANT = 0.672986728863


def nested(payload: object, key: str) -> object:
    value = payload
    for part in key.split("."):
        if isinstance(value, dict):
            value = value[part]
        elif isinstance(value, list):
            value = value[int(part)]
        else:
            raise AssertionError((key, part))
    return value


def matrix_hash(matrix: np.ndarray) -> str:
    encoded = json.dumps(matrix.astype(int).tolist(), separators=(",", ":"))
    return hashlib.sha256(encoded.encode()).hexdigest()


def validate_signing(matrix: np.ndarray) -> None:
    n = len(matrix)
    if matrix.shape != (n, n):
        raise AssertionError(matrix.shape)
    if not np.array_equal(matrix, matrix.T):
        raise AssertionError("matrix is not symmetric")
    if np.any(np.diag(matrix)):
        raise AssertionError("matrix diagonal is not zero")
    off = matrix[~np.eye(n, dtype=bool)]
    if not np.all(np.isin(off, (-1, 1))):
        raise AssertionError("matrix is not a signing")


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins modulo global negation, coordinate zero fixed +1."""
    indices = np.arange(1 << (n - 1), dtype=np.uint64)[:, None]
    bits = (
        (indices >> np.arange(n - 1, dtype=np.uint64)[None, :]) & 1
    ).astype(np.int8)
    return np.concatenate(
        (np.ones((len(indices), 1), dtype=np.int8), 1 - 2 * bits), axis=1
    )


def random_projective_spins(n: int, count: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    tails = rng.integers(0, 2, size=(count, n - 1), dtype=np.int8)
    return np.concatenate(
        (np.ones((count, 1), dtype=np.int8), 1 - 2 * tails), axis=1
    )


def exact_response_fraction(n: int) -> Fraction:
    """Return E ||A X||_1, the same rational number for every signing A."""
    walk_length = n - 1
    numerator = (
        n
        * walk_length
        * math.comb(walk_length - 1, (walk_length - 1) // 2)
    )
    denominator = 1 << (walk_length - 1)
    return Fraction(numerator, denominator)


def fraction_payload(value: Fraction) -> Dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


_SPIN_CACHE: Dict[int, np.ndarray] = {}


def spins_of_order(n: int) -> np.ndarray:
    if n not in _SPIN_CACHE:
        _SPIN_CACHE[n] = projective_spins(n).astype(np.int32)
    return _SPIN_CACHE[n]


def exact_one_sided_caps(block: np.ndarray) -> Tuple[int, int]:
    m = len(block)
    if m <= 1:
        return 0, 0
    spins = spins_of_order(m)
    fields = spins @ block.astype(np.int32)
    energies = np.sum(fields * spins, axis=1, dtype=np.int64)
    return int(energies.max()), int(-energies.min())


def projector_floor_from_eigenvalues(
    eigenvalues: np.ndarray, target_sign: int
) -> float:
    signed = target_sign * eigenvalues
    positive = signed[signed > 1e-9]
    rank = len(positive)
    if rank == 0:
        return 0.0
    mass = float(positive.sum())
    theta = min(1.0, mass / (2 * ARCSIN_KAPPA * rank))
    return max(
        0.0,
        (2 / math.pi)
        * (theta * mass - ARCSIN_KAPPA * theta * theta * rank),
    )


def weighted_projector_floor(matrix: np.ndarray, target_sign: int) -> float:
    """Weighted projector certificate for a symmetric zero-diagonal matrix.

    Unlike ``projector_floor_from_eigenvalues``, the arcsine remainder keeps
    the actual absolute edge weights.
    """
    m = len(matrix)
    if m <= 1:
        return 0.0
    eigenvalues, eigenvectors = np.linalg.eigh(target_sign * matrix.astype(float))
    mask = eigenvalues > 1e-9
    if not np.any(mask):
        return 0.0
    projector = eigenvectors[:, mask] @ eigenvectors[:, mask].T
    mass = float(eigenvalues[mask].sum())
    weights = np.abs(matrix.astype(float))
    off_diagonal_penalty = float(
        np.sum(weights * projector * projector)
        - np.sum(np.diag(weights) * np.diag(projector) ** 2)
    )
    if off_diagonal_penalty <= 1e-15:
        theta = 1.0
    else:
        theta = min(1.0, mass / (2 * ARCSIN_KAPPA * off_diagonal_penalty))
    return max(
        0.0,
        (2 / math.pi)
        * (theta * mass - ARCSIN_KAPPA * theta * theta * off_diagonal_penalty),
    )


def greedy_aligned_floor(
    matrix: np.ndarray, target_sign: int, zero_ties: np.ndarray
) -> Tuple[int, int, int]:
    """A deterministic aligned one-flip local-search witness.

    The last column of ``matrix`` is the collapsed cross field.  Initialize
    the free spins to align its contribution with ``target_sign`` (using the
    supplied original free spin at a zero cross field), fix the collapsed
    coordinate initially to +1, and repeatedly take the largest strictly
    improving single-coordinate flip.  This tie rule is switching covariant.
    The returned energy is a rigorous lower bound on the aligned one-sided
    cap; the flip count is for the complexity audit.
    """
    dimension = len(matrix)
    if dimension <= 1:
        return 0, 0, 0
    cross_field = matrix[:-1, -1]
    if len(zero_ties) != len(cross_field):
        raise AssertionError((len(zero_ties), len(cross_field)))
    aligned_free_spin = np.where(
        target_sign * cross_field > 0,
        1,
        np.where(target_sign * cross_field < 0, -1, zero_ties),
    )
    spin = np.concatenate(
        (
            aligned_free_spin,
            np.ones(1, dtype=np.int32),
        )
    ).astype(np.int32)
    flips = 0
    # The absolute undirected edge weight is O(n^2) for collapsed shores.
    # Every strict integer improvement is at least four, so this guard is a
    # proved polynomial upper bound, not a heuristic iteration cutoff.
    edge_weight = int(np.abs(np.triu(matrix, 1)).sum())
    initial_energy = max(0, int(target_sign * (spin @ matrix @ spin)))
    fields = matrix @ spin
    while flips <= edge_weight:
        gains = -4 * target_sign * spin * fields
        coordinate = int(np.argmax(gains))
        if int(gains[coordinate]) <= 0:
            if not np.array_equal(fields, matrix @ spin):
                raise AssertionError("incremental field update drifted")
            energy = int(target_sign * (spin @ matrix @ spin))
            terminal_margins = target_sign * spin * fields
            if np.any(terminal_margins < 0):
                raise AssertionError("terminal spin is not one-flip stable")
            if energy != int(np.abs(fields).sum()):
                raise AssertionError("terminal l1 identity failed")
            terminal_energy = max(0, energy)
            if terminal_energy < initial_energy:
                raise AssertionError("coordinate ascent decreased objective")
            return initial_energy, terminal_energy, flips
        old_spin = int(spin[coordinate])
        spin[coordinate] = -old_spin
        fields = fields - 2 * old_spin * matrix[:, coordinate]
        flips += 1
    raise AssertionError("strict coordinate ascent exceeded edge-weight bound")


def shore_statistics(block: np.ndarray, include_exact: bool) -> Dict[str, object]:
    m = len(block)
    if m <= 1:
        answer = {
            "order": m,
            "projector_plus": 0.0,
            "projector_minus": 0.0,
            "nuclear": 0.0,
        }
        if include_exact:
            answer.update({"exact_plus": 0, "exact_minus": 0})
        return answer

    eigenvalues = np.linalg.eigvalsh(block.astype(np.float64))
    nuclear_norm = float(np.abs(eigenvalues).sum())
    answer = {
        "order": m,
        "projector_plus": projector_floor_from_eigenvalues(eigenvalues, 1),
        "projector_minus": projector_floor_from_eigenvalues(eigenvalues, -1),
        "nuclear": max(
            0.0,
            nuclear_norm / math.pi - (1 - 2 / math.pi) * m,
        ),
    }
    if include_exact:
        positive, negative = exact_one_sided_caps(block)
        answer.update({"exact_plus": positive, "exact_minus": negative})
    return answer


def key_from_indicator(indicator: np.ndarray) -> bytes:
    return np.packbits(indicator.astype(np.uint8), bitorder="little").tobytes()


def indicator_from_key(key: bytes, n: int) -> np.ndarray:
    raw = np.frombuffer(key, dtype=np.uint8)
    return np.unpackbits(raw, bitorder="little")[:n].astype(bool)


def mean_standard_error(values: np.ndarray, exact: bool) -> Tuple[float, float]:
    mean = float(np.mean(values))
    if exact or len(values) <= 1:
        return mean, 0.0
    return mean, float(np.std(values, ddof=1) / math.sqrt(len(values)))


def exact_quadratic_cap(matrix: np.ndarray) -> int:
    n = len(matrix)
    spins = projective_spins(n).astype(np.int32)
    best = 0
    for start in range(0, len(spins), 32768):
        chunk = spins[start : start + 32768]
        energies = np.sum((chunk @ matrix) * chunk, axis=1, dtype=np.int64)
        best = max(best, int(np.max(np.abs(energies))))
    return best


def audit_law(
    matrix: np.ndarray,
    label: str,
    source: str,
    source_key: str,
    mode: str,
    sample_count: Optional[int],
    seed: int,
    include_exact_face: bool,
    independently_enumerate_cap: bool,
    recorded_quadratic_cap: Optional[int] = None,
) -> Dict[str, object]:
    validate_signing(matrix)
    n = len(matrix)
    population = 1 << (n - 1)
    exact = mode == "exact"
    if exact:
        spins = projective_spins(n)
    else:
        if sample_count is None:
            raise AssertionError("sample count required")
        spins = random_projective_spins(n, sample_count, seed)

    x = spins.astype(np.int32)
    fields = x @ matrix.astype(np.int32)
    # Switching-equivariant Boolean completion at zero fields.
    y = np.where(fields > 0, 1, np.where(fields < 0, -1, x)).astype(np.int32)
    agreement = x == y
    p = np.where(agreement, x, 0)
    q = np.where(~agreement, x, 0)
    energy_i = np.sum((p @ matrix) * p, axis=1, dtype=np.int64)
    energy_j = np.sum((q @ matrix) * q, axis=1, dtype=np.int64)
    response = np.sum(np.abs(fields), axis=1, dtype=np.int64)
    if not np.array_equal(response, energy_i - energy_j):
        raise AssertionError("agreement decomposition failed")

    packed = np.packbits(agreement.astype(np.uint8), axis=1, bitorder="little")
    unique_packed, inverse, multiplicities = np.unique(
        packed, axis=0, return_inverse=True, return_counts=True
    )
    cache: Dict[bytes, Dict[str, object]] = {}

    def get_stats(indicator: np.ndarray) -> Dict[str, object]:
        key = key_from_indicator(indicator)
        if key not in cache:
            vertices = np.flatnonzero(indicator)
            block = matrix[np.ix_(vertices, vertices)]
            cache[key] = shore_statistics(block, include_exact_face)
        return cache[key]

    exact_delta_by_type = (
        np.zeros((len(unique_packed), 2, 2), dtype=np.float64)
        if include_exact_face
        else None
    )
    projector_by_type = np.zeros((len(unique_packed), 2, 2), dtype=np.float64)
    nuclear_by_type = np.zeros((len(unique_packed), 2), dtype=np.float64)
    shore_size_by_type = np.zeros((len(unique_packed), 2), dtype=np.int16)

    # Axis 1 is shore (I,J); axis 2 is sign (-,+) encoded as (0,1).
    for type_index, packed_row in enumerate(unique_packed):
        indicator = indicator_from_key(packed_row.tobytes(), n)
        stats_i = get_stats(indicator)
        stats_j = get_stats(~indicator)
        shore_size_by_type[type_index] = (int(indicator.sum()), int((~indicator).sum()))
        for shore_index, stats in enumerate((stats_i, stats_j)):
            projector_by_type[type_index, shore_index, 0] = float(
                stats["projector_minus"]
            )
            projector_by_type[type_index, shore_index, 1] = float(
                stats["projector_plus"]
            )
            nuclear_by_type[type_index, shore_index] = float(stats["nuclear"])
            if include_exact_face and exact_delta_by_type is not None:
                exact_delta_by_type[type_index, shore_index, 0] = float(
                    stats["exact_minus"]
                )
                exact_delta_by_type[type_index, shore_index, 1] = float(
                    stats["exact_plus"]
                )

    opposite = energy_i * energy_j < 0
    sign_i_index = (energy_i > 0).astype(np.int8)
    sign_j_index = (energy_j > 0).astype(np.int8)
    # delta_J uses the sign of P on J; delta_I uses the sign of R on I.
    projector_j = projector_by_type[inverse, 1, sign_i_index]
    projector_i = projector_by_type[inverse, 0, sign_j_index]
    projector_delta = np.where(
        opposite,
        np.minimum(
            np.maximum(np.abs(energy_j) - projector_j, 0.0),
            np.maximum(np.abs(energy_i) - projector_i, 0.0),
        ),
        0.0,
    )
    nuclear_j = nuclear_by_type[inverse, 1]
    nuclear_i = nuclear_by_type[inverse, 0]
    nuclear_delta = np.where(
        opposite,
        np.minimum(
            np.maximum(np.abs(energy_j) - nuclear_j, 0.0),
            np.maximum(np.abs(energy_i) - nuclear_i, 0.0),
        ),
        0.0,
    )

    exact_delta = None
    if include_exact_face and exact_delta_by_type is not None:
        exact_j = exact_delta_by_type[inverse, 1, sign_i_index]
        exact_i = exact_delta_by_type[inverse, 0, sign_j_index]
        exact_delta = np.where(
            opposite,
            np.minimum(
                np.maximum(np.abs(energy_j) - exact_j, 0.0),
                np.maximum(np.abs(energy_i) - exact_i, 0.0),
            ),
            0.0,
        )

    # Cross-aware recoupling.  Anchor on I: h=A[J,I]p and augment A[J] by
    # one Boolean coordinate coupled to h.  Anchor on J symmetrically.  The
    # cheap witness takes the free shore r=sign(h), using the original free
    # spin at zero.  Exact/projector data depend only on the agreement subset
    # and anchor restriction; tied polynomial witnesses additionally depend
    # on the original free restriction.  Cache those two layers separately.
    cross_base_cache: Dict[Tuple[bytes, bytes, str], Dict[str, object]] = {}
    cross_witness_cache: Dict[
        Tuple[Tuple[bytes, bytes, str], bytes], Dict[str, object]
    ] = {}
    cross_state_order_histogram: Dict[int, int] = {}
    greedy_cached_flip_counts: List[int] = []

    def cross_channel(
        anchor_indicator: np.ndarray,
        anchor_spin: np.ndarray,
        anchor_energy: int,
        tag: str,
    ) -> Dict[str, object]:
        subset_key = key_from_indicator(anchor_indicator)
        anchor_vertices = np.flatnonzero(anchor_indicator)
        free_vertices = np.flatnonzero(~anchor_indicator)
        anchor_values = anchor_spin[anchor_vertices]
        free_values = anchor_spin[free_vertices]
        anchor_key = np.packbits(
            (anchor_values > 0).astype(np.uint8), bitorder="little"
        ).tobytes()
        free_key = np.packbits(
            (free_values > 0).astype(np.uint8), bitorder="little"
        ).tobytes()
        base_key = (subset_key, anchor_key, tag)
        if base_key not in cross_base_cache:
            free_block = matrix[np.ix_(free_vertices, free_vertices)]
            h = matrix[np.ix_(free_vertices, anchor_vertices)] @ anchor_values
            augmented = np.zeros(
                (len(free_vertices) + 1, len(free_vertices) + 1),
                dtype=np.int32,
            )
            augmented[:-1, :-1] = free_block
            augmented[:-1, -1] = h
            augmented[-1, :-1] = h
            augmented_order = len(augmented)
            cross_state_order_histogram[augmented_order] = (
                cross_state_order_histogram.get(augmented_order, 0) + 1
            )
            target_sign = 1 if anchor_energy > 0 else -1
            exact_cap: Optional[int] = None
            if include_exact_face:
                exact_cap = exact_one_sided_caps(augmented)[
                    0 if target_sign > 0 else 1
                ]
            cross_base_cache[base_key] = {
                "free_block": free_block,
                "h": h,
                "augmented": augmented,
                "target_sign": target_sign,
                "exact_cap": exact_cap,
                "projector_cap": weighted_projector_floor(
                    augmented, target_sign
                ),
            }
        base = cross_base_cache[base_key]
        witness_key = (base_key, free_key)
        if witness_key in cross_witness_cache:
            return {**base, **cross_witness_cache[witness_key]}
        free_block = np.asarray(base["free_block"])
        h = np.asarray(base["h"])
        augmented = np.asarray(base["augmented"])
        target_sign = int(base["target_sign"])
        initial_cap, greedy_cap, greedy_flips = greedy_aligned_floor(
            augmented, target_sign, free_values
        )
        exact_cap = base["exact_cap"]
        if exact_cap is not None and greedy_cap > exact_cap:
            raise AssertionError((exact_cap, greedy_cap))
        r = np.where(h > 0, 1, np.where(h < 0, -1, free_values)).astype(
            np.int32
        )
        free_energy = int(r @ free_block @ r) if len(free_vertices) else 0
        cheap_bound = abs(anchor_energy + free_energy) + 2 * int(np.abs(h).sum())
        witness = {
            "initial_aligned_cap": initial_cap,
            "greedy_cap": greedy_cap,
            "greedy_flips": greedy_flips,
            "cheap_bound": float(cheap_bound),
        }
        cross_witness_cache[witness_key] = witness
        greedy_cached_flip_counts.append(greedy_flips)
        return {**base, **witness}

    cross_exact_certificate = np.zeros(len(x), dtype=np.float64)
    cross_projector_certificate = np.zeros(len(x), dtype=np.float64)
    cross_initial_certificate = np.zeros(len(x), dtype=np.float64)
    cross_greedy_certificate = np.zeros(len(x), dtype=np.float64)
    cross_cheap_certificate = np.zeros(len(x), dtype=np.float64)
    for row_index in np.flatnonzero(opposite):
        indicator = agreement[row_index]
        anchor_i = cross_channel(
            indicator, x[row_index], int(energy_i[row_index]), "I"
        )
        anchor_j = cross_channel(
            ~indicator, x[row_index], int(energy_j[row_index]), "J"
        )
        if include_exact_face:
            assert anchor_i["exact_cap"] is not None and anchor_j["exact_cap"] is not None
            cross_exact_certificate[row_index] = max(
                abs(int(energy_i[row_index])) + float(anchor_i["exact_cap"]),
                abs(int(energy_j[row_index])) + float(anchor_j["exact_cap"]),
            )
        cross_projector_certificate[row_index] = max(
            abs(int(energy_i[row_index])) + float(anchor_i["projector_cap"]),
            abs(int(energy_j[row_index])) + float(anchor_j["projector_cap"]),
        )
        cross_initial_certificate[row_index] = max(
            abs(int(energy_i[row_index]))
            + float(anchor_i["initial_aligned_cap"]),
            abs(int(energy_j[row_index]))
            + float(anchor_j["initial_aligned_cap"]),
        )
        cross_greedy_certificate[row_index] = max(
            abs(int(energy_i[row_index])) + float(anchor_i["greedy_cap"]),
            abs(int(energy_j[row_index])) + float(anchor_j["greedy_cap"]),
        )
        cross_cheap_certificate[row_index] = max(
            float(anchor_i["cheap_bound"]), float(anchor_j["cheap_bound"])
        )

    # Same-sign shores are already losslessly recoupled at the response value.
    cross_exact_certificate[~opposite] = response[~opposite]
    cross_projector_certificate[~opposite] = response[~opposite]
    cross_initial_certificate[~opposite] = response[~opposite]
    cross_greedy_certificate[~opposite] = response[~opposite]
    cross_cheap_certificate[~opposite] = response[~opposite]
    cross_exact_defect = (
        np.maximum(response - cross_exact_certificate, 0.0)
        if include_exact_face
        else None
    )
    cross_projector_defect = np.maximum(
        response - cross_projector_certificate, 0.0
    )
    cross_initial_defect = np.maximum(response - cross_initial_certificate, 0.0)
    cross_greedy_defect = np.maximum(response - cross_greedy_certificate, 0.0)
    cross_cheap_defect = np.maximum(response - cross_cheap_certificate, 0.0)
    if include_exact_face:
        for name, candidate_defect in (
            ("weighted projector", cross_projector_defect),
            ("initial aligned", cross_initial_defect),
            ("greedy", cross_greedy_defect),
        ):
            if np.any(candidate_defect < cross_exact_defect - 1e-8):
                raise AssertionError(f"{name} witness exceeded exact augmented cap")

    exact_response = exact_response_fraction(n)
    response_mean, response_se = mean_standard_error(response, exact)
    if exact and Fraction(int(response.sum()), len(response)) != exact_response:
        raise AssertionError((label, response_mean, float(exact_response)))

    def summarize_defect(values: np.ndarray, name: str) -> Dict[str, object]:
        mean, se = mean_standard_error(values, exact)
        certificate = float(exact_response) - mean
        answer = {
            "name": name,
            "mean": mean,
            "standard_error": se,
            "mean_over_n_to_three_halves": mean / n ** 1.5,
            "maximum": float(np.max(values)),
            "positive_fraction": float(np.mean(values > 1e-10)),
            "quantiles": {
                "0.5": float(np.quantile(values, 0.5)),
                "0.9": float(np.quantile(values, 0.9)),
                "0.99": float(np.quantile(values, 0.99)),
            },
            "Q_certificate": certificate,
            "Q_certificate_over_n_to_three_halves": certificate / n ** 1.5,
            "exceeds_current_doubled_constant": (
                certificate / n ** 1.5 > CURRENT_DOUBLED_CONSTANT
            ),
        }
        if exact and np.allclose(values, np.rint(values), atol=1e-10):
            integer_values = np.rint(values).astype(np.int64)
            total = int(integer_values.sum())
            answer["exact_mean"] = fraction_payload(Fraction(total, len(values)))
            distinct, counts = np.unique(integer_values, return_counts=True)
            answer["exact_integer_histogram"] = {
                str(int(value)): int(count)
                for value, count in zip(distinct, counts)
            }
        return answer

    projector_summary = summarize_defect(projector_delta, "projector")
    nuclear_summary = summarize_defect(nuclear_delta, "nuclear")
    exact_summary = (
        summarize_defect(exact_delta, "exact_face")
        if exact_delta is not None
        else None
    )
    cross_exact_summary = (
        summarize_defect(cross_exact_defect, "cross_augmented_exact_face")
        if cross_exact_defect is not None
        else None
    )
    cross_projector_summary = summarize_defect(
        cross_projector_defect, "cross_augmented_weighted_projector"
    )
    cross_initial_summary = summarize_defect(
        cross_initial_defect, "cross_augmented_initial_aligned_sign_h"
    )
    cross_greedy_summary = summarize_defect(
        cross_greedy_defect, "cross_augmented_greedy_coordinate_ascent"
    )
    cross_cheap_summary = summarize_defect(
        cross_cheap_defect, "cross_aware_sign_h_witness"
    )

    quadratic_cap = recorded_quadratic_cap
    cap_classification = "recorded"
    if independently_enumerate_cap:
        quadratic_cap = exact_quadratic_cap(matrix.astype(np.int32))
        cap_classification = "exact exhaustive recomputation"
        if recorded_quadratic_cap is not None and quadratic_cap != recorded_quadratic_cap:
            raise AssertionError((label, quadratic_cap, recorded_quadratic_cap))

    if quadratic_cap is not None:
        samplewise_certificates = {
            "old_projector": response - projector_delta,
            "old_nuclear": response - nuclear_delta,
            "augmented_projector": response - cross_projector_defect,
            "augmented_initial": response - cross_initial_defect,
            "augmented_greedy": response - cross_greedy_defect,
            "cheap_sign_h": response - cross_cheap_defect,
        }
        if exact_delta is not None:
            samplewise_certificates["old_exact"] = response - exact_delta
        if cross_exact_defect is not None:
            samplewise_certificates["augmented_exact"] = (
                response - cross_exact_defect
            )
        for certificate_name, certificate_values in samplewise_certificates.items():
            if float(np.max(certificate_values)) > quadratic_cap + 1e-7:
                raise AssertionError(
                    (label, certificate_name, quadratic_cap, np.max(certificate_values))
                )

    for summary in (
        exact_summary,
        projector_summary,
        nuclear_summary,
        cross_exact_summary,
        cross_projector_summary,
        cross_initial_summary,
        cross_greedy_summary,
        cross_cheap_summary,
    ):
        if summary is not None and quadratic_cap is not None:
            # Numerical spectral floors can overshoot only by roundoff.
            if float(summary["Q_certificate"]) > quadratic_cap + 1e-7:
                raise AssertionError((label, summary["name"], quadratic_cap, summary))

    shore_histogram: Dict[str, int] = {}
    for shore_sizes, count in zip(shore_size_by_type, multiplicities):
        key = f"{int(shore_sizes[0])},{int(shore_sizes[1])}"
        shore_histogram[key] = shore_histogram.get(key, 0) + int(count)

    zero_field_count = int(np.count_nonzero(fields == 0))
    output = {
        "label": label,
        "n": n,
        "source": source,
        "source_key": source_key,
        "matrix_sha256": matrix_hash(matrix),
        "law": "X uniform projectively; Y_i=sign((AX)_i), with Y_i=X_i at zero",
        "mode": mode,
        "seed": None if exact else seed,
        "population_size": population,
        "evaluated_X_count": len(x),
        "sampling_with_replacement": not exact,
        "include_exact_face": include_exact_face,
        "quadratic_cap_Q": quadratic_cap,
        "quadratic_cap_classification": cap_classification if quadratic_cap is not None else None,
        "exact_expected_bilinear_response": fraction_payload(exact_response),
        "sample_or_enumerated_response_mean": response_mean,
        "sample_response_standard_error": response_se,
        "response_over_n_to_three_halves": float(exact_response) / n ** 1.5,
        "zero_field_coordinate_fraction": zero_field_count / fields.size,
        "opposite_shore_sign_fraction": float(np.mean(opposite)),
        "unique_agreement_subsets": len(unique_packed),
        "cached_shore_subsets": len(cache),
        "cached_cross_anchor_states": len(cross_base_cache),
        "cached_cross_tied_witness_states": len(cross_witness_cache),
        "cross_augmented_state_order_histogram": {
            str(order): count
            for order, count in sorted(cross_state_order_histogram.items())
        },
        "cross_augmented_largest_state_order": (
            max(cross_state_order_histogram) if cross_state_order_histogram else 0
        ),
        "cross_augmented_cached_boolean_candidates": (
            sum(
                count * (1 << max(order - 1, 0))
                for order, count in cross_state_order_histogram.items()
            )
            if include_exact_face
            else None
        ),
        "cross_augmented_greedy_cached_flip_count": {
            "maximum": max(greedy_cached_flip_counts) if greedy_cached_flip_counts else 0,
            "mean": (
                float(np.mean(greedy_cached_flip_counts))
                if greedy_cached_flip_counts
                else 0.0
            ),
        },
        "agreement_size_histogram": shore_histogram,
        "exact_face": exact_summary,
        "projector": projector_summary,
        "nuclear": nuclear_summary,
        "cross_augmented_exact_face": cross_exact_summary,
        "cross_augmented_weighted_projector": cross_projector_summary,
        "cross_augmented_initial_aligned_sign_h": cross_initial_summary,
        "cross_augmented_greedy_coordinate_ascent": cross_greedy_summary,
        "cross_aware_sign_h_witness": cross_cheap_summary,
    }
    return output


def load_case(path_text: str, key: str) -> np.ndarray:
    payload = json.loads((ROOT / path_text).read_text())
    return np.asarray(nested(payload, key), dtype=np.int32)


def paley_conference(prime: int) -> np.ndarray:
    """Symmetric Paley conference signing of order prime+1, prime=1 mod 4."""
    if prime % 4 != 1:
        raise ValueError(prime)
    matrix = np.zeros((prime + 1, prime + 1), dtype=np.int32)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for left in range(prime):
        for right in range(prime):
            if left != right:
                residue = (left - right) % prime
                matrix[left + 1, right + 1] = (
                    1
                    if pow(residue, (prime - 1) // 2, prime) == 1
                    else -1
                )
    if not np.array_equal(matrix @ matrix, prime * np.eye(prime + 1, dtype=np.int32)):
        raise AssertionError((prime, "conference identity"))
    return matrix


def random_signing(order: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    upper = np.triu(
        rng.choice(np.asarray((-1, 1), dtype=np.int32), size=(order, order)),
        1,
    )
    return upper + upper.T


def case_specifications() -> List[Dict[str, object]]:
    cases: List[Dict[str, object]] = []
    for n in range(3, 9):
        path = f"computations/results/m{n}_minimizer_orbits.json"
        payload = json.loads((ROOT / path).read_text())
        for index, row in enumerate(payload["classes"]):
            cases.append(
                {
                    "label": f"min-n{n}-orbit{row['class']}",
                    "matrix": np.asarray(row["representative_matrix"], dtype=np.int32),
                    "source": path,
                    "source_key": f"classes.{index}.representative_matrix",
                    "mode": "exact",
                    "sample_count": None,
                    "seed": 1000 + n * 10 + index,
                    "include_exact_face": True,
                    "independently_enumerate_cap": True,
                    "recorded_quadratic_cap": 2 * int(payload["target_cap"]),
                }
            )

    cases.extend(
        [
            {
                "label": "exact-n10",
                "path": "computations/results/exact_m10.json",
                "key": "matrix",
                "mode": "exact",
                "sample_count": None,
                "seed": 1100,
                "include_exact_face": True,
                "independently_enumerate_cap": True,
                "recorded_quadratic_cap": 26,
            },
            {
                "label": "witness-n11",
                "path": "computations/results/nested_10_in_11_cap17.json",
                "key": "matrix",
                "mode": "exact",
                "sample_count": None,
                "seed": 1111,
                "include_exact_face": True,
                "independently_enumerate_cap": True,
                "recorded_quadratic_cap": 34,
            },
            {
                "label": "witness-n12",
                "path": "computations/results/extension_nested_m11_to_12.json",
                "key": "parent_matrix",
                "mode": "exact",
                "sample_count": None,
                "seed": 1212,
                "include_exact_face": True,
                "independently_enumerate_cap": True,
                "recorded_quadratic_cap": 36,
            },
            {
                "label": "witness-n13",
                "path": "computations/results/bridge_6_7_sign1_cap20.json",
                "key": "parent_matrix",
                "mode": "exact",
                "sample_count": None,
                "seed": 1313,
                "include_exact_face": True,
                "independently_enumerate_cap": True,
                "recorded_quadratic_cap": 40,
            },
            {
                "label": "witness-n14",
                "path": "computations/results/heuristic_m14_from_conference.json",
                "key": "matrix",
                "mode": "exact",
                "sample_count": None,
                "seed": 1414,
                "include_exact_face": True,
                "independently_enumerate_cap": True,
                "recorded_quadratic_cap": 42,
            },
        ]
    )

    conference = [
        ("conference-n6", "conference_double_p5.json", 6, "exact", None, True, 10),
        ("conference-n10", "conference_order10_gf9.json", 10, "exact", None, True, 30),
        ("conference-n14", "conference_double_p13.json", 14, "exact", None, True, 42),
        (
            "conference-n18",
            "conference_double_p17.json",
            18,
            "monte_carlo",
            32768,
            True,
            66,
        ),
        (
            "conference-n26",
            "conference_order26_gf25.json",
            26,
            "monte_carlo",
            16384,
            False,
            130,
        ),
        (
            "conference-n98",
            "conference_double_p97.json",
            98,
            "monte_carlo",
            2048,
            False,
            None,
        ),
    ]
    for index, (label, filename, n, mode, count, include_exact, cap) in enumerate(conference):
        cases.append(
            {
                "label": label,
                "path": f"computations/results/{filename}",
                "key": "conference_matrix",
                "mode": mode,
                "sample_count": count,
                "seed": 26081300 + n,
                "include_exact_face": include_exact,
                "independently_enumerate_cap": n <= 18,
                "recorded_quadratic_cap": cap,
            }
        )
    # Intermediate Paley orders make the asymptotic trend falsifiable rather
    # than comparing only n=26 with n=98.  They are generated exactly, while
    # every reported shore-law statistic is reproducible Monte Carlo.
    for prime in (29, 37, 41, 53, 61, 73, 89):
        n = prime + 1
        cases.append(
            {
                "label": f"paley-conference-n{n}",
                "matrix": paley_conference(prime),
                "source": "generated by paley_conference in this auditor",
                "source_key": f"prime={prime}",
                "mode": "monte_carlo",
                "sample_count": 4096,
                "seed": 26081300 + n,
                "include_exact_face": False,
                "independently_enumerate_cap": False,
                "recorded_quadratic_cap": None,
            }
        )
    # Two post-selected finite stress cases record the worst mean greedy
    # defect in the explicit seed sweeps 9300:9308 (n=12) and 9400:9412
    # (n=14).  Larger random cases are held out fixed seeds, not selected.
    for n, seed, label in (
        (12, 9304, "random-sweep-worst-n12-seed9304"),
        (14, 9407, "random-sweep-worst-n14-seed9407"),
    ):
        cases.append(
            {
                "label": label,
                "matrix": random_signing(n, seed),
                "source": "generated by random_signing in this auditor",
                "source_key": f"numpy-default_rng-seed={seed}",
                "mode": "exact",
                "sample_count": None,
                "seed": seed,
                "include_exact_face": True,
                "independently_enumerate_cap": True,
                "recorded_quadratic_cap": None,
            }
        )
    for n, seed, count in (
        (30, 271830, 4096),
        (62, 271862, 4096),
        (98, 271898, 2048),
    ):
        cases.append(
            {
                "label": f"random-heldout-n{n}-seed{seed}",
                "matrix": random_signing(n, seed),
                "source": "generated by random_signing in this auditor",
                "source_key": f"numpy-default_rng-seed={seed}",
                "mode": "monte_carlo",
                "sample_count": count,
                "seed": seed + 1000000,
                "include_exact_face": False,
                "independently_enumerate_cap": False,
                "recorded_quadratic_cap": None,
            }
        )
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--only", action="append", default=[])
    args = parser.parse_args()

    records = []
    for spec in case_specifications():
        label = str(spec["label"])
        if args.only and label not in args.only:
            continue
        if "matrix" in spec:
            matrix = np.asarray(spec["matrix"], dtype=np.int32)
            source = str(spec["source"])
            source_key = str(spec["source_key"])
        else:
            source = str(spec["path"])
            source_key = str(spec["key"])
            matrix = load_case(source, source_key)
        print(f"auditing {label} (n={len(matrix)}, mode={spec['mode']})", flush=True)
        records.append(
            audit_law(
                matrix=matrix,
                label=label,
                source=source,
                source_key=source_key,
                mode=str(spec["mode"]),
                sample_count=spec["sample_count"],
                seed=int(spec["seed"]),
                include_exact_face=bool(spec["include_exact_face"]),
                independently_enumerate_cap=bool(spec["independently_enumerate_cap"]),
                recorded_quadratic_cap=spec["recorded_quadratic_cap"],
            )
        )

    output = {
        "schema": "quadratic-signing-row-sign-recoupling-law-audit-v1",
        "classification": (
            "exact finite enumeration through n=14; exact response and cap but "
            "sampled shore defects at n=18; reproducible Monte Carlo polynomial "
            "shore diagnostics on conference and held-out random families through "
            "n=98; deterministic augmented coordinate ascent is a rigorous "
            "samplewise certificate, while Monte Carlo averages and numerical "
            "spectral values are not asymptotic proofs"
        ),
        "normalization": "doubled Q(A)=max_z |z^T A z|",
        "current_doubled_project_constant": CURRENT_DOUBLED_CONSTANT,
        "zero_field_tie_rule": "Y_i=X_i when (AX)_i=0",
        "cross_field_zero_tie_rule": (
            "use the corresponding original free-shore X coordinate; this "
            "makes augmented initialization switching covariant"
        ),
        "greedy_algorithm": (
            "initialize w=(sign(sigma*h),1), using the declared zero tie; "
            "repeatedly flip the smallest-index coordinate among those with "
            "maximum positive gain -4*sigma*w_i*(Ew)_i"
        ),
        "universal_response_identity": (
            "E[X^T A sign(AX)]=n E|S_{n-1}| for every signing A; "
            "S_m is a length-m Rademacher sum"
        ),
        "records": records,
    }
    canonical = json.dumps(output, sort_keys=True, separators=(",", ":"))
    output["canonical_payload_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "records": len(records),
                "canonical_payload_sha256": output["canonical_payload_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
