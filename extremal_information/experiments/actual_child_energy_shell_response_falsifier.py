#!/usr/bin/env python3
"""Actual finite falsifier for scalar energy data versus shell response.

The two order-eight children are the certified pressure-minimizer classes at
raw temperature t >= 3.  They have the same exact signed energy histogram and
cap.  Attach the unique order-two child, enumerate the complete 8 x 2 bridge
cube at t=3, and compare the combined-energy marginal of the averaged forward
posterior under the negative disorder escort q proportional p^(-1).

The actual-t calculations use long double.  An additional exact rational-
function certificate evaluates the shell responses at z=2, where z=e^t.
Nonidentity there proves that the two response functions in Q(z) are not the
same.  Since e^3 is transcendental, their values at the actual t=3 cannot be
equal.  The KL values themselves are reported numerically.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_radial_ceiling_witness as witness  # noqa: E402


RAW_T = 3.0
DISORDER_EXPONENT = -1.0
EPSILON = 1


def children() -> tuple[list[np.ndarray], np.ndarray]:
    left = [
        np.asarray(witness.A0, dtype=np.int8),
        np.asarray(witness.A1, dtype=np.int8),
    ]
    right = np.asarray(((0, 1), (1, 0)), dtype=np.int8)
    return left, right


def latent_pairs(
    left: np.ndarray, right: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = exact.projective_spins(len(left)).astype(np.int16)
    y = exact.projective_spins(len(right)).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
    words: list[np.ndarray] = []
    shells: list[int] = []
    for xx, hx in zip(x, ex):
        for yy, hy in zip(y, ey):
            words.append((xx[:, None] * yy[None, :]).reshape(-1))
            shells.append(int(hx + EPSILON * hy))
    return np.asarray(words, dtype=np.int16), np.asarray(shells), ex


def bridges(d: int) -> np.ndarray:
    indices = np.arange(1 << d, dtype=np.uint64)
    bits = (
        (indices[:, None] >> np.arange(d, dtype=np.uint64)) & 1
    ).astype(np.int8)
    return (1 - 2 * bits).astype(np.int16)


def actual_shell_response(
    words: np.ndarray, shells: np.ndarray, bridge: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    shell_values = np.unique(shells)
    shell_id = np.searchsorted(shell_values, shells)
    log_weight = np.log(np.cosh(np.longdouble(RAW_T) * shells))
    log_weight -= np.max(log_weight)
    prior_word = np.exp(log_weight)
    prior_word /= np.sum(prior_word)
    prior_shell = np.bincount(
        shell_id, weights=np.asarray(prior_word, dtype=np.float64),
        minlength=len(shell_values),
    )

    output = np.empty(len(bridge), dtype=np.longdouble)
    posterior_shell = np.empty(
        (len(bridge), len(shell_values)), dtype=np.longdouble
    )
    for start in range(0, len(bridge), 2048):
        block = bridge[start : start + 2048]
        correlation = block @ words.T
        terms = np.cosh(np.longdouble(RAW_T) * correlation)
        terms *= prior_word[None, :]
        denominator = np.sum(terms, axis=1)
        output[start : start + len(block)] = denominator
        for shell in range(len(shell_values)):
            posterior_shell[start : start + len(block), shell] = (
                np.sum(terms[:, shell_id == shell], axis=1) / denominator
            )

    escort = output ** np.longdouble(DISORDER_EXPONENT)
    escort /= np.sum(escort)
    averaged = np.asarray(escort @ posterior_shell, dtype=np.float64)
    kl = float(np.dot(averaged, np.log(averaged / prior_shell)))
    return shell_values, prior_shell, averaged, kl


def integer_cosh_at_two(k: int) -> int:
    """Common-scale integer for cosh(k log 2), valid for 0 <= k <= 16."""

    return (1 << (16 + k)) + (1 << (16 - k))


def exact_response_at_z_two(
    words: np.ndarray, shells: np.ndarray, bridge: np.ndarray
) -> tuple[np.ndarray, list[Fraction], int]:
    """Return the exact a=-1 shell response after setting z=e^t=2."""

    shell_values = np.unique(shells)
    shell_id = np.searchsorted(shell_values, shells)
    correlation = np.abs(bridge @ words.T)
    vectors = np.zeros((len(bridge), len(shell_values)), dtype=object)
    for shell, energy in enumerate(shell_values):
        latent = np.flatnonzero(shell_id == shell)
        radial = np.zeros(len(bridge), dtype=np.uint64)
        for value in range(0, words.shape[1] + 1, 2):
            radial += (
                np.count_nonzero(correlation[:, latent] == value, axis=1)
                .astype(np.uint64)
                * np.uint64(integer_cosh_at_two(value))
            )
        coefficient = integer_cosh_at_two(abs(int(energy)))
        vectors[:, shell] = [coefficient * int(item) for item in radial]

    grouped: dict[tuple[int, ...], int] = {}
    for row in vectors:
        key = tuple(int(item) for item in row)
        grouped[key] = grouped.get(key, 0) + 1

    normalizer = Fraction(0)
    numerator = [Fraction(0) for _ in shell_values]
    for row, multiplicity in grouped.items():
        total = sum(row)
        normalizer += Fraction(multiplicity, total)
        for shell, value in enumerate(row):
            numerator[shell] += Fraction(multiplicity * value, total * total)
    return shell_values, [item / normalizer for item in numerator], len(grouped)


def digest_fraction(value: Fraction) -> dict[str, object]:
    payload = f"{value.numerator}/{value.denominator}".encode()
    return {
        "sign": (value > 0) - (value < 0),
        "numerator_decimal_digits": len(str(abs(value.numerator))),
        "denominator_decimal_digits": len(str(value.denominator)),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def main() -> None:
    left_children, right = children()
    bridge = bridges(len(left_children[0]) * len(right))
    records: list[dict[str, object]] = []
    exact_responses: list[list[Fraction]] = []
    common_histogram: dict[int, int] | None = None
    for class_id, left in enumerate(left_children):
        words, shells, energies = latent_pairs(left, right)
        histogram = dict(sorted(Counter(map(int, energies)).items()))
        if common_histogram is None:
            common_histogram = histogram
        else:
            assert histogram == common_histogram
        assert max(map(abs, histogram)) == 10

        shell_values, prior, averaged, kl = actual_shell_response(
            words, shells, bridge
        )
        exact_values, exact_response, group_count = exact_response_at_z_two(
            words, shells, bridge
        )
        assert np.array_equal(shell_values, exact_values)
        exact_responses.append(exact_response)
        records.append(
            {
                "class_id": class_id,
                "shell_values": list(map(int, shell_values)),
                "prior_shell": list(map(float, prior)),
                "negative_path_averaged_shell": list(map(float, averaged)),
                "shell_kl": kl,
                "z_equals_two_exact_group_count": group_count,
                "z_equals_two_shell_minus_11": float(
                    exact_response[int(np.where(shell_values == -11)[0][0])]
                ),
            }
        )

    assert records[0]["prior_shell"] == records[1]["prior_shell"]
    shell_index = records[0]["shell_values"].index(-11)
    exact_difference = (
        exact_responses[0][shell_index] - exact_responses[1][shell_index]
    )
    assert exact_difference != 0

    result = {
        "status": "actual finite scalar-data falsifier",
        "raw_t": RAW_T,
        "scaled_beta_for_parent_order_10": RAW_T * math.sqrt(10.0),
        "negative_disorder_exponent": DISORDER_EXPONENT,
        "orientation": EPSILON,
        "left_order": 8,
        "right_order": 2,
        "bridge_cube_size": len(bridge),
        "common_exact_signed_energy_histogram": common_histogram,
        "records": records,
        "exact_rational_function_nonidentity_at_z_equals_two": {
            "shell": -11,
            **digest_fraction(exact_difference),
        },
        "scope": (
            "The exact certificate proves that the shell-response rational "
            "functions differ.  The t=3 KL values are numerical.  This is a "
            "finite actual-minimizer falsifier, not a scalable small-t theorem."
        ),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
