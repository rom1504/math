#!/usr/bin/env python3
"""Exact Boolean spectral-deficit audit for cyclic two-fiber conferences."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from two_fiber_cyclic_conference import CERTIFICATES, circulant, verify


PALEY_EQUIVALENCES = {
    1: {
        "permutation": [0, 1, 2, 3, 5, 4],
        "switching": [1, 1, 1, 1, 1, -1],
    },
    2: {
        "permutation": [0, 1, 14, 16, 4, 2, 7, 9, 5,
                        13, 3, 10, 11, 15, 6, 17, 8, 12],
        "switching": [1, 1, 1, 1, -1, -1, 1, 1, 1,
                      1, 1, 1, -1, 1, -1, -1, 1, -1],
    },
}


def conference_matrix(k: int) -> np.ndarray:
    data = CERTIFICATES[k]
    verify(k, data["a"], data["c"])
    matrix_a = circulant(data["a"])
    matrix_c = circulant(data["c"])
    return np.block([[matrix_a, matrix_c], [matrix_c.T, -matrix_a]])


def paley_conference(prime: int) -> np.ndarray:
    matrix = np.zeros((prime + 1, prime + 1), dtype=np.int64)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for left in range(prime):
        for right in range(prime):
            if left == right:
                continue
            residue = (left - right) % prime
            matrix[left + 1, right + 1] = (
                1 if pow(residue, (prime - 1) // 2, prime) == 1 else -1
            )
    return matrix


def exact_cap(matrix: np.ndarray) -> tuple[int, int, np.ndarray]:
    order = len(matrix)
    codes = np.arange(1 << (order - 1), dtype=np.uint64)
    positions = np.arange(order - 1, dtype=np.uint64)
    tails = 1 - 2 * ((codes[:, None] >> positions) & 1).astype(np.int8)
    spins = np.column_stack((np.ones(len(codes), dtype=np.int8), tails))
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2
    index = int(np.argmax(np.abs(energies)))
    return abs(int(energies[index])), int(energies[index]), spins[index]


def local_field_record(matrix: np.ndarray, spin: np.ndarray, k: int) -> dict[str, object]:
    energy = int(spin @ matrix @ spin // 2)
    orientation = 1 if energy >= 0 else -1
    fields = orientation * spin * (matrix @ spin)
    penalties = (fields - (2 * k - 1)) * (fields - (2 * k + 1))
    s = 2 * k * k + 1
    deficit = 2 * k * s - abs(energy)
    if np.any(penalties < 0) or int(np.sum(penalties)) != 8 * k * deficit:
        raise AssertionError("local-field deficit identity failed")
    return {
        "exact_energy": energy,
        "absolute_energy": abs(energy),
        "arithmetic_upper_before_parity": 2 * k * s,
        "deficit": deficit,
        "oriented_local_field_histogram": {
            str(key): value for key, value in sorted(Counter(map(int, fields)).items())
        },
        "penalty_sum": int(np.sum(penalties)),
        "positive_bit_hex_little_endian": format(
            sum(1 << index for index, value in enumerate(spin) if value > 0), "x"
        ),
    }


def gauge_delete_negative_graph(matrix: np.ndarray, root: int = 0) -> np.ndarray:
    order = len(matrix)
    switching = np.ones(order, dtype=np.int64)
    other = np.asarray([index for index in range(order) if index != root])
    switching[other] = matrix[root, other]
    gauged = switching[:, None] * matrix * switching[None, :]
    return (gauged[np.ix_(other, other)] < 0).astype(np.int8)


def clique_counts(adjacency: np.ndarray, sizes: tuple[int, ...]) -> dict[str, int]:
    counts = {}
    for size in sizes:
        total = 0
        for vertices in itertools.combinations(range(len(adjacency)), size):
            if all(adjacency[i, j] for i, j in itertools.combinations(vertices, 2)):
                total += 1
        counts[str(size)] = total
    return counts


def periodic_correlation(left: np.ndarray, right: np.ndarray | None = None) -> list[int]:
    if right is None:
        right = left
    return [int(left @ np.roll(right, shift)) for shift in range(len(left))]


def self_indexed_asds_record(k: int) -> dict[str, object]:
    data = CERTIFICATES[k]
    a = np.asarray(data["a"], dtype=np.int64)
    c = np.asarray(data["c"], dtype=np.int64)
    alpha = a.copy()
    alpha[0] = 1
    autocorrelation_sum = np.asarray(periodic_correlation(alpha)) + np.asarray(
        periodic_correlation(c)
    )
    expected = 2 * a
    expected[0] = 2 * len(a)
    if not np.array_equal(autocorrelation_sum, expected):
        raise AssertionError("self-indexed autocorrelation identity failed")
    cross = periodic_correlation(alpha, c)
    cross_symmetric = all(cross[shift] == cross[-shift] for shift in range(len(a)))
    return {
        "k": k,
        "length": len(a),
        "autocorrelation_sum": autocorrelation_sum.tolist(),
        "self_indexed_identity_verified": True,
        "alpha_c_cross_correlation": cross,
        "oqs_cross_symmetry": cross_symmetric,
    }


def gray_four_channel_product(
    alpha: np.ndarray, c: np.ndarray, beta: np.ndarray, d: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    ac = np.multiply.outer(alpha, beta)
    cd = np.multiply.outer(c, d)
    ad = np.multiply.outer(alpha, d)
    cb = np.multiply.outer(c, beta)
    first = (ac - cd - ad - cb) // 2
    second = (ac - cd + ad + cb) // 2
    if not np.all(np.isin(first, (-1, 1))) or not np.all(np.isin(second, (-1, 1))):
        raise AssertionError("Gray four-channel outputs are not Boolean")
    return first, second


def periodic_correlation_2d(array: np.ndarray, shift: tuple[int, int]) -> int:
    return int(np.sum(array * np.roll(array, shift, axis=(0, 1))))


def four_channel_axial_obstruction() -> dict[str, object]:
    pairs = []
    for k in (1, 2):
        data = CERTIFICATES[k]
        a = np.asarray(data["a"], dtype=np.int64)
        alpha = a.copy()
        alpha[0] = 1
        pairs.append((alpha, np.asarray(data["c"], dtype=np.int64)))
    (alpha, c), (beta, d) = pairs
    first, second = gray_four_channel_product(alpha, c, beta, d)
    sums = {}
    for left_shift in range(len(alpha)):
        for right_shift in range(len(beta)):
            value = periodic_correlation_2d(first, (left_shift, right_shift))
            value += periodic_correlation_2d(second, (left_shift, right_shift))
            sums[f"{left_shift},{right_shift}"] = value
            s_left = (
                periodic_correlation(alpha)[left_shift]
                + periodic_correlation(c)[left_shift]
            )
            s_right = (
                periodic_correlation(beta)[right_shift]
                + periodic_correlation(d)[right_shift]
            )
            k_left = (
                periodic_correlation(c, alpha)[left_shift]
                - periodic_correlation(alpha, c)[left_shift]
            )
            k_right = (
                periodic_correlation(d, beta)[right_shift]
                - periodic_correlation(beta, d)[right_shift]
            )
            expected = (s_left * s_right - k_left * k_right) // 2
            if value != expected:
                raise AssertionError("four-channel correlation formula failed")
    axial_value = sums["1,0"]
    if abs(axial_value) != 2 * len(beta):
        raise AssertionError(axial_value)
    return {
        "factor_lengths": [len(alpha), len(beta)],
        "output_shape": list(first.shape),
        "exact_correlation_formula_verified_at_every_shift": True,
        "axial_shift": [1, 0],
        "axial_autocorrelation_sum": axial_value,
        "required_self_indexed_magnitude": 2,
        "scaling_obstruction": "absolute axial magnitude equals 2 times second-factor length",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = []
    for k in (1, 2):
        matrix = conference_matrix(k)
        cap, _, spin = exact_cap(matrix)
        records.append({
            "k": k,
            "classification": "exact exhaustive cap",
            "cap": cap,
            "witness": local_field_record(matrix, spin, k),
        })

    matrix_3 = conference_matrix(3)
    witness_code = int("1e1d4bf05d", 16)
    witness_3 = np.asarray(
        [1 if witness_code & (1 << index) else -1 for index in range(38)],
        dtype=np.int64,
    )
    record_3 = {
        "k": 3,
        "classification": "exact witness lower bound; no exact cap claim",
        "witness": local_field_record(matrix_3, witness_3, 3),
    }
    if record_3["witness"]["exact_energy"] != 109:
        raise AssertionError("saved order-38 witness changed")
    records.append(record_3)

    equivalences = []
    for k, certificate in PALEY_EQUIVALENCES.items():
        matrix = conference_matrix(k)
        paley = paley_conference(4 * k * k + 1)
        permutation = certificate["permutation"]
        switching = np.asarray(certificate["switching"], dtype=np.int64)
        permuted = paley[np.ix_(permutation, permutation)]
        if not np.array_equal(
            matrix, switching[:, None] * permuted * switching[None, :]
        ):
            raise AssertionError("Paley equivalence certificate failed")
        equivalences.append({
            "k": k,
            "paley_prime": 4 * k * k + 1,
            "permutation": permutation,
            "switching": switching.tolist(),
            "verified_exactly": True,
        })

    paley_37 = paley_conference(37)
    nonpaley_cliques = clique_counts(gauge_delete_negative_graph(matrix_3), (4, 5))
    paley_cliques = clique_counts(gauge_delete_negative_graph(paley_37), (4, 5))
    if nonpaley_cliques != {"4": 615, "5": 65}:
        raise AssertionError(nonpaley_cliques)
    if paley_cliques != {"4": 555, "5": 0}:
        raise AssertionError(paley_cliques)

    payload = {
        "schema": "quadratic-signing-two-fiber-boolean-spectral-deficit-v1",
        "classification": (
            "proved universal arithmetic cap inequality, exact finite caps and witness, "
            "and exact equivalence/non-equivalence certificates"
        ),
        "universal_theorem": {
            "order": "N=4*k^2+2=2*s",
            "conference_square": "S^2=(4*k^2+1)I",
            "local_field_identity": (
                "2*k*s-|H|=(1/(8*k))*sum_i "
                "(u_i-(2*k-1))*(u_i-(2*k+1))"
            ),
            "cap_bound_with_parity": "cap(S)<=2*k*s-1",
        },
        "finite_boolean_audits": records,
        "paley_equivalences": equivalences,
        "k3_non_paley_certificate": {
            "invariant": "clique counts in the root-0 gauged negative-edge graph",
            "cyclic_two_fiber": nonpaley_cliques,
            "paley_q37": paley_cliques,
            "note": (
                "Paley switching automorphisms are vertex-transitive, so a root-0 "
                "difference excludes switching/permutation equivalence"
            ),
        },
        "self_indexed_asds_audits": [
            self_indexed_asds_record(k) for k in (1, 2, 3)
        ],
        "gray_four_channel_product_audit": four_channel_axial_obstruction(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
