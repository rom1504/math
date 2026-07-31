#!/usr/bin/env python3
"""Exact finite audit for the edge-twisted two-Hadamard multiplier."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


H = np.asarray(
    [[1, 1, 1, 1], [1, -1, 1, -1],
     [1, 1, -1, -1], [1, -1, -1, 1]],
    dtype=np.int64,
)
D = np.diag(np.asarray([-1, -1, -1, 1], dtype=np.int64))
K = D @ H @ D


def bareiss_determinant(matrix: np.ndarray) -> int:
    values = [[int(entry) for entry in row] for row in matrix]
    order = len(values)
    previous = 1
    sign = 1
    for k in range(order - 1):
        if values[k][k] == 0:
            swap = next((i for i in range(k + 1, order) if values[i][k]), None)
            if swap is None:
                return 0
            values[k], values[swap] = values[swap], values[k]
            sign *= -1
        pivot = values[k][k]
        for i in range(k + 1, order):
            for j in range(k + 1, order):
                values[i][j] = (
                    values[i][j] * pivot - values[i][k] * values[k][j]
                ) // previous
        previous = pivot
    return sign * values[-1][-1]


def characteristic_coefficients(matrix: np.ndarray) -> list[int]:
    """Return det(lambda*I-matrix), coefficients in descending order."""
    order = len(matrix)
    total = [0] * (order + 1)  # ascending powers
    for permutation in itertools.permutations(range(order)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(order) for j in range(i + 1, order)
        )
        polynomial = [1]
        for i, j in enumerate(permutation):
            factor = [-int(matrix[i, j]), 1] if i == j else [-int(matrix[i, j])]
            product = [0] * (len(polynomial) + len(factor) - 1)
            for left, a in enumerate(polynomial):
                for right, b in enumerate(factor):
                    product[left + right] += a * b
            polynomial = product
        sign = -1 if inversions % 2 else 1
        for power, coefficient in enumerate(polynomial):
            total[power] += sign * coefficient
    return list(reversed(total))


def boolean_extremal_vectors(matrix: np.ndarray) -> list[list[int]]:
    result = []
    for values in itertools.product((-1, 1), repeat=len(matrix)):
        vector = np.asarray(values, dtype=np.int64)
        image = matrix @ vector
        root = int(round(len(matrix) ** 0.5))
        if np.array_equal(image, root * vector) or np.array_equal(
            image, -root * vector
        ):
            result.append(vector.tolist())
    return result


def boolean_transport_domain(matrix: np.ndarray) -> list[list[int]]:
    result = []
    for values in itertools.product((-1, 1), repeat=len(matrix)):
        vector = np.asarray(values, dtype=np.int64)
        image = matrix @ vector
        if np.all(np.abs(image) == 2):
            result.append(vector.tolist())
    return result


def triangle_lift() -> np.ndarray:
    """All-negative triangle, b=(0,1,1), tau_ij=b_i*b_j."""
    generators = (H, K)
    bits = (0, 1, 1)
    matrix = np.zeros((12, 12), dtype=np.int64)
    for i, bit in enumerate(bits):
        block = generators[bit].copy()
        np.fill_diagonal(block, 0)
        matrix[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] = block
    for i in range(3):
        for j in range(i + 1, 3):
            block = -generators[bits[i] * bits[j]]
            matrix[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] = block
            matrix[4 * j : 4 * j + 4, 4 * i : 4 * i + 4] = block.T
    return matrix


def build_lift(base: np.ndarray, bits: tuple[int, ...]) -> np.ndarray:
    generators = (H, K)
    order = len(base)
    matrix = np.zeros((4 * order, 4 * order), dtype=np.int64)
    for i, bit in enumerate(bits):
        block = generators[bit].copy()
        np.fill_diagonal(block, 0)
        matrix[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] = block
    for i in range(order):
        for j in range(i + 1, order):
            block = int(base[i, j]) * generators[bits[i] * bits[j]]
            matrix[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] = block
            matrix[4 * j : 4 * j + 4, 4 * i : 4 * i + 4] = block.T
    return matrix


def root_gauged_signing(order: int, code: int) -> np.ndarray:
    matrix = np.ones((order, order), dtype=np.int64)
    np.fill_diagonal(matrix, 0)
    bit = 0
    for i in range(1, order):
        for j in range(i + 1, order):
            if code & (1 << bit):
                matrix[i, j] = matrix[j, i] = -1
            bit += 1
    return matrix


def vectorized_cap(matrix: np.ndarray) -> int:
    order = len(matrix)
    codes = np.arange(1 << (order - 1), dtype=np.uint64)
    bit_positions = np.arange(order - 1, dtype=np.uint64)
    tails = 1 - 2 * ((codes[:, None] >> bit_positions) & 1).astype(np.int8)
    spins = np.column_stack((np.ones(len(codes), dtype=np.int8), tails))
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2
    return int(np.max(np.abs(energies)))


def first_root_gauged_with_cap(order: int, target: int) -> tuple[int, np.ndarray]:
    edge_bits = (order - 1) * (order - 2) // 2
    for code in range(1 << edge_bits):
        matrix = root_gauged_signing(order, code)
        if vectorized_cap(matrix) == target:
            return code, matrix
    raise AssertionError((order, target))


def balanced_split_audit(order: int, target_cap: int) -> dict[str, object]:
    code, base = first_root_gauged_with_cap(order, target_cap)
    one_count = (order + 1) // 2
    records = []
    for ones in itertools.combinations(range(order), one_count):
        bits = tuple(1 if i in ones else 0 for i in range(order))
        lifted_cap = vectorized_cap(build_lift(base, bits))
        records.append({"fiber_bits": list(bits), "lifted_cap": lifted_cap})
    return {
        "order": order,
        "root_gauged_code": code,
        "base_matrix": base.tolist(),
        "base_cap": vectorized_cap(base),
        "balanced_splits": records,
        "lifted_cap_multiset": dict(Counter(row["lifted_cap"] for row in records)),
    }


def spin_from_positive_bit_hex(encoded: str, order: int) -> np.ndarray:
    value = int(encoded, 16)
    return np.asarray(
        [1 if value & (1 << index) else -1 for index in range(order)],
        dtype=np.int64,
    )


def fiber_type_counts(spin: np.ndarray) -> dict[str, int]:
    counts = Counter(
        "".join("+" if value > 0 else "-" for value in spin[start : start + 4])
        for start in range(0, len(spin), 4)
    )
    return dict(sorted(counts.items()))


def projective_spin_from_positive_tail_code(code: int, order: int) -> np.ndarray:
    """First spin +1; a set tail bit means +1 and an unset bit means -1."""
    return np.asarray(
        [1, *(1 if code & (1 << (index - 1)) else -1
              for index in range(1, order))],
        dtype=np.int64,
    )


def energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    return int(spin @ matrix @ spin // 2)


def product_channel_audit(source_path: str) -> dict[str, object]:
    """Audit the recursively aligned v=(---+) product channel exactly."""
    source = json.loads(Path(source_path).read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    order = len(base)
    half = order // 2
    subset = np.arange(half, order)
    restricted = base[np.ix_(subset, subset)]
    records = []
    for code in range(1 << (order - 1)):
        spin = projective_spin_from_positive_tail_code(code, order)
        full_energy = energy(base, spin)
        restricted_energy = energy(restricted, spin[subset])
        limit_coefficient = full_energy - restricted_energy + order // 2
        records.append(
            (abs(limit_coefficient), limit_coefficient, code,
             full_energy, restricted_energy)
        )
    _, limit_coefficient, code, full_energy, restricted_energy = max(records)
    spin = projective_spin_from_positive_tail_code(code, order)
    v = np.asarray([-1, -1, -1, 1], dtype=np.int64)
    first = build_lift(base, tuple(0 if i < half else 1 for i in range(order)))
    first_spin = np.kron(spin, v)
    first_energy = energy(first, first_spin)
    expected_first = 8 * full_energy - 4 * restricted_energy + 3 * order
    if first_energy != expected_first:
        raise AssertionError((source_path, first_energy, expected_first))
    first_subset = np.arange(2 * order, 4 * order)
    first_restricted_energy = energy(
        first[np.ix_(first_subset, first_subset)], first_spin[first_subset]
    )
    expected_first_restricted = 4 * restricted_energy + order
    if first_restricted_energy != expected_first_restricted:
        raise AssertionError(
            (source_path, first_restricted_energy, expected_first_restricted)
        )
    second = build_lift(
        first, tuple(0 if i < 2 * order else 1 for i in range(4 * order))
    )
    second_spin = np.kron(first_spin, v)
    second_energy = energy(second, second_spin)
    expected_second = 64 * full_energy - 48 * restricted_energy + 32 * order
    if second_energy != expected_second:
        raise AssertionError((source_path, second_energy, expected_second))
    return {
        "source": source_path,
        "base_order": order,
        "base_cap": int(source["profile"]["M"]),
        "exhaustive_best_positive_tail_code": code,
        "base_spin": spin.tolist(),
        "base_full_energy_F": full_energy,
        "base_second_half_energy_R0": restricted_energy,
        "limit_coefficient_F_minus_R0_plus_n_over_2": limit_coefficient,
        "limiting_absolute_normalized_energy": (
            abs(limit_coefficient) / order**1.5
        ),
        "depth_1_exact_product_energy": first_energy,
        "depth_1_formula": expected_first,
        "depth_1_second_half_energy": first_restricted_energy,
        "depth_2_exact_product_energy": second_energy,
        "depth_2_formula": expected_second,
    }


def all_balanced_product_partitions_audit(source_path: str) -> dict[str, object]:
    """Exhaust every balanced choice of the K-colored half and every spin."""
    source = json.loads(Path(source_path).read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    order = len(base)
    codes = np.arange(1 << (order - 1), dtype=np.uint64)
    bit_positions = np.arange(order - 1, dtype=np.uint64)
    tails = 2 * ((codes[:, None] >> bit_positions) & 1).astype(np.int8) - 1
    spins = np.column_stack((np.ones(len(codes), dtype=np.int8), tails))
    full = np.einsum("bi,ij,bj->b", spins, base, spins, optimize=True) // 2
    records = []
    for subset_tuple in itertools.combinations(range(order), order // 2):
        subset = np.asarray(subset_tuple, dtype=np.int64)
        restricted = np.einsum(
            "bi,ij,bj->b", spins[:, subset], base[np.ix_(subset, subset)],
            spins[:, subset], optimize=True,
        ) // 2
        limit_values = np.abs(full - restricted + order // 2)
        first_values = np.abs(8 * full - 4 * restricted + 3 * order)
        limit_index = int(np.argmax(limit_values))
        first_index = int(np.argmax(first_values))
        records.append({
            "second_half": list(subset_tuple),
            "maximum_absolute_limit_coefficient": int(limit_values[limit_index]),
            "limit_witness_positive_tail_code": limit_index,
            "maximum_absolute_depth_1_product_energy": int(first_values[first_index]),
            "depth_1_witness_positive_tail_code": first_index,
        })
    minimum_limit = min(
        records,
        key=lambda row: (
            row["maximum_absolute_limit_coefficient"],
            row["maximum_absolute_depth_1_product_energy"],
            row["second_half"],
        ),
    )
    minimum_first = min(
        records,
        key=lambda row: (
            row["maximum_absolute_depth_1_product_energy"],
            row["maximum_absolute_limit_coefficient"],
            row["second_half"],
        ),
    )
    limit_histogram = Counter(
        row["maximum_absolute_limit_coefficient"] for row in records
    )
    first_histogram = Counter(
        row["maximum_absolute_depth_1_product_energy"] for row in records
    )
    return {
        "source": source_path,
        "base_order": order,
        "balanced_partition_count": len(records),
        "minimum_over_partitions_of_worst_limit": minimum_limit,
        "minimum_limiting_absolute_normalized_energy": (
            minimum_limit["maximum_absolute_limit_coefficient"] / order**1.5
        ),
        "minimum_over_partitions_of_worst_depth_1_energy": minimum_first,
        "worst_limit_histogram": dict(sorted(limit_histogram.items())),
        "worst_depth_1_energy_histogram": dict(sorted(first_histogram.items())),
    }


def n14_saved_witness_audit() -> dict[str, object]:
    source_path = "computations/results/heuristic_m14_from_conference.json"
    source = json.loads(Path(source_path).read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    bits = tuple(0 if i < 7 else 1 for i in range(14))
    lifted = build_lift(base, bits)
    encoded = "07810007778967"
    spin = spin_from_positive_bit_hex(encoded, 56)
    exact = energy(lifted, spin)
    if exact != 238:
        raise AssertionError(exact)
    return {
        "source": source_path,
        "order": 56,
        "positive_bit_hex_little_endian": encoded,
        "exact_energy": exact,
        "absolute_normalized_energy": abs(exact) / 56**1.5,
        "fiber_type_counts": fiber_type_counts(spin),
    }


def recursive_n10_witness_audit() -> dict[str, object]:
    source = json.loads(Path("computations/results/exact_m10.json").read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    first_bits = tuple(0 if i < 5 else 1 for i in range(10))
    first = build_lift(base, first_bits)
    first_hex = "e13ae1eeee"
    first_spin = spin_from_positive_bit_hex(first_hex, 40)
    first_energy = int(first_spin @ first @ first_spin // 2)
    if first_energy != -170:
        raise AssertionError(first_energy)

    second_bits = tuple(0 if i < 20 else 1 for i in range(40))
    second = build_lift(first, second_bits)
    second_hex = "7889877787767889788987777888788878887888"
    second_spin = spin_from_positive_bit_hex(second_hex, 160)
    second_energy = int(second_spin @ second @ second_spin // 2)
    if second_energy != 1488:
        raise AssertionError(second_energy)
    return {
        "source": "computations/results/exact_m10.json",
        "base_cap": int(source["profile"]["M"]),
        "first_depth": {
            "order": 40,
            "positive_bit_hex_little_endian": first_hex,
            "exact_energy": first_energy,
            "absolute_normalized_energy": abs(first_energy) / 40**1.5,
            "ratio_to_8_times_base_cap": abs(first_energy) / (8 * 13),
            "fiber_type_counts": fiber_type_counts(first_spin),
        },
        "second_depth": {
            "order": 160,
            "positive_bit_hex_little_endian": second_hex,
            "exact_energy": second_energy,
            "absolute_normalized_energy": abs(second_energy) / 160**1.5,
            "fiber_type_counts": fiber_type_counts(second_spin),
        },
    }


def exact_energy_profile(matrix: np.ndarray) -> dict[str, object]:
    histogram: Counter[int] = Counter()
    witness = None
    cap = -1
    for tail in itertools.product((-1, 1), repeat=len(matrix) - 1):
        spin = np.asarray((1, *tail), dtype=np.int64)
        energy = int(spin @ matrix @ spin // 2)
        histogram[energy] += 1
        if abs(energy) > cap:
            cap = abs(energy)
            witness = spin.tolist()
    return {
        "minimum": min(histogram),
        "maximum": max(histogram),
        "cap": cap,
        "witness": witness,
        "histogram": {str(key): value for key, value in sorted(histogram.items())},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    identity = np.eye(4, dtype=np.int64)
    for matrix in (H, K):
        if not np.array_equal(matrix, matrix.T):
            raise AssertionError("generator is not symmetric")
        if not np.array_equal(matrix @ matrix, 4 * identity):
            raise AssertionError("generator is not Hadamard")
    extremal_h = boolean_extremal_vectors(H)
    extremal_k = boolean_extremal_vectors(K)
    if {tuple(row) for row in extremal_h} & {tuple(row) for row in extremal_k}:
        raise AssertionError("generators have a common Boolean extremal vector")
    domain_h = boolean_transport_domain(H)
    domain_k = boolean_transport_domain(K)
    set_domain_h = {tuple(row) for row in domain_h}
    set_domain_k = {tuple(row) for row in domain_k}
    boolean_cube = set(itertools.product((-1, 1), repeat=4))
    if set_domain_h & set_domain_k or set_domain_h | set_domain_k != boolean_cube:
        raise AssertionError("transport domains do not partition the cube")
    off_domain_bilinear_maxima = {}
    for label, matrix, domain in (("H", H, set_domain_h), ("K", K, set_domain_k)):
        maximum = 0
        for left in boolean_cube:
            for right in boolean_cube:
                if left not in domain or right not in domain:
                    value = abs(int(np.asarray(left) @ matrix @ np.asarray(right)))
                    maximum = max(maximum, value)
        if maximum != 4:
            raise AssertionError((label, maximum))
        off_domain_bilinear_maxima[label] = maximum
    hk = H @ K
    characteristic = characteristic_coefficients(hk)
    if characteristic != [1, -4, 0, -64, 256]:
        raise AssertionError(characteristic)
    plus_boolean = [
        list(values) for values in itertools.product((-1, 1), repeat=4)
        if np.array_equal(hk @ np.asarray(values), 4 * np.asarray(values))
    ]
    minus_boolean = [
        list(values) for values in itertools.product((-1, 1), repeat=4)
        if np.array_equal(hk @ np.asarray(values), -4 * np.asarray(values))
    ]
    if plus_boolean or minus_boolean:
        raise AssertionError("HK has a Boolean vector at eigenvalue plus/minus four")
    determinant_plus = bareiss_determinant(hk + 4 * identity)
    if determinant_plus != 1024:
        raise AssertionError("minus-four exclusion certificate changed")
    signing = triangle_lift()
    off_diagonal = signing[~np.eye(12, dtype=bool)]
    if np.any(np.diag(signing)) or not np.all(np.abs(off_diagonal) == 1):
        raise AssertionError("triangle lift is not a signing")
    profile = exact_energy_profile(signing)
    if profile["cap"] != 20:
        raise AssertionError(profile)
    minimizer_audits = [balanced_split_audit(4, 4), balanced_split_audit(5, 4)]
    if minimizer_audits[0]["lifted_cap_multiset"] != {42: 5, 48: 1}:
        raise AssertionError(minimizer_audits[0])
    if minimizer_audits[1]["lifted_cap_multiset"] != {56: 10}:
        raise AssertionError(minimizer_audits[1])
    recursive_witness = recursive_n10_witness_audit()
    product_channels = [
        product_channel_audit("computations/results/exact_m10.json"),
        product_channel_audit(
            "computations/results/heuristic_m14_from_conference.json"
        ),
    ]
    all_balanced_product_partitions = [
        all_balanced_product_partitions_audit(
            "computations/results/exact_m10.json"
        ),
        all_balanced_product_partitions_audit(
            "computations/results/heuristic_m14_from_conference.json"
        ),
    ]
    n14_saved_witness = n14_saved_witness_audit()
    output = {
        "schema": "quadratic-signing-edge-twisted-hadamard-audit-v1",
        "classification": "exact generator algebra and exhaustive order-12 cap",
        "H": H.tolist(),
        "D_diagonal": np.diag(D).tolist(),
        "K_equals_DHD": K.tolist(),
        "H_boolean_extremal_vectors": extremal_h,
        "K_boolean_extremal_vectors": extremal_k,
        "H_boolean_transport_domain": domain_h,
        "K_boolean_transport_domain": domain_k,
        "transport_domains_partition_boolean_cube": True,
        "off_domain_bilinear_maximum": off_domain_bilinear_maxima,
        "HK": hk.tolist(),
        "HK_characteristic_coefficients": characteristic,
        "HK_characteristic_polynomial": "(lambda-4)^2*(lambda^2+4*lambda+16)",
        "det_HK_plus_4I": determinant_plus,
        "HK_plus_or_minus_4_boolean_vectors": [],
        "triangle": {
            "seed": "all-negative order-three signing",
            "fiber_bits": [0, 1, 1],
            "edge_color_rule": "tau_ij=b_i*b_j",
            "order": 12,
            "matrix": signing.tolist(),
            "exact_energy_profile": profile,
            "normalized_cap": profile["cap"] / 12**1.5,
            "half_threshold_at_order_12": 0.5 * 12**1.5,
        },
        "exact_minimizer_balanced_split_audits": minimizer_audits,
        "recursive_n10_heuristic_search_exact_witnesses": recursive_witness,
        "recursively_aligned_product_channel_audits": product_channels,
        "all_balanced_product_partition_audits": all_balanced_product_partitions,
        "n14_heuristic_search_exact_witness": n14_saved_witness,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
