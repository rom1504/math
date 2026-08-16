#!/usr/bin/env python3
"""Finite checks for geodesic cycle contraction and sharp synchronization.

The exhaustive check covers every spanning support in F_2^w for w <= 3,
every independent subfamily B, the cycle-contraction equivalence, and the
pointwise quotient-section bound.  The sharp check builds the six-channel
bent example at (D,k)=(6,4), verifies every quotient cycle, computes the
whole Cayley metric, and checks the exact distance from linear sections.

No external package or solver is required.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
from collections import Counter, deque
from pathlib import Path


def popcount(x: int) -> int:
    return bin(x).count("1")


def xor_sum(values: list[int] | tuple[int, ...]) -> int:
    total = 0
    for value in values:
        total ^= value
    return total


def rank_binary(vectors: list[int] | tuple[int, ...], dimension: int) -> int:
    pivots = [0] * dimension
    rank = 0
    for vector in vectors:
        x = vector
        while x:
            i = x.bit_length() - 1
            if pivots[i]:
                x ^= pivots[i]
            else:
                pivots[i] = x
                rank += 1
                break
    return rank


def span(vectors: list[int] | tuple[int, ...]) -> set[int]:
    result = {0}
    for vector in vectors:
        result |= {x ^ vector for x in tuple(result)}
    return result


def word_lengths(support: tuple[int, ...], dimension: int) -> list[int]:
    distances = [-1] * (1 << dimension)
    distances[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for s in support:
            y = x ^ s
            if distances[y] < 0:
                distances[y] = distances[x] + 1
                queue.append(y)
    return distances


def basis_weight(x: int, basis: tuple[int, ...]) -> int:
    for mask in range(1 << len(basis)):
        if xor_sum([b for i, b in enumerate(basis) if mask >> i & 1]) == x:
            return popcount(mask)
    raise AssertionError("vector is not in the basis span")


def cycle_condition(support: tuple[int, ...], basis: tuple[int, ...]) -> bool:
    w_span = span(basis)
    outside = tuple(s for s in support if s not in set(basis))
    for mask in range(1 << len(outside)):
        chosen = [s for i, s in enumerate(outside) if mask >> i & 1]
        total = xor_sum(chosen)
        if total in w_span and basis_weight(total, basis) > len(chosen):
            return False
    return True


def quotient_basis_representatives(
    support: tuple[int, ...], basis: tuple[int, ...], dimension: int
) -> tuple[int, ...]:
    selected: list[int] = []
    current_rank = len(basis)
    for s in support:
        new_rank = rank_binary(list(basis) + selected + [s], dimension)
        if new_rank > current_rank:
            selected.append(s)
            current_rank = new_rank
    if current_rank != dimension:
        raise AssertionError("support did not span quotient")
    return tuple(selected)


def check_pointwise_section(
    support: tuple[int, ...], basis: tuple[int, ...], dimension: int
) -> bool:
    reps = quotient_basis_representatives(support, basis, dimension)
    w_span = span(basis)
    for s in support:
        found = False
        for mask in range(1 << len(reps)):
            lift = xor_sum([c for i, c in enumerate(reps) if mask >> i & 1])
            error = s ^ lift
            if error in w_span:
                found = True
                if basis_weight(error, basis) > popcount(mask) + 1:
                    return False
                break
        if not found:
            return False
    return True


def exhaustive_equivalence(max_w: int) -> dict[str, object]:
    result: dict[str, object] = {}
    for w in range(1, max_w + 1):
        universe = tuple(range(1, 1 << w))
        supports = 0
        independent_bases = 0
        geodesics = 0
        section_checks = 0
        size_histogram: Counter[int] = Counter()
        for support_mask in range(1 << len(universe)):
            support = tuple(
                universe[i]
                for i in range(len(universe))
                if support_mask >> i & 1
            )
            if rank_binary(support, w) != w:
                continue
            supports += 1
            distances = word_lengths(support, w)
            for size in range(1, w + 1):
                for basis in itertools.combinations(support, size):
                    if rank_binary(basis, w) != size:
                        continue
                    independent_bases += 1
                    target = xor_sum(basis)
                    shortest = distances[target] == size
                    contracts = cycle_condition(support, basis)
                    if shortest != contracts:
                        raise AssertionError(
                            f"equivalence failed: w={w}, S={support}, B={basis}"
                        )
                    if shortest:
                        geodesics += 1
                        if not check_pointwise_section(support, basis, w):
                            raise AssertionError(
                                f"section bound failed: w={w}, S={support}, B={basis}"
                            )
                        section_checks += 1
                        size_histogram[size] += 1
        result[str(w)] = {
            "spanning_supports": supports,
            "independent_subfamilies_checked": independent_bases,
            "geodesic_subfamilies": geodesics,
            "pointwise_section_checks": section_checks,
            "geodesic_size_histogram": dict(sorted(size_histogram.items())),
        }
    return result


BLOCK_ROWS = (
    (8, 4, 2, 1),
    (12, 8, 1, 3),
    (4, 12, 3, 2),
)


def block_diagonal_rows(rows4: tuple[int, ...], blocks: int) -> tuple[int, ...]:
    rows: list[int] = []
    for block in range(blocks):
        for row in rows4:
            rows.append(row << (4 * block))
    return tuple(rows)


def quadratic_value(q: int, rows: tuple[int, ...]) -> int:
    value = 0
    k = len(rows)
    for i in range(k):
        if not (q >> i & 1):
            continue
        for j in range(i + 1, k):
            if (rows[i] >> j & 1) and (q >> j & 1):
                value ^= 1
    return value


def sharp_f(q: int, forms: tuple[tuple[int, ...], ...]) -> int:
    if q == 0:
        return 0
    c = 1
    value = 0
    for i, rows in enumerate(forms):
        phi = quadratic_value(q, rows)
        value |= phi << (2 * i)
        value |= (phi ^ c) << (2 * i + 1)
    return value


def linear_truth(mask: int, q: int) -> int:
    return popcount(mask & q) & 1


def verify_sharp_example(dimension_w: int = 6, k: int = 4) -> dict[str, object]:
    if dimension_w < 6 or k % 4:
        raise ValueError("sharp construction needs D>=6 and k divisible by four")
    forms = tuple(
        block_diagonal_rows(rows, k // 4) for rows in BLOCK_ROWS
    )
    if any(rank_binary(list(rows), k) != k for rows in forms):
        raise AssertionError("an alternating form is singular")
    if any(
        forms[2][i] != (forms[0][i] ^ forms[1][i]) for i in range(k)
    ):
        raise AssertionError("third polar form is not the sum of the first two")

    nonzero_q = tuple(range(1, 1 << k))
    f_values = {q: sharp_f(q, forms) for q in range(1 << k)}

    zero_sum_cycles = 0
    cycle_size_histogram: Counter[int] = Counter()
    max_cycle_defect = 0
    for mask in range(1 << len(nonzero_q)):
        chosen = [q for i, q in enumerate(nonzero_q) if mask >> i & 1]
        if xor_sum(chosen) != 0:
            continue
        zero_sum_cycles += 1
        defect = xor_sum([f_values[q] for q in chosen])
        defect_weight = popcount(defect)
        if defect_weight > len(chosen):
            raise AssertionError(
                f"cycle contraction failed: R={chosen}, defect={defect:b}"
            )
        max_cycle_defect = max(max_cycle_defect, defect_weight)
        cycle_size_histogram[len(chosen)] += 1

    support = tuple(1 << i for i in range(dimension_w)) + tuple(
        f_values[q] | (q << dimension_w) for q in nonzero_q
    )
    distances = word_lengths(support, dimension_w + k)
    target = (1 << dimension_w) - 1
    if distances[target] != dimension_w or max(distances) != dimension_w:
        raise AssertionError("the constructed geodesic is not diametral")

    coordinate_distances: list[int] = []
    coordinate_linear_masks: list[int] = []
    for output_coordinate in range(6):
        best = len(nonzero_q) + 1
        best_mask = 0
        for linear_mask in range(1 << k):
            distance = sum(
                ((f_values[q] >> output_coordinate) & 1)
                != linear_truth(linear_mask, q)
                for q in nonzero_q
            )
            if distance < best:
                best = distance
                best_mask = linear_mask
        coordinate_distances.append(best)
        coordinate_linear_masks.append(best_mask)
    total_distance = sum(coordinate_distances)
    expected = 3 * ((1 << k) - (1 << (k // 2)) - 1)
    if total_distance != expected:
        raise AssertionError(
            f"linear-section distance {total_distance} != formula {expected}"
        )

    triple_defects = {
        popcount(f_values[x] ^ f_values[y] ^ f_values[x ^ y])
        for x in nonzero_q
        for y in nonzero_q
        if x != y
    }
    if triple_defects != {3}:
        raise AssertionError(f"unexpected triple defects: {triple_defects}")

    synchronized_errors = {
        q: popcount(
            f_values[q]
            ^ sum(
                linear_truth(mask, q) << coordinate
                for coordinate, mask in enumerate(coordinate_linear_masks)
            )
        )
        for q in range(1 << k)
    }
    if max(synchronized_errors.values()) > 9:
        raise AssertionError("uniform BLR radius-nine bound failed")

    # Exhaustive roots for deterministic random raw futures.  This checks
    # the all-context comparison on finite samples, while the draft proof is
    # uniform over every future.
    linear_support = tuple(1 << i for i in range(dimension_w)) + tuple(
        (
            sum(
                linear_truth(mask, q) << coordinate
                for coordinate, mask in enumerate(coordinate_linear_masks)
            )
            | (q << dimension_w)
        )
        for q in nonzero_q
    )
    rng = random.Random(20260816)
    future_checks = []
    universe = tuple(range(1, 1 << (dimension_w + k)))
    for probability in (0.0, 0.01, 0.03, 0.08, 0.2):
        for _ in range(8):
            future = tuple(x for x in universe if rng.random() < probability)
            source_distances = word_lengths(tuple(set(support) | set(future)), dimension_w + k)
            linear_distances = word_lengths(
                tuple(set(linear_support) | set(future)), dimension_w + k
            )
            forward = max(b - a for a, b in zip(source_distances, linear_distances))
            reverse = max(a - b for a, b in zip(source_distances, linear_distances))
            if forward > 11 or reverse > 9:
                raise AssertionError("sampled raw-context comparison failed")
            future_checks.append((forward, reverse))

    return {
        "D": dimension_w,
        "k": k,
        "support_size": len(support),
        "cayley_diameter": max(distances),
        "diametral_target_distance": distances[target],
        "zero_sum_cycles_checked": zero_sum_cycles,
        "cycle_size_histogram": dict(sorted(cycle_size_histogram.items())),
        "maximum_cycle_defect_weight": max_cycle_defect,
        "generic_triple_defect_weights": sorted(triple_defects),
        "coordinate_distances_to_linear": coordinate_distances,
        "coordinate_linear_masks": coordinate_linear_masks,
        "total_distance_to_linear_section": total_distance,
        "theorem_formula": expected,
        "maximum_synchronized_point_error": max(synchronized_errors.values()),
        "raw_future_samples_checked": len(future_checks),
        "largest_sampled_forward_gap": max(x for x, _ in future_checks),
        "largest_sampled_reverse_gap": max(y for _, y in future_checks),
    }


def verify_fibre_stripping_sharp(max_h: int = 4) -> dict[str, object]:
    result: dict[str, object] = {}
    for h in range(1, max_h + 1):
        d = 2 * h
        basis = tuple(1 << i for i in range(d))
        selected = tuple((1 << (d + i)) for i in range(h))
        dropped = tuple(
            selected[i] ^ (1 << (2 * i)) ^ (1 << (2 * i + 1))
            for i in range(h)
        )
        support = basis + selected + dropped
        transversal = basis + selected
        if not cycle_condition(support, basis):
            raise AssertionError(f"sharp fibre example failed contraction at h={h}")
        target = xor_sum(dropped)
        full_distances = word_lengths(support, d + h)
        transversal_distances = word_lengths(transversal, d + h)
        gap = transversal_distances[target] - full_distances[target]
        if full_distances[target] != h or transversal_distances[target] != 3 * h:
            raise AssertionError(f"sharp fibre distances failed at h={h}")
        if gap != 2 * h:
            raise AssertionError(f"sharp fibre gap failed at h={h}")
        result[str(h)] = {
            "ambient_dimension": d + h,
            "full_distance": full_distances[target],
            "transversal_distance": transversal_distances[target],
            "gap": gap,
            "bound": 2 * h,
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-exhaustive-w", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_geodesic_synchronization_results.json"
        ),
    )
    args = parser.parse_args()
    result = {
        "claim": (
            "Cycle contraction is equivalent to geodesicity; it yields "
            "uniform BLR synchronization and an all-raw-future quotient. "
            "The average constant three and the fibre-stripping factor two "
            "have explicit sharp families."
        ),
        "exhaustive_equivalence": exhaustive_equivalence(args.max_exhaustive_w),
        "sharp_bent_example": verify_sharp_example(),
        "sharp_fibre_stripping": verify_fibre_stripping_sharp(),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
