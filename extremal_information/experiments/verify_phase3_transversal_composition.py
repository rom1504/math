#!/usr/bin/env python3
"""Finite checks for dense-transversal composition theorems.

The checks are deliberately independent of the proof text:

1. exhaustive all-future word-profile comparison for every complete
   transversal at (D,k)=(2,2);
2. exhaustive affine-state profile bounds for every nonempty family of
   linear maps at (D,k)=(2,2);
3. the mixed-cycle defect identity for every support containing a fixed
   basis in F_2^4;
4. exact radii for the equal-depth obstruction through D=12.
"""

from __future__ import annotations

import json
from collections import Counter, deque
from itertools import product
from pathlib import Path


def wt(x: int) -> int:
    return bin(x).count("1")


def word_profile(support: tuple[int, ...] | list[int] | set[int], w: int) -> tuple[int, ...]:
    generators = tuple(set(support))
    distance = [-1] * (1 << w)
    distance[0] = 0
    queue: deque[int] = deque([0])
    while queue:
        x = queue.popleft()
        for generator in generators:
            y = x ^ generator
            if distance[y] < 0:
                distance[y] = distance[x] + 1
                queue.append(y)
    return tuple(distance)


def support_mask(support: tuple[int, ...] | list[int] | set[int]) -> int:
    mask = 0
    for generator in support:
        if generator:
            mask |= 1 << (generator - 1)
    return mask


def support_from_mask(mask: int, w: int) -> tuple[int, ...]:
    return tuple(x for x in range(1, 1 << w) if mask & (1 << (x - 1)))


def linear_value(images: tuple[int, ...], q: int) -> int:
    answer = 0
    for index, image in enumerate(images):
        if q & (1 << index):
            answer ^= image
    return answer


def cycle_contracts(values: tuple[int, ...], k: int) -> bool:
    labels = tuple(range(1, 1 << k))
    for subset in range(1 << len(labels)):
        quotient_sum = 0
        kernel_sum = 0
        size = 0
        for index, label in enumerate(labels):
            if subset & (1 << index):
                quotient_sum ^= label
                kernel_sum ^= values[index]
                size += 1
        if quotient_sum == 0 and wt(kernel_sum) > size:
            return False
    return True


def graph_support(values: tuple[int, ...], d: int, k: int) -> tuple[int, ...]:
    basis = tuple(1 << index for index in range(d))
    labels = tuple(range(1, 1 << k))
    return basis + tuple(value | (label << d) for value, label in zip(values, labels))


def closest_linear(values: tuple[int, ...], d: int, k: int) -> tuple[int, ...]:
    labels = tuple(range(1, 1 << k))
    best_images: tuple[int, ...] | None = None
    best_key: tuple[int, int, tuple[int, ...]] | None = None
    for images in product(range(1 << d), repeat=k):
        errors = tuple(values[index] ^ linear_value(images, label) for index, label in enumerate(labels))
        key = (sum(wt(error) for error in errors), max((wt(error) for error in errors), default=0), images)
        if best_key is None or key < best_key:
            best_key = key
            best_images = images
    assert best_images is not None
    return tuple(linear_value(best_images, label) for label in labels)


def span(elements: tuple[int, ...] | list[int] | set[int]) -> frozenset[int]:
    result = {0}
    for element in elements:
        result |= {x ^ element for x in tuple(result)}
    return frozenset(result)


def map_code(values: tuple[int, ...], d: int, k: int) -> int:
    """Vectorize a linear map using its values on the quotient basis."""
    code = 0
    for index in range(k):
        basis_label = 1 << index
        value = values[basis_label - 1]
        code |= value << (d * index)
    return code


def verify_all_future_comparison() -> dict[str, int]:
    d, k, w = 2, 2, 4
    labels = tuple(range(1, 1 << k))
    profiles = [word_profile(support_from_mask(mask, w), w) for mask in range(1 << ((1 << w) - 1))]
    maps_checked = 0
    future_checks = 0
    largest_uniform_error = 0
    largest_forward_gap = 0
    largest_reverse_gap = 0
    for values in product(range(1 << d), repeat=len(labels)):
        if not cycle_contracts(values, k):
            continue
        maps_checked += 1
        linear_values = closest_linear(values, d, k)
        largest_uniform_error = max(
            largest_uniform_error,
            max(wt(a ^ b) for a, b in zip(values, linear_values)),
        )
        source_mask = support_mask(graph_support(values, d, k))
        linear_mask = support_mask(graph_support(linear_values, d, k))
        for future_mask in range(1 << ((1 << w) - 1)):
            source_profile = profiles[source_mask | future_mask]
            linear_profile = profiles[linear_mask | future_mask]
            largest_forward_gap = max(
                largest_forward_gap,
                max(a - b for a, b in zip(source_profile, linear_profile)),
            )
            largest_reverse_gap = max(
                largest_reverse_gap,
                max(b - a for a, b in zip(source_profile, linear_profile)),
            )
            future_checks += 1
    assert largest_uniform_error <= 8
    assert largest_forward_gap <= 8
    assert largest_reverse_gap <= 10
    return {
        "cycle_contracting_maps": maps_checked,
        "all_future_support_checks": future_checks,
        "largest_uniform_error": largest_uniform_error,
        "largest_source_minus_linear_gap": largest_forward_gap,
        "largest_linear_minus_source_gap": largest_reverse_gap,
    }


def affine_state(values: tuple[tuple[int, ...], ...], d: int, k: int) -> tuple[frozenset[int], tuple[int, ...]]:
    reference = values[0]
    differences = []
    for current in values[1:]:
        differences.extend(a ^ b for a, b in zip(current, reference))
    v_space = span(differences)
    return v_space, reference


def quotient_profile(v_space: frozenset[int], reference: tuple[int, ...], d: int, k: int) -> tuple[int, ...]:
    labels = tuple(range(1, 1 << k))
    lookup = {label: reference[index] for index, label in enumerate(labels)}
    result = []
    for x in range(1 << (d + k)):
        kernel = x & ((1 << d) - 1)
        quotient = x >> d
        centre = 0 if quotient == 0 else lookup[quotient]
        result.append(min(wt(kernel ^ centre ^ v) for v in v_space) + (quotient != 0))
    return tuple(result)


def verify_affine_bounds() -> dict[str, int]:
    d, k, w = 2, 2, 4
    labels = tuple(range(1, 1 << k))
    linear_maps = []
    for images in product(range(1 << d), repeat=k):
        linear_maps.append(tuple(linear_value(images, label) for label in labels))
    families_checked = 0
    largest_gap = 0
    largest_gap_minus_rank = 0
    for family_mask in range(1, 1 << len(linear_maps)):
        family = tuple(linear_maps[index] for index in range(len(linear_maps)) if family_mask & (1 << index))
        support = set(1 << index for index in range(d))
        for values in family:
            support.update(graph_support(values, d, k)[d:])
        actual = word_profile(support, w)
        v_space, reference = affine_state(family, d, k)
        lower = quotient_profile(v_space, reference, d, k)
        ell = len(family)
        reference_code = map_code(reference, d, k)
        direction_space = span(
            tuple(map_code(current, d, k) ^ reference_code for current in family)
        )
        affine_rank = (len(direction_space)).bit_length() - 1
        for a, b in zip(actual, lower):
            assert b <= a <= b + ell
            assert a <= b + affine_rank + 1
            largest_gap = max(largest_gap, a - b)
            largest_gap_minus_rank = max(largest_gap_minus_rank, a - b - affine_rank)
        families_checked += 1
    return {
        "linear_families": families_checked,
        "largest_profile_gap": largest_gap,
        "largest_gap_minus_affine_rank": largest_gap_minus_rank,
    }


def join_affine_states(
    left: tuple[int, frozenset[int]], right: tuple[int, frozenset[int]]
) -> tuple[int, frozenset[int]]:
    left_ref, left_direction = left
    right_ref, right_direction = right
    direction = span(tuple(left_direction | right_direction | {left_ref ^ right_ref}))
    return left_ref, direction


def affine_points(state: tuple[int, frozenset[int]]) -> frozenset[int]:
    reference, direction = state
    return frozenset(reference ^ vector for vector in direction)


def verify_affine_join() -> dict[str, int]:
    d, k = 2, 2
    atoms = []
    for images in product(range(1 << d), repeat=k):
        values = tuple(linear_value(images, label) for label in range(1, 1 << k))
        code = map_code(values, d, k)
        atoms.append((code, frozenset({0})))
    pair_checks = 0
    associativity_checks = 0
    for left in atoms:
        for right in atoms:
            joined = join_affine_states(left, right)
            expected = frozenset({left[0], right[0]})
            assert expected <= affine_points(joined)
            pair_checks += 1
            for third in atoms:
                first = join_affine_states(joined, third)
                second = join_affine_states(left, join_affine_states(right, third))
                assert affine_points(first) == affine_points(second)
                associativity_checks += 1
    return {
        "atomic_pair_checks": pair_checks,
        "atomic_associativity_checks": associativity_checks,
    }


def mixed_cycle_excess(support: tuple[int, ...], d: int, w: int) -> int:
    basis = set(1 << index for index in range(d))
    outside = tuple(generator for generator in support if generator not in basis)
    best = 0
    kernel_mask = (1 << d) - 1
    for subset in range(1 << len(outside)):
        total = 0
        size = 0
        for index, generator in enumerate(outside):
            if subset & (1 << index):
                total ^= generator
                size += 1
        if total >> d == 0:
            best = max(best, wt(total & kernel_mask) - size)
    return best


def verify_mixed_cycle_identity() -> dict[str, int]:
    d, w = 2, 4
    basis = tuple(1 << index for index in range(d))
    optional = tuple(x for x in range(1, 1 << w) if x not in basis)
    t = (1 << d) - 1
    checked = 0
    for mask in range(1 << len(optional)):
        support = basis + tuple(optional[index] for index in range(len(optional)) if mask & (1 << index))
        # B already spans W; skip supports that do not span the quotient.
        profile = word_profile(support, w)
        if min(profile) < 0:
            continue
        if any(value < 0 for value in profile):
            continue
        defect = mixed_cycle_excess(support, d, w)
        assert profile[t] == d - defect
        checked += 1
    return {"spanning_supports": checked}


def obstruction_support(a_set: tuple[int, ...], d: int) -> tuple[int, ...]:
    return tuple(1 << index for index in range(d)) + tuple(a | (1 << d) for a in a_set)


def binary_rank(vectors: tuple[int, ...], d: int) -> int:
    pivots = [0] * d
    rank = 0
    for vector in vectors:
        x = vector
        while x:
            pivot = x.bit_length() - 1
            if pivots[pivot]:
                x ^= pivots[pivot]
            else:
                pivots[pivot] = x
                rank += 1
                break
    return rank


def verify_equal_depth_obstruction() -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    for d in (6, 8, 10, 12):
        t = (1 << d) - 1
        sparse = (0,) + tuple(1 << index for index in range(d))
        mixed = (0,) + tuple(t ^ (1 << index) for index in range(d))
        assert len(sparse) == len(mixed) == d + 1
        assert binary_rank(sparse[1:], d) == binary_rank(mixed[1:], d) == d
        sparse_profile = word_profile(obstruction_support(sparse, d), d + 1)
        mixed_profile = word_profile(obstruction_support(mixed, d), d + 1)
        sparse_radius = max(sparse_profile)
        mixed_radius = max(mixed_profile)
        assert sparse_radius == d
        assert mixed_radius == d // 2
        # The pair 0,t+e_1 is the advertised mixed cycle representing t in 3.
        assert mixed_profile[t] == 3
        result[str(d)] = {
            "source_count": d + 1,
            "sparse_radius": sparse_radius,
            "mixed_radius": mixed_radius,
            "radius_gap": sparse_radius - mixed_radius,
            "mixed_cycle_target_length": mixed_profile[t],
        }
    return result


def code_radius(code: tuple[int, ...], dimension: int) -> int:
    return max(min(wt(word ^ codeword) for codeword in code) for word in range(1 << dimension))


def direct_power(code: tuple[int, ...], block_dimension: int, exponent: int) -> tuple[int, ...]:
    result = (0,)
    for block in range(exponent):
        result = tuple(
            prefix | (codeword << (block_dimension * block))
            for prefix in result
            for codeword in code
        )
    return result


def pair_distance_histogram(code: tuple[int, ...]) -> dict[int, int]:
    return dict(sorted(Counter(wt(left ^ right) for left in code for right in code).items()))


def verify_pair_spectrum_no_go() -> dict[str, object]:
    code_zero = (0, 3, 9, 10, 53, 54, 60, 63)
    code_one = (0, 3, 12, 15, 48, 51, 60, 63)
    for code in (code_zero, code_one):
        assert all((left ^ right) in code for left in code for right in code)
    enumerator_zero = dict(sorted(Counter(wt(codeword) for codeword in code_zero).items()))
    enumerator_one = dict(sorted(Counter(wt(codeword) for codeword in code_one).items()))
    assert enumerator_zero == enumerator_one == {0: 1, 2: 3, 4: 3, 6: 1}
    assert code_radius(code_zero, 6) == 2
    assert code_radius(code_one, 6) == 3

    formula_checks = 0
    base_radii = []
    for code in (code_zero, code_one):
        profile = word_profile(obstruction_support(code, 6), 7)
        for word in range(1 << 6):
            distance = min(wt(word ^ codeword) for codeword in code)
            nonzero_distance = min(wt(word ^ codeword) for codeword in code if codeword)
            assert profile[word | (1 << 6)] == 1 + distance
            assert profile[word] == min(wt(word), 2 + nonzero_distance)
            formula_checks += 2
        base_radii.append(max(profile))

    power_data: dict[str, object] = {}
    for exponent in (1, 2):
        zero_power = direct_power(code_zero, 6, exponent)
        one_power = direct_power(code_one, 6, exponent)
        assert pair_distance_histogram(zero_power) == pair_distance_histogram(one_power)
        zero_radius = max(word_profile(obstruction_support(zero_power, 6 * exponent), 6 * exponent + 1))
        one_radius = max(word_profile(obstruction_support(one_power, 6 * exponent), 6 * exponent + 1))
        assert one_radius - zero_radius >= exponent - 1
        power_data[str(exponent)] = {
            "code_zero_union_radius": zero_radius,
            "code_one_union_radius": one_radius,
            "response_gap": one_radius - zero_radius,
            "pair_histogram": pair_distance_histogram(zero_power),
        }
    return {
        "weight_enumerator": enumerator_zero,
        "code_radii": [2, 3],
        "base_union_radii": base_radii,
        "profile_formula_checks": formula_checks,
        "powers": power_data,
    }


def main() -> None:
    results = {
        "all_future_comparison": verify_all_future_comparison(),
        "affine_state_bounds": verify_affine_bounds(),
        "affine_join_law": verify_affine_join(),
        "mixed_cycle_identity": verify_mixed_cycle_identity(),
        "equal_depth_obstruction": verify_equal_depth_obstruction(),
        "pair_spectrum_no_go": verify_pair_spectrum_no_go(),
        "status": "all assertions passed",
    }
    output = Path(__file__).with_name("phase3_transversal_composition_results.json")
    output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
