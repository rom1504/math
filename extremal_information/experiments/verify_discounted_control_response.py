#!/usr/bin/env python3
"""Exact rational checks for the discounted-control response benchmark."""

from __future__ import annotations

import itertools
import json
from fractions import Fraction
from typing import Iterable, List, Optional, Sequence, Tuple


Q = Fraction
Entry = Optional[Fraction]
Kernel = Tuple[Tuple[Entry, ...], ...]
Vector = Tuple[Fraction, ...]


def apply_block(kernel: Kernel, duration: int, value: Vector, lam: Fraction) -> Vector:
    """Apply U_(K,duration), using None for an impossible transition."""

    out = []
    discount = lam ** duration
    for row in kernel:
        candidates = [
            reward + discount * value[j]
            for j, reward in enumerate(row)
            if reward is not None
        ]
        assert candidates
        out.append(max(candidates))
    return tuple(out)


def compose(
    left: Kernel,
    left_duration: int,
    right: Kernel,
    lam: Fraction,
) -> Kernel:
    """Return the discounted max-plus kernel left odot right."""

    size = len(left)
    discount = lam ** left_duration
    rows: List[Tuple[Entry, ...]] = []
    for i in range(size):
        row: List[Entry] = []
        for k in range(size):
            candidates = [
                left[i][j] + discount * right[j][k]
                for j in range(size)
                if left[i][j] is not None and right[j][k] is not None
            ]
            row.append(max(candidates) if candidates else None)
        rows.append(tuple(row))
    return tuple(rows)


def finite_row_kernels(size: int, alphabet: Sequence[Entry]) -> Iterable[Kernel]:
    """Enumerate square kernels with at least one finite entry in every row."""

    for entries in itertools.product(alphabet, repeat=size * size):
        kernel = tuple(
            tuple(entries[size * i + j] for j in range(size))
            for i in range(size)
        )
        if all(any(entry is not None for entry in row) for row in kernel):
            yield kernel


def verify_composition() -> int:
    lam = Q(1, 2)
    kernels = tuple(finite_row_kernels(2, (None, Q(-1), Q(1))))
    values = tuple(
        tuple(Q(x) for x in raw)
        for raw in itertools.product((-2, 0, 3), repeat=2)
    )
    checks = 0
    for left in kernels:
        for right in kernels:
            combined = compose(left, 1, right, lam)
            for value in values:
                direct = apply_block(
                    combined, 3, value, lam
                )
                serial = apply_block(
                    left, 1, apply_block(right, 2, value, lam), lam
                )
                assert direct == serial
                checks += 1
    return checks


def hard_kernel(size: int, target: int) -> Kernel:
    """Every input state follows one zero-reward path to target."""

    return tuple(
        tuple(Q(0) if j == target else None for j in range(size))
        for _ in range(size)
    )


def oscillation(vector: Vector) -> Fraction:
    return max(vector) - min(vector)


def verify_response_metrics() -> Tuple[int, int]:
    values = tuple(
        tuple(Q(x) for x in raw)
        for raw in itertools.product((-2, -1, 1, 3), repeat=3)
    )
    absolute_checks = 0
    projective_checks = 0
    for lam in (Q(1, 2), Q(2, 3)):
        for depth in range(5):
            scale = lam ** depth
            for first_index, first in enumerate(values):
                for second in values[first_index + 1 :]:
                    difference = tuple(a - b for a, b in zip(first, second))
                    predicted_absolute = scale * max(abs(x) for x in difference)
                    exposed_absolute = max(
                        max(
                            abs(a - b)
                            for a, b in zip(
                                apply_block(hard_kernel(3, target), depth, first, lam),
                                apply_block(hard_kernel(3, target), depth, second, lam),
                            )
                        )
                        for target in range(3)
                    )
                    assert exposed_absolute == predicted_absolute
                    absolute_checks += 1

                    response_differences = tuple(scale * x for x in difference)
                    exposed_projective = oscillation(response_differences) / 2
                    predicted_projective = scale * oscillation(difference) / 2
                    assert exposed_projective == predicted_projective
                    projective_checks += 1
    return absolute_checks, projective_checks


def verify_bellman_contraction() -> int:
    """Exhaust max-choice kernels and check the sharp discounted upper bound."""

    kernels = tuple(finite_row_kernels(2, (None, Q(-1), Q(1))))
    values = tuple(
        tuple(Q(x) for x in raw)
        for raw in itertools.product((-2, 0, 3), repeat=2)
    )
    checks = 0
    for lam in (Q(1, 2), Q(2, 3)):
        for depth in (1, 2, 3):
            for kernel in kernels:
                for first_index, first in enumerate(values):
                    for second in values[first_index + 1 :]:
                        image_first = apply_block(kernel, depth, first, lam)
                        image_second = apply_block(kernel, depth, second, lam)
                        after = max(
                            abs(a - b) for a, b in zip(image_first, image_second)
                        )
                        before = max(abs(a - b) for a, b in zip(first, second))
                        assert after <= lam ** depth * before
                        checks += 1
    return checks


def geometric_sum(lam: Fraction, horizon: int) -> Fraction:
    return (1 - lam ** horizon) / (1 - lam)


def verify_value_cubes_and_fixed_points() -> Tuple[int, int]:
    horizon_checks = 0
    fixed_point_checks = 0
    reward_bound = Q(3)
    for lam in (Q(1, 2), Q(2, 3), Q(3, 4)):
        for horizon in range(1, 7):
            bound = reward_bound * geometric_sum(lam, horizon)
            candidates = (-bound, -bound / 3, Q(0), bound / 2, bound)
            for value in itertools.product(candidates, repeat=3):
                reward = tuple(
                    coordinate * (1 - lam) / (1 - lam ** horizon)
                    for coordinate in value
                )
                assert all(abs(x) <= reward_bound for x in reward)
                reconstructed = tuple(
                    x * geometric_sum(lam, horizon) for x in reward
                )
                assert reconstructed == value
                horizon_checks += 1

        rewards = tuple(
            tuple(Q(x) for x in raw)
            for raw in itertools.product((-3, -1, 0, 2, 3), repeat=3)
        )
        for reward in rewards:
            fixed = tuple(x / (1 - lam) for x in reward)
            assert tuple(x + lam * v for x, v in zip(reward, fixed)) == fixed
            fixed_point_checks += 1
        for first_index, first in enumerate(rewards):
            for second in rewards[first_index + 1 :]:
                reward_distance = max(abs(a - b) for a, b in zip(first, second))
                first_fixed = tuple(x / (1 - lam) for x in first)
                second_fixed = tuple(x / (1 - lam) for x in second)
                value_distance = max(
                    abs(a - b) for a, b in zip(first_fixed, second_fixed)
                )
                assert value_distance == reward_distance / (1 - lam)
                fixed_point_checks += 1
    return horizon_checks, fixed_point_checks


def verify_repeated_error_bound() -> int:
    checks = 0
    for lam in (Q(1, 3), Q(1, 2), Q(3, 4)):
        for delta in (Q(1, 7), Q(2, 5)):
            approximate = Q(0)
            for step in range(21):
                predicted = delta * (1 - lam ** step) / (1 - lam)
                assert approximate == predicted
                assert approximate <= delta / (1 - lam)
                checks += 1
                approximate = lam * approximate + delta
    return checks


def verify_finite_packing() -> Tuple[int, int]:
    """Check one absolute grid packing and one matching coordinate cover."""

    lam = Q(1, 2)
    depth = 2
    epsilon = Q(1, 4)
    bound = Q(4)
    input_radius = epsilon / (lam ** depth)

    packing_axis = tuple(
        -bound + 3 * input_radius * index
        for index in range(3)
    )
    packing = tuple(itertools.product(packing_axis, repeat=2))
    packing_checks = 0
    for index, first in enumerate(packing):
        for second in packing[index + 1 :]:
            distance = lam ** depth * max(
                abs(a - b) for a, b in zip(first, second)
            )
            assert distance > 2 * epsilon
            packing_checks += 1

    cover_axis = tuple(
        -bound + input_radius * index
        for index in range(1 + int(2 * bound / input_radius))
    )
    cover_checks = 0
    test_axis = tuple(-bound + Q(index, 4) for index in range(33))
    for value in itertools.product(test_axis, repeat=2):
        nearest = tuple(min(cover_axis, key=lambda point: abs(point - x)) for x in value)
        error = lam ** depth * max(abs(a - b) for a, b in zip(value, nearest))
        assert error <= epsilon
        cover_checks += 1
    return packing_checks, cover_checks


def main() -> None:
    absolute_checks, projective_checks = verify_response_metrics()
    horizon_checks, fixed_point_checks = verify_value_cubes_and_fixed_points()
    packing_checks, cover_checks = verify_finite_packing()
    print(
        json.dumps(
            {
                "bellman_contraction_checks": verify_bellman_contraction(),
                "composition_checks": verify_composition(),
                "depth_response_metric_checks": absolute_checks,
                "fixed_point_and_sharp_perturbation_checks": fixed_point_checks,
                "horizon_cube_realization_checks": horizon_checks,
                "projective_response_metric_checks": projective_checks,
                "response_cover_checks": cover_checks,
                "response_packing_pair_checks": packing_checks,
                "repeated_error_recurrence_checks": verify_repeated_error_bound(),
                "status": "all discounted-control response checks passed",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
