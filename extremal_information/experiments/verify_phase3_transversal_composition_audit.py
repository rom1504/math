#!/usr/bin/env python3
"""Independent finite audit for the dense-transversal composition note.

This script focuses on claims not isolated by the originating verifier:

* well-definedness, commutativity, associativity, and idempotence of the
  affine quotient product;
* the affine decoder after arbitrary appended futures;
* exact-support separation by the all-future radius query class; and
* the equal-depth obstruction, including the exceptional small order D=4.
"""

from __future__ import annotations

import json
from collections import deque
from itertools import product
from pathlib import Path


def wt(x: int) -> int:
    return bin(x).count("1")


def span(elements: tuple[int, ...] | list[int] | set[int]) -> frozenset[int]:
    result = {0}
    for element in elements:
        result |= {x ^ element for x in tuple(result)}
    return frozenset(result)


def subspaces(d: int) -> tuple[frozenset[int], ...]:
    spaces = {span(tuple(x for x in range(1 << d) if mask & (1 << x)))
              for mask in range(1 << (1 << d))}
    return tuple(sorted(spaces, key=lambda space: (len(space), tuple(space))))


def linear_value(images: tuple[int, ...], q: int) -> int:
    answer = 0
    for index, image in enumerate(images):
        if q & (1 << index):
            answer ^= image
    return answer


def image_space(images: tuple[int, ...]) -> frozenset[int]:
    return span(images)


def space_sum(*spaces: frozenset[int]) -> frozenset[int]:
    return span(tuple(x for space in spaces for x in space))


def coset_representative(x: int, space: frozenset[int]) -> int:
    return min(x ^ v for v in space)


State = tuple[tuple[int, ...], frozenset[int]]


def canonical_state(images: tuple[int, ...], space: frozenset[int]) -> State:
    return tuple(coset_representative(x, space) for x in images), space


def state_product(left: State, right: State) -> State:
    l_images, v_space = left
    k_images, z_space = right
    difference = tuple(a ^ b for a, b in zip(l_images, k_images))
    total_space = space_sum(v_space, z_space, image_space(difference))
    return canonical_state(l_images, total_space)


AffineState = tuple[int, frozenset[int]]


def affine_state_key(reference: int, direction: frozenset[int]) -> AffineState:
    return min(reference ^ vector for vector in direction), direction


def affine_join(left: AffineState, right: AffineState) -> AffineState:
    left_reference, left_direction = left
    right_reference, right_direction = right
    direction = span(tuple(left_direction | right_direction | {left_reference ^ right_reference}))
    return affine_state_key(left_reference, direction)


def verify_affine_product() -> dict[str, int]:
    d, k = 2, 2
    spaces = subspaces(d)
    raw_maps = tuple(product(range(1 << d), repeat=k))
    states = sorted(
        {canonical_state(tuple(images), space) for space in spaces for images in raw_maps},
        key=lambda state: (len(state[1]), state[0], tuple(state[1])),
    )

    representative_checks = 0
    for images, v_space in states:
        perturbations = tuple(product(tuple(v_space), repeat=k))
        for other_images, z_space in states:
            other_perturbations = tuple(product(tuple(z_space), repeat=k))
            expected = state_product((images, v_space), (other_images, z_space))
            for a in perturbations:
                for b in other_perturbations:
                    changed_left = tuple(x ^ y for x, y in zip(images, a))
                    changed_right = tuple(x ^ y for x, y in zip(other_images, b))
                    assert state_product(
                        canonical_state(changed_left, v_space),
                        canonical_state(changed_right, z_space),
                    ) == expected
                    representative_checks += 1

    law_checks = 0
    for a in states:
        assert state_product(a, a) == a
        for b in states:
            assert state_product(a, b) == state_product(b, a)
            for c in states:
                assert state_product(state_product(a, b), c) == state_product(
                    a, state_product(b, c)
                )
                law_checks += 1

    return {
        "canonical_states": len(states),
        "representative_checks": representative_checks,
        "associativity_triples": law_checks,
    }


def verify_affine_subspace_join() -> dict[str, int]:
    # This algebra is independent of the Cayley realization.  Dimension three
    # keeps all representative changes and associativity triples exhaustive.
    dimension = 3
    spaces = subspaces(dimension)
    states = sorted(
        {affine_state_key(reference, direction)
         for direction in spaces for reference in range(1 << dimension)},
        key=lambda state: (len(state[1]), state[0], tuple(state[1])),
    )
    representative_checks = 0
    associativity_checks = 0
    for left in states:
        for right in states:
            expected = affine_join(left, right)
            for left_shift in left[1]:
                for right_shift in right[1]:
                    assert affine_join(
                        affine_state_key(left[0] ^ left_shift, left[1]),
                        affine_state_key(right[0] ^ right_shift, right[1]),
                    ) == expected
                    representative_checks += 1
            assert affine_join(left, right) == affine_join(right, left)
            for third in states:
                assert affine_join(affine_join(left, right), third) == affine_join(
                    left, affine_join(right, third)
                )
                associativity_checks += 1
        assert affine_join(left, left) == left
    return {
        "affine_subspace_states": len(states),
        "representative_checks": representative_checks,
        "associativity_triples": associativity_checks,
    }


def word_profile(support: set[int] | tuple[int, ...], w: int) -> tuple[int, ...]:
    distances = [-1] * (1 << w)
    distances[0] = 0
    queue: deque[int] = deque([0])
    while queue:
        x = queue.popleft()
        for generator in support:
            y = x ^ generator
            if distances[y] < 0:
                distances[y] = distances[x] + 1
                queue.append(y)
    return tuple(distances)


def graph_support(images: tuple[int, ...], d: int, k: int) -> set[int]:
    support = {1 << i for i in range(d)}
    for q in range(1, 1 << k):
        support.add(linear_value(images, q) | (q << d))
    return support


def delta_profile(state: State, d: int, k: int) -> tuple[int, ...]:
    images, space = state
    answer = []
    for x in range(1 << (d + k)):
        kernel = x & ((1 << d) - 1)
        q = x >> d
        centre = linear_value(images, q)
        answer.append(min(wt(kernel ^ centre ^ v) for v in space) + (q != 0))
    return tuple(answer)


def minplus(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        min(left[y] + right[x ^ y] for y in range(len(left)) if right[x ^ y] >= 0)
        for x in range(len(left))
    )


def verify_future_decoder() -> dict[str, int]:
    d, k, w = 2, 1, 3
    maps = tuple(product(range(1 << d), repeat=k))
    futures = []
    for mask in range(1 << ((1 << w) - 1)):
        support = {x for x in range(1, 1 << w) if mask & (1 << (x - 1))}
        futures.append((support, word_profile(support, w)))

    checks = 0
    largest_gap = 0
    for family_mask in range(1, 1 << len(maps)):
        family = tuple(maps[index] for index in range(len(maps)) if family_mask & (1 << index))
        union: set[int] = set()
        for images in family:
            union |= graph_support(tuple(images), d, k)
        reference = tuple(family[0])
        difference_images = tuple(
            a ^ b for current in family[1:] for a, b in zip(current, reference)
        )
        direction_space = span(difference_images)
        state = canonical_state(reference, direction_space)
        delta = delta_profile(state, d, k)
        ell = len(family)
        affine_rank = len(direction_space).bit_length() - 1
        for future, future_profile in futures:
            actual = word_profile(union | future, w)
            decoded = minplus(delta, future_profile)
            for a, b in zip(actual, decoded):
                assert b <= a <= b + ell
                assert a <= b + affine_rank + 1
                largest_gap = max(largest_gap, a - b)
            checks += 1
    return {"family_future_checks": checks, "largest_decoder_gap": largest_gap}


def verify_query_separation() -> dict[str, int]:
    w = 4
    universe = set(range(1, 1 << w))
    checks = 0
    for mask in range(1 << len(universe)):
        support = {x for x in universe if mask & (1 << (x - 1))}
        for x in universe:
            future = universe - {x}
            radius = max(word_profile(support | future, w))
            assert radius == (1 if x in support else 2)
            checks += 1
    return {"support_membership_queries": checks}


def obstruction_radius(d: int, mixed: bool) -> int:
    t = (1 << d) - 1
    labels = [0]
    labels.extend((t ^ (1 << i)) if mixed else (1 << i) for i in range(d))
    support = {1 << i for i in range(d)}
    support |= {a | (1 << d) for a in labels}
    return max(word_profile(support, d + 1))


def verify_equal_depth_scope() -> dict[str, dict[str, int]]:
    answer: dict[str, dict[str, int]] = {}
    for d in (2, 4, 6, 8, 10, 12):
        sparse = obstruction_radius(d, mixed=False)
        mixed = obstruction_radius(d, mixed=True)
        if d >= 6:
            assert sparse == d and mixed == d // 2
        answer[str(d)] = {"sparse_radius": sparse, "mixed_radius": mixed}
    assert answer["4"]["mixed_radius"] == 3  # explains the theorem's D>=6 scope
    return answer


def main() -> None:
    results = {
        "affine_product": verify_affine_product(),
        "affine_subspace_join": verify_affine_subspace_join(),
        "all_future_affine_decoder": verify_future_decoder(),
        "all_future_query_separation": verify_query_separation(),
        "equal_depth_scope": verify_equal_depth_scope(),
        "status": "all assertions passed",
    }
    output = Path(__file__).with_name("phase3_transversal_composition_audit_results.json")
    output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
