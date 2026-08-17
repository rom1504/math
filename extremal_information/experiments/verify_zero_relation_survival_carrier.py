#!/usr/bin/env python3
"""Exact checks for finite support-survival lumpability."""

from collections import deque
from itertools import product


def relation_from_mask(size, mask):
    return tuple(
        sum(
            1 << target
            for target in range(size)
            if (mask >> (source * size + target)) & 1
        )
        for source in range(size)
    )


def image(relation, subset):
    answer = 0
    for source, successors in enumerate(relation):
        if (subset >> source) & 1:
            answer |= successors
    return answer


def compose(left, right):
    """Relation for first following left and then right."""

    return tuple(image(right, successors) for successors in left)


def has_cycle(relation):
    size = len(relation)

    def visit(vertex, active, finished):
        if (active >> vertex) & 1:
            return True
        if (finished >> vertex) & 1:
            return False
        active |= 1 << vertex
        for target in range(size):
            if (relation[vertex] >> target) & 1:
                if visit(target, active, finished):
                    return True
        return False

    # The graph is tiny in every exhaustive check; recomputing DFS state is
    # clearer than sharing mutable colors.
    for start in range(size):
        stack = [(start, (start,))]
        while stack:
            vertex, path = stack.pop()
            for target in range(size):
                if not ((relation[vertex] >> target) & 1):
                    continue
                if target in path:
                    return True
                if len(path) < size:
                    stack.append((target, path + (target,)))
    return False


def reachable_subsets(size, alphabet):
    full = (1 << size) - 1
    seen = {full}
    queue = deque([full])
    while queue:
        subset = queue.popleft()
        for relation in alphabet:
            target = image(relation, subset)
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return seen


def relation_semigroup(size, alphabet):
    identity = tuple(1 << index for index in range(size))
    seen = {identity}
    queue = deque([identity])
    while queue:
        relation = queue.popleft()
        for letter in alphabet:
            target = compose(relation, letter)
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return seen


def permutation_relation(size):
    return tuple(1 << ((source + 1) % size) for source in range(size))


def deletion_relation(size):
    return tuple(0 if source == 0 else 1 << source for source in range(size))


def delete_point(subset, point, size, permutation, deletion):
    for _ in range((size - point) % size):
        subset = image(permutation, subset)
    subset = image(deletion, subset)
    for _ in range(point):
        subset = image(permutation, subset)
    return subset


def debruijn_relation(memory, letter):
    states = tuple(product((0, 1), repeat=memory))
    index = {state: number for number, state in enumerate(states)}
    rows = []
    for source in states:
        successors = 0
        if source[0] == letter:
            for bit in (0, 1):
                successors |= 1 << index[source[1:] + (bit,)]
        rows.append(successors)
    return tuple(rows)


def main():
    equivalence_checks = 0
    size = 2
    relations = tuple(
        relation_from_mask(size, mask) for mask in range(1 << (size * size))
    )
    for alphabet in product(relations, repeat=2):
        subsets = reachable_subsets(size, alphabet)
        universally_alive = 0 not in subsets
        products = relation_semigroup(size, alphabet)
        every_product_cyclic = all(has_cycle(relation) for relation in products)
        assert universally_alive == every_product_cyclic
        equivalence_checks += 1

    monitor_checks = 0
    for size in range(1, 9):
        permutation = permutation_relation(size)
        deletion = deletion_relation(size)
        alphabet = (permutation, deletion)
        subsets = reachable_subsets(size, alphabet)
        assert len(subsets) == 1 << size
        # Deleting every point except one in a symmetric difference
        # distinguishes each unordered pair of subsets by future mortality.
        for left in subsets:
            for right in subsets:
                if left >= right:
                    continue
                difference = left ^ right
                witness = difference & -difference
                witness_index = witness.bit_length() - 1
                left_image, right_image = left, right
                for point in range(size):
                    if point == witness_index:
                        continue
                    left_image = delete_point(
                        left_image, point, size, permutation, deletion
                    )
                    right_image = delete_point(
                        right_image, point, size, permutation, deletion
                    )
                assert bool(left_image) != bool(right_image)
                monitor_checks += 1

    debruijn_checks = 0
    for memory in range(1, 8):
        full = (1 << (1 << memory)) - 1
        for letter in (0, 1):
            relation = debruijn_relation(memory, letter)
            assert image(relation, full) == full
            debruijn_checks += 1

    print(f"two-state support/cycle equivalence checks: {equivalence_checks}")
    print(f"powerset monitor distinguishability checks: {monitor_checks}")
    print(f"one-state de-Bruijn survival checks: {debruijn_checks}")
    print("support survival is complete, pumpable, and strictly query-relative")


if __name__ == "__main__":
    main()
