#!/usr/bin/env python3
"""Exact checks for the query-local Walsh amalgamation theorem."""

from __future__ import annotations

from itertools import product
from math import sqrt


def dot(x: int, y: int) -> int:
    return bin(x & y).count("1") & 1


def span_value(labels: tuple[int, ...], mask: int) -> int:
    value = 0
    for i, label in enumerate(labels):
        if (mask >> i) & 1:
            value ^= label
    return value


def rooted_state(labels: tuple[int, ...], m: int) -> tuple[object, ...]:
    k = len(labels)
    gram = tuple(dot(x, y) for x in labels for y in labels)
    relations = tuple(
        mask for mask in range(1 << k) if span_value(labels, mask) == 0
    )
    omega = (1 << m) - 1
    roots = tuple(
        mask for mask in range(1 << k) if span_value(labels, mask) == omega
    )
    return gram, relations, roots


def cross_form(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(dot(x, y) for x in a for y in b)


def mixed_fibres(
    a: tuple[int, ...], b: tuple[int, ...], m: int
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    omega = (1 << m) - 1
    coincidences = []
    roots = []
    for ca in range(1 << len(a)):
        ua = span_value(a, ca)
        for cb in range(1 << len(b)):
            ub = span_value(b, cb)
            if ua == ub:
                coincidences.append((ca, cb))
            if (ua ^ ub) == omega:
                roots.append((ca, cb))
    return tuple(coincidences), tuple(roots)


def restrict_tuple(labels: tuple[int, ...], indices: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(labels[i] for i in indices)


def verify_local_reconstruction() -> int:
    """Restrictions reproduce the direct local state on random-looking cases."""
    checks = 0
    m = 5
    tuples_a = ((1, 2, 3), (3, 5, 6), (1, 1, 7))
    tuples_b = ((4, 6, 7), (2, 3, 5), (7, 7, 1))
    supports = (
        ((0,), (0, 1)),
        ((0, 2), (1,)),
        ((1, 2), (0, 2)),
        ((0, 1, 2), (0, 1, 2)),
    )
    for a in tuples_a:
        for b in tuples_b:
            for ia, ib in supports:
                local_a = restrict_tuple(a, ia)
                local_b = restrict_tuple(b, ib)
                direct = rooted_state(local_a + local_b, m)
                gram = (
                    rooted_state(local_a, m)[0],
                    rooted_state(local_b, m)[0],
                    cross_form(local_a, local_b),
                )
                # The three blocks determine exactly the direct Gram matrix.
                rebuilt_gram = []
                ka, kb = len(local_a), len(local_b)
                for i in range(ka + kb):
                    for j in range(ka + kb):
                        if i < ka and j < ka:
                            rebuilt_gram.append(gram[0][i * ka + j])
                        elif i >= ka and j >= ka:
                            rebuilt_gram.append(
                                gram[1][(i - ka) * kb + (j - ka)]
                            )
                        elif i < ka:
                            rebuilt_gram.append(gram[2][i * kb + (j - ka)])
                        else:
                            rebuilt_gram.append(gram[2][j * kb + (i - ka)])
                assert tuple(rebuilt_gram) == direct[0]
                coincidences, roots = mixed_fibres(local_a, local_b, m)
                direct_relations = tuple(
                    sorted(ca | (cb << ka) for ca, cb in coincidences)
                )
                direct_roots = tuple(
                    sorted(ca | (cb << ka) for ca, cb in roots)
                )
                assert direct_relations == direct[1]
                assert direct_roots == direct[2]
                checks += 1
    return checks


def verify_semantic_packing(max_h: int = 6) -> int:
    checks = 0
    for h in range(1, max_h + 1):
        unused = 1 if h % 2 == 0 else 2
        m = 3 * h + unused
        assert m % 2 == 1
        a_labels = []
        c_labels = []
        for i in range(h):
            base = 3 * i
            a_labels.append((1 << base) | (1 << (base + 1)) | (1 << (base + 2)))
            c_labels.append(1 << (base + 2))

        piece_a = tuple(x for a in a_labels for x in (a, a))
        reference_b_state = None
        reference_cross = None
        observed_j = set()
        for sigma in product((0, 1), repeat=h):
            piece_b = tuple(
                c_labels[i] if sigma[i] else a_labels[i] for i in range(h)
            )
            state_b = rooted_state(piece_b, m)
            cross = cross_form(piece_a, piece_b)
            coincidences, roots = mixed_fibres(piece_a, piece_b, m)
            if reference_b_state is None:
                reference_b_state = state_b
                reference_cross = cross
            assert state_b == reference_b_state
            assert cross == reference_cross
            assert roots == ()

            # Record, for each middle label, whether its singleton span meets
            # the endpoint span nontrivially.  This is precisely sigma_i=0.
            visible_bits = []
            endpoint_span = {
                span_value(piece_a, mask) for mask in range(1 << len(piece_a))
            }
            for label in piece_b:
                visible_bits.append(int(label in endpoint_span))
            assert tuple(visible_bits) == tuple(1 - bit for bit in sigma)
            observed_j.add(tuple(coincidences))
            checks += 1
        assert len(observed_j) == 1 << h
        checks += 1

    delta = (7 - 3 * sqrt(3)) / 2
    assert delta > 0
    assert 3.5 - 1.5 * sqrt(3) == delta
    checks += 2
    return checks


def main() -> None:
    checks = verify_local_reconstruction() + verify_semantic_packing()
    print(f"query-local Walsh amalgamation checks passed: {checks}")


if __name__ == "__main__":
    main()
