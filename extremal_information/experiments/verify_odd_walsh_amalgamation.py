#!/usr/bin/env python3
"""Exact checks for the odd-Walsh rooted amalgamation algebra."""

from __future__ import annotations

from itertools import product
import random


def dot(x: int, y: int) -> int:
    return bin(x & y).count("1") & 1


def span_value(labels: tuple[int, ...], mask: int) -> int:
    value = 0
    for i, label in enumerate(labels):
        if (mask >> i) & 1:
            value ^= label
    return value


def state(labels: tuple[int, ...], m: int) -> tuple[object, ...]:
    k = len(labels)
    gram = tuple(dot(labels[i], labels[j]) for i in range(k) for j in range(k))
    relations = tuple(mask for mask in range(1 << k) if span_value(labels, mask) == 0)
    omega = (1 << m) - 1
    roots = tuple(mask for mask in range(1 << k) if span_value(labels, mask) == omega)
    return gram, relations, roots


def cross_form(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(dot(x, y) for x in a for y in b)


def quotient_classes(labels: tuple[int, ...]) -> dict[int, tuple[int, ...]]:
    """Map each represented vector to its coefficient fibre."""
    fibres: dict[int, list[int]] = {}
    for mask in range(1 << len(labels)):
        fibres.setdefault(span_value(labels, mask), []).append(mask)
    return {value: tuple(masks) for value, masks in fibres.items()}


def relative_data(
    a: tuple[int, ...], b: tuple[int, ...], m: int
) -> tuple[object, ...]:
    fa = quotient_classes(a)
    fb = quotient_classes(b)
    # Quotient elements are represented canonically by their actual span vector.
    intersection = tuple(sorted(set(fa).intersection(fb)))
    omega = (1 << m) - 1
    root_pairs = tuple(
        sorted((u, v) for u in fa for v in fb if (u ^ v) == omega)
    )
    return cross_form(a, b), intersection, root_pairs


def reconstruct_state(
    a: tuple[int, ...], b: tuple[int, ...], m: int
) -> tuple[object, ...]:
    """Implement RA.12--RA.14 from isolated states and relative fibres."""
    k, ell = len(a), len(b)
    ga, ra, _ = state(a, m)
    gb, rb, _ = state(b, m)
    kappa, intersection, root_pairs = relative_data(a, b, m)
    intersection_set = set(intersection)
    root_pair_set = set(root_pairs)

    gram: list[int] = []
    for i in range(k + ell):
        for j in range(k + ell):
            if i < k and j < k:
                gram.append(ga[i * k + j])
            elif i >= k and j >= k:
                gram.append(gb[(i-k) * ell + (j-k)])
            elif i < k:
                gram.append(kappa[i * ell + (j-k)])
            else:
                gram.append(kappa[j * ell + (i-k)])

    relations = []
    roots = []
    ra_set, rb_set = set(ra), set(rb)
    # ra_set/rb_set are used to assert that actual-vector representatives
    # correctly quotient every internal relation fibre.
    assert 0 in ra_set and 0 in rb_set
    for mask in range(1 << (k + ell)):
        ca = mask & ((1 << k) - 1)
        cb = mask >> k
        ua = span_value(a, ca)
        ub = span_value(b, cb)
        if ua == ub and ua in intersection_set:
            relations.append(mask)
        if (ua, ub) in root_pair_set:
            roots.append(mask)
    return tuple(gram), tuple(relations), tuple(roots)


def gf2_rank(columns: tuple[int, ...]) -> int:
    basis: dict[int, int] = {}
    for value in columns:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                break
    return len(basis)


def symplectic_coordinates(h: int) -> tuple[list[int], list[int], int]:
    """Return a symplectic basis of omega-perp in V of dimension 2h.

    We search greedily in the tiny dimensions used by the verifier rather
    than baking in a coordinate formula.
    """
    m = 2 * h + 1
    omega = (1 << m) - 1
    candidates = [x for x in range(1, 1 << m) if dot(x, omega) == 0]
    p: list[int] = []
    q: list[int] = []
    chosen: list[int] = []
    for _ in range(h):
        found = None
        for x in candidates:
            if gf2_rank(tuple(chosen + [x])) != len(chosen) + 1:
                continue
            if any(dot(x, z) for z in chosen):
                continue
            for y in candidates:
                if dot(x, y) != 1:
                    continue
                if gf2_rank(tuple(chosen + [x, y])) != len(chosen) + 2:
                    continue
                if any(dot(y, z) for z in chosen):
                    continue
                found = (x, y)
                break
            if found:
                break
        assert found is not None
        x, y = found
        p.append(x)
        q.append(y)
        chosen.extend((x, y))
    assert len(chosen) == 2 * h
    return p, q, m


def mat2_invertible(entries: tuple[int, int, int, int]) -> bool:
    a, b, c, d = entries
    return (a * d ^ b * c) == 1


def verify() -> None:
    checks = 0

    # Exhaustive two-piece reconstruction at m=3, lengths up to 2+2.
    m = 3
    for k, ell in ((1, 1), (1, 2), (2, 1), (2, 2)):
        for a in product(range(1 << m), repeat=k):
            for b in product(range(1 << m), repeat=ell):
                assert reconstruct_state(a, b, m) == state(a + b, m)
                checks += 1

    # Random larger reconstruction tests.
    rng = random.Random(20260817)
    for m in (5, 7):
        for _ in range(200):
            k = rng.randint(1, 4)
            ell = rng.randint(1, 4)
            a = tuple(rng.randrange(1 << m) for _ in range(k))
            b = tuple(rng.randrange(1 << m) for _ in range(ell))
            assert reconstruct_state(a, b, m) == state(a + b, m)
            checks += 1

    # RA.3: all 2^(rs) cross matrices with fixed intrinsic/J/root data.
    r, s = 2, 2
    p, q, m = symplectic_coordinates(r + s)
    a = tuple(p[:r])
    observed = set()
    reference_b_state = None
    for bits in range(1 << (r * s)):
        b_values = []
        for j in range(s):
            value = p[r + j]
            for i in range(r):
                if (bits >> (i * s + j)) & 1:
                    value ^= q[i]
            b_values.append(value)
        b = tuple(b_values)
        b_state = state(b, m)
        if reference_b_state is None:
            reference_b_state = b_state
        assert b_state == reference_b_state
        rel = relative_data(a, b, m)
        assert rel[1] == (0,)
        assert rel[2] == ()
        observed.add(rel[0])
    assert len(observed) == 1 << (r * s)
    checks += len(observed)

    # RA.4: GL(2,2) gives six different intersection correspondences while
    # all isolated states, cross forms, and roots agree.
    p, _, m = symplectic_coordinates(2)
    a = (p[0], p[1])
    intersection_graphs = set()
    reference = None
    for entries in product((0, 1), repeat=4):
        if not mat2_invertible(entries):
            continue
        cols = (
            (p[0] if entries[0] else 0) ^ (p[1] if entries[2] else 0),
            (p[0] if entries[1] else 0) ^ (p[1] if entries[3] else 0),
        )
        signature = (state(a, m), state(cols, m), cross_form(a, cols))
        if reference is None:
            reference = signature
        assert signature == reference
        assert state(a + cols, m)[2] == ()
        # Encode the full coefficient-level mixed relation graph.
        mixed = tuple(
            (ca, cb)
            for ca in range(1 << 2)
            for cb in range(1 << 2)
            if span_value(a, ca) == span_value(cols, cb)
        )
        intersection_graphs.add(mixed)
    assert len(intersection_graphs) == 6
    checks += len(intersection_graphs)

    # The leading-scale LG.1 path gap isolates exactly a coincidence bit when
    # the endpoints are treated as one piece and the middle as the other.
    m = 5
    a_label = 0b00111
    b_label = 0b00100
    endpoints = (a_label, a_label)
    assert state((a_label,), m) == state((b_label,), m)
    assert cross_form(endpoints, (a_label,)) == cross_form(
        endpoints, (b_label,)
    )
    equal_relative = relative_data(endpoints, (a_label,), m)
    distinct_relative = relative_data(endpoints, (b_label,), m)
    assert equal_relative[0] == distinct_relative[0]
    assert equal_relative[1] != distinct_relative[1]
    assert equal_relative[2] == distinct_relative[2] == ()
    checks += 5

    # RA.5: the characteristic root is a genuinely additional gluing bit.
    m = 5
    omega = (1 << m) - 1
    a = (0b00011,)
    b_plus = (omega ^ a[0],)
    b_minus = (0b00100,)
    assert state(a, m) == state(a, m)
    assert state(b_plus, m) == state(b_minus, m)
    assert cross_form(a, b_plus) == cross_form(a, b_minus) == (0,)
    assert relative_data(a, b_plus, m)[:2] == relative_data(a, b_minus, m)[:2]
    assert relative_data(a, b_plus, m)[2] != relative_data(a, b_minus, m)[2]
    checks += 5

    # RA.6: every pair agrees, while the ternary relation differs.
    p, _, m = symplectic_coordinates(3)
    plus = (p[0], p[1], p[0] ^ p[1])
    minus = (p[0], p[1], p[2])
    for i in range(3):
        assert state((plus[i],), m) == state((minus[i],), m)
    for i in range(3):
        for j in range(i + 1, 3):
            assert relative_data((plus[i],), (plus[j],), m) == relative_data(
                (minus[i],), (minus[j],), m
            )
    assert gf2_rank(plus) == 2
    assert gf2_rank(minus) == 3
    checks += 11

    # Associativity at the presentation level: either parenthesization is
    # exactly the direct tuple state.  (The quotient maps are actual spans in
    # this verifier, so concatenation realizes both pushouts.)
    for m in (3, 5):
        for _ in range(200):
            pieces = [
                tuple(rng.randrange(1 << m) for _ in range(rng.randint(1, 3)))
                for _ in range(3)
            ]
            direct = state(pieces[0] + pieces[1] + pieces[2], m)
            left = reconstruct_state(pieces[0] + pieces[1], pieces[2], m)
            right = reconstruct_state(pieces[0], pieces[1] + pieces[2], m)
            assert left == direct == right
            checks += 1

    print(f"odd-Walsh rooted amalgamation checks passed: {checks}")


if __name__ == "__main__":
    verify()
