#!/usr/bin/env python3
"""Exact checks for drafts/boolean_port_product_algebra_closure.md."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import floor
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_bcx_two_port_holonomy import regular_hadamard  # noqa: E402


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def selector(a: tuple[int, ...] | np.ndarray) -> int:
    total = int(sum(a))
    if total > 0:
        return 1
    if total < 0:
        return -1
    return int(a[0])


def subsets(p: int):
    for size in range(p + 1):
        yield from combinations(range(p), size)


def fourier_numerator(p: int, subset: tuple[int, ...]) -> int:
    total = 0
    for a in product((-1, 1), repeat=p):
        character = 1
        for i in subset:
            character *= a[i]
        total += selector(a) * character
    return total


def active_sets(p: int) -> list[tuple[int, ...]]:
    return [s for s in subsets(p) if fourier_numerator(p, s) != 0]


def port_product(ports: np.ndarray, subset: tuple[int, ...]) -> np.ndarray:
    if not subset:
        return np.ones(ports.shape[1], dtype=np.int64)
    return np.prod(ports[list(subset)], axis=0)


def selector_witness(ports: np.ndarray, eps: tuple[int, ...]) -> np.ndarray:
    signed = np.asarray(eps, dtype=np.int64)[:, None] * ports
    return np.asarray(
        [selector(tuple(signed[:, j])) for j in range(ports.shape[1])],
        dtype=np.int64,
    )


def assert_closed(h: np.ndarray, r: int, ports: np.ndarray,
                  active: list[tuple[int, ...]]) -> None:
    for w in ports:
        assert np.array_equal(h @ w, r * w)
    for s in active:
        w = port_product(ports, s)
        assert np.array_equal(h @ w, r * w), s


def projective_histogram(ports: np.ndarray) -> dict[tuple[int, ...], int]:
    histogram: dict[tuple[int, ...], int] = {}
    for row in ports.T:
        gauged = tuple((int(row[0]) * row).tolist())
        histogram[gauged] = histogram.get(gauged, 0) + 1
    return histogram


def histogram_convolution(
    left: dict[tuple[int, ...], int],
    right: dict[tuple[int, ...], int],
) -> dict[tuple[int, ...], int]:
    result: dict[tuple[int, ...], int] = {}
    for s, count_s in left.items():
        for t, count_t in right.items():
            product_type = tuple(a * b for a, b in zip(s, t))
            result[product_type] = result.get(product_type, 0) + count_s * count_t
    return result


def check_selector_fourier() -> int:
    checks = 0
    for p in range(1, 7):
        active = active_sets(p)
        assert all(len(s) % 2 == 1 for s in active)
        denominator = 2**p
        for a in product((-1, 1), repeat=p):
            reconstructed = Fraction(0)
            for s in active:
                character = 1
                for i in s:
                    character *= a[i]
                reconstructed += Fraction(fourier_numerator(p, s), denominator) * character
            assert reconstructed == selector(a)
            assert selector(tuple(-x for x in a)) == -selector(a)
            if sum(a):
                assert selector(a) == (1 if sum(a) > 0 else -1)
            checks += 1
    return checks


def antipodal_benchmark(p: int = 4) -> int:
    points = np.asarray(list(product((-1, 1), repeat=p)), dtype=np.int64)
    n = len(points)
    index = {tuple(a): i for i, a in enumerate(points)}
    antipode = np.zeros((n, n), dtype=np.int64)
    for i, a in enumerate(points):
        antipode[i, index[tuple(-a)]] = 1
    h = -antipode
    r = 1
    ports = points.T.copy()
    active = active_sets(p)
    assert np.array_equal(h, h.T)
    assert np.array_equal(h @ h, np.eye(n, dtype=np.int64))
    assert int(np.trace(h)) == 0
    assert_closed(h, r, ports, active)

    xs = cube(n)
    quadratic = np.einsum("bi,ij,bj->b", xs, h, xs, optimize=True) // 2
    m = 2
    checks = 0
    for eps in product((-1, 1), repeat=p):
        z = np.asarray(eps, dtype=np.int64) @ ports
        witness = selector_witness(ports, eps)
        assert np.array_equal(h @ witness, witness)
        assert int(z @ witness) == int(np.abs(z).sum())
        exact = int(np.max(np.abs(quadratic) + m * (xs @ z)))
        predicted = n // 2 + m * int(np.abs(z).sum())
        assert exact == predicted
        assert np.array_equal(h @ z, z)
        checks += 1

    # Two copies retain closure and add their labelled field responses.
    h_block = np.block([[h, np.zeros_like(h)], [np.zeros_like(h), h]])
    ports_block = np.concatenate([ports, ports], axis=1)
    assert_closed(h_block, r, ports_block, active)
    for eps in product((-1, 1), repeat=p):
        z = np.asarray(eps, dtype=np.int64) @ ports
        z_block = np.asarray(eps, dtype=np.int64) @ ports_block
        assert int(np.abs(z_block).sum()) == 2 * int(np.abs(z).sum())
        checks += 1
    return checks


def h16_seed() -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    r, h, _ = regular_hadamard(2)
    xs = cube(len(h))
    tops = xs[np.all(xs @ h.T == r * xs, axis=1)]
    reps = tops[tops[:, 0] == 1]
    ports = reps[[0, 1, 6]]
    triple = np.prod(ports, axis=0)
    assert np.array_equal(h @ triple, r * triple)
    return r, h, ports, tops


def check_hadamard_seed_and_tensor() -> int:
    r, h, ports, tops = h16_seed()
    n = len(h)
    active3 = active_sets(3)
    assert_closed(h, r, ports, active3)

    xs = cube(n)
    q = np.einsum("bi,ij,bj->b", xs, h, xs, optimize=True) // 2
    checks = 0
    for eps in product((-1, 1), repeat=3):
        z = np.asarray(eps, dtype=np.int64) @ ports
        witness = selector_witness(ports, eps)
        assert np.array_equal(h @ witness, r * witness)
        assert int(z @ witness) == int(np.abs(z).sum())
        exact = int(np.max(np.abs(q) + (xs @ z)))
        assert exact == r * n // 2 + int(np.abs(z).sum())
        checks += 1

    # The base odd products form an affine rank-two multiplicative coset.
    a, b, c = ports
    c1 = a * b
    c2 = a * c
    group0 = [np.ones(n, dtype=np.int64), c1, c2, c1 * c2]
    affine0 = {tuple((a * g).tolist()) for g in group0}
    odd0 = {
        tuple(port_product(ports, s).tolist())
        for s in subsets(3) if len(s) % 2 == 1
    }
    assert affine0 == odd0 and len(affine0) == 4

    # Coordinatewise tensor of two closed triples is again closed.
    h2 = np.kron(h, h)
    tensor_ports = np.asarray([np.kron(w, w) for w in ports])
    assert_closed(h2, r * r, tensor_ports, active3)
    histogram0 = projective_histogram(ports)
    assert projective_histogram(tensor_ports) == histogram_convolution(
        histogram0, histogram0
    )

    # Common-pole amplification is convolution with a point mass, hence
    # simply multiplies every histogram count by the tail order.
    one = np.ones(n, dtype=np.int64)
    pole_ports = np.asarray([np.kron(w, one) for w in ports])
    assert_closed(h2, r * r, pole_ports, active3)
    assert projective_histogram(pole_ports) == {
        row_type: n * count for row_type, count in histogram0.items()
    }

    # Growing family at j=2: base pole plus four relative generator ports.
    a2 = np.kron(a, a)
    generators = [
        np.kron(c1, one), np.kron(c2, one),
        np.kron(one, c1), np.kron(one, c2),
    ]
    ports2 = np.asarray([a2] + [a2 * g for g in generators])
    p2 = len(ports2)
    assert p2 == 5

    group2 = set()
    for bits in product((0, 1), repeat=4):
        g = np.ones(n * n, dtype=np.int64)
        for bit, generator in zip(bits, generators):
            if bit:
                g *= generator
        group2.add(tuple(g.tolist()))
    assert len(group2) == 16
    affine2 = {tuple((a2 * np.asarray(g, dtype=np.int64)).tolist()) for g in group2}
    odd2 = {
        tuple(port_product(ports2, s).tolist())
        for s in subsets(p2) if len(s) % 2 == 1
    }
    assert affine2 == odd2 and len(odd2) == 16
    for w in odd2:
        vector = np.asarray(w, dtype=np.int64)
        assert np.array_equal(h2 @ vector, (r * r) * vector)

    active5 = active_sets(5)
    assert_closed(h2, r * r, ports2, active5)
    for eps in product((-1, 1), repeat=p2):
        z = np.asarray(eps, dtype=np.int64) @ ports2
        witness = selector_witness(ports2, eps)
        assert np.array_equal(h2 @ witness, (r * r) * witness)
        assert int(z @ witness) == int(np.abs(z).sum())
        checks += 1

    histogram = projective_histogram(ports2)
    assert len(histogram) <= 2 ** (p2 - 1) == 16
    assert sum(histogram.values()) == n * n
    m2 = floor((r * r) / p2)
    assert m2 == 3 and m2 * p2 <= r * r
    checks += 20

    # The common-Hadamard four-port collision is outside tau closure.
    plus = tops[[0, 1, 2, 4]]
    minus = tops[[0, 1, 2, 5]]
    active4 = active_sets(4)
    assert all(len(s) in (1, 3) for s in active4)
    failures = []
    for tuple_ports in (plus, minus):
        failures.append(
            [s for s in active4
             if not np.array_equal(h @ port_product(tuple_ports, s),
                                   r * port_product(tuple_ports, s))]
        )
    assert failures[0] and failures[1]
    assert any(len(s) == 3 for s in failures[0])
    assert any(len(s) == 3 for s in failures[1])
    checks += 2
    return checks


def main() -> None:
    selector_checks = check_selector_fourier()
    antipodal_checks = antipodal_benchmark()
    hadamard_checks = check_hadamard_seed_and_tensor()
    print(
        "Boolean product-algebra closure verified: "
        f"selector={selector_checks}, antipodal={antipodal_checks}, "
        f"Hadamard/tensor={hadamard_checks}"
    )


if __name__ == "__main__":
    main()
