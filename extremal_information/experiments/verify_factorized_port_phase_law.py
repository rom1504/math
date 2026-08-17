#!/usr/bin/env python3
"""Finite checks for drafts/factorized_port_phase_law.md."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import sqrt
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_bcx_two_port_holonomy import regular_hadamard  # noqa: E402


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def seed() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r, h, _ = regular_hadamard(2)
    xs = cube(len(h))
    tops = xs[np.all(xs @ h.T == r * xs, axis=1)]
    reps = tops[tops[:, 0] == 1]
    a, b, c = reps[[0, 1, 6]]
    assert r == 4 and len(h) == 16
    assert int(np.trace(h)) == 0
    assert np.array_equal(h @ h, 16 * np.eye(16, dtype=np.int64))
    for w in (a, b, c, a * b * c):
        assert np.array_equal(h @ w, 4 * w)
    assert int(a @ b) == 8
    assert int(a @ c) == 0
    return h, a, b, c


def check_cartesian_closure() -> int:
    h, a, b, c = seed()
    one = np.ones(16, dtype=np.int64)
    gx, gy = a * b, a * c
    base = np.kron(a, a)
    generators = (
        np.kron(gx, one),
        np.kron(gy, one),
        np.kron(one, gx),
        np.kron(one, gy),
    )
    ports = (base,) + tuple(base * g for g in generators)
    h2 = np.kron(h, h)
    checks = 0
    # Every odd port product is one Cartesian affine pole and is +16.
    for bits in product((0, 1), repeat=len(ports)):
        if sum(bits) % 2 != 1:
            continue
        w = np.ones(256, dtype=np.int64)
        for bit, port in zip(bits, ports):
            if bit:
                w *= port
        assert np.array_equal(h2 @ w, 16 * w)
        checks += 1
    assert checks == 16
    return checks


def exact_block_law() -> list[tuple[int, int, Fraction]]:
    _, a, b, c = seed()
    gx, gy = a * b, a * c
    counts: dict[tuple[int, int], int] = {}
    for x, y in zip(gx.tolist(), gy.tolist()):
        counts[(x, y)] = counts.get((x, y), 0) + 1
    assert counts == {(1, 1): 4, (1, -1): 8, (-1, 1): 4}
    assert sum(x * count for (x, _), count in counts.items()) == 8
    assert sum(y * count for (_, y), count in counts.items()) == 0
    return [(x, y, Fraction(count, 16)) for (x, y), count in counts.items()]


def support_for_blocks(
    blocks: list[list[tuple[tuple[int, ...], Fraction]]],
) -> tuple[Fraction, int]:
    """Exhaust exact max E|eps0+sum eps.block| for small blocks."""
    q_total = sum(len(block[0][0]) for block in blocks)
    best = Fraction(-1)
    for eps in product((-1, 1), repeat=q_total + 1):
        eps0 = eps[0]
        split: list[tuple[int, ...]] = []
        at = 1
        for block in blocks:
            q = len(block[0][0])
            split.append(eps[at:at + q])
            at += q
        expectation = Fraction(0)
        for outcomes in product(*blocks):
            total = eps0
            probability = Fraction(1)
            for signs_and_probability, local_eps in zip(outcomes, split):
                signs, local_probability = signs_and_probability
                total += sum(e * x for e, x in zip(local_eps, signs))
                probability *= local_probability
            expectation += probability * abs(total)
        best = max(best, expectation)
    return best, q_total + 1


def check_first_moment_bound() -> int:
    pair_law = exact_block_law()
    block = [((x, y), probability) for x, y, probability in pair_law]
    checks = 0
    for length in range(1, 5):
        maximum, p = support_for_blocks([block] * length)
        first_moment = Fraction(1) + Fraction(length, 2)
        variance_roof = 2 * sqrt(length)  # sqrt(sum q_t^2), q_t=2.
        assert maximum >= first_moment
        assert float(maximum - first_moment) <= variance_roof + 1e-12
        theta = first_moment / p
        assert abs(float(maximum / p - theta)) <= variance_roof / p + 1e-12
        checks += 1
    return checks


def one_bit_law(kind: str) -> list[tuple[tuple[int, ...], Fraction]]:
    if kind == "X":
        return [((1,), Fraction(3, 4)), ((-1,), Fraction(1, 4))]
    if kind == "Y":
        return [((1,), Fraction(1, 2)), ((-1,), Fraction(1, 2))]
    raise ValueError(kind)


def check_mixed_factor_phase() -> int:
    checks = 0
    for word in ("X", "Y", "XY", "XXY", "XYYX"):
        blocks = [one_bit_law(kind) for kind in word]
        maximum, p = support_for_blocks(blocks)
        k = word.count("X")
        first_moment = Fraction(1) + Fraction(k, 2)
        assert maximum >= first_moment
        assert float(maximum - first_moment) <= sqrt(len(word)) + 1e-12
        assert p == len(word) + 1
        checks += 1
    return checks


def check_normalization_and_blocks() -> int:
    checks = 0
    # Repeated two-generator seed: the explicit limiting phase is 1/4 and
    # the completed-parent normalization tends to 3/4.
    for length in (10, 100, 1000):
        p = 2 * length + 1
        n = 16**length
        r = 4**length
        m = r // p
        c_mass = Fraction(m * p, r)
        theta = (Fraction(1) + Fraction(length, 2)) / p
        predicted = Fraction(1, 2) + c_mass * theta
        assert abs(float(theta) - 0.25) <= 1 / p
        assert abs(float(predicted) - 0.75) <= 2 / p
        auxiliary = m * p
        assert auxiliary <= r
        assert Fraction(auxiliary * auxiliary, r * n) <= Fraction(1, r)
        checks += 1

    # Alternating dominant blocks force X densities near one and zero.
    total = 1
    x_count = 1
    for k in range(2, 18):
        block = k * k * total
        if k % 2 == 1:  # X block
            x_count += block
        total += block
        density = Fraction(x_count, total)
        if k % 2 == 1:
            assert density >= Fraction(k * k, k * k + 1)
        else:
            assert density <= Fraction(1, k * k + 1)
        theta = (Fraction(1) + Fraction(x_count, 2)) / (total + 1)
        if k % 2 == 1:
            assert float(theta) > 0.5 - 1 / k
        else:
            assert float(theta) < 1 / k
        checks += 1
    return checks


def main() -> None:
    closure = check_cartesian_closure()
    moment = check_first_moment_bound()
    mixed = check_mixed_factor_phase()
    asymptotic = check_normalization_and_blocks()
    total = closure + moment + mixed + asymptotic
    print(
        "factorized port phase verification passed: "
        f"closure={closure}, moment={moment}, mixed={mixed}, "
        f"asymptotic={asymptotic}, total={total}"
    )


if __name__ == "__main__":
    main()
