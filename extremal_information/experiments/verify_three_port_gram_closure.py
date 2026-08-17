#!/usr/bin/env python3
"""Exact checks for drafts/three_port_gram_closure.md."""

from __future__ import annotations

from itertools import combinations, product
from math import isclose, sqrt
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_bcx_two_port_holonomy import regular_hadamard  # noqa: E402


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def word(text: str) -> np.ndarray:
    return np.asarray([1 if c == "+" else -1 for c in text], dtype=np.int64)


def boolean_response(h: np.ndarray, ports: np.ndarray, m: int,
                     xs: np.ndarray) -> int:
    quadratics = np.abs(
        np.einsum("bi,ij,bj->b", xs, h, xs, optimize=True) // 2
    )
    fields = np.abs(xs @ ports.T).sum(axis=1)
    return int(np.max(quadratics + m * fields))


def endpoint_supports(ports: np.ndarray) -> list[int]:
    p = len(ports)
    return sorted(
        int(np.abs(np.asarray(eps, dtype=np.int64) @ ports).sum())
        for eps in product((-1, 1), repeat=p)
    )


def top_vectors(h: np.ndarray, r: int, xs: np.ndarray) -> np.ndarray:
    return xs[np.all(xs @ h.T == r * xs, axis=1)]


def check_two_port_formula(h: np.ndarray, r: int, xs: np.ndarray) -> int:
    n = len(h)
    tops = top_vectors(h, r, xs)
    reps = tops[tops[:, 0] == 1]
    checks = 0
    for ia, ib in combinations(range(len(reps)), 2):
        ports = reps[[ia, ib]]
        dot = int(ports[0] @ ports[1])
        for m in (0, 1, 2):
            predicted = r * n // 2 + m * (n + abs(dot))
            assert boolean_response(h, ports, m, xs) == predicted
            spherical = r * n / 2 + m * n * sqrt(2 * (1 + abs(dot) / n))
            z = ports[0] + (1 if dot >= 0 else -1) * ports[1]
            assert np.array_equal(h @ z, r * z)
            assert isclose(spherical, r * n / 2 + m * sqrt(n) * np.linalg.norm(z))
            checks += 1
    return checks


def check_three_port_closure(h: np.ndarray, r: int, xs: np.ndarray) -> int:
    n = len(h)
    tops = top_vectors(h, r, xs)
    reps = tops[tops[:, 0] == 1]
    ports = reps[[0, 1, 6]]
    triple = np.prod(ports, axis=0)
    assert np.array_equal(h @ triple, r * triple)
    assert "".join("+" if a == 1 else "-" for a in triple) == "+++++--++--+++++"

    gram = ports @ ports.T
    assert np.array_equal(
        gram,
        np.asarray([[16, 8, 0], [8, 16, -8], [0, -8, 16]], dtype=np.int64),
    )
    pair_max = max(
        sum(gram[i, j] * eps[i] * eps[j]
            for i in range(3) for j in range(i + 1, 3))
        for eps in product((-1, 1), repeat=3)
    )
    assert pair_max == n  # T(G)=1.

    checks = 0
    for eps in product((-1, 1), repeat=3):
        epsv = np.asarray(eps, dtype=np.int64)
        z = epsv @ ports
        majority = (
            eps[0] * ports[0] + eps[1] * ports[1] + eps[2] * ports[2]
            - eps[0] * eps[1] * eps[2] * triple
        ) // 2
        assert set(majority.tolist()) <= {-1, 1}
        assert np.array_equal(majority, np.sign(z))
        assert np.array_equal(h @ majority, r * majority)
        pair_sum = sum(
            gram[i, j] * eps[i] * eps[j]
            for i in range(3) for j in range(i + 1, 3)
        )
        assert int(np.abs(z).sum()) == (3 * n + pair_sum) // 2
        assert int(z @ z) == 3 * n + 2 * pair_sum
        checks += 1

    for m in (0, 1, 2, 4):
        predicted = r * n // 2 + 2 * m * n
        assert boolean_response(h, ports, m, xs) == predicted
        spherical = r * n / 2 + m * n * sqrt(5)
        assert spherical >= predicted
        checks += 1

    # The tensor lift preserves product closure and normalized Gram data.
    h2 = np.kron(h, h)
    one = np.ones(n, dtype=np.int64)
    ports2 = np.asarray([np.kron(w, one) for w in ports])
    triple2 = np.prod(ports2, axis=0)
    assert np.array_equal(h2 @ ports2.T, (r * r) * ports2.T)
    assert np.array_equal(h2 @ triple2, (r * r) * triple2)
    assert np.array_equal(ports2 @ ports2.T, n * gram)
    checks += 3
    return checks


def check_four_port_collision(h: np.ndarray, r: int, xs: np.ndarray) -> int:
    n = len(h)
    tops = top_vectors(h, r, xs)  # Full lexicographic list, including antipodes.
    assert len(tops) == 20
    plus = tops[[0, 1, 2, 4]]
    minus = tops[[0, 1, 2, 5]]

    expected_words = (
        "----------------",
        "-----++--++-----",
        "---+--+-+-++-+++",
        "---+++-+-+---+++",
        "--+----+-++++-++",
    )
    actual_words = tuple(
        "".join("+" if a == 1 else "-" for a in w)
        for w in (tops[0], tops[1], tops[2], tops[4], tops[5])
    )
    assert actual_words == expected_words

    raw_gram = np.asarray(
        [[16, 8, 0, 0], [8, 16, 0, 0],
         [0, 0, 16, 0], [0, 0, 0, 16]], dtype=np.int64
    )
    assert np.array_equal(plus @ plus.T, raw_gram)
    assert np.array_equal(minus @ minus.T, raw_gram)
    # Every port is a +r eigenvector, so raw Rayleigh and Gram tables agree.
    assert np.array_equal(plus @ h @ plus.T, r * raw_gram)
    assert np.array_equal(minus @ h @ minus.T, r * raw_gram)

    supports_plus = endpoint_supports(plus)
    supports_minus = endpoint_supports(minus)
    assert supports_plus == [16] * 4 + [24] * 8 + [32] * 4
    assert supports_minus == [20] * 8 + [28] * 8

    # The Euclidean-maximizing channels have eps_0 eps_1=+1, hence norm
    # sqrt(5n).  Their best l1 values expose the extra flatness coordinate.
    euclidean_plus = []
    euclidean_minus = []
    for eps in product((-1, 1), repeat=4):
        if eps[0] == eps[1]:
            epsv = np.asarray(eps, dtype=np.int64)
            zp = epsv @ plus
            zm = epsv @ minus
            assert int(zp @ zp) == 5 * n
            assert int(zm @ zm) == 5 * n
            euclidean_plus.append(int(np.abs(zp).sum()))
            euclidean_minus.append(int(np.abs(zm).sum()))
    assert max(euclidean_plus) == 2 * n
    assert set(euclidean_minus) == {7 * n // 4}

    exact = {1: (56, 56), 2: (88, 82), 4: (152, 138)}
    for m, values in exact.items():
        assert boolean_response(h, plus, m, xs) == values[0]
        assert boolean_response(h, minus, m, xs) == values[1]

    witness = word("------+---+-----")
    assert int(np.abs(plus @ witness).sum()) == 32
    assert int(abs(witness @ h @ witness) // 2) == 24

    # Scalable tensor certificate at j=2.  No 2^256 enumeration is used.
    tail_n = n
    h2 = np.kron(h, h)
    one = np.ones(tail_n, dtype=np.int64)
    plus2 = np.asarray([np.kron(w, one) for w in plus])
    minus2 = np.asarray([np.kron(w, one) for w in minus])
    n2 = n * tail_n
    r2 = r * r
    assert np.array_equal(h2 @ plus2.T, r2 * plus2.T)
    assert np.array_equal(h2 @ minus2.T, r2 * minus2.T)
    assert np.array_equal(plus2 @ plus2.T, tail_n * raw_gram)
    assert np.array_equal(minus2 @ minus2.T, tail_n * raw_gram)
    assert max(endpoint_supports(plus2)) == 2 * n2
    assert max(endpoint_supports(minus2)) == 7 * n2 // 4

    witness2 = np.kron(witness, one)
    assert int(np.abs(plus2 @ witness2).sum()) == 2 * n2
    assert int(abs(witness2 @ h2 @ witness2) // 2) == 3 * r2 * n2 // 8
    lower_plus = 3 * r2 * n2 // 8 + r2 * 2 * n2
    upper_minus = r2 * n2 // 2 + r2 * 7 * n2 // 4
    assert lower_plus - upper_minus == r2 * n2 // 8
    return 25


def main() -> None:
    r, h, _ = regular_hadamard(2)
    n = len(h)
    xs = cube(n)
    assert r == 4 and n == 16
    assert np.array_equal(h @ h, n * np.eye(n, dtype=np.int64))
    assert np.array_equal(h @ np.ones(n, dtype=np.int64), r * np.ones(n, dtype=np.int64))
    assert int(np.trace(h)) == 0

    two = check_two_port_formula(h, r, xs)
    three = check_three_port_closure(h, r, xs)
    four = check_four_port_collision(h, r, xs)
    print(
        "three-port Gram closure verified: "
        f"two-port checks={two}, three-port checks={three}, "
        f"four-port/tensor checks={four}"
    )


if __name__ == "__main__":
    main()
