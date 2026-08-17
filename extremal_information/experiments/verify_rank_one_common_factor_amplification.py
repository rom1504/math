#!/usr/bin/env python3
"""Exact diagnostics for the rank-one common-factor amplification theorem."""

from __future__ import annotations

from itertools import product
from math import log, sqrt

import numpy as np

from verify_bcx_two_port_holonomy import regular_hadamard


def word(text: str) -> np.ndarray:
    assert len(text) == 16 and set(text) <= {"+", "-"}
    return np.asarray([1 if symbol == "+" else -1 for symbol in text], dtype=np.int64)


W_MINUS = word("----+--++--+----")
W_PLUS = word("+--+---++-------")
TOP_ORTHOGONAL = word("---+--+-+-++-+++")


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def response(h: np.ndarray, w: np.ndarray, multiplicity: int, spins: np.ndarray) -> int:
    quadratic = np.abs(np.einsum("bi,ij,bj->b", spins, h, spins)) // 2
    field = np.abs(spins @ w)
    return int(np.max(quadratic + multiplicity * field))


def verify_seed() -> tuple[np.ndarray, np.ndarray, int]:
    r, h, _ = regular_hadamard(2)
    spins = cube(16)
    assert r == 4
    assert np.array_equal(h @ h, 16 * np.eye(16, dtype=np.int64))
    assert np.array_equal(h @ np.ones(16, dtype=np.int64), 4 * np.ones(16, dtype=np.int64))
    assert int(np.trace(h)) == 0
    assert int(W_MINUS @ h @ W_MINUS) == int(W_PLUS @ h @ W_PLUS) == 0
    assert response(h, W_MINUS, r, spins) == 64
    assert response(h, W_PLUS, r, spins) == 78
    return h, spins, r


def verify_exact_contraction(h: np.ndarray, spins: np.ndarray) -> int:
    rng = np.random.default_rng(1729)
    checks = 0
    for k in (1, 2, 3, 5, 8):
        ones = np.ones(k, dtype=np.int64)
        j = np.ones((k, k), dtype=np.int64)
        hk = np.kron(h, j)
        for w in (W_MINUS, W_PLUS):
            wk = np.kron(w, ones)
            for _ in range(20):
                x = rng.choice((-1, 1), size=16 * k).astype(np.int64)
                fibre_sums = x.reshape(16, k).sum(axis=1)
                y = fibre_sums / k
                left = abs(int(x @ hk @ x)) / 2 + 4 * k * abs(int(wk @ x))
                right = k * k * (abs(float(y @ h @ y)) / 2 + 4 * abs(float(w @ y)))
                assert abs(left - right) < 1e-9
                checks += 1

        # Every base Boolean vertex occurs in the fibre grid.
        base_values_plus = (
            np.abs(np.einsum("bi,ij,bj->b", spins, h, spins)) // 2
            + 4 * np.abs(spins @ W_PLUS)
        )
        maximizer = spins[int(np.argmax(base_values_plus))]
        lifted = np.kron(maximizer, ones)
        assert abs(int(lifted @ hk @ lifted)) // 2 + 4 * k * abs(
            int(np.kron(W_PLUS, ones) @ lifted)
        ) == 78 * k * k

        # The equal state and the analytic 6 k^2 separation certificate.
        assert int(np.kron(W_MINUS, ones) @ np.kron(W_MINUS, ones)) == 16 * k
        assert int(np.kron(W_PLUS, ones) @ np.kron(W_PLUS, ones)) == 16 * k
        assert int(np.kron(W_MINUS, ones) @ hk @ np.kron(W_MINUS, ones)) == 0
        assert int(np.kron(W_PLUS, ones) @ hk @ np.kron(W_PLUS, ones)) == 0
        assert 78 * k * k - (64 + 8) * k * k == 6 * k * k
        assert 6 * k * k * 32 == 3 * (4 * k) * (16 * k)
        checks += 7
    return checks


def verify_rounding_identity(h: np.ndarray) -> int:
    """Check the channel identity behind the deterministic upper bound."""

    rng = np.random.default_rng(2718)
    checks = 0
    assert int(np.trace(h)) == 0
    for _ in range(200):
        y = rng.uniform(-1, 1, size=16)
        for sigma in (-1, 1):
            for epsilon in (-1, 1):
                # E[Z_i Z_j]=y_i y_j off the diagonal and E[Z_i^2]=1.
                expected_quadratic = float(y @ h @ y) + float(
                    np.sum(np.diag(h) * (1 - y * y))
                )
                expected_channel = (
                    sigma * expected_quadratic / 2 + 4 * epsilon * float(W_MINUS @ y)
                )
                channel_at_mean = sigma * float(y @ h @ y) / 2 + 4 * epsilon * float(
                    W_MINUS @ y
                )
                correction = sigma * float(np.sum(np.diag(h) * y * y)) / 2
                assert abs(channel_at_mean - (expected_channel + correction)) < 1e-10
                assert abs(correction) <= 8 + 1e-12
                checks += 2
    return checks


def verify_spherical_value(h: np.ndarray) -> int:
    checks = 0
    j = h / 4
    for w in (W_MINUS, W_PLUS):
        z = j @ w
        assert abs(float(w @ z)) < 1e-12
        assert abs(float(z @ z) - 16) < 1e-12
        for sigma in (-1, 1):
            u = (sqrt(3) * w + sigma * z) / 2
            assert abs(float(u @ u) - 16) < 1e-10
            value = sigma * float(u @ h @ u) / 2 + 4 * float(w @ u)
            assert abs(value - 3 * sqrt(3) * 4 * 16 / 4) < 1e-9
            checks += 3
    return checks


def verify_hidden_mode(h: np.ndarray) -> int:
    # The first strict regular-Hadamard common factor already contains an
    # orthogonal Boolean top pole.
    k = h
    a = np.ones(16, dtype=np.int64)
    b = TOP_ORTHOGONAL
    assert np.array_equal(k @ a, 4 * a)
    assert np.array_equal(k @ b, 4 * b)
    assert int(a @ b) == 0

    tensor = np.kron(h, k)
    port = np.kron(W_MINUS, a)
    x_one = np.kron(np.ones(16, dtype=np.int64), b)
    x_zero = np.kron(W_MINUS, b)
    for x in (x_one, x_zero):
        fibre_magnetization = x.reshape(16, 16) @ a
        assert np.array_equal(fibre_magnetization, np.zeros(16, dtype=np.int64))
        assert int(port @ x) == 0
    assert abs(int(x_one @ tensor @ x_one)) // 2 == (16 * 256) // 2
    assert int(x_zero @ tensor @ x_zero) // 2 == 0
    return 9


def verify_completion_arithmetic() -> int:
    checks = 0
    for k in (64, 100, 256, 1000):
        q_bound = 2 * (4 * k) ** 1.5
        retained = 6 * k * k - 2 * q_bound
        assert retained > 0
        log_probability_bound = (4 * k + 1) * log(2) - 16 * k
        assert log_probability_bound < 0
        checks += 2
    return checks


def main() -> None:
    h, spins, _ = verify_seed()
    checks = 8
    checks += verify_exact_contraction(h, spins)
    checks += verify_rounding_identity(h)
    checks += verify_spherical_value(h)
    checks += verify_hidden_mode(h)
    checks += verify_completion_arithmetic()
    print(f"rank-one amplification checks passed: {checks}")


if __name__ == "__main__":
    main()
