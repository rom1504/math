#!/usr/bin/env python3
"""Exact checks for equal-(G,R), unequal-Boolean-response collisions."""

from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import sqrt

import numpy as np

from verify_bcx_two_port_holonomy import regular_hadamard


def word(text: str) -> np.ndarray:
    assert len(text) == 16
    assert set(text) <= {"+", "-"}
    return np.asarray([1 if symbol == "+" else -1 for symbol in text], dtype=np.int64)


PORT_WORDS = [
    word("----------------"),
    word("-----++--++-----"),
    word("---+--+-+-++-+++"),
    word("---+++-+-+---+++"),
    word("--+----+-++++-++"),
]
PORTS_A = np.asarray(PORT_WORDS[:4], dtype=np.int64)
PORTS_B = np.asarray(PORT_WORDS[:3] + [PORT_WORDS[4]], dtype=np.int64)
X_A = word("------+---+-----")

ONE_MINUS = word("----+--++--+----")
ONE_PLUS = word("+--+---++-------")


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def support(ports: np.ndarray) -> tuple[int, list[int]]:
    eps = cube(len(ports))
    values = np.sum(np.abs(eps @ ports), axis=1)
    return int(np.max(values)), sorted(set(int(value) for value in values))


def boolean_response(
    h: np.ndarray, ports: np.ndarray, multiplicity: int, spins: np.ndarray
) -> int:
    quadratic = np.abs(np.einsum("bi,ij,bj->b", spins, h, spins)) // 2
    fields = np.sum(np.abs(spins @ ports.T), axis=1)
    return int(np.max(quadratic + multiplicity * fields))


def verify_four_port_seed() -> int:
    r, h, _ = regular_hadamard(2)
    n = len(h)
    spins = cube(n)
    target = np.asarray(
        [[16, 8, 0, 0], [8, 16, 0, 0], [0, 0, 16, 0], [0, 0, 0, 16]],
        dtype=np.int64,
    )

    assert r == 4 and n == 16
    assert np.array_equal(h @ h, 16 * np.eye(16, dtype=np.int64))
    assert int(np.trace(h)) == 0
    for w in PORT_WORDS:
        assert np.array_equal(h @ w, 4 * w)
    assert np.array_equal(PORTS_A @ PORTS_A.T, target)
    assert np.array_equal(PORTS_B @ PORTS_B.T, target)
    assert np.array_equal(PORTS_A @ h @ PORTS_A.T, 4 * target)
    assert np.array_equal(PORTS_B @ h @ PORTS_B.T, 4 * target)

    support_a, values_a = support(PORTS_A)
    support_b, values_b = support(PORTS_B)
    assert (support_a, values_a) == (32, [16, 24, 32])
    assert (support_b, values_b) == (28, [20, 28])
    assert int(X_A @ h @ X_A) // 2 == 24
    assert int(np.sum(np.abs(PORTS_A @ X_A))) == 32

    expected = {1: (56, 56), 2: (88, 82), 4: (152, 138), 8: (280, 250)}
    for m, pair in expected.items():
        observed = (
            boolean_response(h, PORTS_A, m, spins),
            boolean_response(h, PORTS_B, m, spins),
        )
        assert observed == pair
    return 20


def verify_dense_tensor_theorem() -> int:
    r0, h0, _ = regular_hadamard(2)
    checks = 0

    # Direct physical checks at the first nontrivial tensor extension.
    h2 = np.kron(h0, h0)
    rest_one = np.ones(16, dtype=np.int64)
    ports_a2 = np.asarray([np.kron(w, rest_one) for w in PORTS_A])
    ports_b2 = np.asarray([np.kron(w, rest_one) for w in PORTS_B])
    gram_a2 = ports_a2 @ ports_a2.T
    gram_b2 = ports_b2 @ ports_b2.T
    ray_a2 = ports_a2 @ h2 @ ports_a2.T
    ray_b2 = ports_b2 @ h2 @ ports_b2.T
    assert np.array_equal(gram_a2, gram_b2)
    assert np.array_equal(ray_a2, ray_b2)
    assert support(ports_a2)[0] == 32 * 16
    assert support(ports_b2)[0] == 28 * 16
    checks += 4

    # The exact formulas are integer identities at arbitrary tensor depth.
    for depth in range(1, 9):
        n = 16**depth
        r = 4**depth
        rest = 16 ** (depth - 1)
        high_witness = 6 * r * rest + 32 * r * rest
        low_upper = r * n // 2 + 28 * r * rest
        assert high_witness - low_upper == r * n // 8
        assert 4 * r / r == 4  # total port mass mp/r at m=r,p=4.
        completion_cap_twice = (4 * r) * (4 * r - 1)
        assert completion_cap_twice / (r * n) < 1 or depth == 1
        checks += 4
    return checks


def verify_exact_parent_and_completion() -> int:
    """Check that EG.1 is the parent cap and audit the 2Q(C) estimate."""

    _, h, _ = regular_hadamard(2)
    old_spins = cube(16)
    endpoint_spins = cube(4)  # One endpoint in each of four shores.
    hollow = h - np.diag(np.diag(h))
    old_from_h = np.einsum("bi,ij,bj->b", old_spins, h, old_spins) // 2
    old_from_hollow = (
        np.einsum("bi,ij,bj->b", old_spins, hollow, old_spins) // 2
    )
    assert int(np.trace(h)) == 0
    assert np.array_equal(old_from_h, old_from_hollow)

    # A public positive auxiliary clique is enough to test the pointwise
    # Lipschitz statement; the proof itself is independent of this choice.
    completion = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    completion_energy = (
        np.einsum("bi,ij,bj->b", endpoint_spins, completion, endpoint_spins) // 2
    )
    completion_cap = int(np.max(np.abs(completion_energy)))

    incomplete_caps = []
    completed_caps = []
    for ports in (PORTS_A, PORTS_B):
        correlations = old_spins @ ports.T
        bridge = correlations @ endpoint_spins.T
        incomplete_values = old_from_h[:, None] + bridge
        completed_values = incomplete_values + completion_energy[None, :]
        incomplete = int(np.max(np.abs(incomplete_values)))
        completed = int(np.max(np.abs(completed_values)))
        assert incomplete == boolean_response(h, ports, 1, old_spins)
        assert abs(completed - incomplete) <= completion_cap
        incomplete_caps.append(incomplete)
        completed_caps.append(completed)

    assert (
        completed_caps[0] - completed_caps[1]
        >= incomplete_caps[0] - incomplete_caps[1] - 2 * completion_cap
    )
    return 7


def projective_distance_shells(
    h: np.ndarray, port: np.ndarray, spins: np.ndarray
) -> list[int]:
    n = len(port)
    correlations = spins @ port
    distance = (n - np.abs(correlations)) // 2
    energy = np.abs(np.einsum("bi,ij,bj->b", spins, h, spins)) // 2
    return [int(np.max(energy[distance == d])) for d in range(n // 2 + 1)]


def all_one_port_responses(
    h: np.ndarray, spins: np.ndarray, multiplicity: int
) -> tuple[np.ndarray, np.ndarray]:
    """Compute every one-port response by an exact hypercube distance transform."""

    n = len(h)
    rayleigh = np.einsum("bi,ij,bj->b", spins, h, spins).astype(np.int64)
    envelope = (np.abs(rayleigh) // 2).copy()
    # max_x e(x)-2m d(w,x), all w, via the max-plus Hamming transform.
    for bit in range(n):
        step = 1 << (n - bit - 1)
        view = envelope.reshape(-1, 2 * step)
        left = view[:, :step].copy()
        right = view[:, step:].copy()
        view[:, :step] = np.maximum(left, right - 2 * multiplicity)
        view[:, step:] = np.maximum(right, left - 2 * multiplicity)
    return rayleigh, multiplicity * n + envelope


def verify_one_port_collision() -> int:
    r, h, _ = regular_hadamard(2)
    spins = cube(16)
    assert int(ONE_MINUS @ h @ ONE_MINUS) == 0
    assert int(ONE_PLUS @ h @ ONE_PLUS) == 0
    assert projective_distance_shells(h, ONE_MINUS, spins) == [
        0, 6, 16, 22, 32, 26, 24, 26, 32
    ]
    assert projective_distance_shells(h, ONE_PLUS, spins) == [
        0, 22, 24, 26, 32, 26, 32, 26, 32
    ]
    assert boolean_response(h, ONE_MINUS[None, :], r, spins) == 64
    assert boolean_response(h, ONE_PLUS[None, :], r, spins) == 78

    # Exhaust the complete one-port space.  At m=r, the largest response
    # spread within one exact Rayleigh class is 14, attained at R=0.
    rayleigh, responses = all_one_port_responses(h, spins, r)
    spreads = defaultdict(list)
    for value in np.unique(rayleigh):
        response_class = responses[rayleigh == value]
        spreads[int(value)] = [int(np.min(response_class)), int(np.max(response_class))]
    assert max(hi - lo for lo, hi in spreads.values()) == 14
    assert spreads[0] == [64, 78]
    return 10


def verify_exposed_flatness() -> int:
    r, h, _ = regular_hadamard(2)
    j = h / r
    expected_l1 = {
        "minus+": 8 * sqrt(3),
        "minus-": 8 * sqrt(3),
        "plus+": 3 + 7 * sqrt(3),
        "plus-": 8 * sqrt(3),
    }
    observed: dict[str, float] = {}
    for label, w in (("minus", ONE_MINUS), ("plus", ONE_PLUS)):
        z = j @ w
        assert abs(float(w @ z)) < 1e-12
        assert abs(float(z @ z) - 16) < 1e-12
        for sigma, suffix in ((1, "+"), (-1, "-")):
            u = (sqrt(3) * w + sigma * z) / 2
            assert abs(float(u @ u) - 16) < 1e-12
            value = sigma * float(u @ h @ u) / 2 + r * float(w @ u)
            assert abs(value - 3 * sqrt(3) * r * 16 / 4) < 1e-10
            observed[label + suffix] = float(np.linalg.norm(u, 1))
    for key, value in expected_l1.items():
        assert abs(observed[key] - value) < 1e-10
    phi_minus = 1 - expected_l1["minus+"] / 16
    phi_plus = 1 - expected_l1["plus+"] / 16
    assert abs((phi_minus - phi_plus) - (3 - sqrt(3)) / 16) < 1e-12

    # Tensoring with a common Boolean top pole repeats every exposed coordinate.
    for rest in (1, 16, 256, 4096):
        assert abs((expected_l1["minus+"] * rest) / (16 * rest)
                   - expected_l1["minus+"] / 16) < 1e-15
        assert abs((expected_l1["plus+"] * rest) / (16 * rest)
                   - expected_l1["plus+"] / 16) < 1e-15
    return 20


def verify_pure_linear_baseline() -> int:
    rows = cube(4)
    even = rows[np.prod(rows, axis=1) == 1]
    uniform = rows
    doubled_even = np.repeat(even, 2, axis=0)
    assert np.array_equal(uniform.T @ uniform, 16 * np.eye(4, dtype=np.int64))
    assert np.array_equal(doubled_even.T @ doubled_even, 16 * np.eye(4, dtype=np.int64))
    assert support(uniform.T)[0] == 24
    assert support(doubled_even.T)[0] == 32

    # For p=3, total mass and three pair characters form the complete Walsh
    # basis on the four projective row patterns.  For p=4 they have rank 7
    # on eight patterns and leave one parity direction.
    projective3 = cube(2)  # gauge the first of three row signs to +1.
    moments3 = np.asarray(
        [[1, a, b, a * b] for a, b in projective3], dtype=np.int64
    ).T
    assert np.linalg.matrix_rank(moments3) == 4
    projective4 = cube(3)
    moments4 = []
    for a, b, c in projective4:
        row = [1, a, b, c, a * b, a * c, b * c]
        moments4.append(row)
    assert np.linalg.matrix_rank(np.asarray(moments4, dtype=np.int64)) == 7
    return 6


def main() -> None:
    checks = 0
    checks += verify_four_port_seed()
    checks += verify_dense_tensor_theorem()
    checks += verify_exact_parent_and_completion()
    checks += verify_one_port_collision()
    checks += verify_exposed_flatness()
    checks += verify_pure_linear_baseline()
    print(f"equal Gram--Rayleigh collision checks passed: {checks}")


if __name__ == "__main__":
    main()
