#!/usr/bin/env python3
"""Small exact/numerical checks for approximate Boolean recovery.

The convergence proofs are topological.  This script checks the finite
algebra used by the lift/retraction certificate and the draft's falsifiers.
"""

from __future__ import annotations

from itertools import product
from math import log, sin

import numpy as np


def boolean_vectors(n: int) -> list[np.ndarray]:
    return [np.asarray(x, dtype=np.float64) for x in product((-1, 1), repeat=n)]


def boolean_bilinear_norm(matrix: np.ndarray) -> float:
    vectors = boolean_vectors(len(matrix))
    return max(abs(float(x @ matrix @ y)) for x in vectors for y in vectors)


def quadratic_extrema(matrix: np.ndarray) -> tuple[float, float]:
    values = [float(x @ matrix @ x) for x in boolean_vectors(len(matrix))]
    return max(values), max(abs(value) for value in values)


def block_sum(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    answer = np.zeros(
        (len(left) + len(right), len(left) + len(right)), dtype=np.float64
    )
    answer[: len(left), : len(left)] = left
    answer[len(left) :, len(left) :] = right
    return answer


def verify() -> None:
    checks = 0

    # Exact regular-Hadamard compression is the zero-defect case of AR.27.
    h = np.asarray(
        [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]],
        dtype=np.float64,
    )
    u = np.asarray([1, 1, 1, -1], dtype=np.float64)
    a = np.asarray([[2, -1], [-1, 0]], dtype=np.float64)
    n = len(a)
    t = np.kron(np.eye(n), u.reshape(-1, 1))
    amplified = np.kron(a, h)
    c0 = a / n**1.5
    c1 = amplified / (4 * n) ** 1.5
    d_up = t.T @ c1 @ t - c0
    assert np.allclose(d_up, 0)
    assert boolean_bilinear_norm(d_up) < 1e-12
    checks += 2

    # A symmetric non-tensor perturbation: the compressed-kernel norm is
    # exactly the exhaustive all-pairs lift distortion and obeys the simple
    # operator certificate from AR.28.
    raw = np.arange(1, (4 * n) ** 2 + 1, dtype=np.float64).reshape(4 * n, 4 * n)
    error = 1e-3 * (raw + raw.T)
    perturbed = amplified + error
    c1_perturbed = perturbed / (4 * n) ** 1.5
    d_up = t.T @ c1_perturbed @ t - c0
    matrix_distortion = boolean_bilinear_norm(d_up)
    direct_distortion = 0.0
    for x in boolean_vectors(n):
        for y in boolean_vectors(n):
            direct_distortion = max(
                direct_distortion,
                abs(float((t @ x) @ c1_perturbed @ (t @ y) - x @ c0 @ y)),
            )
    assert abs(matrix_distortion - direct_distortion) < 1e-12
    assert matrix_distortion <= np.linalg.norm(error, 2) / np.sqrt(4 * n) + 1e-12
    assert matrix_distortion <= np.sum(np.abs(d_up)) + 1e-12
    assert matrix_distortion <= n * np.linalg.norm(d_up, 2) + 1e-12
    checks += 4

    # A signed coordinate selector is a literal left inverse of replication;
    # its kernel defect exactly equals the reverse all-pairs distortion.
    s = np.zeros((n, 4 * n), dtype=np.float64)
    for i in range(n):
        s[i, 4 * i] = u[0]
    assert np.array_equal(s @ t, np.eye(n))
    d_down = s.T @ c0 @ s - c1_perturbed
    reverse_matrix_distortion = boolean_bilinear_norm(d_down)
    reverse_direct_distortion = 0.0
    for z in boolean_vectors(4 * n):
        sz = s @ z
        for w in boolean_vectors(4 * n):
            reverse_direct_distortion = max(
                reverse_direct_distortion,
                abs(float(sz @ c0 @ (s @ w) - z @ c1_perturbed @ w)),
            )
    assert abs(reverse_matrix_distortion - reverse_direct_distortion) < 1e-12
    checks += 2

    # AR.32--AR.33: equality of every self-quadratic does not control pairs.
    zero = np.zeros((2, 2), dtype=np.float64)
    trace_zero_diagonal = np.diag([1.0, -1.0])
    for x in boolean_vectors(2):
        assert float(x @ zero @ x) == float(x @ trace_zero_diagonal @ x) == 0.0
        checks += 1
    x = np.asarray([1.0, 1.0])
    y = np.asarray([1.0, -1.0])
    assert float(x @ zero @ y) == 0.0
    assert float(x @ trace_zero_diagonal @ y) == 2.0
    checks += 2

    # AR.34--AR.35: exact forward lifts allow arbitrarily slow innovation.
    current = np.zeros((1, 1), dtype=np.float64)
    partial_sum = 0.0
    for b in (0.25, 0.125, 0.0625):
        cell = b * np.asarray([[-1, 1], [1, -1]], dtype=np.float64)
        next_matrix = block_sum(current, cell)
        append = np.asarray([1.0, 1.0])
        assert float(append @ cell @ append) == 0.0
        for old_x in boolean_vectors(len(current)):
            lifted_x = np.concatenate((old_x, append))
            for old_y in boolean_vectors(len(current)):
                lifted_y = np.concatenate((old_y, append))
                assert float(lifted_x @ next_matrix @ lifted_y) == float(
                    old_x @ current @ old_y
                )
                checks += 1
        current = next_matrix
        partial_sum += b
        self_values = [float(z @ current @ z) for z in boolean_vectors(len(current))]
        assert min(self_values) == -4 * partial_sum
        assert max(self_values) == 0.0
        checks += 2

    # AR.30--AR.31: growing-size duplication has vanishing distortion while
    # its scalar carrier follows the prescribed nonconvergent coefficient.
    for r in range(6):
        n_r = 2**r
        n_next = 2 * n_r
        c_r = sin(log(r + 2))
        c_next = sin(log(r + 3))
        a_r = c_r * np.eye(n_r) / n_r  # normalized A_r/n_r^(3/2)
        a_next = c_next * np.eye(n_next) / n_next
        duplicate = np.kron(np.eye(n_r), np.ones((2, 1)))
        drift = duplicate.T @ a_next @ duplicate - a_r
        assert np.allclose(drift, ((c_next - c_r) / n_r) * np.eye(n_r))
        assert abs(np.sum(np.abs(drift)) - abs(c_next - c_r)) < 1e-12
        assert abs(c_next - c_r) <= log((r + 3) / (r + 2)) + 1e-12
        # The all-ones witness displays the scalar self-response c_r; the
        # identity for every Boolean vector follows algebraically from x^Tx=n.
        witness = np.ones(n_r)
        unnormalized = c_r * np.sqrt(n_r) * np.eye(n_r)
        assert abs(float(witness @ unnormalized @ witness) / n_r**1.5 - c_r) < 1e-12
        checks += 4

    # Numerically sample the two analytic subsequences used to prove that
    # sin(log(r+2)) has limits +1 and -1 along rounded exponential phases.
    for k in range(1, 5):
        positive_index = round(np.exp(np.pi / 2 + 2 * np.pi * k))
        negative_index = round(np.exp(3 * np.pi / 2 + 2 * np.pi * k))
        assert sin(log(positive_index)) > 1 - 1e-7
        assert sin(log(negative_index)) < -1 + 1e-7
        checks += 2

    # AR.29d--AR.29h: flipping a symmetric perfect matching is a genuinely
    # non-tensor full-sign perturbation of operator norm two.  Its normalized
    # maxima obey the quasi-monotone bound, and hollowing costs only trace/2.
    sign_matrices = [np.asarray([[1.0]])]
    matching_errors: list[np.ndarray] = []
    for _ in range(2):
        exact_tensor = np.kron(sign_matrices[-1], h)
        next_sign = exact_tensor.copy()
        for i in range(0, len(next_sign), 2):
            next_sign[i, i + 1] *= -1
            next_sign[i + 1, i] *= -1
        matching_error = next_sign - exact_tensor
        assert np.array_equal(next_sign, next_sign.T)
        assert set(next_sign.ravel()) == {-1.0, 1.0}
        assert abs(np.linalg.norm(matching_error, 2) - 2.0) < 1e-12
        sign_matrices.append(next_sign)
        matching_errors.append(matching_error)
        checks += 3

    sign_extrema = [quadratic_extrema(matrix) for matrix in sign_matrices]
    for r, matching_error in enumerate(matching_errors):
        n0 = len(sign_matrices[r])
        n1 = len(sign_matrices[r + 1])
        delta = np.linalg.norm(matching_error, 2) / np.sqrt(n1)
        for coordinate in (0, 1):
            p0 = sign_extrema[r][coordinate] / (2 * n0**1.5)
            p1 = sign_extrema[r + 1][coordinate] / (2 * n1**1.5)
            assert p1 + delta / 2 + 1e-12 >= p0
            checks += 1

    final_full = sign_matrices[-1]
    final_hollow = final_full.copy()
    np.fill_diagonal(final_hollow, 0)
    _, full_absolute = quadratic_extrema(final_full)
    _, hollow_absolute = quadratic_extrema(final_hollow)
    n_final = len(final_full)
    p_full = full_absolute / (2 * n_final**1.5)
    q_hollow = hollow_absolute / (2 * n_final**1.5)
    assert abs(q_hollow - p_full) <= 1 / (2 * np.sqrt(n_final)) + 1e-12
    for x in boolean_vectors(n_final):
        assert float(x @ final_full @ x - x @ final_hollow @ x) == float(
            np.trace(final_full)
        )
    checks += 1 + 2**n_final

    # Absolute response cannot orient the signed carrier.
    assert abs((-1.0) ** 4) == abs((-1.0) ** 5) == 1.0
    assert (-1.0) ** 4 != (-1.0) ** 5
    checks += 2

    print(f"approximate Boolean recovery checks passed: {checks}")


if __name__ == "__main__":
    verify()
