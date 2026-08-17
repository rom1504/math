#!/usr/bin/env python3
"""Independent checks for the Boolean--spherical regular-Hadamard gap.

Adds a direct SA.3/GE normalization check and exact completed-signing
Lipschitz checks at the order-16 base instance.
"""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np

from verify_bounded_cap_contextual_metric_compiler import build


V0 = np.asarray(
    (
        -1, -1, -1, 1,
        -1, -1, 1, -1,
        1, -1, 1, 1,
        -1, 1, 1, 1,
    ),
    dtype=np.int64,
)


def hollow_energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    return sum(
        int(matrix[i, j]) * int(spin[i]) * int(spin[j])
        for i in range(len(spin))
        for j in range(i + 1, len(spin))
    )


def auxiliary_table(completion: np.ndarray, m: int):
    rows = []
    for word in product((-1, 1), repeat=2 * m):
        spin = np.asarray(word, dtype=np.int64)
        rows.append(
            (
                int(np.sum(spin[:m])),
                int(np.sum(spin[m:])),
                hollow_energy(completion, spin),
            )
        )
    return rows


def completed_boolean(
    h: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    m: int,
    completion: np.ndarray,
) -> int:
    table = auxiliary_table(completion, m)
    best = 0
    trace = int(np.trace(h))
    for old_word in product((-1, 1), repeat=len(h)):
        old = np.asarray(old_word, dtype=np.int64)
        child = (int(old @ h @ old) - trace) // 2
        field_a = int(a @ old)
        field_b = int(b @ old)
        for sum_a, sum_b, auxiliary in table:
            best = max(best, abs(child + sum_a * field_a + sum_b * field_b + auxiliary))
    return best


def completed_spherical(
    n: int, r: int, m: int, completion: np.ndarray
) -> float:
    """Exact sphere optimum because every auxiliary field is in E_+(H)."""

    best = 0.0
    for sum_a, sum_b, auxiliary in auxiliary_table(completion, m):
        field_norm = sqrt(n * (sum_a * sum_a + sum_b * sum_b))

        # Positive outer channel: a convex function of ||u_+||, so one of
        # the two endpoints is optimal.
        positive = max(
            r * n / 2 + field_norm * sqrt(n) + auxiliary,
            -r * n / 2 + auxiliary,
        )

        # Negative outer channel: a concave quadratic in s=||u_+||.
        s = min(sqrt(n), field_norm / (2 * r))
        negative = -r * s * s + r * n / 2 + field_norm * s - auxiliary
        best = max(best, positive, negative)
    return best


def completion_cap(completion: np.ndarray) -> int:
    return max(
        abs(hollow_energy(completion, np.asarray(word, dtype=np.int64)))
        for word in product((-1, 1), repeat=len(completion))
    )


def check_sa_ge_normalization() -> int:
    checks = 0
    for mu in (0.0, 0.25, 0.5, 1.0, 1.7):
        # p=2, G=R=I.  In the positive channel g=h=2, hence GE has
        # a=g+h=4, b=0 and A=4 mu^2.
        spherical_normalized = 0.5 + sqrt((4 * mu * mu) / 2)
        boolean_normalized = 0.5 + mu
        assert abs(spherical_normalized - (0.5 + sqrt(2) * mu)) < 1e-12
        assert abs(
            spherical_normalized - boolean_normalized - (sqrt(2) - 1) * mu
        ) < 1e-12
        c = 2 * mu
        assert abs(
            spherical_normalized - boolean_normalized - c * (sqrt(2) - 1) / 2
        ) < 1e-12
        checks += 1
    return checks


def check_exact_completions() -> int:
    r, n, h_list, _ = build(2)
    h = np.asarray(h_list, dtype=np.int64)
    a = np.ones(n, dtype=np.int64)
    b = V0
    m = r // 2
    boolean_zero = r * n / 2 + m * n
    spherical_zero = r * n / 2 + sqrt(2) * m * n

    completions = []
    all_positive = np.ones((2 * m, 2 * m), dtype=np.int64) - np.eye(2 * m, dtype=np.int64)
    completions.append(all_positive)
    patterned = np.zeros((2 * m, 2 * m), dtype=np.int64)
    for i in range(2 * m):
        for j in range(i + 1, 2 * m):
            patterned[i, j] = patterned[j, i] = 1 if (i + 2 * j) % 3 else -1
    completions.append(patterned)

    checks = 0
    for completion in completions:
        q = completion_cap(completion)
        boolean = completed_boolean(h, a, b, m, completion)
        spherical = completed_spherical(n, r, m, completion)
        assert abs(boolean - boolean_zero) <= q
        assert abs(spherical - spherical_zero) <= q + 1e-10
        assert spherical - boolean >= spherical_zero - boolean_zero - 2 * q - 1e-10
        assert q <= (2 * m) * (2 * m - 1) // 2
        checks += 1
    return checks


def check_tensor_without_large_enumeration() -> int:
    r0, n0, h_list, _ = build(2)
    h0 = np.asarray(h_list, dtype=np.int64)
    h2 = np.kron(h0, h0)
    n, r = n0**2, r0**2
    a = np.ones(n, dtype=np.int64)
    b = np.kron(V0, np.ones(n0, dtype=np.int64))
    assert np.trace(h2) == 0
    assert np.array_equal(h2 @ h2, r * r * np.eye(n, dtype=np.int64))
    assert np.array_equal(h2 @ a, r * a)
    assert np.array_equal(h2 @ b, r * b)
    assert int(a @ b) == 0
    m = r // 2
    gap = (sqrt(2) - 1) * m * n
    assert abs(gap / n**1.5 - (sqrt(2) - 1) / 2) < 1e-12
    return 6


def main() -> None:
    normalization = check_sa_ge_normalization()
    completions = check_exact_completions()
    tensor = check_tensor_without_large_enumeration()
    print(
        "regular-Hadamard Boolean--spherical independent audit: PASS",
        f"normalization={normalization}",
        f"completions={completions}",
        f"tensor={tensor}",
    )


if __name__ == "__main__":
    main()
