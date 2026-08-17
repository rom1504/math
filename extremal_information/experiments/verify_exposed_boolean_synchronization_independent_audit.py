#!/usr/bin/env python3
"""Independent checks for exposed Boolean synchronization.

Adds exact correlated-port/GE normalization checks, the rounding theorem on
the exposed Walsh optimizer, and completed-signing checks at order 16.
"""

from __future__ import annotations

from itertools import product
from math import gamma, pi, sqrt

import numpy as np

from verify_exposed_boolean_synchronization import construction


def hollow_energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    return sum(
        int(matrix[i, j]) * int(spin[i]) * int(spin[j])
        for i in range(len(spin))
        for j in range(i + 1, len(spin))
    )


def auxiliary_table(completion: np.ndarray, width: int):
    table = []
    for word in product((-1, 1), repeat=2 * width):
        spin = np.asarray(word, dtype=np.int64)
        table.append(
            (
                int(np.sum(spin[:width])),
                int(np.sum(spin[width:])),
                hollow_energy(completion, spin),
            )
        )
    return table


def completed_boolean(
    h: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    width: int,
    completion: np.ndarray,
) -> int:
    table = auxiliary_table(completion, width)
    trace = int(np.trace(h))
    best = 0
    for word in product((-1, 1), repeat=len(h)):
        old = np.asarray(word, dtype=np.int64)
        child = (int(old @ h @ old) - trace) // 2
        field_a, field_b = int(a @ old), int(b @ old)
        for sum_a, sum_b, auxiliary in table:
            best = max(best, abs(child + sum_a * field_a + sum_b * field_b + auxiliary))
    return best


def completed_spherical(
    n: int,
    r: int,
    rho: float,
    width: int,
    completion: np.ndarray,
) -> float:
    best = 0.0
    for sum_a, sum_b, auxiliary in auxiliary_table(completion, width):
        squared_field = n * (
            sum_a * sum_a + sum_b * sum_b + 2 * rho * sum_a * sum_b
        )
        field_norm = sqrt(max(0.0, squared_field))
        positive = max(
            r * n / 2 + field_norm * sqrt(n) + auxiliary,
            -r * n / 2 + auxiliary,
        )
        positive_radius = min(sqrt(n), field_norm / (2 * r))
        negative = (
            -r * positive_radius**2
            + r * n / 2
            + field_norm * positive_radius
            - auxiliary
        )
        best = max(best, positive, negative)
    return best


def completion_cap(completion: np.ndarray) -> int:
    return max(
        abs(hollow_energy(completion, np.asarray(word, dtype=np.int64)))
        for word in product((-1, 1), repeat=len(completion))
    )


def check_correlated_formula_and_ge() -> int:
    q, _, h, _, _, _, a, b = construction(2)
    n, r, width = len(a), q, q // 2
    rho = float(a @ b) / n
    assert abs(rho - 0.5) < 1e-12

    boolean = r * n / 2 + width * n * (1 + rho)
    spherical = r * n / 2 + width * n * sqrt(2 * (1 + rho))

    # GE positive channel: g=h=2(1+rho), so a_GE=4(1+rho), b_GE=0.
    mu = width / r
    ge_normalized = 0.5 + sqrt(mu * mu * 4 * (1 + rho) / 2)
    assert abs(ge_normalized - spherical / (r * n)) < 1e-12

    # Exact finite Boolean enumeration, independent of the closed formula.
    trace = int(np.trace(h))
    enumerated = 0
    for word in product((-1, 1), repeat=n):
        old = np.asarray(word, dtype=np.int64)
        child = abs((int(old @ h @ old) - trace) // 2)
        field = width * (abs(int(a @ old)) + abs(int(b @ old)))
        enumerated = max(enumerated, child + field)
    assert enumerated == boolean
    return 4


def check_rounding_bound_on_exposed_optimizer() -> int:
    q, _, h, _, _, _, a, b = construction(2)
    n, r, width = len(a), q, q // 2
    rho = float(a @ b) / n
    u = (a + b) / sqrt(2 * (1 + rho))
    x = np.where(u >= 0, 1.0, -1.0)
    phi = 1 - float(np.sum(np.abs(u))) / n
    assert abs(float(np.sum((x - u) ** 2)) / n - 2 * phi) < 1e-12

    sphere = float(u @ h @ u) / 2 + width * float((a + b) @ u)
    cube_same_channel = float(x @ h @ x) / 2 + width * float((a + b) @ x)
    c = 2 * width / r
    local_loss = sphere - cube_same_channel
    assert local_loss <= r * n * (1 + c) * sqrt(2 * phi) + 1e-9

    exact_boolean = r * n / 2 + width * n * (1 + rho)
    exact_spherical = r * n / 2 + width * n * sqrt(2 * (1 + rho))
    assert exact_spherical - exact_boolean <= r * n * (1 + c) * sqrt(2 * phi)
    return 3


def check_completed_family() -> int:
    q, _, h, _, _, _, a, b = construction(2)
    n, r, width = len(a), q, q // 2
    rho = float(a @ b) / n
    boolean_zero = r * n / 2 + width * n * (1 + rho)
    spherical_zero = r * n / 2 + width * n * sqrt(2 * (1 + rho))

    completions = []
    positive = np.ones((2 * width, 2 * width), dtype=np.int64)
    np.fill_diagonal(positive, 0)
    completions.append(positive)
    patterned = np.zeros_like(positive)
    for i in range(len(patterned)):
        for j in range(i + 1, len(patterned)):
            patterned[i, j] = patterned[j, i] = -1 if (i + j) % 3 == 0 else 1
    completions.append(patterned)

    checks = 0
    for completion in completions:
        cap = completion_cap(completion)
        boolean = completed_boolean(h, a, b, width, completion)
        spherical = completed_spherical(n, r, rho, width, completion)
        assert abs(boolean - boolean_zero) <= cap
        assert abs(spherical - spherical_zero) <= cap + 1e-10
        assert 0 <= spherical - boolean
        assert spherical - boolean <= spherical_zero - boolean_zero + 2 * cap + 1e-10
        checks += 1
    return checks


def check_gamma_constants() -> int:
    gamma_2 = sqrt(2) * gamma(1) / (sqrt(pi) * gamma(1.5))
    assert abs(gamma_2 - 2 * sqrt(2) / pi) < 1e-12
    assert gamma_2 < 1
    distance = sqrt(2 * (1 - gamma_2))
    assert distance > 0
    return 3


def check_common_pole_corollary() -> int:
    q, _, _, _, _, _, a, b = construction(2)
    n, r, width, p = len(a), q, q // 2, 2
    ports = (a, b)
    deficit = 1 - sum(abs(int(port @ a)) for port in ports) / (p * n)
    rho = float(a @ b) / n
    assert abs(deficit - 1 / q) < 1e-12
    c = width * p / r
    spherical = r * n / 2 + width * n * sqrt(2 * (1 + rho))
    boolean = r * n / 2 + width * n * (1 + rho)
    assert spherical <= r * n / 2 + width * p * n + 1e-12
    assert boolean >= r * n / 2 + width * p * n * (1 - deficit) - 1e-12
    assert spherical - boolean <= c * deficit * r * n + 1e-12
    return 4


def main() -> None:
    correlated = check_correlated_formula_and_ge()
    rounding = check_rounding_bound_on_exposed_optimizer()
    completion = check_completed_family()
    gamma_checks = check_gamma_constants()
    common_pole = check_common_pole_corollary()
    print(
        "exposed Boolean synchronization independent audit: PASS",
        f"correlated={correlated}",
        f"rounding={rounding}",
        f"completion={completion}",
        f"gamma={gamma_checks}",
        f"common_pole={common_pole}",
    )


if __name__ == "__main__":
    main()
