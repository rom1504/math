"""Exact finite audit of the range-localization composition.

All integer computations; no numerical Grothendieck or asymptotic claims.
Run from repository root with .venv/bin/python.
"""

from itertools import product
import json
from pathlib import Path
import random


def profile(matrix):
    n = len(matrix)
    if n == 0:
        return 0, 0, 0
    minimum = 10**9
    maximum = -10**9
    beta = 0
    for x in product((-1, 1), repeat=n):
        ax = [sum(matrix[i][j] * x[j] for j in range(n))
              for i in range(n)]
        doubled = sum(x[i] * ax[i] for i in range(n))
        assert doubled % 2 == 0
        energy = doubled // 2
        minimum = min(minimum, energy)
        maximum = max(maximum, energy)
        beta = max(beta, sum(abs(v) for v in ax))
    assert minimum <= 0 <= maximum
    assert beta <= 2 * (maximum - minimum)  # beta <= 4 half-range
    return minimum, maximum, beta


def make_matrix(n, signs):
    matrix = [[0] * n for _ in range(n)]
    for value, (i, j) in zip(signs,
                            ((i, j) for i in range(n)
                             for j in range(i + 1, n))):
        matrix[i][j] = matrix[j][i] = value
    return matrix


def audit_matrix(matrix):
    n = len(matrix)
    profiles = {}
    full = (1 << n) - 1
    for mask in range(1 << n):
        indices = [i for i in range(n) if mask >> i & 1]
        sub = [[matrix[i][j] for j in indices] for i in indices]
        profiles[mask] = profile(sub)
    low, high, _ = profiles[full]
    checks = 0
    for mask in range(1 << n):
        low_s, high_s, _ = profiles[mask]
        low_t, high_t, _ = profiles[full ^ mask]
        # These two one-sided statements follow separately by bridge reversal.
        assert high >= high_s + high_t
        assert low <= low_s + low_t
        assert high - low >= high_s - low_s + high_t - low_t
        assert 2 * max(high, -low) >= high_s - low_s + high_t - low_t
        checks += 1
    return checks, len(profiles)


def main():
    rng = random.Random(202609172103)
    matrices = partitions = profiles = 0
    by_order = {}
    for n in range(2, 6):
        count = 0
        for signs in product((-1, 1), repeat=n * (n - 1) // 2):
            p, s = audit_matrix(make_matrix(n, signs))
            matrices += 1
            partitions += p
            profiles += s
            count += 1
        by_order[n] = {"full_signings": count, "exhaustive": True}
    for n in range(6, 11):
        for _ in range(8):
            signs = [rng.choice((-1, 1)) for _ in range(n * (n - 1) // 2)]
            p, s = audit_matrix(make_matrix(n, signs))
            matrices += 1
            partitions += p
            profiles += s
        by_order[n] = {"full_signings": 8, "exhaustive": False}
    result = {
        "status": "PASS",
        "full_signings": matrices,
        "principal_profiles_beta_half_range_checks": profiles,
        "partition_range_superadditivity_checks": partitions,
        "by_order": by_order,
        "scope": "Exact finite normalization audit; no numerical proof of Grothendieck or asymptotic lower bound.",
    }
    output = Path("tmp/paper_portfolio_2026_09_17/bernoulli/range_localization_audit.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
