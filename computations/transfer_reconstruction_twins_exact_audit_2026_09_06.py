"""Independent integer/rational regression audit of the twin-seed obstruction.

Enumerates all sign diagonals in the stated finite cases. Precision and
determinant inequalities are exponentiated and checked over Fractions, with
no floating eigenvalues, logarithms, or Cholesky factors.
"""

from fractions import Fraction as F
from itertools import combinations, product
import json
import random


def determinant(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    sign = previous = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def product_entries(values):
    result = 1
    for value in values:
        result *= value
    return result


def quadratic(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(i + 1, len(a)))


def symmetric_signing(rng, n):
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = rng.choice((-1, 1))
    return a


def image(a, x):
    return [sum(entry * value for entry, value in zip(row, x)) for row in a]


def check_row_expectation(a):
    n = len(a)
    p = F(3, 8)
    images = [list(map(abs, image(a, (1,) + suffix))) for suffix in product((-1, 1), repeat=n - 1)]
    bilinear = max(map(sum, images))
    row_sum_by_size = [0] * (n + 1)
    for size in range(n + 1):
        for selected in combinations(range(n), size):
            row_sum_by_size[size] += max(sum(values[i] for i in selected) for values in images)
    expected_row_norm = sum(
        total * p**size * (1 - p) ** (n - size)
        for size, total in enumerate(row_sum_by_size)
    )
    signed_norm_by_size = [0] * (n + 1)
    for word in product((-1, 0, 1), repeat=n):
        size = sum(value != 0 for value in word)
        signed_norm_by_size[size] += sum(map(abs, image(a, word)))
    expected_signed_norm = sum(
        total * (p / 2) ** size * (1 - p) ** (n - size)
        for size, total in enumerate(signed_norm_by_size)
    )
    centered = expected_row_norm - p * bilinear
    assert 0 <= centered <= 2 * expected_signed_norm
    assert expected_signed_norm**2 <= n * n * p * n
    assert centered**2 <= 4 * n * n * p * n
    return {
        "order": n,
        "row_subsets": 2**n,
        "signed_bernoulli_words": 3**n,
        "bernoulli_mean": str(p),
        "expected_row_norm": str(expected_row_norm),
        "expected_signed_norm": str(expected_signed_norm),
        "full_bilinear_norm": bilinear,
    }


def pivots(matrix):
    a = [[F(entry) for entry in row] for row in matrix]
    out = []
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot > 0
        out.append(pivot)
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                value = a[i][j] - a[i][k] * a[k][j] / pivot
                a[i][j] = a[j][i] = value
    return out


def check_case(rng, n, count):
    a = symmetric_signing(rng, n)
    spins = [(1,) + suffix for suffix in product((-1, 1), repeat=n - 1)]
    abs_images = [list(map(abs, image(a, x))) for x in spins]
    row_norm, selected = min(
        (max(sum(values[i] for i in rows) for values in abs_images), rows)
        for rows in combinations(range(n), 2 * count)
    )
    pairs = list(zip(selected[::2], selected[1::2]))
    edge_signs = [rng.choice((-1, 1)) for _ in pairs]
    representative = list(range(n))
    for r, s in pairs:
        representative[s] = r
    b = [[a[representative[i]][representative[j]] for j in range(n)] for i in range(n)]
    for (r, s), sign in zip(pairs, edge_signs):
        b[r][s] = b[s][r] = sign
    assert all(b[i][i] == 0 for i in range(n))
    assert all(abs(b[i][j]) == 1 for i in range(n) for j in range(n) if i != j)
    cap_a = cap_b = 0
    for x in spins:
        y = [0] * n
        for i, value in enumerate(x):
            y[representative[i]] += value
        z = [yy - xx for yy, xx in zip(y, x)]
        pa, pb, py = quadratic(a, x), quadratic(b, x), quadratic(a, y)
        assert pb == py + sum(sign * x[r] * x[s] for (r, s), sign in zip(pairs, edge_signs))
        assert 2 * abs(py - pa) <= 3 * row_norm
        assert abs(sum(z[i] * value for i, value in enumerate(image(a, x)))) <= row_norm
        assert abs(sum(z[i] * value for i, value in enumerate(image(a, z)))) <= row_norm
        cap_a, cap_b = max(cap_a, abs(pa)), max(cap_b, abs(pb))
    assert 2 * cap_b <= 2 * cap_a + 3 * row_norm + 2 * count
    gram_det_bound = F(4, n) ** count * F(n, n - count) ** (n - count)
    completions = nonsingular = precision_checks = 0
    largest_determinant_ratio = F(0)
    for diagonal in product((-1, 1), repeat=n):
        full = [row[:] for row in b]
        for i, value in enumerate(diagonal):
            full[i][i] = value
        for r, s in pairs:
            difference = [full[i][r] - full[i][s] for i in range(n)]
            assert all(value == 0 for i, value in enumerate(difference) if i not in (r, s))
            assert sum(value * value for value in difference) <= 8
        det_b = determinant(full)
        det_g = F(det_b * det_b, n**n)
        assert det_g <= gram_det_bound
        largest_determinant_ratio = max(largest_determinant_ratio, det_g / gram_det_bound)
        completions += 1
        if not det_b:
            continue
        nonsingular += 1
        if nonsingular > 8:
            continue
        t = F(4)
        precisions = [t * F(i + 1, n + 1) for i in range(n)]
        for order in (list(range(n)), list(reversed(range(n))), rng.sample(range(n), n)):
            matrix = [
                [sum(precisions[k] * full[k][i] * full[k][j] for k in range(n)) / n for j in order]
                for i in order
            ]
            delta = pivots(matrix)
            assert all(0 < value <= t for value in delta)
            assert product_entries(delta) == det_g * product_entries(precisions)
            prefactor_fourth_power = product_entries(
                value / t * (2 - value / t) for value in precisions
            ) / product_entries(value / t * (2 - value / t) for value in delta)
            assert prefactor_fourth_power >= 1 / (det_g * 2**n)
            precision_checks += 1
    return {
        "order": n,
        "twin_pairs": count,
        "all_spin_checks_up_to_global_reversal": len(spins),
        "cap_before": cap_a,
        "cap_after": cap_b,
        "chosen_row_norm": row_norm,
        "all_sign_diagonal_completions": completions,
        "nonsingular_completions": nonsingular,
        "rational_precision_order_checks": precision_checks,
        "largest_determinant_to_bound_ratio": str(largest_determinant_ratio),
    }


def main():
    rng = random.Random(2026090602)
    results = {
        "row_restriction_expectation": check_row_expectation(symmetric_signing(rng, 8)),
        "twin_cases": [check_case(rng, n, m) for n, m in ((6, 1), (8, 2), (10, 3))],
        "arithmetic": "integer and rational only",
        "scope": "finite regression checks; asymptotic theorem requires its written proof",
    }
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
