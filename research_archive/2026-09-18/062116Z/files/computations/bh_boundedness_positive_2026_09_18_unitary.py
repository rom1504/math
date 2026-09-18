"""Exact tangent-gap certificate for the five-bit unimodular seed.

No files are written.  All assertions use integer or Fraction arithmetic.
The matrix uses h=c+i*sqrt(3)*d so that its entries are rational; the
coefficient norm is then ||c||_2^2+3||d||_2^2.
"""

from fractions import Fraction
import json


ROOTS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
MASKS = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
POWERS = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
AFFINE = [0, 1, 2, 4, 8, 16]


def sign(mask, word):
    return (-1) ** bin(mask & word).count("1")


def value(word):
    a = sum(sign(mask, word) * ROOTS[power][0] for mask, power in zip(MASKS, POWERS))
    b = sum(sign(mask, word) * ROOTS[power][1] for mask, power in zip(MASKS, POWERS))
    assert a % 4 == b % 4 == 0
    out = a // 4, b // 4
    assert out in ROOTS
    return out


def positive_pivots(matrix):
    a = [list(row) for row in matrix]
    pivots = []
    for i in range(len(a)):
        pivot = a[i][i]
        assert pivot > 0
        pivots.append(pivot)
        for j in range(i + 1, len(a)):
            for k in range(i + 1, len(a)):
                a[j][k] -= a[j][i] * a[i][k] / pivot
    return pivots


def main():
    rows = []
    for word in range(32):
        a, b = value(word)
        characters = [sign(mask, word) for mask in AFFINE]
        # Re(f)=a+b/2; sqrt(3) Im(f)=3b/2.
        rows.append([Fraction(2 * a + b, 2) * z for z in characters]
                    + [Fraction(3 * b, 2) * z for z in characters])
    gram = [[sum(row[i] * row[j] for row in rows) / 32
             for j in range(12)] for i in range(12)]
    c = [[int(gram[i][j + 6] * Fraction(16, 3)) for j in range(6)] for i in range(6)]
    for i in range(6):
        assert sum(abs(v) for v in c[i]) == 2
        for j in range(6):
            assert c[i][j] == c[j][i]
            assert gram[i][j] == Fraction(7 if i == j else 0, 16)
            assert gram[i + 6][j + 6] == Fraction(27 if i == j else 0, 16)
            assert gram[i][j + 6] == Fraction(3 * c[i][j], 16)
    metric = [1] * 6 + [3] * 6
    shifted = [[gram[i][j] - (Fraction(metric[i], 4) if i == j else 0)
                for j in range(12)] for i in range(12)]
    pivots = positive_pivots(shifted)
    print(json.dumps({
        "status": "PASS", "base_variables": 5, "base_degree": 2,
        "base_fourier_support": 16,
        "affine_tangent_gap_squared_lower_bound": "1/4",
        "uniform_extra_variable_rigidity_radius": "1/4 (strict)",
        "cross_block_integer_matrix": c,
        "shifted_gram_exact_positive_pivots": [str(v) for v in pivots],
    }, indent=2))


if __name__ == "__main__":
    main()
