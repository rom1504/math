"""Exact Gaussian-integer certificate for the sharp unitary bilinear seed.

P has eight coefficients of modulus one and |P|^2=8.  Consequently
P/sqrt(8) is unimodular, degree two, with BH ratio 8^(1/4).
No floating-point arithmetic or external packages are used.
"""

import json


def multiply(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def main():
    coefficients = {5: (1, 0), 6: (1, 0), 9: (0, 1), 10: (0, 1),
                    17: (1, 0), 18: (-1, 0), 33: (0, 1), 34: (0, -1)}
    square = {}
    for left_mask, left in coefficients.items():
        for right_mask, right in coefficients.items():
            value = multiply(left, (right[0], -right[1]))
            mask = left_mask ^ right_mask
            old = square.get(mask, (0, 0))
            square[mask] = (old[0] + value[0], old[1] + value[1])
    square = {mask: value for mask, value in square.items() if value != (0, 0)}
    assert square == {0: (8, 0)}
    for word in range(64):
        real = imag = 0
        for mask, value in coefficients.items():
            sign = (-1) ** bin(word & mask).count("1")
            real += sign * value[0]
            imag += sign * value[1]
        assert real * real + imag * imag == 8
    assert all(bin(mask).count("1") == 2 for mask in coefficients)
    print(json.dumps({"exact_certificate": "PASS", "variables": 6,
                      "degree": 2, "support": 8, "vertices_checked": 64,
                      "unnormalized_modulus_squared": 8,
                      "critical_ratio": "2^(3/4)"}))


if __name__ == "__main__":
    main()
