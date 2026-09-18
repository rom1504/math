"""Exact sharpness certificates for unimodular block-affine polynomials.

For k=1,...,5, construct two disjoint Boolean address polynomials G,H
with denominator N=2^(k-1).  The Gaussian-integer numerator G+iH has
constant squared modulus 2*N^2 and exactly 2*N^2 Fourier coefficients.
After division by sqrt(2)*N this attains ratio 2^(1-1/(2k)).
"""

import json


def product(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def certificate(k):
    controls = k - 1
    addresses = 2 ** controls
    half_variables = controls + addresses
    coefficients = {}
    block_of_bit = {}
    for copy in range(2):
        offset = copy * half_variables
        for bit in range(controls):
            block_of_bit[offset + bit] = bit
        for address in range(addresses):
            block_of_bit[offset + controls + address] = controls
            for subset in range(addresses):
                mask = (subset << offset) | (1 << (offset + controls + address))
                sign = (-1) ** bin(subset & address).count("1")
                coefficients[mask] = (sign, 0) if copy == 0 else (0, sign)
    assert len(coefficients) == 2 * addresses ** 2
    for mask in coefficients:
        blocks = [block_of_bit[bit] for bit in block_of_bit if mask & (1 << bit)]
        assert len(blocks) == len(set(blocks))
    assert max(bin(mask).count("1") for mask in coefficients) == k
    square = {}
    for left_mask, left in coefficients.items():
        for right_mask, right in coefficients.items():
            value = product(left, (right[0], -right[1]))
            mask = left_mask ^ right_mask
            old = square.get(mask, (0, 0))
            square[mask] = (old[0] + value[0], old[1] + value[1])
    square = {mask: value for mask, value in square.items() if value != (0, 0)}
    assert square == {0: (2 * addresses ** 2, 0)}
    return {"exact_certificate": "PASS", "blocks": k,
            "variables": 2 * half_variables, "degree": k,
            "support": len(coefficients),
            "unnormalized_modulus_squared": 2 * addresses ** 2,
            "critical_ratio": "2^(1-1/(2*%d))" % k}


if __name__ == "__main__":
    for block_count in range(1, 6):
        print(json.dumps(certificate(block_count)), flush=True)
