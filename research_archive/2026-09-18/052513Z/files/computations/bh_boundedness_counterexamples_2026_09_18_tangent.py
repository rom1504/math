#!/usr/bin/env python3
"""Exact Eisenstein certificate for a complex Boolean BH ratio >2.

All theorem assertions below use integer/Fraction arithmetic. The final
decimal norm is explicitly diagnostic.
"""
from collections import Counter
from fractions import Fraction
import json
import math


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0] + z[1] * w[1]


def conj(z):
    return z[0] + z[1], -z[1]


def norm2(z):
    return z[0] ** 2 + z[0] * z[1] + z[1] ** 2


def convolution(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i ^ j] = add(out.get(i ^ j, (0, 0)), mul(a, b))
    return {i: a for i, a in out.items() if a != (0, 0)}


def evaluate(p, x):
    out = (0, 0)
    for i, z in p.items():
        sign = (-1) ** bin(x & i).count("1")
        out = add(out, (sign * z[0], sign * z[1]))
    return out


def main():
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    a = {0: (1, -1), 4: (0, 1), 8: (1, 0)}
    b16 = convolution(convolution(f4, f4), {i: conj(z) for i, z in a.items()})
    assert all(z[0] % 4 == z[1] % 4 == 0 for z in b16.values())
    b4 = {i: (z[0] // 4, z[1] // 4) for i, z in b16.items()}
    assert len(b4) == 12
    assert max(bin(i).count("1") for i in b4) == 3
    expected = {0: (1, -1), 3: (-1, -1), 4: (-1, 0), 7: (1, -2),
                8: (0, -1), 11: (0, -3), 17: (-1, -1), 18: (-3, 3),
                21: (3, 0), 22: (1, -2), 25: (2, -1), 26: (-2, 1)}
    assert b4 == expected
    F16 = {i | (j << 5): mul(v, w) for i, v in f4.items() for j, w in f4.items()}
    h4 = {}
    for i, v in a.items():
        for j, w in b4.items():
            key = i | (j << 5)
            h4[key] = add(h4.get(key, (0, 0)), mul(v, w))
    for i, v in b4.items():
        for j, w in a.items():
            key = i | (j << 5)
            z = mul(v, w)
            h4[key] = add(h4.get(key, (0, 0)), (-z[0], -z[1]))
    h4 = {i: z for i, z in h4.items() if z != (0, 0)}
    assert len(F16) == 256 and len(h4) == 60
    assert max(bin(i).count("1") for i in h4) == 4
    outside = {i: z for i, z in h4.items() if i not in F16}
    assert len(outside) == 36
    assert Counter(norm2(z) for z in outside.values()) == {3: 24, 9: 12}
    assert Fraction(sum(norm2(z) for z in outside.values()), 16) == Fraction(45, 4)
    assert h4[21 << 5] == (3, -3)
    gradient2 = sum(2 * mul(conj(F16.get(i, (0, 0))), z)[0]
                    + mul(conj(F16.get(i, (0, 0))), z)[1]
                    for i, z in h4.items())
    assert gradient2 == 0
    histogram = Counter()
    for x in range(32):
        fraw = evaluate(f4, x)
        assert fraw[0] % 4 == fraw[1] % 4 == 0
        fx = (fraw[0] // 4, fraw[1] // 4)
        assert fx in roots
        for y in range(32):
            word = x | (y << 5)
            Fraw, hraw = evaluate(F16, word), evaluate(h4, word)
            assert Fraw[0] % 16 == Fraw[1] % 16 == 0
            Fvalue = (Fraw[0] // 16, Fraw[1] // 16)
            assert hraw[0] % 4 == hraw[1] % 4 == 0
            hvalue = (hraw[0] // 4, hraw[1] // 4)
            quotient = mul(hvalue, conj(Fvalue))
            assert quotient in [(0, 0), (-4, 8), (4, -8)]
            histogram[quotient] += 1
    assert histogram == {(0, 0): 640, (-4, 8): 192, (4, -8): 192}
    # Exact radial normalization exits degree four: C=F |h/F|^2
    # =-h^2 conjugate(F), with coefficient denominator 256.
    C256 = convolution(convolution(h4, h4),
                       {i: conj(z) for i, z in F16.items()})
    C256 = {i: (-z[0], -z[1]) for i, z in C256.items()}
    assert len(C256) == 751
    assert max(bin(i).count("1") for i in C256) == 8
    assert bin(503).count("1") == 8
    assert C256[503] == (192, 0)
    epsilon = Fraction(1, 1024)
    denom_increase = Fraction(4, 5) * 48 * epsilon ** 2
    assert denom_increase == Fraction(3, 81920)
    gain = (Fraction(45, 262144) - 4 * denom_increase) / (1 + denom_increase)
    assert gain == Fraction(33, 16 * 81923) > Fraction(1, 65536)
    # A=2^(8/5)<7/2 is certified by fifth powers.
    assert 2 ** 8 < Fraction(7, 2) ** 5
    assert Fraction(7, 2) + Fraction(1, 65536) < 4
    # On(0,4], derivative( t^(5/8)) >1/4, since4^(3/8)<2.
    assert 4 ** 3 < 2 ** 8

    def complex_value(z):
        return complex(z[0] + z[1] / 2.0, math.sqrt(3) * z[1] / 2.0)

    q = 8 / 5.0
    total = sum(abs(complex_value(F16.get(i, (0, 0))) / 16
                    + float(epsilon) * complex_value(h4.get(i, (0, 0))) / 4) ** q
                for i in set(F16) | set(h4))
    diagnostic = total ** (1 / q) / math.sqrt(1 + 48 * float(epsilon) ** 2)
    print(json.dumps({"status": "PASS", "exact_lower_bound": "2 + 2^(-18)",
                      "degree": 4, "variables": 10, "epsilon": str(epsilon),
                      "old_support": len(F16), "new_support": len(outside),
                      "outside_squared_mass": "45/4", "rational_q_mass_gain": str(gain),
                      "radial_normalization_degree": 8,
                      "radial_normalization_top_coefficient": {"503": "3/4"},
                      "quotient_histogram": {str(k): v for k, v in histogram.items()},
                      "diagnostic_ratio_only": diagnostic}, indent=2))


if __name__ == "__main__":
    main()
