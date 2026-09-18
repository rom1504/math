#!/usr/bin/env python3
"""Exact independent physical-cube audit of the complex BH lower seed.

No imports of the other certificates. Work in Z[sqrt(-3)], evaluate every
physical vertex first, then invert by an unnormalized Walsh transform.
All assertions use integers or Fraction; no floating-point evidence.
"""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path


def mul(z, w):
    return z[0]*w[0]-3*z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conjugate(z):
    return z[0], -z[1]


def norm(z):
    return z[0]*z[0]+3*z[1]*z[1]


def subtract(z, w):
    return z[0]-w[0], z[1]-w[1]


def transform(values):
    result = list(values)
    stride = 1
    while stride < len(result):
        for start in range(0, len(result), 2*stride):
            for offset in range(stride):
                a, b = result[start+offset], result[start+stride+offset]
                result[start+offset] = a[0]+b[0], a[1]+b[1]
                result[start+stride+offset] = subtract(a, b)
        stride *= 2
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    # Twice the six powers of zeta=(1+sqrt(-3))/2.
    roots2 = [(2, 0), (1, 1), (-1, 1), (-2, 0), (-1, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    phases = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f8, a2, b128 = [], [], []
    for word in range(32):
        value = [0, 0]
        for mask, phase in zip(masks, phases):
            sign = (-1)**bin(mask & word).count("1")
            value[0] += sign*roots2[phase][0]
            value[1] += sign*roots2[phase][1]
        f = tuple(value)
        assert norm(f) == 64
        x3, x4 = 1-2*((word >> 2) & 1), 1-2*((word >> 3) & 1)
        a = (1+x3+2*x4, -1+x3)
        b = mul(mul(f, f), conjugate(a))
        f8.append(f)
        a2.append(a)
        b128.append(b)
    assert max(bin(i).count("1") for i, z in enumerate(transform(b128)) if z != (0, 0)) == 3

    F64, h256, C4194304 = [], [], []
    quotient_counts = Counter()
    for word in range(1024):
        x, y = word & 31, word >> 5
        F = mul(f8[x], f8[y])
        h = subtract(mul(a2[x], b128[y]), mul(b128[x], a2[y]))
        quotient = mul(h, conjugate(F))
        assert quotient[0] == 0
        assert quotient[1] in (0, -65536, 65536)
        quotient_counts[quotient[1] // 16384] += 1
        assert norm(h) in (0, 48*256**2)
        F64.append(F)
        h256.append(h)
        C4194304.append((F[0]*norm(h), F[1]*norm(h)))
    assert quotient_counts == {0: 640, -4: 192, 4: 192}
    Fhat, hhat, Chat = transform(F64), transform(h256), transform(C4194304)
    old = {i for i, z in enumerate(Fhat) if z != (0, 0)}
    support_h = {i for i, z in enumerate(hhat) if z != (0, 0)}
    outside = support_h-old
    assert len(old) == 256 and len(support_h) == 60 and len(outside) == 36
    assert all(Fraction(norm(Fhat[i]), (64*1024)**2) == Fraction(1, 256) for i in old)
    assert max(bin(i).count("1") for i in old | support_h) == 4
    outside_norms = Counter(Fraction(norm(hhat[i]), (256*1024)**2) for i in outside)
    assert outside_norms == {Fraction(3, 16): 24, Fraction(9, 16): 12}
    assert sum(outside_norms[x]*x for x in outside_norms) == Fraction(45, 4)
    assert sum(Fraction(norm(z), (256*1024)**2) for z in hhat) == 18
    assert max(bin(i).count("1") for i, z in enumerate(Chat) if z != (0, 0)) == 8
    assert Fraction(Chat[503][0], 4194304*1024) == Fraction(3, 4)
    assert Chat[503][1] == 0
    assert Fraction(Chat[503][0], 48*4194304*1024) == Fraction(1, 64)
    # FR=C/48 has the exact degree needed by the joint normalization.
    assert all(norm(h) == 0 or norm(h) == 48*256**2 for h in h256)

    epsilon = Fraction(1, 1024)
    denominator_cost = Fraction(4, 5)*48*epsilon**2
    rational_gain = (Fraction(45, 4)*Fraction(1, 65536)-4*denominator_cost)/(1+denominator_cost)
    assert rational_gain == Fraction(33, 16*81923) > Fraction(1, 65536)
    assert 2**8 < Fraction(7, 2)**5
    assert Fraction(7, 2)+Fraction(1, 65536) < 4
    assert 4**3 < 2**8
    result = {
        "status": "PASS", "evidence": "integer and rational arithmetic only",
        "algorithm": "physical cube evaluation then inverse Walsh transform",
        "vertices": 1024, "old_support": 256, "tangent_support": 60,
        "new_support": 36, "new_squared_mass": "45/4", "tangent_degree": 4,
        "normalization_term_degree": 8, "unnormalized_C_top_coefficient": "3/4",
        "FR_top_coefficient": "1/64",
        "quotient_over_sqrt_minus3_histogram": dict(quotient_counts),
        "rational_q_mass_gain": str(rational_gain),
        "finite_lower_bound": "B_4 > 2 + 2^(-18)",
        "scope": "finite certificate; all-order normalization also needs the analytic proof",
    }
    payload = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
