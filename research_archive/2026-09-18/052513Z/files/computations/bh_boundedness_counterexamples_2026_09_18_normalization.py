#!/usr/bin/env python3
"""Replay degree-efficient normalization.

The finite rational assertions audit the constants in Section 8. The
500-digit entropy values are numerical diagnostics, not theorem inputs.
"""
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import mpmath as mp


def main():
    source = Path(__file__).with_name("bh_boundedness_counterexamples_2026_09_18_tangent.py")
    spec = importlib.util.spec_from_file_location("tangent", str(source))
    t = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t)
    assert Fraction(17, 4) ** 2 > 18
    assert 30 ** 2 > 864
    assert 24 * Fraction(17, 4) / 1024 + 15 < 16
    assert Fraction(17, 4) + Fraction(16, 1024) < Fraction(9, 2)
    assert Fraction(1, 8) + Fraction(9, 2048) < Fraction(17, 128)
    assert Fraction(7 * 81, 4) + 256 * Fraction(153, 256) ** 2 < 256
    assert (Fraction(1, 16) - Fraction(9, 2048)) ** 2 > Fraction(1, 512)
    assert 1 - Fraction(128, 1024) > Fraction(1, 2)
    assert 4 + Fraction(384, 1024 ** 2) < 5
    assert 402 * Fraction(2, 3) - 256 == 12
    assert Fraction(12, 10) > 1
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    a = {0: (1, -1), 4: (0, 1), 8: (1, 0)}
    raw = t.convolution(t.convolution(f4, f4), {i: t.conj(z) for i, z in a.items()})
    b4 = {i: (z[0] // 4, z[1] // 4) for i, z in raw.items()}
    F16 = {i | (j << 5): t.mul(v, w) for i, v in f4.items() for j, w in f4.items()}
    h4 = {}
    for i, v in a.items():
        for j, w in b4.items():
            key = i | (j << 5)
            h4[key] = t.add(h4.get(key, (0, 0)), t.mul(v, w))
    for i, v in b4.items():
        for j, w in a.items():
            key = i | (j << 5)
            z = t.mul(v, w)
            h4[key] = t.add(h4.get(key, (0, 0)), (-z[0], -z[1]))
    h4 = {i: z for i, z in h4.items() if z != (0, 0)}
    C256 = t.convolution(t.convolution(h4, h4), {i: t.conj(z) for i, z in F16.items()})
    C256 = {i: (-z[0], -z[1]) for i, z in C256.items()}
    RR = t.convolution(h4, {i: t.conj(z) for i, z in h4.items()})
    assert max(bin(i).count("1") for i in RR) == 8
    mp.mp.dps = 500
    zeta = mp.mpc(mp.mpf(".5"), mp.sqrt(3) / 2)

    def value(pair):
        return pair[0] + pair[1] * zeta

    reports = []
    for exponent in [10, 20, 64, 104]:
        epsilon = mp.mpf(2) ** -exponent
        u = 1 / mp.sqrt(1 + 48 * epsilon ** 2)
        probabilities = []
        for i in range(1024):
            coefficient = (value(F16.get(i, (0, 0))) / 16
                           + u * epsilon * value(h4.get(i, (0, 0))) / 4
                           + (u - 1) * value(C256.get(i, (0, 0))) / (256 * 48))
            probabilities.append(abs(coefficient) ** 2)
        assert abs(sum(probabilities) - 1) < mp.mpf("1e-490")
        entropy = -sum(p * mp.log(p) for p in probabilities if p)
        lower = 8 * mp.log(2) + epsilon ** 2 * (mp.mpf(45) / 4 * mp.log(1 / epsilon) - 256)
        assert entropy >= lower
        lam = mp.log(1 + 48 * epsilon ** 2) / 2
        d = 4 + 16 * lam
        scaled_gap = (entropy - 2 * d * mp.log(2)) / epsilon ** 2
        if exponent == 104:
            assert scaled_gap > 12
            assert mp.exp(entropy / (2 * d)) > 2 + mp.mpf(2) ** -207
        reports.append({"epsilon_exponent": exponent,
                        "entropy_minus_degree_budget_over_epsilon_squared_diagnostic": mp.nstr(scaled_gap, 18)})
    print(json.dumps({"status": "PASS", "rational_constant_checks": 11,
                      "R_degree_exact": 8, "entropy_diagnostics": reports,
                      "theorem": "liminf B_m > 2 + 2^(-207)"}, indent=2))


if __name__ == "__main__":
    main()
