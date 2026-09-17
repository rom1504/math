"""Replay the scalar zonotope obstruction; no large vectors are materialized.

Run from the repository root with .venv/bin/python.
Analytic proof: artifacts/paper_bernoulli_2026_09_17.md, Section 8.
The numerical checks verify finite witnesses, not the universal theorem.
"""

import json
from fractions import Fraction

import mpmath as mp


def normal_hinge(a):
    """E (|G|-a)_+ for a>=0, evaluated at high precision."""
    a = mp.mpf(a)
    density = mp.exp(-a * a / 2) / mp.sqrt(2 * mp.pi)
    upper_tail = mp.erfc(a / mp.sqrt(2)) / 2
    return 2 * (density - a * upper_tail)


def verify_lattice(m):
    assert m >= 2 and m % 2 == 0
    tiny_count = m ** 3
    tiny_numerators = range(-tiny_count, tiny_count + 1, 2)
    actual = set()
    for first_sign in (-1, 1):
        actual.update(first_sign * m ** 2 + s for s in tiny_numerators)
    expected = set(range(-(m + 1) * m ** 2, (m + 1) * m ** 2 + 1, 2))
    assert actual == expected
    assert Fraction(tiny_count, m ** 4) == Fraction(1, m)
    return len(actual)


def finite_witness(c):
    c = mp.mpf(c)
    h = normal_hinge(2 * c)
    # Deliberately generous, explicit even integer; source dimension m^3+1.
    cutoff = max(4 * c + 2, 16 * (c + 1) / h)
    m = 2 * int(mp.ceil(cutoff / 2))
    assert m > 2 * c
    continuous_left = h - normal_hinge(m)
    quantized_left_lower = continuous_left - mp.mpf(2) / m ** 2
    reference_right_upper = c / (4 * m)
    assert quantized_left_lower > reference_right_upper
    return {
        "C": str(c),
        "m": m,
        "source_dimension": m ** 3 + 1,
        "AA_star": "1 + 1/m",
        "continuous_left": mp.nstr(continuous_left, 20),
        "quantized_left_lower": mp.nstr(quantized_left_lower, 20),
        "reference_right_upper": mp.nstr(reference_right_upper, 20),
        "ratio_lower": mp.nstr(quantized_left_lower / reference_right_upper, 20),
    }


def main():
    mp.mp.dps = 80
    lattice = {m: verify_lattice(m) for m in (2, 4, 6, 8, 10)}
    # Check the elementary hinge/variance inequality on exact rationals.
    rational_checks = 0
    for denominator in range(1, 21):
        for numerator in range(0, 20 * denominator + 1):
            s = Fraction(numerator, denominator)
            assert max(Fraction(0), s - 1) <= s * s / 4
            rational_checks += 1
    report = {
        "status": "PASS",
        "lattice_image_sizes": lattice,
        "exact_hinge_checks": rational_checks,
        "finite_witnesses": [finite_witness(c) for c in ("0.5", "1", "2", "3")],
        "scope": "Scalar zonotope proxy fails; not an ambient-subGaussian sign-lift counterexample.",
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
