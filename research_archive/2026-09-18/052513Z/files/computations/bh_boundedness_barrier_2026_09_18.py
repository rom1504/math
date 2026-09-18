#!/usr/bin/env python3
"""Finite checks for scoped Boolean BH weighted-functional barriers.

Exact Fraction checks test symmetric Walsh recurrences and cube enumeration.
Floating-point CLT/contraction values are diagnostics, not asymptotic proofs.
The proofs are in artifacts/bh_boundedness_barrier_2026_09_18.md.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln


def symmetric_chebyshev(degree: int, dimension: int) -> list[Fraction]:
    """Common coefficient on each Walsh level of T_degree(sum(x)/dimension)."""
    old = [Fraction(0) for _ in range(degree + 2)]
    old[0] = Fraction(1)
    if degree == 0:
        return old[:1]
    cur = [Fraction(0) for _ in old]
    cur[1] = Fraction(1, dimension)
    for _ in range(1, degree):
        new = [Fraction(0) for _ in old]
        for r in range(min(degree, dimension) + 1):
            down = r * cur[r - 1] if r else Fraction(0)
            up = (dimension - r) * cur[r + 1]
            new[r] = 2 * (down + up) / dimension - old[r]
        old, cur = cur, new
    return cur[:degree + 1]


def chebyshev_fraction(k: int, z: Fraction) -> Fraction:
    if k == 0:
        return Fraction(1)
    prev, cur = Fraction(1), z
    for _ in range(1, k):
        prev, cur = cur, 2 * z * cur - prev
    return cur


def exact_checks() -> dict:
    assertions = 0
    enumerations = []
    for n in range(1, 10):
        for k in range(1, min(n, 7) + 1):
            coeff = symmetric_chebyshev(k, n)
            direct = [Fraction(0) for _ in range(k + 1)]
            for x in itertools.product((-1, 1), repeat=n):
                value = chebyshev_fraction(k, Fraction(sum(x), n))
                assert abs(value) <= 1
                assertions += 1
                for r in range(k + 1):
                    direct[r] += value * math.prod(x[:r]) / 2**n
            assert coeff == direct
            assertions += 1
            assert coeff[k] == Fraction(2**(k - 1) * math.factorial(k), n**k)
            assertions += 1
            assert all(coeff[r] == 0 for r in range(k + 1) if (r-k) % 2)
            assertions += 1
            enumerations.append([n, k])
    small_clt = []
    for k in (1, 3, 5, 7, 9, 15, 21, 31):
        coeff = symmetric_chebyshev(k, k*k)
        value = abs(k * coeff[1])
        q = 2*k/(k+1)
        full_norm = sum(math.comb(k*k, r)*abs(float(a))**q
                        for r, a in enumerate(coeff))**(1/q)
        small_clt.append({"k": k, "sqrt_n_abs_a1": float(value),
                          "full_q_k_norm_diagnostic": full_norm,
                          "full_q_k_norm_limit": math.sqrt((1-math.exp(-2))/2),
                          "exact_numerator": str(value.numerator),
                          "exact_denominator": str(value.denominator)})
        assert value <= 1
        assertions += 1
    return {"assertions": assertions, "enumeration_cases": enumerations,
            "exact_chebyshev_clt_diagnostics": small_clt}


def first_level_diagnostic(k: int, lam: float = 1.0) -> dict:
    n = round(lam * k*k)
    center = n / 2
    radius = math.ceil(8 * math.sqrt(n))
    minus = np.arange(max(0, math.floor(center-radius)),
                      min(n, math.ceil(center+radius)) + 1, dtype=np.float64)
    log_probability = (gammaln(n+1) - gammaln(minus+1)
                       - gammaln(n-minus+1) - n*math.log(2))
    probability = np.exp(log_probability)
    z = (n - 2*minus) / math.sqrt(n)
    h = np.sin(k * np.arcsin(z / math.sqrt(n)))
    value = float(np.sum(probability*z*h))
    return {"k": k, "n": n, "lambda": lam,
            "retained_probability": float(np.sum(probability)),
            "sqrt_n_abs_a1": abs(value),
            "limit": lam**(-0.5)*math.exp(-0.5/lam)}


def damping(d: int, e: int) -> float:
    return math.exp(d*e/(2*(d+e)) *
                    (math.log1p(2/(d-1)) + math.log1p(2/(e-1))))


def contraction_diagnostics() -> dict:
    central = []
    for d in (2, 3, 5, 10, 100, 1000, 100000):
        central.append({"d": d, "D_d_d": damping(d, d),
                        "s0_factor": damping(d, d),
                        "s1_factor": damping(d, d)/2,
                        "s_critical_factor": damping(d, d)/math.e})
    minimum_margin = math.inf
    for d in range(2, 251):
        for e in range(2, 251):
            theta = d/(d+e)
            entropy = -theta*math.log(theta)-(1-theta)*math.log1p(-theta)
            factor = damping(d, e)*math.exp(-entropy/math.log(2))
            minimum_margin = min(minimum_margin, factor-1)
            assert factor > 1
    return {"critical_s": 1/math.log(2), "central": central,
            "grid_minimum_critical_factor_minus_1": minimum_margin,
            "multi_color_thresholds": [
                {"colors": j, "threshold": (j-1)/math.log(j)}
                for j in range(2, 11)]}


def interval_scale(lo: Fraction, hi: Fraction, scale: Fraction):
    values = (scale*lo, scale*hi)
    return min(values), max(values)


def distance_to_zero(lo: Fraction, hi: Fraction) -> Fraction:
    return max(Fraction(0), lo, -hi)


def atan_interval(reciprocal: int, terms: int = 24):
    value = sum((Fraction((-1)**j, (2*j+1)*reciprocal**(2*j+1))
                 for j in range(terms)), Fraction(0))
    next_value = value + Fraction((-1)**terms,
                                 (2*terms+1)*reciprocal**(2*terms+1))
    return min(value, next_value), max(value, next_value)


def torsion_nondensity_exact() -> dict:
    # Components represent 15*z: real + i*(a*sqrt(6)+b*sqrt(2)+c*sqrt(11)).
    components = [(15, 0, 0, 0), (15, 0, 0, 0),
                  (-3, 6, 0, 0), (-3, -6, 0, 0),
                  (-5, 0, 10, 0), (-5, 0, -10, 0),
                  (-7, 0, 0, 4), (-7, 0, 0, -4)]
    roots = [(6, 2449489742), (2, 1414213562), (11, 3316624790)]
    intervals = []
    for integer, numerator in roots:
        lo, hi = Fraction(numerator, 10**9), Fraction(numerator+1, 10**9)
        assert lo*lo < integer < hi*hi
        intervals.append((lo, hi))
    assert [sum(row[j] for row in components) for j in range(4)] == [0]*4
    attained = []
    for mask in range(1, 255):
        total = [sum(components[j][a] for j in range(8) if mask >> j & 1)
                 for a in range(4)]
        assert total != [0]*4
        lo = hi = Fraction(0)
        for coefficient, bounds in zip(total[1:], intervals):
            lower, upper = interval_scale(*bounds, Fraction(coefficient, 15))
            lo += lower
            hi += upper
        real = Fraction(total[0], 15)
        lower_squared_norm = real*real + distance_to_zero(lo, hi)**2
        assert lower_squared_norm >= Fraction(1, 225)
        if lower_squared_norm == Fraction(1, 225):
            attained.append(mask)
    # Exact Machin enclosure of pi, with alternating Taylor remainders.
    a5, a239 = atan_interval(5), atan_interval(239)
    pi_lo, pi_hi = 16*a5[0]-4*a239[1], 16*a5[1]-4*a239[0]
    # Outward rounding keeps subsequent exact Taylor fractions manageable.
    denominator = 10**20
    lower_scaled, upper_scaled = denominator*pi_lo, denominator*pi_hi
    pi_lo = Fraction(lower_scaled.numerator//lower_scaled.denominator, denominator)
    pi_hi = Fraction(-((-upper_scaled.numerator)//upper_scaled.denominator), denominator)
    assert Fraction(3) < pi_lo < pi_hi < Fraction(4)
    target_real = Fraction(-1, 5)
    target_imag = interval_scale(*intervals[0], Fraction(2, 5))
    closest_lower_squared = Fraction(100)
    for j in range(210):
        reduced = j if j <= 105 else j - 210
        theta_lo, theta_hi = interval_scale(pi_lo, pi_hi, Fraction(reduced, 105))
        center = (theta_lo+theta_hi)/2
        radius = (theta_hi-theta_lo)/2
        cosine = sum(((-1)**a * center**(2*a)/math.factorial(2*a)
                      for a in range(23)), Fraction(0))
        sine = sum(((-1)**a * center**(2*a+1)/math.factorial(2*a+1)
                    for a in range(23)), Fraction(0))
        cosine_error = radius + Fraction(4**46, math.factorial(46))
        sine_error = radius + Fraction(4**47, math.factorial(47))
        re_lo = cosine-cosine_error-target_real
        re_hi = cosine+cosine_error-target_real
        im_lo = sine-sine_error-target_imag[1]
        im_hi = sine+sine_error-target_imag[0]
        lower = distance_to_zero(re_lo, re_hi)**2 + distance_to_zero(im_lo, im_hi)**2
        assert lower > Fraction(1, 40000)
        closest_lower_squared = min(closest_lower_squared, lower)
    return {"proper_subsets_exactly_checked": 254,
            "proper_subsum_modulus_lower_bound": "1/15",
            "minimum_subsum_attaining_masks": attained,
            "torsion_ratios_exactly_checked": 210,
            "distance_of_z3_to_mu210_lower_bound": "1/200",
            "uniform_distance_from_degree2_torsion_polynomials_lower_bound": "1/400",
            "diagnostic_minimum_root_distance_lower_bound":
                math.sqrt(float(closest_lower_squared))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"scope": "Exact finite identities; asymptotics require artifact proofs",
              "exact": exact_checks(),
              "clt_diagnostics": [first_level_diagnostic(k, lam)
                                  for lam in (0.5, 1.0, 2.0)
                                  for k in (31, 101, 301, 1001)],
              "contraction_diagnostics": contraction_diagnostics(),
              "torsion_nondensity_exact": torsion_nondensity_exact()}
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload)


if __name__ == "__main__":
    main()
