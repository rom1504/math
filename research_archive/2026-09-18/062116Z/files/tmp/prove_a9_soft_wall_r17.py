#!/usr/bin/env python3
"""Exact polynomial certificate for the A9 fixed-temperature cavity wall."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, sqrt

from selected_child_r17 import A9


def energies(vertices):
    vertices = tuple(vertices)
    n = len(vertices)
    out = []
    for tail in product((-1, 1), repeat=n - 1):
        x = (1,) + tail
        e = sum(2 * A9[vertices[i]][vertices[j]] * x[i] * x[j]
                for i in range(n) for j in range(i + 1, n))
        out.extend((e, -e))  # the two orientations
    return out


def deficit_poly(vertices):
    es = energies(vertices)
    q = max(es)
    counts = Counter(q - e for e in es)
    return q, counts


def bernstein_coefficients(poly, degree):
    """Power polynomial coefficients converted to degree-n Bernstein basis."""
    a = [poly.get(k, 0) for k in range(degree + 1)]
    # x^k = sum_{j=k}^n C(j,k)/C(n,k) B_{j,n}(x)
    return [sum(Fraction(a[k] * comb(j, k), comb(degree, k))
                for k in range(j + 1)) for j in range(degree + 1)]


def pretty(poly):
    return " + ".join(f"{poly[d]}s^{d}" for d in sorted(poly))


def main():
    root = tuple(range(9))
    q9, z9 = deficit_poly(root)
    assert q9 == 24
    print("Z_A9(s)=", pretty(z9))
    for missing in root:
        child = tuple(i for i in root if i != missing)
        qc, zc = deficit_poly(child)
        assert qc == 24
        # kappa < 3 is exactly s^3 Z_parent < 2 Z_child, since the
        # normalized state counts differ by a factor two.
        p = Counter({d: 2 * c for d, c in zc.items()})
        for d, c in z9.items():
            p[d + 3] -= c
        degree = max(p)
        b = bernstein_coefficients(p, degree)
        print("missing", missing, "Z_child=", pretty(zc))
        print("  2Zc-s^3Z9=", pretty(p))
        print("  degree", degree, "min Bernstein coeff", min(b),
              "all nonnegative", all(x >= 0 for x in b))
        # Nonnegative Bernstein coefficients prove p(s)>=0 on [0,1].
        # The coefficient of B_{0,degree} is 4, so p(s)>0 for 0<=s<1.
        assert all(x >= 0 for x in b) and b[0] > 0
    # The alpha-centered zero-error threshold exceeds three exactly:
    # alpha=24/9^(3/2)=8/9 and Delta r^(3/2)=27-16sqrt(2).
    threshold = Fraction(8, 9) * 27 - Fraction(128, 9) * sqrt(2)
    assert threshold > 3
    print("alpha threshold=24-(128/9)sqrt(2) ~", threshold, ">3")
    print("PASS: every A9 cavity reward is <3 for every beta>0")


if __name__ == "__main__":
    main()
