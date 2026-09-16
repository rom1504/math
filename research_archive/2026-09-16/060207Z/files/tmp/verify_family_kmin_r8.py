#!/usr/bin/env python3
"""Exact scalar audit of K_min for the family (10.416)."""

from fractions import Fraction

from kmin_carleson_r8 import structured_family_k


def closed_form(m):
    if m == 4:
        return Fraction(12, 5)
    return Fraction(12 * m * m, 5 * m * m + 3 * m - 8)


for m in (4, 16, 64, 256, 1024, 4096):
    exact, service = structured_family_k(m)
    assert exact == closed_form(m)
    if m >= 16:
        slope = Fraction(2 * (m * m - 4), 3)
        q = exact / 2
        child_imbalance = m * (m - 2)
        # The minimizing root edge is in the unsaturated linear regime both
        # for its two equal local buckets and for its clique child.
        assert q < 2
        assert slope * q < child_imbalance
        assert service(q) == slope * q
    print("m", m, "K_min", exact, "float", float(exact))

print("limit 12/5")
