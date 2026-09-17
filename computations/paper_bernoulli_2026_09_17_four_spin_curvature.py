"""Exact finite identity replay for every full order-four child."""

import itertools
import json
from collections import Counter
from fractions import Fraction

edges = list(itertools.combinations(range(4), 2))
spins = list(itertools.product([-1, 1], repeat=4))
patterns = list(itertools.product([-1, 1], repeat=3))
caps = Counter()
identity_checks = 0
sign_checks = 0
for signs in itertools.product([-1, 1], repeat=6):
    def energy(y):
        return sum(a * y[i] * y[j] for (i, j), a in zip(edges, signs))

    caps[max(abs(energy(y)) for y in spins)] += 1
    row_one_count = 0
    for z in patterns:
        fields = [4 * t for t in z]
        weights = {w: abs(energy((1,) + w)) for w in patterns}
        ap = max(weights[w] + sum(a * b for a, b in zip(fields, w)) for w in patterns)
        am = max(weights[w] - sum(a * b for a, b in zip(fields, w)) for w in patterns)
        offset = Fraction(ap - am, 2)
        row = sum(signs[j] * z[j] for j in range(3))
        h = energy((1,) + z) - row
        expected_offset = Fraction(abs(h + row) - abs(h - row), 2)
        assert offset == expected_offset
        assert abs(offset) in [1, 3]
        row_one_count += abs(row) == 1
        sign_checks += 1
        for g1 in [-10, -4, -1, 0, 1, 4, 10]:
            field = [g1] + fields
            absolute_maximum = max(abs(energy(y) + sum(a * b for a, b in zip(field, y))) for y in spins)
            formula = Fraction(ap + am, 2) + abs(g1 + offset)
            assert absolute_maximum == formula
            identity_checks += 1
    assert row_one_count == 6

curvature_lower = 16 * (Fraction(3, 4) * Fraction(116, 125) * Fraction(3, 80)
                        - Fraction(9, 125) * Fraction(1, 10))
assert curvature_lower == Fraction(189, 625)
print(json.dumps({
    "full_children": 64,
    "full_cube_cap_histogram": dict(sorted(caps.items())),
    "sign_pattern_checks": sign_checks,
    "full_absolute_identity_checks": identity_checks,
    "strict_curvature_lower_bound": str(curvature_lower),
    "all_pass": True
}, indent=2))
