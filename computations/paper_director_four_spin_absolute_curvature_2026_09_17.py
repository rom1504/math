"""Exact finite combinatorics for the analytic absolute-curvature falsifier.

Enumerates all 64 full signings and 8 external sign patterns. Gaussian
tail/density inequalities are proved in the accompanying artifact, not
inferred from this finite computation.
"""
import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

edges = list(itertools.combinations(range(4), 2))
words = list(itertools.product((-1, 1), repeat=4))
patterns = list(itertools.product((-1, 1), repeat=3))
caps = []
checks = 0
for signs in itertools.product((-1, 1), repeat=6):
    def energy(word):
        return sum(a * word[i] * word[j] for a, (i, j) in zip(signs, edges))

    caps.append(max(abs(energy(word)) for word in words))
    weight = {z: abs(energy((1,) + z)) for z in patterns}
    unit_offsets = 0
    for z in patterns:
        field = tuple((4 + i) * z[i] for i in range(3))
        plus = max(weight[t] + sum(a*b for a, b in zip(field, t)) for t in patterns)
        minus = max(weight[t] - sum(a*b for a, b in zip(field, t)) for t in patterns)
        offset = Fraction(plus - minus, 2)
        negative = tuple(-a for a in z)
        assert offset == Fraction(weight[z] - weight[negative], 2)
        assert abs(offset) in (1, 3)
        unit_offsets += (abs(offset) == 1)
        for first in (Fraction(-7, 2), Fraction(-1, 4), Fraction(0), Fraction(5, 2)):
            full = max(abs(energy(word) + first*word[0]
                           + sum(a*b for a, b in zip(field, word[1:]))) for word in words)
            assert full == Fraction(plus + minus, 2) + abs(first + offset)
        checks += 1
    assert unit_offsets >= 6

numerator = Fraction(3, 4)*Fraction(116, 125)*Fraction(1, 20)*Fraction(3, 4)
numerator -= Fraction(9, 125)*Fraction(1, 10)
assert numerator == Fraction(189, 10000)
curvature = numerator * 16
assert curvature == Fraction(189, 625)
assert min(caps) == 4
result = dict(status="exact rational/enumeration PASS", full_signings=len(caps),
              external_sign_patterns_checked=checks, exact_M4=min(caps),
              algebraic_curvature_lower=str(curvature),
              scope="full absolute four-spin cap; unequal independent Gaussian fields")
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path)
args = parser.parse_args()
encoded = json.dumps(result, indent=2) + "\n"
if args.output:
    args.output.write_text(encoded)
print(encoded, end="")
