"""Exact checks for the transitive spin-amplification obstruction.

This script writes no files.  The ternary-tree recurrences are exact
Fractions; small cubes independently check degree, influence, and the
orbit-quotient commutator identity by enumeration.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json


def ternary_value(x):
    if len(x) == 1:
        return x[0]
    block = len(x) // 3
    a, b, c = [ternary_value(x[j * block:(j + 1) * block]) for j in range(3)]
    return (a * b + a * c + b * c - 1) // 2


def orbit_key(x):
    if len(x) == 1:
        return x[0]
    block = len(x) // 3
    return tuple(sorted(orbit_key(x[j * block:(j + 1) * block]) for j in range(3)))


def walsh_integers(values):
    out = list(values)
    half = 1
    while half < len(out):
        for start in range(0, len(out), 2 * half):
            for offset in range(half):
                a, b = out[start + offset], out[start + offset + half]
                out[start + offset] = a + b
                out[start + offset + half] = a - b
        half *= 2
    return out


def enumerate_depth(depth):
    n = 3 ** depth
    cube = list(product((-1, 1), repeat=n))
    values = [ternary_value(x) for x in cube]
    transformed = walsh_integers(values)
    denominator = 2 ** (2 * n)
    influence = sum(Fraction(bin(mask).count("1") * coefficient ** 2, denominator)
                    for mask, coefficient in enumerate(transformed))
    degree = max(bin(mask).count("1") for mask, coefficient in enumerate(transformed) if coefficient)
    groups = defaultdict(list)
    for x in cube:
        groups[orbit_key(x)].append(x)
    keys = sorted(groups)
    edge_counts = {}
    for key, points in groups.items():
        reference = None
        for x in points:
            neighbors = Counter()
            for i in range(n):
                y = x[:i] + (-x[i],) + x[i + 1:]
                neighbors[orbit_key(y)] += 1
            if reference is None:
                reference = neighbors
            assert neighbors == reference
        edge_counts[key] = reference
    column_energy = Fraction(0)
    edge_energy = Fraction(0)
    edge_influence = Fraction(0)
    minimum_amplification = None
    for i, left in enumerate(keys):
        w_left = Fraction(len(groups[left]), 2 ** n)
        f_left = ternary_value(groups[left][0])
        for right in keys:
            a = edge_counts[left][right]
            b = edge_counts[right][left]
            delta = ternary_value(groups[right][0]) - f_left
            column_energy += w_left * a * b * delta ** 2 / 4
        for right in keys[i + 1:]:
            a = edge_counts[left][right]
            b = edge_counts[right][left]
            if not a:
                continue
            assert w_left * a == Fraction(len(groups[right]), 2 ** n) * b
            delta = ternary_value(groups[right][0]) - f_left
            flow = w_left * a
            edge_energy += flow * (a + b) * delta ** 2 / 4
            edge_influence += flow * delta ** 2 / 2
            amplification = Fraction(a + b, 2)
            minimum_amplification = min(minimum_amplification or amplification, amplification)
    assert column_energy == edge_energy
    assert edge_influence == influence
    assert edge_energy <= degree ** 2
    assert degree == 2 ** depth
    return {
        "depth": depth, "variables": n, "degree": degree,
        "vertex_orbits": len(groups), "influence": str(influence),
        "orbit_commutator_energy": str(edge_energy),
        "minimum_orbit_edge_amplification": str(minimum_amplification),
    }


def main():
    mean, influence = Fraction(0), Fraction(1)
    rows = []
    for depth in range(13):
        n, degree = 3 ** depth, 2 ** depth
        ratio = n * influence / degree ** 2
        assert ratio >= Fraction(9, 8) ** depth
        rows.append({"depth": depth, "variables": n, "degree": degree,
                     "mean": float(mean), "influence": float(influence),
                     "n_influence_over_degree_squared": float(ratio)})
        influence *= Fraction(3, 2) * (1 + mean ** 2)
        mean = (3 * mean ** 2 - 1) / 2
    print(json.dumps({"exact_small_cube_checks": [enumerate_depth(t) for t in range(3)],
                      "exact_fraction_recurrence_diagnostics": rows}, indent=2))


if __name__ == "__main__":
    main()
