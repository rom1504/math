"""Independently rebuild the finite eight-block dual certificate.

Only the 71 proposed integer coefficients are imported. All patterns,
affine symmetries, orbit averages, masses, and rational log enclosures
are rebuilt without NumPy, sparse matrices, or optimizer arithmetic.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import factorial

from continued_feedback_eight_block_profile_lp_2026_09_06 import CERTIFIED_DUAL


def parity(x):
    return bin(x).count('1') % 2


FREQUENCIES = (0, 1, 2, 4)


def projection(x):
    return tuple(sum((-1) ** parity(f & u) * x[u] for u in range(8))
                 for f in FREQUENCIES)


@lru_cache(None)
def log_interval(x):
    power = 0
    while x < 1:
        x *= 2
        power -= 1
    while x >= 2:
        x /= 2
        power += 1

    def atanh_interval(value):
        z = (value - 1) / (value + 1)
        lower = 2 * sum((z ** (2*j+1) / (2*j+1) for j in range(35)), F(0))
        return lower, lower + 2*z**71 / (71*(1-z*z))

    lower, upper = atanh_interval(x)
    l2, u2 = atanh_interval(F(2))
    if power < 0:
        return lower + power*u2, upper + power*l2
    return lower + power*l2, upper + power*u2


def main():
    patterns = tuple(product((-1, 0, 1), repeat=8))
    masses = Counter()
    for x in patterns:
        masses[projection(x)] += 2**x.count(0) * 15**(8-x.count(0))
    assert len(masses) == 1537
    assert sum(masses.values()) == 32**8

    # Every permutation of the four affine-basis frequencies determines
    # exactly one affine automorphism preserving this tetrahedron.
    actions = set()
    for images in permutations(FREQUENCIES):
        offset = images[0]
        columns = tuple(images[j] ^ offset for j in range(1, 4))
        coordinate_permutation = tuple(
            sum(parity(c & u) << j for j, c in enumerate(columns))
            for u in range(8))
        assert set(coordinate_permutation) == set(range(8))
        for global_sign in (-1, 1):
            signs = tuple(global_sign * (-1)**parity(offset & u) for u in range(8))
            actions.add((coordinate_permutation, signs))
    assert len(actions) == 48

    basis = []
    for degree in range(5):
        for coordinates in combinations(range(8), degree):
            for symbols in product((-1, 1), repeat=degree):
                basis.append(tuple(zip(coordinates, symbols)))
    assert len(basis) == 1697

    terms = [(coefficient, basis[index]) for index, coefficient in CERTIFIED_DUAL.items()]

    def integer_polynomial(x):
        return sum(coefficient for coefficient, conditions in terms
                   if all(x[i] == symbol for i, symbol in conditions))

    seen = set()
    orbit_count = 0
    minimum_margin = None
    for x in patterns:
        if x in seen:
            continue
        orbit = {tuple(signs[u] * x[permutation[u]] for u in range(8))
                 for permutation, signs in actions}
        assert not (orbit & seen)
        assert all(z.count(0) == x.count(0) for z in orbit)
        output_masses = {masses[projection(z)] for z in orbit}
        assert len(output_masses) == 1
        mass = output_masses.pop()
        polynomial_average = F(sum(integer_polynomial(z) for z in orbit),
                               len(orbit) * 10**10)
        _, log_upper = log_interval(F(mass, 32**8))
        margin = polynomial_average - log_upper
        assert margin >= 0
        minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
        seen.update(orbit)
        orbit_count += 1
    assert len(seen) == 6561 and orbit_count == 255

    objective = sum((F(coefficient, 10**10) * F(15, 32)**len(conditions)
                     for coefficient, conditions in terms), F(0))
    assert objective == F(-65899554490727779, 10485760000000000)
    _, log2_upper = log_interval(F(2))
    sqrt15_lower = F(3872983346207416, 10**15)
    assert sqrt15_lower**2 < 15
    assert sum((F(16)**j / factorial(j) for j in range(50)), F(0)) > 8000000
    tilted_upper = F(15, 16)*log2_upper + objective/8 + F(1, 16000000) + 4 - sqrt15_lower
    assert tilted_upper < F(-874170918, 10**11)
    print(json.dumps({
        'status': 'independent exact reconstruction passed',
        'patterns': len(patterns), 'outputs': len(masses),
        'actions': len(actions), 'orbits': orbit_count, 'basis_rows': len(basis),
        'nonzero_coefficients': len(terms), 'objective_exact': str(objective),
        'minimum_pointwise_orbit_margin_lower': float(minimum_margin),
        'tilted_upper': float(tilted_upper),
        'scope': 'four-wise marginal dual and near-flat class only; no full cap claim'
    }, indent=2))


if __name__ == '__main__':
    main()
