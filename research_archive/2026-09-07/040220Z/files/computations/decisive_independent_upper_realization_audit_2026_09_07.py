"""Independent integer weave and rational final-cap checks. No file writes."""

from fractions import Fraction as F
from itertools import product
from math import isqrt
import json
import random


def logarithm_bounds(x):
    x = F(x)
    exponent = 0
    while x < 1:
        x *= 2
        exponent -= 1
    while x > 2:
        x /= 2
        exponent += 1

    def series(y):
        w = (y - 1) / (y + 1)
        value = 2 * sum(w ** (2*j+1) / (2*j+1) for j in range(64))
        tail = 2*w**129 / (129*(1-w*w))
        return value, value + tail

    low, high = series(x)
    low2, high2 = series(F(2))
    if exponent >= 0:
        return low + exponent*low2, high + exponent*high2
    return low + exponent*high2, high + exponent*low2


def square_root_bounds(x):
    x = F(x)
    scale = 10**18
    integer = isqrt(x.numerator*scale**2 // x.denominator)
    low, high = F(integer, scale), F(integer+1, scale)
    assert low*low <= x <= high*high
    return low, high


def cap_check():
    p = F(31, 32)
    sl, su = square_root_bounds(257)
    rl, ru = (sl-1)/16, (su-1)/16
    gl = -4*(1-rl) + logarithm_bounds(1-ru*ru)[0]/4
    gu = -4*(1-ru) + logarithm_bounds(1-rl*rl)[1]/4
    pl, pu = square_root_bounds(p)
    ll, lu = logarithm_bounds(2)
    nl, nu = 4+p*ll+gl, 4+p*lu+gu
    assert nl > 0
    low, high = nl/(8*pu), nu/(8*pl)
    published = F(7787631971809, 15748015748016)
    decimal_cap = F(494515125, 10**9)
    assert high < published < decimal_cap
    return {
        "independent_lower_float": float(low),
        "independent_upper_float": float(high),
        "published_upper": str(published),
        "exact_published_decimal_gap": str(decimal_cap-published),
    }


def transpose(a):
    return [list(row) for row in zip(*a)]


def check_hadamard(h):
    m = len(h)
    assert all(abs(a) == 1 for row in h for a in row)
    assert all(sum(a*b for a, b in zip(h[i], h[j])) == m*(i == j)
               for i in range(m) for j in range(m))


def h12():
    squares = {i*i % 11 for i in range(1, 11)}
    def character(i):
        return 0 if i % 11 == 0 else (1 if i % 11 in squares else -1)
    h = [[1]*12]
    h += [[-1]+[int(i == j)+character(i-j) for j in range(11)]
          for i in range(11)]
    check_hadamard(h)
    assert h != transpose(h)
    return h


def permuted(h, rng):
    m = len(h)
    indices = rng.sample(range(m), m)
    signs = [rng.choice((-1, 1)) for _ in range(m)]
    return [[row[indices[j]]*signs[j] for j in range(m)] for row in h]


def weave_check(basis, retained, rng, exhaustive):
    m = len(basis)
    # Physical bases are transposes of the recursive input transformations.
    bases = [transpose(permuted(basis, rng)) for _ in range(m)]
    selectors = [rng.sample(range(m), retained) for _ in range(m)]
    outer = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i, m):
            outer[i][j] = outer[j][i] = rng.choice((-1, 1))
    labels = [(i, a) for i in range(m) for a in selectors[i]]
    matrix = [[outer[i][j]*bases[i][a][j]*bases[j][b][i]
               for j, b in labels] for i, a in labels]
    n = len(labels)
    assert matrix == transpose(matrix)
    assert all(abs(a) == 1 for row in matrix for a in row)
    diagonal = sum(matrix[i][i] for i in range(n))
    assert diagonal == retained*sum(outer[i][i] for i in range(m))
    words = (product((-1, 1), repeat=n) if exhaustive else
             ([rng.choice((-1, 1)) for _ in range(n)] for _ in range(128)))
    count = 0
    for x in words:
        h = [[sum(bases[i][a][j]*x[i*retained+u]
                  for u, a in enumerate(selectors[i])) for j in range(m)]
             for i in range(m)]
        assert sum(a*a for row in h for a in row) == m*m*retained
        energy = sum(matrix[i][j]*x[i]*x[j] for i in range(n) for j in range(n))
        assert energy == sum(outer[i][j]*h[i][j]*h[j][i]
                             for i in range(m) for j in range(m))
        for sigma in (-1, 1):
            defect = sum((h[i][j]-sigma*outer[i][j]*h[j][i])**2
                         for i in range(m) for j in range(m))
            assert defect == 2*(m*m*retained-sigma*energy)
            for i in range(m):
                for j in range(i+1, m):
                    # This exact doubling fixes the exponential kernel scale.
                    a, b, s = h[i][j], h[j][i], sigma*outer[i][j]
                    assert (a-s*b)**2+(b-s*a)**2 == 2*(a-s*b)**2
        assert abs(energy-diagonal) <= abs(energy)+n
        count += 1
    return {"m": m, "k": retained, "checked_words": count,
            "exhaustive": exhaustive}


def main():
    rng = random.Random(730091)
    small = h12()
    left, right = permuted(small, rng), permuted(small, rng)
    # Exact unnormalized version of a single binary recursive gate.
    recursive = [row+row for row in left]+[row+[-a for a in row] for row in right]
    recursive = permuted(recursive, rng)
    check_hadamard(recursive)
    result = {
        "status": "independent exact checks passed",
        "cap": cap_check(),
        "nonsymmetric_terminal": weave_check(small, 1, rng, True),
        "recursive_nonsymmetric_terminal": weave_check(recursive, 2, rng, False),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
