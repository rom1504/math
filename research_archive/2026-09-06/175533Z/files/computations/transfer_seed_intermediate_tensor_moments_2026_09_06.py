"""Exact polynomial and finite-moment checks for random tensor restrictions.

No input/output files. Rational interpolation supplies a global polynomial
minorant by the separate Hermite-remainder proof, not by a grid assertion.
"""
from fractions import Fraction as F
import itertools
import json
import math


NODES = [F(x, 10000) for x in (7991, 16067, 24324, 32891, 41962, 51901, 63639)]


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return out


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def evaluate(a, t):
    result = F(0)
    for coefficient in reversed(a):
        result = result * t + coefficient
    return result


def polynomial():
    """q(0)=0, q(s_j^2)=s_j, q'(s_j^2)=1/(2s_j), deg(q)<=14."""
    squares = [x * x for x in NODES]
    out = [F(0)]
    for j, (s, t) in enumerate(zip(NODES, squares)):
        lagrange = [F(1)]
        logarithmic_derivative = F(0)
        for ell, other in enumerate(squares):
            if ell != j:
                lagrange = multiply(lagrange, [-other / (t - other), 1 / (t - other)])
                logarithmic_derivative += 1 / (t - other)
        slope = -1 / (2 * s) - 2 * s * logarithmic_derivative
        local = multiply([F(0), 1 / t], multiply(lagrange, lagrange))
        local = multiply(local, [s - t * slope, slope])
        out = add(out, local)
    assert len(out) == 15 and out[0] == 0
    derivative = [k * out[k] for k in range(1, len(out))]
    for s in NODES:
        assert evaluate(out, s * s) == s
        assert evaluate(derivative, s * s) == 1 / (2 * s)
    return out


def gaussian_moment(k):
    return math.factorial(2 * k) // (2 ** k * math.factorial(k))


def scalar_moment_checks():
    # Directly check the absolute-column-correlation estimate eta_q <= theta.
    seeds = [[[1, 1], [1, -1]], [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]]
    result = []
    for b in seeds:
        d = len(b)
        beta = max(sum(abs(sum(b[i][j] * x[i] for i in range(d))) for j in range(d))
                   for x in itertools.product((-1, 1), repeat=d))
        theta = F(beta, d * d)
        values = {}
        for q in (2, 4, 6):
            total = 0
            for colors in itertools.product(range(d), repeat=q):
                total += abs(sum(math.prod(b[a][z] for z in colors) for a in range(d)))
            eta = F(total, d ** (q + 1))
            assert eta <= theta
            values[str(q)] = str(eta)
        result.append({'d': d, 'theta': str(theta), 'eta_even': values})
    return result


def adaptive_moment_checks():
    b = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
    theta = F(5, 9)
    checks = 0
    for depth in (1, 2):
        words = list(itertools.product(range(3), repeat=depth))
        matrix = [[math.prod(b[x][y] for x, y in zip(u, v)) for v in words] for u in words]
        ambient = len(words)
        count = 3
        parity_counts = {}
        iid_moments = {}
        for degree in (2, 4, 6):
            counts = [0] * (1 << count)
            for indices in itertools.product(range(count), repeat=degree):
                mask = 0
                for index in indices:
                    mask ^= 1 << index
                counts[mask] += 1
            parity_counts[degree] = counts
            iid_moments[degree] = F(sum(sum(x) ** degree for x in itertools.product((-1, 1), repeat=count)),
                                     2 ** count * count ** (degree // 2))
        sum_errors = {degree: F(0) for degree in parity_counts}
        for history in itertools.product(range(ambient), repeat=count):
            signs = [1] * count
            for i in range(1, count):
                field = sum(signs[j] * matrix[history[i]][history[j]] for j in range(i))
                signs[i] = 1 if field >= 0 else -1
            correlations = [F(0)] * (1 << count)
            for mask in range(1, 1 << count):
                correlations[mask] = F(sum(math.prod(matrix[w][history[j]] for j in range(count)
                                                      if mask & (1 << j)) for w in range(ambient)), ambient)
            fields = [sum(signs[j] * matrix[w][history[j]] for j in range(count)) for w in range(ambient)]
            for degree, counts in parity_counts.items():
                normalized = F(sum(z ** degree for z in fields), ambient * count ** (degree // 2))
                error = sum((counts[mask] * abs(correlations[mask]) for mask in range(1, len(counts))), F(0))
                error /= count ** (degree // 2)
                assert abs(normalized - iid_moments[degree]) <= error
                sum_errors[degree] += error
                checks += 1
        for degree, total in sum_errors.items():
            assert total / ambient ** count <= count ** (degree // 2) * theta ** depth
    return checks


def main():
    q = polynomial()
    normal_expectation = sum((a * gaussian_moment(k) for k, a in enumerate(q)), F(0))
    assert normal_expectation > F(377, 500)
    assert sum(map(abs, q)) < 9
    # Sanity checks only; global validity follows from the Hermite remainder.
    for numerator in range(501):
        s = F(numerator, 50)
        assert evaluate(q, s * s) <= s
    display_scale = 10 ** 12
    outward_lower = F(normal_expectation.numerator * display_scale // normal_expectation.denominator,
                      display_scale)
    assert outward_lower > F(377, 500)
    print(json.dumps({'nodes': list(map(str, NODES)), 'degree_in_x': 28,
                      'gaussian_expectation_downward_1e12': str(outward_lower),
                      'gaussian_expectation_display': float(normal_expectation),
                      'certified_greedy_limit_lower': str(F(377, 750)),
                      'sum_absolute_coefficients_display': float(sum(map(abs, q))),
                      'scalar_correlation_checks': scalar_moment_checks(),
                      'adaptive_conditional_moment_checks': adaptive_moment_checks()}, indent=2))


if __name__ == '__main__':
    main()
