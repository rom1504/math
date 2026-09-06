"""Exact checks for exceptional and random very thin tensor restrictions.

No runtime input files or output files. All finite data are constructed here.
"""
import itertools
import json
import math
from fractions import Fraction as F

import numpy as np


def full_seed(d, code):
    b = np.ones((d, d), dtype=np.int64)
    for bit, (i, j) in enumerate(itertools.combinations_with_replacement(range(d), 2)):
        b[i, j] = b[j, i] = 1 if code & (1 << bit) else -1
    return b


def hollow_target(n, code):
    a = np.zeros((n, n), dtype=np.int64)
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        a[i, j] = a[j, i] = 1 if code & (1 << bit) else -1
    return a


def anchored(b):
    return b * b[:, 0, None] * b[None, 0, :] * b[0, 0]


def encode(b, a):
    h = anchored(b)
    pairs = np.argwhere(h == -1)
    if not len(pairs):
        return None
    left, right = map(int, pairs[0])
    n = len(a)
    words = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        words[i, i] = left
        for j in range(i):
            prior = int(np.prod(h[words[j, :j], words[i, :j]]))
            words[i, j] = right if prior != a[j, i] else 0
    selected = np.ones((n, n), dtype=np.int64)
    for k in range(n):
        selected *= b[words[:, k, None], words[None, :, k]]
    gauge = np.prod(b[words, 0], axis=1)
    switched = selected * gauge[:, None] * gauge[None, :] * int(b[0, 0]) ** n
    np.fill_diagonal(switched, 0)
    assert len(set(map(tuple, words))) == n
    assert np.array_equal(switched, a)
    return words


def beta(b):
    spins = np.array(list(itertools.product((-1, 1), repeat=len(b))), dtype=np.int64)
    return int(np.max(np.sum(np.abs(spins @ b), axis=1)))


def graph_fourier(b, n):
    edges = list(itertools.combinations(range(n), 2))
    law = [F(0) for _ in range(1 << len(edges))]
    for colors in itertools.product(range(len(b)), repeat=n):
        code = sum((b[colors[i], colors[j]] == -1) << k for k, (i, j) in enumerate(edges))
        law[code] += F(1, len(b) ** n)
    coefficients = []
    for mask in range(len(law)):
        coefficients.append(sum((prob * (-1 if bin(mask & code).count('1') % 2 else 1)
                                 for code, prob in enumerate(law)), F(0)))
    theta = F(beta(b), len(b) ** 2)
    assert coefficients[0] == 1
    assert all(abs(value) <= theta for value in coefficients[1:])
    # Exact two-coordinate XOR convolution checks tensorization and Parseval.
    convolution = [sum((law[a] * law[a ^ c] for a in range(len(law))), F(0))
                   for c in range(len(law))]
    for mask, coefficient in enumerate(coefficients):
        got = sum((prob * (-1 if bin(mask & code).count('1') % 2 else 1)
                   for code, prob in enumerate(convolution)), F(0))
        assert got == coefficient ** 2
    tv = sum((abs(prob - F(1, len(law))) for prob in convolution), F(0)) / 2
    assert tv ** 2 <= sum((x ** 4 for x in coefficients[1:]), F(0)) / 4
    return theta, tv


def greedy(a):
    x = np.ones(len(a), dtype=np.int64)
    increments = []
    for i in range(1, len(a)):
        field = int(a[i, :i] @ x[:i])
        x[i] = 1 if field >= 0 else -1
        increments.append(abs(field))
    assert sum(increments) == int(x @ a @ x) // 2
    return tuple(increments)


def main():
    encodings = 0
    rank_one = 0
    for d in (2, 3):
        for code in range(1 << (d * (d + 1) // 2)):
            b = full_seed(d, code)
            if np.all(anchored(b) == 1):
                assert beta(b) == d * d
                rank_one += 1
                continue
            assert beta(b) < d * d
            for target in range(64):
                assert encode(b, hollow_target(4, target)) is not None
                encodings += 1
    rng = np.random.default_rng(20260906)
    for d in (4, 6, 16):
        for _ in range(30):
            b = rng.choice((-1, 1), size=(d, d)).astype(np.int64)
            b = np.triu(b) + np.triu(b, 1).T
            n = 32
            a = rng.choice((-1, 1), size=(n, n)).astype(np.int64)
            a = np.triu(a, 1) + np.triu(a, 1).T
            assert encode(b, a) is not None
            encodings += 1
    fourier = []
    for b in (np.array([[1, 1], [1, -1]]), np.ones((3, 3), dtype=int) - 2 * np.eye(3, dtype=int)):
        theta, tv = graph_fourier(b, 4)
        fourier.append({'d': len(b), 'theta': str(theta), 'r2_tv': str(tv)})
    n = 5
    energies = [sum(greedy(hollow_target(n, code))) for code in range(1 << (n * (n - 1) // 2))]
    mean = F(sum(energies), len(energies))
    mu = sum((F(m * math.comb(m - 1, (m - 1) // 2), 2 ** (m - 1)) for m in range(1, n)), F(0))
    assert mean == mu
    variance = sum(((F(x) - mean) ** 2 for x in energies), F(0)) / len(energies)
    expected_variance = sum((F(m) - F(m * math.comb(m - 1, (m - 1) // 2), 2 ** (m - 1)) ** 2
                             for m in range(1, n)), F(0))
    assert variance == expected_variance
    print(json.dumps({'encodings_checked': encodings, 'rank_one_seeds_checked': rank_one,
                      'fourier_checks': fourier, 'greedy_exact_n': n,
                      'greedy_mean': str(mean), 'greedy_variance': str(variance),
                      'greedy_limit': (2 / 3) * math.sqrt(2 / math.pi)}, indent=2))


if __name__ == '__main__':
    main()
