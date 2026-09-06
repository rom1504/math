#!/usr/bin/env python3
"""Checks the Wave 32 row-noise identities and Paley block wall."""

from itertools import combinations, product
from math import comb

import numpy as np


def row_cost(a, x):
    y = a @ x
    return float(y @ y)


def random_signing(n, rng):
    u = rng.choice((-1.0, 1.0), size=(n, n))
    a = np.triu(u, 1)
    return a + a.T


def check_product_noise():
    rng = np.random.default_rng(3204)
    n = 7
    a = random_signing(n, rng)
    b = a @ a
    x = rng.choice((-1.0, 1.0), size=n)
    p = 0.27
    theta = 1 - 2 * p
    values = []
    weights = []
    for eps0 in product((-1.0, 1.0), repeat=n):
        eps = np.asarray(eps0)
        neg = int(np.count_nonzero(eps < 0))
        weights.append(p**neg * (1-p)**(n-neg))
        values.append(row_cost(a, x * eps))
    values = np.asarray(values)
    weights = np.asarray(weights)
    mean = float(weights @ values)
    var = float(weights @ ((values - mean) ** 2))
    target_mean = theta**2 * row_cost(a, x) + (1-theta**2) * n * (n-1)
    assert abs(mean-target_mean) < 2e-10

    v = 1-theta**2
    q = 2 * b * np.outer(x, x)
    np.fill_diagonal(q, 0)
    ell = q.sum(axis=1)
    target_var = v * theta**2 * float(ell @ ell)
    target_var += v**2 * sum(q[i, j]**2 for i in range(n) for j in range(i+1, n))
    assert abs(var-target_var) < 2e-8


def check_fixed_spheres():
    rng = np.random.default_rng(3214)
    n = 9
    a = random_signing(n, rng)
    x = rng.choice((-1.0, 1.0), size=n)
    r0 = row_cost(a, x)
    for d in range(n+1):
        vals = []
        for f in combinations(range(n), d):
            y = x.copy()
            y[list(f)] *= -1
            vals.append(row_cost(a, y))
        mean = sum(vals) / comb(n, d)
        gamma = ((n-2*d)**2-n)/(n*(n-1))
        target = gamma*r0+(1-gamma)*n*(n-1)
        assert abs(mean-target) < 2e-10


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1


def paley_conference(q):
    assert q % 4 == 1
    g = q + 1
    c = np.zeros((g, g), dtype=float)
    c[0, 1:] = 1
    c[1:, 0] = 1
    for i in range(q):
        for j in range(q):
            if i != j:
                c[i+1, j+1] = legendre(i-j, q)
    assert np.array_equal(c, c.T)
    assert np.array_equal(c @ c, q*np.eye(g))
    return c


def check_paley_block():
    for q in (5, 13):
        c = paley_conference(q)
        g = q+1
        j = np.ones((g, g))
        a = np.kron(c, j)+np.kron(np.eye(g), j-np.eye(g))
        n = g*g
        assert np.all(np.diag(a) == 0)
        assert set(np.unique(a)) == {-1.0, 0.0, 1.0}
        op = float(np.linalg.norm(a, 2))
        assert abs(op-(g*np.sqrt(q)+g-1)) < 2e-9
        one = np.ones(n)
        direct = row_cost(a, one)
        formula = g*np.linalg.norm(g*(c@np.ones(g))+(g-1)*np.ones(g))**2
        assert abs(direct-formula) < 2e-7
        assert direct > 0.03*n**2.5
        d = max(1, int(n**0.6))
        gap = np.sqrt(direct)-2*op*np.sqrt(d)
        assert gap > 0
        lower = gap**2
        # The smallest test order has poor finite-size constants; the point
        # of (N7)--(N8) is the positive asymptotic gap for d=o(N).
        assert lower > 1e-5*n**2.5


if __name__ == "__main__":
    check_product_noise()
    check_fixed_spheres()
    check_paley_block()
    print("PASS row_noise_proximity_r32")
