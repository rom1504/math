#!/usr/bin/env python3
"""Exact audit of liberal one-path harvesting for the (10.409) matrix."""

from functools import lru_cache
from itertools import product

A = (
    (0, -1, 1, -1, -1, -1),
    (-1, 0, 1, -1, 1, 1),
    (1, 1, 0, -1, 1, -1),
    (-1, -1, -1, 0, 1, -1),
    (-1, 1, 1, 1, 0, -1),
    (-1, 1, -1, -1, -1, 0),
)


def energy(vertices, spin):
    return sum(
        A[i][j] * spin[r] * spin[s]
        for r, i in enumerate(vertices)
        for s, j in enumerate(vertices)
        if r != s
    )


@lru_cache(None)
def endpoints(vertices):
    spins = tuple(product((-1, 1), repeat=len(vertices)))
    vals = tuple(energy(vertices, x) for x in spins)
    hi, lo = max(vals), min(vals)
    ps = tuple(x for x, v in zip(spins, vals) if v == hi)
    ns = tuple(x for x, v in zip(spins, vals) if v == lo)
    return hi, -lo, ps, ns


def cross_l1(vertices, p, child_positions):
    child = set(child_positions)
    other = [j for j in range(len(vertices)) if j not in child]
    return sum(
        abs(sum(A[vertices[i]][vertices[j]] * p[i] for i in child))
        for j in other
    )


@lru_cache(None)
def best_path(vertices):
    """Maximum sum of positive mu along one path; optimize every endpoint tie."""
    if len(vertices) <= 1:
        return 0, ()
    P, N, ps, ns = endpoints(vertices)
    Q = max(P, N)
    best = (-10**9, ())
    for p in ps:
        for n in ns:
            sides = {}
            for bit in (-1, 1):
                positions = tuple(i for i in range(len(vertices)) if p[i] * n[i] == bit)
                if positions:
                    sides[positions] = None
            if len(sides) != 2:
                continue
            for positions in sides:
                child_vertices = tuple(vertices[i] for i in positions)
                h = energy(child_vertices, tuple(p[i] for i in positions))
                L = cross_l1(vertices, p, positions)
                child_value, child_trace = best_path(child_vertices)
                for sigma in (-1, 1):
                    mu = 2 * L - (max(endpoints(child_vertices)[:2]) - sigma * h)
                    value = max(mu, 0) + child_value
                    trace = ((vertices, child_vertices, sigma, mu, L, h),) + child_trace
                    if value > best[0]:
                        best = (value, trace)
    return best


root = tuple(range(len(A)))
P, N, _, _ = endpoints(root)
value, trace = best_path(root)
print("P,N,R,best_path", P, N, P + N, value)
for row in trace:
    print(row)

# The positive triangle is already an asymmetric, fixed-orientation version:
# P=6>N=2, Range=8, but its most liberal one-path total is 6.
A = (
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0),
)
endpoints.cache_clear()
best_path.cache_clear()
root = tuple(range(len(A)))
P, N, _, ns = endpoints(root)
value, trace = best_path(root)
negative_endpoint = ns[0]
print("triangle P,N,R,best_path", P, N, P + N, value, "n", negative_endpoint)
for row in trace:
    print(row)
for k in range(1, 8):
    vals = []
    for dspin in product((-1, 1), repeat=k):
        m = sum(dspin)
        for y in product((-1, 1), repeat=len(root)):
            vals.append(m * m - k + energy(root, y) + 2 * m * sum(
                a * b for a, b in zip(negative_endpoint, y)
            ))
    q = max(max(vals), -min(vals))
    target = k * k - k - N + 2 * k * len(root)
    print("triangle parent", k, "Q", q, "target", target, "ground", target == q)

# Embed the negative endpoint as the positive restriction of a larger parent:
# a positive clique D of size k, with every D-to-root row equal to n.
negative_endpoint = endpoints(root)[3][0]
for k in range(1, 9):
    vals = []
    witnesses = []
    for dspin in product((-1, 1), repeat=k):
        m = sum(dspin)
        for y in product((-1, 1), repeat=len(root)):
            val = m * m - k + energy(root, y) + 2 * m * sum(
                a * b for a, b in zip(negative_endpoint, y)
            )
            vals.append(val)
            witnesses.append((dspin, y, val))
    q = max(max(vals), -min(vals))
    target = k * k - k - N + 2 * k * len(root)
    target_is_positive_ground = target == q
    print("parent", k, "Q", q, "target", target, "ground", target_is_positive_ground)

# A random order-eight endpoint cut where both child restriction energies lie
# strictly between their asymmetric child-range midpoints.
A = (
    (0, 1, -1, 1, -1, -1, 1, 1),
    (1, 0, 1, 1, -1, -1, 1, 1),
    (-1, 1, 0, -1, 1, -1, 1, 1),
    (1, 1, -1, 0, 1, -1, -1, -1),
    (-1, -1, 1, 1, 0, 1, 1, 1),
    (-1, -1, -1, -1, 1, 0, -1, 1),
    (1, 1, 1, -1, 1, -1, 0, 1),
    (1, 1, 1, -1, 1, 1, 1, 0),
)
endpoints.cache_clear()
best_path.cache_clear()
root = tuple(range(len(A)))
P, N, _, _ = endpoints(root)
value, trace = best_path(root)
print("order8 P,N,R,best_path", P, N, P + N, value)
for row in trace:
    print(row)
