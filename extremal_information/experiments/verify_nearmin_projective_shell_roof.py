#!/usr/bin/env python3
"""Finite regression checks for drafts/nearmin_projective_shell_roof.md.

The proof is exact and does not rely on this script.  We exhaust all hollow
signings through order four, exhaust all order-five signings against a
deterministic query bank, and sample larger orders.  The checks cover:

* maximal projective packings really cover the positive shell;
* every covered atom has the claimed +/- signed local-chart representation;
* the sorting evaluator equals direct subset maximization; and
* the response error never exceeds 2*k*(k-1).
"""

from __future__ import annotations

import itertools
import random


def spins(n: int):
    """One representative of each global-spin class."""
    return [(1,) + tail for tail in itertools.product((-1, 1), repeat=n - 1)]


def edge_list(n: int):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def cut_word(x, edges):
    return tuple(x[i] * x[j] for i, j in edges)


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def d_projective(a, b):
    d = hamming(a, b)
    return min(d, len(a) - d)


def energy(a, c):
    return sum(x * y for x, y in zip(a, c))


def local_fields(a, sigma, u, edges):
    ell = [0] * len(u)
    for ae, (i, j) in zip(a, edges):
        s = sigma * ae * u[i] * u[j]
        ell[i] += s
        ell[j] += s
    return tuple(ell)


def signed_atoms(xs, cuts):
    out = []
    for x, c in zip(xs, cuts):
        for sigma in (-1, 1):
            out.append((sigma, x, tuple(sigma * t for t in c)))
    return out


def maximal_net(shell, radius):
    centres = []
    for atom in shell:
        if all(d_projective(atom[2], c[2]) > radius for c in centres):
            centres.append(atom)
    return centres


def k_for_radius(n, radius):
    return max(d for d in range(n // 2 + 1) if d * (n - d) <= radius)


def atlas_data(a, centres, edges, cuts_by_x):
    out = []
    for sigma, u, _ in centres:
        h = sigma * energy(a, cuts_by_x[u])
        ell = local_fields(a, sigma, u, edges)
        assert sum(ell) == 2 * h
        out.append((sigma, u, h, ell))
    return out


def atlas_sort(data, g, k):
    best = -10**30
    for _, u, h, ell in data:
        gu = sum(gi * ui for gi, ui in zip(g, u))
        for eta in (-1, 1):
            for tau in (-1, 1):
                increments = sorted(
                    (-2 * eta * li - 2 * tau * gi * ui
                     for li, gi, ui in zip(ell, g, u)),
                    reverse=True,
                )
                gain = sum(v for v in increments[:k] if v > 0)
                best = max(best, eta * h + tau * gu + gain)
    return best


def atlas_subsets(data, g, k):
    best = -10**30
    n = len(g)
    for _, u, h, ell in data:
        gu = sum(gi * ui for gi, ui in zip(g, u))
        for eta in (-1, 1):
            for tau in (-1, 1):
                for size in range(k + 1):
                    for s in itertools.combinations(range(n), size):
                        val = eta * h + tau * gu
                        val -= 2 * eta * sum(ell[i] for i in s)
                        val -= 2 * tau * sum(g[i] * u[i] for i in s)
                        best = max(best, val)
    return best


def boolean_response(a, xs, cuts, g):
    return max(
        sigma * energy(a, c) + abs(sum(gi * xi for gi, xi in zip(g, x)))
        for x, c in zip(xs, cuts)
        for sigma in (-1, 1)
    )


def check_instance(a, n, g, G, radius, do_subset_check=False):
    edges = edge_list(n)
    E = len(edges)
    assert radius < E - (n * n // 4)
    xs = spins(n)
    cuts = [cut_word(x, edges) for x in xs]
    cuts_by_x = dict(zip(xs, cuts))
    atoms = signed_atoms(xs, cuts)
    Q = max(abs(energy(a, c)) for c in cuts)
    shell = [atom for atom in atoms if energy(a, atom[2]) >= Q - G]
    centres = maximal_net(shell, radius)
    assert centres
    assert all(
        any(d_projective(z, c[2]) <= radius for c in centres)
        for _, _, z in shell
    )
    assert all(
        d_projective(centres[i][2], centres[j][2]) > radius
        for i in range(len(centres)) for j in range(i)
    )

    k = k_for_radius(n, radius)
    # Explicitly verify the signed-chart lemma for every shell/centre cover.
    for sigma, x, z in shell:
        centre = next(c for c in centres if d_projective(z, c[2]) <= radius)
        sigma0, u, z0 = centre
        eta = 1 if hamming(z, z0) <= hamming(z, tuple(-v for v in z0)) else -1
        assert sigma == eta * sigma0
        d = min(hamming(x, u), n - hamming(x, u))
        assert d <= k

    data = atlas_data(a, centres, edges, cuts_by_x)
    bhat = atlas_sort(data, g, k)
    if do_subset_check:
        assert abs(bhat - atlas_subsets(data, g, k)) <= 1e-9
    b = boolean_response(a, xs, cuts, g)
    assert abs(bhat - b) <= 2 * k * (k - 1) + 1e-9, (
        n, a, g, G, radius, k, len(centres), b, bhat
    )
    return len(centres)


def deterministic_queries(n, G):
    bank = {(0,) * n}
    for i in range(n):
        for s in (-1, 1):
            g = [0] * n
            g[i] = s * min(2, G)
            bank.add(tuple(g))
    # Mixed fields exercise the absolute-field sign split.
    for signs in itertools.product((-1, 1), repeat=min(n, 4)):
        g = list(signs) + [0] * (n - len(signs))
        if sum(abs(v) for v in g) <= G:
            bank.add(tuple(g))
    return sorted(bank)


def main():
    checks = 0
    max_centres = 0

    # Exhaustive matrices and a broad query grid through order four.
    for n in range(2, 5):
        edges = edge_list(n)
        E = len(edges)
        max_radius = E - (n * n // 4) - 1
        for a in itertools.product((-1, 1), repeat=E):
            for G in (0, 2, 4):
                fields = [
                    g for g in itertools.product((-2, -1, 0, 1, 2), repeat=n)
                    if sum(abs(v) for v in g) <= G
                ]
                for radius in range(max_radius + 1):
                    for g in fields:
                        l = check_instance(
                            a, n, g, G, radius,
                            do_subset_check=(checks % 37 == 0),
                        )
                        max_centres = max(max_centres, l)
                        checks += 1

    # Every order-five signing, with deterministic mixed queries.
    n = 5
    edges = edge_list(n)
    E = len(edges)
    for a in itertools.product((-1, 1), repeat=E):
        for G in (0, 3):
            for radius in range(E - (n * n // 4)):
                for g in deterministic_queries(n, G):
                    l = check_instance(a, n, g, G, radius)
                    max_centres = max(max_centres, l)
                    checks += 1

    # Reproducible larger random stress tests, including real-valued fields.
    rng = random.Random(20260817)
    for n in range(6, 10):
        edges = edge_list(n)
        E = len(edges)
        max_radius = E - (n * n // 4) - 1
        for _ in range(80):
            a = tuple(rng.choice((-1, 1)) for _ in edges)
            G = rng.uniform(0.0, 8.0)
            radius = rng.randrange(max_radius + 1)
            raw = [rng.uniform(-1.0, 1.0) for _ in range(n)]
            norm = sum(abs(v) for v in raw)
            scale = 0.0 if norm == 0 else rng.random() * G / norm
            g = tuple(scale * v for v in raw)
            l = check_instance(
                a, n, g, G, radius,
                do_subset_check=(_ == 0),
            )
            max_centres = max(max_centres, l)
            checks += 1

    print(
        "PASS: projective shell roof; "
        f"{checks} matrix/query/radius checks; max centres {max_centres}"
    )


if __name__ == "__main__":
    main()
