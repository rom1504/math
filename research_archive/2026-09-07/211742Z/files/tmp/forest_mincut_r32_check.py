#!/usr/bin/env python3
"""Independent finite checks for Wave 32's completion-forest formula."""

from itertools import combinations, product
import random


def dh(x, y):
    return sum(a != b for a, b in zip(x, y))


def dpr(x, y):
    d = dh(x, y)
    return min(d, len(x) - d)


def adjacency(num_vertices, edges):
    adj = [[] for _ in range(num_vertices)]
    for e, (u, v) in enumerate(edges):
        adj[u].append((v, e))
        adj[v].append((u, e))
    return adj


def coordinate_dp(edges, pins):
    """Minimum tree total variation with binary pins; tree is rooted at 0."""
    nvert = len(edges) + 1
    adj = adjacency(nvert, edges)
    inf = 10**9

    def rec(v, parent):
        child_tables = [rec(u, v) for u, _ in adj[v] if u != parent]
        out = {}
        for a in (-1, 1):
            if v in pins and pins[v] != a:
                out[a] = inf
            else:
                out[a] = sum(min(tab[c] + (a != c) for c in (-1, 1))
                             for tab in child_tables)
        return out

    tab = rec(0, -1)
    return min(tab.values())


def coordinate_separator(edges, pins):
    """Brute minimum edge set separating all oppositely pinned terminals."""
    nvert = len(edges) + 1
    best = len(edges) + 1
    for mask in range(1 << len(edges)):
        used = bin(mask).count("1")
        if used >= best:
            continue
        parent = list(range(nvert))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            a, b = find(a), find(b)
            if a != b:
                parent[b] = a

        for e, (u, v) in enumerate(edges):
            if not (mask >> e) & 1:
                union(u, v)
        seen = {}
        valid = True
        for v, sign in pins.items():
            comp = find(v)
            if comp in seen and seen[comp] != sign:
                valid = False
                break
            seen[comp] = sign
        if valid:
            best = used
    return best


def tree_path_mask(edges, start, target):
    adj = adjacency(len(edges) + 1, edges)
    stack = [(start, -1, 0)]
    while stack:
        v, parent, mask = stack.pop()
        if v == target:
            return mask
        for u, e in adj[v]:
            if u != parent:
                stack.append((u, v, mask | (1 << e)))
    raise AssertionError("tree path not found")


def coordinate_path_packing(edges, pins):
    """Brute maximum number of edge-disjoint opposite-terminal paths."""
    plus = [v for v, a in pins.items() if a == 1]
    minus = [v for v, a in pins.items() if a == -1]
    paths = [tree_path_mask(edges, u, v) for u in plus for v in minus]
    best = 0
    for choice in range(1 << len(paths)):
        if bin(choice).count("1") <= best:
            continue
        occupied = 0
        valid = True
        for j, path in enumerate(paths):
            if (choice >> j) & 1:
                if occupied & path:
                    valid = False
                    break
                occupied |= path
        if valid:
            best = bin(choice).count("1")
    return best


def formula_cost(ncoord, edges, z, masks, words):
    q = len(masks)
    best = 10**9
    for sigma in product((-1, 1), repeat=q):
        total = 0
        for i in range(ncoord):
            pins = {0: z[i]}
            for v in range(q):
                if i in masks[v]:
                    pins[v + 1] = sigma[v] * words[v][i]
            d = coordinate_dp(edges, pins)
            assert d == coordinate_separator(edges, pins)
            assert d == coordinate_path_packing(edges, pins)
            total += d
        best = min(best, total)
    return best


def projective_candidates(ncoord, mask, word):
    ans = []
    for x in product((-1, 1), repeat=ncoord):
        restriction = tuple(x[i] for i in sorted(mask))
        target = tuple(word[i] for i in sorted(mask))
        if restriction == target or restriction == tuple(-a for a in target):
            ans.append(x)
    return ans


def brute_projective_cost(ncoord, edges, z, masks, words):
    candidates = [projective_candidates(ncoord, s, y)
                  for s, y in zip(masks, words)]
    best = 10**9
    for chosen in product(*candidates):
        labels = (tuple(z),) + chosen
        cost = sum(dpr(labels[u], labels[v]) for u, v in edges)
        best = min(best, cost)
    return best


def random_formula_checks():
    rng = random.Random(3201)
    for _ in range(160):
        ncoord = rng.randint(1, 4)
        q = rng.randint(1, 3)
        # The index order makes this a rooted random recursive tree.
        edges = [(v, rng.randrange(v)) for v in range(1, q + 1)]
        z = tuple(rng.choice((-1, 1)) for _ in range(ncoord))
        masks = []
        words = []
        for _v in range(q):
            mask = {i for i in range(ncoord) if rng.randrange(2)}
            if not mask:
                mask = {rng.randrange(ncoord)}
            word = tuple(rng.choice((-1, 1)) for _ in range(ncoord))
            masks.append(mask)
            words.append(word)
        exact = formula_cost(ncoord, edges, z, masks, words)
        brute = brute_projective_cost(ncoord, edges, z, masks, words)
        assert exact == brute, (ncoord, edges, z, masks, words, exact, brute)


def orientation_and_hidden_conflict_checks():
    # A node orientation is shared across coordinates: the answer is 1, not 0.
    ncoord = 2
    edges = [(0, 1)]
    z = (1, 1)
    masks = [{0, 1}]
    words = [(1, -1)]
    assert formula_cost(ncoord, edges, z, masks, words) == 1
    false_coordinatewise = sum(
        min(coordinate_dp(edges, {0: z[i], 1: s * words[0][i]})
            for s in (-1, 1))
        for i in range(ncoord)
    )
    assert false_coordinatewise == 0

    # Adjacent selector-overlap checks see no conflict through an unpinned node.
    edges = [(0, 1), (1, 2)]
    pins = {0: 1, 2: -1}
    assert coordinate_dp(edges, pins) == 1
    assert coordinate_separator(edges, pins) == 1
    assert coordinate_path_packing(edges, pins) == 1


def weighted_gadget_check():
    n = 5
    A = [[0] * n for _ in range(n)]
    weights = {
        (0, 1): 10, (0, 2): 10, (1, 2): 20, (3, 4): 20,
        (0, 3): 1, (0, 4): 1, (1, 3): -1, (1, 4): -1,
    }
    for (i, j), a in weights.items():
        A[i][j] = A[j][i] = a

    def matvec(x):
        return tuple(sum(A[i][j] * x[j] for j in range(n)) for i in range(n))

    def energy(x):
        ax = matvec(x)
        return sum(x[i] * ax[i] for i in range(n))

    def row_square(x):
        ax = matvec(x)
        return sum(a * a for a in ax)

    all_spins = list(product((-1, 1), repeat=n))
    maximum = max(map(energy, all_spins))
    grounds = [x for x in all_spins if energy(x) == maximum]
    x = (1, 1, 1, 1, 1)
    F = {0, 1, 2}
    xf = tuple(-a if i in F else a for i, a in enumerate(x))
    assert maximum == 120
    assert x in grounds and xf in grounds
    assert len(grounds) == 4  # two projective ground classes
    assert dpr(x, xf) == 2
    assert row_square(x) == 2968
    assert row_square(xf) == 3048
    boundary = tuple(sum(A[i][j] for j in range(n)
                         if (i in F) != (j in F)) for i in range(n))
    assert boundary == (2, -2, 0, 0, 0)
    assert energy(x) - energy(xf) == 0


if __name__ == "__main__":
    random_formula_checks()
    orientation_and_hidden_conflict_checks()
    weighted_gadget_check()
    print("forest_mincut_r32_check: PASS")
