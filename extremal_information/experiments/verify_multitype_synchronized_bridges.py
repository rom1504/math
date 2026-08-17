#!/usr/bin/env python3
"""Exact small checks of Theorem MT.1 and its typewise holonomy witness."""

from itertools import product
import json


def spins(n):
    return list(product((-1, 1), repeat=n))


def counts(x, cells):
    return tuple(sum(x[i] == 1 for i in cell) for cell in cells)


def micro_energy(xs, cells, alpha, beta, hs):
    out = sum(hs[a][counts(xs[a], cells)] for a in range(len(xs)))
    for a in range(len(xs)):
        for b in range(a + 1, len(xs)):
            for c, cell in enumerate(cells):
                out += alpha[a, b, c] * sum(xs[a][i] * xs[b][i] for i in cell)
            sa = [sum(xs[a][i] for i in cell) for cell in cells]
            sb = [sum(xs[b][i] for i in cell) for cell in cells]
            for c in range(len(cells)):
                for d in range(len(cells)):
                    out += beta[a, b, c, d] * sa[c] * sb[d]
    return out


def quotient_energy(ks, cells, alpha, beta, hs):
    out = sum(hs[a][ks[a]] for a in range(len(ks)))
    ns = [len(cell) for cell in cells]
    for a in range(len(ks)):
        for b in range(a + 1, len(ks)):
            out += sum(
                alpha[a, b, c] * (ns[c] - 2 * abs(ks[a][c] - ks[b][c]))
                for c in range(len(cells))
            )
            for c in range(len(cells)):
                for d in range(len(cells)):
                    sa = 2 * ks[a][c] - ns[c]
                    sb = 2 * ks[b][d] - ns[d]
                    out += beta[a, b, c, d] * sa * sb
    return out


def check_instance(n, cells, blocks):
    alpha = {}
    beta = {}
    for a in range(blocks):
        for b in range(a + 1, blocks):
            for c in range(len(cells)):
                alpha[a, b, c] = 1 + ((a + 2 * b + c) % 3)
            for c in range(len(cells)):
                for d in range(len(cells)):
                    beta[a, b, c, d] = ((a + b + 2 * c + d) % 3) - 1
    count_space = list(product(*(range(len(cell) + 1) for cell in cells)))
    hs = []
    for a in range(blocks):
        hs.append({k: ((a + 1) * sum((c + 1) * v for c, v in enumerate(k))) % 5 - 2
                   for k in count_space})

    raw = spins(n)
    micro = max(
        micro_energy(xs, cells, alpha, beta, hs)
        for xs in product(raw, repeat=blocks)
    )
    quot = max(
        quotient_energy(ks, cells, alpha, beta, hs)
        for ks in product(count_space, repeat=blocks)
    )
    assert micro == quot
    return micro


def check_unbalanced_triangle(cell_size):
    assert cell_size % 2 == 0
    raw = [x for x in spins(cell_size) if sum(x) == 0]
    true_opt = max(
        sum(x[i] * y[i] + y[i] * z[i] - x[i] * z[i]
            for i in range(cell_size))
        for x, y, z in product(raw, repeat=3)
    )
    assert true_opt == cell_size
    assert 3 * cell_size - true_opt == 2 * cell_size


def main():
    checks = []
    for n, cells, blocks in [
        (3, ((0,), (1, 2)), 2),
        (4, ((0, 1), (2, 3)), 2),
        (4, ((0,), (1,), (2, 3)), 2),
        (3, ((0,), (1, 2)), 3),
    ]:
        checks.append({"n": n, "types": len(cells), "blocks": blocks,
                       "optimum": check_instance(n, cells, blocks)})
    for size in (2, 4):
        check_unbalanced_triangle(size)
    print(json.dumps({
        "exact_micro_vs_multitype_quotient": checks,
        "unbalanced_type_cycles_checked": 2,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
