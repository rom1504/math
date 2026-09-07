#!/usr/bin/env python3
"""Independent pure-Python audit of the Wave 12 replacement certificates."""

from itertools import product


A7 = (
    (0, 1, 1, 1, 1, 1, 1),
    (1, 0, 1, 1, -1, -1, 1),
    (1, 1, 0, 1, -1, 1, -1),
    (1, 1, 1, 0, 1, -1, -1),
    (1, -1, -1, 1, 0, -1, -1),
    (1, -1, 1, -1, -1, 0, -1),
    (1, 1, -1, -1, -1, -1, 0),
)


def spin_vectors(n):
    return ((1,) + tail for tail in product((-1, 1), repeat=n - 1))


def all_spin_vectors(n):
    return product((-1, 1), repeat=n)


def energy(a, x):
    return 2 * sum(a[i][j] * x[i] * x[j]
                   for i in range(len(a)) for j in range(i + 1, len(a)))


def q(a):
    return max(abs(energy(a, x)) for x in spin_vectors(len(a)))


def gauge_matrix(n, bits):
    a = [[0] * n for _ in range(n)]
    for j in range(1, n):
        a[0][j] = a[j][0] = 1
    pos = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = 1 if (bits >> pos) & 1 else -1
            pos += 1
    return tuple(tuple(row) for row in a)


def gauge_minima(n):
    best = n * (n - 1)
    mins = []
    for bits in range(1 << ((n - 1) * (n - 2) // 2)):
        a = gauge_matrix(n, bits)
        qa = q(a)
        if qa < best:
            best, mins = qa, [a]
        elif qa == best:
            mins.append(a)
    return best, mins


def switch(a, s):
    return tuple(tuple(a[i][j] * s[i] * s[j] for j in range(len(a)))
                 for i in range(len(a)))


def all_minima(n):
    qn, gauges = gauge_minima(n)
    out = {}
    for a in gauges:
        for tail in product((-1, 1), repeat=n - 1):
            s = (1,) + tail
            b = switch(a, s)
            out[b] = b
    return qn, tuple(out.values()), len(gauges)


def principal(a, vertices):
    return tuple(tuple(a[i][j] for j in vertices) for i in vertices)


def cross(a, heads, tails):
    return tuple(tuple(a[i][j] for j in tails) for i in heads)


def layer(c, d):
    qd = q(d)
    value = 0
    for y in spin_vectors(len(d)):
        exposure = 2 * sum(abs(sum(row[j] * y[j] for j in range(len(y))))
                           for row in c)
        value = max(value, exposure - qd + abs(energy(d, y)))
    return value


def block_matrix(ds, crosses):
    sizes = [len(d) for d in ds]
    starts = [0]
    for size in sizes:
        starts.append(starts[-1] + size)
    n = starts[-1]
    a = [[0] * n for _ in range(n)]
    for b, d in enumerate(ds):
        off = starts[b]
        for i in range(len(d)):
            for j in range(len(d)):
                a[off + i][off + j] = d[i][j]
    for (b, c), block in crosses.items():
        for i in range(sizes[b]):
            for j in range(sizes[c]):
                a[starts[b] + i][starts[c] + j] = block[i][j]
                a[starts[c] + j][starts[b] + i] = block[i][j]
    return tuple(tuple(row) for row in a)


def direct_phi(ds, crosses):
    a = block_matrix(ds, crosses)
    baseline = sum(q(d) for d in ds)
    # Algebraic definition and direct common-mosaic maximization.
    algebraic = q(a) - baseline
    sizes = [len(d) for d in ds]
    best = -10**9
    for local in product(*[tuple(all_spin_vectors(s)) for s in sizes]):
        for sigma in (-1, 1):
            cross_term = 0
            for (b, c), block in crosses.items():
                cross_term += 2 * sigma * sum(
                    block[i][j] * local[b][i] * local[c][j]
                    for i in range(sizes[b]) for j in range(sizes[c])
                )
            deficits = sum(q(ds[b]) - sigma * energy(ds[b], local[b])
                           for b in range(len(ds)))
            best = max(best, cross_term - deficits)
    assert best == algebraic
    return algebraic, q(a)


def main():
    q5, mins5, g5 = all_minima(5)
    q6, mins6, g6 = all_minima(6)
    q7, gauges7 = gauge_minima(7)
    g7 = len(gauges7)
    assert (q5, q6, q7) == (8, 10, 18)
    assert (g5, g6, g7) == (12, 12, 3240)
    assert q(A7) == 18

    tails6 = (0, 1, 2, 4, 5, 6)
    heads6 = (3,)
    d6, c6 = principal(A7, tails6), cross(A7, heads6, tails6)
    assert (q(d6), layer(c6, d6), min(layer(c6, g) for g in mins6)) == (18, 0, 8)

    tails5 = (0, 1, 2, 3, 4)
    heads5 = (5, 6)
    d5, h2 = principal(A7, tails5), principal(A7, heads5)
    c5 = cross(A7, heads5, tails5)
    assert (q(d5), layer(c5, d5), min(layer(c5, g) for g in mins5)) == (12, 4, 12)

    # R12.18 on the displayed 2+5 split and one minimizing replacement.
    crosses = {(0, 1): tuple(tuple(c5[j][i] for j in range(2)) for i in range(5))}
    phi_d, qa = direct_phi((d5, h2), crosses)
    phi_g, qtilde = direct_phi((mins5[0], h2), crosses)
    assert qa == q7
    assert phi_g - phi_d == (q(d5) - q5) + (q(h2) - q(h2)) + (qtilde - q7)

    print("independent q/minimizer counts and A7 boundary claims: PASS")
    print("common-mosaic signed response identity: PASS", phi_d, phi_g, qtilde)


if __name__ == "__main__":
    main()
