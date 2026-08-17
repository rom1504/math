#!/usr/bin/env python3
"""Exact width-two Ising and weighted-automaton switching checks."""

from fractions import Fraction as Q
from itertools import product
from random import Random


SPINS = tuple(product((-1, 1), repeat=2))


def hilbert(v, w):
    d = [x - y for x, y in zip(v, w)]
    return (max(d) - min(d)) / 2


def normalize(v):
    m = max(v)
    return tuple(x - m for x in v)


def transfer(v, fields, vertical, horizontal):
    out = []
    for y in SPINS:
        local = fields[0] * y[0] + fields[1] * y[1] + vertical * y[0] * y[1]
        out.append(local + max(
            v[k] + horizontal[0] * x[0] * y[0]
                 + horizontal[1] * x[1] * y[1]
            for k, x in enumerate(SPINS)
        ))
    return normalize(tuple(out))


def nearest_multiple(x, mesh):
    z = x / mesh
    q = z.numerator // z.denominator
    if 2 * (z - q) >= 1:
        q += 1
    return q * mesh


def quantize(v, mesh):
    return normalize(tuple(nearest_multiple(x, mesh) for x in normalize(v)))


def ising_checks():
    rng = Random(20260817)
    checks = 0
    for _ in range(1200):
        v = normalize(tuple(Q(rng.randrange(-16, 17), 8) for _ in SPINS))
        w = normalize(tuple(Q(rng.randrange(-16, 17), 8) for _ in SPINS))
        j = (Q(rng.randrange(-8, 9), 16), Q(rng.randrange(-8, 9), 16))
        fields = (Q(rng.randrange(-4, 5), 8), Q(rng.randrange(-4, 5), 8))
        vertical = Q(rng.randrange(-4, 5), 8)
        fv = transfer(v, fields, vertical, j)
        fw = transfer(w, fields, vertical, j)
        assert hilbert(fv, fw) <= 2 * sum(abs(x) for x in j)
        checks += 1

    mesh = Q(1, 64)
    eta = mesh / 2
    gap = 4
    weak = (Q(1, 32), Q(-1, 32))
    reset_diameter = 2 * sum(abs(x) for x in weak)
    bound = reset_diameter + gap * eta
    for _ in range(120):
        exact = normalize(tuple(Q(rng.randrange(-16, 17), 8) for _ in SPINS))
        approx = quantize(exact, mesh)
        for t in range(48):
            horizontal = weak if t % gap == 0 else (
                Q(rng.randrange(-8, 9), 16), Q(rng.randrange(-8, 9), 16)
            )
            fields = (Q(rng.randrange(-4, 5), 8), Q(rng.randrange(-4, 5), 8))
            vertical = Q(rng.randrange(-4, 5), 8)
            exact = transfer(exact, fields, vertical, horizontal)
            approx = quantize(transfer(approx, fields, vertical, horizontal), mesh)
            assert hilbert(exact, approx) <= bound
            checks += 1
    return checks


NEG = None


def block_max(row, block):
    vals = [row[j] for j in block if row[j] is not NEG]
    return max(vals) if vals else NEG


def refine_partition(partition, matrices, terminal):
    while True:
        signature = {}
        for old_id, block in enumerate(partition):
            for i in block:
                signature[i] = (terminal[i], old_id, tuple(
                    block_max(matrix[i], target)
                    for matrix in matrices for target in partition
                ))
        groups = {}
        for i, sig in signature.items():
            groups.setdefault(sig, []).append(i)
        new = tuple(tuple(v) for _, v in sorted(groups.items(), key=lambda z: min(z[1])))
        if new == partition:
            return new
        partition = new


def maxplus_mv(matrix, v):
    out = []
    for i, row in enumerate(matrix):
        vals = [row[j] + v[j] for j in range(len(v)) if row[j] is not NEG]
        out.append(max(vals))
    return tuple(out)


def automaton_checks():
    # The quotient {0,1}|{2,3} is discovered from terminal and block signatures.
    matrix = (
        (Q(0), Q(-2), Q(1), Q(0)),
        (Q(-1), Q(0), Q(0), Q(1)),
        (Q(2), Q(0), Q(0), Q(-2)),
        (Q(0), Q(2), Q(-1), Q(0)),
    )
    terminal = (Q(0), Q(0), Q(1), Q(1))
    found = refine_partition(((0, 1), (2, 3)), (matrix,), terminal)
    assert found == ((0, 1), (2, 3))

    # A forced one-block quotient fails by a pumpable projective toll.
    delta = Q(1, 7)
    defective = ((Q(0), NEG), (NEG, delta))
    split = refine_partition(((0, 1),), (defective,), (Q(0), Q(0)))
    assert split == ((0,), (1,))
    v = (Q(0), Q(0))
    for t in range(1, 31):
        v = maxplus_mv(defective, v)
        assert (max(v) - min(v)) / 2 == t * delta / 2
    return 32


def main():
    print(f"width-two Ising switching checks: {ising_checks()}")
    print(f"weighted-automaton lumpability checks: {automaton_checks()}")


if __name__ == "__main__":
    main()
