#!/usr/bin/env python3
"""Exact tiny checks for kernel partitions and paired selector cycles."""

from fractions import Fraction as Q
from itertools import product


def compose(rho, sigma):
    return tuple(rho[sigma[i]] for i in range(len(rho)))


def kernel(f):
    return tuple(tuple(f[i] == f[j] for j in range(len(f)))
                 for i in range(len(f)))


def pullback(partition, sigma):
    return tuple(tuple(partition[sigma[i]][sigma[j]]
                       for j in range(len(sigma)))
                 for i in range(len(sigma)))


def clip(z, lo=Q(0), hi=Q(1)):
    return min(max(z, lo), hi)


def main():
    checked = 0
    for r in (2, 3):
        maps = tuple(product(range(r), repeat=r))
        ident = tuple(range(r))
        for depth in range(1, 5):
            for word in product(maps, repeat=depth):
                rho = ident
                partitions = kernel(rho)
                has_reset = False
                for sigma in word:
                    rho = compose(rho, sigma)
                    partitions = pullback(partitions, sigma)
                    assert partitions == kernel(rho)
                    has_reset |= len(set(rho)) == 1
                # Constants are a two-sided ideal: once reset, total is reset.
                assert has_reset == (len(set(rho)) == 1)
                checked += 1

    # Different selectors force the cross carrier.
    for z in map(Q, range(-5, 6)):
        x = (Q(0), z)
        y = x
        sigma = (0, 1)
        tau = (1, 0)
        old_diag = tuple(y[i] - x[i] for i in range(2))
        new_diag = tuple(y[tau[i]] - x[sigma[i]] for i in range(2))
        assert old_diag == (0, 0)
        assert new_diag == (z, -z)
        checked += 1

    # A local middle-cell self-loop is not pumpable.
    delta = Q(1, 4)
    for seed in (Q(1, 2), Q(3, 4), Q(5, 4)):
        z = seed
        middle_steps = 0
        for _ in range(12):
            if delta < z < 1 + delta:
                middle_steps += 1
            z = clip(z - delta)
        assert middle_steps <= 4
        checked += 1

    # Two bounded involutions can make a drifting mixed word.
    a = lambda z: -z
    b = lambda z: 1 - z
    for z in map(Q, range(-3, 4)):
        assert a(a(z)) == z and b(b(z)) == z
        assert b(a(z)) == z + 1
        checked += 1

    print(f"selector-partition and paired-cycle checks: {checked}")


if __name__ == "__main__":
    main()
