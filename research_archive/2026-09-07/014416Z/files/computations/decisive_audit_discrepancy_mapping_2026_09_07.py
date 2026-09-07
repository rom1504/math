"""Exact finite check of induced and balanced graph-discrepancy mappings."""
from fractions import Fraction
from itertools import combinations, product
import random


def check(n, edges):
    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    total = sum(edges)
    egraph = sum(a == 1 for a in edges)
    degree = [0] * n
    for (i, j), a in zip(pairs, edges):
        degree[i] += a
        degree[j] += a
    ninduced = nbalanced = 0
    for x in product((-1, 1), repeat=n):
        size = sum(v == 1 for v in x)
        energy = sum(a*x[i]*x[j] for (i, j), a in zip(pairs, edges))
        ein = sum(a == 1 and x[i] == x[j] == 1
                  for (i, j), a in zip(pairs, edges))
        d = ein-Fraction(egraph, m)*(size*(size-1)//2)
        rhs = Fraction(energy+sum(a*b for a, b in zip(degree, x))+total, 8)
        rhs -= Fraction(total, 2*m)*(size*(size-1)//2)
        assert d == rhs
        ninduced += 1
        if size == n//2:
            h = n*(n-2)//4
            overlap = sum(a == 1 and x[i] == x[j]
                          for (i, j), a in zip(pairs, edges))
            discrepancy = overlap-Fraction(egraph*h, m)
            assert discrepancy == Fraction(energy, 4)+Fraction(total, 4*(n-1))
            nbalanced += 1
    return ninduced, nbalanced


def main():
    rng = random.Random(20260907)
    ni = nb = matrices = 0
    for n in (4, 6, 8):
        m = n*(n-1)//2
        rows = product((-1, 1), repeat=m) if n == 4 else (
            [rng.choice((-1, 1)) for _ in range(m)] for _ in range(40))
        for edges in rows:
            a, b = check(n, edges)
            ni += a
            nb += b
            matrices += 1
    print("PASS: %d actual signings, %d induced identities, %d balanced identities"
          % (matrices, ni, nb))


if __name__ == "__main__":
    main()
