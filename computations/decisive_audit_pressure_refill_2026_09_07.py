"""Finite exact-state checks of the pressure refill inequality (zero noise).

The finite signing average is exhaustive. Floating arithmetic is only a
diagnostic for log-sum-exp; the analytic proof is in the paired artifact.
"""
from itertools import combinations, product
import math
import random
import numpy as np


def pressure(a, beta):
    n = len(a)
    x = np.asarray(list(product((-1, 1), repeat=n)), dtype=float)
    h = np.einsum("ij,jk,ik->i", x, a, x)/2
    z = np.concatenate((h, -h))*beta/math.sqrt(n)
    peak = float(z.max())
    return (peak+math.log(float(np.exp(z-peak).sum())))/(beta*n)


def main():
    rng = random.Random(20260907)
    cases = refills = 0
    worst_ratio = 0.0
    for n, deleted in ((3, 1), (4, 1), (4, 2), (5, 1), (6, 1)):
        pairs = list(combinations(range(n), 2))
        missing = [(i, j) for i, j in pairs if j >= n-deleted]
        for beta in (0.2, 0.7, 1.3):
            a = np.zeros((n, n))
            for i, j in pairs:
                a[i, j] = a[j, i] = rng.choice((-1, 1))
            a0 = a.copy()
            for i, j in missing:
                a0[i, j] = a0[j, i] = 0
            base = pressure(a0, beta)
            assert pressure(a, beta) >= base-1e-12
            vals = []
            for signs in product((-1, 1), repeat=len(missing)):
                b = a0.copy()
                for (i, j), sign in zip(missing, signs):
                    b[i, j] = b[j, i] = sign
                value = pressure(b, beta)
                assert value >= base-1e-12
                vals.append(value)
                refills += 1
            increase = sum(vals)/len(vals)-base
            bound = len(missing)*math.log(math.cosh(beta/math.sqrt(n)))/(beta*n)
            assert increase <= bound+1e-12
            assert bound <= beta*len(missing)/(2*n*n)+1e-12
            worst_ratio = max(worst_ratio, increase/bound)
            cases += 1
    print("PASS: %d parameter/core cases, %d exhaustive refills; max ratio %.12f"
          % (cases, refills, worst_ratio))


if __name__ == "__main__":
    main()
