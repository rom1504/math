#!/usr/bin/env python3
"""Finite checks for the exact random-coset identities in the Wave 28 memo."""

from itertools import product
import random
import numpy as np


def qform(h, x):
    x = np.asarray(x, dtype=int)
    return int(x @ h @ x)


def check(seed=2801):
    rng = random.Random(seed)
    for n in range(4, 9):
        for trial in range(20):
            a = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    a[i, j] = a[j, i] = rng.choice((-1, 1))
            m = rng.randrange(2, n)
            s = set(rng.sample(range(n), m))
            p2_num, p2_den = m * (m - 1), n * (n - 1)
            # Work with p2_den * H to keep all arithmetic integral.
            h = -p2_num * a
            for i in s:
                for j in s:
                    h[i, j] += p2_den * a[i, j]

            # A random nonempty labeled partition.
            k = rng.randrange(1, min(4, n) + 1)
            labels = list(range(k)) + [rng.randrange(k) for _ in range(n - k)]
            rng.shuffle(labels)
            subgroup = []
            for z in product((-1, 1), repeat=k):
                subgroup.append(tuple(z[labels[i]] for i in range(n)))

            vals = []
            spins = list(product((-1, 1), repeat=n))
            for x in spins:
                vals.append(abs(qform(h, x)))
            threshold = rng.choice(vals)
            e = {x for x, val in zip(spins, vals) if val >= threshold}

            hit_count = 0
            saturation = set()
            for eps in spins:
                hit = False
                for b in subgroup:
                    xb = tuple(eps[i] * b[i] for i in range(n))
                    saturation.add(xb) if xb in e else None
                    if xb in e:
                        hit = True
                hit_count += hit

            # Build E K directly, rather than just the E points seen above.
            ek = {
                tuple(x[i] * b[i] for i in range(n))
                for x in e
                for b in subgroup
            }
            assert hit_count == len(ek)
            assert len(e) <= hit_count <= min(2**n, 2 ** (k - 1) * len(e))

            hf2 = int(np.sum(h * h))
            expected = (
                p2_num * (p2_den - p2_num) * p2_den
            )
            assert hf2 == expected, (n, m, hf2, expected)
    print("PASS: quotient saturation, amplification bounds, and Frobenius identity")


if __name__ == "__main__":
    check()
