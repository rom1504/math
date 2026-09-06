#!/usr/bin/env python3
"""Exact small-order audit of the augmented-cut boundary-state identity."""

from itertools import product
from random import Random


def hamming(a, b):
    return bin(a ^ b).count("1")


def instance(m, ell, exhaustive=False, samples=5000):
    n = m + ell
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    internal_ids = tuple(k for k, (i, j) in enumerate(edges)
                         if (i < m) == (j < m))
    cross_ids = tuple(k for k in range(len(edges)) if k not in internal_ids)

    def restrict(word, ids):
        return sum(((word >> k) & 1) << q for q, k in enumerate(ids))

    # Gauge z_0=0.  The n remaining bits are t,z_1,...,z_{n-1}.
    code = set()
    for t in (0, 1):
        for tail in product((0, 1), repeat=n - 1):
            z = (0,) + tail
            word = sum((t ^ z[i] ^ z[j]) << k
                       for k, (i, j) in enumerate(edges))
            code.add(word)
    assert len(code) == 2 ** n

    ci = {restrict(c, internal_ids) for c in code
          if restrict(c, cross_ids) == 0}
    cx = {restrict(c, cross_ids) for c in code
          if restrict(c, internal_ids) == 0}
    assert len(ci) == len(cx) == 2

    # Projection to X modulo the cross-supported subcode.
    projection = {restrict(c, cross_ids) for c in code}
    unused = set(projection)
    classes = []
    while unused:
        r = min(unused)
        cls = {r ^ k for k in cx}
        assert cls <= projection
        unused -= cls
        classes.append(cls)
    assert len(classes) == 2 ** (n - 2)

    # A lift for each state class, and the corresponding left fiber.
    lifts = []
    fibers = []
    for cls in classes:
        r = min(cls)
        candidates = [c for c in code if restrict(c, cross_ids) == r]
        assert candidates
        lifts.append((restrict(candidates[0], internal_ids), r))
        fibers.append({restrict(c, internal_ids) for c in candidates})

    def W(xi, xx):
        full = sum(((xi >> q) & 1) << k
                   for q, k in enumerate(internal_ids))
        full |= sum(((xx >> q) & 1) << k
                    for q, k in enumerate(cross_ids))
        return min(hamming(full, c) for c in code)

    if exhaustive:
        inputs = product(
            range(1 << len(internal_ids)),
            range(1 << len(cross_ids)),
        )
        audit_label = "exhaustive"
    else:
        rng = Random(20260729 + 100 * m + ell)
        inputs = (
            (
                rng.randrange(1 << len(internal_ids)),
                rng.randrange(1 << len(cross_ids)),
            )
            for _ in range(samples)
        )
        audit_label = f"{samples} sampled entries"

    # Door identity, exhaustive at 3+3 and sampled at 3+4.
    for xi, xx in inputs:
        rhs = min(
            min(hamming(xi, f) for f in fibers[q])
            + min(hamming(xx, r) for r in classes[q])
            for q in range(len(classes))
        )
        assert W(xi, xx) == rhs

    # The transversal block has zero diagonal and positive off diagonal.
    for i, (xi, _) in enumerate(lifts):
        for j, (_, xx) in enumerate(lifts):
            if i == j:
                assert W(xi, xx) == 0
            else:
                assert W(xi, xx) >= 1

    print(
        f"{m}+{ell}: dim C={n}, dim C_I=1, dim C_X=1, "
        f"state dimension={n-2}, states={len(classes)}; "
        f"door ({audit_label})/block PASS"
    )


if __name__ == "__main__":
    instance(3, 3, exhaustive=True)
    instance(3, 4, samples=5000)
