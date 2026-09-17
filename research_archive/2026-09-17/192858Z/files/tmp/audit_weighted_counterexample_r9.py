#!/usr/bin/env python3
"""Exact audit of the 5+4 weighted same-positive counterexample."""

from fractions import Fraction as F
from itertools import product


U = tuple(range(3))
V = tuple(range(3, 5))
C = tuple(range(5, 8))
D = (8,)
BLOCKS = (U, V, C, D)


def put(w, A, B, value, same=False):
    for i in A:
        for j in B:
            if i < j and (not same or i != j):
                w[i, j] = value


def weights():
    w = {}
    put(w, U, U, F(1, 8), True)
    put(w, U, V, F(-1, 12))
    put(w, V, V, F(1, 6), True)
    put(w, C, C, F(1, 8), True)
    put(w, C, D, F(1, 24))
    put(w, U, C, F(1, 24))
    put(w, U, D, F(1, 24))
    put(w, V, C, F(1, 12))
    put(w, V, D, F(0))
    assert len(w) == 36
    return w


def strict_weights():
    """Mix 49/50 of the witness with 1/50 of K_{5,4}/20."""
    w = {e: F(49, 50) * a for e, a in weights().items()}
    for i in U + V:
        for j in C + D:
            w[i, j] += F(1, 1000)
    return w


def cut(w, z, vertices=tuple(range(9))):
    X = set(vertices)
    return sum(a for (i, j), a in w.items()
               if i in X and j in X and z[i] != z[j])


def main():
    for label, w in (("raw", weights()), ("strict", strict_weights())):
      print("--", label)
      all_states = [(0,) + z for z in product((0, 1), repeat=8)]
      vals = [cut(w, z) for z in all_states]
      S, T = U + V, C + D
      s_states = [z + (0,) * 4 for z in product((0, 1), repeat=5)]
      t_states = [(0,) * 5 + z for z in product((0, 1), repeat=4)]
      svals = [cut(w, z, S) for z in s_states]
      tvals = [cut(w, z, T) for z in t_states]
      aS = sum(a for (i, j), a in w.items() if i in S and j in S)
      aT = sum(a for (i, j), a in w.items() if i in T and j in T)
      c = sum(a for (i, j), a in w.items() if i in S and j in T)
      nontriv = [v for v, z in zip(vals, all_states) if 0 < sum(z) < 9]
      endpoint = tuple(int(i in T) for i in range(9))
      nonend = [v for v, z in zip(vals, all_states)
                if z not in (endpoint, tuple(1-b for b in endpoint))]
      print("global min/max", min(vals), max(vals), "counts", vals.count(min(vals)), vals.count(max(vals)))
      print("strict lower/upper gaps", min(nontriv), c-max(nonend))
      print("S a,m,M,dominance rhs", aS, min(svals), max(svals), aS-min(svals))
      print("T a,m,M,dominance rhs", aT, min(tvals), max(tvals), aT-min(tvals))
      print("c, target lhs, violation", c, aS+aT-min(svals)-min(tvals), aS+aT-min(svals)-min(tvals)-c)
      if label == "raw":
          print("full value histogram", sorted((v, vals.count(v)) for v in set(vals)))


if __name__ == "__main__":
    main()
