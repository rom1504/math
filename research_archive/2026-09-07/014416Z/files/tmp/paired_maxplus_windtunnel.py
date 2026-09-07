#!/usr/bin/env python3
"""Exact tiny wind tunnel for paired 2-state max-plus maps.

Matrices use rows as input coordinates and columns as output coordinates.
None denotes -infinity.  All arithmetic is rational.
"""

from fractions import Fraction as Q
from itertools import product


def q(x):
    return x if isinstance(x, Q) else Q(x)


def apply_matrix(S, z):
    """Projective scalar z=x_2-x_1 -> z'=y_2-y_1."""
    x = (Q(0), q(z))
    out = []
    for j in range(2):
        vals = [x[i] + S[i][j] for i in range(2) if S[i][j] is not None]
        if not vals:
            raise ValueError("empty output column")
        out.append(max(vals))
    return out[1] - out[0]


def active_sets(S, z):
    x = (Q(0), q(z))
    ans = []
    for j in range(2):
        vals = [(x[i] + S[i][j], i) for i in range(2) if S[i][j] is not None]
        m = max(v for v, _ in vals)
        ans.append(tuple(i for v, i in vals if v == m))
    return tuple(ans)


def hilbert_error(z, w):
    return abs(q(z) - q(w)) / 2


def enumerate_words(raw, pert, depth, seed=Q(0)):
    labels = tuple(raw)
    states = {"": (seed, seed)}
    rows = []
    for t in range(1, depth + 1):
        nxt = {}
        for word, (x, y) in states.items():
            for letter in labels:
                nxt[word + letter] = (
                    apply_matrix(raw[letter], x),
                    apply_matrix(pert[letter], y),
                )
        states = nxt
        winner = max(states.items(), key=lambda kv: hilbert_error(*kv[1]))
        rows.append((t, winner[0], winner[1], hilbert_error(*winner[1])))
    return rows


def tangent_at_tie(S, z, e):
    """Exact one-sided directional derivative for an unperturbed kernel."""
    active = active_sets(S, z)
    ev = (Q(0), q(e))
    out = [max(ev[i] for i in I) for I in active]
    return out[1] - out[0]


def fmt(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main():
    neg = None

    # Unique-selector translation: exact linear drift.
    ident = ((Q(0), neg), (neg, Q(0)))
    translate = ((Q(0), neg), (neg, Q(1)))
    z = w = Q(0)
    translation_rows = []
    for t in range(1, 9):
        z, w = apply_matrix(ident, z), apply_matrix(translate, w)
        translation_rows.append(fmt(hilbert_error(z, w)))
    print("unique-selector translation errors:", " ".join(translation_rows))

    # Two reflections.  Each letter alone has period two, but alternating
    # them creates a translation.
    anti0 = ((neg, Q(0)), (Q(0), neg))       # z -> -z
    anti1 = ((neg, Q(1)), (Q(0), neg))       # z -> 1-z
    raw = {"a": anti0, "b": anti0}
    pert = {"a": anti1, "b": anti0}
    print("reflection-switching maxima (depth, word, raw, pert, error):")
    for t, word, (x, y), err in enumerate_words(raw, pert, 10):
        print(t, word, fmt(x), fmt(y), fmt(err))
    for letter in "ab":
        x = y = Q(0)
        errs = []
        for _ in range(8):
            x = apply_matrix(raw[letter], x)
            y = apply_matrix(pert[letter], y)
            errs.append(hilbert_error(x, y))
        assert max(errs) <= Q(1, 2)
    x = y = Q(0)
    alt = []
    for letter in "ab" * 6:
        x = apply_matrix(raw[letter], x)
        y = apply_matrix(pert[letter], y)
        alt.append(hilbert_error(x, y))
    assert alt[-1] == 3
    print("alternating a,b errors:", " ".join(fmt(v) for v in alt))

    # Bounded clips with a symbolic middle-cell self-loop carrying drift.
    # S_c realizes z -> clip(z+c,0,1).
    def clip_matrix(c):
        c = q(c)
        return ((Q(0), Q(0)), (c - 1, c))

    S0, Sc = clip_matrix(0), clip_matrix(Q(1, 4))
    print("clip-pair errors under repetition:")
    for seed in (Q(-2), Q(0), Q(1, 5), Q(4, 5), Q(2)):
        x = y = seed
        errs = []
        for _ in range(12):
            x, y = apply_matrix(S0, x), apply_matrix(Sc, y)
            errs.append(hilbert_error(x, y))
        assert max(errs) <= Q(1, 2)
        print("seed", fmt(seed), ":", " ".join(fmt(v) for v in errs))

    # A tie face whose true directional map is a reset.  Arbitrary selector
    # branching would incorrectly permit identity or reflection derivatives.
    tie_reset = ((Q(0), Q(1)), (Q(-1), Q(0)))  # every z maps to 1
    assert apply_matrix(tie_reset, Q(1)) == 1
    assert active_sets(tie_reset, Q(1)) == ((0, 1), (0, 1))
    for e in range(-5, 6):
        assert tangent_at_tie(tie_reset, Q(1), Q(e)) == 0
    print("tie-reset active sets:", active_sets(tie_reset, Q(1)))
    print("tie-reset true directional outputs: all zero for e=-5,...,5")


if __name__ == "__main__":
    main()
