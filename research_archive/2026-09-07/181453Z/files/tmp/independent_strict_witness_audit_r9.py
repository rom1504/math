#!/usr/bin/env python3
"""Independent exact audit of the strict 4+5 weighted witness."""

from fractions import Fraction as F
from itertools import product


S = tuple(range(4))
T = tuple(range(4, 9))
N = 9


def build_weights():
    w = {}
    for i in range(N):
        for j in range(i + 1, N):
            if i in S and j in S:
                z = 54
            elif i in S and j in (4, 5, 6):
                z = 28
            elif i in S and j == 7:
                z = 55
            elif i in S and j == 8:
                z = 29
            elif i in (4, 5, 6) and j in (4, 5, 6):
                z = 83
            elif i in (4, 5, 6) and j in (7, 8):
                z = -55
            elif (i, j) == (7, 8):
                z = 113
            else:
                raise AssertionError((i, j))
            w[i, j] = F(z, 672)
    return w


W = build_weights()


def cut(mask, vertices=tuple(range(N))):
    v = set(vertices)
    return sum(
        a
        for (i, j), a in W.items()
        if i in v and j in v and bool(mask >> i & 1) != bool(mask >> j & 1)
    )


def total(vertices):
    v = set(vertices)
    return sum(a for (i, j), a in W.items() if i in v and j in v)


def child_data(vertices):
    masks = []
    first = vertices[0]
    rest = vertices[1:]
    for bits in product((0, 1), repeat=len(rest)):
        mask = sum(b << i for b, i in zip(bits, rest))
        assert not (mask >> first & 1)
        masks.append(mask)
    vals = {mask: cut(mask, vertices) for mask in masks}
    a = total(vertices)
    m = min(vals.values())
    M = max(vals.values())
    h = 2 * a
    P = h - 4 * m
    NN = 4 * M - h
    return a, m, M, h, P, NN, [q for q, value in vals.items() if value == M]


def spin(mask, i):
    return -1 if mask >> i & 1 else 1


def cross_bilinear(ms, mt):
    return sum(
        a * spin(ms, i) * spin(mt, j)
        for (i, j), a in W.items()
        if i in S and j in T
    )


def main():
    states = [
        sum(bit << i for i, bit in enumerate(bits, start=1))
        for bits in product((0, 1), repeat=N - 1)
    ]
    values = {mask: cut(mask) for mask in states}
    endpoint = sum(1 << i for i in T)
    interior = [v for mask, v in values.items() if mask not in (0, endpoint)]
    assert values[0] == 0
    assert values[endpoint] == 1
    assert min(interior) == F(1, 112)
    assert max(interior) == F(111, 112)
    lower_pair_constant = min(
        value / (bin(mask).count("1") * (N - bin(mask).count("1")))
        for mask, value in values.items() if mask != 0
    )
    upper_pair_constant = min(
        (1 - value)
        / (bin(mask ^ endpoint).count("1")
           * (N - bin(mask ^ endpoint).count("1")))
        for mask, value in values.items() if mask != endpoint
    )

    sd = child_data(S)
    td = child_data(T)
    assert sd[:6] == (
        F(27, 56), F(0), F(9, 28), F(27, 28), F(27, 28), F(9, 28)
    )
    assert td[:6] == (
        F(1, 21), F(-55, 112), F(19, 112), F(2, 21), F(173, 84), F(49, 84)
    )

    a_s, m_s, M_s, h_s, p_s, n_s, max_s = sd
    a_t, m_t, M_t, h_t, p_t, n_t, max_t = td
    c = values[endpoint]
    q_s, q_t = max(p_s, n_s), max(p_t, n_t)
    b_s = 2 * c - q_s + abs(h_s)
    b_t = 2 * c - q_t + abs(h_t)
    imbalance = 2 * abs(h_s + h_t)
    stopping_gap = imbalance - b_s - b_t
    strong_gap = q_s + q_t - (4 * c - abs(h_s + h_t))
    cut_gap = a_s + a_t - m_s - m_t - c
    assert (b_s, b_t, imbalance, stopping_gap) == (
        F(2), F(1, 28), F(89, 42), F(1, 12)
    )
    assert strong_gap == F(1, 12)
    assert cut_gap == F(1, 48)

    zvals = [abs(cross_bilinear(ms, mt)) for ms in max_s for mt in max_t]
    max_z = max(zvals)
    averaging_lhs = 2 * (M_s + M_t) + max_z
    averaging_rhs = a_s + a_t - m_s - m_t
    assert averaging_lhs < averaging_rhs

    print("parent nonendpoint interval", min(interior), max(interior))
    print("lower/upper gap per macro disagreement", lower_pair_constant,
          upper_pair_constant)
    print("S a,m,M,h,P,N,Q,b", *sd[:6], q_s, b_s)
    print("T a,m,M,h,P,N,Q,b", *td[:6], q_t, b_t)
    print("I, I-B, strong gap, cut gap", imbalance, stopping_gap, strong_gap, cut_gap)
    print("maximum-cut mask counts", len(max_s), len(max_t))
    print("max |z|, 10.442 lhs/rhs/gap", max_z, averaging_lhs, averaging_rhs,
          averaging_rhs - averaging_lhs)


if __name__ == "__main__":
    main()
