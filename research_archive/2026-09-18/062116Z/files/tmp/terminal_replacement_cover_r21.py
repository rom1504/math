#!/usr/bin/env python3
"""Exact finite audit for the Wave 21 terminal translate-cover memo."""

from __future__ import annotations

import itertools

import numpy as np

from check_response_dual_r16 import A9
from puncture_cycles_r20 import A6


def main() -> None:
    A = np.asarray(A9, dtype=int)
    C = np.asarray(A6, dtype=int)
    terminal = list(range(3, 9))
    outside = list(range(3))
    m = len(terminal)
    B = A[np.ix_(terminal, terminal)]

    states = [
        (sigma, np.asarray((1,) + tail, dtype=int))
        for sigma in (-1, 1)
        for tail in itertools.product((-1, 1), repeat=m - 1)
    ]
    assert len(states) == 2**m
    index = {(sigma, tuple(map(int, x))): k
             for k, (sigma, x) in enumerate(states)}

    external = []
    energy_b = []
    energy_c = []
    for sigma, x_t in states:
        extensions = []
        for x_u_tuple in itertools.product((-1, 1), repeat=len(outside)):
            x = np.empty(len(A), dtype=int)
            x[terminal] = x_t
            x[outside] = x_u_tuple
            extensions.append(
                sigma * int(x @ A @ x) - sigma * int(x_t @ B @ x_t)
            )
        external.append(max(extensions))
        energy_b.append(sigma * int(x_t @ B @ x_t))
        energy_c.append(sigma * int(x_t @ C @ x_t))

    q, q_b, q_c = 24, max(energy_b), max(energy_c)
    assert (q_b, q_c) == (14, 10)
    epsilon = q_b - q_c
    surplus = np.asarray(external) - (q - q_b)
    deficit_b = q_b - np.asarray(energy_b)
    deficit_c = q_c - np.asarray(energy_c)
    assert np.all(surplus <= deficit_b)
    assert int(np.max(surplus - deficit_b)) == 0
    assert int(np.sum(surplus == deficit_b)) == 17
    parent_deficits = sorted(map(int, deficit_b - surplus))
    assert {value: parent_deficits.count(value) for value in set(parent_deficits)} == {
        0: 17, 8: 20, 16: 18, 24: 8, 32: 1
    }
    assert (int(np.min(surplus)), int(np.max(surplus))) == (-4, 20)

    translate_margins = []
    for h_sigma, h_x in states:
        values = []
        for d_index, (sigma, x) in enumerate(states):
            u_index = index[(h_sigma * sigma, tuple(map(int, h_x * x)))]
            values.append(int(surplus[d_index] - deficit_c[u_index]))
        translate_margins.append(max(values))
    assert min(translate_margins) == 8

    layer_data = []
    cover_sum = 0
    for level in sorted(set(map(int, deficit_c))):
        near = int(np.sum(deficit_c <= level))
        resource = int(np.sum(surplus >= epsilon + level))
        layer_data.append((level, near, resource))
        cover_sum += near * resource
    assert layer_data == [(0, 12, 35), (4, 32, 17),
                          (16, 52, 1), (20, 64, 0)]
    assert cover_sum >= len(states)

    print({"q": q, "q_B": q_b, "q_C": q_c, "epsilon": epsilon,
           "min_replacement_margin": min(translate_margins),
           "layers": layer_data})
    print("PASS: exact terminal translate-cover and A9 compensation wall")


if __name__ == "__main__":
    main()
