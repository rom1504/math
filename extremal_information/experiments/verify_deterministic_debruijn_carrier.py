#!/usr/bin/env python3
"""Exact finite checks for deterministic de Bruijn support carriers."""

from __future__ import annotations

from itertools import product
import json


def words(q, length):
    return list(product(range(q), repeat=length))


def shift(state, letter):
    return state[1:] + (letter,)


def apply_word(state, word):
    for e in word:
        state = shift(state, e)
    return state


def suffix_support(states, suffix):
    if not suffix:
        return frozenset(states)
    return frozenset(s for s in states if s[-len(suffix):] == suffix)


def verify():
    fixed_point_checks = 0
    singleton_image_checks = 0
    carrier_edge_checks = 0
    rows = []
    for q in (2, 3):
        for m in range(1, 6):
            states = words(q, m)
            for length in range(1, m + 3):
                for word in words(q, length):
                    fixed = [s for s in states if apply_word(s, word) == s]
                    assert fixed
                    fixed_point_checks += 1
            for word in words(q, m):
                image = {apply_word(s, word) for s in states}
                assert image == {word}
                singleton_image_checks += 1
            for L in range(m):
                supports = {
                    u: suffix_support(states, u)
                    for depth in range(L + 1)
                    for u in words(q, depth)
                }
                for u, support in supports.items():
                    for e in range(q):
                        if len(u) < L:
                            target = supports[u + (e,)]
                            image = frozenset(shift(s, e) for s in support)
                            assert target.issubset(image)
                            shortfall = 0
                        else:
                            target = supports[()]
                            shortfall = 1
                        beta = 1 / (L + 1)
                        psi_source = -len(u) / (L + 1)
                        psi_target = -(0 if len(u) == L else len(u) + 1) / (L + 1)
                        assert abs(shortfall - (beta + psi_target - psi_source)) < 1e-12
                        carrier_edge_checks += 1
                rows.append({"q": q, "m": m, "L": L,
                             "carrier_states": len(supports)})
    return {
        "fixed_point_checks": fixed_point_checks,
        "length_m_singleton_images": singleton_image_checks,
        "carrier_edges_checked": carrier_edge_checks,
        "sample_constructions": rows,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
