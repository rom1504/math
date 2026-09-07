#!/usr/bin/env python3
"""Exact audit of the recursive complement/half-stopping K<=4 construction.

Scratch only.  Uses a fixed endpoint tree from minimax_allocation_r8.
All allocation arithmetic and all checks are rational/integer exact.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, "/home/math/quadra/tmp")
import minimax_allocation_r8 as m  # noqa: E402
from verify_kmin_witnesses_r8 import A7, A8  # noqa: E402


def complement_slack(v):
    """B + (dS+dT)/2 - dU, for the best one-orientation buckets."""
    if not v.children:
        return Fraction(0)
    B = sum(max(caps) for caps in v.caps)
    D = sum(ch.d for ch in v.children)
    return Fraction(B) + Fraction(D, 2) - v.d


def theorem_allocation(root):
    """Return allocations realizing the theorem's root K<=4 certificate."""
    alloc = defaultdict(lambda: Fraction(0))
    obligation = {root.label: Fraction(root.R)}

    if root.R == 0:
        return Fraction(0), Fraction(0), alloc, obligation, m.relabel(root)

    C = sum(sum(caps) for caps in root.caps)
    D = sum(ch.d for ch in root.children)
    assert C + D >= root.R
    t0 = Fraction(root.R, C + D)
    for j, caps in enumerate(root.caps):
        for s, cap in enumerate(caps):
            alloc[root.label, j, s] = t0 * cap
        obligation[root.children[j].label] = t0 * root.children[j].d

    def service(v, t):
        w = obligation.get(v.label, Fraction(0))
        assert 0 <= t <= 1 and 0 <= w <= t * v.d
        if not v.children:
            assert w == 0
            return

        # Available slots: one best parent bucket on each outgoing edge at
        # load t, followed by child obligations at normalized level t/2.
        slots = []
        for j, caps in enumerate(v.caps):
            s = max(range(2), key=lambda k: caps[k])
            slots.append(("a", j, s, t * caps[s]))
        for j, ch in enumerate(v.children):
            slots.append(("w", j, None, t * ch.d / 2))

        assert sum(slot[3] for slot in slots) >= w, (
            v.ids, w, t, complement_slack(v), slots
        )
        rem = w
        for kind, j, s, cap in slots:
            take = min(rem, cap)
            if kind == "a":
                alloc[v.label, j, s] += take
            else:
                obligation[v.children[j].label] = take
            rem -= take
        assert rem == 0
        for ch in v.children:
            service(ch, t / 2)

    for ch in root.children:
        service(ch, t0)

    # Exact conservation and capacity checks.
    nodes = m.relabel(root)
    for v in nodes:
        incoming = obligation.get(v.label, Fraction(0))
        if v is root:
            incoming = Fraction(root.R)
        outgoing = Fraction(0)
        for j, ch in enumerate(v.children):
            outgoing += obligation.get(ch.label, Fraction(0))
            for s, cap in enumerate(v.caps[j]):
                a = alloc[v.label, j, s]
                assert 0 <= a <= cap
                outgoing += a
        assert incoming == outgoing, (v.ids, incoming, outgoing)

    def theta(v):
        if not v.children:
            return Fraction(0)
        ans = Fraction(0)
        for j, ch in enumerate(v.children):
            load = sum(
                alloc[v.label, j, s] / cap
                for s, cap in enumerate(v.caps[j]) if cap
            )
            ans += max(load, theta(ch))
        return ans

    return t0, theta(root), alloc, obligation, nodes


def audit(name, A):
    root = m.make_tree(A)
    m.relabel(root)
    slacks = [complement_slack(v) for v in m.relabel(root) if v.children]
    assert min(slacks) >= 0
    t, th, _, _, nodes = theorem_allocation(root)
    assert th <= 4 * t <= 4
    print(name, "R", root.R, "nodes", len(nodes), "t", t,
          "theta", th, "4t", 4*t, "min complement slack", min(slacks))


def main():
    audit("triangle", m.all_sign(3))
    audit("A37/(10.420)", m.A37)
    audit("A7", A7)
    audit("A8", A8)
    audit("reset12", m.RESET12)
    audit("rank2 m=4", m.family_A(4))


if __name__ == "__main__":
    main()
