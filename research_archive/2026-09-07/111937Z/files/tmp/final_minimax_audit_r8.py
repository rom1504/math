#!/usr/bin/env python3
"""Final exact rational certificate audit for Wave 8 allocation minimax."""

from fractions import Fraction
from itertools import product
import sys

sys.path.insert(0, "/home/math/quadra/tmp")
import minimax_allocation_r8 as m  # noqa: E402


def clique_sig(n):
    if n == 1:
        return (0, ())
    child = clique_sig(n // 2)
    return (n * (n - 2), (((n, n * n // 2), child), ((n, n * n // 2), child)))


def am_tree(n):
    C = n * n + n
    sig = (0, (((C, C), clique_sig(n)), ((C, C), clique_sig(n))))
    root = m.node_from_signature(sig)
    root.P = root.N = 2 * n * n
    return root


def root_signatures_with_fixed_positive(a, p):
    P, N, ps, ns = m.extrema(a)
    assert p in ps
    out = set()
    for n in ns:
        ratio = tuple(x * y for x, y in zip(p, n))
        sides = [tuple(i for i, z in enumerate(ratio) if z == s) for s in (-1, 1)]
        branches = []
        for side in sides:
            other = tuple(i for i in range(len(a)) if i not in side)
            child = m.induced(a, side)
            pc, nc, _, _ = m.extrema(child)
            q = max(pc, nc)
            h = m.energy(child, tuple(p[i] for i in side))
            L = sum(abs(sum(a[j][i] * p[i] for i in side)) for j in other)
            caps = tuple(sorted((max(0, 2 * L - (q - h)), max(0, 2 * L - (q + h)))))
            branches.append(tuple((caps, sig) for sig in m.tree_signatures(child)))
        for pair in product(*branches):
            out.add((abs(P - N), tuple(sorted(pair, key=repr))))
    return out


def exact_value(root, W):
    p = m.exact_primal_certificate(root, W)[0]
    d = m.exact_dual_certificate(root, W)[0]
    assert p == d
    return p


def main():
    examples = (
        ("triangle", m.make_tree(m.all_sign(3)), 8, Fraction(2)),
        ("A37", m.make_tree(m.A37), 32, Fraction(5, 2)),
        ("A_m, m=4", am_tree(4), 64, Fraction(12, 5)),
        ("A_m, m=16", am_tree(16), 1024, Fraction(128, 55)),
        ("reset lex, W=84", m.make_tree(m.RESET12), 84, Fraction(22, 7)),
        ("reset lex, W=40", m.make_tree(m.RESET12), 40, Fraction(58, 51)),
    )
    for name, root, W, expected in examples:
        value = exact_value(root, W)
        assert value == expected
        print(name, value)

    P, N, _, _ = m.extrema(m.RESET12)
    all_duals = []
    all_gap_duals = []
    best_sig = None
    best_gap_sig = None
    for sig in m.tree_signatures(m.RESET12):
        root = m.node_from_signature(sig)
        root.P, root.N = P, N
        value = m.exact_dual_certificate(root, 84)[0]
        gap_value = m.exact_dual_certificate(root, 40)[0]
        all_duals.append(value)
        all_gap_duals.append(gap_value)
        if best_sig is None or value < best_sig[0]:
            best_sig = (value, sig)
        if best_gap_sig is None or gap_value < best_gap_sig[0]:
            best_gap_sig = (gap_value, sig)
    root = m.node_from_signature(best_sig[1])
    root.P, root.N = P, N
    assert min(all_duals) == m.exact_primal_certificate(root, 84)[0] == Fraction(28, 13)
    print("reset all ties, W=84", len(all_duals), "signatures, min", min(all_duals))
    root = m.node_from_signature(best_gap_sig[1])
    root.P, root.N = P, N
    assert min(all_gap_duals) == m.exact_primal_certificate(root, 40)[0] == Fraction(19, 23)
    print("reset all ties, W=40", len(all_gap_duals), "signatures, min", min(all_gap_duals))

    fresh = (-1, -1) + (1,) * 10
    fixed_sigs = root_signatures_with_fixed_positive(m.RESET12, fresh)
    fixed_duals = []
    fixed_range_duals = []
    fixed_best = None
    fixed_range_best = None
    for sig in fixed_sigs:
        root = m.node_from_signature(sig)
        root.P, root.N = P, N
        value = m.exact_dual_certificate(root, 40)[0]
        range_value = m.exact_dual_certificate(root, 84)[0]
        fixed_duals.append(value)
        fixed_range_duals.append(range_value)
        if fixed_best is None or value < fixed_best[0]:
            fixed_best = (value, sig)
        if fixed_range_best is None or range_value < fixed_range_best[0]:
            fixed_range_best = (range_value, sig)
    root = m.node_from_signature(fixed_best[1])
    root.P, root.N = P, N
    assert min(fixed_duals) == m.exact_primal_certificate(root, 40)[0] == Fraction(104, 119)
    print("reset prescribed fresh p, W=40", len(fixed_duals), "signatures, min", min(fixed_duals))
    root = m.node_from_signature(fixed_range_best[1])
    root.P, root.N = P, N
    assert min(fixed_range_duals) == m.exact_primal_certificate(root, 84)[0] == Fraction(22, 9)
    print("reset prescribed fresh p, W=84", len(fixed_range_duals), "signatures, min", min(fixed_range_duals))


if __name__ == "__main__":
    main()
