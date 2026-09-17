#!/usr/bin/env python3
"""Independent exhaustive order-six audit of the K<=4 construction."""

from fractions import Fraction
import sys

sys.path.insert(0, "/home/math/quadra/tmp")
import half_stopping_k4_r10 as hs  # noqa: E402
import minimax_allocation_r8 as mm  # noqa: E402


def main():
    n = 6
    count = 1 << ((n - 1) * (n - 2) // 2)
    worst = (Fraction(-1), None, None, None)
    minimum_slack = None
    for code in range(count):
        root = mm.make_tree(mm.matrix_code(n, code))
        mm.relabel(root)
        slacks = [
            hs.complement_slack(v)
            for v in mm.relabel(root)
            if v.children
        ]
        assert not slacks or min(slacks) >= 0
        if slacks:
            minimum_slack = min(slacks) if minimum_slack is None else min(
                minimum_slack, min(slacks)
            )
        t, theta, _, _, _ = hs.theorem_allocation(root)
        assert theta <= 4 * t <= 4
        if theta > worst[0]:
            worst = (theta, code, t, root.R)
    print("PASS matrices", count, "minimum complement slack", minimum_slack)
    print("maximum constructed theta/code/t/R", worst)


if __name__ == "__main__":
    main()
