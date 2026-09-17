"""Exact verifier for an order-12 tie-preserving reset-chain counterexample.

It disproves the candidate bound sum(reset ordinary deficits) <= 2(P+N)
even when every orientation change is forced by strict one-sided dominance.
"""

import sys

sys.path.insert(0, "/home/math/quadra/tmp")
from reset_bellman_r6 import Instance, edge_list  # noqa: E402


A = [
    [0, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1],
    [-1, 0, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1],
    [-1, -1, 0, 1, -1, 1, 1, -1, 1, 1, 1, -1],
    [-1, -1, 1, 0, -1, 1, 1, 1, 1, 1, 1, 1],
    [-1, -1, -1, -1, 0, 1, 1, 1, 1, -1, 1, 1],
    [-1, -1, 1, 1, 1, 0, 1, -1, -1, -1, -1, -1],
    [-1, 1, 1, 1, 1, 1, 0, -1, 1, -1, -1, -1],
    [1, -1, -1, 1, 1, -1, -1, 0, -1, 1, 1, 1],
    [-1, 1, 1, 1, 1, -1, 1, -1, 0, 1, 1, 1],
    [1, -1, 1, 1, -1, -1, -1, 1, 1, 0, 1, -1],
    [-1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 0, -1],
    [-1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 0],
]
BLOCKS = (2, 3, 2, 2, 2, 1)


def main():
    n = len(A)
    assert all(A[i][j] == A[j][i] for i in range(n) for j in range(n))
    signs = tuple(A[i][j] for i, j in edge_list(n))
    inst = Instance(n, signs)
    starts = [0]
    for d in BLOCKS:
        starts.append(starts[-1] + d)

    rows = []
    rho = -1
    total = 0
    for t in range(len(BLOCKS) - 1):
        mask = sum(1 << i for i in range(starts[t], n))
        child = sum(1 << i for i in range(starts[t + 1], n))
        tau = 1 if t % 2 == 0 else -1
        metric = inst.metric(mask)
        old = metric.p if rho == 1 else metric.n
        new = metric.p if tau == 1 else metric.n
        assert tau == -rho
        assert new == metric.q and new > old
        assert (tau, child) in metric.transitions
        a = metric.q + rho * metric.h
        # Direct cut formula in the inherited-one gauge.
        cross = sum(
            A[i][j]
            for i in range(starts[t], starts[t + 1])
            for j in range(starts[t + 1], n)
        )
        assert a == -4 * tau * cross
        rows.append((t, metric.p, metric.n, metric.h, rho, tau, a))
        total += a
        rho = tau

    root = inst.metric((1 << n) - 1)
    assert (root.p, root.n, root.r) == (44, 40, 84)
    assert [r[-1] for r in rows] == [48, 60, 32, 24, 8]
    assert total == 172 > 2 * root.r == 168
    print("rows (t,P,N,H,rho,tau,a):", rows)
    print("root range", root.r, "sum reset a", total, "ratio", total / root.r)


if __name__ == "__main__":
    main()
