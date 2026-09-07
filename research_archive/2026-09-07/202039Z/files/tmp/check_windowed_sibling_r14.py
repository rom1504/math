"""Exact checks for the Wave-14 disjoint-sibling certificate.

Scratch only.  All matrix calculations are integer calculations.  Radical
signs used for the temporal examples are certified by squaring positive
quantities, not by floating-point comparisons.
"""

from itertools import product


A5 = (
    (0, -1, 1, -1, 1),
    (-1, 0, -1, 1, 1),
    (1, -1, 0, 1, 1),
    (-1, 1, 1, 0, 1),
    (1, 1, 1, 1, 0),
)

A7 = (
    (0, 1, 1, 1, 1, 1, 1),
    (1, 0, 1, 1, -1, -1, 1),
    (1, 1, 0, 1, -1, 1, -1),
    (1, 1, 1, 0, 1, -1, -1),
    (1, -1, -1, 1, 0, -1, -1),
    (1, -1, 1, -1, -1, 0, -1),
    (1, 1, -1, -1, -1, -1, 0),
)

A8 = (
    (0, 1, 1, 1, 1, 1, 1, 1),
    (1, 0, 1, -1, 1, 1, -1, -1),
    (1, 1, 0, 1, -1, 1, -1, -1),
    (1, -1, 1, 0, -1, -1, -1, 1),
    (1, 1, -1, -1, 0, -1, 1, -1),
    (1, 1, 1, -1, -1, 0, 1, 1),
    (1, -1, -1, -1, 1, 1, 0, 1),
    (1, -1, -1, 1, -1, 1, 1, 0),
)


def energy(A, x, vertices=None):
    if vertices is None:
        vertices = tuple(range(len(A)))
    return sum(
        A[i][j] * x[ii] * x[jj]
        for ii, i in enumerate(vertices)
        for jj, j in enumerate(vertices)
    )


def stats(A, vertices=None):
    if vertices is None:
        vertices = tuple(range(len(A)))
    spins = [(1,) + tail for tail in product((-1, 1), repeat=len(vertices) - 1)]
    scored = [(energy(A, x, vertices), x) for x in spins]
    P = max(e for e, _ in scored)
    N = -min(e for e, _ in scored)
    return P, N, max(P, N), [x for e, x in scored if e == P], [x for e, x in scored if e == -N]


def bucket(A, U, X, p, sigma):
    """Return exact (c, partial, r, e, sibling energy, u) data."""
    U = tuple(U)
    X = tuple(X)
    D = tuple(i for i in U if i not in X)
    at = {v: k for k, v in enumerate(U)}
    y = tuple(p[at[i]] for i in X)
    Q_U = stats(A, U)[2]
    Q_X = stats(A, X)[2]
    h = energy(A, y, X)
    fields = tuple(sum(A[i][j] * y[k] for k, j in enumerate(X)) for i in D)
    L = sum(abs(v) for v in fields)
    c = max(0, 2 * L - (Q_X - sigma * h))
    partial = Q_U - Q_X
    r = max(0, c - partial)
    # Either sign is allowed at a zero field; choose +1.
    u = tuple(1 if sigma * v >= 0 else -1 for v in fields)
    z = tuple(y[X.index(i)] if i in X else u[D.index(i)] for i in U)
    completion_deficit = Q_U - sigma * energy(A, z, U)
    sibling_energy = energy(A, u, D)
    if r > 0:
        assert completion_deficit >= 0
        assert -sigma * sibling_energy == r + completion_deficit
    return c, partial, r, completion_deficit, sibling_energy, u, D


def root_residual_totals(A):
    U = tuple(range(len(A)))
    P, N, Q, ps, ns = stats(A, U)
    totals = []
    for p in ps:
        for n in ns:
            shores = [tuple(i for i in U if p[i] * n[i] == val) for val in (-1, 1)]
            if not all(shores):
                continue
            total = 0
            for X in shores:
                for sigma in (1, -1):
                    total += bucket(A, U, X, p, sigma)[2]
            totals.append(total)
    return (P, N, Q), totals


# All A5 and A8 root endpoint pairs reproduce the exact residual census.
assert root_residual_totals(A5)[0] == (8, 8, 8)
totals5 = root_residual_totals(A5)[1]
assert totals5.count(0) == 5 and totals5.count(4) == 20 and len(totals5) == 25
assert root_residual_totals(A8)[0] == (20, 20, 20)
assert set(root_residual_totals(A8)[1]) == {0}


# A5: a zero-completion-deficit certificate accompanies positive temporal
# demand.  It is also a positive-probability field-proportional outcome.
p5 = (1, -1, -1, -1, -1)
n5 = (1, -1, -1, 1, 1)
X5 = (0, 1, 2)
assert tuple(p5[i] * n5[i] for i in range(5)) == (1, 1, 1, -1, -1)
b5 = bucket(A5, range(5), X5, p5, -1)
assert b5[:5] == (4, 2, 2, 0, 2)
gauged_fields5 = tuple(p5[i] * sum(A5[i][j] * p5[j] for j in range(5)) for i in range(5))
assert gauged_fields5 == (0, 2, 0, 4, 2)
# With eligible set D={3,4}, deletion probabilities are (1,1/2), so the
# outcome D has probability 1/2.  Its demand is 6-24 sqrt(15)/25 > 0.
assert 150**2 > (24**2) * 15
# The other outcome {3} also has positive demand 8-64 sqrt(5)/25 > 0.
assert 200**2 > (64**2) * 5


# A concrete two-level A7 chain exercises the mixed-orientation aggregation.
# The first bucket has sigma=+ and the second sigma=-; sibling supports are
# disjoint.  Direct enumeration checks S_+<=N(W), S_-<=P(W), hence their sum
# is at most the range P(W)+N(W).
U0 = tuple(range(7))
p0 = (1, 1, 1, 1, -1, -1, 1)
b0 = bucket(A7, U0, (1, 2, 3, 6), p0, 1)
U1 = (1, 2, 3, 6)
p1 = (1, 1, 1, -1)
b1 = bucket(A7, U1, (1, 6), p1, -1)
assert b0[:5] == (12, 10, 2, 4, -6)
assert b1[:5] == (8, 6, 2, 0, 2)
assert set(b0[6]).isdisjoint(b1[6])
W = b0[6] + b1[6]
P_W, N_W, _, _, _ = stats(A7, W)
S_plus = b0[2] + b0[3]
S_minus = b1[2] + b1[3]
assert S_plus <= N_W and S_minus <= P_W
assert S_plus + S_minus <= P_W + N_W


# The A7 proper-suffix wall has a matching endpoint pair at its order-six
# child, but both residuals for retaining the named order-five child are zero.
B6 = tuple(tuple(A7[i + 1][j + 1] for j in range(6)) for i in range(6))
p6 = (1, -1, -1, -1, -1, 1)
n6 = (1, -1, -1, -1, -1, -1)
assert stats(B6)[2] == 10 and stats(B6, (0, 1, 2, 3, 4))[2] == 8
assert tuple(p6[i] * n6[i] for i in range(6)) == (1, 1, 1, 1, 1, -1)
assert bucket(B6, range(6), range(5), p6, 1)[2] == 0
assert bucket(B6, range(6), range(5), p6, -1)[2] == 0


print("PASS: exact completion identities, mixed-chain range bound, and A5/A7/A8 walls")
