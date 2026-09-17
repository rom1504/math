#!/usr/bin/env python3
"""Exact primal/dual audit for the corrected alternate t=1/2 macro tree."""

from fractions import Fraction as F


# node -> (imbalance cap, ((child, (orientation capacities)), ...))
TREE = {
    0: (F(25480), ((1, (0, 480)), (12, (480, 24000)))),
    1: (F(17640), ((2, (6860, 15680)), (7, (6860, 15680)))),
    2: (F(1960), ((3, (3430, 3430)), (4, (1470, 3430)))),
    3: (F(0), ()),
    4: (F(0), ((5, (980, 980)), (6, (980, 980)))),
    5: (F(0), ()),
    6: (F(0), ()),
    7: (F(1960), ((8, (3430, 3430)), (9, (1470, 3430)))),
    8: (F(0), ()),
    9: (F(0), ((10, (980, 980)), (11, (980, 980)))),
    10: (F(0), ()),
    11: (F(0), ()),
    12: (F(7840), ((13, (5880, 7840)), (16, (1960, 7840)))),
    13: (F(0), ((14, (980, 980)), (15, (980, 980)))),
    14: (F(0), ()),
    15: (F(0), ()),
    16: (F(0), ((17, (2940, 2940)), (18, (2940, 2940)))),
    17: (F(0), ()),
    18: (F(0), ()),
}
W = F(48000)
K = F(3575, 1959)


# Nonzero primal variables.  Missing variables are zero.
A = {
    (0, 0, 1): F(258560, 653),
    (0, 1, 1): F(24000),
    (1, 0, 1): F(8960),
    (1, 1, 1): F(7786240, 1959),
    (2, 0, 0): F(1960),
    (7, 0, 0): F(1703240, 1959),
    (12, 0, 1): F(7840),
}
WIN = {
    1: F(10293920, 653),
    2: F(1960),
    7: F(1703240, 1959),
    12: F(7840),
}
THETA = {
    0: K,
    1: F(1616, 1959),
    2: F(4, 7),
    7: F(3476, 13713),
    12: F(1),
}
PATH_Z = {
    (0, 0): F(1616, 1959),
    (0, 1): F(1),
    (1, 0): F(4, 7),
    (1, 1): F(3476, 13713),
    (2, 0): F(4, 7),
    (7, 0): F(3476, 13713),
    (12, 0): F(1),
}


# Nonzero dual variables.  Missing variables are zero.
Q = {
    (0, 0): F(16, 653),
    (0, 1): F(1175, 1959),
    (1, 0): F(1568, 1959),
    (1, 1): F(1568, 1959),
    (2, 0): F(343, 1959),
    (2, 1): F(343, 1959),
    (7, 0): F(343, 1959),
    (7, 1): F(343, 1959),
    (12, 0): F(784, 1959),
    (12, 1): F(784, 1959),
}
Y = {v: F(1, 19590) for v in (0, 1, 2, 3, 7, 8, 12)}
BUCKET_SLACK = {(0, 1, 1): F(49, 1880640)}
OBLIGATION_SLACK = {v: F(1, 19590) for v in (4, 9, 13, 16)}


def paths(v=0, prefix=()):
    children = TREE[v][1]
    if not children:
        return (prefix,)
    ans = ()
    for j, (child, _) in enumerate(children):
        ans += paths(child, prefix + ((v, j),))
    return ans


def audit_primal():
    for v, (d, children) in TREE.items():
        incoming = W if v == 0 else WIN.get(v, F(0))
        if v != 0:
            assert 0 <= incoming <= d
        outgoing = F(0)
        for j, (child, caps) in enumerate(children):
            outgoing += WIN.get(child, F(0))
            load = F(0)
            for s, cap in enumerate(caps):
                amount = A.get((v, j, s), F(0))
                assert 0 <= amount <= cap
                if cap:
                    load += amount / cap
                else:
                    assert amount == 0
                outgoing += amount
            zz = PATH_Z.get((v, j), F(0))
            assert zz >= load
            assert zz >= THETA.get(child, F(0))
        assert incoming == outgoing
        if children:
            assert THETA.get(v, F(0)) >= sum(
                (PATH_Z.get((v, j), F(0)) for j in range(len(children))), F(0)
            )
        else:
            assert THETA.get(v, F(0)) == 0
    assert THETA[0] == K


def audit_dual():
    for path in paths():
        assert sum((Q.get(e, F(0)) for e in path), F(0)) <= 1

    parent = {}
    for v, (_, children) in TREE.items():
        for j, (child, caps) in enumerate(children):
            parent[child] = v
            q = Q.get((v, j), F(0))
            for s, cap in enumerate(caps):
                if cap:
                    assert (
                        Y.get(v, F(0)) - BUCKET_SLACK.get((v, j, s), F(0))
                        <= q / cap
                    )

    for child, par in parent.items():
        assert (
            Y.get(par, F(0)) - Y.get(child, F(0))
            - OBLIGATION_SLACK.get(child, F(0)) <= 0
        )

    objective = W * Y[0]
    for (v, j, s), z in BUCKET_SLACK.items():
        objective -= F(TREE[v][1][j][1][s]) * z
    for v, s in OBLIGATION_SLACK.items():
        objective -= TREE[v][0] * s
    assert objective == K


def show_nonzero(label, values):
    print(label)
    for key in sorted(values, key=repr):
        print(" ", key, "=", values[key])


if __name__ == "__main__":
    audit_primal()
    audit_dual()
    show_nonzero("primal a", A)
    show_nonzero("primal w", WIN)
    show_nonzero("primal theta", THETA)
    show_nonzero("primal path z", PATH_Z)
    show_nonzero("dual q", Q)
    show_nonzero("dual y", Y)
    show_nonzero("dual bucket slack", BUCKET_SLACK)
    show_nonzero("dual obligation slack", OBLIGATION_SLACK)
    print("exact primal value = exact dual value =", K)
