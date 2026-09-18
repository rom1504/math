#!/usr/bin/env python3
"""Independent exact-data audit of the n=7 and n=8 K_min witnesses."""

from fractions import Fraction

from kmin_carleson_r8 import EndpointSystem, Model


A7 = (
    (0,1,1,1,1,1,1),
    (1,0,1,-1,1,1,-1),
    (1,1,0,1,1,1,-1),
    (1,-1,1,0,1,1,1),
    (1,1,1,1,0,-1,1),
    (1,1,1,1,-1,0,-1),
    (1,-1,-1,1,1,-1,0),
)

A8 = (
    (0,1,1,1,1,1,1,1),
    (1,0,1,-1,1,1,-1,-1),
    (1,1,0,1,-1,1,-1,-1),
    (1,-1,1,0,-1,-1,-1,1),
    (1,1,-1,-1,0,-1,1,-1),
    (1,1,1,-1,-1,0,1,1),
    (1,-1,-1,-1,1,1,0,1),
    (1,-1,-1,1,-1,1,1,0),
)


def root_data(A):
    S = EndpointSystem(A)
    U = (1 << len(A)) - 1
    P, N, ps, ns = S.endpoints(U)
    acts = S.actions(U)
    return S, U, P, N, ps, ns, acts


def audit7():
    S, U, P, N, ps, ns, acts = root_data(A7)
    assert (P, N, len(ps), len(ns), len(acts)) == (22, 18, 2, 2, 4)
    for act in acts:
        rows = sorted(act.shores, key=lambda z: z.imbalance, reverse=True)
        a, b = rows
        assert (a.caps, a.imbalance, a.P, a.N, a.h, a.L, bin(a.child).count("1")) == ((16,12),4,2,6,2,10,3)
        assert (b.caps, b.imbalance, b.P, b.N, b.h, b.L, bin(b.child).count("1")) == ((12,12),0,8,8,0,10,4)
        # The passed obligation four is served by one child capacity four.
        child_actions = S.actions(a.child)
        assert all(any(4 in shore.caps for shore in ca.shores) for ca in child_actions)

    # Exact lower certificate.  If x<=4 is passed, local mass is >=40-x>=36.
    # For edge maxima q1,q2, local mass <=(12q1+4)+12q2.
    lower = Fraction(36 - 4, 12)
    assert lower == Fraction(8, 3)
    # Exact primal: pass 4; local allocations 16 and 20.  Child theta=1.
    upper = max(Fraction(16,16), Fraction(1)) + Fraction(20,12)
    assert upper == lower

    M = Model(S)
    res = M.solve(time_limit=120)
    assert res.status == 0 and abs(res.fun - float(lower)) < 1e-8
    print("n=7: P,N,R", P, N, P+N, "endpoint actions", len(acts), "K_min", lower)


def audit8():
    S, U, P, N, ps, ns, acts = root_data(A8)
    assert (P, N, len(ps), len(ns), len(acts)) == (20, 20, 4, 4, 16)
    for act in acts:
        for shore in act.shores:
            assert (shore.caps, shore.imbalance, shore.P, shore.N, shore.h, shore.L, bin(shore.child).count("1")) == ((12,12),0,8,8,0,10,4)
    # No mass can pass.  Every allocation coordinate has denominator 12,
    # hence theta_root=sum a/12=40/12 for every feasible allocation.
    exact = Fraction(40,12)
    assert exact == Fraction(10,3)
    M = Model(S)
    res = M.solve(time_limit=120)
    assert res.status == 0 and abs(res.fun - float(exact)) < 1e-8
    print("n=8: P,N,R", P, N, P+N, "endpoint actions", len(acts), "K_min", exact)


if __name__ == "__main__":
    audit7()
    audit8()
