#!/usr/bin/env python3
"""Independent exact audit of a four-reset compatible n=9 path."""

from itertools import product

A = [
 [0,-1, 1, 1, 1, 1,-1, 1, 1],
 [-1,0, 1, 1,-1, 1, 1, 1, 1],
 [1, 1, 0,-1,-1,-1,-1, 1,-1],
 [1, 1,-1, 0,-1,-1,-1,-1,-1],
 [1,-1,-1,-1, 0,-1, 1, 1, 1],
 [1, 1,-1,-1,-1, 0, 1, 1, 1],
 [-1,1,-1,-1, 1, 1, 0,-1,-1],
 [1, 1, 1,-1, 1, 1,-1, 0, 1],
 [1, 1,-1,-1, 1, 1,-1, 1, 0],
]
parts = (2, 2, 2, 1, 2)


def H(inds, x):
    return 2*sum(A[i][j]*x[a]*x[b]
                 for a,i in enumerate(inds)
                 for b,j in enumerate(inds) if a < b)


start = 0
tau = 1
total = 0
rows = []
for shore in parts[:-1]:
    inds = tuple(range(start, 9))
    vals = [(H(inds,x),x) for x in product((-1,1), repeat=len(inds))]
    P = max(v for v,x in vals)
    N = -min(v for v,x in vals)
    y = (1,)*len(inds)
    h = tau*H(inds,y)
    p = P if tau == 1 else N
    q = N if tau == 1 else P
    assert q > p
    x = (-1,)*shore + (1,)*(len(inds)-shore)
    tau_new = -tau
    assert tau_new*H(inds,x) == q
    g = q-h
    a = q+h
    assert a == 2*q-g
    # Ordinary cut identity in the fresh orientation.
    D = inds[:shore]; E = inds[shore:]
    cross = sum(A[i][j] for i in D for j in E)
    assert a == -4*tau_new*cross
    rows.append((inds, tau, p, q, h, g, a, cross))
    total += a
    tau = tau_new
    start += shore

P0 = max(H(tuple(range(9)),x) for x in product((-1,1),repeat=9))
N0 = -min(H(tuple(range(9)),x) for x in product((-1,1),repeat=9))
assert (P0,N0,total) == (24,32,104)
print("P0,N0,Q0,R0,sum_a", P0,N0,max(P0,N0),P0+N0,total)
for row in rows:
    print(row)
print("sum_a/Q0", f"{total}/{max(P0,N0)} = {total/max(P0,N0):.6f}")
