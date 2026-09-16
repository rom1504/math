import itertools
import numpy as np


A9 = np.array([
    [0,1,1,-1,1,1,-1,1,1],
    [1,0,1,1,1,1,1,-1,1],
    [1,1,0,-1,1,-1,1,-1,-1],
    [-1,1,-1,0,1,-1,1,-1,1],
    [1,1,1,1,0,-1,1,1,1],
    [1,1,-1,-1,-1,0,1,1,1],
    [-1,1,1,1,1,1,0,1,-1],
    [1,-1,-1,-1,1,1,1,0,-1],
    [1,1,-1,1,1,1,-1,-1,0],
], dtype=int)


def spins(n):
    for bits in itertools.product((1, -1), repeat=n):
        yield np.array(bits, dtype=int)


def energy(A, x):
    return int(x @ A @ x // 2)


def exact_rank(rows):
    a = [list(map(int, row)) for row in rows]
    if not a:
        return 0
    m, n = len(a), len(a[0])
    rank = 0
    from fractions import Fraction
    a = [[Fraction(x) for x in row] for row in a]
    for col in range(n):
        pivot = next((i for i in range(rank, m) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        q = a[rank][col]
        a[rank] = [x/q for x in a[rank]]
        for i in range(m):
            if i != rank and a[i][col]:
                q = a[i][col]
                a[i] = [x-q*y for x,y in zip(a[i], a[rank])]
        rank += 1
    return rank


def endpoint_data(E, sign=1):
    xs = list(spins(len(E)))
    vals = [energy(E, x) for x in xs]
    endpoint = max(vals) if sign == 1 else min(vals)
    deficits = [(endpoint-v if sign == 1 else v-endpoint) for v in vals]
    levels = sorted(set(deficits))
    out = []
    for t in levels:
        cloud = [x for x,d in zip(xs, deficits) if d <= t]
        R = sum((np.outer(x,x) for x in cloud), np.zeros_like(E, dtype=float))/len(cloud)
        ev = np.linalg.eigvalsh(R)
        pos = ev[ev > 1e-8]
        out.append((t, len(cloud), exact_rank(cloud), float(pos.min()) if len(pos) else 0.0))
    return endpoint, out


z = np.ones(9, dtype=int)
for i in range(9):
    if (50 >> i) & 1:
        z[i] = -1
T = -np.diag(z) @ A9 @ np.diag(z)
U = [i for i in range(9) if (99 >> i) & 1]
V = [i for i in range(9) if i not in U]
B = T[np.ix_(U,U)]
D = T[np.ix_(V,V)]
C = T[np.ix_(U,V)]

print("U,V", U,V)
print("all-one traffic", int(np.ones(len(U), dtype=int) @ C @ np.ones(len(V), dtype=int)))
print("cross rank", exact_rank(C.tolist()), "op", np.linalg.svd(C, compute_uv=False)[0])
for name,E in (("B",B),("D",D)):
    p, layers = endpoint_data(E, 1)
    print(name, "p", p, "positive layers (threshold,count,rank,kappa)")
    print(layers)
    mn, nlayers = endpoint_data(E, -1)
    print(name, "min", mn, "negative layers", nlayers)

vals = [energy(T,x) for x in spins(9)]
print("full endpoints", min(vals), max(vals))
print("p excess", max(vals)-endpoint_data(B,1)[0]-endpoint_data(D,1)[0])
