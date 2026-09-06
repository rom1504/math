#!/usr/bin/env python3
"""MILP search for prescribed forced alternating singleton reset chains.

The inherited state is all-one on every suffix R_t={t,...,n-1}.
At t the endpoint flips only t and has fresh orientation theta_t, with
theta_t alternating.  Endpoint and strict-dominance constraints are imposed
against every Boolean state of the suffix.  Edge variables are exactly ±1.
"""

from itertools import product
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix


def solve_partition(parts, C=2.0, time_limit=300):
    """parts=(|D_0|,...,|D_{k-1}|,|terminal core|)."""
    n = sum(parts)
    k = len(parts) - 1
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ei = {e: z for z, e in enumerate(edges)}
    m = len(edges)
    # Last variable p is the subdominant root endpoint magnitude.
    nv = m + 1
    rows, cols, vals, lo, hi = [], [], [], [], []

    def ecoeff(inds, spin):
        # Affine doubled energy for a_ij=2*b_ij-1:
        # H=4 sum b_ij x_i x_j - 2 sum x_i x_j.
        pos = dict(zip(inds, spin))
        out = np.zeros(m)
        const = 0
        for ii, i in enumerate(inds):
            for j in inds[ii + 1:]:
                prod = pos[i] * pos[j]
                out[ei[(i, j)]] = 4 * prod
                const -= 2 * prod
        return out, const

    def add(co, const=0, lower=-np.inf, upper=np.inf):
        r = len(lo)
        nz = np.flatnonzero(co)
        rows.extend([r] * len(nz)); cols.extend(nz.tolist())
        vals.extend(co[nz].tolist()); lo.append(lower-const); hi.append(upper-const)

    score = np.zeros(nv)
    score_const = 0
    q0co = None
    q0const = None
    theta0 = -1
    starts = [0]
    for s in parts[:-1]:
        starts.append(starts[-1] + s)
    for t in range(k):
        start = starts[t]
        shore = parts[t]
        inds = tuple(range(start, n))
        theta = theta0 * ((-1) ** t)
        one = (1,) * len(inds)
        x = (-1,) * shore + (1,) * (len(inds) - shore)
        Hx, cx = ecoeff(inds, x)
        H1, c1 = ecoeff(inds, one)
        if t == 0:
            q0co = theta * Hx
            q0const = theta * cx
        # Ordinary reset cut layer a=theta(Hx-H1).
        score[:m] += theta * (Hx - H1)
        score_const += theta * (cx-c1)
        for z in product((-1, 1), repeat=len(inds)):
            Hz, cz = ecoeff(inds, z)
            # theta Hx >= theta Hz (fresh-channel maximum).
            co = np.zeros(nv); co[:m] = theta * (Hx - Hz)
            add(co, const=theta*(cx-cz), lower=0)
            # Strict dominance over old channel: theta Hx >= -theta Hz + 4.
            co = np.zeros(nv); co[:m] = theta * (Hx + Hz)
            add(co, const=theta*(cx+cz), lower=4)

    # p0 >= -theta0 H(z) for all root states, and p0 is minimized by objective.
    inds = tuple(range(n))
    for z in product((-1, 1), repeat=n):
        Hz, cz = ecoeff(inds, z)
        co = np.zeros(nv); co[:m] = theta0 * Hz; co[m] = 1
        # p + theta0 H(z) >= 0, i.e. p >= -theta0 H(z).
        add(co, const=theta0*cz, lower=0)

    # Maximize score - C(q0+p0).
    obj = score.copy()
    obj[:m] -= C * q0co
    obj[m] -= C
    A = coo_matrix((vals, (rows, cols)), shape=(len(lo), nv)).tocsr()
    c = -obj
    integrality = np.zeros(nv, dtype=int); integrality[:m] = 1
    lower = np.r_[np.zeros(m), 0.0]
    upper = np.r_[np.ones(m), np.inf]
    res = milp(c, integrality=integrality, bounds=Bounds(lower, upper),
               constraints=LinearConstraint(A, np.array(lo), np.array(hi)),
               options={"time_limit": time_limit, "mip_rel_gap": 0.0})
    if not res.success:
        return res, None
    bvec = np.rint(res.x[:m]).astype(int)
    avec = 2*bvec-1
    S = int(round(score[:m] @ bvec + score_const))
    q0 = int(round(q0co @ bvec + q0const))
    p0 = int(round(res.x[m]))
    mat = [[0] * n for _ in range(n)]
    for (i, j), a in zip(edges, avec):
        mat[i][j] = mat[j][i] = int(a)
    return res, (S, q0, p0, mat)


def solve(n, k, C=2.0, time_limit=300):
    return solve_partition((1,) * k + (n-k,), C=C, time_limit=time_limit)


if __name__ == "__main__":
    for n in range(4, 13):
        k = n - 2
        res, data = solve(n, k, C=2.0, time_limit=120)
        print("n,k,status", n, k, res.message)
        if data:
            S, q, p, A = data
            print("S,q,p,S/R,S/Q,obj", S, q, p, S/(q+p), S/q,
                  S-2*(q+p))
            if S > 2*(q+p):
                print("A=", A)
