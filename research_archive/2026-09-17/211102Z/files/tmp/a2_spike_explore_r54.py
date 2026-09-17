#!/usr/bin/env python3
"""Exploratory exact integer statistics for Wave 54 A^2 box cancellation."""

from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A as A10
from envelope_block_cover_r27 import A8, A9


def spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def records(a: np.ndarray, m: int):
    n = len(a)
    a2 = a @ a
    ans = []
    for ss in itertools.combinations(range(n), m):
        s = np.asarray(ss, dtype=int)
        t = np.asarray([i for i in range(n) if i not in ss], dtype=int)
        p = a[np.ix_(s, s)]
        vals = [(y, int(y @ p @ y)) for y in spins(m)]
        qs = max(abs(v) for _, v in vals)
        for y, raw in vals:
            if abs(raw) != qs:
                continue
            sig = 1 if raw > 0 else -1
            r = sig * y * (p @ y)
            u = a[np.ix_(t, s)] @ y
            I = int(r @ r)
            X = int(u @ u)
            d = a2[np.ix_(t, s)] @ y
            H = a2[np.ix_(t, t)] - (n - 1) * np.eye(n-m, dtype=np.int64)
            vals_w = []
            qvals = []
            lvals = []
            for ww in itertools.product((-1, 1), repeat=n-m):
                w = np.asarray(ww, dtype=np.int64)
                vals_w.append(I + X + (n-m)*(n-1) + int(2*d@w + w@H@w))
                qvals.append(int(w@H@w))
                lvals.append(int(d@w))
            lam = float(np.linalg.eigvalsh(H.astype(float))[0]) if len(t) else 0.0
            ans.append(dict(S=ss, y=tuple(map(int,y)), qs=qs, r=tuple(map(int,r)),
                            rmax=int(max(r)), I=I, X=X, d2=int(d@d),
                            d1=int(np.abs(d).sum()), Hf2=int((H*H).sum()),
                            lam=lam, qmin=min(qvals), lmax=max(map(abs,lvals)),
                            mu=I+X+(n-m)*(n-1), V=min(vals_w)))
    return ans


def report(name: str, a: np.ndarray, m: int):
    rr = records(a,m)
    print(f"\n{name} m={m} records={len(rr)}")
    for key, fun in [
        ("spike", lambda z:(z['rmax'],z['I'],-z['X'])),
        ("I", lambda z:(z['I'],-z['X'])),
        ("lowX",lambda z:(-z['X'],z['I'])),
        ("worst residual",lambda z:(z['V']-z['X']-(len(a)-m)*(len(a)-1), z['I'])),
        ("best cancel",lambda z:(z['mu']-z['V'],z['I'])),
    ]:
        z=max(rr,key=fun)
        print(key, {k:z[k] for k in ('S','y','qs','r','rmax','I','X','d2','d1','Hf2','lam','qmin','lmax','mu','V')})
    # Report extrema of proposed scale-free ratios (zero denominators omitted).
    print("max (V-X-k(n-1))/I", max((z['V']-z['X']-(len(a)-m)*(len(a)-1))/z['I'] for z in rr))
    print("min cancellation/I", min((z['mu']-z['V'])/z['I'] for z in rr))


def main():
    for name,a,ms in [('A8',A8,(4,5,6,7)),('A9',A9,(5,6,7,8)),('A10',A10,(5,6,7,8,9))]:
        for m in ms:
            report(name,a,m)


if __name__ == '__main__':
    main()
