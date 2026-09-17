#!/usr/bin/env python3
"""Exact finite diagnostics and Euler identities for Wave 44."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from anchored_conflict_r40_check import A5
from envelope_block_cover_r27 import A6, A8, A9


def spins(n):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qvalue(a):
    return max(abs(int(x @ a @ x)) for x in spins(len(a)))


def analyze(a, m, name):
    n = len(a)
    q = qvalue(a)
    sels = list(itertools.combinations(range(n), m))
    xs = list(spins(n))
    cols = []
    lookup = {}
    for xi, x in enumerate(xs):
        row = int((a @ x) @ (a @ x))
        raw = int(x @ a @ x)
        for sigma in (-1, 1):
            e = sigma * raw
            inc = []
            for si, ss in enumerate(sels):
                idx = np.asarray(ss)
                c = sigma * int(x[idx] @ a[np.ix_(idx, idx)] @ x[idx])
                if 2 * c - e >= q:
                    inc.append(si)
            alpha = len(inc) / len(sels)
            sid = (sigma, tuple(int(v) for v in x))
            lookup[sid] = len(cols)
            cols.append(dict(sigma=sigma, x=x, E=e, row=row, inc=tuple(inc),
                             alpha=alpha,
                             h=(-math.log(alpha) if alpha else math.inf)))

    # Verify exact single-vertex row and incidence-score change identities.
    for d in cols:
        sigma, x = d["sigma"], d["x"]
        sgn = sigma * a * np.outer(x, x)
        r = sgn.sum(axis=1)
        assert int(r.sum()) == d["E"]
        assert int(r @ r) == d["row"]
        total_delta = 0
        for i in range(n):
            xp = x.copy(); xp[i] *= -1
            if xp[0] == -1: xp *= -1
            dp = cols[lookup[(sigma, tuple(int(v) for v in xp))]]
            t_i = int(sgn[i] @ r)
            delta = 4 * ((n - 1) - t_i)
            assert dp["row"] - d["row"] == delta
            total_delta += delta
            for si, ss in enumerate(sels):
                idx = np.asarray(ss)
                c = sigma * int(x[idx] @ a[np.ix_(idx, idx)] @ x[idx])
                cp = sigma * int(xp[idx] @ a[np.ix_(idx, idx)] @ xp[idx])
                L = 2*c-d["E"]
                Lp = 2*cp-dp["E"]
                a_i = int(sgn[i, idx].sum()) if i in ss else 0
                expected = L + 4*int(r[i]) - 8*a_i
                assert Lp == expected
        assert total_delta == 4 * (n*(n-1)-d["row"])

    active = [d for d in cols if d["alpha"]]
    pareto = [d for d in active if not any(
        (e["h"] <= d["h"] + 1e-12 and e["row"] <= d["row"] and
         (e["h"] < d["h"] - 1e-12 or e["row"] < d["row"]))
        for e in active)]
    pairs = sorted(set((d["row"], len(d["inc"]), d["E"]) for d in pareto))

    # Representatives of each Pareto type, including the exact one-flip
    # incidence boundary.  This diagnoses whether scalar Euler descent has
    # any automatic incidence retention.
    boundary = []
    for typ in pairs:
        d = next(z for z in pareto
                 if (z["row"], len(z["inc"]), z["E"]) == typ)
        sigma, x = d["sigma"], d["x"]
        desc = []
        for i in range(n):
            xp=x.copy(); xp[i]*=-1
            if xp[0] == -1: xp*=-1
            dp=cols[lookup[(sigma,tuple(int(v) for v in xp))]]
            if dp["row"] < d["row"]:
                desc.append((i,dp["row"]-d["row"],len(dp["inc"])))
        incidences=[]
        for si in d["inc"]:
            ss=sels[si]; idx=np.asarray(ss)
            b=-a.copy()
            b[np.ix_(idx,idx)]=a[np.ix_(idx,idx)]
            vals=[]
            for y in xs:
                rawb=int(y @ b @ y)
                vals.extend((rawb,-rawb))
            M=max(vals)
            L=sigma*int(x @ b @ x)
            retained=0
            for i in range(n):
                xp=x.copy(); xp[i]*=-1
                if sigma*int(xp @ b @ xp) >= q:
                    retained+=1
            fiber=sum(v>=q for v in vals)
            incidences.append((si,L-q,M-q,retained,fiber))
        boundary.append((typ,desc,incidences))

    # Exact hypergeometric variance formula for c_S.
    p2 = m*(m-1)/(n*(n-1))
    p3 = m*(m-1)*(m-2)/(n*(n-1)*(n-2)) if n >= 3 else 0
    p4 = (m*(m-1)*(m-2)*(m-3)/(n*(n-1)*(n-2)*(n-3))
          if n >= 4 else 0)
    for d in cols:
        vals=[]
        sigma,x=d["sigma"],d["x"]
        for ss in sels:
            idx=np.asarray(ss)
            vals.append(sigma*int(x[idx] @ a[np.ix_(idx,idx)] @ x[idx]))
        mean=sum(vals)/len(vals)
        var=sum((v-mean)**2 for v in vals)/len(vals)
        E,R=d["E"],d["row"]
        # c=2 sum selected unordered signed edges.
        # Explanation: 4[p2*N + 2p3*A + 2p4*D]-(p2E)^2,
        # with N=n(n-1)/2, A=(R-n(n-1))/2, and
        # D=(E^2/4-N)/2-A.
        formula2=(4*(p2*n*(n-1)/2
                     +2*p3*(R-n*(n-1))/2
                     +2*p4*((E*E/4-n*(n-1)/2)/2
                            -(R-n*(n-1))/2))
                  -(p2*E)**2)
        assert abs(var-formula2) < 1e-8, (name,var,formula2)

    row_counts=Counter(d["row"] for d in cols)
    best=max(active,key=lambda d:(d["alpha"],-d["row"]))
    low=min(active,key=lambda d:(d["row"],d["h"]))
    hard_cliffs=[]
    for d in active:
        descending=[]
        for i in range(n):
            xp=d["x"].copy(); xp[i]*=-1
            if xp[0] == -1: xp*=-1
            dp=cols[lookup[(d["sigma"],tuple(int(v) for v in xp))]]
            if dp["row"] < d["row"]:
                descending.append(dp)
        if d["row"] > n*(n-1) and descending and all(not z["inc"] for z in descending):
            hard_cliffs.append(d)
    cliff_state=(max(hard_cliffs,key=lambda d:(len(d["inc"]),d["row"]),
                     default=None))
    cliff=None
    if cliff_state is not None:
        descent=[]
        for i in range(n):
            xp=cliff_state["x"].copy(); xp[i]*=-1
            if xp[0] == -1: xp*=-1
            dp=cols[lookup[(cliff_state["sigma"],tuple(int(v) for v in xp))]]
            if dp["row"] < cliff_state["row"]:
                descent.append((i,dp["row"]-cliff_state["row"],len(dp["inc"])))
        cliff={"row":cliff_state["row"],"degree":len(cliff_state["inc"]),
               "E":cliff_state["E"],"sigma":cliff_state["sigma"],
               "x":tuple(int(v) for v in cliff_state["x"]),
               "selectors":tuple(sels[i] for i in cliff_state["inc"]),
               "descent":tuple(descent)}
    print({"name":name,"n":n,"m":m,"q":q,"active":len(active),
           "max_alpha":best["alpha"],"max_alpha_row":best["row"],
           "min_active_row":low["row"],"min_active_alpha":low["alpha"],
           "pareto":pairs,"boundary":boundary,
           "high_row_all_descent_cliff":cliff,
           "row_hist":sorted(row_counts.items())})
    return cols,sels


def main():
    analyze(A5,4,"A5")
    analyze(A6,5,"A6")
    analyze(A8,6,"A8")
    analyze(A9,7,"A9")
    print("PASS r44_direct_incidence_check")


if __name__ == "__main__":
    main()
