#!/usr/bin/env python3
"""Exact finite audit for Wave 45 scalar partition and block orbits."""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from anchored_conflict_r40_check import A5
from envelope_block_cover_r27 import A6, A8, A9


def spins(n):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qvalue(a):
    return max(abs(int(x @ a @ x)) for x in spins(len(a)))


def audit(a, m, name):
    n=len(a); q=qvalue(a)
    sels=list(itertools.combinations(range(n),m))
    xs=list(spins(n))
    states=[]; index={}
    edge_total=0
    for x in xs:
        raw=int(x @ a @ x)
        row=int((a @ x) @ (a @ x))
        for sigma in (-1,1):
            inc=[]
            for si,S in enumerate(sels):
                idx=np.asarray(S)
                c=sigma*int(x[idx] @ a[np.ix_(idx,idx)] @ x[idx])
                if 2*c-sigma*raw >= q:
                    inc.append(si)
            d={"x":x,"sigma":sigma,"row":row,"inc":tuple(inc),
               "alpha":Fraction(len(inc),len(sels))}
            index[(sigma,tuple(int(v) for v in x))]=len(states)
            states.append(d); edge_total += len(inc)

    D=len(states); N=len(sels)
    z0=sum((d["alpha"] for d in states),Fraction())/D
    assert z0 == Fraction(edge_total,D*N)
    # Direct pair averaging identity for Z_0.
    pair_count=0
    for d in states:
        for si in range(N):
            pair_count += int(si in d["inc"])
    assert Fraction(pair_count,D*N)==z0

    lam=n**-1.5
    active=[d for d in states if d["alpha"]]
    F=lambda d:-math.log(float(d["alpha"]))+lam*d["row"]
    best=min(active,key=F)
    zlam=sum(float(d["alpha"])*math.exp(-lam*d["row"])
             for d in states)/D
    assert zlam <= math.exp(-F(best))+1e-14
    assert zlam+1e-14 >= math.exp(-F(best))/D

    # Central fixed-cardinality block orbit around the scalar optimizer.
    k=n//2; alphas=[]; rows=[]
    for U in itertools.combinations(range(n),k):
        xp=best["x"].copy(); xp[list(U)]*=-1
        if xp[0] == -1: xp*=-1
        du=states[index[(best["sigma"],tuple(int(v) for v in xp))]]
        alphas.append(du["alpha"]); rows.append(du["row"])
    theta=Fraction(n*(n-1)-4*k*(n-k),n*(n-1))
    expected_row=Fraction(n*(n-1),1)+theta*(best["row"]-n*(n-1))
    assert Fraction(sum(rows),len(rows))==expected_row

    # Exact fixed-block energy mean and flipped-row average for every d,S.
    p2=Fraction(m*(m-1),n*(n-1))
    p3=Fraction(m*(m-1)*(m-2),n*(n-1)*(n-2))
    coeff=1-4*(p2-p3)
    for d in states:
        sigma=d["sigma"]; x=d["x"]
        rb=[]
        for S in sels:
            idx=np.asarray(S)
            b=-a.copy(); b[np.ix_(idx,idx)]=a[np.ix_(idx,idx)]
            rb.append(int((b @ x) @ (b @ x)))
            base=sigma*int(x @ b @ x)
            orb=[]
            for U in itertools.combinations(range(n),k):
                xp=x.copy(); xp[list(U)]*=-1
                orb.append(sigma*int(xp @ b @ xp))
            assert Fraction(sum(orb),len(orb))==theta*base
        expected_rb=coeff*d["row"]+4*(p2-p3)*n*(n-1)
        assert Fraction(sum(rb),N)==expected_rb

    mean_alpha=sum(alphas,Fraction())/len(alphas)
    empty=sum(not aa for aa in alphas)
    finite_h=[-math.log(float(aa)) for aa in alphas if aa]
    print({"name":name,"n":n,"m":m,"q":q,"Z0":str(z0),
           "minus_log_Z0":-math.log(float(z0)),"lambda":lam,
           "Zlambda":zlam,"scalar_min":F(best),
           "best":(best["row"],len(best["inc"]),best["sigma"],
                   tuple(int(v) for v in best["x"])),
           "central_k":k,"central_mean_alpha":str(mean_alpha),
           "central_empty":(empty,len(alphas)),
           "central_finite_mean_h":(sum(finite_h)/len(finite_h)
                                     if finite_h else math.inf)})


def main():
    audit(A5,4,"A5")
    audit(A6,5,"A6")
    audit(A8,6,"A8")
    audit(A9,7,"A9")
    print("PASS r45_scalar_partition_check")


if __name__=="__main__":
    main()
