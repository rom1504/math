#!/usr/bin/env python3
"""Checks for the Wave 30 matched parent common-mode memo."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, oriented_cuts, restriction_key
from parent_cut_entropy_r26 import Endpoint


def harmonic_endpoint(A: np.ndarray, m: int, beta: float):
    """Directly enumerate the all-k harmonic likelihood identity."""
    n=len(A); k=n-m; cuts=oriented_cuts(n)
    energy=np.einsum("bij,ij->b",cuts,A,optimize=True).astype(float)
    logZA=float(logsumexp(beta*energy))
    nu=np.exp(beta*energy-logZA)
    selectors=list(itertools.combinations(range(n),m))
    reciprocal=np.empty((len(selectors),len(cuts)),dtype=float)
    logZS=[]
    for si,S in enumerate(selectors):
      child=np.asarray([float(np.sum(A[np.ix_(S,S)]*d[np.ix_(S,S)])) for d in cuts])
      groups={}
      for di,d in enumerate(cuts): groups.setdefault(restriction_key(d,S),[]).append(di)
      assert len(groups)==2**m and {len(g) for g in groups.values()}=={2**k}
      child_values=[]
      for group in groups.values():
        assert np.max(child[group])-np.min(child[group])<1e-12
        c=child[group[0]]; child_values.append(c)
        external=energy[group]-c
        # Uniform outside completion has exactly zero mean external energy.
        assert abs(float(np.mean(external)))<1e-12
        logK=float(logsumexp(beta*external))
        assert logK+1e-12>=k*math.log(2)
        reciprocal[si,group]=math.exp(k*math.log(2)-logK)
      logZS.append(float(logsumexp(beta*np.asarray(child_values))))
    logZ0=float(logsumexp(logZS)-math.log(len(selectors)))
    C=math.exp(logZA-logZ0)
    direct_f=(C/(2**k))*np.mean(reciprocal,axis=0)
    U=np.mean(reciprocal,axis=0)
    normalized=U/float(np.dot(nu,U))
    assert np.allclose(direct_f,normalized,rtol=2e-12,atol=2e-12)
    assert abs(float(np.dot(nu,normalized))-1)<2e-12
    entropy=float(np.dot(nu*normalized,np.log(normalized)))
    renyi2=math.log(float(np.dot(nu,U*U))/float(np.dot(nu,U))**2)
    assert entropy<=renyi2+1e-12
    return cuts,nu,normalized,U,entropy


def row_squares(A: np.ndarray,cuts: np.ndarray) -> np.ndarray:
    return np.asarray([int(np.sum((A*d).sum(axis=1)**2)) for d in cuts],dtype=float)


def one_deletion_checks(A: np.ndarray,beta: float):
    n=len(A); cuts,nu,f,U,entropy=harmonic_endpoint(A,n-1,beta)
    fields=np.asarray([[(A*d).sum(axis=1)[i] for i in range(n)] for d in cuts],dtype=float)
    usech=np.mean(1/np.cosh(2*beta*fields),axis=1)
    assert np.allclose(U,usech,rtol=2e-12,atol=2e-12)
    R=np.sum(fields*fields,axis=1)
    meanR=float(np.dot(nu,R))
    linear=2*beta*math.sqrt(meanR/n)
    quadratic=2*beta*beta*meanR/n
    spectral=2*beta*float(np.linalg.norm(A,2))
    assert entropy<=linear+2e-12 and entropy<=quadratic+2e-12
    assert linear<=spectral+2e-12
    return {"entropy":entropy,"mean_R2":meanR,"linear_bound":linear,
            "quadratic_bound":quadratic,"spectral_bound":spectral,
            "U_values":sorted(set(np.round(U,14)))}


def j_minus_i_coefficients(n: int):
    A=np.ones((n,n),dtype=np.int64)-np.eye(n,dtype=np.int64)
    cuts=oriented_cuts(n); selectors=[tuple(j for j in range(n) if j!=i) for i in range(n)]
    g=np.empty((n,len(cuts)),dtype=np.int64)
    for di,d in enumerate(cuts):
      # Recover an arbitrary spin representative; squares erase orientation.
      sigma=int(d[0,1]*d[1,2]*d[2,0]); x=np.ones(n,dtype=np.int64)
      x[1:]=d[0,1:]//sigma
      for outside,S in enumerate(selectors):
        b=int(np.dot(A[outside,list(S)],x[list(S)]))
        g[outside,di]=(n-1)-b*b
    gbar=np.mean(g,axis=0)
    parent=2*float(np.mean(gbar*gbar))
    hell=.5*float(np.mean([(g[i]-g[j])@(g[i]-g[j])/len(cuts)
                            for i in range(n) for j in range(n) if i!=j]))
    return parent,hell


def main():
    # General-k identity, not just one deletion.
    harmonic_endpoint(A6,3,0.37)

    a4=one_deletion_checks(A4,0.5)
    assert abs(a4["entropy"]-0.012957183462625588)<2e-13
    assert abs(float(np.linalg.norm(A4,2))-math.sqrt(5))<1e-12

    # A6 is a symmetric conference matrix: centered R2 is identically zero,
    # while the matched endpoint likelihood is nonconstant.
    assert np.array_equal(A6@A6,5*np.eye(6,dtype=np.int64))
    a6=one_deletion_checks(A6,0.5)
    cuts=oriented_cuts(6); R=row_squares(A6,cuts)
    assert set(R)=={30.0}
    ug=(5/math.cosh(1)+1/math.cosh(5))/6
    ue=(1/math.cosh(1)+1/math.cosh(3))/2
    assert len(a6["U_values"])==2
    assert abs(a6["U_values"][0]-min(ug,ue))<2e-13
    assert abs(a6["U_values"][1]-max(ug,ue))<2e-13
    assert abs(a6["entropy"]-0.008845026487227174)<2e-13

    # Reproduce the exact leading-coefficient J-I wall without using a
    # Taylor remainder at growing n.
    n=8; parent,hell=j_minus_i_coefficients(n)
    assert abs(parent-4*(n-1)*(n-2)**2/n)<1e-12
    assert abs(hell-4*(n-2))<1e-12
    assert abs(parent/hell-(n-1)*(n-2)/n)<1e-12

    print("A4",a4)
    print("A6",a6)
    print("J-I n=8 coefficients",parent,hell,parent/hell)
    print("PASS parent_common_mode_r30")


if __name__=="__main__": main()
