"""Exact finite deleted-slice tests and actual full-sign escape diagnostics.

The asymptotic simultaneous-existence theorem is analytical. Monte Carlo
entries below are labeled diagnostics, not certificates of limiting caps.
"""
from __future__ import annotations

from fractions import Fraction
import itertools
import json
from math import comb, sqrt, pi
from pathlib import Path

import numpy as np


def sample_blocks(m: int, ell: int, rng):
    a = Fraction(m*(m-2), ell*(ell-2))
    upper = np.triu(np.where(rng.random((ell,ell)) < (1-float(a))/2, 1, -1), 1)
    bulk = (upper+upper.T).astype(np.int64)
    cross = rng.choice((-1,1), size=(m,ell)).astype(np.int64)
    return a, bulk, cross


def greedy(bulk):
    ell = len(bulk)
    d = bulk[::2,::2]-bulk[::2,1::2]-bulk[1::2,::2]+bulk[1::2,1::2]
    pairs = ell//2
    left = pairs//2
    field = d[:left,left:].sum(axis=0)
    pair_signs = np.concatenate((np.ones(left,dtype=np.int64), np.where(field>=0,1,-1)))
    y = np.empty(ell,dtype=np.int64)
    y[::2],y[1::2] = pair_signs,-pair_signs
    assert y.sum() == 0
    actual = int(y@bulk@y)//2
    pair_energy = sum(int(d[i,j]*pair_signs[i]*pair_signs[j]) for i in range(pairs) for j in range(i+1,pairs))
    within = -int(bulk[np.arange(0,ell,2),np.arange(1,ell,2)].sum())
    assert actual == pair_energy+within
    return y, actual, int(np.abs(field).sum())


def main():
    rng = np.random.default_rng(271820260917)
    exact_response_checks = 0
    exact_parent_caps = []
    for m,ell in ((4,6),(4,8),(6,8),(6,10)):
        bulk_words = np.asarray(list(itertools.product((-1,1),repeat=ell)),dtype=np.int64)
        balanced = bulk_words[bulk_words.sum(axis=1)==0]
        r=ell//2
        f=Fraction(4*r*comb(r-1,r//2)**2,comb(2*r,r))
        for _ in range(3):
            a,bulk,cross=sample_blocks(m,ell,rng)
            q0=Fraction(m*(m-1),2)+a*ell/2
            eb=np.einsum("bi,ij,bj->b",balanced,bulk,balanced)//2
            field=balanced@cross.sum(axis=0)
            eplus=m*(m-1)//2+eb+field
            eminus=m*(m-1)//2+eb-field
            # Select an actual symmetric slice, using a finite tolerance.
            deviation=np.maximum(np.abs(eplus-float(q0)),np.abs(eminus-float(q0)))
            tolerance=float(np.quantile(deviation,.8))+1e-9
            good=balanced[deviation<=tolerance]
            assert len(good)>0
            missing=Fraction(len(balanced)-len(good),len(balanced))
            error=missing*ell/(1-missing)
            for h in bulk_words:
                response=Fraction(int(np.abs(good@h).sum()),len(good))
                assert response >= f*(1-Fraction(abs(int(h.sum())),ell))-error
                exact_response_checks+=1
            core=np.ones((m,m),dtype=np.int64)-np.eye(m,dtype=np.int64)
            full=np.block([[core,cross],[cross.T,bulk]])
            assert np.all(np.diag(full)==0)
            assert np.all(np.abs(full+np.eye(m+ell,dtype=np.int64))==1)
            words=np.asarray(list(itertools.product((-1,1),repeat=m+ell)),dtype=np.int64)
            energies=np.einsum("bi,ij,bj->b",words,full,words)//2
            cap=int(np.abs(energies).max())
            y, bulk_greedy, cross_gain=greedy(bulk)
            actual_witness=m*(m-1)//2+bulk_greedy+int((cross@y).sum())
            assert cap>=actual_witness
            exact_parent_caps.append({"n":m+ell,"m":m,"cap":cap,"mean_cap":str(q0),
                                      "retained_balanced_fraction":float(1-missing)})
    diagnostics=[]
    for n in (128,256,512,1024):
        m=2*int(n**.75/2)
        ell=n-m
        cross_values=[]
        bulk_values=[]
        for _ in range(8):
            a,bulk,cross=sample_blocks(m,ell,rng)
            y,energy,gain=greedy(bulk)
            cross_values.append(gain/ell**1.5)
            bulk_values.append(energy/ell**1.5)
        diagnostics.append({"n":n,"m":m,"ell":ell,
                            "mean_pair_cross_coefficient":float(np.mean(cross_values)),
                            "mean_bulk_witness_coefficient":float(np.mean(bulk_values)),
                            "asymptotic_target":sqrt(2/pi)/4})
    result={"status":"PASS","exact_deleted_slice_response_checks":exact_response_checks,
            "exact_full_parent_cases":len(exact_parent_caps),"finite_parents":exact_parent_caps,
            "monte_carlo_greedy_diagnostics_not_proof":diagnostics,
            "scope":"actual high-energy counterexample, explicitly NOTground"}
    target=Path(__file__).resolve().parents[1]/"tmp/paper_portfolio_2026_09_17/bernoulli/actual_high_energy_obstruction.json"
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
