"""Exact all-flip verification of anchors, ties, and conditioning.

No numerical probability approximation is used. This checks the finite
identities underlying the actual-sign neighborhood/path theorem, not the
asymptotic chaining theorem or its universal constants.
"""
from fractions import Fraction
from itertools import combinations, product
from math import comb
import json
import numpy as np


def run(n):
    edges = list(combinations(range(n), 2))
    d = len(edges)
    spins = np.array(list(product((-1, 1), repeat=n)), dtype=np.int64)
    features = np.array([spins[:,i]*spins[:,j] for i,j in edges],dtype=np.int64).T
    # Fixed actual full signing; no optimum claim is attached to it.
    signs = np.array([1 if (i+2*j)%4 else -1 for i,j in edges],dtype=np.int64)
    old = features@signs
    cap = int(np.abs(old).max())
    anchor_index = int(np.argmax(np.abs(old)))
    anchor_sign = 1 if old[anchor_index] >= 0 else -1
    signed_features = np.concatenate((features,-features),axis=0)
    signed_old = np.concatenate((old,-old))
    deficits = cap-signed_old
    bitrows = np.array(list(product((0,1),repeat=d)),dtype=np.int64)
    flips = bitrows.sum(axis=1)
    energies = ((1-2*bitrows)*signs)@signed_features.T
    anchor_energy = ((1-2*bitrows)*signs)@(anchor_sign*features[anchor_index])
    new_caps = energies.max(axis=1)
    tests = []
    for tolerance in [0,1,2,4]:
        outside = deficits > tolerance
        if not outside.any():
            continue
        escaping = (energies[:,outside] == new_caps[:,None]).any(axis=1)
        # rho=1/8; F is multiplied by 2*8, exactly.
        numerator, denominator = 1,8
        f_scaled = np.maximum(0,np.max(
            2*denominator*(energies[:,outside]-anchor_energy[:,None])
            +(denominator-2*numerator)*deficits[outside][None,:],axis=1))
        assert np.all(4*f_scaled[escaping] >= 2*denominator*tolerance)
        # Uniform k-subset law = independent Bernoulli(k/d) conditioned
        # on its modal count k. Check each conditional event exactly.
        for k in range(d+1):
            total_k = comb(d,k)
            failures_k = int(escaping[flips==k].sum())
            assert int((flips==k).sum())==total_k
            if k in (0,d):
                mode_probability = Fraction(1)
                bernoulli_failure = Fraction(failures_k,total_k)
            else:
                rho = Fraction(k,d)
                mode_probability = total_k*rho**k*(1-rho)**(d-k)
                bernoulli_failure = sum(
                    Fraction(int(escaping[flips==j].sum()))*rho**j*(1-rho)**(d-j)
                    for j in range(d+1))
            assert mode_probability >= Fraction(1,d+1)
            assert Fraction(failures_k,total_k) <= (d+1)*bernoulli_failure
        tests.append({"window":tolerance,"escaping_flip_sets":int(escaping.sum()),
                      "all_flip_sets":len(bitrows),"conditioning_counts_checked":d+1})
    return {"n":n,"edges":d,"actual_cap":cap,"tests":tests}


if __name__ == "__main__":
    results=[run(n) for n in (4,5,6)]
    print(json.dumps({"status":"PASS exact full-sign anchors and modal conditioning",
                      "results":results},indent=2))
