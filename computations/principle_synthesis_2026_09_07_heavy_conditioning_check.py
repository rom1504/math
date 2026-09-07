#!/usr/bin/env python3
"""Finite regressions for descending matching and tagged-zero conditioning."""
import collections
import itertools
import json
import math

import numpy as np


def tag_check(d, h):
    counts = collections.Counter()
    targets = list(range(d-h,d))
    for permutation in itertools.permutations(range(d)):
        word = list(permutation)
        for tag, target in enumerate(targets):
            current = word.index(tag)
            word[current], word[target] = word[target], word[current]
        assert all(word[target] == tag for tag,target in enumerate(targets))
        counts[tuple(word)] += 1
    assert len(counts) == math.factorial(d-h)
    assert len(set(counts.values())) == 1
    assert next(iter(counts.values())) == math.factorial(d)//math.factorial(d-h)
    return dict(d=d,h=h,outputs=len(counts),preimages_per_output=next(iter(counts.values())),
                status="PASS")


def matching_check(n, seed):
    rng = np.random.default_rng(seed)
    h = np.ones((1,1),dtype=np.int64)
    while len(h)<n:
        h=np.block([[h,h],[h,-h]])
    s=h.copy()
    np.fill_diagonal(s,0)
    # Q(S)<=n*(sqrt(n)+1)/2 by this exact Hadamard operator witness.
    root_upper=math.isqrt(n)
    if root_upper*root_upper<n:
        root_upper+=1
    cap_upper=n*(root_upper+1)//2
    c=2
    lower=[4,5,6,7]
    stubs=[]
    for i in range(n):
        row=[]
        budget=c*n
        for _ in range(int(rng.integers(1,7))):
            a=int(rng.choice(lower))
            if a*a<=budget:
                row.append((a,int(rng.choice([-1,1]))))
                budget-=a*a
        stubs.append(row)
    used=set()
    details=[]
    for a in reversed(lower):
        residual=[[sign for value,sign in row if value==a] for row in stubs]
        matched=[]
        while True:
            found=None
            for i in range(n):
                for j in range(i+1,n):
                    if (i,j) in used:
                        continue
                    for si in residual[i]:
                        for sj in residual[j]:
                            if int(s[i,j])==si*sj:
                                found=(i,j,si,sj)
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found is None:
                break
            i,j,si,sj=found
            residual[i].remove(si)
            residual[j].remove(sj)
            used.add((i,j))
            matched.append(found)
        vertices=[i for i in range(n) if residual[i]]
        sign={i:residual[i][0] for i in vertices}
        energy=sum(int(s[i,j])*sign[i]*sign[j]
                   for i,j in itertools.combinations(vertices,2))
        degrees=[sum(i in edge for edge in used) for i in range(n)]
        assert max(degrees,default=0)<=c*n/(a*a)+1e-12
        for i,j in itertools.combinations(vertices,2):
            if (i,j) not in used:
                assert int(s[i,j])*sign[i]*sign[j]==-1
        upper=2*c*n/(a*a)+1+math.sqrt(2*cap_upper)
        assert len(vertices)<=upper+1e-12
        assert energy<=len(vertices)*c*n/(a*a)-math.comb(len(vertices),2)+1e-12
        details.append(dict(a=a,matched_edges=len(matched),
                            residual_vertices=len(vertices),residual_energy=energy,
                            residual_vertex_upper=upper))
    return dict(n=n,seed=seed,cap_upper=cap_upper,bins=details,status="PASS")


def main():
    output=dict(status="PASS",scope="finite regression, not proof of asymptotic comparison",
                tags=[tag_check(d,h) for d in range(3,8) for h in range(1,min(3,d)+1)],
                matchings=[matching_check(n,20260907+j) for j,n in enumerate([16,32,64])])
    print(json.dumps(output,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
