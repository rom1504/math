"""Exact finite Palm-switch and light-variance checks for signed matchings."""

import itertools
import json

import numpy as np

from paper_discrepancy_2026_09_17_involution_ground import matchings


rng = np.random.default_rng(2026091720)
reports = []
for p in (6,8,10):
    edges = list(itertools.combinations(range(p),2))
    index = {e:j for j,e in enumerate(edges)}
    mm = list(matchings(tuple(range(p))))
    match_edges = np.array([[index[tuple(sorted(e))] for e in pairs] for pairs in mm])
    indicators = np.zeros((len(mm),len(edges)),dtype=np.int64)
    for row,idx in zip(indicators,match_edges):
        row[idx]=1
    partners = np.empty((len(mm),p),dtype=np.int64)
    for row,pairs in zip(partners,mm):
        for a,b in pairs:
            row[a],row[b]=b,a
    lookup = {tuple(row):j for j,row in enumerate(indicators)}
    switches=[]
    for ei,(a,b) in enumerate(edges):
        conditional=[]
        for mi in range(len(mm)):
            row=indicators[mi].copy()
            if not row[ei]:
                u,w=int(partners[mi,a]),int(partners[mi,b])
                row[index[tuple(sorted((a,u)))]]=0
                row[index[tuple(sorted((b,w)))]]=0
                row[ei]=1
                row[index[tuple(sorted((u,w)))]]=1
            conditional.append(lookup[tuple(row)])
        conditional=np.array(conditional)
        target_count=np.bincount(conditional,minlength=len(mm))
        assert np.array_equal(target_count,(p-1)*indicators[:,ei])
        switches.append(conditional)
    weighted_checks=0
    variance_checks=0
    max_weighted_ratio=0.0
    max_variance_ratio=0.0
    for trial in range(40):
        c=rng.normal(size=len(edges))/np.sqrt(p)
        if trial%4==0:
            c*=rng.random(len(edges))<0.2
        degrees=np.zeros(p)
        for (a,b),cc in zip(edges,c):
            degrees[a]+=cc*cc
            degrees[b]+=cc*cc
        d=float(degrees.max())
        for tau in (0.1,0.25,0.5,0.8):
            heavy=np.abs(c)>tau
            light=~heavy
            vm=indicators@(c*c*light)
            vl=float(vm.mean())
            variance=float(vm.var())
            variance_bound=tau*tau*vl+2*vl*vl/(p-3)
            assert variance<=variance_bound+1e-12
            if variance_bound:
                max_variance_ratio=max(max_variance_ratio,variance/variance_bound)
            variance_checks+=1
            total_cost=0.0
            for ei in np.flatnonzero(heavy):
                reduced=indicators[switches[ei]].copy()
                reduced[:,ei]-=1
                cost=np.abs(indicators-reduced)@(np.abs(c)*heavy)
                exact=float(cost.mean())
                bound=abs(c[ei])/(p-1)+3*d/(tau*(p-3))
                assert exact<=bound+1e-12
                total_cost+=exact/(p-1)
                weighted_checks+=1
            final_bound=(p*d/(2*tau*(p-1)**2)
                         +3*p*d*d/(2*tau**3*(p-1)*(p-3)))
            assert total_cost<=final_bound+1e-12
            if final_bound:
                max_weighted_ratio=max(max_weighted_ratio,total_cost/final_bound)
    reports.append(dict(p=p,matchings=len(mm),exact_palm_edges=len(edges),
                        weighted_switch_checks=weighted_checks,
                        exact_light_variance_checks=variance_checks,
                        maximum_weighted_bound_ratio=max_weighted_ratio,
                        maximum_variance_bound_ratio=max_variance_ratio))
print(json.dumps(dict(status='PASS',seed=2026091720,cases=reports,
    assertions=['each conditional matching has exactly p-1 preimages',
                'exact light variance bounded by tau^2 v+2v^2/(p-3)',
                'weighted reduced-Palm switch bound',
                'finite total weighted Poisson error bound']),indent=2))
