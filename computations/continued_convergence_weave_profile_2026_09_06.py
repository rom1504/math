"""Positive-DP evaluation of the one-row soft Finner bound (not a cap).

Each permanent is divided by (m-1)!, so the recursion is an expectation
over a uniform permutation and involves only nonnegative terms.
"""
import argparse
from collections import Counter
from functools import lru_cache
import json
import math
from pathlib import Path
import numpy as np
from continued_convergence_restricted_weave_2026_09_06 import walsh

def profiles(m,k,seed):
    rng=np.random.default_rng(seed); rows=rng.choice(m,k,replace=False)
    a=walsh(m)[rows]
    counts=Counter()
    # Global spin reversal has the same absolute spectrum; count it twice.
    for start in range(0,1<<(k-1),1<<15):
        ids=np.arange(start,min(1<<(k-1),start+(1<<15)),dtype=np.uint64)
        x=np.ones((len(ids),k),dtype=np.int64)
        x[:,1:]=1-2*((ids[:,None]>>np.arange(k-1,dtype=np.uint64))&1).astype(np.int64)
        spectra=np.abs(x@a)
        for spectrum in spectra:counts[tuple(sorted(map(int,spectrum)))]+=2
    return rows,counts

def permanent_expectation(a,k,t):
    values=tuple(sorted(set(a))); multiplicities=tuple(a.count(v) for v in values)
    rows=tuple(sorted(a)); n=len(rows)
    kernel=np.array([[.5*(math.exp(-t*(v-w)**2/k)+math.exp(-t*(v+w)**2/k)) for w in values] for v in rows])
    @lru_cache(None)
    def rec(left):
        remaining=sum(left)
        if not remaining:return 1.0
        i=n-remaining;ans=0.0
        for j,c in enumerate(left):
            if c:
                nxt=list(left);nxt[j]-=1
                ans+=(c/remaining)*kernel[i,j]*rec(tuple(nxt))
        return ans
    return rec(multiplicities)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--m',type=int,default=16)
    ap.add_argument('--k',type=int,default=14);ap.add_argument('--seed',type=int,default=9061001)
    ap.add_argument('--times',default='.25,.5,1,2,4,8,16,32,64,128')
    ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    rows,counts=profiles(args.m,args.k,args.seed);records=[]
    eta=1-math.sqrt(args.k/args.m)
    for t in map(float,args.times.split(',')):
        entries=[]
        for a,count in counts.items():
            options=[]
            for v in set(a):
                shortened=list(a);shortened.remove(v)
                options.append((permanent_expectation(shortened,args.k,t),v))
            expectation,removed=max(options)
            term=count*math.sqrt(expectation)
            entries.append({'profile':a,'count':count,'removed':removed,'L':math.sqrt(expectation),'contribution':term})
        z=sum(e['contribution'] for e in entries)
        records.append({'t':t,'Z':z,'log_Z_over_m':math.log(z)/args.m,
          'first_moment_exponent_per_m2_at_half':math.log(z)/args.m+t*eta,
          'largest_profiles':sorted(entries,key=lambda e:e['contribution'],reverse=True)[:8]})
        print(json.dumps({q:records[-1][q] for q in ['t','Z','log_Z_over_m','first_moment_exponent_per_m2_at_half']}),flush=True)
    args.output.write_text(json.dumps({'m':args.m,'k':args.k,'seed':args.seed,'rows':rows.tolist(),
       'number_of_profiles':len(counts),'total_spin_count':sum(counts.values()),'eta_at_half':eta,
       'records':records,'status':'finite sufficient-bound evaluation; positive exponent is inconclusive'},indent=2)+'\n')

if __name__=='__main__':main()
