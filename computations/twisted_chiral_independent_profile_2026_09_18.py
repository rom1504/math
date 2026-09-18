#!/usr/bin/env python3
"""Independent direct spin enumeration of the fixed-child family for small r.

Unlike the search engine, constructs all parent energy coefficient vectors
with straightforward matrix multiplication; computes every matching cap in
batches. This deliberately does not use Gray updates or chiral reduction.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--orders',type=int,nargs='+',default=[3,4,5,6])
    parser.add_argument('--random-count',type=int,default=20)
    parser.add_argument('--max-exhaustive',type=int,default=6)
    parser.add_argument('--classes',type=int,nargs='+')
    parser.add_argument('--tag',default='independent')
    args=parser.parse_args()
    from twisted_chiral_driver_2026_09_18 import load_seeds
    records=[]
    for r in args.orders:
        for label,a,source in load_seeds(r):
            if args.classes is not None and int(label.split('class')[-1]) not in args.classes:continue
            started=time.time();a=np.asarray(a,dtype=np.int16)
            codes=np.arange(1<<r,dtype=np.uint64)
            spins=(1-2*((codes[:,None]>>np.arange(r,dtype=np.uint64))&1).astype(np.int16))
            projective=spins[::2]
            q=np.einsum('bi,ij,bj->b',spins,a,spins)//2
            # All x_0=+1, all y: no use of chirality.
            base=q[::2,None]-q[None,:]
            t=projective[:,None,:]*spins[None,:,:]
            matching=np.einsum('xyi,di->xyd',t,spins)
            seen=set();hist=Counter();core_hist=Counter();width_hist=Counter();best=10000;witness=None
            rng=np.random.default_rng(20260918+r)
            if r<=args.max_exhaustive:
                params=((p,s) for p in itertools.permutations(range(r)) for s in projective)
                exhaustive=True
            else:
                params=((rng.permutation(r),projective[rng.integers(len(projective))]) for _ in range(args.random_count))
                exhaustive=False
            for p,s in params:
                b=s[:,None]*a[np.ix_(p,p)]*s[None,:]
                key=b.tobytes()
                if key in seen:continue
                seen.add(key)
                ebase=base+projective@b@spins.T
                core_hist[int(np.max(np.abs(ebase)))]+=1
                index=np.arange(len(projective))
                by_t=ebase[index[:,None],2*(index[:,None]^index[None,:])]
                width_hist[int(np.max(np.max(by_t,axis=0)-np.min(by_t,axis=0)))//2]+=1
                caps=np.max(np.abs(ebase[:,:,None]+matching),axis=(0,1))
                cap=int(caps.min());hist[cap]+=1
                if cap<best:
                    best=cap;didx=int(np.argmin(caps))
                    witness=dict(p=list(map(int,p)),s=s.tolist(),d=spins[didx].tolist())
            record=dict(label=label,r=r,source=source,exhaustive=exhaustive,
                best=best,distinct_B=len(seen),best_matching_cap_histogram=dict(hist),
                core_cap_histogram=dict(core_hist),unconstrained_profile_center_radius_histogram=dict(width_hist),
                witness=witness,elapsed=time.time()-started,method='all x0=+1, all y, all d by direct integer matrix multiplication')
            records.append(record);print(json.dumps(record),flush=True)
            (ROOT/f'computations/results/twisted_chiral_2026_09_18_{args.tag}.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':main()
