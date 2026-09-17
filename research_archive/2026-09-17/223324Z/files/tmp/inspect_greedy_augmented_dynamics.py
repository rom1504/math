#!/usr/bin/env python3
"""Inspect best-improvement dynamics for one row-sign collapsed branch."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
ROOT=Path('/home/math/quadra')
CASES=[('conference_double_p17.json',32768),('conference_order26_gf25.json',8192),('conference_double_p97.json',4096)]
rng=np.random.default_rng(260813)

def ascent(D,h,sigma):
 r=np.where(sigma*h>=0,1,-1).astype(np.int64);t=1; counts=np.zeros(len(r),int);tf=0; gains=[]
 while True:
  L=D@r+t*h
  gr=-4*sigma*r*L
  gt=-4*sigma*t*(h@r)
  idx=int(np.argmax(gr)) if len(r) else -1
  bg=int(gr[idx]) if len(r) else -10**9
  if gt>bg and gt>0: t=-t;tf+=1;gains.append(int(gt))
  elif bg>0: r[idx]=-r[idx];counts[idx]+=1;gains.append(bg)
  else:break
  if len(gains)>10000: raise RuntimeError
 return int(sigma*(r@D@r+2*t*(h@r))),counts,tf,gains,r,t

for fn,N in CASES:
 C=np.array(json.load(open(ROOT/'computations/results'/fn))['conference_matrix'],int);n=len(C)
 if N==(1<<(n-1)):
  codes=np.arange(N,dtype=np.uint64)[:,None];bits=((codes>>np.arange(n-1,dtype=np.uint64))&1).astype(int);XX=np.c_[np.ones(N,int),1-2*bits]
 else:XX=rng.choice([-1,1],size=(N,n));XX[:,0]=1
 agg={k:[] for k in ['loss','steps','repeats','maxcount','tflips','initial','target','sumgain','j','R','H','h2']}
 for X in XX:
  B=C*X[:,None]*X[None,:];ell=B@np.ones(n,int);I=np.flatnonzero(ell>0);J=np.flatnonzero(ell<0);D=B[np.ix_(J,J)];h=B[np.ix_(J,I)]@np.ones(len(I),int);P=int(np.ones(len(I),int)@B[np.ix_(I,I)]@np.ones(len(I),int));R=int(np.ones(len(J),int)@D@np.ones(len(J),int));sigma=1 if P>=0 else -1
  initr=np.where(sigma*h>=0,1,-1);init=int(sigma*(initr@D@initr+2*(h@initr)));val,c,tf,g,rr,tt=ascent(D,h,sigma);target=-sigma*R
  vals=dict(loss=max(0,target-val),steps=len(g),repeats=int(c.sum()-(c>0).sum()),maxcount=int(c.max(initial=0)),tflips=tf,initial=init,target=target,sumgain=sum(g),j=len(J),R=R,H=int(h.sum()),h2=int(h@h))
  for k,v in vals.items():agg[k].append(v)
 out={'file':fn,'n':n,'N':N,'stats':{k:{'mean':float(np.mean(v)),'max':int(np.max(v)),'p99':float(np.quantile(v,.99))} for k,v in agg.items()},'loss_hist':dict(zip(*[a.tolist() for a in np.unique(agg['loss'],return_counts=True)]))}
 print(json.dumps(out,sort_keys=True),flush=True)
