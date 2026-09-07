#!/usr/bin/env python3
import itertools,sys,numpy as np
sys.path.insert(0,'/home/math/quadra/tmp')
from envelope_block_cover_r27 import A9

def spins(n):
 for t in itertools.product((-1,1),repeat=n-1): yield np.array((1,)+t,dtype=int)
def grounds(B):
 z=[(x,int(x@B@x)) for x in spins(len(B))]; q=max(abs(v) for x,v in z)
 return q,[(1 if v>0 else -1,x) for x,v in z if abs(v)==q]
def all_unions(B,sig,y):
 m=len(B); edges=[(i,j) for i in range(m) for j in range(i+1,m)]
 idx={e:k for k,e in enumerate(edges)}
 c={(i,j):sig*y[i]*B[i,j]*y[j] for i,j in edges}
 r=sig*y*(B@y)
 masks={0}
 for i in range(m):
  pos=[]
  for j in range(m):
   if i==j: continue
   e=tuple(sorted((i,j)))
   if c[e]==1: pos.append(idx[e])
  choices=[]
  for cc in itertools.combinations(pos,int(r[i])):
   mask=sum(1<<k for k in cc); choices.append(mask)
  masks={u|v for u in masks for v in choices}
 return edges,c,masks

A=A9;n=9;m=7;x=np.array([1,1,-1,1,-1,-1,1,-1,1]);tau=1
data=[]
for S0 in itertools.combinations(range(n),m):
 S=np.array(S0);B=A[np.ix_(S,S)];q,cgs=grounds(B)
 opts=[]
 for sig,y in cgs:
  edges,c,masks=all_unions(B,sig,y)
  opts.append((sig,y,edges,masks))
 data.append((S0,S,B,q,cgs,opts))

def robust_for(tau,x):
 robust=[]
 for S0,S,B,q,cgs,opts in data:
  local=tau*int(x[S]@B@x[S]);delta=q-local
  best=-2; worst=2; nmasks=0
  for sig,y,edges,masks in opts:
   nmasks+=len(masks)
   for mask in masks:
    vals=[tau*A[S[i],S[j]]*x[S[i]]*x[S[j]] for h,(i,j) in enumerate(edges) if mask>>h&1]
    ratio=sum(vals)/len(vals)
    best=max(best,ratio);worst=min(worst,ratio)
  robust.append((best<=.5,worst<=.5,best,worst,delta,len(cgs),nmasks,S0))
 return robust
robust=robust_for(tau,x)
print('robust escape all choices',sum(a for a,*_ in robust),'some choice',sum(b for a,b,*_ in robust))
print('robust unfavorable escape',sum(a and d>=4 for a,b,best,worst,d,*_ in robust))
for z in robust: print(z)

q,pg=grounds(A)
scores=[]
for tau,x in pg:
 rr=robust_for(tau,x)
 scores.append((sum(a for a,*_ in rr),sum(a and z[4]>=4 for z in rr for a in [z[0]]),tau,''.join('+' if u==1 else '-' for u in x)))
print('all ground robust scores')
print(sorted(scores,reverse=True))
