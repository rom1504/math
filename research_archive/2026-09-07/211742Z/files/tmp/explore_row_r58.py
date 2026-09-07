#!/usr/bin/env python3
import itertools, math, sys
import numpy as np
sys.path.insert(0, '/home/math/quadra/tmp')
from envelope_block_cover_r27 import A6,A8,A9,projective_spins

def qnorm(a): return max(abs(int(x@a@x)) for x in projective_spins(len(a)))

def records(a,m,t=0):
 n=len(a); q=qnorm(a); sels=list(itertools.combinations(range(n),m))
 qs=[qnorm(a[np.ix_(s,s)]) for s in sels]
 p2=m*(m-1)/(n*(n-1)); B=((m/n)**1.5-p2)*q
 out=[]
 for xi,x in enumerate(projective_spins(n)):
  raw=int(x@a@x); row=int((a@x)@(a@x))
  for tau in (-1,1):
   deficit=q-tau*raw; good=[]
   for s,qq in zip(sels,qs):
    loc=tau*int(x[list(s)]@a[np.ix_(s,s)]@x[list(s)])
    good.append(qq-loc-p2*deficit-B<=t+1e-12)
   u=sum(good)/len(sels)
   if u:
    fields=tau*x*(a@x)
    out.append(dict(x=x.copy(),xi=xi,tau=tau,u=u,h=-math.log(u),R=row,D=deficit,
                    E=q-deficit,fields=fields,good=np.array(good),sels=sels,B=B,p2=p2))
 return out

for name,a in [('A6',A6),('A8',A8),('A9',A9)]:
 for m in range(2,len(a)):
  rr=records(a,m)
  pareto=[]
  for z in rr:
   if not any(w['h']<=z['h']+1e-12 and w['R']<=z['R'] and
              (w['h']<z['h']-1e-12 or w['R']<z['R']) for w in rr):
    pareto.append(z)
  types=sorted(set((round(z['u'],8),z['R'],z['D'], int(np.maximum(z['fields'],0).sum())) for z in pareto))
  correct=[x for x in types if x[2] <= qnorm(a)]
  high=[x for x in correct if x[3]>=2]
  if high:
   print(name,m,'pareto',types,'CORRECT HIGH',high)

# Detailed A8 exposed high-row correct-sector column and possible responses.
a=A8; m=6; rr=records(a,m); lam=0.001
vals=np.array([z['h']+lam*z['R'] for z in rr])
v=vals.min(); opt=[z for z,val in zip(rr,vals) if abs(val-v)<1e-12]
print('A8m6 lambda',lam,'opt count',len(opt),'types',sorted(set((z['u'],z['R'],z['D']) for z in opt)), 'nextgap', min(val-v for val in vals if val>v+1e-12))
d=next(z for z in opt)
print('d',d['xi'],d['tau'],'x',d['x'].tolist(),'fields',d['fields'].tolist())
pos=[]
for i in range(len(a)):
 for j in range(i+1,len(a)):
  s=d['tau']*a[i,j]*d['x'][i]*d['x'][j]
  if s==1: pos.append((i,j))
grounds=[z for z in rr if z['D']==0]
for ecount in range(2,8):
 for E in itertools.combinations(pos,ecount):
  for w in grounds:
   if w is d: continue
   rev=sum(w['tau']*a[i,j]*w['x'][i]*w['x'][j]==-1 for i,j in E)
   if rev>ecount/4:
    # relative D and exact identities
    kappa=d['tau']*w['tau']; y=d['x']*w['x']
    D=[(i,j) for i in range(len(a)) for j in range(i+1,len(a)) if kappa*y[i]*y[j]==-1]
    sd={(i,j):int(d['tau']*a[i,j]*d['x'][i]*d['x'][j]) for i in range(len(a)) for j in range(i+1,len(a))}
    b=np.array([sum(sd[min(i,j),max(i,j)] for j in range(len(a)) if j!=i and (min(i,j),max(i,j)) in D) for i in range(len(a))])
    print('E',E,'response',w['xi'],w['tau'],'x',w['x'].tolist(),'fields',w['fields'].tolist(),'uRD',w['u'],w['R'],w['D'],'rev',rev,'Dset',D,'b',b.tolist(),'W',sum(sd[e] for e in D),'Rdiff',w['R']-d['R'],'formula',4*int(b@b-d['fields']@b),'objdiff',w['h']+lam*w['R']-(d['h']+lam*d['R']))
    print('fibres counts d,w,inter,onlyd,onlyw',int(d['good'].sum()),int(w['good'].sum()),int((d['good']&w['good']).sum()),int((d['good']&~w['good']).sum()),int((~d['good']&w['good']).sum()))
    raise SystemExit
