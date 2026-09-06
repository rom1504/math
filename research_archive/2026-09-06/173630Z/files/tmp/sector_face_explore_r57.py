#!/usr/bin/env python3
import itertools, math, sys
import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A8,A9
from box_triple_r50_check import A as A10

def spins(n):
    for tail in itertools.product((-1,1), repeat=n-1):
        yield np.array((1,)+tail,dtype=int)

def capgrounds(B):
    vals=[(x,int(x@B@x)) for x in spins(len(B))]
    q=max(abs(v) for x,v in vals)
    gs=[]
    for x,v in vals:
        if abs(v)==q:
            sig=1 if v>0 else -1
            gs.append((sig,x.copy()))
    return q,gs

def edgeblock(A,S,sig,y):
    # Wave-56 all-vertex construction, deterministic first r_i positives.
    B=A[np.ix_(S,S)]
    r=sig*y*(B@y)
    E=set()
    for ii,i in enumerate(S):
        pos=[]
        for jj,j in enumerate(S):
            if ii!=jj and sig*y[ii]*B[ii,jj]*y[jj]==1:
                pos.append((int(i),int(j)))
        assert len(pos)==(len(S)-1+int(r[ii]))//2
        for i0,j0 in pos[:int(r[ii])]: E.add(tuple(sorted((i0,j0))))
    return sorted(E),r

def audit(A,m,name):
    n=len(A); q,pg=capgrounds(A)
    rec=[]
    for S0 in itertools.combinations(range(n),m):
        S=np.array(S0,dtype=int); qs,cg=capgrounds(A[np.ix_(S,S)])
        # lexicographic first oriented child ground
        sig,y=cg[0]
        E,r=edgeblock(A,S,sig,y)
        rec.append((S,qs,sig,y,E,r))
    p=m/n; p2=m*(m-1)/(n*(n-1)); Bconst=(p**1.5-p2)*q
    print(name,"n,m,q,#grounds,#S,B",n,m,q,len(pg),len(rec),Bconst)
    rows=[]
    for tau,x in pg:
        esc=fav=escfav=near=0
        vals=[]
        for S,qs,sig,y,E,r in rec:
            SE=sum(tau*A[i,j]*x[i]*x[j] for i,j in E)
            isesc=SE<=len(E)/2
            local=tau*int(x[S]@A[np.ix_(S,S)]@x[S])
            delta=qs-local
            what=delta-Bconst # Delta zero
            isfav=what<=0
            esc+=isesc; fav+=isfav; escfav+=isesc and isfav
            if isesc: vals.append((what,delta,SE,len(E),tuple(S)))
        rows.append((esc,fav,escfav,tau,''.join('+' if z==1 else '-' for z in x), vals))
    rows.sort(reverse=True)
    for row in rows[:8]:
        esc,fav,ef,tau,xstr,vals=row
        print(" state",tau,xstr,"esc/fav/ef",esc,fav,ef,"esc what range", (min(v[0] for v in vals),max(v[0] for v in vals)) if vals else None)
    # Best escape with worst intersection.
    candidates=[r for r in rows if r[0]>=len(rec)/3]
    print(">=third count",len(candidates),"min escfav",min((r[2] for r in candidates),default=None))

    if name=='A9' and m==7:
        hit=next(r for r in rows if r[3]==-1 and r[4]=='+-++-----')
        esc,fav,ef,tau,xstr,vals=hit
        x=np.array([1 if c=='+' else -1 for c in xstr],dtype=int)
        detail=[]
        for S,qs,sig,y,E,r in rec:
            SE=sum(tau*A[i,j]*x[i]*x[j] for i,j in E)
            local=tau*int(x[S]@A[np.ix_(S,S)]@x[S])
            delta=qs-local
            if SE<=len(E)/2:
                z=x[S]*y; kappa=tau*sig
                cut=sum(1 for i,j in E if z[list(S).index(i)]!=z[list(S).index(j)])
                detail.append((sig,kappa,qs,delta,SE,len(E),cut,tuple(S)))
        print(' obstruction detail full_energy,row,total_s',int(tau*x@A@x),int((A@x)@(A@x)),sum(tau*A[i,j]*x[i]*x[j] for i in range(n) for j in range(i+1,n)))
        print(' orientation counts', {k:sum(d[1]==k for d in detail) for k in (-1,1)})
        print(' delta values', {k:sum(d[3]==k for d in detail) for k in sorted(set(d[3] for d in detail))})
        print(' block sizes', {k:sum(d[5]==k for d in detail) for k in sorted(set(d[5] for d in detail))})
        print(' first obstruction rows',detail[:8])

for A,name in ((A8,'A8'),(A9,'A9'),(A10,'A10')):
  for m in range(max(3,len(A)-3),len(A)):
    audit(A,m,name)
