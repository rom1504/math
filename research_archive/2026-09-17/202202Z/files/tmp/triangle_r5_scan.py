from itertools import product, combinations
import numpy as np

A5 = np.array([
[0,-1,1,-1,1],
[-1,0,-1,1,1],
[1,-1,0,1,1],
[-1,1,1,0,1],
[1,1,1,1,0]], dtype=int)

A9 = np.array([
[0,1,1,-1,1,1,-1,1,1],
[1,0,1,1,1,1,1,-1,1],
[1,1,0,-1,1,-1,1,-1,-1],
[-1,1,-1,0,1,-1,1,-1,1],
[1,1,1,1,0,-1,1,1,1],
[1,1,-1,-1,-1,0,1,1,1],
[-1,1,1,1,1,1,0,1,-1],
[1,-1,-1,-1,1,1,1,0,-1],
[1,1,-1,1,1,1,-1,-1,0]], dtype=int)

def states(n):
    for bits in range(1 << (n-1)):
        x=np.ones(n,dtype=int)
        for i in range(n-1):
            if bits>>i&1:x[i]=-1
        yield x

def Q(A):
    return max(abs(int(x@A@x)) for x in states(len(A)))

def grounds(A):
    q=Q(A)
    out=[]
    for x in states(len(A)):
        e=int(x@A@x)
        for sig in [1,-1]:
            if sig*e==q:
                out.append((sig,x,sig*(x[:,None]*A*x[None,:])))
    return q,out

def cut_parity_nonzero(F, verts):
    # F set of sorted pairs, noncut iff an odd triangle exists
    for tri in combinations(verts,3):
        if sum(tuple(sorted(e)) in F for e in combinations(tri,2))%2:
            return True
    return False

def scan(A):
    n=len(A); q, gs=grounds(A)
    records=[]
    for gi,(sig,x,C) in enumerate(gs):
      assert int(np.ones(n,dtype=int)@C@np.ones(n,dtype=int))==q
      for mask in range(1,1<<(n-1)): # root n-1 in T, nonempty S
        S=[i for i in range(n) if mask>>i&1]; T=[i for i in range(n) if i not in S]
        D=C[np.ix_(S,S)]; E=C[np.ix_(T,T)]; B=C[np.ix_(S,T)]
        oneS=np.ones(len(S),dtype=int);oneT=np.ones(len(T),dtype=int)
        HD=int(oneS@D@oneS);HE=int(oneT@E@oneT);c=int(oneS@B@oneT)
        g=4*c
        if g<0: raise AssertionError((g,S))
        qD=Q(D) if S else 0;qE=Q(E) if T else 0
        gd=qD-HD;ge=qE-HE;z=qD+qE-q
        LF=int(np.abs(B@oneT).sum()); LB=int(np.abs(B.T@oneS).sum())
        ellF=2*LF-ge;ellB=2*LB-gd;f=max(ellF,0)
        h=max(g//2-z-f,0)
        alt=max(g-gd-max(ge,2*LF),0)
        assert h==alt
        Nrow=int(np.maximum(-(B@oneT),0).sum())
        # signed weight maximum over genuinely non-cut internal subsets
        ed=[(S[ii],S[jj]) for ii in range(len(S)) for jj in range(ii+1,len(S))]
        P={e for e in ed if C[e[0],e[1]]>0}
        if len(S)<3:
          best=-10**9;bestF=None
        elif cut_parity_nonzero(P,S):
          best=len(P);bestF=P
        else:
          # Any one-edge toggle leaves the cut space, at signed cost one.
          e=ed[0]
          best=len(P)-1;bestF=P.symmetric_difference({e})
        records.append(dict(gi=gi,S=S,T=T,h=h,g=g,c=c,gd=gd,ge=ge,z=z,
                            LF=LF,LB=LB,ellF=ellF,ellB=ellB,Nrow=Nrow,
                            best_res=best,bestF=bestF))
    return q,records

def residual_depths(A, gi, S):
    q,gs=grounds(A); C=gs[gi][2]; S=list(S)
    ed=[(S[ii],S[jj]) for ii in range(len(S)) for jj in range(ii+1,len(S))]
    layer=[]
    n=len(A)
    for x in states(n):
      base=int(x@C@x)
      for sig in [1,-1]:
        score=sig*base; delta=q-score
        z={(i,j):sig*x[i]*x[j] for i in range(n) for j in range(i+1,n)}
        layer.append((delta,z))
    ans=[]
    for em in range(1<<len(ed)):
      R={ed[k] for k in range(len(ed)) if em>>k&1}
      if not cut_parity_nonzero(R,S):continue
      cert=[]
      for delta,z in layer:
        Y=-4*sum(C[i,j]*z[(i,j)] for i,j in R)
        if Y>=delta: cert.append((delta,Y,z))
      if not cert: raise AssertionError('global certificate absent')
      md=min(x[0] for x in cert)
      minpay=min(x[1] for x in cert if x[0]==md)
      ans.append((md,minpay,R,len(cert)))
    ans.sort(key=lambda z:(z[0],z[1]),reverse=True)
    return C,ans

def exact_report(A, gi, S):
    C,ans=residual_depths(A,gi,S)
    print('EXACT gi,S',gi,S,'ground tuple sig,x=',grounds(A)[1][gi][:2])
    print(C)
    print('residual best forced depth/pay/R/#cert',ans[:10])

def audit_s3(A):
  q,gs=grounds(A); n=len(A); out=[]
  xlist=list(states(n))
  for gi,(_,_,C) in enumerate(gs):
    layers=[]
    for x in xlist:
      e=int(x@C@x)
      for sig in [1,-1]:
        layers.append((q-sig*e, sig*(x[:,None]*x[None,:])))
    for S in combinations(range(n-1),3):
      T=[i for i in range(n) if i not in S]
      D=C[np.ix_(S,S)];E=C[np.ix_(T,T)];B=C[np.ix_(S,T)]
      oS=np.ones(3,dtype=int);oT=np.ones(len(T),dtype=int)
      HD=int(oS@D@oS);HE=int(oT@E@oT);c=int(oS@B@oT);g=4*c
      gd=Q(D)-HD;ge=Q(E)-HE;z=Q(D)+Q(E)-q
      LF=int(np.abs(B@oT).sum());ellF=2*LF-ge
      h=max(g//2-z-max(ellF,0),0)
      if h<=0:continue
      ed=list(combinations(S,2)); stats=[]
      for em in range(1<<3):
        R={ed[k] for k in range(3) if em>>k&1}
        if not cut_parity_nonzero(R,S):continue
        cert=[]
        for delta,Z in layers:
          Y=-4*sum(C[i,j]*Z[i,j] for i,j in R)
          if Y>=delta: cert.append((delta,Y))
        stats.append((min(d for d,y in cert),min(y for d,y in cert),R))
      maxdepth=max(x[0] for x in stats)
      max_minpay=max(x[1] for x in stats)
      key=(len(S),h,g,c,gd,ge,z,LF,ellF)
      out.append((h,maxdepth,max_minpay,gi,S,stats,key))
  from collections import Counter
  print('S3 AUDIT count',len(out),'maxdepth dist',Counter(z[1] for z in out),
        'min max_minpay/h',min((z[2]/z[0],z[:5]) for z in out))
  for z in sorted(out,key=lambda z:(z[2]/z[0],z[0]))[:10]: print('S3 worst',z[:5],z[2]/z[0])
  groups={}
  for z in out:groups.setdefault(z[6],[]).append(z)
  varied=[(k,min(v[2] for v in vs),max(v[2] for v in vs),len(vs)) for k,vs in groups.items() if len({v[2] for v in vs})>1]
  print('S3 same scalar varied payment',varied[:20])
  return out

for name,A in [('n5',A5),('n9',A9)]:
    q,rr=scan(A)
    pos=[r for r in rr if r['h']>0]
    print(name,'q',q,'grounds',len(grounds(A)[1]),'nodes',len(rr),'hpos',len(pos),'maxh',max((r['h'] for r in pos),default=0))
    from collections import Counter
    print('sizes',Counter(len(r['S']) for r in pos),'hvals',Counter(r['h'] for r in pos))
    print('no internal residual',sum(r['bestF'] is None for r in pos))
    for r in sorted(pos,key=lambda z:(-z['h'],len(z['S'])))[:30]:
       print({k:r[k] for k in ['gi','S','h','g','c','gd','ge','z','LF','LB','ellF','ellB','Nrow','best_res','bestF']})
    print()

exact_report(A5,1,[3])
exact_report(A9,0,[0,1,4])
exact_report(A9,0,[0,1,4,6])
audit_s3(A5)
audit_s3(A9)
