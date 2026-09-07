"""Independent directed-interval ternary upper reconstruction.

Does not import the proposing verifier. Avoids its closed-form gain branches:
evaluate the exact rational optimizer at conservative A/B endpoints directly.
Uses v=(p/A-1/B)/(1-p) for gradient interval enclosure. Saves every accepted
rectangle and its rational upper bound in a gzip JSON evidence file.
"""
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import time
import mpmath as mp

mp.iv.dps=60
p,t,threshold=F(24,25),F(97,20),-F(5151,6250)
prefix=Path(__file__).resolve().parent/'results'/Path(__file__).stem

def rat(raw):
    s,m,e,_=raw
    return F((-1 if s else 1)*m)*(F(2)**e)

def ends(v):return tuple(rat(r) for r in v._mpi_)
def iv(q):
    q=F(q)
    return mp.iv.mpf(q.numerator)/q.denominator
def span(a,b):return mp.iv.mpf([iv(a).a,iv(b).b])
P,T=iv(p),iv(t)

def AB(l,z):
    a=mp.iv.exp(l*z*z/P)-1
    h=l*z/P
    b=(mp.iv.exp(h)-mp.iv.exp(-h))**2/2
    return a,b

def conservative_gain(A,B):
    a,b=ends(A)[0],ends(B)[1]
    if a==0 or b==0:
        v=F(1) if a==0 and b>0 else F(0)
    else:
        v=min(F(1),max(F(0),(p/a-1/b)/(1-p)))
    return P*mp.iv.ln(1+iv(v*b))-mp.iv.ln(1+iv(v*a))

def constant(l):return mp.iv.ln(l*(2*T-l)/(T*T))/4-l

def crude(box):
    l,u,a,b=map(iv,box)
    Amin=mp.iv.exp(l*a*a/P)-1
    h=u*b/P
    Bmax=(mp.iv.exp(h)-mp.iv.exp(-h))**2/2
    return ends(constant(l)+conservative_gain(Amin,Bmax))[1]

def mvt(box):
    l,u,a,b=box
    if a==0:return None
    L,Z=span(l,u),span(a,b)
    A,B=AB(L,Z)
    if ends(A)[0]<=0 or ends(B)[0]<=0:return None
    vr=(P/A-1/B)/(1-P)
    vlo,vhi=ends(vr)
    clamp=lambda q:min(F(1),max(F(0),q))
    V=span(clamp(vlo),clamp(vhi))
    h=2*L*Z/P
    S=(mp.iv.exp(h)-mp.iv.exp(-h))/2
    R=mp.iv.exp(L*Z*Z/P)
    # Different grouping: factor the common V out before differentiating.
    gradL=(1/L-1/(2*T-L))/4-1+V*(2*Z*S/(1+V*B)-Z*Z*R/(P*(1+V*A)))
    gradZ=2*V*L*(S/(1+V*B)-Z*R/(P*(1+V*A)))
    boundL=max(abs(x) for x in ends(gradL))
    boundZ=max(abs(x) for x in ends(gradZ))
    midL,midZ=iv((l+u)/2),iv((a+b)/2)
    ca,cb=AB(midL,midZ)
    center=ends(constant(midL)+conservative_gain(ca,cb))[1]
    return center+boundL*(u-l)/2+boundZ*(b-a)/2

pending=[(F(1,2),t,F(0),F(1))]
leaves=[]
visited=meanvalue=0
start=time.monotonic()
while pending:
    box=pending.pop()
    visited+=1
    upper=crude(box)
    method='monotone'
    if upper>=threshold:
        alternative=mvt(box)
        if alternative is not None and alternative<upper:
            upper=alternative
            method='meanvalue'
    if upper<threshold:
        leaves.append([*map(str,box),str(upper),method])
        meanvalue+=method=='meanvalue'
    else:
        l,u,a,b=box
        if (u-l)>=(t-F(1,2))*(b-a):
            mid=(l+u)/2
            pending.extend([(l,mid,a,b),(mid,u,a,b)])
        else:
            mid=(a+b)/2
            pending.extend([(l,u,a,mid),(l,u,mid,b)])
    if visited%10000==0:
        print(json.dumps(dict(visited=visited,leaves=len(leaves),pending=len(pending),seconds=time.monotonic()-start)),flush=True)

rho=(mp.iv.sqrt(1+16*T*T)-1)/(4*T)
gaussian=-T*(1-rho)+mp.iv.ln(1-rho*rho)/4
assert ends(gaussian)[1]<threshold
cap=(T+P*mp.iv.ln(2)+iv(threshold))/(2*T*mp.iv.sqrt(P))
assert ends(cap)[1]<F(493608094,10**9)
# Independent coverage bookkeeping: binary-tree node identity and exact area.
assert visited==2*len(leaves)-1
area=sum((F(u)-F(l))*(F(b)-F(a)) for l,u,a,b,_,_ in leaves)
assert area==t-F(1,2)
result=dict(status='independent directed replay PASS',p=str(p),t=str(t),target=str(threshold),
            visited=visited,leaves=len(leaves),meanvalue_leaves=meanvalue,
            exact_covered_area=str(area),worst_upper=str(max(F(r[4]) for r in leaves)),
            gaussian=list(map(str,ends(gaussian))),cap=list(map(str,ends(cap))),
            seconds=time.monotonic()-start)
prefix.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
with gzip.open(str(prefix)+'_leaves.json.gz','wt') as output:json.dump(leaves,output)
print(json.dumps(result,indent=2),flush=True)
