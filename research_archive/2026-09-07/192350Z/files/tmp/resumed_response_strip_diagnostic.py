"""Numerical-only full Gaussian-direction check for strip feedback."""
import sys, json, math
from pathlib import Path
import numpy as np
from collections import Counter
from collections import defaultdict
from itertools import product
from scipy.special import gammaln
from scipy.special import ndtr
from scipy.integrate import quad_vec
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'computations'))
from fresh_finite_anchor_fixed_point_diagnostic import evaluate

P0=1/math.sqrt(2*math.pi)
def phi(x): return P0*np.exp(-np.asarray(x)**2/2)

raw=json.loads((Path(__file__).resolve().parents[1]/'computations/results/fresh_finite_anchor_fixed_point_certificate.json').read_text())
from fractions import Fraction
base=evaluate(float(Fraction(raw['alpha'])),np.array([float(Fraction(r)) for r in raw['rho']]),raw['anchor_children'],200,True,fixed_a=float(Fraction(raw['resolvent_a'])))
alpha,p0,c0=base['alpha'],base['mass'],base['covariance']
s0=math.sqrt(1-c0*c0/p0)
nodes,weights=leggauss(128)

def inverse_decomposition():
    rho=np.array(base['rho']); a=base['a']; n2=1-rho@rho
    hn=[1.,alpha]
    for r in range(1,200):hn.append(alpha/math.sqrt(r+1)*hn[-1]-math.sqrt(r/(r+1))*hn[-2])
    beta=np.zeros(201);beta[0]=p0
    for r in range(2,201,2): beta[r]=-2*phi(alpha)*hn[r-1]/math.sqrt(r)
    ac=[]
    for child in base['anchor_children']:
        mul=Counter(child);r=len(child)
        cc=beta[r]*math.sqrt(math.factorial(r)/math.prod(math.factorial(v) for v in mul.values()))
        for j,m in mul.items():cc*=rho[j]**m
        ac.append(cc)
    ac=np.array(ac)
    weight=np.zeros(201);weight[0]=p0*p0
    for r in range(2,201,2):
        l=np.arange(r+1)
        weight[:r+1]+=beta[r]**2*np.exp(gammaln(r+1)-gammaln(l+1)-gammaln(r-l+1)+l*math.log(n2)+(r-l)*math.log(1-n2))
    weight[0]-=ac@ac
    norm=np.sum(weight/(a+np.arange(201))**2)
    scale=math.sqrt(n2/norm)
    delta=rho-scale/a*ac
    return rho,beta,ac,norm,scale,delta

def correction_fourth(delta):
    indices=[]
    for child in base['anchor_children']:
        ct=Counter(child); indices.append(tuple(ct[i] for i in range(9)))
    terms=defaultdict(float)
    for a in range(1,len(delta)):
        for b in range(1,len(delta)):
            pairs=[]
            for m,n in zip(indices[a],indices[b]):
                choices=[]
                for k in range(min(m,n)+1):
                    j=m+n-2*k
                    cc=math.factorial(k)*math.comb(m,k)*math.comb(n,k)*math.sqrt(math.factorial(j)/(math.factorial(m)*math.factorial(n)))
                    choices.append((j,cc))
                pairs.append(choices)
            for items in product(*pairs):
                idx=tuple(x[0] for x in items)
                terms[idx]+=delta[a]*delta[b]*math.prod(x[1] for x in items)
    return sum(v*v for v in terms.values()),len(terms)

def moments(variance,covariance,inner,outer):
    s=math.sqrt(1-covariance*covariance/variance)
    vals=np.zeros(7)
    for lo,hi in ((0,inner),(inner,outer)):
        if lo==hi: continue
        w=lo+(nodes+1)*(hi-lo)/2
        x=w/math.sqrt(variance)
        m=covariance/variance*w
        l=(-alpha-m)/s
        u=(alpha-m)/s
        t0=ndtr(l)+ndtr(-u)
        t1=phi(u)-phi(l)
        t2=t0+u*phi(u)-l*phi(l)
        q=np.ones_like(w) if hi==inner else (outer-w)/(outer-inner)
        dp=0 if hi==inner else 1/(outer-inner)**2
        dens=phi(x)/math.sqrt(variance)
        blocks=np.array([t0*q*q,t0*w*q,t0*dp,t0*x*x*dp,t1*x*dp,t2*dp,t0*q])
        vals+=2*(hi-lo)/2*(blocks*dens)@weights
    lam=max(vals[2],np.linalg.eigvalsh([[vals[3],vals[4]],[vals[4],vals[5]]])[-1])
    rem=math.sqrt(variance-covariance*covariance)
    tail=(2*covariance*phi(alpha)*(2*ndtr(covariance*alpha/rem)-1)
          +4*math.sqrt(variance)*P0*ndtr(-alpha*math.sqrt(variance)/rem))
    return vals[0],lam,tail-vals[1],vals[-1]

if __name__=='__main__':
    print(json.dumps(base,indent=2))
    print('s0',s0)
    rho,beta,ac,norm,scale,delta=inverse_decomposition()
    print('inverse correction scale',scale,'norm',np.linalg.norm(delta),'coefficients',delta)
    for outer in [.01,.02,.03,.04,.06,.08,.1,.12,.16,.2,.3]:
        inner=outer/100
        start=moments(p0,c0,inner,outer)
        d=math.sqrt(start[0])
        for k in range(100):
            worst=moments(p0+d*d,c0-s0*d,inner,outer)
            dnew=math.sqrt(worst[0])
            if abs(d-dnew)<1e-13: break
            d=dnew
        lam=0; low=1; mass=0
        for t in np.linspace(0,d,21):
            for c in np.linspace(c0-s0*t,c0+s0*t,41):
                ms=moments(p0+t*t,c,inner,outer)
                lam=max(lam,ms[1]); low=min(low,ms[2]); mass=max(mass,ms[0])
        print({'outer':outer,'d':d,'initial':start,'worst_lambda':lam,'worst_J':low,'invariant_ratio':mass/d**2},flush=True)
