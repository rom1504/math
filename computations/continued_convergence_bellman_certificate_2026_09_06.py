"""Rational Bellman lower policies with independently checked interval costs.

--from-policy creates a certificate by rounding an optimizer diagnostic.
Without it, verification uses only integer arithmetic and mpmath intervals;
it does not optimize and does not import scipy or numpy.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from mpmath import iv


def classes(n):
    return [(i,j,s) for i in range(n) for j in range(i,n)
            for s in ([1] if i==0 else [1,-1])]


def children(q,n):
    rows=[0]*n;plus=[0]*(2*n-1);minus=plus.copy()
    for mass,(i,j,s) in zip(q,classes(n)):
        rows[i]+=mass;rows[j]+=mass
        plus[abs(i+s*j)]+=mass;minus[abs(i-s*j)]+=mass
    return rows,plus,minus


def create(source,denominator):
    import importlib.util
    import numpy as np
    from scipy.special import logsumexp
    raw=json.loads(source.read_text())
    record=max(raw['records'],key=lambda r:r['primal_leaf_value'])
    p=Fraction(str(raw['p']));depth=raw['depth'];t=raw['t'];D=denominator
    assert D%p.denominator==0
    profiles=[[D-int(D*p),int(D*p)]];nodes=[];offset=0
    for d in range(depth):
        output=[]
        for profile in profiles:
            n=len(profile);cs=classes(n);target=np.array(record['policy'][offset:offset+len(cs)])
            offset+=len(cs);row=np.zeros(n)
            for mass,(i,j,s) in zip(target,cs):row[i]+=mass/2;row[j]+=mass/2
            scale=min([1.]+[profile[i]/(D*row[i]) for i in range(n) if row[i]>0])
            masses={};offsum=[0]*n
            for i in range(n):
                for j in range(i+1,n):
                    indices=[c for c,z in enumerate(cs) if z[0]==i and z[1]==j]
                    total=2*int(max(0,float(target[indices].sum()*scale*D/2)))
                    masses[i,j]=total;offsum[i]+=total//2;offsum[j]+=total//2
            for i in range(n):masses[i,i]=profile[i]-offsum[i]
            q=[]
            for i in range(n):
                for j in range(i,n):
                    total=masses[i,j]
                    if i==0:q.append(total);continue
                    inds=[c for c,z in enumerate(cs) if z[0]==i and z[1]==j]
                    ratio=float(target[inds[0]]/target[inds].sum()) if target[inds].sum()>0 else .5
                    first=max(0,min(total,round(total*ratio)))
                    q.extend([first,total-first])
            assert min(q)>=0
            rows,plus,minus=children(q,n)
            assert rows==[2*z for z in profile]
            nodes.append(q);output.extend([plus,minus])
        profiles=output
    J=D*2**20;couplings=[]
    for profile in profiles:
        n=len(profile);probs=np.array(profile)/D;a=np.arange(n)/np.sqrt(float(p)*2**depth)
        logk=logsumexp(np.stack([-t*(a[:,None]-a[None,:])**2,
                                -t*(a[:,None]+a[None,:])**2]),axis=0)-np.log(2)
        keep=probs>0;ids=np.flatnonzero(keep);lp=np.log(probs[keep]);lk=logk[np.ix_(keep,keep)]
        u=np.zeros(len(ids))
        for it in range(2000):
            nxt=.5*(u-logsumexp(lk+lp[None,:]+u[None,:],axis=1))
            if max(abs(nxt-u))<1e-14:u=nxt;break
            u=nxt
        approx=np.exp(lp[:,None]+lp[None,:]+u[:,None]+u[None,:]+lk)
        off=approx.copy();np.fill_diagonal(off,0)
        scale=min([1.]+[probs[i]/off[k].sum() for k,i in enumerate(ids) if off[k].sum()>0])
        q=[[0]*n for _ in range(n)]
        for k,i in enumerate(ids):
            for l,j in enumerate(ids):
                if i<j:q[i][j]=q[j][i]=int(max(0,scale*approx[k,l]*J))
        for i in range(n):q[i][i]=profile[i]*(J//D)-sum(q[i])
        assert all(q[i][i]>=0 for i in range(n))
        couplings.append(q)
    return {'kind':'rational asymmetric Bellman LOWER certificate','p':[p.numerator,p.denominator],
      't':str(t),'depth':depth,'denominator':D,'node_class_masses':nodes,
      'coupling_denominator':J,'leaf_self_couplings':couplings,'source':str(source)}


def verify(cert,precision):
    iv.dps=precision
    D=cert['denominator'];J=cert['coupling_denominator'];depth=cert['depth']
    pnum,pden=cert['p'];p=iv.mpf(pnum)/pden;t=iv.mpf(cert['t']);ln2=iv.log(2)
    assert D%pden==0 and J%D==0
    def ent(masses,denom,multiplicities):
        result=iv.mpf(0)
        for z,mult in zip(masses,multiplicities):
            assert isinstance(z,int) and z>=0
            if z:result-=iv.mpf(z)/denom*iv.log(iv.mpf(z)/(denom*mult))
        return result
    value=iv.mpf(0);profiles=[[D-D*pnum//pden,D*pnum//pden]];index=0
    for d in range(depth):
        output=[]
        for profile in profiles:
            q=cert['node_class_masses'][index];index+=1;n=len(profile);cs=classes(n)
            assert len(q)==len(cs) and sum(q)==D
            rows,plus,minus=children(q,n)
            assert rows==[2*z for z in profile]
            mult=[1 if j==0 else (2 if i==j else 4) for i,j,s in cs]
            value+=(ent(q,D,mult)/2-ent(profile,D,[1]+[2]*(n-1)))/(2**d)
            output.extend([plus,minus])
        profiles=output
    assert index==len(cert['node_class_masses'])
    assert len(profiles)==len(cert['leaf_self_couplings'])
    for profile,q in zip(profiles,cert['leaf_self_couplings']):
        n=len(profile);assert len(q)==n
        assert all(len(row)==n for row in q)
        for i in range(n):
            assert sum(q[i])==profile[i]*(J//D)
            for j in range(n):assert isinstance(q[i][j],int) and q[i][j]>=0 and q[i][j]==q[j][i]
        cost=iv.mpf(0)
        for i in range(n):
            for j in range(n):
                if not q[i][j]:continue
                assert profile[i]>0 and profile[j]>0
                mass=iv.mpf(q[i][j])/J
                logk=-t*(i-j)**2/(p*2**depth)+iv.log(1+iv.exp(-4*t*i*j/(p*2**depth)))-ln2
                cost+=mass*(iv.log(iv.mpf(q[i][j])*D*D/(J*profile[i]*profile[j]))-logk)
        value-=cost/(2**(depth+1))
    exponent=p*ln2+value+t*(1-iv.sqrt(p))
    return {'Bellman_lower_interval':str(value),'certificate_exponent_lower_interval':str(exponent),
      'strictly_positive':bool(exponent.a>0),'precision':precision,'depth':depth,
      'all_integer_constraints':'PASS','meaning':'positive excludes this depth, p, t from the proposed strict upper certificate'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path)
    ap.add_argument('--from-policy',type=Path);ap.add_argument('--denominator',type=int,default=2**30)
    ap.add_argument('--precision',type=int,default=60);args=ap.parse_args()
    if args.from_policy:
        cert=create(args.from_policy,args.denominator)
        cert['verification']=verify(cert,args.precision)
        args.certificate.write_text(json.dumps(cert,indent=2)+'\n')
    else:cert=json.loads(args.certificate.read_text())
    print(json.dumps(verify(cert,args.precision),indent=2))


if __name__=='__main__':main()
