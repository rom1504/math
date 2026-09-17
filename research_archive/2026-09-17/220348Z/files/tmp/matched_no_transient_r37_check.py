#!/usr/bin/env python3
"""Wave 37 finite audit of matched exclusion/context migration splitting."""

from __future__ import annotations

import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp, expit

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, A8, A9
from harmonic_resonance_r34_check import Endpoint


def kl(p,q):
    p=np.asarray(p);q=np.asarray(q); mask=p>0
    return float(np.sum(p[mask]*np.log(p[mask]/q[mask])))


def dynamic_selector_audit(H: Endpoint, t: float):
    def law(s):
        z=H.lognu+s*H.g; z-=float(logsumexp(z)); return np.exp(z)
    mt,m0,mp=law(t),law(0),law(t-1)
    max_omit=max_chain=max_ineq=0.0
    for bit,edges in enumerate(H.edges[1:],start=1):
        for x,y in edges:
            M=mt[x]+mt[y]
            bit_t=np.array([mt[x],mt[y]])/M
            bit_0=np.array([m0[x],m0[y]])/(m0[x]+m0[y])
            bit_prev=np.array([mp[x],mp[y]])/(mp[x]+mp[y])
            joint=np.stack([mt[x]*H.posterior[:,x],mt[y]*H.posterior[:,y]],axis=1)/M
            ps=np.sum(joint,axis=1)
            # Every omitted component has the previous-temperature bit law.
            for si,S in enumerate(H.selectors):
                if bit not in S and ps[si]>1e-300:
                    max_omit=max(max_omit,float(np.max(np.abs(joint[si]/ps[si]-bit_prev))))
            J=sum(ps[si]*kl(joint[si]/ps[si],bit_prev)
                  for si in range(len(ps)) if ps[si]>0)
            I=0.0
            for si in range(len(ps)):
                for b in range(2):
                    if joint[si,b]>0:
                        I += joint[si,b]*np.log(joint[si,b]/(ps[si]*bit_t[b]))
            qshift=kl(bit_t,bit_prev)
            kbase=kl(bit_t,bit_0)
            max_chain=max(max_chain,abs(J-qshift-I))
            max_ineq=max(max_ineq,kbase-t*qshift)
    assert max_omit<3e-8 and max_chain<3e-8 and max_ineq<3e-9
    return max_omit,max_chain,max_ineq


def split_state(H: Endpoint, t: float):
    logmu = H.lognu + t * H.g
    logmu -= float(logsumexp(logmu))
    mu = np.exp(logmu)
    out = []
    for bit, edges in enumerate(H.edges):
        x, y = edges[:, 0], edges[:, 1]
        M = mu[x] + mu[y]
        u0 = H.lognu[y] - H.lognu[x]
        ut = u0 + t * (H.g[y] - H.g[x])
        p = expit(ut)
        Aprime = (1-p)*H.g[x] + p*H.g[y]
        kval = p*(ut-u0) - (np.logaddexp(0,ut)-np.logaddexp(0,u0))
        def cov(a,b):
            return float(M @ ((a-M@a)*(b-M@b)))
        total = cov(Aprime,kval)
        if bit == 0:
            out.append((total, np.nan, np.nan, np.nan))
            continue
        # Exact endpoint omission aggregate b_i(e)=f(d) r_i(d), invariant on edge.
        bx = H.f[x] * H.excl[bit,x]
        by = H.f[y] * H.excl[bit,y]
        assert np.max(np.abs(bx-by)) < 2e-9
        logb = np.log(bx)
        hprime = Aprime-logb
        cb = cov(logb,kval)
        cr = cov(hprime,kval)
        assert abs(total-cb-cr) < 2e-9
        out.append((total,cb,cr,float(np.min(H.excl[bit]))))
    return np.asarray(out)


def audit(A,beta,m,order=96):
    H=Endpoint(A,beta,m)
    nodes,weights=leggauss(order)
    signed=np.zeros((H.n,3)); adverse=np.zeros((H.n,3))
    mins=np.full((H.n,3),np.inf); maxs=np.full((H.n,3),-np.inf)
    for node,w0 in zip(nodes,weights):
        t=(node+1)/2; w=w0/2
        z=split_state(H,float(t))[:,:3]
        signed += w*z
        adverse += w*np.maximum(-z,0)
        mins=np.minimum(mins,z); maxs=np.maximum(maxs,z)
    return H,signed,adverse,mins,maxs


def flux_audit(H: Endpoint,order=96):
    nodes,weights=leggauss(order); flux=np.zeros(H.n); actual=np.zeros(H.n)
    for node,w0 in zip(nodes,weights):
        t=float((node+1)/2); w=w0/2
        def law(s):
            z=H.lognu+s*H.g;z-=float(logsumexp(z));return np.exp(z)
        mt,mp=law(t),law(t-1); mean=float(mt@H.g)
        zsplit=split_state(H,t)
        actual += w*np.maximum(-zsplit[:,0],0)
        for bit,edges in enumerate(H.edges):
            for x,y in edges:
                massprime=mt[x]*(H.g[x]-mean)+mt[y]*(H.g[y]-mean)
                if massprime>=0: continue
                qt=np.array([mt[x],mt[y]])/(mt[x]+mt[y])
                qp=np.array([mp[x],mp[y]])/(mp[x]+mp[y])
                flux[bit] += w*t*(-massprime)*kl(qt,qp)
    assert np.max(actual-flux)<2e-8
    return actual,flux


def main():
    for name,A,beta,m in [
        ("A4",A4,.5,3),("A4",A4,2,3),
        ("A6",A6,.5,3),("A6",A6,2,3),
        ("A8",A8,.5,4),("A8",A8,2,5),("A8",A8,2,7),
        ("A9",A9,.5,3),("A9",A9,2,3),("A9",A9,2,7)]:
        H,s,a,lo,hi=audit(A,beta,m)
        if (name,beta,m) in (("A6",.5,3),("A9",2,7)):
            dyn=[dynamic_selector_audit(H,t) for t in (.2,.5,1.0)]
            af,fl=flux_audit(H)
        else:
            dyn="skipped"
            af=fl="skipped"
        # columns total, omitted-aggregate log b, exclusion-shape
        print(name,beta,m,
              "signed vertex",np.sum(s[1:],axis=0),
              "adverse vertex",np.sum(a[1:],axis=0),
              "min",np.min(lo[1:],axis=0),"max",np.max(hi[1:],axis=0),
              "orientation signed/adverse",s[0,0],a[0,0])
        print(" dynamic max errors",dyn)
        print(" actual/one-step flux",af,fl)
        assert np.max(np.abs(s[1:,0]-s[1:,1]-s[1:,2]))<2e-8
    print("PASS matched_no_transient_r37_check")


if __name__ == "__main__": main()
