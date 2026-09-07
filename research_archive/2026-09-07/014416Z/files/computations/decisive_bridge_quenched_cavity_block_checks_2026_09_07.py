"""Actual globally reoptimized quenched orientation objective; exploratory.

One-dimensional Gaussian quadrature is repeated at two orders. All signings are
enumerated up to switching. Reports the precise direct derivative as well as
the proposed first-plus-second cavity expansion on a block variance path.
"""
import argparse
import itertools
import json
import numpy as np
from scipy.special import logsumexp, roots_hermitenorm
from decisive_bridge_actual_cavity_psd_checks_2026_09_07 import objects

def logcosh(x):
    return np.logaddexp(x,-x)-np.log(2)

def check(n,beta,u,h,order):
    edges,chars,signs=objects(n)
    p=n//2
    r=p/n
    d=np.array([1/r-1 if i<p and j<p else 1/(1-r)-1 if i>=p and j>=p
                else -1 for i,j in edges])
    w=np.sqrt(1+u*d)
    lam=beta/np.sqrt(n)
    E=lam*(signs*w)@chars.T
    lp=logsumexp(E,axis=1)
    lm=logsumexp(-E,axis=1)
    gap=(lp-lm)/2
    g,gw=roots_hermitenorm(order)
    gw/=np.sqrt(2*np.pi)
    # Group equal branch free energies to avoid a huge quadrature array.
    pair=np.round(np.column_stack(((lp+lm)/2,np.abs(gap))),12)
    unique,inv=np.unique(pair,axis=0,return_inverse=True)
    cost=unique[:,0]+logcosh(unique[:,1,None]+h*g)@gw
    vals=cost[inv]
    ids=np.flatnonzero(vals<=vals.min()+2e-10)
    worst=-float('inf')
    best=float('inf')
    direct=[]
    first=[]
    second=[]
    for ix in ids:
        pp=np.exp(E[ix]-lp[ix]); pm=np.exp(-E[ix]-lm[ix])
        rp=pp@chars; rm=pm@chars
        orient=np.tanh(gap[ix]+h*g)
        mean=(rp-rm)/2+(rp+rm)/2*(orient@gw)
        direct.append(float(lam/2*np.sum(d/w*signs[ix]*mean)/n))
        M=[]; V=[]
        for e in range(len(edges)):
            cav=E[ix]-lam*w[e]*signs[ix,e]*chars[:,e]
            cp=logsumexp(cav); cm=logsumexp(-cav)
            xp=np.exp(cav-cp)@chars[:,e]
            xm=np.exp(-cav-cm)@chars[:,e]
            rr=(xp-xm)/2+(xp+xm)/2*np.tanh((cp-cm)/2+h*g)
            M.append(abs(float(rr@gw)))
            V.append(float(rr**2@gw))
        M=np.array(M);V=np.array(V)
        first.append(float(-lam/(2*n)*np.sum(d/w*M)))
        second.append(float(-lam/(2*n)*np.sum(d/w*np.tanh(lam*w)*V)))
        combined=first[-1]+second[-1]
        worst=max(worst,combined);best=min(best,combined)
    baseline=float(lam/(2*n)*np.sum(d/w*np.tanh(lam*w)))
    return dict(n=n,beta=beta,u=u,h=h,quadrature=order,minimizers=len(ids),
                objective=float(vals.min()),direct_derivative=[min(direct),max(direct)],
                first_range=[min(first),max(first)],second_range=[min(second),max(second)],
                combined_range=[best,worst],baseline=baseline)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--n',type=int,default=6)
    args=ap.parse_args()
    for beta in (1.,2.,4.):
        for h in (1.,3.):
            for u in (.25,.5,.75):
                for order in (128,256):
                    print(json.dumps(check(args.n,beta,u,h,order)),flush=True)
