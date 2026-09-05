#!/usr/bin/env python3
"""Finite countertests of the exact smooth-mask optimizer, rho=.94 alpha=.81.

The s value lies in the separately certified interval. Floating point is
used only for Monte Carlo diagnostics, not the rigorous lower certificate.
"""
import argparse
import json
import math

import numpy as np
from scipy.special import ndtr

from fresh_custom_response_falsify import models


RHO = 47/50
ALPHA = 81/100
NOISE = math.sqrt(1-RHO*RHO)
PROB = math.erf(ALPHA/math.sqrt(2))
SD_K = 0.3637488664617
TARGET = 0.385785876909


def hstar(g):
    k = ndtr((ALPHA-RHO*g)/NOISE)-ndtr((-ALPHA-RHO*g)/NOISE)
    return RHO + NOISE*(k-PROB)/SD_K


def test(label,a,rng,samples,batch):
    n=len(a)
    bmat=a/math.sqrt(n-1)
    eigs=np.linalg.eigvalsh(bmat)
    acc=np.zeros(7)
    sq=np.zeros(7)
    for start in range(0,samples,batch):
        count=min(batch,samples-start)
        spins=rng.choice((-1.,1.),size=(count,n))
        g=spins@bmat
        h=hstar(g)
        y=(spins*h)@bmat
        yp=np.where(ALPHA*spins+y>=0,1.,-1.)
        ym=np.where(ALPHA*spins-y>=0,1.,-1.)
        ep=np.sum(yp*(yp@bmat),axis=1)/(2*n)
        em=np.sum(ym*(ym@bmat),axis=1)/(2*n)
        vals=np.stack(((ep-em)/2,(ep+em)/2,h.mean(1),(h*h).mean(1),
                       (y*y).mean(1),(y**4).mean(1),(y*g).mean(1)),axis=1)
        acc+=vals.sum(0)
        sq+=(vals*vals).sum(0)
    mu=acc/samples
    se=np.sqrt(np.maximum(0,sq/samples-mu*mu)/(samples-1))
    keys=('oriented_energy','unoriented_energy','h_mean','h_second',
          'rooted_second','rooted_fourth','rooted_cov_first')
    return {'model':label,'order':n,'samples':samples,'rho':RHO,'alpha':ALPHA,
            'predicted_oriented_energy':TARGET,
            'op_B':float(np.max(np.abs(eigs))),
            'spectral_third':float(np.mean(eigs**3)),
            'estimates':{k:{'mean':float(m),'sample_se':float(s)} for k,m,s in zip(keys,mu,se)},
            'classification':'finite Monte Carlo stress test, not proof'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--orders',type=int,nargs='+',default=(128,256,512))
    p.add_argument('--samples',type=int,default=8192)
    p.add_argument('--batch',type=int,default=256)
    p.add_argument('--seed',type=int,default=2026090524)
    p.add_argument('--quick',action='store_true')
    args=p.parse_args()
    rng=np.random.default_rng(args.seed)
    for n in args.orders:
        for label,a in models(n,rng,args.quick):
            print(json.dumps(test(label,a,rng,args.samples,args.batch)),flush=True)


if __name__=='__main__':
    main()
