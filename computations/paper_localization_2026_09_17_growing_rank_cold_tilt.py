"""Exact finite-cube and Gaussian quadrature checks for cold feature tilts."""
import itertools
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tmp/paper_portfolio_2026_09_17/localization/growing_rank_cold_tilt.json"
RNG = np.random.default_rng(170921)


def cube(n):
    z = np.arange(1 << n, dtype=np.uint64)
    return (1 - 2*((z[:,None] >> np.arange(n,dtype=np.uint64)) & 1).astype(np.int8)).astype(float)


def hadamard(n):
    h = np.ones((1,1))
    while len(h)<n:
        h = np.block([[h,h],[h,-h]])
    return h


def quadrature(r, order=50):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    indices = np.array(list(itertools.product(range(order), repeat=r)))
    return math.sqrt(2)*nodes[indices], np.prod(weights[indices]/math.sqrt(math.pi),axis=1)


rows=[]
fourier_checks=0
for n,r in ((8,1),(8,2),(16,1),(16,2),(16,3),(16,4)):
    h=cube(n)
    u=hadamard(n)[:,:r]/math.sqrt(n)
    pp=u@u.T
    for v in (.3,.6,.85):
        beta=1/v-1
        rr=np.sum((h@u)**2,axis=1)
        raw=np.exp(-beta*rr/2)
        zz=raw.mean()
        law=raw/zz
        cov=h.T@(law[:,None]*h)/len(h)
        ref=np.eye(n)+(v-1)*(pp-np.diag(np.diag(pp)))
        cov_error=float(np.linalg.norm(cov-ref))
        max_response_error=0.
        for j in RNG.integers(0,len(h),size=16):
            a=u.T@h[j]/math.sqrt(n)
            target=math.sqrt(2/math.pi)*math.sqrt(1+(v-1)*float(a@a))
            response=float(np.mean(law*np.abs(h@h[j]))/math.sqrt(n))
            max_response_error=max(max_response_error,abs(response-target))
        rows.append(dict(n=n,rank=r,v=v,partition_relative_error=abs(zz/v**(r/2)-1),
                         covariance_frobenius_error=cov_error,
                         maximum_test_response_error=max_response_error))
        if r<=2:
            g,weights=quadrature(r)
            z=math.sqrt(beta)*(g@u.T)
            cc=np.cos(z)
            qz=float(weights@np.prod(cc,axis=1))
            assert abs(qz-zz)<2e-10,(n,r,v,qz,zz)
            fourier_checks+=1
            for i,j in ((0,1),(0,n-1),(n//2,n-1)):
                others=[k for k in range(n) if k not in (i,j)]
                # No tan/cos ratios: safe even at roots of individual cosines.
                numerator=-float(weights@(np.sin(z[:,i])*np.sin(z[:,j])*np.prod(cc[:,others],axis=1)))
                assert abs(numerator/zz-cov[i,j])<2e-10
                fourier_checks+=1
            for j in (0,len(h)//3,len(h)-1):
                for t in (.3,1.,2.):
                    cf=float(np.mean(law*np.cos(t*(h@h[j])/math.sqrt(n))))
                    qcf=float(weights@np.prod(np.cos(z+t*h[j]/math.sqrt(n)),axis=1))/zz
                    assert abs(cf-qcf)<3e-10,(n,r,v,t,cf,qcf)
                    fourier_checks+=1

assert all(math.isfinite(row["covariance_frobenius_error"]) for row in rows)
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(dict(status="PASS",fourier_identity_checks=fourier_checks,cases=rows),indent=2)+"\n")
print(json.dumps(dict(status="PASS",cases=len(rows),fourier_identity_checks=fourier_checks)))
