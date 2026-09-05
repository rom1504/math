#!/usr/bin/env python3
"""Exact finite spin averages for custom-response orientation differences.

All signings are covered modulo diagonal switching by fixing first-row+1;
global spin reversal halves the spin enumeration. Field thresholds use exact
integer comparisons. The second cutoff floor(a sqrt(n-1)) is 0 for n<=4
and 1 for5<=n<=12, with a=2phi(7/8); these wide intervals are not close to
an integer. A finite negative example would not alone be asymptotic evidence.
"""
import argparse
import json
import numpy as np


def test(n, batch):
    assert 2 <= n <= 12
    m = n-1
    edges = [(i,j) for i in range(1,n) for j in range(i+1,n)]
    count = 1 << len(edges)
    spin_codes = np.arange(1 << (n-1), dtype=np.int64)
    s = np.ones((len(spin_codes), n), dtype=np.int64)
    s[:, 1:] = 1-2*((spin_codes[:,None] >> np.arange(n-1)) & 1)
    cutoff = 0 if n <= 4 else 1
    best = None
    maxval = None
    witness = None
    for start in range(0, count, batch):
        codes = np.arange(start, min(count,start+batch), dtype=np.int64)
        a = np.ones((len(codes),n,n),dtype=np.int64)
        a[:,np.arange(n),np.arange(n)]=0
        for k,(i,j) in enumerate(edges):
            a[:,i,j]=a[:,j,i]=1-2*((codes>>k)&1)
        g = np.einsum('si,aij->asj',s,a,optimize=True)
        v = s[None,:,:]*(64*g*g <= 49*m)
        yr = np.einsum('asi,aij->asj',v,a,optimize=True)
        large = np.abs(yr)>cutoff
        yp = np.where(large,np.sign(yr),s[None,:,:])
        ym = np.where(large,-np.sign(yr),s[None,:,:])
        ep = np.einsum('asi,aij,asj->as',yp,a,yp,optimize=True)
        em = np.einsum('asi,aij,asj->as',ym,a,ym,optimize=True)
        totals = (ep-em).sum(axis=1)
        idx=int(totals.argmin())
        local=int(totals[idx])
        if best is None or local<best:
            best=local
            witness=a[idx].tolist()
        maxval=max(maxval or 0,int(totals.max()))
    return {"order":n,"switching_representatives":count,
            "projective_spins":len(s),"minimum_doubled_energy_difference_sum":best,
            "maximum_doubled_energy_difference_sum":maxval,
            "minimum_e_y_times_sqrt_n_minus_1":f"{best}/{4*n*len(s)}",
            "second_integer_cutoff":cutoff,"minimizing_matrix":witness,
            "classification":"exact integer finite enumeration"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--orders',type=int,nargs='+',default=(3,4,5,6,7))
    p.add_argument('--batch',type=int,default=128)
    args=p.parse_args()
    for n in args.orders:
        print(json.dumps(test(n,args.batch)),flush=True)


if __name__=='__main__':
    main()
