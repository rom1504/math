"""Exhaustive actual global absolute-partition minimizers, n<=6 by default.

Switching fixes the first row of signs, without changing any partition value.
Checks ALL minimizing switching classes and all two-block tests for homogeneous
weights; for inhomogeneous weights checks the interpolation's fixed partition.
Numerical exploration only: minimizer ties use a declared tolerance.
"""
import argparse
import itertools
import json
import numpy as np

def objects(n):
    edges=list(itertools.combinations(range(n),2))
    spins=np.array([(1,)+x for x in itertools.product((-1.,1.),repeat=n-1)])
    chars=np.array([spins[:,i]*spins[:,j] for i,j in edges]).T
    free=[q for q,(i,j) in enumerate(edges) if i!=0]
    signs=np.ones((2**len(free),len(edges)))
    signs[:,free]=np.array(list(itertools.product((-1.,1.),repeat=len(free))))
    return edges,chars,signs

def evaluate(n,beta,u,parts,objects_n):
    edges,chars,signs=objects_n
    r=parts/n
    w=np.array([np.sqrt(1+u*(1/r-1)) if i<parts and j<parts
                else np.sqrt(1+u*(1/(1-r)-1)) if i>=parts and j>=parts
                else np.sqrt(1-u) for i,j in edges])
    lam=beta/np.sqrt(n)
    energies=lam*(signs*w)@chars.T
    absenergies=np.abs(energies)
    shift=absenergies.max(axis=1)
    ec=np.exp(absenergies-shift[:,None])
    ep=np.exp(-absenergies-shift[:,None])
    cosine=(ec+ep)/2
    sine=np.sign(energies)*(ec-ep)/2
    logz=shift+np.log(cosine.sum(axis=1))
    best=logz.min()
    ids=np.flatnonzero(logz<=best+2e-11)
    proj=np.eye(n)-np.ones((n,n))/n
    worst_eig=0.
    worst_block=-float('inf')
    worst_sign=None
    for ix in ids:
        full=(sine[ix]@chars)/cosine[ix].sum()
        z=signs[ix]*np.tanh(lam*w)
        cavity=(full-z)/(1-z*full)
        assert np.max(signs[ix]*cavity)<=2e-8
        K=np.zeros((n,n))
        for e,(i,j) in enumerate(edges):
            K[i,j]=K[j,i]=abs(cavity[e])
        eig=float(np.linalg.eigvalsh(proj@K@proj)[0])
        worst_eig=min(worst_eig,eig/np.sqrt(n))
        blocks=itertools.chain.from_iterable(itertools.combinations(range(n),a)
                      for a in range(1,n//2+1)) if u==0 else [range(parts)]
        for block in blocks:
            a=len(block)
            q=np.full(n,-np.sqrt(a/(n-a)))
            q[list(block)]=np.sqrt((n-a)/a)
            defect=float(-.5*q@K@q/n**1.5)
            if defect>worst_block:
                worst_block=defect
                worst_sign=int(ix)
    return dict(n=n,beta=beta,u=u,parts=parts,minimizers=len(ids),
                conditional_eig_over_sqrtn=worst_eig,
                worst_block_defect_over_n32=worst_block,witness=worst_sign)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-n',type=int,default=6)
    args=parser.parse_args()
    for n in range(3,args.max_n+1):
        ob=objects(n)
        for beta in (.2,1.,2.,4.,8.):
            for u in (0.,.25,.5,.75,1.):
                print(json.dumps(evaluate(n,beta,u,n//2,ob)),flush=True)
