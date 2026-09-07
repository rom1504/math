"""Bounded diagnostic for broader Rademacher-convolution reference laws.

Uses exact Gaussian integration and floating finite discrete expectations.
This is not a classification of physical Walsh spectra or a rigorous
maximization over all mixtures. Explicit two-Rademacher counterexamples
also have elementary rational lower proofs in the companion artifact.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


def source(coefficients):
    if len(set(coefficients))==1:
        length=len(coefficients)
        return ((2*np.arange(length+1)-length)/math.sqrt(length),
                np.array([math.comb(length,j)/2**length for j in range(length+1)]))
    a=np.array(coefficients,dtype=float); a/=np.linalg.norm(a)
    points=np.array([np.dot(a,s) for s in itertools.product((-1,1),repeat=len(a))])
    values,counts=np.unique(np.round(points,13),return_counts=True)
    return values,counts/len(points)


def run():
    coefficients=[(1,),(1,1),(1,2),(1,3),(1,1,1),(1,1,2),(1,2,3),
                  (1,2,4),(1,1,1,1),(1,1,1,2),(1,2,3,4),
                  (1,)*8,(1,)*16,(1,)*32]
    bases=[source(c) for c in coefficients]
    kinds=[(i,a) for i in range(len(bases)) for a in (.5,.8,1.)]
    kinds.insert(0,(0,0.))
    t=2*math.sqrt(15); lam=7.5
    prepared={}
    for i in range(len(bases)):
        for j in range(i,len(bases)):
            vi,wi=bases[i]; vj,wj=bases[j]
            inds=np.array(list(itertools.product(range(len(vi)),range(len(vi)),
                                                range(len(vj)),range(len(vj)))),dtype=np.int16)
            a,b=vi[inds[:,0]],vi[inds[:,1]]; c,d=vj[inds[:,2]],vj[inds[:,3]]
            squares=a*a+b*b+c*c+d*d
            cross=math.sqrt(2)*(a*c+a*d+b*c-b*d)
            weights=np.log(wi[inds[:,:2]]).sum(axis=1)+np.log(wj[inds[:,2:]]).sum(axis=1)
            prepared[i,j]=(squares,cross,weights)
    size=len(kinds); matrix=np.zeros((size,size))
    for i,(bi,alpha) in enumerate(kinds):
        for j in range(i,size):
            bj,beta=kinds[j]; p=alpha*beta; den=16-15*p*p
            sq,cross,weight=prepared[min(bi,bj),max(bi,bj)]
            value=-math.log(den)+logsumexp(weight-lam*p*p/den*sq+t*p/den*cross)
            matrix[i,j]=matrix[j,i]=value
    rng=np.random.default_rng(9070905); mixtures=[]
    for start in [np.eye(size)[i] for i in range(size)]+[rng.dirichlet(np.ones(size)) for _ in range(20)]:
        def objective(w): return -float(w@matrix@w),-2*(matrix@w)
        fit=minimize(objective,start,jac=True,method='SLSQP',bounds=[(0,1)]*size,
                     constraints=[dict(type='eq',fun=lambda w:w.sum()-1,jac=lambda w:np.ones(size))],
                     options=dict(maxiter=1000,ftol=1e-12))
        val=float(fit.x@matrix@fit.x)
        mixtures.append(dict(log_kernel_average=val,weights=fit.x.tolist(),success=bool(fit.success),
                             full_entropy_cap_expression=(math.log(2)+lam+val/4)/(2*t)))
    best=max(mixtures,key=lambda r:r['log_kernel_average'])
    compatible=[]
    for label,allowed in [('even_Walsh',[0,8,12]),('odd_Walsh',[1,11,13])]:
        selected=[i for i,(base,alpha) in enumerate(kinds) if base in allowed or alpha==0]
        small=matrix[np.ix_(selected,selected)]; dimension=len(selected); choices=[]
        for start in [np.eye(dimension)[i] for i in range(dimension)]+[rng.dirichlet(np.ones(dimension)) for _ in range(20)]:
            fit=minimize(lambda w:(-float(w@small@w),-2*small@w),start,jac=True,
                         method='SLSQP',bounds=[(0,1)]*dimension,
                         constraints=[dict(type='eq',fun=lambda w:w.sum()-1,jac=lambda w:np.ones(dimension))],
                         options=dict(maxiter=1000,ftol=1e-12))
            val=float(fit.x@small@fit.x)
            choices.append(dict(log_kernel_average=val,weights=fit.x.tolist(),
                                full_entropy_cap_expression=(math.log(2)+lam+val/4)/(2*t)))
        winner=max(choices,key=lambda r:r['log_kernel_average'])
        compatible.append(dict(label=label,selected_indices=selected,best=winner))
    flat=kinds.index((0,1.)); semi=kinds.index((1,1.))
    out=dict(status='BOUNDED FLOATING KERNEL AND MIXTURE DIAGNOSTIC',
             coefficients=[list(c) for c in coefficients],
             kinds=[dict(base=i,alpha=a) for i,a in kinds],log_kernel=matrix.tolist(),
             flat_semibent_kernel=float(math.exp(matrix[flat,semi])),
             semibent_self_kernel=float(math.exp(matrix[semi,semi])),mixtures=mixtures,best=best,
             parity_compatible_mixtures=compatible)
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_rademacher_convolution_kernel_extended.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(flat_semibent=out['flat_semibent_kernel'],semibent_self=out['semibent_self_kernel'],
                         best_expression=best['full_entropy_cap_expression'],
                         support=[(kinds[i],float(w)) for i,w in enumerate(best['weights']) if w>1e-5],
                         compatible=[(r['label'],r['best']['full_entropy_cap_expression']) for r in compatible])),flush=True)


if __name__=='__main__': run()
