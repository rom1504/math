"""Actual noisy-bent source family: numerical homogeneous transport diagnostic.

Physical row entropy is h(epsilon), not an arbitrary source entropy.
Gauss-Hermite discretization and floating convex optimization are diagnostics,
not a certificate for a growing-order upper or lower construction.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from numpy.polynomial.hermite import hermgauss


def source(alpha,order):
    if alpha==1: return np.array([1.]),np.array([1.])
    z,w=hermgauss(order); values=np.abs(alpha+math.sqrt(2*(1-alpha*alpha))*z)
    weights=w/math.sqrt(math.pi)
    if alpha==0:
        values=values[order//2:]; weights=2*weights[order//2:]
    return values,weights


def solve(values,weights,t,start=None):
    n=len(values)
    inds=np.array(list(itertools.product(range(n),repeat=4)),dtype=np.int16)
    a,b,c,d=[values[inds[:,q]] for q in range(4)]
    energies=[]
    for e,f,g,h in itertools.product((-1,1),repeat=4):
        energies.append(math.sqrt(2)*(a*c*e*g+a*d*e*h+b*c*f*g-b*d*f*h))
    logkernel=logsumexp(t*np.array(energies),axis=0)-math.log(16)
    base=np.log(weights[inds]).sum(axis=1)+logkernel

    def objective(z):
        phi=np.r_[z,0.]
        logs=base+phi[inds].sum(axis=1); normal=logsumexp(logs)
        probability=np.exp(logs-normal)
        gradient=-4*weights.copy()
        for column in range(4):
            gradient+=np.bincount(inds[:,column],weights=probability,minlength=n)
        return float(normal-4*weights@phi),gradient[:-1]

    if n==1: return objective(np.array([]))[0],0.,[],True
    guess=np.zeros(n-1) if start is None else np.array(start)
    fit=minimize(objective,guess,jac=True,method='L-BFGS-B',
                 options=dict(gtol=2e-9,ftol=1e-13,maxiter=600))
    value,gradient=objective(fit.x)
    return value,float(np.max(np.abs(gradient))),fit.x.tolist(),bool(fit.success)


def run():
    records=[]; times=[1.,2.,4.,8.,16.]
    for order in (12,20):
        for alpha in (0.,.2,.4,.6,.8,.95,1.):
            values,weights=source(alpha,order); epsilon=(1-alpha)/2
            entropy=0. if epsilon==0 else -epsilon*math.log(epsilon)-(1-epsilon)*math.log(1-epsilon)
            entries=[]; start=None
            for t in times:
                value,residual,start,success=solve(values,weights,t,start)
                entries.append(dict(t=t,transport=value,gradient_residual=residual,success=success,
                                    cap_expression=(entropy+value/4)/(2*t),potential=start))
            row=dict(order=order,alpha=alpha,epsilon=epsilon,physical_entropy=entropy,
                     support=values.tolist(),weights=weights.tolist(),values=entries,
                     best_expression=min(v['cap_expression'] for v in entries))
            records.append(row)
            print(json.dumps({k:row[k] for k in ('order','alpha','best_expression')}),flush=True)
    out=dict(status='FLOATING SOURCE/TRANSPORT DIAGNOSTIC ONLY',records=records,
             gaussian_exact_best=math.sqrt(15)/8)
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_noisy_bent_transport.json'
    target.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__': run()
