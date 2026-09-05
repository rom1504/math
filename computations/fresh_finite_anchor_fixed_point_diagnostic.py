#!/usr/bin/env python3
"""Finite ancestor anchors, exploratory Gaussian calculation only.

No matrix simulation and no claim of numerical certification.  Every anchor
is indexed by an even multiset of earlier children, exactly as in U.
"""
import argparse
from collections import Counter
import json
import math
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import gammaln, ndtr


def tree_family(max_size):
    children, sizes = [()], [1]
    for target in range(3, max_size+1, 2):
        candidates = []
        def visit(start, remaining, current):
            if remaining == 0:
                if len(current) % 2 == 0:
                    candidates.append(tuple(current))
                return
            for j in range(start, len(sizes)):
                if sizes[j] <= remaining:
                    visit(j, remaining-sizes[j], current+[j])
        visit(0, target-1, [])
        for item in candidates:
            children.append(item)
            sizes.append(target)
    return children, sizes


def evaluate(alpha, rho, children, degree=200, details=False, fixed_a=None):
    total = float(np.dot(rho,rho))
    if not (.15 < alpha < 1.5 and .01 < total < .99999):
        return 10.+total
    ph = math.exp(-alpha*alpha/2)/math.sqrt(2*math.pi)
    p = 2*ndtr(alpha)-1
    hn = [1., alpha]
    for r in range(1,degree):
        hn.append(alpha/math.sqrt(r+1)*hn[-1]-math.sqrt(r/(r+1))*hn[-2])
    c = np.zeros(degree+1)
    c[0] = p
    for r in range(2,degree+1,2):
        c[r] = -2*ph*hn[r-1]/math.sqrt(r)
    noise2 = 1-total
    weights = np.zeros(degree+1)
    weights[0] = p*p
    for r in range(2,degree+1,2):
        ell = np.arange(r+1)
        logpmf = (gammaln(r+1)-gammaln(ell+1)-gammaln(r-ell+1)
                  +ell*math.log(noise2)+(r-ell)*math.log(total))
        weights[:r+1] += c[r]**2*np.exp(logpmf)
    anchor_coeff = []
    for child in children:
        d = len(child)
        mult = Counter(child)
        coeff = c[d]*math.sqrt(math.factorial(d)/math.prod(
            math.factorial(v) for v in mult.values()))
        for j,m in mult.items():
            coeff *= rho[j]**m
        anchor_coeff.append(coeff)
    weights[0] -= np.dot(anchor_coeff,anchor_coeff)
    if weights[0] <= 1e-16:
        return 10.
    ell = np.arange(degree+1)
    def ratio(a):
        terms = weights/(a+ell)**2
        return np.dot(ell,terms)/np.sum(terms)
    if fixed_a is not None:
        a = fixed_a
    elif np.dot(ell,weights)/np.sum(weights) < .999:
        a = 1e12
    else:
        a = brentq(lambda z:ratio(z)-.999,1e-10,1e10)
    norm = np.sum(weights/(a+ell)**2)
    covariance = (np.dot(rho,anchor_coeff)+math.sqrt(noise2)
                  *np.sum(weights/(a+ell))/math.sqrt(norm))
    if not 0 < covariance < math.sqrt(p):
        return 10.
    sigma = math.sqrt(p-covariance*covariance)
    value = (2*covariance*ph*(2*ndtr(covariance*alpha/sigma)-1)
             +2*math.sqrt(2*p/math.pi)*ndtr(-alpha*math.sqrt(p)/sigma))
    if details:
        return {'bound':value,'alpha':alpha,'rho':list(rho),'a':a,
                'conditional_derivative_energy':ratio(a),'mass':p,
                'covariance':covariance,'degree':degree,
                'anchor_children':[list(t) for t in children],
                'certified':False}
    return -value


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-size',type=int,default=7)
    parser.add_argument('--degree',type=int,default=200)
    args=parser.parse_args()
    children,sizes=tree_family(args.max_size)
    rho=np.zeros(len(children))
    rho[:2]=[.8201,-.3210]
    # Learn the extra anchor directions from the threshold's own Hermite
    # coefficients; this initializes, it does not impose a fixed point.
    initial=np.concatenate(([.727],rho))
    result=minimize(lambda z:evaluate(z[0],z[1:],children,args.degree),
                    initial,method='BFGS',
                    options={'maxiter':250,'gtol':2e-7})
    report=evaluate(result.x[0],result.x[1:],children,args.degree,True)
    report['optimizer_success']=bool(result.success)
    report['optimizer_message']=str(result.message)
    report['anchor_sizes']=sizes
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
