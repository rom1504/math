#!/usr/bin/env python3
"""Finite two-anchor Gaussian fixed point; exploratory, not a certificate.

Anchors are G0=U1 and G1=U He2(G0)/sqrt(2).  The innovation is orthogonal
to both anchors.  This tests a genuinely larger universal tree response.
"""
import argparse
import json
import math

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import gammaln, ndtr


def evaluate(alpha, rho0, rho1, degree=200, details=False, fixed_a=None):
    total = rho0*rho0 + rho1*rho1
    if not (0.05 < alpha < 2 and 0 < rho0 < .9999 and .01 < total < .9999):
        return 1.
    ph = math.exp(-alpha*alpha/2)/math.sqrt(2*math.pi)
    p = 2*ndtr(alpha)-1
    hn = [1., alpha]
    for r in range(1, degree):
        hn.append(alpha/math.sqrt(r+1)*hn[-1]-math.sqrt(r/(r+1))*hn[-2])
    weights = np.zeros(degree+1)
    noise2 = 1-total
    for r in range(2, degree+1, 2):
        c2 = 4*ph*ph*hn[r-1]**2/r
        ell = np.arange(r+1)
        logpmf = (gammaln(r+1)-gammaln(ell+1)-gammaln(r-ell+1)
                  + ell*math.log(noise2)+(r-ell)*math.log(total))
        weights[:r+1] += c2*np.exp(logpmf)
    # Remove the He2(G0)/sqrt(2) coefficient as well as the constant,
    # which is absent already because the total degree starts at two.
    c_star = -math.sqrt(2)*alpha*ph*rho0*rho0
    weights[0] -= c_star*c_star
    if min(weights) < -1e-12 or weights[0] <= 0:
        return 1.
    ell = np.arange(degree+1)

    def ratio(a):
        terms = weights/(a+ell)**2
        return np.dot(ell, terms)/np.sum(terms)

    if fixed_a is not None:
        a = fixed_a
    elif np.dot(ell, weights)/np.sum(weights) <= .999:
        a = 1e12
    else:
        a = brentq(lambda t: ratio(t)-.999, 1e-10, 1e10)
    norm = np.sum(weights/(a+ell)**2)
    cov = np.sum(weights/(a+ell))/math.sqrt(norm)
    w = rho0*p + rho1*c_star + math.sqrt(noise2)*cov
    if not 0 < w*w < p:
        return 1.
    sigma = math.sqrt(p-w*w)
    value = (2*w*ph*(2*ndtr(w*alpha/sigma)-1)
             +2*math.sqrt(2*p/math.pi)*ndtr(-alpha*math.sqrt(p)/sigma))
    if details:
        return {'bound':value, 'alpha':alpha, 'rho0':rho0, 'rho1':rho1,
                'a':a, 'derivative_energy':ratio(a), 'mass':p,
                'excluded_star_coefficient':c_star,
                'covariance':w, 'degree':degree, 'certified':False}
    return -value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--degree', type=int, default=200)
    args = parser.parse_args()
    reports = []
    for initial in ((.73,.8,-.2),(.7,.75,-.35),(.75,.9,-.2),(.7,.6,-.1)):
        result = minimize(lambda x:evaluate(*x,args.degree), initial,
                          method='Nelder-Mead',
                          options={'maxiter':600,'xatol':1e-8,'fatol':1e-10})
        report = evaluate(*result.x,args.degree,details=True)
        reports.append(report)
        print(json.dumps(report), flush=True)
    print(json.dumps({'best':max(reports,key=lambda r:r['bound']),
                      'certified':False}), flush=True)


if __name__ == '__main__':
    main()
