#!/usr/bin/env python3
"""Independent floating-point anchored fixed-point optimization.

All Hermite sums are finite. This is a diagnostic, not an interval proof.
"""
import argparse
import json
import math

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import gammaln, ndtr


def objective(alpha, rho, degree=200, details=False, fixed_a=None):
    if not (0.05 < alpha < 2.5 and 0.01 < rho < .9999):
        return 1.
    ph = math.exp(-alpha*alpha/2)/math.sqrt(2*math.pi)
    p = 2*ndtr(alpha)-1
    hn = [1., alpha]
    for r in range(1, degree):
        hn.append(alpha/math.sqrt(r+1)*hn[-1]-math.sqrt(r/(r+1))*hn[-2])
    weights = np.zeros(degree+1)
    noise2 = 1-rho*rho
    for r in range(2, degree+1, 2):
        beta2 = 4*ph*ph*hn[r-1]**2/r
        ell = np.arange(r+1)
        logpmf = (gammaln(r+1)-gammaln(ell+1)-gammaln(r-ell+1)
                  + ell*math.log(noise2)+(r-ell)*math.log(rho*rho))
        weights[:r+1] += beta2*np.exp(logpmf)
    ell = np.arange(degree+1)

    def ratio(a):
        terms = weights/(a+ell)**2
        return np.dot(ell, terms)/np.sum(terms)

    if fixed_a is not None:
        a = fixed_a
        d = ratio(a)
        norm = np.sum(weights/(a+ell)**2)
        covariance = np.sum(weights/(a+ell))/math.sqrt(norm)
    elif np.dot(ell, weights)/np.sum(weights) <= 1:
        a = None
        d = np.dot(ell, weights)/np.sum(weights)
        covariance = math.sqrt(np.sum(weights))
    else:
        a = brentq(lambda a: ratio(a)-1, 1e-9, 1e7)
        d = ratio(a)
        norm = np.sum(weights/(a+ell)**2)
        covariance = np.sum(weights/(a+ell))/math.sqrt(norm)
    w = rho*p+math.sqrt(noise2)*covariance
    sigma = math.sqrt(max(0, p-w*w))
    value = (2*w*ph*(2*ndtr(w*alpha/sigma)-1)
             +2*math.sqrt(2*p/math.pi)*ndtr(-alpha*math.sqrt(p)/sigma))
    if details:
        return {'bound': value, 'alpha': alpha, 'rho': rho, 'a': a,
                'conditional_derivative_energy': d, 'mask_mass': p,
                'centered_covariance': covariance, 'total_covariance': w,
                'maximum_total_Hermite_degree': degree}
    return -value


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--degree', type=int, default=200)
    args = p.parse_args()
    best = None
    for initial in ((.7,.7),(.7,.85),(.8,.95),(.6,.98)):
        result = minimize(lambda x: objective(*x,args.degree), initial,
                          method='Nelder-Mead',
                          options={'xatol':1e-8,'fatol':1e-10,'maxiter':400})
        record = objective(*result.x,args.degree,details=True)
        print(json.dumps(record), flush=True)
        if best is None or record['bound'] > best['bound']:
            best = record
    print(json.dumps({'best': best, 'certified': False}), flush=True)


if __name__ == '__main__':
    main()
