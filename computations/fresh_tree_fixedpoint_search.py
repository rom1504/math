"""Numerical search for a subcritical Gaussian tree fixed-point certificate.

These numbers are diagnostics; a separate interval calculation is needed.
The first candidate is a normalized OU-smoothed threshold mask. The second
is the exact finite-Hermite solution of the two quadratic constraints.
"""
import json
import numpy as np
from scipy.special import ndtr
from scipy.optimize import differential_evolution, minimize, brentq


def phi(x):
    return np.exp(-np.asarray(x)**2/2)/np.sqrt(2*np.pi)


def coefficients(alpha, cutoff=120):
    h = [1., alpha]
    for k in range(1, 2*cutoff):
        h.append(alpha*h[-1]/np.sqrt(k+1)-np.sqrt(k/(k+1))*h[-2])
    return np.r_[2*ndtr(alpha)-1,
                 [-2*phi(alpha)*h[2*j-1]/np.sqrt(2*j)
                  for j in range(1, cutoff+1)]]


def tail_value(alpha, w):
    p = 2*ndtr(alpha)-1
    s = np.sqrt(max(p-w*w, 1e-30))
    return float(2*w*phi(alpha)*(2*ndtr(w*alpha/s)-1)
                 +4*np.sqrt(p)*phi(0)*ndtr(-alpha*np.sqrt(p)/s))


def ou_candidate(theta, detail=False):
    rho, alpha = theta
    coeff = coefficients(alpha)
    powers = np.arange(len(coeff))*2
    k = coeff[1:]*rho**powers[1:]
    sd = np.linalg.norm(k)
    g = np.r_[rho, np.sqrt(1-rho*rho)*k/sd]
    derivative = float(powers@(g*g))
    w = float(g@coeff)
    val = tail_value(alpha, w)
    if detail:
        return dict(rho=float(rho),alpha=float(alpha),derivative=derivative,
                    overlap=w,value=val,mean=float(g[0]))
    return val-100*max(0,derivative-.999)


def resolvent_candidate(alpha, cutoff=120, target=.999, detail=False):
    coeff = coefficients(alpha,cutoff)
    degree = np.arange(len(coeff))*2
    def calc(a):
        g = coeff/(a+degree)
        g /= np.linalg.norm(g)
        return g, float(degree@(g*g))
    # The ratio rises monotonically from zero toward the threshold influence.
    if calc(1e6)[1] <= target:
        a = float('inf')
        g = coeff/np.linalg.norm(coeff)
        d = float(degree@(g*g))
    else:
        a = brentq(lambda x: calc(x)[1]-target, 1e-9, 1e6)
        g,d = calc(a)
    w = float(g@coeff)
    val = tail_value(alpha,w)
    if detail:
        return dict(alpha=float(alpha),cutoff=cutoff,resolvent=a,
                    derivative=d,overlap=w,value=val,
                    Hermite_coefficients=list(map(float,g)))
    return val


def main():
    ou = differential_evolution(lambda x:-ou_candidate(x),
        [(.02,.999),(.1,2)],seed=20260905,maxiter=180,tol=1e-10)
    out = {"evidence":"Gaussian numerical diagnostic, not certificate",
           "OU_candidate":ou_candidate(ou.x,True),"finite_resolvents":[]}
    for degree in (2,4,8,16,32,64,120):
        opt = minimize(lambda x:-resolvent_candidate(float(x[0]),degree),
            [.8],bounds=[(.05,2)],method='Nelder-Mead',
            options={'xatol':1e-10,'fatol':1e-12})
        out['finite_resolvents'].append(resolvent_candidate(
            float(opt.x[0]),degree,detail=True))
    print(json.dumps(out,indent=2),flush=True)


if __name__ == '__main__':
    main()
