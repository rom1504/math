"""Exact-formula floating checks at exhaustive global width minima, n<=6."""
import json
import numpy as np
from decisive_bridge_actual_cavity_psd_checks_2026_09_07 import objects


def logsumexp(x, axis=-1):
    m = np.max(x, axis=axis, keepdims=True)
    return (m + np.log(np.exp(x-m).sum(axis=axis, keepdims=True))).squeeze(axis)


def check(n, beta, u):
    edges, chars, signs = objects(n)
    m = n // 2
    z = np.r_[np.full(m, np.sqrt((n-m)/m)),
              np.full(n-m, -np.sqrt(m/(n-m)))]
    direction = np.array([z[i]*z[j] for i,j in edges])
    tau = np.log(np.cosh(2*beta/np.sqrt(n)))
    lam = np.arccosh(np.exp(tau*(1+u*direction)))/2
    energy = (signs*lam) @ chars.T
    values = (logsumexp(energy)+logsumexp(-energy))/2
    ix = int(np.argmin(values))
    a = signs[ix]
    e = energy[ix]
    pp = np.exp(e-logsumexp(e))
    pm = np.exp(-e-logsumexp(-e))
    cp, cm = pp @ chars, pm @ chars
    product = cp*cm
    response = a*(cp-cm)/2
    slack = 1-product-2*response/np.tanh(2*lam)
    assert slack.min() > -2e-9
    t = np.tanh(lam)
    cp0 = (cp-a*t)/(1-a*t*cp)
    cm0 = (cm+a*t)/(1+a*t*cm)
    b = np.abs(cp0-cm0)
    assert np.max(a*(cp0-cm0)) < 2e-9
    q = 1-t*b-t*t*cp0*cm0
    cavity_slack = b*(1-t*t)**2/(2*t*q)
    assert np.max(np.abs(slack-cavity_slack)) < 2e-8
    costs = []
    for j in range(len(edges)):
        enew = e-2*lam[j]*a[j]*chars[:,j]
        costs.append((logsumexp(enew)+logsumexp(-enew))/2-values[ix])
    cost_slack = np.expm1(2*np.array(costs))/np.sinh(2*lam)**2
    assert np.max(np.abs(slack-cost_slack)) < 2e-8
    pmat = np.eye(n)
    for j,(v,w) in enumerate(edges):
        pmat[v,w] = pmat[w,v] = product[j]
    assert np.linalg.eigvalsh(pmat).min() > -2e-10
    psd = -tau*(z @ pmat @ z)/8
    flip = -tau*(direction @ slack)/4
    deriv = (tau*direction/(2*np.tanh(2*lam))) @ response
    assert abs(deriv-psd-flip) < 2e-9
    return dict(n=n,beta=beta,u=u,derivative=float(deriv),psd=float(psd),
                weighted_flip=float(flip),checks="PASS")


if __name__ == "__main__":
    count=0
    largest_positive=-float("inf")
    witness=None
    for n in [3,4,5,6]:
        for beta in [.2,1.,2.,4.]:
            for u in [0.,.25,.6,.9]:
                result=check(n,beta,u)
                count+=1
                if result["derivative"]>largest_positive:
                    largest_positive=result["derivative"]
                    witness=result
    print(json.dumps(dict(total_checks=count,largest_derivative_witness=witness,
                          checks="PASS")))
