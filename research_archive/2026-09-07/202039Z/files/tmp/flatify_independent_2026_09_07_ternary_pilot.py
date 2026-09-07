"""Numerical-only independent exploration of the exact ternary E formula."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution, minimize


def gain(lam, z, p):
    A = math.expm1(lam*z*z/p)
    B = 2*math.sinh(lam*z/p)**2
    if B == 0 or p*B <= A:
        return 0.
    if A == 0 or p*B-A >= (1-p)*A*B:
        return p*math.log1p(B)-math.log1p(A)
    return (p*math.log(p)+(1-p)*math.log1p(-p)-p*math.log(A)
            +math.log(B)-(1-p)*math.log(B-A))


def phi(t):
    rho=4*t/(math.sqrt(1+16*t*t)+1)
    return -t*(1-rho)+math.log1p(-rho*rho)/4


def evaluate(p,t,seed=73207):
    def f(v):
        lam,z=v
        return .25*math.log(lam*(2*t-lam)/(t*t))-lam+gain(lam,z,p)
    result=differential_evolution(lambda v:-f(v),[(.5,t),(0,1)],seed=seed,
                                  tol=1e-10,popsize=20,maxiter=1000,polish=True)
    high=-result.fun
    E=max(phi(t),high)
    return {'p':float(p),'t':float(t),'E':E,'gaussian':phi(t),'high':high,
            'argmax':list(map(float,result.x)),
            'upper':(t+p*math.log(2)+E)/(2*t*math.sqrt(p))}


def fast_evaluate(p,t):
    def f(v):
        lam,z=v
        return .25*math.log(lam*(2*t-lam)/(t*t))-lam+gain(lam,z,p)
    results=[minimize(lambda v:-f(v),start,method='Nelder-Mead',
                      bounds=[(.5,t),(0,1)],
                      options={'xatol':1e-11,'fatol':1e-13,'maxiter':1000})
             for start in [[.65*t,p],[.88*t,.99]]]
    high=max(-result.fun for result in results)
    E=max(phi(t),high)
    return {'p':float(p),'t':float(t),'E':E,'gaussian':phi(t),'high':high,
            'branches':[{'value':-r.fun,'point':r.x.tolist()} for r in results],
            'upper':(t+p*math.log(2)+E)/(2*t*math.sqrt(p))}


if __name__=='__main__':
    records=[]
    for p in [31/32,.96,.95,.94,.93,.97,.98,.99]:
        for t in [4.,4.2,4.4,4.6,4.8,5.]:
            result=evaluate(p,t)
            records.append(result)
            print(json.dumps(result),flush=True)
    print(json.dumps({'best':sorted(records,key=lambda r:r['upper'])[:10]}),flush=True)
    Path('computations/results/flatify_independent_2026_09_07_ternary_pilot.json').write_text(
        json.dumps({'status':'NUMERICAL ONLY; DE CAN MISS A LOCAL MAXIMUM',
                    'seed':73207,'records':records},indent=2)+'\n')
