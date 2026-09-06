"""Analytic one-segment concave hull for the ternary latent envelope.

Floating values are diagnostics, not an interval certificate. The shape
and unique common-tangent characterization have a separate exact proof.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def entropy(z):
    return -sum(x*math.log(x) for x in (z,1-z) if x>0)


def xcoth_minus_one(a):
    if abs(a)<.05:
        x=a*a
        return x*(1/3+x*(-1/45+x*(2/945+x*(-1/4725+x*(2/93555)))))
    return a/math.tanh(a)-1


def threshold(a):
    return 1.5 if a==0 else a*a/(2*xcoth_minus_one(a))


def magnetization(kappa,z):
    coupling=2*kappa*z
    if coupling<=1:return 0.
    return brentq(lambda a:xcoth_minus_one(a)-(coupling-1),0,coupling,xtol=1e-13)


def g(kappa,z):
    a=magnetization(kappa,z)
    if a==0:return entropy(z)+z*math.log(2)
    log2cosh=float(np.logaddexp(a,-a))
    return entropy(z)+z*log2cosh-a*a/(4*kappa)


def tangent(kappa):
    if kappa<=1.5:return None
    maximum=magnetization(kappa,1.)
    inflection=brentq(lambda a:threshold(a)-kappa,0,maximum,xtol=1e-13)
    def data(a):
        e=4*kappa*math.sinh(a)/a-2*math.cosh(a)
        return e,math.log((e+2*math.cosh(a))/(e+2))-a*a/(4*kappa)
    a=brentq(lambda a:data(a)[1],inflection,maximum,xtol=2e-13)
    e,_=data(a)
    low=2/(e+2);high=a/(2*kappa*math.tanh(a));slope=math.log(e)
    return {'magnetization_field':a,'lower_endpoint':low,'upper_endpoint':high,
       'slope':slope,'grand_pressure':math.log1p(2/e),'inflection_field':inflection}


def hull(kappa,p):
    data=tangent(kappa)
    if data is None or p<data['lower_endpoint'] or p>data['upper_endpoint']:
        return g(kappa,p),data
    return data['grand_pressure']+p*data['slope'],data


def envelope_value(p,t,kappa):
    if kappa<=0:return -math.inf
    r=1-p*kappa/t
    if not 0<=r<1:return -math.inf
    value,_=hull(kappa,p)
    return -entropy(p)+t*(1-math.sqrt(p))-p*kappa+.25*math.log1p(-r*r)+value


def diagnostic(p,t,grid):
    rho=(math.sqrt(1+16*t*t)-1)/(4*t)
    gaussian=p*math.log(2)+t*(1-math.sqrt(p))-t*(1-rho)+.25*math.log1p(-rho*rho)
    xs=np.linspace(.5,t/p,grid+1)
    vals=[envelope_value(p,t,float(x)) for x in xs]
    candidates=[(gaussian,t*(1-rho)/p)]
    for i in range(1,len(xs)-1):
        if vals[i]>=vals[i-1] and vals[i]>=vals[i+1]:
            result=minimize_scalar(lambda x:-envelope_value(p,t,x),bounds=(xs[i-1],xs[i+1]),
                method='bounded',options={'xatol':1e-12})
            candidates.append((-float(result.fun),float(result.x)))
    candidates.extend([(vals[0],float(xs[0])),(vals[-1],float(xs[-1]))])
    value,kappa=max(candidates)
    _,data=hull(kappa,p)
    return {'p':p,'t':t,'maximum_diagnostic':value,'kappa':kappa,
      'gaussian_candidate':gaussian,'candidates':candidates,'hull_data':data,
      'status':'exact analytic hull shape, floating outer optimization only; NOT a Bellman upper certificate'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--p',type=float,default=31/32)
    parser.add_argument('--t',type=float,default=4.)
    parser.add_argument('--grid',type=int,default=300)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=diagnostic(args.p,args.t,args.grid)
    print(json.dumps(result),flush=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
