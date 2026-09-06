"""Exploratory LOWER witnesses for the deep Bellman latent envelope.

The proved two-posterior-orbit reduction is in the companion artifact.
Numerical optimization here is not a global upper certificate.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import xlogy


def h(x):
    return -xlogy(x, x)-xlogy(1-x, 1-x)


def gaussian_phi(u):
    u=np.asarray(u)
    rho=4*u/(1+np.sqrt(1+16*u*u))
    return -u*(1-rho)+.25*np.log1p(-rho*rho)


def channel(p,t,x):
    z1=p*x[0];z2=p+(1-p)*x[1];s1=x[2];s2=x[3]
    lam=(z2-p)/(z2-z1) if z2>z1 else .5
    conditional_entropy=lam*(h(z1)+z1*h((1+s1)/2))+(1-lam)*(h(z2)+z2*h((1+s2)/2))
    variance=max(0.,1-(lam*z1*z1*s1*s1+(1-lam)*z2*z2*s2*s2)/p)
    info=h(p)+p*np.log(2)-conditional_entropy
    bound=p*np.log(2)-info+gaussian_phi(t*variance)+t*(1-np.sqrt(p))
    return {'exponent_lower_diagnostic':float(bound),'conditional_variance':float(variance),
            'mutual_information':float(info),'orbit_weights':[float(lam),float(1-lam)],
            'posterior_z':[float(z1),float(z2)],'posterior_bias':[float(s1),float(s2)]}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--p',type=float,default=15/16)
    ap.add_argument('--tilts',default='4,4.75,5,6');ap.add_argument('--iterations',type=int,default=700)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    records=[]
    for t in map(float,args.tilts.split(',')):
        for seed in [1,7]:
            answer=differential_evolution(lambda x:-channel(args.p,t,x)['exponent_lower_diagnostic'],
                [(0,1)]*4,seed=seed,popsize=22,maxiter=args.iterations,tol=1e-11,polish=True)
            result=channel(args.p,t,answer.x)
            result.update(p=args.p,t=t,seed=seed,iterations=answer.nit,success=bool(answer.success))
            records.append(result);print(json.dumps(result),flush=True)
    args.output.write_text(json.dumps({'records':records,'status':'selected-channel LOWER diagnostics only'},indent=2)+'\n')


if __name__=='__main__':main()
