"""Necessary supersolution test at a binary source's temperature transition.

Uses the independently proved ternary posterior hull for child LOWER values
and an independent rate-distortion dual for the source UPPER diagnostic.
"""
import argparse
import heapq
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar


def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


hull=load('hull','continued_feedback_ternary_concave_envelope_2026_09_06.py')
dual=load('dual','continued_convergence_latent_supersolution_test_2026_09_06.py')


def child(p,t):
    if p==0:return 0.,None
    result=hull.diagnostic(p,t,55)
    candidates=[(result['maximum_diagnostic'],result['kappa'])]
    maximum=t/p
    if maximum>.5:
        # The rare-information branch may peak arbitrarily close to the
        # endpoint and need not be a local maximum on a uniform grid.
        answer=minimize_scalar(lambda k:-hull.envelope_value(p,t,k),
          bounds=(max(.5,maximum-.15),maximum),method='bounded',options={'xatol':2e-13})
        candidates.append((-float(answer.fun),float(answer.x)))
    value,kappa=max(candidates)
    actual=value-p*np.log(2)-t*(1-np.sqrt(p))
    return actual,kappa


def test(t,grid):
    def constant(lam):return .25*np.log(lam*(2*t-lam)/(t*t))
    def cost(lam):return np.log(2)+lam-hull.g(lam,1.)
    lower=max(dual.g(t),-np.log(2))
    if t>.5:
        candidate=minimize_scalar(lambda k:-(constant(k)-cost(k)),bounds=(max(.5,t-.15),t),
          method='bounded',options={'xatol':2e-13})
        lower=max(lower,-float(candidate.fun))
    queue=[(0.,0.,t)];count=0
    while queue and count<1200:
        negative,a,b=heapq.heappop(queue)
        if -negative<=lower+2e-6:heapq.heappush(queue,(negative,a,b));break
        mid=(a+b)/2;lower=max(lower,constant(mid)-cost(mid))
        for left,right in [(a,mid),(mid,b)]:
            upper=constant(right)-(cost(left) if left else 0)
            heapq.heappush(queue,(-upper,left,right))
        count+=1
    source={'lower':lower,'upper':-queue[0][0],'intervals':count,
      'qualification':'exact binary Curie-Weiss cost, floating monotone interval cover; not interval arithmetic'}
    def pair(e):
        plus,kp=child(e,2*e*t);minus,km=child(1-e,2*(1-e)*t)
        return .5*(plus+minus)-.5*(np.log(2)-hull.entropy(e)),kp,km
    xs=np.unique(np.concatenate([np.linspace(0,.5,grid+1),np.geomspace(1e-12,.05,90)]))
    records=[(pair(float(e))[0],float(e)) for e in xs]
    for i in range(1,len(xs)-1):
        if records[i][0]>=records[i-1][0] and records[i][0]>=records[i+1][0]:
            answer=minimize_scalar(lambda e:-pair(e)[0],bounds=(xs[i-1],xs[i+1]),
              method='bounded',options={'xatol':1e-9})
            records.append((-float(answer.fun),float(answer.x)))
    value,e=max(records);_,kp,km=pair(e)
    return {'t':t,'source_envelope':source,'B_E_lower':value,'equal_sign_fraction':e,
      'child_kappas':[kp,km],'gap_lower_diagnostic':value-source['upper'],
      'grid_records':records,'scope':'positive gap would require interval replay; negative gap is not supersolution proof'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tilts',default='2.5,2.75,3,3.25,4')
    ap.add_argument('--grid',type=int,default=30);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();records=[]
    for t in map(float,args.tilts.split(',')):
        record=test(t,args.grid);records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k!='grid_records'}),flush=True)
    args.output.write_text(json.dumps({'records':records},indent=2)+'\n')


if __name__=='__main__':main()
