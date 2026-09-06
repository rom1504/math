"""Bounded falsifier search for B E <= E, where E is the latent envelope.

The source-envelope approximation has a finite-grid rate-distortion dual
and an analytic reproduction-quantization allowance. Float outputs are
diagnostics; positive candidates require interval replay before banking.
"""
import argparse
import heapq
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution


LN2=np.log(2)


def ent(w):
    w=np.asarray(w);p=w[w>0];return float(-p@np.log(p))


def g(tvar):
    if tvar<=0:return 0.
    rho=4*tvar/(np.sqrt(1+16*tvar*tvar)+1)
    return -tvar*(1-rho)+.25*np.log1p(-rho*rho)


class SourceEnvelope:
    def __init__(self,absvalues,absweights,t,grid=301):
        self.x=np.array([-v for v in absvalues]+list(absvalues))
        self.p=np.array(list(absweights)+list(absweights))/2
        self.t=t;self.y=np.linspace(-max(absvalues),max(absvalues),grid)
        self.h=self.y[1]-self.y[0];self.cache={}

    def cost(self,lam):
        if lam in self.cache:return self.cache[lam]
        K=np.exp(-lam*(self.x[:,None]-self.y[None,:])**2)
        q=np.full(len(self.y),1/len(self.y))
        for step in range(3000):
            Z=K@q;scores=(self.p/Z)@K
            primal=-float(self.p@np.log(Z));dual=primal-np.log(max(scores))
            if primal-dual<2e-8:break
            q*=scores;q/=q.sum()
        # Unbiased rounding of an optimal reconstruction to neighboring
        # grid points increases MSE by at most h^2/4, and contracts I.
        answer={'lower':dual-lam*self.h*self.h/4,'grid_primal':primal,
          'grid_dual':dual,'dual_gap':primal-dual,'iterations':step+1}
        self.cache[lam]=answer;return answer

    def constant(self,lam):
        return .25*np.log(lam*(2*self.t-lam)/(self.t*self.t))

    def interval_upper(self,a,b):
        if a==0:return self.constant(b)
        return self.constant(b)-self.cost(a)['lower']

    def upper(self,tolerance=.001,max_intervals=400):
        # E>=Gaussian and >=-source entropy provide incumbent lower values.
        variance=float(self.p@(self.x*self.x))
        lower=max(g(self.t*variance),-ent(self.p))
        queue=[(-self.interval_upper(0,self.t),0.,self.t)];count=0
        while queue and count<max_intervals:
            neg,a,b=heapq.heappop(queue);upper=-neg
            if upper<=lower+tolerance:
                heapq.heappush(queue,(neg,a,b));break
            mid=(a+b)/2
            cost=self.cost(mid)
            lower=max(lower,self.constant(mid)-cost['grid_primal'])
            for left,right in [(a,mid),(mid,b)]:
                ub=self.interval_upper(left,right)
                heapq.heappush(queue,(-ub,left,right))
            count+=1
        return {'lower':lower,'upper':-queue[0][0],'intervals':count,
          'maximum_dual_gap':max(z['dual_gap'] for z in self.cache.values()),
          'grid_spacing':self.h,'qualification':'finite dual+quantization upper, floating arithmetic not interval certified'}


def child_lower(values,weights,t,rho,hfrac):
    values=np.array(values);weights=np.array(weights)
    E=float(weights@(values*values));meanabs=float(weights@values)
    H=ent(weights)+(1-weights[values==0].sum())*LN2
    lam=t*(1-rho);h=hfrac*meanabs
    field=2*lam*h*values
    logcosh=np.logaddexp(field,-field)-LN2
    binary=.25*np.log1p(-rho*rho)-lam*(E+h*h)+float(weights@logcosh)
    answers=[g(t*E),-H,binary]
    # Explicit boundary channels prevent a smooth optimizer from missing
    # coarse-sign concentration. Small magnitudes are either given a fair
    # binary label or assigned their own uninformative zero label.
    for threshold in sorted(set(values)):
        keep=(values>=threshold)&(values>0)
        q=float(weights[keep].sum());m=float(weights[keep]@values[keep])
        if q:
            answers.append(g(t*max(0,E-m*m))-q*LN2)
            answers.append(g(t*max(0,E-m*m/q))-ent([q,1-q])-q*LN2)
    return max(answers)


def pair_test(values,weights,t,iterations,seed):
    a,b=values;p,q=weights
    def objective(z):
        u=z[0];phase=z[1:4]
        pairmasses=[p-u,2*u,q-u]
        pairs=[(a,a),(a,b),(b,b)]
        outplus={};outminus={};D=0.
        theta=np.array([p-u,u,u,q-u]);product=np.array([p*p,p*q,p*q,q*q]);keep=theta>0
        D=float(np.sum(theta[keep]*np.log(theta[keep]/product[keep])))
        for mass,al,(x,y) in zip(pairmasses,phase,pairs):
            D+=mass*(LN2-ent([al,1-al]))
            add=round((x+y)/np.sqrt(2),12);sub=round(abs(x-y)/np.sqrt(2),12)
            for out,m1,m2 in [(outplus,mass*al,mass*(1-al)),(outminus,mass*(1-al),mass*al)]:
                out[add]=out.get(add,0)+m1;out[sub]=out.get(sub,0)+m2
        outputs=[]
        for out in [outplus,outminus]:
            vals=sorted(out);ws=[out[v] for v in vals];outputs.append((vals,ws))
        lp=child_lower(*outputs[0],t,z[4],z[5]);lm=child_lower(*outputs[1],t,z[6],z[7])
        return -.5*(lp+lm)+.5*D
    answer=differential_evolution(objective,[(0,min(p,q))]+[(0,1)]*3+
      [(0,1-1e-10),(0,1)]*2,popsize=12,maxiter=iterations,seed=seed,tol=1e-8,polish=True)
    candidates=[(float(answer.fun),answer.x.tolist(),'optimizer')]
    for u in [0,p*q,min(p,q)]:
        for phase in itertools.product([0,.5,1],repeat=3):
            z=[u,*phase,0,0,0,0]
            candidates.append((objective(z),z,'explicit boundary'))
    value,z,kind=min(candidates,key=lambda record:record[0])
    return {'lower':-float(value),'parameters':z,'iterations':answer.nit,'kind':kind}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--values',default='.99,1.01')
    ap.add_argument('--weight',type=float,default=.5);ap.add_argument('--t',type=float,default=4)
    ap.add_argument('--iterations',type=int,default=180);ap.add_argument('--grid',type=int,default=301)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    values=list(map(float,args.values.split(',')));weights=[args.weight,1-args.weight]
    envelope=SourceEnvelope(values,weights,args.t,args.grid).upper()
    print(json.dumps({'source_envelope':envelope}),flush=True)
    pair=pair_test(values,weights,args.t,args.iterations,9061225)
    result={'values':values,'weights':weights,'t':args.t,'source_envelope':envelope,
      'B_E_selected_lower':pair,'gap_lower_diagnostic':pair['lower']-envelope['upper'],
      'status':'bounded necessary-supersolution falsifier test, not a proof if gap nonpositive'}
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
