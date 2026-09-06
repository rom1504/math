"""Finite full asymmetric Bellman policies. Every reported optimum is a LOWER
diagnostic, not an upper certificate. Global spin reversal and input swap are
averaged; independent input flips are deliberately NOT averaged.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, LinearConstraint
from scipy.special import logsumexp


def entropy(p, multiplicity=None):
    p=np.asarray(p)
    m=np.ones_like(p) if multiplicity is None else np.asarray(multiplicity)
    return -float(np.sum(p*np.log(np.maximum(p,1e-300)/m)))


def folded_phi(p, amplitude, t):
    """Value and simplex gradient; also a feasible-coupling primal value."""
    p=np.maximum(p,0)
    a=np.arange(len(p))*amplitude
    logk=logsumexp(np.stack([-t*(a[:,None]-a[None,:])**2,
                            -t*(a[:,None]+a[None,:])**2]),axis=0)-np.log(2)
    keep=p>1e-16
    logp=np.log(p[keep]); lk=logk[np.ix_(keep,keep)]
    u=np.zeros(keep.sum())
    for it in range(500):
        nxt=.5*(u-logsumexp(lk+logp[None,:]+u[None,:],axis=1))
        if np.max(abs(nxt-u))<3e-12:
            u=nxt;break
        u=nxt
    all_u=-logsumexp(logk[:,keep]+logp[None,:]+u[None,:],axis=1)
    all_u[keep]=u
    joint=np.exp(logp[:,None]+logp[None,:]+u[:,None]+u[None,:]+lk)
    rows=joint.sum(axis=1)
    scale=min(1.,float(np.min(p[keep]/rows)))
    feasible=scale*joint
    feasible[np.diag_indices_from(feasible)]+=np.maximum(p[keep]-scale*rows,0)
    cost=float(np.sum(feasible*(np.log(np.maximum(feasible,1e-300))
                               -logp[:,None]-logp[None,:]-lk)))
    return -float(p[keep]@u),-all_u,-cost/2,float(np.max(abs(rows-p[keep])))


def node_maps(n):
    classes=[]
    # A class records magnitudes i<=j and relative sign (+ or -).
    # Relative signs are indistinguishable if at least one input is zero.
    for i in range(n):
        for j in range(i,n):
            for phase in ([1] if i==0 else [1,-1]):
                classes.append((i,j,phase))
    count=len(classes); marginal=np.zeros((n,count))
    plus=np.zeros((2*n-1,count));minus=plus.copy(); mult=[]
    for c,(i,j,phase) in enumerate(classes):
        marginal[i,c]+=.5;marginal[j,c]+=.5
        plus[abs(i+phase*j),c]=1
        minus[abs(i-phase*j),c]=1
        mult.append(1 if j==0 else (2 if i==j else 4))
    return classes,marginal,plus,minus,np.array(mult)


class Problem:
    def __init__(self,p,t,depth):
        self.p=p;self.t=t;self.depth=depth;self.nodes=[];self.size=0
        self.layers=[]
        for d in range(depth):
            layer=[]
            for j in range(2**d):
                maps=node_maps(2**d+1); count=len(maps[0])
                node={'depth':d,'index':j,'slice':slice(self.size,self.size+count),'maps':maps}
                self.size+=count;layer.append(len(self.nodes));self.nodes.append(node)
            self.layers.append(layer)
        equations=[];rhs=[]
        for ix,node in enumerate(self.nodes):
            d=node['depth'];j=node['index'];M=node['maps'][1]
            block=np.zeros((M.shape[0],self.size));block[:,node['slice']]=M
            if d==0: target=np.array([1-p,p])
            else:
                parent=self.nodes[self.layers[d-1][j//2]]
                block[:,parent['slice']]-=parent['maps'][2+j%2]
                target=np.zeros(M.shape[0])
            equations.extend(block);rhs.extend(target)
        self.constraints=np.array(equations);self.rhs=np.array(rhs)

    def initial(self,phase=.5,diagonal_fraction=0.):
        v=np.zeros(self.size);profiles={(0,0):np.array([1-self.p,self.p])}
        for node in self.nodes:
            d=node['depth'];j=node['index'];p=profiles[d,j]
            classes,M,P,Q,mult=node['maps']; mass=[]
            theta=(1-diagonal_fraction)*np.outer(p,p)+diagonal_fraction*np.diag(p)
            for i,k,s in classes:
                base=theta[i,k]*(1 if i==k else 2)
                if i>0:base*=phase if s==1 else 1-phase
                mass.append(base)
            mass=np.array(mass);v[node['slice']]=mass
            profiles[d+1,2*j]=P@mass;profiles[d+1,2*j+1]=Q@mass
        return v

    def objective(self,v,primal=False):
        value=0.;grad=np.zeros(self.size);maxerr=0.
        for node in self.nodes:
            d=node['depth'];sl=node['slice'];classes,M,P,Q,mult=node['maps']
            q=np.maximum(v[sl],0);p=M@q;factor=2.**(-d)
            signs=np.full(len(p),2.);signs[0]=1.
            value+=factor*(.5*entropy(q,mult)-entropy(p,signs))
            grad[sl]+=factor*(.5*(-np.log(np.maximum(q,1e-300)/mult)-1)
                              +M.T@(np.log(np.maximum(p,1e-300)/signs)+1))
            if d==self.depth-1:
                for out in [P,Q]:
                    phi,g,phi_primal,err=folded_phi(out@q,1/np.sqrt(self.p*2**self.depth),self.t)
                    value+=factor*.5*(phi_primal if primal else phi)
                    grad[sl]+=factor*.5*(out.T@g)
                    maxerr=max(maxerr,err)
        return value,grad,maxerr

    def run(self,iterations,starts,warm=None):
        records=[]
        constraint=LinearConstraint(self.constraints,self.rhs,self.rhs)
        initializations=[(phase,diag,self.initial(phase,diag)) for phase,diag in starts]
        if warm is not None:initializations.insert(0,(None,None,np.array(warm)))
        for phase,diag,initial in initializations:
            answer=minimize(lambda v:(-self.objective(v)[0],-self.objective(v)[1]),initial,
                jac=True,method='SLSQP',bounds=[(0,1)]*self.size,constraints=[constraint],
                options={'ftol':2e-10,'maxiter':iterations})
            value,_,err=self.objective(answer.x)
            primal,_,_=self.objective(answer.x,primal=True)
            record={'initial_phase':phase,'initial_diagonal_fraction':diag,
              'value':value,'primal_leaf_value':primal,
              'exponent':self.p*np.log(2)+primal+self.t*(1-np.sqrt(self.p)),
              'constraint_error':float(np.max(abs(self.constraints@answer.x-self.rhs))),
              'sinkhorn_error':err,'success':bool(answer.success),'message':answer.message,
              'iterations':int(answer.nit),'policy':answer.x.tolist()}
            records.append(record)
            print(json.dumps({k:v for k,v in record.items() if k!='policy'}),flush=True)
        return {'p':self.p,'t':self.t,'depth':self.depth,'variables':self.size,
          'records':records,'status':'admissible asymmetric policy diagnostics, LOWER bounds only'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--p',type=float,default=15/16)
    ap.add_argument('--t',type=float,default=4);ap.add_argument('--depth',type=int,default=3)
    ap.add_argument('--iterations',type=int,default=250);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--starts',type=str,default='.5,0;.8,.4;1,1')
    ap.add_argument('--initial-policy',type=Path)
    args=ap.parse_args();problem=Problem(args.p,args.t,args.depth)
    starts=[tuple(map(float,item.split(','))) for item in args.starts.split(';') if item]
    warm=None
    if args.initial_policy:
        raw=json.loads(args.initial_policy.read_text())
        assert raw['p']==args.p and raw['depth']==args.depth
        warm=max(raw['records'],key=lambda r:r['primal_leaf_value'])['policy']
    result=problem.run(args.iterations,starts,warm)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
